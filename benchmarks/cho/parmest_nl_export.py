"""
Simplified CHO Parameter Estimation with NL File Export
========================================================

Hey John! Welcome to the CHO cell culture parameter estimation script.
This is a self-contained Pyomo DAE example that Victor put together so you
can kick the tires on different NLP solvers using the exported .nl file.
No external local imports needed — just Pyomo, NumPy, SciPy, and friends.

What this script does (the elevator pitch):
    1. Generates 10 synthetic CHO cell culture batches with randomized
       initial conditions and lognormal measurement noise (because real
       bioprocess data is never clean — surprise!).
    2. Builds ONE large monolithic Pyomo model with shared kinetic parameters
       across all batches. Each batch lives in its own Pyomo Block.
    3. Discretizes using orthogonal collocation on finite elements
       (60 FE x 3 CP x 10 batches x 12 state+derivative vars = ~21,700 variables).
       Yes, it's big. That's the point.
    4. Exports a .nl file so you can throw your favorite solver at it.
       IPOPT, KNITRO, CONOPT, or whatever you've got lying around in
       your solver drawer. We won't judge.
    5. Solves with IPOPT (ma57 linear solver) and prints a parameter
       recovery table comparing true vs estimated values.
    6. Generates two figures:
       - Training fit (10 batches used in estimation)
       - Validation fit (5 NEW held-out batches, forward-simulated with
         the estimated parameters — the real test of whether we actually
         learned anything or just overfit like a neural network at a
         buffet).

The CHO Model (6-state batch culture):
    Xv  — Viable cell density [10^6 cells/mL]
    GLC — Glucose [mM]  (the cells' favorite snack)
    GLN — Glutamine [mM] (the other favorite snack)
    LAC — Lactate [mM]   (metabolic waste — the cells' regret)
    AMM — Ammonia [mM]   (more regret)
    B1  — Product/mAb [a.u.] (what we actually care about)

Kinetics follow Monod growth with substrate inhibition (Xing 2010,
Kornecki & Strube 2018). Death rate is driven by toxic byproducts
because even cells have consequences for their lifestyle choices.

Usage:
    python parmest_nl_export.py

The .nl file lands in nl_export_results/cho_parmest.nl.
IPOPT log saved to nl_export_results/ipopt.log.

P.S. — If IPOPT struggles, remember: it's not a bug, it's a
nonconvex NLP. Welcome to our world.

References
----------
- Kornecki & Strube (2018), Bioengineering, 5(1), 25.
- Xing et al. (2010), Biotechnol. Prog., 26(1), 282-291.
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend: save to disk, no GUI windows
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

import pyomo.environ as pyo
import pyomo.dae as dae
from scipy.integrate import solve_ivp


# ══════════════════════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════════════════════
# (Tweak these knobs to your heart's content, John.)

IPOPT_PATH = None         # Use system default ipopt (None = whatever's on PATH)
T_FINAL = 240.0           # Simulation horizon [hours] (10 days of cell culture)
N_BATCHES = 10            # Training batches for parameter estimation
N_VALIDATION = 5          # Held-out batches for validation (the honesty check)
NFE = 60                  # Finite elements per batch (controls discretization granularity)
NCP = 3                   # Collocation points per element (Lagrange-Radau)
N_MEAS = 30               # Measurement time points per batch (sparse sampling)
NOISE_LEVEL = 0.03        # Lognormal noise sigma (3% CV — optimistic for real bioprocess data)
RNG_SEED = 42             # The answer to life, the universe, and reproducibility
OUTDIR = "nl_export_results"

# State variable names (order matters — matches the ODE system)
STATES = ["Xv", "GLC", "GLN", "LAC", "AMM", "B1"]

# Parameters to estimate (the 12 unknowns we're hunting for)
THETA_NAMES = [
    "K_glc", "K_gln", "KI_amm", "KI_lac", "KD_amm", "KD_lac",
    "m_glc", "a1", "a2", "d_gln", "r_amm", "Q_B1",
]

# Bounds on estimated parameters — wide enough to be fair, tight enough
# to keep IPOPT from wandering off into nonsense-land.
THETA_BOUNDS = {
    "K_glc":  (0.01,  5.0),    # [mM] Monod half-saturation, glucose
    "K_gln":  (0.001, 0.5),    # [mM] Monod half-saturation, glutamine
    "KI_amm": (1.0,   20.0),   # [mM] Ammonia inhibition constant
    "KI_lac": (10.0,  100.0),  # [mM] Lactate inhibition constant
    "KD_amm": (1.0,   20.0),   # [mM] Ammonia death constant
    "KD_lac": (10.0,  100.0),  # [mM] Lactate death constant
    "m_glc":  (1e-6,  1e-3),   # [mmol/(10^6 cells·h)] Glucose maintenance
    "a1":     (1e-7,  1e-4),   # [mmol/(10^6 cells·h)] Max glutamine maintenance
    "a2":     (0.1,   10.0),   # [mM] Glutamine maintenance half-saturation
    "d_gln":  (1e-4,  0.1),    # [1/h] Glutamine spontaneous degradation
    "r_amm":  (1e-7,  1e-4),   # [mmol/(10^6 cells·h)] Ammonia consumption rate
    "Q_B1":   (1e-8,  1e-3),   # [a.u./(10^6 cells·h)] Specific mAb productivity
}


# ══════════════════════════════════════════════════════════════════════════════
# True Parameters
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Params:
    """
    Container for CHO cell culture kinetic parameters.

    These are the "ground truth" values used to generate synthetic data.
    The parameter estimation will try to recover them from noisy measurements.
    (Spoiler: it usually does pretty well, unless you crank up the noise.)

    All units are consistent with Xv expressed as 10^6 cells/mL.

    Parameters
    ----------
    mu_max : float
        Maximum specific growth rate [1/h]. How fast the cells can possibly
        grow when everything is perfect (which it never is).
    k_d : float
        Base death rate coefficient [1/h]. Even in paradise, some cells
        don't make it.
    Y_Xv_glc : float
        Biomass yield on glucose [(10^6 cells/mL)/mmol]. How many cells
        you get per unit of sugar consumed.
    Y_Xv_gln : float
        Biomass yield on glutamine [(10^6 cells/mL)/mmol]. Same idea,
        different substrate.
    Y_lac_glc : float
        Lactate yield on glucose [mol/mol]. The price of glycolysis.
    Y_amm_gln : float
        Ammonia yield on glutamine [mol/mol]. The price of glutaminolysis.
    K_glc : float
        Monod half-saturation constant for glucose [mM]. Concentration at
        which growth is half its maximum w.r.t. glucose.
    K_gln : float
        Monod half-saturation constant for glutamine [mM].
    KI_amm : float
        Ammonia inhibition constant [mM]. Higher = more tolerant cells.
    KI_lac : float
        Lactate inhibition constant [mM]. Higher = more tolerant cells.
    KD_amm : float
        Ammonia death rate constant [mM]. Drives toxicity-induced death.
    KD_lac : float
        Lactate death rate constant [mM]. Drives toxicity-induced death.
    a2 : float
        Half-saturation for glutamine maintenance term [mM].
    d_gln : float
        Spontaneous glutamine degradation rate [1/h]. Glutamine is unstable
        in aqueous solution — it doesn't need cells to disappear.
    m_glc : float
        Maintenance glucose consumption [mmol/(10^6 cells·h)]. Cells gotta
        keep the lights on even when they're not dividing.
    a1 : float
        Maximum glutamine maintenance rate [mmol/(10^6 cells·h)].
    r_amm : float
        Ammonia consumption/detoxification rate [mmol/(10^6 cells·h)].
        The cells' attempt at cleaning up after themselves.
    Q_B1 : float
        Specific mAb productivity [a.u./(10^6 cells·h)]. The whole reason
        we're growing these cells in the first place.
    """
    # Growth / death [1/h]
    mu_max:    float = 0.039
    k_d:       float = 0.004

    # Biomass yields [(10^6 cells/mL) per mmol substrate]
    Y_Xv_glc:  float = 357.0
    Y_Xv_gln:  float = 974.0

    # By-product yields [mol/mol]
    Y_lac_glc: float = 0.70
    Y_amm_gln: float = 0.67

    # Monod / inhibition / death constants [mM]
    K_glc:     float = 1.00
    K_gln:     float = 0.047
    KI_amm:    float = 6.51
    KI_lac:    float = 43.0
    KD_amm:    float = 6.51
    KD_lac:    float = 45.8

    # Maintenance & degradation
    a2:        float = 2.1       # [mM]
    d_gln:     float = 7.2e-3   # [1/h]
    m_glc:     float = 6.92e-5  # [mmol/(10^6 cells·h)]
    a1:        float = 3.2e-6   # [mmol/(10^6 cells·h)]
    r_amm:     float = 6.3e-6   # [mmol/(10^6 cells·h)]
    Q_B1:      float = 4.0e-6   # [a.u./(10^6 cells·h)]


# ══════════════════════════════════════════════════════════════════════════════
# ODE System
# ══════════════════════════════════════════════════════════════════════════════

def rhs(t: float, y: np.ndarray, p: Params) -> np.ndarray:
    """
    Right-hand side of the 6-state CHO batch culture ODE system.

    This is the "ground truth" continuous-time model used to generate
    synthetic data via scipy's solve_ivp. The same equations appear
    in the Pyomo model as algebraic constraints after collocation.

    The growth rate mu follows double-Monod kinetics (glucose + glutamine)
    with product inhibition from lactate and ammonia. Death rate mu_d is
    driven by toxic byproduct accumulation. Think of it as the cells
    partying too hard and suffering the consequences.

    Parameters
    ----------
    t : float
        Current time [h]. Not explicitly used (autonomous system), but
        required by solve_ivp's interface. The cells don't know what
        time it is either.
    y : np.ndarray, shape (6,)
        State vector: [Xv, GLC, GLN, LAC, AMM, B1].
    p : Params
        Kinetic parameter set.

    Returns
    -------
    np.ndarray, shape (6,)
        Time derivatives: [dXv/dt, dGLC/dt, dGLN/dt, dLAC/dt, dAMM/dt, dB1/dt].

    Notes
    -----
    - Glucose maintenance has a smooth availability limiter (GLC/(K_glc+GLC))
      to prevent negative glucose consumption when GLC ~ 0.
    - mAb production (dB1/dt) is inversely related to growth rate, because
      cells that are busy dividing are too distracted to make antibodies.
      (Relatable, honestly.)
    """
    Xv, GLC, GLN, LAC, AMM, B1 = y

    # Specific growth rate: double-Monod with inhibition
    mu = (p.mu_max
          * (GLC / (p.K_glc + GLC))
          * (GLN / (p.K_gln + GLN))
          * (p.KI_lac / (p.KI_lac + LAC))
          * (p.KI_amm / (p.KI_amm + AMM)))

    # Specific death rate: driven by toxic byproducts
    mu_d = p.k_d * (LAC / (p.KD_lac + LAC)) * (AMM / (p.KD_amm + AMM))

    # Smooth glucose availability (prevents negative consumption near depletion)
    glc_avail = GLC / (p.K_glc + GLC)

    # Glutamine maintenance (saturable)
    m_gln = p.a1 * GLN / (p.a2 + GLN)

    # Net growth
    net = mu - mu_d

    # ODEs
    dXv  = net * Xv
    dGLC = -((net / p.Y_Xv_glc) + p.m_glc * glc_avail) * Xv
    dGLN = -((net / p.Y_Xv_gln) + m_gln) * Xv - p.d_gln * GLN
    dLAC = p.Y_lac_glc * ((net / p.Y_Xv_glc) + p.m_glc * glc_avail) * Xv
    dAMM = p.Y_amm_gln * (net / p.Y_Xv_gln) * Xv - p.r_amm * Xv + p.d_gln * GLN
    dB1  = p.Q_B1 * Xv * (1.0 - mu / p.mu_max)

    return np.array([dXv, dGLC, dGLN, dLAC, dAMM, dB1])


def simulate(p: Params, tgrid: np.ndarray, y0: np.ndarray) -> pd.DataFrame:
    """
    Forward-simulate the CHO ODE system using scipy's solve_ivp.

    Uses LSODA (auto-switching stiff/non-stiff) with tight tolerances
    and dense output for accurate interpolation at arbitrary time points.
    Basically, we let the ODE solver do the hard work so we can focus
    on the fun part (parameter estimation).

    Parameters
    ----------
    p : Params
        Kinetic parameters for the simulation.
    tgrid : np.ndarray, shape (N,)
        Time points at which to evaluate the solution [h].
        Must be monotonically increasing.
    y0 : np.ndarray, shape (6,)
        Initial conditions: [Xv0, GLC0, GLN0, LAC0, AMM0, B10].

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ["Xv", "GLC", "GLN", "LAC", "AMM", "B1", "time"],
        evaluated at each point in tgrid.

    Raises
    ------
    RuntimeError
        If the ODE solver fails to integrate (e.g., stiffness issues,
        which shouldn't happen with LSODA unless the parameters are
        truly bonkers).
    """
    sol = solve_ivp(lambda t, y: rhs(t, y, p), (tgrid[0], tgrid[-1]),
                    y0, dense_output=True, method="LSODA",
                    rtol=1e-10, atol=1e-12, max_step=0.5)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = sol.sol(tgrid).T
    return pd.DataFrame(Y, columns=STATES).assign(time=tgrid)


# ══════════════════════════════════════════════════════════════════════════════
# Synthetic Data Generation
# ══════════════════════════════════════════════════════════════════════════════

def generate_batches(
    p_true: Params,
    n_batches: int,
    n_meas: int,
    noise_sigma: float,
    seed: int,
) -> List[Dict]:
    """
    Generate synthetic CHO batch culture datasets with randomized IVPs and noise.

    Each batch gets a different initial condition drawn uniformly around a
    nominal operating point. The "measurements" are the true trajectories
    corrupted with multiplicative lognormal noise — which guarantees
    nonnegativity (because negative cell counts would be... concerning).

    John, if you want to make the estimation harder, crank up noise_sigma.
    If you want to make it trivially easy, set it to 0. We won't tell anyone.

    Parameters
    ----------
    p_true : Params
        Ground-truth parameters used to simulate the "real" cell culture.
    n_batches : int
        Number of synthetic batches to generate.
    n_meas : int
        Number of uniformly-spaced measurement time points per batch.
    noise_sigma : float
        Standard deviation of the lognormal noise. A value of 0.03 means
        ~3% coefficient of variation. Real bioprocess data is often worse,
        but let's be kind to the optimizer.
    seed : int
        Random seed for reproducibility. Change this to get different
        random batches (but the same seed always gives the same data).

    Returns
    -------
    list of dict
        Each dict contains:
        - "id" : int — Batch index.
        - "y0" : np.ndarray, shape (6,) — Initial conditions used.
        - "truth" : pd.DataFrame — Noise-free simulated trajectory.
        - "meas" : pd.DataFrame — Noisy "measured" trajectory.

    Notes
    -----
    Nominal IVPs (the "average" bioreactor inoculation):
        Xv=4.0, GLC=35.0, GLN=4.0, LAC=0.5, AMM=0.3, B1=0.0

    The spread on each state is ±[2.0, 10.0, 2.0, 0.5, 0.3, 0.0], so
    you get a decent variety of starting conditions. B1 always starts at
    zero because you can't have product before the cells start working.
    That would violate thermodynamics AND common sense.
    """
    rng = np.random.default_rng(seed)
    tgrid = np.linspace(0, T_FINAL, n_meas)

    # Nominal initial conditions and their half-ranges
    nominal = np.array([4.0, 35.0, 4.0, 0.5, 0.3, 0.0])
    spreads = np.array([2.0, 10.0, 2.0, 0.5, 0.3, 0.0])

    batches = []
    for i in range(n_batches):
        # Random perturbation around nominal
        y0 = nominal + rng.uniform(-1, 1, size=6) * spreads
        y0 = np.maximum(y0, 0.01)  # Concentrations must be positive
        y0[5] = 0.0                # B1 always starts at zero (no free antibodies)

        truth = simulate(p_true, tgrid, y0)

        # Apply multiplicative lognormal noise: y_meas = y_true * exp(N(0, sigma))
        meas = truth.copy()
        for s in STATES:
            eps = rng.normal(0, noise_sigma, size=len(tgrid))
            meas[s] = truth[s].values * np.exp(eps)
        meas[STATES] = meas[STATES].clip(lower=0.0)

        batches.append({"id": i, "y0": y0, "truth": truth, "meas": meas})

    return batches


# ══════════════════════════════════════════════════════════════════════════════
# Pyomo Model Construction
# ══════════════════════════════════════════════════════════════════════════════

def build_full_model(
    batches: List[Dict],
    nfe: int = NFE,
    ncp: int = NCP,
) -> pyo.ConcreteModel:
    """
    Build a monolithic Pyomo DAE model for multi-batch parameter estimation.

    This is where the magic happens. We create ONE big optimization problem
    where 12 kinetic parameters are shared across all batches, and each
    batch has its own set of state trajectories discretized via orthogonal
    collocation on finite elements.

    The resulting NLP has ~21,700 variables and ~21,660 constraints.
    It's big enough to be interesting but small enough that IPOPT handles
    it in reasonable time. Perfect for benchmarking your favorite solver,
    John. (No pressure.)

    Architecture:
        m (ConcreteModel)
        ├── Shared Vars: K_glc, K_gln, KI_amm, ... (12 estimated parameters)
        ├── Fixed Params: mu_max, k_d, Y_Xv_glc, ... (6 known constants)
        └── m.batch[0..9] (Blocks)
            ├── ContinuousSet: b.t
            ├── State Vars: b.Xv[t], b.GLC[t], ..., b.B1[t]
            ├── DerivativeVars: b.dXv[t], b.dGLC[t], ...
            ├── Expressions: b.mu[t], b.mu_d[t], b.glc_avail[t], b.m_gln_expr[t]
            └── ODE Constraints: b.ode_Xv[t], b.ode_GLC[t], ...

    Objective: minimize SSE = sum over batches, measurement times, and states
    of (y_model - y_meas)^2. Classic least squares, no weights. If you want
    weighted LSQ, that's a homework exercise.

    Parameters
    ----------
    batches : list of dict
        Synthetic batch data from ``generate_batches()``. Each dict must
        contain keys: "meas" (noisy measurements), "truth" (for initialization),
        and "y0" (initial conditions).
    nfe : int, default=60
        Number of finite elements for collocation discretization per batch.
        More elements = finer time resolution = more variables = more fun
        for the solver.
    ncp : int, default=3
        Number of collocation points per finite element (Lagrange-Radau).
        3 is standard. If you set this to 1, things get less accurate.
        If you set it to 6, things get less... stable. 3 is the sweet spot.

    Returns
    -------
    pyo.ConcreteModel
        Fully constructed, discretized Pyomo model ready for .write("*.nl")
        or direct solving. States are initialized from truth trajectories
        (via linear interpolation onto the collocation grid) to give the
        solver a fighting chance.

    Notes
    -----
    - Estimated parameters are initialized at 80% of their true values.
      This simulates the common scenario of "we have a rough idea but
      need data to nail it down." If you initialize at the true values,
      of course it converges instantly. That's cheating. Don't do that.
      (Unless you're debugging. Then it's "verification.")
    - The collocation discretization uses Pyomo's dae.collocation
      transformation (Lagrange-Radau polynomials). This converts each
      ODE into a set of algebraic constraints — the "simultaneous"
      approach to dynamic optimization. No sequential integration needed.
    """
    m = pyo.ConcreteModel("CHO_parmest_10batches")

    # ── Shared estimated parameters (the unknowns) ──────────────────────────
    # These live on the top-level model so they're shared across all batches.
    # Initialized at 80% of truth because in real life you never know the
    # exact values. (If you did, you wouldn't need parameter estimation.)
    for name in THETA_NAMES:
        true_val = getattr(Params(), name)
        setattr(m, name, pyo.Var(
            bounds=THETA_BOUNDS[name],
            initialize=true_val * 0.8,
            domain=pyo.NonNegativeReals,
        ))

    # ── Fixed (known) parameters ────────────────────────────────────────────
    # These are NOT estimated — they're assumed known from prior literature.
    # In practice, you might want to estimate some of these too, but we're
    # keeping the problem tractable. You're welcome.
    p = Params()
    m.mu_max    = pyo.Param(initialize=p.mu_max)     # [1/h]
    m.k_d       = pyo.Param(initialize=p.k_d)        # [1/h]
    m.Y_Xv_glc  = pyo.Param(initialize=p.Y_Xv_glc)  # [(10^6 cells/mL)/mmol]
    m.Y_Xv_gln  = pyo.Param(initialize=p.Y_Xv_gln)  # [(10^6 cells/mL)/mmol]
    m.Y_lac_glc = pyo.Param(initialize=p.Y_lac_glc)  # [mol/mol]
    m.Y_amm_gln = pyo.Param(initialize=p.Y_amm_gln)  # [mol/mol]

    # ── One Block per batch ─────────────────────────────────────────────────
    # Each batch is an independent experiment with its own IVP and time
    # discretization, but sharing the same kinetic parameters. This is
    # the multi-experiment structure that makes parameter estimation robust.
    m.BATCHES = pyo.Set(initialize=range(len(batches)))
    m.batch = pyo.Block(m.BATCHES)

    for idx in m.BATCHES:
        b = m.batch[idx]
        batch = batches[idx]
        meas = batch["meas"]
        truth = batch["truth"]
        y0 = batch["y0"]
        timepoints = list(meas["time"].values)

        # ── Time set ────────────────────────────────────────────────────────
        # ContinuousSet is initialized with measurement times but will be
        # augmented with collocation points after discretization.
        b.t = dae.ContinuousSet(initialize=timepoints)

        # ── State variables ─────────────────────────────────────────────────
        # All states are nonneg (you can't have negative cells or negative
        # glucose... unless your model is very wrong).
        b.Xv  = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[0])
        b.GLC = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[1])
        b.GLN = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[2])
        b.LAC = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[3])
        b.AMM = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[4])
        b.B1  = pyo.Var(b.t, domain=pyo.NonNegativeReals, initialize=y0[5])

        # ── Derivative variables ────────────────────────────────────────────
        # Pyomo.DAE creates these as algebraic variables. After collocation,
        # dx/dt = f(x, theta) becomes a big system of algebraic equations.
        b.dXv  = dae.DerivativeVar(b.Xv,  wrt=b.t)
        b.dGLC = dae.DerivativeVar(b.GLC, wrt=b.t)
        b.dGLN = dae.DerivativeVar(b.GLN, wrt=b.t)
        b.dLAC = dae.DerivativeVar(b.LAC, wrt=b.t)
        b.dAMM = dae.DerivativeVar(b.AMM, wrt=b.t)
        b.dB1  = dae.DerivativeVar(b.B1,  wrt=b.t)

        # ── Fix initial conditions ──────────────────────────────────────────
        # These are "known" from the experiment — we measured the inoculum.
        t0 = timepoints[0]
        b.Xv[t0].fix(float(y0[0]))
        b.GLC[t0].fix(float(y0[1]))
        b.GLN[t0].fix(float(y0[2]))
        b.LAC[t0].fix(float(y0[3]))
        b.AMM[t0].fix(float(y0[4]))
        b.B1[t0].fix(float(y0[5]))

        # ── Kinetic expressions ─────────────────────────────────────────────
        # These reference the SHARED parameters on `m` (not `b`), which is
        # the key trick for multi-batch estimation.

        # Growth rate: double Monod × double inhibition
        # mu = mu_max * (GLC/(K_glc+GLC)) * (GLN/(K_gln+GLN))
        #            * (KI_lac/(KI_lac+LAC)) * (KI_amm/(KI_amm+AMM))
        b.mu = pyo.Expression(b.t, rule=lambda b, t:
            m.mu_max
            * (b.GLC[t] / (m.K_glc + b.GLC[t]))
            * (b.GLN[t] / (m.K_gln + b.GLN[t]))
            * (m.KI_lac / (m.KI_lac + b.LAC[t]))
            * (m.KI_amm / (m.KI_amm + b.AMM[t]))
        )

        # Death rate: driven by toxic byproduct accumulation
        # mu_d = k_d * (LAC/(KD_lac+LAC)) * (AMM/(KD_amm+AMM))
        b.mu_d = pyo.Expression(b.t, rule=lambda b, t:
            m.k_d
            * (b.LAC[t] / (m.KD_lac + b.LAC[t]))
            * (b.AMM[t] / (m.KD_amm + b.AMM[t]))
        )

        # Glucose availability limiter (smooth, ∈ [0,1])
        # Prevents the maintenance term from consuming glucose that isn't there.
        b.glc_avail = pyo.Expression(b.t, rule=lambda b, t:
            b.GLC[t] / (m.K_glc + b.GLC[t])
        )

        # Glutamine maintenance (saturable Michaelis-Menten-like)
        b.m_gln_expr = pyo.Expression(b.t, rule=lambda b, t:
            m.a1 * (b.GLN[t] / (m.a2 + b.GLN[t]))
        )

        # ── ODE constraints ─────────────────────────────────────────────────
        # Each ODE is written as: dState/dt == f(states, params)
        # After collocation, these become algebraic equality constraints.

        @b.Constraint(b.t)
        def ode_Xv(b, t):
            """Cell growth: dXv/dt = (mu - mu_d) * Xv"""
            return b.dXv[t] == (b.mu[t] - b.mu_d[t]) * b.Xv[t]

        @b.Constraint(b.t)
        def ode_GLC(b, t):
            """Glucose consumption: growth-associated + maintenance"""
            return b.dGLC[t] == -(((b.mu[t] - b.mu_d[t]) / m.Y_Xv_glc)
                                  + m.m_glc * b.glc_avail[t]) * b.Xv[t]

        @b.Constraint(b.t)
        def ode_GLN(b, t):
            """Glutamine consumption + spontaneous degradation"""
            return b.dGLN[t] == -(((b.mu[t] - b.mu_d[t]) / m.Y_Xv_gln)
                                  + b.m_gln_expr[t]) * b.Xv[t] - m.d_gln * b.GLN[t]

        @b.Constraint(b.t)
        def ode_LAC(b, t):
            """Lactate production: stoichiometric to glucose consumption"""
            return b.dLAC[t] == m.Y_lac_glc * (((b.mu[t] - b.mu_d[t]) / m.Y_Xv_glc)
                                                + m.m_glc * b.glc_avail[t]) * b.Xv[t]

        @b.Constraint(b.t)
        def ode_AMM(b, t):
            """Ammonia: production from glutamine + degradation - consumption"""
            return b.dAMM[t] == (m.Y_amm_gln * ((b.mu[t] - b.mu_d[t]) / m.Y_Xv_gln) * b.Xv[t]
                                 - m.r_amm * b.Xv[t] + m.d_gln * b.GLN[t])

        @b.Constraint(b.t)
        def ode_B1(b, t):
            """mAb production: inversely related to growth (non-growth associated)"""
            return b.dB1[t] == m.Q_B1 * b.Xv[t] * (1.0 - b.mu[t] / m.mu_max)

        # ── Collocation discretization ──────────────────────────────────────
        # This is the Biegler-style simultaneous approach: convert the DAE
        # into a (very large) NLP by representing state trajectories as
        # piecewise polynomials on each finite element.
        disc = pyo.TransformationFactory("dae.collocation")
        disc.apply_to(b, nfe=nfe, ncp=ncp, wrt=b.t)

        # ── Initialize states from truth ────────────────────────────────────
        # Linear interpolation of the true trajectory onto the collocation
        # grid. This gives IPOPT a warm start so it doesn't wander around
        # in the dark. (Good initialization is half the battle in NLP.)
        t_grid_model = sorted(pyo.value(t) for t in b.t)
        t_src = truth["time"].values
        for s in STATES:
            y_src = truth[s].values
            vals = np.interp(t_grid_model, t_src, y_src)
            vals = np.maximum(vals, 0.0)
            for t_val, v in zip(t_grid_model, vals):
                getattr(b, s)[t_val].set_value(float(v))

    # ── SSE Objective ───────────────────────────────────────────────────────
    # Sum of squared errors across all batches, all measurement times,
    # and all 6 states. No weighting (all states treated equally).
    # If your states have wildly different magnitudes, you might want to
    # scale this. But for this example, it works fine.
    obj_expr = 0.0
    for idx in m.BATCHES:
        b = m.batch[idx]
        meas = batches[idx]["meas"]
        for _, row in meas.iterrows():
            t = float(row["time"])
            # Match measurement time to nearest collocation point
            t_set = sorted(pyo.value(tt) for tt in b.t)
            t_closest = min(t_set, key=lambda tc: abs(tc - t))
            if abs(t_closest - t) < 1e-6:
                for s in STATES:
                    y_meas = float(row[s])
                    y_model = getattr(b, s)[t_closest]
                    obj_expr += (y_model - y_meas) ** 2

    m.obj = pyo.Objective(expr=obj_expr, sense=pyo.minimize)

    return m


# ══════════════════════════════════════════════════════════════════════════════
# Main Pipeline
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """
    Run the full parameter estimation pipeline.

    This is the main entry point. It orchestrates everything:
    data generation, model construction, NL export, solving, parameter
    comparison, training plots, and validation plots.

    If you're reading this, John, here's the TL;DR of what happens
    when you hit "python parmest_nl_export.py":

    Step 1: Generate 10 synthetic batches (different IVPs, 3% noise)
    Step 2: Build a ~21,700-variable Pyomo model
    Step 3: Export .nl file (for your solver benchmarking pleasure)
    Step 4: Solve with IPOPT (MA57 linear solver)
    Step 5: Print parameter recovery table (true vs estimated)
    Step 6: Plot training fit (10 batches)
    Step 7: Generate 5 NEW validation batches, forward-simulate with
            estimated params, and plot (the moment of truth)

    Outputs (all in nl_export_results/):
        - cho_parmest.nl   : The NL file. Feed this to any AMPL solver.
        - ipopt.log        : Full IPOPT iteration log.
        - parmest_fit.png  : Training fit plot.
        - validation_fit.png : Validation plot (held-out batches).
    """
    os.makedirs(OUTDIR, exist_ok=True)
    p_true = Params()

    # ── Step 1: Generate synthetic data ─────────────────────────────────────
    print(f"Generating {N_BATCHES} synthetic batches with {NOISE_LEVEL*100:.0f}% noise...")
    batches = generate_batches(p_true, N_BATCHES, N_MEAS, NOISE_LEVEL, RNG_SEED)

    # ── Step 2: Build Pyomo model ───────────────────────────────────────────
    print(f"Building Pyomo model ({NFE} FE, {NCP} CP, {N_BATCHES} batches)...")
    model = build_full_model(batches, nfe=NFE, ncp=NCP)

    # Report size (this is the part where you go "wow, that's a lot of variables")
    n_vars = len(list(model.component_data_objects(pyo.Var, active=True)))
    n_cons = len(list(model.component_data_objects(pyo.Constraint, active=True)))
    print(f"Model size: {n_vars:,} variables, {n_cons:,} constraints")

    # ── Step 3: Export NL file ──────────────────────────────────────────────
    # This is why we're here, John. Take this .nl file and throw your best
    # solver at it. May the best algorithm win.
    nl_path = os.path.join(OUTDIR, "cho_parmest.nl")
    model.write(nl_path, format="nl")
    print(f"NL file written: {os.path.abspath(nl_path)}")

    # ── Step 4: Solve with IPOPT ────────────────────────────────────────────
    print("Solving with IPOPT (grab a coffee, this might take a minute)...")
    solver = pyo.SolverFactory("ipopt")
    solver.options["max_iter"] = 5000
    solver.options["tol"] = 1e-6
    solver.options["linear_solver"] = "ma57"
    solver.options["print_level"] = 5
    solver.options["output_file"] = os.path.join(OUTDIR, "ipopt.log")
    result = solver.solve(model, tee=True)

    print(f"\nSolver status: {result.solver.status}")
    print(f"Termination:   {result.solver.termination_condition}")
    print(f"Objective:     {pyo.value(model.obj):.6e}\n")

    # ── Step 5: Parameter recovery table ────────────────────────────────────
    # The moment of truth: did we get the parameters back?
    print(f"{'Parameter':<12} {'True':>12} {'Estimated':>12} {'Rel. Error':>12}")
    print("-" * 50)
    for name in THETA_NAMES:
        true_val = getattr(p_true, name)
        est_val = pyo.value(getattr(model, name))
        rel_err = abs(est_val - true_val) / (abs(true_val) + 1e-30)
        print(f"{name:<12} {true_val:>12.4e} {est_val:>12.4e} {rel_err:>11.2%}")

    # ── Step 6: Training fit plot ───────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.ravel()
    colors = plt.cm.tab10(np.linspace(0, 1, N_BATCHES))

    for k, s in enumerate(STATES):
        ax = axes[k]
        for idx in range(N_BATCHES):
            batch = batches[idx]
            c = colors[idx]
            # Noisy measurements (what the optimizer saw)
            ax.scatter(batch["meas"]["time"], batch["meas"][s],
                       color=c, s=12, alpha=0.5, zorder=2)
            # Ground truth (what we wish we had)
            ax.plot(batch["truth"]["time"], batch["truth"][s],
                    color=c, ls="--", lw=0.8, alpha=0.5, zorder=1)
            # Pyomo collocation solution (what the optimizer found)
            b = model.batch[idx]
            t_model = sorted(pyo.value(t) for t in b.t)
            y_model = [pyo.value(getattr(b, s)[t]) for t in t_model]
            ax.plot(t_model, y_model, color=c, lw=1.2, zorder=3)

        ax.set_title(s)
        ax.set_xlabel("Time (h)")
        ax.grid(True, alpha=0.3)

    fig.suptitle("Parameter Estimation: 10 Training Batches\n"
                 "(dashed=truth, dots=noisy data, solid=estimated model)",
                 fontsize=12)
    fig.tight_layout()
    fig_path = os.path.join(OUTDIR, "parmest_fit.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"\nFigure saved: {fig_path}")
    plt.close(fig)

    # ── Step 7: Validation on held-out batches ──────────────────────────────
    # This is the real test. New IVPs the optimizer has NEVER seen.
    # If the model generalizes, the estimated-parameter trajectories should
    # track the truth. If not, well... back to the drawing board.
    # (Or blame the noise. That's what we usually do.)
    print(f"\n{'='*50}")
    print(f"Validation: {N_VALIDATION} held-out batches (new IVPs)")
    print(f"{'='*50}")

    # Build a Params object with the estimated values
    p_est = Params()
    for name in THETA_NAMES:
        setattr(p_est, name, pyo.value(getattr(model, name)))

    # Generate validation batches with a completely different seed
    val_batches = generate_batches(p_true, N_VALIDATION, N_MEAS,
                                   NOISE_LEVEL, seed=RNG_SEED + 999)

    # Forward-simulate each validation batch with both true and estimated params
    tgrid_val = np.linspace(0, T_FINAL, 500)
    fig_v, axes_v = plt.subplots(2, 3, figsize=(14, 8))
    axes_v = axes_v.ravel()
    val_colors = plt.cm.Set2(np.linspace(0, 1, N_VALIDATION))

    for k, s in enumerate(STATES):
        ax = axes_v[k]
        for vi, vbatch in enumerate(val_batches):
            c = val_colors[vi]
            y0 = vbatch["y0"]

            # Truth: forward sim with the actual parameters
            sim_true = simulate(p_true, tgrid_val, y0)
            # Prediction: forward sim with what we estimated
            sim_est = simulate(p_est, tgrid_val, y0)

            # Noisy measurements (for visual context)
            ax.scatter(vbatch["meas"]["time"], vbatch["meas"][s],
                       color=c, s=15, alpha=0.6, zorder=2,
                       label=f"Val {vi+1} data" if k == 0 else None)
            # True trajectory
            ax.plot(sim_true["time"], sim_true[s],
                    color=c, ls="--", lw=1.0, alpha=0.6, zorder=1)
            # Estimated-parameter trajectory
            ax.plot(sim_est["time"], sim_est[s],
                    color=c, lw=1.5, zorder=3)

        ax.set_title(s)
        ax.set_xlabel("Time (h)")
        ax.grid(True, alpha=0.3)

    fig_v.suptitle(f"Validation: {N_VALIDATION} Held-Out Batches\n"
                   "(dashed=truth, dots=noisy data, solid=estimated params)",
                   fontsize=12)
    fig_v.tight_layout()
    fig_v_path = os.path.join(OUTDIR, "validation_fit.png")
    fig_v.savefig(fig_v_path, dpi=200, bbox_inches="tight")
    print(f"Validation figure saved: {fig_v_path}")
    plt.close(fig_v)


if __name__ == "__main__":
    main()
