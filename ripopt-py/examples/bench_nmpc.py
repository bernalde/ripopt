"""NMPC-style reuse benchmark for ripopt-py.

Compares the wall time of solving the *same* NLP structure 100 times in a
row through:

  1. ``ripopt.minimize`` — single-shot scipy-style API, re-traces and
     re-jits the Lagrangian Hessian on every call.
  2. ``ripopt.Problem`` — persistent NLP object that JIT-compiles every
     callback exactly once and reuses the XLA cache across all solves.

On a real closed-loop eNMPC workload (see saudi-aramco/CSTR_pyomo/jax_ripopt),
switching from (1) to (2) brings wall time within ~1.5× of cyipopt's
reference implementation. This script reproduces the effect on a smaller,
self-contained LQR-style tracking problem so the perf gap is visible in
this repo without any external CUTEst/CSTR setup.

Run::

    python examples/bench_nmpc.py

Expected output: Problem.solve path is several× faster than re-entering
minimize, and the gap is dominated by JAX tracing + XLA lowering of
``jax.hessian(lagrangian)`` on every minimize call.
"""

from __future__ import annotations

import time

import jax
import jax.numpy as jnp
import numpy as np

from ripopt import Problem, minimize

jax.config.update("jax_enable_x64", True)


# ------------------- problem definition -------------------

# LQR-style tracking: scalar state, scalar control, simple LTI dynamics.
# Decision vars laid out as [x_1, ..., x_N, u_0, ..., u_{N-1}], with the
# current state x_0 threaded in as a `params` PyTree.

N = 30                 # horizon length → n = 2N = 60 decision variables
A = 0.9
B = 0.5
R = 1.0                # reference state
Q = 1.0                # state cost weight
W = 0.01               # control cost weight
U_LB = -2.0
U_UB = 2.0


def _split(z):
    x = z[:N]
    u = z[N:]
    return x, u


def stage_cost(z, x0):
    """Sum of (x_k - R)^2 and W*u_k^2 over the horizon."""
    del x0
    x, u = _split(z)
    return jnp.sum(Q * (x - R) ** 2) + jnp.sum(W * u ** 2)


def dynamics(z, x0):
    """Equality constraints g_k = x_{k+1} - A*x_k - B*u_k, 0 <= k < N."""
    x, u = _split(z)
    x_prev = jnp.concatenate([jnp.atleast_1d(x0[0]), x[:-1]])
    return x - A * x_prev - B * u


def make_bounds():
    lb = np.concatenate([np.full(N, -np.inf), np.full(N, U_LB)])
    ub = np.concatenate([np.full(N, np.inf), np.full(N, U_UB)])
    return lb, ub


def initial_guess():
    return np.zeros(2 * N)


# ------------------- benchmarks -------------------

N_STEPS = 100
RNG = np.random.default_rng(42)
# Closed-loop NMPC simulation: at each step we solve for optimal u_0, apply
# it to the plant, integrate one step forward, and repeat. This gives a
# smoothly-evolving initial state so consecutive solves share structure —
# the exact regime in which dual warm-starting is meant to help.
_X0_INITIAL = 4.0
_DISTURBANCE = 0.05 * RNG.standard_normal(N_STEPS)


def run_minimize_path() -> tuple[float, float, int]:
    """Baseline: re-enter minimize() every step with x0 pinned via bounds.

    Because minimize() has no `params` support via closure capture without
    retracing, we embed x0 as an extra decision variable pinned by equal
    bounds (the same workaround used in the real CSTR port).
    """

    def f(z):
        # z = [x_0_pinned, x_1..x_N, u_0..u_{N-1}]  -> length 2N+1
        x0 = z[0:1]
        zz = z[1:]
        return stage_cost(zz, x0)

    def g(z):
        x0 = z[0:1]
        zz = z[1:]
        return dynamics(zz, x0)

    lb_z, ub_z = make_bounds()
    plant_x = _X0_INITIAL
    z0 = np.concatenate([[plant_x], initial_guess()])
    x_opt = z0.copy()

    t_total_start = time.perf_counter()
    total_iters = 0
    first_solve = 0.0

    for step in range(N_STEPS):
        lb = np.concatenate([[plant_x], lb_z])
        ub = np.concatenate([[plant_x], ub_z])
        x_opt[0] = plant_x
        t0 = time.perf_counter()
        res = minimize(
            f,
            x0=x_opt,
            bounds=(lb, ub),
            constraints={"fun": g, "lb": 0.0, "ub": 0.0},
            options={"print_level": 0, "tol": 1e-6, "max_iter": 200},
        )
        dt = time.perf_counter() - t0
        if step == 0:
            first_solve = dt
        assert res.success, f"minimize path failed at step {step}: {res.status}"
        x_opt = res.x
        total_iters += res.iterations
        # Apply u_0 to the plant and evolve one step (plus disturbance).
        u0 = float(res.x[N + 1])  # after [x0_pinned, x_1..x_N, u_0..u_{N-1}]
        plant_x = A * plant_x + B * u0 + _DISTURBANCE[step]

    return time.perf_counter() - t_total_start, first_solve, total_iters


def run_problem_path(*, warm_start: bool) -> tuple[float, float, int]:
    """Persistent Problem: construct once, update_parameters each step.

    When ``warm_start`` is True, each solve past the first reuses the
    previous solve's dual iterates via the Ipopt-style warm_start hook.
    """

    lb_z, ub_z = make_bounds()
    plant_x = _X0_INITIAL
    prob = Problem(
        stage_cost,
        x0=initial_guess(),
        bounds=(lb_z, ub_z),
        constraints={"fun": dynamics, "lb": 0.0, "ub": 0.0},
        options={"print_level": 0, "tol": 1e-6, "max_iter": 200},
        params=jnp.array([plant_x]),
        jac_mode="reverse",
    )

    t_total_start = time.perf_counter()
    total_iters = 0
    first_solve = 0.0

    for step in range(N_STEPS):
        prob.update_parameters(jnp.array([plant_x]))
        t0 = time.perf_counter()
        res = prob.solve(warm_start=warm_start)
        dt = time.perf_counter() - t0
        if step == 0:
            first_solve = dt
        assert res.success, f"Problem path failed at step {step}: {res.status}"
        total_iters += res.iterations
        u0 = float(res.x[N])  # [x_1..x_N, u_0..u_{N-1}] → u_0 is index N
        plant_x = A * plant_x + B * u0 + _DISTURBANCE[step]

    return time.perf_counter() - t_total_start, first_solve, total_iters


def main() -> None:
    # Warm-up JAX and Python import caches (no wall-time effect on either
    # path since both hit the same JAX global cache, but makes the first
    # measured step representative of steady-state).
    _ = stage_cost(initial_guess(), jnp.array([0.0]))
    _ = dynamics(initial_guess(), jnp.array([0.0]))

    print(f"NMPC-style benchmark: N={N}, n={2*N}, m={N}, steps={N_STEPS}")
    print()

    t_min, first_min, iters_min = run_minimize_path()
    print(
        f"  minimize() path         : total {t_min:7.3f}s   "
        f"(first solve {first_min*1000:7.1f}ms,  "
        f"per-step avg {t_min/N_STEPS*1000:6.1f}ms,  "
        f"total iters {iters_min})"
    )

    t_prob_cold, first_cold, iters_cold = run_problem_path(warm_start=False)
    print(
        f"  Problem cold            : total {t_prob_cold:7.3f}s   "
        f"(first solve {first_cold*1000:7.1f}ms,  "
        f"per-step avg {t_prob_cold/N_STEPS*1000:6.1f}ms,  "
        f"total iters {iters_cold})"
    )

    t_prob_warm, first_warm, iters_warm = run_problem_path(warm_start=True)
    print(
        f"  Problem warm_start=True : total {t_prob_warm:7.3f}s   "
        f"(first solve {first_warm*1000:7.1f}ms,  "
        f"per-step avg {t_prob_warm/N_STEPS*1000:6.1f}ms,  "
        f"total iters {iters_warm})"
    )

    print()
    print(
        f"  speedup vs minimize()   : "
        f"{t_min / t_prob_cold:.2f}× (Problem cold) / "
        f"{t_min / t_prob_warm:.2f}× (Problem warm)"
    )


if __name__ == "__main__":
    main()
