"""Regression guard for the L1-penalty NMPC stall.

Reproduces the failure observed in the CSTR eNMPC port (see
saudi-aramco/CSTR_pyomo/jax_ripopt). The NLP is a closed-loop eNMPC step
with a Lagrange-Radau collocation discretization, soft state/control
bounds, and an L1 endpoint penalty (``rho * sum(p + n)`` with
``rho = 1e8``) driving the terminal state to a precomputed economic
steady state. cyipopt/IPOPT solves step 1 in ~21 iterations to
``fun ≈ -118.8``; ripopt currently either reports ``NumericalError`` or
stalls and returns ``Optimal`` with a wildly-suboptimal iterate via the
"acceptable convergence" fallback.

This test is marked ``xfail(strict=True)`` so it will flip to a visible
failure the moment the underlying ripopt filter/line-search fix lands —
serving as both a regression guard and a wake-up signal.

The pathology (captured in the diagnostic session that motivated the
fix):

  1. The NLP scaling floor (was 1e-2, now 1e-8 — Layer 1 fix) used to
     clamp ``obj_scaling`` and leave the scaled L1 gradient at ~1e6,
     freezing ``inf_du`` at 2.4e-1 from iter 0.
  2. With the floor lowered, ripopt reaches ``inf_pr ≈ 5e-10`` by iter
     10 but then takes ``alpha_pr = 0`` every subsequent iteration
     (filter rejects all 12 backtrack trials), ``mu`` stalls at ~1e-4,
     restoration never fires (``theta ≈ 0`` → nothing to restore), and
     eventually ripopt hits max_iter and returns the best-so-far iterate
     which is 5000× worse than IPOPT's optimum.

The IPOPT reference solution satisfies ``res.fun < -100`` on step 1
(IPOPT: -118.80). We assert a generous bound so the test flips the
moment ripopt gets within ~20% of IPOPT's objective.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from ripopt import Problem, minimize

jax.config.update("jax_enable_x64", True)


# ---- CSTR constants (copied from CSTR_pyomo/jax_ripopt/cstr_jax.py) -----
V = 10.0
K = 1.2
CB_LB, CB_UB = 0.45, 1.0
Q_LB, Q_UB = 10.0, 20.0


def _radau_ncp3():
    tau_col = np.array([0.155051025721682, 0.644948974278318, 1.0])
    tau = np.concatenate([[0.0], tau_col])
    n = len(tau)
    D = np.zeros((n - 1, n))
    for j in range(1, n):
        for ell in range(n):
            d = 0.0
            for mm in range(n):
                if mm == ell:
                    continue
                term = 1.0 / (tau[ell] - tau[mm])
                for k in range(n):
                    if k == ell or k == mm:
                        continue
                    term *= (tau[j] - tau[k]) / (tau[ell] - tau[k])
                d += term
            D[j - 1, ell] = d
    return jnp.asarray(D)


def _solve_steady_state():
    def f(z):
        return -z[2] * (2.0 * z[1] - 0.5)

    def g(z):
        ca, cb, q = z[0], z[1], z[2]
        return jnp.array([
            q / V * (1.0 - ca) - K * ca,
            q / V * (-cb) + K * ca,
        ])

    res = minimize(
        f,
        x0=np.array([0.5, 0.5, 12.0]),
        bounds=(
            np.array([0.0, CB_LB, Q_LB]),
            np.array([1e20, CB_UB, Q_UB]),
        ),
        constraints={"fun": g, "lb": [0.0, 0.0], "ub": [0.0, 0.0]},
        options={
            "tol": 1e-8,
            "print_level": 0,
            "hessian_approximation": "limited-memory",
        },
    )
    assert res.success, f"steady state sub-solve failed: {res.status}"
    return float(res.x[0]), float(res.x[1]), float(res.x[2])


def _build_cstr_step1_problem(rho: float = 1.0e8):
    N, ncp, h, caf = 10, 3, 1.0, 1.0
    D = _radau_ncp3()
    ss_ca, ss_cb, ss_q = _solve_steady_state()
    delta = 0.5

    n_ca = N * ncp
    n_cb = N * ncp
    n_q = N
    n_slack_soft = 4 * N
    n_end = 6
    off_ca = 0
    off_cb = off_ca + n_ca
    off_q = off_cb + n_cb
    off_rcb_lb = off_q + n_q
    off_rcb_ub = off_rcb_lb + N
    off_rq_lb = off_rcb_ub + N
    off_rq_ub = off_rq_lb + N
    off_end = off_rq_ub + N
    n = n_ca + n_cb + n_q + n_slack_soft + n_end

    m_coll = 2 * N * ncp
    m_soft = 4 * N
    m_end = 3
    m_stab = 1
    m = m_coll + m_soft + m_end + m_stab

    def unpack(z):
        ca = jnp.reshape(z[off_ca: off_ca + N * ncp], (N, ncp))
        cb = jnp.reshape(z[off_cb: off_cb + N * ncp], (N, ncp))
        q = z[off_q: off_q + N]
        rcb_lb = z[off_rcb_lb: off_rcb_lb + N]
        rcb_ub = z[off_rcb_ub: off_rcb_ub + N]
        rq_lb = z[off_rq_lb: off_rq_lb + N]
        rq_ub = z[off_rq_ub: off_rq_ub + N]
        p1, p2, p3 = z[off_end], z[off_end + 1], z[off_end + 2]
        n1, n2, n3 = z[off_end + 3], z[off_end + 4], z[off_end + 5]
        return ca, cb, q, rcb_lb, rcb_ub, rq_lb, rq_ub, p1, p2, p3, n1, n2, n3

    def boundaries(ca, cb, ca_ic, cb_ic):
        ca_b = jnp.concatenate([jnp.atleast_1d(ca_ic), ca[:, -1]])
        cb_b = jnp.concatenate([jnp.atleast_1d(cb_ic), cb[:, -1]])
        return ca_b, cb_b

    def coll_residuals(ca, cb, q, ca_ic, cb_ic):
        ca_starts = jnp.concatenate([jnp.atleast_1d(ca_ic), ca[:-1, -1]])
        cb_starts = jnp.concatenate([jnp.atleast_1d(cb_ic), cb[:-1, -1]])
        ca_full = jnp.concatenate([ca_starts[:, None], ca], axis=1)
        cb_full = jnp.concatenate([cb_starts[:, None], cb], axis=1)
        dca_dt = (ca_full @ D.T) / h
        dcb_dt = (cb_full @ D.T) / h
        qb = q[:, None]
        rhs_ca = qb / V * (caf - ca) - K * ca
        rhs_cb = qb / V * (-cb) + K * ca
        return dca_dt - rhs_ca, dcb_dt - rhs_cb

    def val_function(ca_b, cb_b):
        return jnp.sum(
            (ca_b[:-1] - ca_b[1:]) ** 2 + (cb_b[:-1] - cb_b[1:]) ** 2
        )

    def fun(z, params):
        ca_ic, cb_ic, _pre_val, _pre_first = params[0], params[1], params[2], params[3]
        (ca, cb, q,
         rcb_lb, rcb_ub, rq_lb, rq_ub,
         p1, p2, p3, n1_, n2_, n3_) = unpack(z)
        ca_b, _ = boundaries(ca, cb, ca_ic, cb_ic)
        _, cb_b = boundaries(ca, cb, ca_ic, cb_ic)
        eco = jnp.sum(-q * (2.0 * cb_b[:-1] - 0.5))
        pen = jnp.sum(rcb_lb + rcb_ub + rq_lb + rq_ub)
        l1 = (p1 + n1_) + (p2 + n2_) + (p3 + n3_)
        return eco + rho * pen + rho * l1

    def g_all(z, params):
        ca_ic, cb_ic, pre_val, pre_first = params[0], params[1], params[2], params[3]
        (ca, cb, q,
         rcb_lb, rcb_ub, rq_lb, rq_ub,
         p1, p2, p3, n1_, n2_, n3_) = unpack(z)
        res_ca, res_cb = coll_residuals(ca, cb, q, ca_ic, cb_ic)
        ca_b, cb_b = boundaries(ca, cb, ca_ic, cb_ic)
        cb_end_fe = cb_b[1:]
        g_cb_lb = cb_end_fe + rcb_lb - CB_LB
        g_cb_ub = CB_UB + rcb_ub - cb_end_fe
        g_q_lb = q + rq_lb - Q_LB
        g_q_ub = Q_UB + rq_ub - q
        g_p1 = ca_b[-1] - ss_ca - p1 + n1_
        g_p2 = cb_b[-1] - ss_cb - p2 + n2_
        g_p3 = q[-1] - ss_q - p3 + n3_
        vfun = val_function(ca_b, cb_b)
        g_stab = vfun - pre_val + delta * pre_first
        return jnp.concatenate([
            res_ca.flatten(), res_cb.flatten(),
            g_cb_lb, g_cb_ub, g_q_lb, g_q_ub,
            jnp.array([g_p1, g_p2, g_p3]),
            jnp.array([g_stab]),
        ])

    cl = np.zeros(m)
    cu = np.zeros(m)
    s = m_coll
    cu[s: s + m_soft] = 1e20
    s += m_soft
    s += m_end
    cl[s: s + m_stab] = -1e20

    lb = np.zeros(n)
    ub = np.full(n, 1e20)

    z0 = np.zeros(n)
    z0[off_ca: off_ca + n_ca] = 0.5
    z0[off_cb: off_cb + n_cb] = 0.5
    z0[off_q: off_q + n_q] = 12.0

    params0 = jnp.array([0.1, 1.0, 1.0e6, 0.0])

    return Problem(
        fun,
        x0=z0,
        bounds=(lb, ub),
        constraints={"fun": g_all, "lb": cl, "ub": cu},
        options={"tol": 1e-6, "max_iter": 3000, "print_level": 0},
        params=params0,
        jac_mode="reverse",
        sparsity="dense",
    )


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Layer-2 IPM stall: filter/line-search gets alpha_pr=0 every "
        "iteration after primal becomes feasible at iter ~10. See "
        "test_l1_penalty_stall.py module docstring for the full trace."
    ),
)
def test_cstr_l1_step1_matches_ipopt():
    """NMPC step 1 must reach an IPOPT-quality objective.

    IPOPT (via cyipopt) converges to ``fun ≈ -118.80`` in 21 iterations.
    ripopt currently returns ``fun ≈ 5.7e5`` (via "acceptable
    convergence" fallback) or ``NumericalError``. We accept any ripopt
    solution with ``res.fun < -100`` as "close enough to IPOPT".
    """
    prob = _build_cstr_step1_problem(rho=1.0e8)
    res = prob.solve()
    assert res.success, f"solve failed: status={res.status}"
    assert res.fun < -100.0, (
        f"ripopt returned a suboptimal iterate: fun={res.fun:.3e} "
        f"(IPOPT reference: -118.80; anything ≥ -100 means the L1 slacks "
        f"were not driven to the bound)"
    )
