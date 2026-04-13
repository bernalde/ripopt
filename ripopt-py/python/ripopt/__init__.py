"""Direct Python interface to ripopt, with JAX autodiff for derivatives.

Users provide only ``fun`` (and optionally constraint functions); gradients,
Jacobians, and the Lagrangian Hessian are built automatically via ``jax.grad``,
``jax.jacfwd``, and ``jax.hessian``.

Example
-------
>>> import jax.numpy as jnp
>>> from ripopt import minimize
>>> def f(x): return (x[0] - 1.0)**2 + (x[1] - 2.5)**2
>>> res = minimize(f, x0=[0.0, 0.0])
>>> float(res.x[0]), float(res.x[1])
(1.0, 2.5)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional, Sequence, Union

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

from ._ripopt import _solve  # type: ignore[attr-defined]

__all__ = ["minimize", "OptimizeResult"]


@dataclass
class OptimizeResult:
    x: np.ndarray
    fun: float
    status: str
    success: bool
    iterations: int
    constraint_multipliers: np.ndarray
    bound_multipliers_lower: np.ndarray
    bound_multipliers_upper: np.ndarray
    constraint_values: np.ndarray
    wall_time_secs: float
    final_primal_inf: float
    final_dual_inf: float
    final_compl: float


Bounds = Union[None, tuple, Sequence]
Constraint = dict
Constraints = Union[None, Constraint, Sequence[Constraint]]


def _normalize_bounds(bounds: Bounds, n: int) -> tuple[np.ndarray, np.ndarray]:
    if bounds is None:
        return (
            np.full(n, -np.inf, dtype=np.float64),
            np.full(n, np.inf, dtype=np.float64),
        )

    if isinstance(bounds, tuple) and len(bounds) == 2 and not isinstance(bounds[0], (list, tuple)):
        lb, ub = bounds
        lb_arr = np.broadcast_to(np.asarray(lb, dtype=np.float64), (n,)).astype(np.float64)
        ub_arr = np.broadcast_to(np.asarray(ub, dtype=np.float64), (n,)).astype(np.float64)
        return np.ascontiguousarray(lb_arr), np.ascontiguousarray(ub_arr)

    # scipy-style: sequence of (lb_i, ub_i) pairs with None for missing sides
    pairs = list(bounds)
    if len(pairs) != n:
        raise ValueError(f"bounds has {len(pairs)} entries, expected {n}")
    lb = np.array(
        [(-np.inf if p[0] is None else p[0]) for p in pairs], dtype=np.float64
    )
    ub = np.array(
        [(np.inf if p[1] is None else p[1]) for p in pairs], dtype=np.float64
    )
    return lb, ub


def _normalize_constraints(
    constraints: Constraints,
    x0_jnp,
    n: int,
) -> tuple[Optional[Callable], np.ndarray, np.ndarray, int]:
    if constraints is None:
        return None, np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.float64), 0

    items = [constraints] if isinstance(constraints, dict) else list(constraints)

    wrapped = []
    lbs, ubs, sizes = [], [], []
    for c in items:
        if not isinstance(c, dict) or "fun" not in c:
            raise ValueError("each constraint must be a dict containing 'fun'")
        raw = c["fun"]

        def make_wrapper(f):
            return lambda x: jnp.atleast_1d(f(x))

        wf = make_wrapper(raw)
        out = np.asarray(wf(x0_jnp))
        k = int(out.size)
        sizes.append(k)
        wrapped.append(wf)

        lb = c.get("lb", 0.0)
        ub = c.get("ub", 0.0)
        lbs.append(
            np.broadcast_to(np.asarray(lb, dtype=np.float64), (k,)).astype(np.float64)
        )
        ubs.append(
            np.broadcast_to(np.asarray(ub, dtype=np.float64), (k,)).astype(np.float64)
        )

    m = int(sum(sizes))
    if len(wrapped) == 1:
        g_base = wrapped[0]
    else:
        _funs = tuple(wrapped)

        def g_base(x):
            return jnp.concatenate([f(x) for f in _funs])

    g_l = np.concatenate(lbs).astype(np.float64)
    g_u = np.concatenate(ubs).astype(np.float64)
    return g_base, g_l, g_u, m


def _to_scalar(fn):
    def w(x):
        return float(fn(x))
    return w


def _to_1d(fn):
    def w(x):
        return np.ascontiguousarray(np.asarray(fn(x), dtype=np.float64))
    return w


def _to_2d_hess(fn):
    def w(x, sigma, lam):
        return np.ascontiguousarray(
            np.asarray(fn(x, sigma, lam), dtype=np.float64)
        )
    return w


def _build_sparsity(
    mode: str,
    jac_g_np: Optional[Callable],
    hess_l_np: Optional[Callable],
    n: int,
    m: int,
    x0: np.ndarray,
):
    """Return ``(jac_rows, jac_cols, hes_rows, hes_cols)`` triplet lists.

    ``'dense'`` declares every entry structurally nonzero.

    ``'detect'`` evaluates the Jacobian at two points and the Lagrangian
    Hessian with two random ``(sigma, lambda)`` probes, then takes the
    union of nonzero entries. With random multipliers, cancellations in
    ``sigma*H_f + sum_i lambda_i * H_{g_i}`` occur with probability zero,
    so any structurally nonzero entry shows up.
    """
    if mode == "dense":
        if m > 0:
            jac_rows = np.repeat(np.arange(m, dtype=np.int64), n).tolist()
            jac_cols = np.tile(np.arange(n, dtype=np.int64), m).tolist()
        else:
            jac_rows, jac_cols = [], []
        if hess_l_np is not None:
            hes_rows, hes_cols = [], []
            for i in range(n):
                for j in range(i + 1):
                    hes_rows.append(i)
                    hes_cols.append(j)
        else:
            hes_rows, hes_cols = [], []
        return jac_rows, jac_cols, hes_rows, hes_cols

    if mode != "detect":
        raise ValueError(f"sparsity must be 'dense' or 'detect', got {mode!r}")

    rng = np.random.default_rng(0)
    perturb = 0.01 * rng.standard_normal(n)

    if m > 0 and jac_g_np is not None:
        J1 = np.asarray(jac_g_np(x0))
        J2 = np.asarray(jac_g_np(x0 + perturb))
        mask_j = (J1 != 0.0) | (J2 != 0.0)
        rows_arr, cols_arr = np.nonzero(mask_j)
        jac_rows = rows_arr.astype(np.int64).tolist()
        jac_cols = cols_arr.astype(np.int64).tolist()
    else:
        jac_rows, jac_cols = [], []

    if hess_l_np is not None:
        lam_a = rng.standard_normal(m) if m > 0 else np.zeros(0)
        lam_b = rng.standard_normal(m) if m > 0 else np.zeros(0)
        H1 = np.asarray(hess_l_np(x0, 1.0, lam_a))
        H2 = np.asarray(hess_l_np(x0 + perturb, 1.0, lam_b))
        mask_h = (H1 != 0.0) | (H2 != 0.0)
        # Symmetrize — if either (i,j) or (j,i) is nonzero, mark it.
        mask_h = mask_h | mask_h.T
        hes_rows, hes_cols = [], []
        for i in range(n):
            for j in range(i + 1):
                if mask_h[i, j]:
                    hes_rows.append(i)
                    hes_cols.append(j)
    else:
        hes_rows, hes_cols = [], []

    return jac_rows, jac_cols, hes_rows, hes_cols


def minimize(
    fun: Callable,
    x0,
    *,
    bounds: Bounds = None,
    constraints: Constraints = None,
    options: Optional[dict] = None,
    sparsity: str = "dense",
) -> OptimizeResult:
    """Solve ``min fun(x)`` subject to bounds and (optional) constraints.

    Parameters
    ----------
    fun
        Scalar-valued objective, JAX-traceable: ``fun(x) -> scalar``.
    x0
        Initial point, array-like of length ``n``.
    bounds
        ``None``, ``(lb, ub)`` where each is a scalar or length-``n`` array,
        or a length-``n`` sequence of ``(lb_i, ub_i)`` pairs (scipy-style;
        ``None`` in a pair means no bound).
    constraints
        A dict ``{'fun': g, 'lb': ..., 'ub': ...}`` or a list of such dicts.
        Each ``g`` must be JAX-traceable and return a scalar or 1-D array.
        Use ``lb == ub`` for equality constraints.
    options
        Optional solver options dict. Recognized keys: ``tol``, ``max_iter``,
        ``print_level``, ``constr_viol_tol``, ``dual_inf_tol``,
        ``compl_inf_tol``, ``max_wall_time``, ``mu_init``,
        ``hessian_approximation`` (``'exact'`` or ``'limited-memory'``).
    sparsity
        Sparsity strategy for Jacobian and Hessian structure:

        * ``'dense'`` (default) — every entry is declared a structural nonzero.
          Safe and simple; best for small problems.
        * ``'detect'`` — probe the Jacobian at two points and the Lagrangian
          Hessian with two random ``(sigma, lambda)`` samples, take the union
          of nonzero entries. Does not reduce JAX work, but shrinks the KKT
          system ripopt factors. For larger problems with genuine sparsity,
          this can substantially reduce iteration cost.
    """
    x0_np = np.ascontiguousarray(np.asarray(x0, dtype=np.float64))
    if x0_np.ndim != 1:
        raise ValueError("x0 must be 1-D")
    n = int(x0_np.size)
    x0_jnp = jnp.asarray(x0_np)

    x_l, x_u = _normalize_bounds(bounds, n)
    g_base, g_l, g_u, m = _normalize_constraints(constraints, x0_jnp, n)

    f_jit = jax.jit(fun)
    grad_f_jit = jax.jit(jax.grad(fun))

    if g_base is not None:
        g_jit = jax.jit(g_base)
        jac_g_jit = jax.jit(jax.jacfwd(g_base))
    else:
        g_jit = None
        jac_g_jit = None

    opts = dict(options or {})
    use_lbfgs = opts.get("hessian_approximation") == "limited-memory"

    if use_lbfgs:
        hess_l_jit = None
    else:
        if g_base is None:
            def lagrangian(x, sigma, lam):
                return sigma * fun(x)
        else:
            def lagrangian(x, sigma, lam):
                return sigma * fun(x) + jnp.dot(lam, g_base(x))
        hess_l_jit = jax.jit(jax.hessian(lagrangian))

    f_np = _to_scalar(f_jit)
    grad_f_np = _to_1d(grad_f_jit)
    g_np = _to_1d(g_jit) if g_jit is not None else None
    jac_g_np = _to_1d(jac_g_jit) if jac_g_jit is not None else None
    hess_l_np = _to_2d_hess(hess_l_jit) if hess_l_jit is not None else None

    # Pre-trace each callable so the first solver iteration isn't inflated
    # by JAX tracing/compilation time.
    _ = f_np(x0_np)
    _ = grad_f_np(x0_np)
    if g_np is not None:
        _ = g_np(x0_np)
        _ = jac_g_np(x0_np)
    if hess_l_np is not None:
        _ = hess_l_np(x0_np, 1.0, np.zeros(m, dtype=np.float64))

    jac_rows, jac_cols, hes_rows, hes_cols = _build_sparsity(
        sparsity, jac_g_np, hess_l_np, n, m, x0_np
    )

    raw = _solve(
        n, m,
        x0_np, x_l, x_u, g_l, g_u,
        f_np, grad_f_np, g_np, jac_g_np, hess_l_np,
        jac_rows, jac_cols, hes_rows, hes_cols,
        opts,
    )

    return OptimizeResult(
        x=np.asarray(raw["x"], dtype=np.float64),
        fun=float(raw["fun"]),
        status=str(raw["status"]),
        success=bool(raw["success"]),
        iterations=int(raw["iterations"]),
        constraint_multipliers=np.asarray(
            raw["constraint_multipliers"], dtype=np.float64
        ),
        bound_multipliers_lower=np.asarray(
            raw["bound_multipliers_lower"], dtype=np.float64
        ),
        bound_multipliers_upper=np.asarray(
            raw["bound_multipliers_upper"], dtype=np.float64
        ),
        constraint_values=np.asarray(raw["constraint_values"], dtype=np.float64),
        wall_time_secs=float(raw["wall_time_secs"]),
        final_primal_inf=float(raw["final_primal_inf"]),
        final_dual_inf=float(raw["final_dual_inf"]),
        final_compl=float(raw["final_compl"]),
    )
