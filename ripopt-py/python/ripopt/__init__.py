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

__all__ = ["minimize", "Problem", "OptimizeResult"]


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
    params: Any = None,
) -> tuple[Optional[Callable], np.ndarray, np.ndarray, int]:
    if constraints is None:
        return None, np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.float64), 0

    items = [constraints] if isinstance(constraints, dict) else list(constraints)
    has_params = params is not None

    wrapped = []
    lbs, ubs, sizes = [], [], []
    for c in items:
        if not isinstance(c, dict) or "fun" not in c:
            raise ValueError("each constraint must be a dict containing 'fun'")
        raw = c["fun"]

        if has_params:
            def make_wrapper(f):
                return lambda x, p: jnp.atleast_1d(f(x, p))
            wf = make_wrapper(raw)
            out = np.asarray(wf(x0_jnp, params))
        else:
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
        if has_params:
            def g_base(x, p):
                return jnp.concatenate([f(x, p) for f in _funs])
        else:
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


def _raw_to_result(raw) -> OptimizeResult:
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


class Problem:
    """Persistent NLP object — builds JAX-autodiff callbacks once, reuses them
    across many ``solve()`` calls.

    This is the right shape for workloads (e.g. closed-loop NMPC) that solve
    the *same* NLP structure repeatedly with only a handful of numeric
    parameters changing between calls. The ``minimize`` convenience function
    constructs a fresh ``Problem`` on every call and is therefore the wrong
    choice for such workloads: every call would retrace and re-lower the
    Lagrangian Hessian.

    A ``Problem`` instead:

    * JIT-compiles ``fun``, ``grad(fun)``, ``g``, ``jac(g)`` and
      ``hessian(lagrangian)`` exactly **once** in ``__init__``. The
      Lagrangian is defined outside any per-solve scope so its XLA lowering
      is cached for the lifetime of the object.
    * Builds the sparsity pattern and bounds arrays once.
    * Warm-starts each solve from the previous solution by default.
    * Optionally threads a user-supplied ``params`` PyTree through ``fun``
      and each constraint as a second positional argument; mutate it with
      :meth:`update_parameters` between solves so the compiled HLO is
      reused without retracing.

    Parameters
    ----------
    fun
        Scalar objective. JAX-traceable. Signature is ``fun(x)`` when
        ``params is None``; otherwise ``fun(x, params)``.
    x0
        Initial point, 1-D array-like of length ``n``.
    bounds, constraints, options, sparsity
        Same semantics as :func:`minimize`. When ``params`` is supplied,
        every constraint's ``fun`` must take ``(x, params)`` as well.
    params
        Optional PyTree of extra inputs (typically jnp arrays or Python
        scalars). Threaded as a second positional arg into ``fun``, each
        constraint, and the Lagrangian Hessian. Mutate via
        :meth:`update_parameters`; as long as the PyTree structure and each
        leaf's (shape, dtype) are stable across calls, JAX reuses the
        compiled HLO.
    jac_mode
        ``'forward'`` (default) uses :func:`jax.jacfwd`; ``'reverse'`` uses
        :func:`jax.jacrev`. Reverse mode is typically cheaper when
        ``m < n`` (wide problems with fewer constraints than variables).
    """

    def __init__(
        self,
        fun: Callable,
        x0,
        *,
        bounds: Bounds = None,
        constraints: Constraints = None,
        options: Optional[dict] = None,
        sparsity: str = "dense",
        params: Any = None,
        jac_mode: str = "forward",
    ) -> None:
        x0_np = np.ascontiguousarray(np.asarray(x0, dtype=np.float64))
        if x0_np.ndim != 1:
            raise ValueError("x0 must be 1-D")
        n = int(x0_np.size)
        x0_jnp = jnp.asarray(x0_np)

        if jac_mode not in ("forward", "reverse"):
            raise ValueError(
                f"jac_mode must be 'forward' or 'reverse', got {jac_mode!r}"
            )

        self._n = n
        self._params = params
        pass_params = params is not None
        self._pass_params = pass_params

        x_l, x_u = _normalize_bounds(bounds, n)
        g_base, g_l, g_u, m = _normalize_constraints(
            constraints, x0_jnp, n, params=params
        )

        # Build JIT'd callables ONCE. With pass_params=True, params is a
        # traced argument so mutating self._params does not invalidate the
        # XLA cache as long as (pytree structure, leaf shape/dtype) is stable.
        if pass_params:
            f_jit = jax.jit(fun)
            grad_f_jit = jax.jit(jax.grad(fun, argnums=0))
        else:
            f_jit = jax.jit(fun)
            grad_f_jit = jax.jit(jax.grad(fun))

        jac_fn = jax.jacfwd if jac_mode == "forward" else jax.jacrev
        if g_base is not None:
            g_jit = jax.jit(g_base)
            jac_g_jit = jax.jit(jac_fn(g_base, argnums=0))
        else:
            g_jit = None
            jac_g_jit = None

        opts = dict(options or {})
        use_lbfgs = opts.get("hessian_approximation") == "limited-memory"

        # Lagrangian is defined ONCE — NOT inside any per-solve closure — so
        # jax.hessian + jax.jit trace and lower it exactly one time for the
        # lifetime of this Problem. This is the single biggest perf win on
        # NMPC workloads.
        if use_lbfgs:
            hess_l_jit = None
        else:
            if g_base is None:
                if pass_params:
                    def lagrangian(x, sigma, lam, p):
                        return sigma * fun(x, p)
                else:
                    def lagrangian(x, sigma, lam):
                        return sigma * fun(x)
            else:
                if pass_params:
                    def lagrangian(x, sigma, lam, p):
                        return sigma * fun(x, p) + jnp.dot(lam, g_base(x, p))
                else:
                    def lagrangian(x, sigma, lam):
                        return sigma * fun(x) + jnp.dot(lam, g_base(x))
            hess_l_jit = jax.jit(jax.hessian(lagrangian, argnums=0))

        # numpy-facing wrappers. When params is threaded, these re-read
        # self._params on every call so update_parameters takes effect
        # without rebuilding any jit'd function.
        if pass_params:
            problem_self = self  # for closure readability

            def f_np(x):
                return float(f_jit(x, problem_self._params))

            def grad_f_np(x):
                return np.ascontiguousarray(
                    np.asarray(
                        grad_f_jit(x, problem_self._params), dtype=np.float64
                    )
                )

            if g_jit is not None:
                def g_np(x):
                    return np.ascontiguousarray(
                        np.asarray(g_jit(x, problem_self._params), dtype=np.float64)
                    )

                def jac_g_np(x):
                    return np.ascontiguousarray(
                        np.asarray(
                            jac_g_jit(x, problem_self._params), dtype=np.float64
                        )
                    )
            else:
                g_np = None
                jac_g_np = None

            if hess_l_jit is not None:
                def hess_l_np(x, sigma, lam):
                    return np.ascontiguousarray(
                        np.asarray(
                            hess_l_jit(x, sigma, lam, problem_self._params),
                            dtype=np.float64,
                        )
                    )
            else:
                hess_l_np = None
        else:
            f_np = _to_scalar(f_jit)
            grad_f_np = _to_1d(grad_f_jit)
            g_np = _to_1d(g_jit) if g_jit is not None else None
            jac_g_np = _to_1d(jac_g_jit) if jac_g_jit is not None else None
            hess_l_np = _to_2d_hess(hess_l_jit) if hess_l_jit is not None else None

        # Pre-trace / warm the XLA cache ONCE at construction. Subsequent
        # solves pay zero compile cost for these callables.
        _ = f_np(x0_np)
        _ = grad_f_np(x0_np)
        if g_np is not None:
            _ = g_np(x0_np)
            _ = jac_g_np(x0_np)
        if hess_l_np is not None:
            _ = hess_l_np(x0_np, 1.0, np.zeros(m, dtype=np.float64))

        # Sparsity ONCE. Structural nonzeros don't depend on the numeric
        # value of params for typical well-posed problems.
        jac_rows, jac_cols, hes_rows, hes_cols = _build_sparsity(
            sparsity, jac_g_np, hess_l_np, n, m, x0_np
        )

        self._m = m
        self._x_current = x0_np.copy()
        self._x_l = x_l
        self._x_u = x_u
        self._g_l = g_l
        self._g_u = g_u
        self._last_lam_g: Optional[np.ndarray] = None
        self._last_z_l: Optional[np.ndarray] = None
        self._last_z_u: Optional[np.ndarray] = None
        self._f_np = f_np
        self._grad_f_np = grad_f_np
        self._g_np = g_np
        self._jac_g_np = jac_g_np
        self._hess_l_np = hess_l_np
        self._jac_rows = jac_rows
        self._jac_cols = jac_cols
        self._hes_rows = hes_rows
        self._hes_cols = hes_cols
        self._options = opts

    @property
    def n(self) -> int:
        return self._n

    @property
    def m(self) -> int:
        return self._m

    @property
    def x(self) -> np.ndarray:
        """Most recent primal iterate (initial point before the first solve)."""
        return self._x_current

    def update_parameters(self, params: Any) -> None:
        """Replace the ``params`` PyTree. JAX re-uses the compiled HLO as
        long as the PyTree structure and each leaf's (shape, dtype) match.
        """
        if not self._pass_params and params is not None:
            raise ValueError(
                "Problem was constructed without params=...; "
                "cannot set params on this Problem"
            )
        self._params = params

    def update_bounds(
        self,
        lb: Optional[Any] = None,
        ub: Optional[Any] = None,
    ) -> None:
        """Fast in-place bounds update. Scalar or length-n array accepted."""
        n = self._n
        if lb is not None:
            self._x_l = np.ascontiguousarray(
                np.broadcast_to(np.asarray(lb, dtype=np.float64), (n,)).astype(
                    np.float64
                )
            )
        if ub is not None:
            self._x_u = np.ascontiguousarray(
                np.broadcast_to(np.asarray(ub, dtype=np.float64), (n,)).astype(
                    np.float64
                )
            )

    def solve(
        self,
        x0: Any = None,
        *,
        bounds: Bounds = None,
        options_overrides: Optional[dict] = None,
        warm_start: bool = False,
        lam0: Any = None,
        z_l0: Any = None,
        z_u0: Any = None,
    ) -> OptimizeResult:
        """Re-enter the solver with the cached callbacks.

        Parameters
        ----------
        x0
            Initial primal iterate. Defaults to the previous solution (warm
            start of primal), or the construction-time ``x0`` on the first call.
        bounds
            Optional bounds override; same forms as in ``__init__``. When
            omitted, the current cached bounds are used.
        options_overrides
            Optional dict merged on top of the construction-time options.
        warm_start
            When ``True`` and this Problem has already been solved at least
            once, reuse the previous solve's constraint multipliers and bound
            multipliers to warm-start the dual iterates on this solve. This
            is the ripopt equivalent of Ipopt's ``warm_start_init_point=yes``
            and typically cuts the NMPC iteration count by 3-4× on step 2+.
            Ignored on the very first solve, where the solver falls back to
            its own least-squares multiplier estimate.
        lam0, z_l0, z_u0
            Explicit warm-start dual iterates. Must be supplied as a
            complete triple (constraint multipliers, lower-bound multipliers,
            upper-bound multipliers); partial specs raise ``ValueError``.
            When supplied, ``warm_start=True`` is implied.
        """
        n = self._n
        if x0 is None:
            x0_np = self._x_current
        else:
            x0_np = np.ascontiguousarray(np.asarray(x0, dtype=np.float64))
            if x0_np.ndim != 1 or x0_np.size != n:
                raise ValueError(f"x0 must be 1-D of length {n}")

        if bounds is not None:
            x_l, x_u = _normalize_bounds(bounds, n)
        else:
            x_l, x_u = self._x_l, self._x_u

        # Resolve warm-start multiplier triple. Explicit lam0/z_l0/z_u0
        # always take precedence over stored previous-solve values.
        explicit_count = sum(v is not None for v in (lam0, z_l0, z_u0))
        if 0 < explicit_count < 3:
            raise ValueError(
                "lam0, z_l0, z_u0 must all be provided together or not at all"
            )
        if explicit_count == 3:
            init_lam = np.ascontiguousarray(np.asarray(lam0, dtype=np.float64))
            init_zl = np.ascontiguousarray(np.asarray(z_l0, dtype=np.float64))
            init_zu = np.ascontiguousarray(np.asarray(z_u0, dtype=np.float64))
        elif warm_start and self._last_lam_g is not None:
            init_lam = self._last_lam_g
            init_zl = self._last_z_l
            init_zu = self._last_z_u
        else:
            init_lam = None
            init_zl = None
            init_zu = None

        if init_lam is not None:
            if init_lam.shape != (self._m,):
                raise ValueError(
                    f"lam0 must have shape ({self._m},), got {init_lam.shape}"
                )
            if init_zl.shape != (n,) or init_zu.shape != (n,):
                raise ValueError(
                    f"z_l0, z_u0 must have shape ({n},); got "
                    f"{init_zl.shape}, {init_zu.shape}"
                )

        if options_overrides:
            opts = dict(self._options)
            opts.update(options_overrides)
        else:
            opts = self._options

        raw = _solve(
            n,
            self._m,
            x0_np,
            x_l,
            x_u,
            self._g_l,
            self._g_u,
            self._f_np,
            self._grad_f_np,
            self._g_np,
            self._jac_g_np,
            self._hess_l_np,
            self._jac_rows,
            self._jac_cols,
            self._hes_rows,
            self._hes_cols,
            opts,
            init_lam,
            init_zl,
            init_zu,
        )

        result = _raw_to_result(raw)
        self._x_current = result.x
        self._last_lam_g = result.constraint_multipliers
        self._last_z_l = result.bound_multipliers_lower
        self._last_z_u = result.bound_multipliers_upper
        return result


def minimize(
    fun: Callable,
    x0,
    *,
    bounds: Bounds = None,
    constraints: Constraints = None,
    options: Optional[dict] = None,
    sparsity: str = "dense",
    params: Any = None,
    jac_mode: str = "forward",
) -> OptimizeResult:
    """Solve ``min fun(x)`` subject to bounds and (optional) constraints.

    This is a thin wrapper around :class:`Problem`: it constructs a
    ``Problem``, calls :meth:`Problem.solve` once, and returns the result.
    For workloads that solve the *same* NLP structure many times in a row
    (e.g. closed-loop NMPC), prefer using :class:`Problem` directly so JAX
    tracing and XLA lowering are amortized across solves.

    Parameters
    ----------
    fun
        Scalar-valued objective, JAX-traceable: ``fun(x) -> scalar``, or
        ``fun(x, params) -> scalar`` when ``params`` is supplied.
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
        ``mu_strategy`` (``'adaptive'`` or ``'monotone'``; default adaptive),
        ``bound_push``, ``bound_frac``, ``sb`` (``'yes'``/``'no'``; accepted
        for cyipopt compatibility), ``warm_start`` / ``warm_start_init_point``
        (aliases), ``hessian_approximation`` (``'exact'`` or
        ``'limited-memory'``). Several Ipopt options common in cyipopt code
        (``linear_solver``, ``nlp_scaling_method``, ``acceptable_*``,
        ``limited_memory_max_history``) are accepted but ignored with a
        stderr warning so cyipopt option dicts can be reused verbatim.
    sparsity
        Sparsity strategy for Jacobian and Hessian structure:

        * ``'dense'`` (default) — every entry is declared a structural nonzero.
          Safe and simple; best for small problems.
        * ``'detect'`` — probe the Jacobian at two points and the Lagrangian
          Hessian with two random ``(sigma, lambda)`` samples, take the union
          of nonzero entries. Does not reduce JAX work, but shrinks the KKT
          system ripopt factors. For larger problems with genuine sparsity,
          this can substantially reduce iteration cost.
    params
        Optional PyTree threaded as a second positional argument into
        ``fun`` and each constraint. See :class:`Problem` for details.
    jac_mode
        ``'forward'`` (default) or ``'reverse'``. See :class:`Problem`.
    """
    problem = Problem(
        fun,
        x0,
        bounds=bounds,
        constraints=constraints,
        options=options,
        sparsity=sparsity,
        params=params,
        jac_mode=jac_mode,
    )
    return problem.solve()
