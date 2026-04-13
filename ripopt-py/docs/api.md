# API Reference

## `ripopt.minimize`

```python
ripopt.minimize(
    fun: Callable,
    x0: ArrayLike,
    *,
    bounds: Bounds = None,
    constraints: Constraints = None,
    options: dict | None = None,
    sparsity: str = "dense",
) -> OptimizeResult
```

Solve

```
min   f(x)
 x
s.t.  g_l ≤ g(x) ≤ g_u
      x_l ≤   x  ≤ x_u
```

All derivatives are built from `fun` and the constraint functions via JAX autodiff; the user never hand-writes a gradient, Jacobian, or Hessian.

### Parameters

| name | type | description |
|---|---|---|
| `fun` | `Callable[[jax.Array], scalar]` | Objective. Must be JAX-traceable: written with `jax.numpy` ops on the input `x`. |
| `x0` | array-like of length `n` | Initial primal point. Cast to `float64`. |
| `bounds` | see below | Variable bounds `x_l ≤ x ≤ x_u`. |
| `constraints` | see below | Nonlinear constraints on `g(x)`. |
| `options` | `dict` | Solver options; see **Options** below. |
| `sparsity` | `"dense"` \| `"detect"` | Sparsity strategy for the Jacobian and Lagrangian Hessian. |

### `bounds`

One of:

- `None` — no bounds (both sides `±inf`).
- `(lb, ub)` — two scalars or two length-`n` arrays. Broadcast if scalar. Use `np.inf` / `-np.inf` for missing sides.
- A length-`n` sequence of `(lb_i, ub_i)` pairs (scipy-style). `None` inside a pair means no bound on that side.

### `constraints`

One of:

- `None` — unconstrained (beyond bounds).
- A **single dict** `{"fun": g, "lb": ..., "ub": ...}` where `g(x)` returns a scalar or a 1-D array. `lb` and `ub` may be scalars or arrays matching `g(x)`. Use `lb == ub` for equality constraints.
- A **list of such dicts**, each defining a block of constraints. The blocks are concatenated into a single `g` vector internally.

Example — one inequality and one equality:

```python
def g(x):
    return jnp.array([
        x[0] * x[1] * x[2] * x[3],              # inequality
        x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,  # equality
    ])

constraints = {"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]}
```

### `options`

Recognized keys:

| key | type | default | meaning |
|---|---|---|---|
| `tol` | float | `1e-8` | Overall convergence tolerance. |
| `max_iter` | int | `3000` | Maximum IPM iterations. |
| `print_level` | int (0–5) | `0` | Verbosity; `0` is silent. |
| `constr_viol_tol` | float | | Primal infeasibility tolerance. |
| `dual_inf_tol` | float | | Dual infeasibility tolerance. |
| `compl_inf_tol` | float | | Complementarity tolerance. |
| `max_wall_time` | float | `0.0` | Time budget in seconds; `0` = none. |
| `mu_init` | float | | Initial barrier parameter. |
| `hessian_approximation` | `"exact"` \| `"limited-memory"` | `"exact"` | `"limited-memory"` skips Hessian construction entirely. |

Unrecognized keys raise `ValueError`.

### `sparsity`

- `"dense"` (default) — every entry of the `(m, n)` Jacobian and lower triangle of the `(n, n)` Hessian is declared a structural nonzero. Safe and simple.
- `"detect"` — the wrapper probes the Jacobian at `x0` and a random perturbation, and the Lagrangian Hessian at two random `(σ, λ)` samples, and takes the union of nonzero entries. With random multipliers the probability of missing a structural nonzero is zero in theory; in practice, pathological cancellations are rare.

Sparsity detection shrinks the KKT system ripopt factors; it does not reduce the per-callback JAX work (JAX still computes a dense Jacobian/Hessian).

### Returns

An `OptimizeResult` — see below.

### Raises

- `ValueError` for malformed inputs or unknown option keys.
- `KeyboardInterrupt` if `Ctrl-C` is hit during the solve (stashed in a callback and re-raised on return).
- Any `PyErr` raised inside the user's `fun` or constraint callables is captured and re-raised after the solve.

---

## `ripopt.OptimizeResult`

`@dataclass`. Returned by `minimize`.

| field | type | meaning |
|---|---|---|
| `x` | `np.ndarray` shape `(n,)` | Primal optimum `x*`. |
| `fun` | `float` | `f(x*)`. |
| `status` | `str` | `"Optimal"`, `"Infeasible"`, `"LocalInfeasibility"`, `"MaxIterations"`, `"NumericalError"`, `"Unbounded"`, `"RestorationFailed"`, `"InternalError"`. |
| `success` | `bool` | `status == "Optimal"`. |
| `iterations` | `int` | IPM iterations used. |
| `constraint_multipliers` | `np.ndarray` shape `(m,)` | Lagrange multipliers for `g(x) = g_l` / `g_u`. |
| `bound_multipliers_lower` | `np.ndarray` shape `(n,)` | Multipliers for active lower variable bounds. |
| `bound_multipliers_upper` | `np.ndarray` shape `(n,)` | Multipliers for active upper variable bounds. |
| `constraint_values` | `np.ndarray` shape `(m,)` | `g(x*)`. |
| `wall_time_secs` | `float` | Total ripopt wall time, including time spent inside callbacks. |
| `final_primal_inf` | `float` | Primal infeasibility at the final iterate. |
| `final_dual_inf` | `float` | Dual infeasibility at the final iterate. |
| `final_compl` | `float` | Complementarity slack at the final iterate. |
