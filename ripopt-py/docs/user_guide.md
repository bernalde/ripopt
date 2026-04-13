# User Guide

`ripopt` is a direct Python interface to the Rust interior-point NLP solver of the same name. You write only the **objective** and the **constraint functions**, both as JAX-traceable Python; gradients, Jacobians, and the Lagrangian Hessian are generated automatically via `jax.grad`, `jax.jacfwd`, and `jax.hessian`.

This guide is task-oriented. For a strict reference, see [`api.md`](api.md); for runnable examples, see [`../examples/tutorial.ipynb`](../examples/tutorial.ipynb).

---

## Installation

```bash
cd ripopt-py
maturin develop --release
```

The release build ships ripopt with MUMPS-style sparse LDLᵀ as the default linear solver and is what the numbers in the README are measured against. A debug build works but is several times slower.

Dependencies: `numpy`, `jax`, `jaxlib`. JAX is configured to use `float64` on import.

## The problem

`minimize` solves

```
min   f(x)
 x
s.t.  g_l ≤ g(x) ≤ g_u        (constraints, optional)
      x_l ≤   x  ≤ x_u        (bounds, optional)
```

with `x ∈ ℝⁿ`, `g: ℝⁿ → ℝᵐ`. Equality constraints are encoded as `lb == ub`.

## Writing the objective

The objective is any Python callable `f(x)` that returns a scalar and is traceable by JAX:

```python
import jax.numpy as jnp

def rosenbrock(x):
    return (1.0 - x[0])**2 + 100.0 * (x[1] - x[0]**2)**2
```

**Rules of thumb**:

- Use `jax.numpy` (`jnp`) operations on `x`, not `numpy` — JAX needs to trace through arithmetic to differentiate.
- Avoid Python `if` on tracer values. Use `jnp.where` / `jax.lax.cond` instead.
- Indexing (`x[0]`, `x[:3]`) is fine.
- Returning `0-D` (scalar) is required.

If you see `TracerBoolConversionError` or similar, you've used a regular `if` on a JAX tracer — replace it with `jnp.where`.

## Starting point

`x0` is a plain Python list or NumPy array of length `n`. It's cast to contiguous `float64` internally. Choose a point that is at least well-defined (no divisions by zero, no `log` of non-positive, etc.) because the wrapper probes `fun` at `x0` during setup.

## Bounds

Three forms:

```python
# No bounds
minimize(f, x0)

# Scalar broadcast
minimize(f, x0, bounds=(-1.0, 1.0))

# Per-variable
minimize(f, x0, bounds=(np.zeros(n), np.full(n, np.inf)))

# scipy-style list of (lb, ub) pairs, with None for missing sides
minimize(f, x0, bounds=[(0.0, None), (-1.0, 1.0), (None, 5.0)])
```

Use `np.inf` / `-np.inf` (or `None` in the list form) for missing sides. ripopt treats any bound with magnitude greater than `1e19` as "no bound".

## Constraints

A constraint block is a dict:

```python
{"fun": g, "lb": ..., "ub": ...}
```

where `g(x)` returns a scalar or 1-D array and `lb` / `ub` are scalars or arrays matching the length of `g(x)`.

**Single block** — the common case:

```python
def g(x):
    return jnp.array([
        x[0] * x[1] * x[2] * x[3],             # >= 25
        x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2, # == 40
    ])

constraints = {"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]}
```

Equality constraints are written `lb == ub` — no separate equality argument.

**Multiple blocks** — useful when constraints come from different physics / layers of the model. The wrapper concatenates them internally:

```python
constraints = [
    {"fun": power_balance, "lb": 0.0, "ub": 0.0},          # equality
    {"fun": thermal_limits, "lb": -np.inf, "ub": 1.0},     # inequality
    {"fun": ramp_rates,     "lb": -0.1,    "ub": 0.1},
]
```

This is equivalent to writing one `g` that stacks the three blocks, just more readable.

## Options

Pass an options dict to `minimize`:

```python
minimize(f, x0, options={
    "tol": 1e-10,
    "max_iter": 200,
    "print_level": 5,                     # per-iteration log to stderr
    "hessian_approximation": "exact",     # or "limited-memory"
})
```

The full list of keys is in [`api.md`](api.md). Unknown keys raise `ValueError`.

For a first-pass diagnosis of a failing solve, set `print_level=5` and watch the dual / primal infeasibility columns. If dual infeas plateaus above tolerance while primal infeas is fine, your Hessian is likely ill-conditioned near the optimum and L-BFGS mode may help.

## Sparsity modes

By default (`sparsity="dense"`) the wrapper tells ripopt that every Jacobian and Hessian entry is structurally nonzero. For small-to-medium problems this is fine — the JAX work is identical either way, and ripopt's KKT factorization is fast.

For genuinely sparse problems, pass `sparsity="detect"`:

```python
minimize(f, x0, constraints=..., sparsity="detect")
```

The wrapper evaluates the Jacobian at `x0` and at `x0 + perturb`, and the Lagrangian Hessian with two random `(σ, λ)` samples, then unions the nonzero patterns. Random multipliers make structural cancellations probability-zero events in theory.

**When this helps**: problems where the Jacobian or Hessian has genuine structural zeros — e.g., a chain of link-local constraints, a Lagrangian whose Hessian is a banded or block-diagonal matrix.

**When this does not help**: small problems (HS071 and smaller), or problems where the JAX callback cost dominates the KKT factorization cost. You can run with both modes and compare `res.wall_time_secs`.

**Caveat**: the JAX Jacobian / Hessian are still computed densely each callback. The sparsity reduction only shrinks what ripopt stores and factors. A future version will plumb `jax.experimental.sparse` through to JAX as well.

## L-BFGS Hessian approximation

When the Hessian is expensive — either because `n` is large or because the Lagrangian contains a lot of nonlinear coupling — skip it:

```python
minimize(f, x0, options={"hessian_approximation": "limited-memory"})
```

With this setting the wrapper never constructs `jax.hessian(lagrangian)`, and ripopt uses an L-BFGS approximation internally. You only need `f` and the constraint functions; gradients and the Jacobian are still computed exactly via JAX.

L-BFGS is typically 2–5× more iterations than exact-Hessian mode, but each iteration is much cheaper. For problems with `n ≥ ~500` or Hessians that take longer than a few milliseconds, it's usually the right choice.

## The result

`minimize` returns an `OptimizeResult` dataclass. See [`api.md`](api.md) for the full field list.

Common patterns:

```python
res = minimize(...)

if not res.success:
    raise RuntimeError(f"solve failed: {res.status}")

x_opt = res.x
lam   = res.constraint_multipliers       # g multipliers
z_l   = res.bound_multipliers_lower      # x >= x_l multipliers
z_u   = res.bound_multipliers_upper      # x <= x_u multipliers

# Which bounds are active?
active_lower = z_l > 1e-6
active_upper = z_u > 1e-6
```

The `wall_time_secs` field is ripopt's internal measurement — it includes the time spent inside your JAX callbacks.

## Troubleshooting

**`TracerBoolConversionError`** — you have a Python `if` on a JAX tracer. Replace with `jnp.where(cond, a, b)`.

**`MaxIterations` status** — check `res.final_primal_inf` and `res.final_dual_inf`. If dual infeas is the holdout, try `options={"hessian_approximation": "limited-memory"}`. If primal infeas is the holdout, your constraints may be ill-posed (redundant, nearly-infeasible, or badly scaled).

**`NumericalError` status** — usually a singular KKT matrix. Check that your Hessian is correct (evaluate it at `x0` and look at its eigenvalues). Falling back to L-BFGS mode often works.

**Slow per-iteration** — JAX callback dispatch is ~100–300 µs per call. ripopt calls callbacks roughly 5× per iteration. For an 18-iteration HS071 solve this is ~30 ms of overhead. If your problem has hundreds of iterations, this dominates; consider `sparsity="detect"` (if applicable) and `hessian_approximation="limited-memory"` to reduce the number of callback invocations.

**Different answers from cyipopt** — ripopt's defaults match ipopt's where possible, but not everywhere. Set `options={"tol": 1e-10}` first; if they still differ significantly, the problem may have multiple local minima.

**`KeyboardInterrupt` not firing during a long solve** — make sure you rebuilt the extension after signal-handling was added. Each callback polls `PyErr_CheckSignals`; if you built before that change, the extension won't check.

## Relationship to other Python interfaces

| package | interface | user writes |
|---|---|---|
| `cyipopt` | direct, callback-based | `f`, `∇f`, `g`, `∇g`, `∇²L` |
| `pyomo-ripopt` | symbolic (Pyomo) | an algebraic model; derivatives via Pyomo's AD |
| `ripopt` (this package) | direct, callback-based + autodiff | `f`, `g` |

Use `pyomo-ripopt` when you already have a Pyomo model. Use this package when you want to stay in a JAX-native workflow (gradients for use in larger JAX pipelines, mixing optimization with neural nets, JIT-compiled models).
