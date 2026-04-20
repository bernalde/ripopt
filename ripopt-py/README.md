# ripopt (Python)

Direct Python bindings for [ripopt](https://github.com/jkitchin/ripopt), the Rust interior-point NLP solver. Derivatives are built automatically from the user's objective and constraints via [JAX](https://jax.readthedocs.io/) autodiff, so you only have to write `f(x)` and `g(x)`.

This is the analogue of `cyipopt` for ripopt. For symbolic modeling via Pyomo, see the separate `pyomo-ripopt` package.

## Install

```bash
pip install maturin jax jaxlib numpy
maturin develop --release
```

## Example: HS071

```python
import jax.numpy as jnp
import numpy as np
from ripopt import minimize

def f(x):
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

def g(x):
    return jnp.array([
        x[0] * x[1] * x[2] * x[3],
        x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
    ])

res = minimize(
    f,
    x0=[1.0, 5.0, 5.0, 1.0],
    bounds=(1.0, 5.0),
    constraints={"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]},
    options={"tol": 1e-8},
)
print(res.x, res.fun)
```

## Reusing a problem across many solves

`minimize(...)` is a single-shot scipy-style API: every call re-traces and re-JITs `jax.hessian(lagrangian)`, which is the single most expensive step for a medium-sized problem. When you solve the *same* NLP structure many times in a row (e.g. closed-loop NMPC), use the persistent `Problem` class instead — it builds every JAX callback exactly once and reuses the XLA cache across all solves:

```python
import jax.numpy as jnp
from ripopt import Problem

def stage_cost(z, x0):
    x, u = z[:N], z[N:]
    return jnp.sum((x - r) ** 2) + 0.01 * jnp.sum(u ** 2)

def dynamics(z, x0):
    x, u = z[:N], z[N:]
    x_prev = jnp.concatenate([jnp.atleast_1d(x0[0]), x[:-1]])
    return x - A * x_prev - B * u

prob = Problem(
    stage_cost,
    x0=np.zeros(2 * N),
    bounds=(lb, ub),
    constraints={"fun": dynamics, "lb": 0.0, "ub": 0.0},
    params=jnp.array([x0_current]),   # extra traced arg threaded into fun and g
    jac_mode="reverse",                # jacrev is typically cheaper when m ≤ n
    options={"tol": 1e-6},
)

for step in range(100):
    prob.update_parameters(jnp.array([x0_current]))
    res = prob.solve()                 # warm-starts from previous x* by default
```

`examples/bench_nmpc.py` runs this same pattern 100 times and compares `minimize()` vs `Problem.solve()` — the persistent object is typically several× faster because JAX tracing and XLA lowering of the Lagrangian Hessian happen exactly once for the Problem's lifetime, not once per call.

### Warm-starting from an external solution

`Problem.solve()` accepts an explicit dual warm-start triple `(lam0, z_l0, z_u0)` — the constraint multipliers and upper/lower bound multipliers. This is the ripopt equivalent of seeding Ipopt with a cyipopt solution and lets you verify "is the reported x, λ, z_L, z_U a KKT point of this NLP?" independently of the solver's own initial-point and LS-estimate paths:

```python
# e.g. a cyipopt solution you want ripopt to verify
x_star, lam_star, z_l_star, z_u_star = cyipopt_result

res = prob.solve(
    x0=x_star,
    lam0=lam_star,
    z_l0=z_l_star,
    z_u0=z_u_star,
)
```

All three dual arrays must be supplied together (partial specs raise `ValueError`); passing any of them implies `warm_start=True`. On plain NMPC-style restarts where you simply want to reuse the previous solve's multipliers, `prob.solve(warm_start=True)` is the shorter form.

## API

`minimize(fun, x0, *, bounds=None, constraints=None, options=None, sparsity='dense', params=None, jac_mode='forward') -> OptimizeResult`

- `fun(x)` — scalar, JAX-traceable.
- `bounds` — `None`, `(lb, ub)` with scalars or arrays of length n, or a scipy-style list of `(lb_i, ub_i)` pairs.
- `constraints` — a dict `{"fun": g, "lb": ..., "ub": ...}` or a list of such dicts. Use `lb == ub` for equality.
- `options` — recognized keys: `tol`, `max_iter`, `print_level`, `constr_viol_tol`, `dual_inf_tol`, `compl_inf_tol`, `max_wall_time`, `mu_init`, `hessian_approximation` (`"exact"` or `"limited-memory"`).
- `sparsity` — `'dense'` (default) declares every entry structurally nonzero. `'detect'` probes the Jacobian and Lagrangian Hessian with two random `(x, sigma, lambda)` samples and takes the union of nonzero entries; with random multipliers, structural nonzeros are identified almost surely.

`hessian_approximation="limited-memory"` skips building the Hessian entirely and uses ripopt's L-BFGS path — useful when `jax.hessian` becomes expensive for larger problems.

`Ctrl-C` is honored: each callback polls `PyErr_CheckSignals` and the stashed `KeyboardInterrupt` is re-raised when the solve returns.

## Performance

Rough numbers for HS071 on an M-series Mac, release build, warm JIT cache:

| configuration              | total wall | solver only | iterations |
| -------------------------- | ---------- | ----------- | ---------- |
| `sparsity='dense'`         | ~155 ms    | ~32 ms      | 18         |
| `sparsity='detect'`        | ~155 ms    | ~32 ms      | 18         |

Raw ripopt (Rust) on the same problem is ~3 ms. The gap is per-callback JAX dispatch (~300 µs × 5 callbacks × 18 iterations ≈ 27 ms); for problems where ripopt dominates the wall clock this gap is amortized. Run `python examples/bench_hs071.py` to reproduce.

## Documentation

- [`examples/tutorial.ipynb`](examples/tutorial.ipynb) — runnable notebook walking through unconstrained, box-constrained, HS071, options, L-BFGS, sparsity detection, and multi-block constraints
- [`docs/user_guide.md`](docs/user_guide.md) — task-oriented prose guide (JAX tips, constraint patterns, troubleshooting)
- [`docs/api.md`](docs/api.md) — strict reference for `minimize`, `OptimizeResult`, and the options dict

## Status

Early prototype. Sparsity detection uses a two-probe random-multiplier union — safe for structurally-honest code but defeated by rare numerical cancellations. A wheel is built by `maturin build --release`; cross-platform wheels, PyPI release, and a symbolic sparsity path are future work.
