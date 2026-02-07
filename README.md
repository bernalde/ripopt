# ripopt

A memory-safe interior point optimizer written in Rust, inspired by [Ipopt](https://github.com/coin-or/Ipopt).

## What is this?

ripopt solves nonlinear programming (NLP) problems of the form:

```
min  f(x)
s.t. g_l <= g(x) <= g_u
     x_l <= x    <= x_u
```

It implements a primal-dual interior point method with a barrier formulation, similar to the algorithm described in the Ipopt papers. The solver is written entirely in Rust with no external C/Fortran dependencies.

## Features

- Primal-dual interior point method with logarithmic barrier
- Dense LDL^T factorization via Bunch-Kaufman pivoting with inertia detection
- Filter line search with switching condition and Armijo criterion
- Second-order corrections (SOC) for improved step acceptance
- Adaptive and monotone barrier parameter strategies
- Fraction-to-boundary rule for primal and dual step sizes
- Support for equality constraints, inequality constraints, and variable bounds
- Warm-start initialization
- Gauss-Newton restoration phase with Levenberg-Marquardt regularization

## Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
ripopt = { path = "." }
```

### Defining a Problem

Implement the `NlpProblem` trait:

```rust
use ripopt::NlpProblem;

struct Rosenbrock;

impl NlpProblem for Rosenbrock {
    fn num_variables(&self) -> usize { 2 }
    fn num_constraints(&self) -> usize { 0 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        // Unconstrained: use infinity bounds
        for i in 0..2 {
            x_l[i] = f64::NEG_INFINITY;
            x_u[i] = f64::INFINITY;
        }
    }

    fn constraint_bounds(&self, _g_l: &mut [f64], _g_u: &mut [f64]) {}
    fn initial_point(&self, x0: &mut [f64]) { x0[0] = -1.0; x0[1] = 1.0; }

    fn objective(&self, x: &[f64]) -> f64 {
        100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2)
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
        grad[1] = 200.0 * (x[1] - x[0] * x[0]);
    }

    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {}
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) { (vec![], vec![]) }
    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {}

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0, 1, 1], vec![0, 0, 1])  // lower triangle
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        vals[0] = obj_factor * (-400.0 * (x[1] - 3.0 * x[0] * x[0]) + 2.0);
        vals[1] = obj_factor * (-400.0 * x[0]);
        vals[2] = obj_factor * 200.0;
    }
}
```

### Solving

```rust
use ripopt::{SolverOptions, solve};

let problem = Rosenbrock;
let options = SolverOptions::default();
let result = solve(&problem, &options);

println!("Status: {:?}", result.status);
println!("Objective: {:.6e}", result.objective);
println!("Solution: {:?}", result.x);
println!("Iterations: {}", result.iterations);
```

### Solver Options

Key options (all have Ipopt-matching defaults):

| Option | Default | Description |
|--------|---------|-------------|
| `tol` | 1e-8 | Convergence tolerance |
| `max_iter` | 3000 | Maximum iterations |
| `acceptable_tol` | 1e-6 | Acceptable (less strict) tolerance |
| `acceptable_iter` | 15 | Consecutive acceptable iterations needed |
| `mu_init` | 0.1 | Initial barrier parameter |
| `print_level` | 5 | Output verbosity (0=silent, 5=verbose) |
| `mu_strategy_adaptive` | true | Adaptive vs monotone barrier update |
| `max_soc` | 4 | Maximum second-order correction steps |

### Result

`SolveResult` contains:

- `x` -- optimal primal variables
- `objective` -- optimal objective value f(x*)
- `constraint_multipliers` -- Lagrange multipliers for constraints (y)
- `bound_multipliers_lower` / `bound_multipliers_upper` -- bound multipliers (z_L, z_U)
- `constraint_values` -- constraint values g(x*)
- `status` -- one of: `Optimal`, `Acceptable`, `Infeasible`, `MaxIterations`, `NumericalError`, `Unbounded`, `RestorationFailed`, `InternalError`
- `iterations` -- number of IPM iterations

## Examples

```bash
# Rosenbrock function (unconstrained with bounds)
cargo run --example rosenbrock

# HS071 (constrained NLP with inequalities)
cargo run --example hs071

# Benchmark timing across 5 problems
cargo run --release --example benchmark
```

## Tests

```bash
cargo test
```

23 tests total:
- **17 unit tests**: Dense LDL factorization (7), convergence checking (4), filter line search (5), fraction-to-boundary (1)
- **6 integration tests**: Rosenbrock, SimpleQP, HS071, HS035, PureBoundConstrained, MultipleEqualityConstraints

## Validation Against Ipopt (cyipopt)

Validated against the **120-problem Hock-Schittkowski test suite**, comparing ripopt with cyipopt (Python wrapper for Ipopt with MUMPS):

| Metric | ripopt | cyipopt |
|--------|--------|---------|
| Problems solved | **119/120 (99.2%)** | 118/120 (98.3%) |
| Matching solutions (rel diff < 1e-4) | 107/117 | — |
| Mean iterations (where both solve) | 14.5 | 14.5 |

ripopt solves one more problem than cyipopt (TP116, TP223) while failing only on TP374 (trigonometric Chebyshev problem with no variable bounds).

See `hs_suite/HS_VALIDATION_REPORT.md` for the full comparison report.

## Architecture

```
src/
  lib.rs              Public API (solve function, re-exports)
  problem.rs          NlpProblem trait definition
  options.rs          SolverOptions with Ipopt-matching defaults
  result.rs           SolveResult and SolveStatus
  ipm.rs              Main IPM loop, barrier updates, line search
  kkt.rs              KKT system assembly and solution
  convergence.rs      Convergence checking (primal/dual/complementarity)
  filter.rs           Filter line search mechanism
  linear_solver/
    mod.rs            LinearSolver trait, SymmetricMatrix
    dense.rs          Dense LDL^T (Bunch-Kaufman) factorization
  restoration.rs      Gauss-Newton restoration phase
  warmstart.rs        Warm-start initialization

tests/
  correctness.rs      Integration tests (6 NLP problems)

examples/
  rosenbrock.rs       Unconstrained optimization
  hs071.rs            Constrained NLP
  benchmark.rs        Timing benchmark across 5 problems

validation/
  compare_cyipopt.py  Solve same problems with cyipopt for comparison
  generate_report.py  Generate markdown validation report
```

## Implementation Status

**Version: 0.1.0**

### What works

- Unconstrained, equality-constrained, inequality-constrained, and bound-constrained NLP
- **119/120 Hock-Schittkowski problems solved** (99.2%), surpassing cyipopt (118/120)
- Iterative refinement for KKT solves
- Gauss-Newton restoration phase with gradient descent fallback
- Dense problems up to a few hundred variables

### Known limitations

- **Dense linear solver only** -- no sparse support, which limits practical problem size
- **No NLP scaling** -- gradient-based or user-provided scaling not implemented
- **No quasi-Newton** (L-BFGS) Hessian approximation -- user must supply exact Hessians

### Not yet implemented (present in Ipopt)

- Sparse linear solvers (MA27, MA57, MUMPS, Pardiso)
- NLP scaling
- Mehrotra predictor-corrector
- Quasi-Newton Hessian approximation (L-BFGS)

## Sign Convention

ripopt uses the Ipopt convention where the Lagrangian is:

```
L = f(x) + y^T g(x)
```

For inequality constraints `g(x) >= g_l`, the multiplier `y` is negative at optimality. For equality constraints, `y` can be positive or negative.

## License

[EPL-2.0](LICENSE) (Eclipse Public License 2.0), consistent with Ipopt.
