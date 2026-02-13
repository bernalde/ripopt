# ripopt

[![Tests](https://github.com/jkitchin/ripopt/actions/workflows/test.yml/badge.svg)](https://github.com/jkitchin/ripopt/actions/workflows/test.yml)

![img](./ipopt-rust.png)

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
- Gauss-Newton restoration phase with adaptive Levenberg-Marquardt regularization
- Watchdog strategy for escaping narrow feasible corridors
- Automatic NE-to-LS reformulation for overdetermined nonlinear equation systems
- NLP scaling (gradient-based objective and constraint scaling)
- Condensed KKT system for problems with many more constraints than variables
- Local infeasibility detection for inconsistent constraint systems

## Benchmarks

### Hock-Schittkowski Test Suite (120 problems)

| Metric | ripopt | cyipopt (Ipopt + MUMPS) |
|--------|--------|-------------------------|
| Problems solved | **119/120 (99.2%)** | 118/120 (98.3%) |
| Matching solutions (rel diff < 1e-4) | 107/117 | -- |
| Mean iterations (where both solve) | 14.5 | 14.5 |
| Geometric mean speedup | **116x faster** | -- |

ripopt solves one more problem than cyipopt (TP116, TP223) while failing only on TP374 (trigonometric Chebyshev problem). See `hs_suite/HS_VALIDATION_REPORT.md` for the full comparison.

### CUTEst Test Set (727 problems)

| Metric | ripopt | Ipopt (C++ with MUMPS) |
|--------|--------|------------------------|
| Problems solved | **552/727 (75.9%)** | 557/727 (76.6%) |
| Both solved | 519/727 | -- |
| Matching solutions (rel diff < 1e-4) | 431/519 | -- |

ripopt is within 5 problems of Ipopt on the full CUTEst benchmark. 99 overdetermined nonlinear equation problems are correctly classified as locally infeasible via automatic NE-to-LS reformulation. See `cutest_suite/CUTEST_REPORT.md` for the full comparison.

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
| `acceptable_tol` | 1e-4 | Acceptable (less strict) tolerance |
| `acceptable_iter` | 10 | Consecutive acceptable iterations needed |
| `mu_init` | 0.1 | Initial barrier parameter |
| `print_level` | 5 | Output verbosity (0=silent, 5=verbose) |
| `mu_strategy_adaptive` | true | Adaptive vs monotone barrier update |
| `max_soc` | 4 | Maximum second-order correction steps |
| `max_wall_time` | 0.0 | Wall-clock time limit in seconds (0=no limit) |
| `warm_start` | false | Enable warm-start initialization |
| `constr_viol_tol` | 1e-4 | Constraint violation tolerance |
| `dual_inf_tol` | 1.0 | Dual infeasibility tolerance |

### Result

`SolveResult` contains:

- `x` -- optimal primal variables
- `objective` -- optimal objective value f(x*)
- `constraint_multipliers` -- Lagrange multipliers for constraints (y)
- `bound_multipliers_lower` / `bound_multipliers_upper` -- bound multipliers (z_L, z_U)
- `constraint_values` -- constraint values g(x*)
- `status` -- one of: `Optimal`, `Acceptable`, `Infeasible`, `LocalInfeasibility`, `MaxIterations`, `NumericalError`, `Unbounded`, `RestorationFailed`, `InternalError`
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

113 tests total:
- **77 unit tests**: Dense LDL factorization, convergence checking, filter line search, fraction-to-boundary, KKT assembly, restoration
- **21 integration tests**: Rosenbrock, SimpleQP, HS071, HS035, PureBoundConstrained, MultipleEqualityConstraints, NE-to-LS reformulation, and more
- **15 HS regression tests**: Selected Hock-Schittkowski problems for regression checking

## Architecture

```
src/
  lib.rs              Public API (solve function, re-exports)
  problem.rs          NlpProblem trait definition
  options.rs          SolverOptions with Ipopt-matching defaults
  result.rs           SolveResult and SolveStatus
  ipm.rs              Main IPM loop, barrier updates, line search, NE-to-LS detection, NLP scaling
  kkt.rs              KKT system assembly, solution, and inertia correction
  convergence.rs      Convergence checking (primal/dual/complementarity)
  filter.rs           Filter line search mechanism
  linear_solver/
    mod.rs            LinearSolver trait, SymmetricMatrix
    dense.rs          Dense LDL^T (Bunch-Kaufman) factorization
  restoration.rs      Gauss-Newton restoration phase with adaptive LM regularization
  warmstart.rs        Warm-start initialization

tests/
  correctness.rs      Integration tests (21 NLP problems)
  hs_regression.rs    HS suite regression tests (15 problems)

examples/
  rosenbrock.rs       Unconstrained optimization
  hs071.rs            Constrained NLP
  benchmark.rs        Timing benchmark
```

## Algorithm Details

### Core Interior Point Method

The solver follows the primal-dual barrier method from the Ipopt papers (Wachter & Biegler, 2006). At each iteration it:

1. Assembles and factors the KKT system using dense LDL^T with Bunch-Kaufman pivoting
2. Computes inertia of the factorization and applies regularization if needed
3. Computes search directions with iterative refinement (up to 3 rounds)
4. Applies second-order corrections (SOC) when the initial step is rejected
5. Uses a filter line search with backtracking to ensure sufficient progress
6. Updates the barrier parameter adaptively based on complementarity

### Condensed KKT System

For problems where the number of constraints m exceeds 2n, the solver automatically uses a condensed (Schur complement) formulation. This reduces the factorization cost from O((n+m)^3) to O(n^2 m + n^3), enabling efficient handling of problems with many constraints and few variables.

### NE-to-LS Reformulation

When the solver detects an overdetermined nonlinear equation system (m > n, f(x) = 0, all equality constraints, starting point not already feasible), it automatically reformulates the problem as unconstrained least-squares minimization:

```
min  0.5 * ||g(x) - target||^2
```

using a Gauss-Newton Hessian approximation (H = J^T J). If the residual is small at the solution, the original system is consistent and `Optimal` is reported. Otherwise, `LocalInfeasibility` is reported with the best least-squares solution.

### Restoration Phase

When the filter line search fails, the solver enters a Gauss-Newton restoration phase that minimizes constraint violation ||g(x) - target||^2. This uses:

- Adaptive Levenberg-Marquardt regularization (scaling from 1e-8 up to 1e-2)
- Penalty-regularized fallback when the standard GN step fails
- Gradient descent fallback when GN line search fails

### Watchdog Strategy

The solver implements a watchdog mechanism that temporarily relaxes the filter acceptance criteria when progress stalls due to shortened steps. This helps escape narrow feasible corridors where strict Armijo conditions are too conservative.

## Implementation Status

**Version: 0.2.0-dev**

### What works

- Unconstrained, equality-constrained, inequality-constrained, and bound-constrained NLP
- **119/120 Hock-Schittkowski problems solved** (99.2%), surpassing cyipopt (118/120)
- **552/727 CUTEst problems solved** (75.9%), within 5 of Ipopt (557/727)
- Automatic detection and reformulation of overdetermined nonlinear equations
- NLP scaling (gradient-based objective and constraint scaling)
- Condensed KKT for problems with many more constraints than variables
- Iterative refinement for KKT solves
- Gauss-Newton restoration phase with gradient descent fallback
- Dense problems up to a few hundred variables

### Known limitations

- **Dense linear solver only** -- no sparse support, which limits practical problem size to ~500 variables
- **No quasi-Newton** (L-BFGS) Hessian approximation -- user must supply exact Hessians
- **No Mehrotra predictor-corrector** -- single centrality correction only

### Not yet implemented (present in Ipopt)

- Sparse linear solvers (MA27, MA57, MUMPS, Pardiso)
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
