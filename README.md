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

It implements a primal-dual interior point method with a barrier formulation, similar to the algorithm described in the Ipopt papers. The solver is written entirely in Rust (~9,800 lines) with no external C/Fortran dependencies.

## Features

- Primal-dual interior point method with logarithmic barrier
- Dense LDL^T factorization via Bunch-Kaufman pivoting with inertia detection
- Sparse LDL^T factorization (via `faer`) for larger problems (n+m >= 100)
- Filter line search with switching condition and Armijo criterion
- Second-order corrections (SOC) for improved step acceptance
- Adaptive and monotone barrier parameter strategies
- Fraction-to-boundary rule for primal and dual step sizes
- Support for equality constraints, inequality constraints, and variable bounds
- Warm-start initialization
- Two-phase restoration: fast Gauss-Newton + full NLP restoration subproblem
- Multi-attempt recovery with systematic barrier landscape perturbation
- Watchdog strategy for escaping narrow feasible corridors
- Automatic NE-to-LS reformulation for overdetermined nonlinear equation systems
- NLP scaling (gradient-based objective and constraint scaling)
- Condensed KKT system for problems with many more constraints than variables
- Local infeasibility detection for inconsistent constraint systems
- **Multi-solver fallback architecture**: L-BFGS, Augmented Lagrangian, and explicit slack reformulation

## Benchmarks

### Hock-Schittkowski Test Suite (120 problems)

| Metric          | ripopt             | Ipopt (native, MUMPS) |
|-----------------|--------------------|-----------------------|
| Problems solved | **120/120 (100%)** | 118/120 (98.3%)       |
| Optimal         | 97                 | 118                   |
| Acceptable      | 23                 | 0                     |

ripopt solves all 120 problems. Ipopt fails on TP214 (`InvalidNumberDetected`) and TP223 (declared infeasible despite being feasible -- ripopt solves it in 4 iterations).

### CUTEst Benchmark Suite (727 problems)

| Metric        | ripopt              | Ipopt (C++ with MUMPS) |
|---------------|---------------------|------------------------|
| Total solved  | **599/727 (82.4%)** | 556/727 (76.5%)        |
| Constrained   | **377/493**         | 340/493                |
| Unconstrained | **222/234**         | 216/234                |
| Both solve    | 556                 | 556                    |
| ripopt only   | **43**              | --                     |
| Ipopt only    | --                  | **0**                  |
| Both fail     | 128                 | 128                    |

ripopt solves **43 more problems** than Ipopt, and solves everything Ipopt solves.

#### Speed comparison

**Methodology.** Both solvers are timed on the same CUTEst problems using the same evaluation callbacks. For Ipopt, only the `IpoptSolve` call is timed (excluding `CreateIpoptProblem`, option setting, and `FreeIpoptProblem`). For ripopt, the full `solve()` API is timed (including internal scaling, NE-to-LS detection, and all fallbacks). Each problem runs 3 times with the best time reported. Both use adaptive mu, tol=1e-8, max_iter=3000.

**On 272 problems where both solvers take >= 0.1ms** (excluding trivially fast sub-0.1ms solves where measurement noise dominates):

| Metric                          | Value            |
|---------------------------------|------------------|
| Median speedup (ripopt faster)  | **3.3x**         |
| Geometric mean speedup          | **2.2x**         |
| Problems where ripopt is faster | 186/272 (68%)    |
| Problems where Ipopt is faster  | 86/272 (32%)     |
| ripopt 10x+ faster              | 76/272 (28%)     |

Speed by problem size (problems >= 0.1ms):

| Problem size          | Count | Median speedup |
|-----------------------|-------|----------------|
| Small (n <= 10)       | 176   | **4.3x**       |
| Medium (10 < n <= 50) | 64    | **3.7x**       |
| Large (n > 50)        | 32    | **1.6x**       |

**Interpreting the speed numbers.** The raw median across all 556 mutually-solved problems is 30x, but this is misleading. Most CUTEst problems are small (n < 10) and solve in microseconds for ripopt, while Ipopt has a ~1-3ms floor from internal initialization (MA27 symbolic analysis, TNLP setup) that dominates for trivial problems. The 3.3x median on non-trivial problems is a more meaningful comparison.

The speed advantage comes from three sources:

1. **Lower per-iteration overhead.** ripopt's dense Bunch-Kaufman factorization avoids sparse symbolic analysis and has minimal allocation. For small-to-medium problems (n < 50), this gives 2-5x per-iteration speedup.
2. **Condensed KKT for over-constrained problems.** When m >> n, ripopt reduces an (n+m)x(n+m) factorization to nxn, giving dramatic speedups on problems like PT (67x, n=2, m=501) and SIPOW (5-8x, n=2-4, m=2000).
3. **Fewer iterations on some problems.** NE-to-LS reformulation, better restoration recovery, and adaptive strategies sometimes converge in fewer iterations.

Where Ipopt is faster:

1. **Large sparse problems.** Ipopt's MUMPS/MA27 scales better than ripopt's dense solver for n > 50, and better than ripopt's faer-based sparse solver for some sparsity patterns.
2. **Problems requiring many iterations.** When both take 100+ iterations, Ipopt's mature sparse linear algebra can win on per-iteration cost for larger systems.
3. **Some difficult nonlinear problems.** Ipopt's extensive tuning of barrier parameter updates and restoration gives it an edge on specific hard problems.

#### Large problems (n+m >= 100)

On 48 problems with n+m >= 100 (exercising the sparse LDL solver):

| Metric         | Value        |
|----------------|--------------|
| Both solve     | 48/48 (100%) |
| ripopt faster  | 29/48 (60%)  |
| Median speedup | 1.6x         |

#### Attribution of ripopt-only solves (43 problems)

- **NE-to-LS reformulation** (~10): Nonlinear equation systems where Ipopt returns status -10
- **Two-phase restoration** (~10): GN phase or multi-attempt recovery succeeds where Ipopt's restoration fails
- **Pragmatic inertia correction** (~5): ripopt proceeds with approximate factorization instead of failing
- **Explicit slack fallback** (~3): Stabilizes multiplier oscillation at degenerate points
- **L-BFGS / AL fallback** (~5): L-BFGS for unconstrained, augmented Lagrangian for constrained
- **Best-du tracking** (~5): Recovers acceptable solutions from cycling/stalling at max iterations
- **Other algorithmic differences** (~5): Barrier parameter, filter, convergence criteria

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

138 tests total:
- **102 unit tests**: Dense LDL factorization, convergence checking, filter line search, fraction-to-boundary, KKT assembly, restoration, linear solver
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
  restoration.rs      Gauss-Newton restoration phase with adaptive LM regularization
  restoration_nlp.rs  Full NLP restoration subproblem (Phase 2)
  lbfgs.rs            L-BFGS solver for unconstrained/bound-constrained problems
  augmented_lagrangian.rs  Augmented Lagrangian fallback for constrained problems
  slack_formulation.rs     Explicit slack reformulation fallback
  warmstart.rs        Warm-start initialization
  linear_solver/
    mod.rs            LinearSolver trait, SymmetricMatrix
    dense.rs          Dense LDL^T (Bunch-Kaufman) factorization

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

1. Assembles and factors the KKT system using dense LDL^T (Bunch-Kaufman) or sparse LDL^T (faer)
2. Computes inertia of the factorization and applies regularization if needed
3. Computes search directions with iterative refinement (up to 3 rounds)
4. Applies second-order corrections (SOC) when the initial step is rejected
5. Uses a filter line search with backtracking to ensure sufficient progress
6. Updates the barrier parameter adaptively based on complementarity

### Multi-Solver Fallback Architecture

When the primary IPM fails, ripopt automatically tries alternative solvers:

1. **L-BFGS**: Tried first for unconstrained problems (m=0, no bounds); used as fallback for bound-constrained problems
2. **Augmented Lagrangian**: PHR penalty method for constrained problems, with the IPM solving each AL subproblem
3. **Explicit slack reformulation**: Converts g(x) to g(x)-s=0 with bounds on s, stabilizing multiplier oscillation at degenerate points
4. **Best-du tracking**: Throughout the solve, tracks the iterate with lowest dual infeasibility and recovers it at max iterations

### Condensed KKT System

For problems where the number of constraints m exceeds 2n, the solver automatically uses a condensed (Schur complement) formulation. This reduces the factorization cost from O((n+m)^3) to O(n^2 m + n^3), enabling efficient handling of problems with many constraints and few variables.

### NE-to-LS Reformulation

When the solver detects an overdetermined nonlinear equation system (m >= n, f(x) = 0, all equality constraints, starting point not already feasible), it automatically reformulates the problem as unconstrained least-squares minimization:

```
min  0.5 * ||g(x) - target||^2
```

using a full Hessian (J^T J + sum of r_i * nabla^2 g_i). If the residual is small at the solution, the original system is consistent and `Optimal` is reported. Otherwise, `LocalInfeasibility` is reported with the best least-squares solution.

### Two-Phase Restoration

When the filter line search fails:

1. **Phase 1 (Gauss-Newton)**: Fast feasibility solver minimizing ||violation||^2 with gradient descent fallback
2. **Phase 2 (NLP restoration)**: Full barrier subproblem with p/n slack variables (Ipopt formulation)
3. **Multi-attempt recovery**: Up to 6 attempts cycling barrier parameter perturbations [10x, 0.1x, 100x, 0.01x, 1000x, 0.001x] with x perturbation

### Watchdog Strategy

The solver implements a watchdog mechanism that temporarily relaxes the filter acceptance criteria when progress stalls due to shortened steps. This helps escape narrow feasible corridors where strict Armijo conditions are too conservative.

## Profiling

### Built-in Phase Timing

ripopt includes per-iteration phase timing instrumentation. When `print_level >= 5` (the default), a summary table is printed at the end of each solve showing where CPU time is spent:

```
Phase breakdown (47 iterations):
  Problem eval           0.234s (45.2%)
  KKT assembly           0.089s (17.2%)
  Factorization          0.156s (30.1%)
  Direction solve        0.012s  (2.3%)
  Line search            0.021s  (4.1%)
  Other                  0.006s  (1.1%)
  Total                  0.518s
```

To suppress timing output, set `print_level: 0` in `SolverOptions`.

### External Profiling with samply

Release builds include debug symbols (`debug = true` in `[profile.release]`), so external profilers can show function names. [samply](https://github.com/mstange/samply) provides flamegraph visualization on macOS and Linux:

```bash
cargo install samply
cargo build --release --bin hs_suite
samply record target/release/hs_suite
```

This opens a Firefox Profiler UI in the browser with a full call tree and flamegraph. Look for wide bars under `solve_ipm` to identify dominant functions.

On macOS, Instruments (Xcode) also works without any additional setup:

```bash
cargo build --release --bin hs_suite
xcrun xctrace record --template "Time Profiler" --launch target/release/hs_suite
```

## Sign Convention

ripopt uses the Ipopt convention where the Lagrangian is:

```
L = f(x) + y^T g(x)
```

For inequality constraints `g(x) >= g_l`, the multiplier `y` is negative at optimality. For equality constraints, `y` can be positive or negative.

## License

[EPL-2.0](LICENSE) (Eclipse Public License 2.0), consistent with Ipopt.
