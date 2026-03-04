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

It implements a primal-dual interior point method with a barrier formulation, similar to the algorithm described in the Ipopt papers. The solver is written entirely in Rust (~14,200 lines) with no external C/Fortran dependencies.

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
- **Preprocessing**: Automatic elimination of fixed variables, redundant constraints, and bound tightening from single-variable linear constraints
- **Near-linear constraint detection**: Automatically identifies linear constraints and skips their Hessian contribution
- **Multi-solver fallback architecture**: L-BFGS, Augmented Lagrangian, SQP, and explicit slack reformulation
- **C API** mirroring the Ipopt C interface for direct linking from C/C++/Python/Julia
- **AMPL NL interface** with Pyomo integration via `SolverFactory('ripopt')`

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
| Total solved  | **598/727 (82.3%)** | 555/727 (76.3%)        |
| Constrained   | **376/493**         | 339/493                |
| Unconstrained | **222/234**         | 216/234                |
| Both solve    | 555                 | 555                    |
| ripopt only   | **43**              | --                     |
| Ipopt only    | --                  | **0**                  |
| Both fail     | 129                 | 129                    |

ripopt solves **43 more problems** than Ipopt, and solves everything Ipopt solves. (Counts vary ±2 between runs due to timing sensitivity on borderline problems.)

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

**Interpreting the speed numbers.** The raw median across all 555 mutually-solved problems is 30x, but this is misleading. Most CUTEst problems are small (n < 10) and solve in microseconds for ripopt, while Ipopt has a ~1-3ms floor from internal initialization (MA27 symbolic analysis, TNLP setup) that dominates for trivial problems. The 3.3x median on non-trivial problems is a more meaningful comparison.

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
| `enable_preprocessing` | true | Eliminate fixed variables and redundant constraints |
| `detect_linear_constraints` | true | Skip Hessian for linear constraints |
| `enable_sqp_fallback` | true | SQP fallback for constrained problems |

### Result

`SolveResult` contains:

- `x` -- optimal primal variables
- `objective` -- optimal objective value f(x*)
- `constraint_multipliers` -- Lagrange multipliers for constraints (y)
- `bound_multipliers_lower` / `bound_multipliers_upper` -- bound multipliers (z_L, z_U)
- `constraint_values` -- constraint values g(x*)
- `status` -- one of: `Optimal`, `Acceptable`, `Infeasible`, `LocalInfeasibility`, `MaxIterations`, `NumericalError`, `Unbounded`, `RestorationFailed`, `InternalError`
- `iterations` -- number of IPM iterations

## C API

ripopt exposes a C API that mirrors the [Ipopt C interface](https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_C), enabling direct linking from C, C++, Python (`ctypes`/`cffi`), Julia, and any language with C FFI support — without the subprocess/file overhead of the NL interface. If you have existing Ipopt C code, migrating to ripopt requires only header/function renaming; the callback signatures are identical.

### Build the shared library

```bash
cargo build --release
# produces target/release/libripopt.dylib (macOS) or libripopt.so (Linux)
```

### C header

Include `ripopt.h` (repo root) in your C project. It defines version macros, callback typedefs, return status codes, and all public functions:

```c
#include "ripopt.h"

// Check version at compile time
printf("ripopt %s\n", RIPOPT_VERSION);  // "0.2.0"
```

### Callback signatures

The five callback types are identical to the Ipopt C interface. All callbacks return `1` on success, `0` on error (the solver will abort if a callback returns `0`):

```c
typedef int (*Eval_F_CB)   (int n, const double *x, int new_x,
                             double *obj_value, void *user_data);
typedef int (*Eval_Grad_F_CB)(int n, const double *x, int new_x,
                              double *grad_f, void *user_data);
typedef int (*Eval_G_CB)   (int n, const double *x, int new_x,
                             int m, double *g, void *user_data);
typedef int (*Eval_Jac_G_CB)(int n, const double *x, int new_x,
                              int m, int nele_jac,
                              int *iRow, int *jCol, double *values,
                              void *user_data);
typedef int (*Eval_H_CB)   (int n, const double *x, int new_x,
                             double obj_factor,
                             int m, const double *lambda, int new_lambda,
                             int nele_hess,
                             int *iRow, int *jCol, double *values,
                             void *user_data);
```

**Two-call protocol for Jacobian and Hessian:** When `values == NULL`, fill `iRow`/`jCol` with the sparsity pattern (0-based indexing); when `values != NULL`, fill numerical values in the same element order as the pattern. The Hessian uses the **lower triangle** only.

**Sign convention:** ripopt uses the Ipopt convention L = f(x) + y^T g(x). The Hessian callback receives `obj_factor` and `lambda` and should compute `obj_factor * ∇²f + Σ lambda[i] * ∇²g_i`.

### Lifecycle

```c
#include "ripopt.h"

// 1. Create handle
RipoptProblem nlp = ripopt_create(
    n, x_l, x_u,          // variable bounds (use ±1e30 for ±∞)
    m, g_l, g_u,           // constraint bounds (g_l == g_u for equality)
    nele_jac, nele_hess,   // number of nonzeros
    eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h);

// 2. Set options (Ipopt-compatible key names)
ripopt_add_int_option(nlp, "print_level", 5);
ripopt_add_num_option(nlp, "tol",         1e-8);
ripopt_add_str_option(nlp, "mu_strategy", "adaptive");

// 3. Solve  (x: in = initial point, out = solution)
double obj_val;
int status = ripopt_solve(nlp, x, NULL, &obj_val,
                          NULL, NULL, NULL, NULL);
// status == 0  →  RIPOPT_SOLVE_SUCCEEDED

// 4. Free
ripopt_free(nlp);
```

For unconstrained problems, pass `m=0` and `NULL` for `g_l`/`g_u`.

**Infinity bounds:** Use `HUGE_VAL` (from `<math.h>`) for "no bound". Internally, any value beyond `±1e19` is treated as unbounded. Avoid using finite large values like `1e30` — they may cause numerical issues.

### Extracting multipliers

All output pointers except `x` are optional (pass `NULL` to skip). Here is how to extract the full solution including Lagrange multipliers and bound multipliers:

```c
double x[4]      = {1.0, 5.0, 5.0, 1.0};  // initial point
double obj_val   = 0.0;
double g[2]      = {0.0, 0.0};              // constraint values at solution
double mult_g[2] = {0.0, 0.0};              // constraint multipliers (lambda)
double mult_xl[4]= {0.0, 0.0, 0.0, 0.0};   // lower bound multipliers (z_L)
double mult_xu[4]= {0.0, 0.0, 0.0, 0.0};   // upper bound multipliers (z_U)

int status = ripopt_solve(nlp, x, g, &obj_val,
                          mult_g, mult_xl, mult_xu,
                          NULL);  // user_data

// At the solution:
// - x[]       contains the optimal primal variables
// - obj_val   is f(x*)
// - g[]       contains g(x*) — verify constraints are satisfied
// - mult_g[]  contains the Lagrange multipliers for constraints
//             (nonzero for active constraints)
// - mult_xl[] contains z_L (positive when x is at its lower bound)
// - mult_xu[] contains z_U (positive when x is at its upper bound)
```

The `user_data` pointer is forwarded to every callback unchanged — use it to pass problem-specific data (e.g., model parameters) without globals.

### Return status

| Code | Enum constant | Meaning |
|------|---------------|---------|
| 0  | `RIPOPT_SOLVE_SUCCEEDED` | Converged to optimal solution |
| 1  | `RIPOPT_ACCEPTABLE_LEVEL` | Converged to acceptable (less strict) tolerance |
| 2  | `RIPOPT_INFEASIBLE_PROBLEM` | Problem is locally infeasible |
| 5  | `RIPOPT_MAXITER_EXCEEDED` | Reached iteration limit |
| 6  | `RIPOPT_RESTORATION_FAILED` | Feasibility restoration failed |
| 7  | `RIPOPT_ERROR_IN_STEP_COMPUTATION` | Numerical difficulties |
| 10 | `RIPOPT_NOT_ENOUGH_DEGREES_OF_FREEDOM` | Problem has too few free variables |
| 11 | `RIPOPT_INVALID_PROBLEM_DEFINITION` | Problem appears unbounded |
| -1 | `RIPOPT_INTERNAL_ERROR` | Internal error |

Status 0 and 1 indicate a successful solve. All others indicate failure — check your problem formulation, initial point, or try adjusting options.

### Options reference

Option-setting functions return `1` on success, `0` if the keyword is unknown. All option keywords match Ipopt naming conventions.

**Numeric options** (`ripopt_add_num_option`):

| Option | Default | Description |
|--------|---------|-------------|
| `tol` | 1e-8 | Convergence tolerance |
| `acceptable_tol` | 1e-4 | Acceptable convergence tolerance |
| `acceptable_constr_viol_tol` | 1e-2 | Acceptable constraint violation |
| `acceptable_dual_inf_tol` | 1e10 | Acceptable dual infeasibility |
| `acceptable_compl_inf_tol` | 1e-2 | Acceptable complementarity |
| `mu_init` | 0.1 | Initial barrier parameter |
| `mu_min` | 1e-11 | Minimum barrier parameter |
| `bound_push` | 1e-2 | Initial bound push |
| `bound_frac` | 1e-2 | Initial bound fraction |
| `constr_viol_tol` | 1e-4 | Constraint violation tolerance |
| `dual_inf_tol` | 100.0 | Dual infeasibility tolerance |
| `compl_inf_tol` | 1e-4 | Complementarity tolerance |
| `max_wall_time` | 0.0 | Wall-clock time limit in seconds (0 = no limit) |
| `warm_start_bound_push` | 1e-3 | Warm-start bound push |
| `warm_start_bound_frac` | 1e-3 | Warm-start bound fraction |
| `warm_start_mult_bound_push` | 1e-3 | Warm-start multiplier push |
| `nlp_lower_bound_inf` | -1e19 | Threshold for -infinity bounds |
| `nlp_upper_bound_inf` | 1e19 | Threshold for +infinity bounds |
| `kappa` | 10.0 | Adaptive mu divisor |
| `constr_mult_init_max` | 1000.0 | Max initial constraint multiplier |
| `barrier_tol_factor` | 10.0 | Barrier tolerance factor |

**Integer options** (`ripopt_add_int_option`):

| Option | Default | Description |
|--------|---------|-------------|
| `max_iter` | 3000 | Maximum iterations |
| `print_level` | 5 | Output verbosity (0 = silent, 5 = verbose, 12 = debug) |
| `acceptable_iter` | 10 | Consecutive acceptable iterations for convergence |
| `max_soc` | 4 | Maximum second-order correction steps |
| `sparse_threshold` | 110 | KKT dimension threshold for sparse solver |
| `restoration_max_iter` | 200 | Max iterations in NLP restoration subproblem |

**String options** (`ripopt_add_str_option`):

| Option | Default | Values | Description |
|--------|---------|--------|-------------|
| `mu_strategy` | `"adaptive"` | `"adaptive"`, `"monotone"` | Barrier parameter update strategy |
| `warm_start_init_point` | `"no"` | `"yes"`, `"no"` | Enable warm-start initialization |
| `mu_allow_increase` | `"yes"` | `"yes"`, `"no"` | Allow barrier parameter increase |
| `least_squares_mult_init` | `"yes"` | `"yes"`, `"no"` | LS estimate for initial multipliers |
| `enable_slack_fallback` | `"yes"` | `"yes"`, `"no"` | Slack reformulation fallback |
| `enable_lbfgs_fallback` | `"yes"` | `"yes"`, `"no"` | L-BFGS fallback for unconstrained |
| `enable_al_fallback` | `"yes"` | `"yes"`, `"no"` | Augmented Lagrangian fallback |
| `enable_preprocessing` | `"yes"` | `"yes"`, `"no"` | Preprocessing (fixed vars, redundant constraints) |
| `detect_linear_constraints` | `"yes"` | `"yes"`, `"no"` | Skip Hessian for linear constraints |
| `enable_sqp_fallback` | `"yes"` | `"yes"`, `"no"` | SQP fallback for constrained problems |

### Error handling

- **Callback errors:** If any callback returns `0`, the solver aborts and returns `RIPOPT_ERROR_IN_STEP_COMPUTATION` (7). Always return `1` from callbacks unless you detect a problem (e.g., NaN in inputs).
- **Unknown options:** `ripopt_add_*_option` returns `0` for unrecognized keywords. Check the return value if you want to detect typos.
- **NULL safety:** `ripopt_free(NULL)` is a no-op (safe to call). All output pointers in `ripopt_solve` except `x` may be `NULL`.
- **Memory:** The problem handle owns all internal memory. Call `ripopt_free()` once when done. Do not use the handle after freeing.

### Migrating from Ipopt

If you have existing Ipopt C code, the migration is straightforward:

1. **Header:** `#include "IpStdCInterface.h"` → `#include "ripopt.h"`
2. **Handle type:** `IpoptProblem` → `RipoptProblem` (both are `void*`)
3. **Functions:** Rename `CreateIpoptProblem` → `ripopt_create`, `FreeIpoptProblem` → `ripopt_free`, `AddIpoptNumOption` → `ripopt_add_num_option`, etc.
4. **Callbacks:** No changes required — signatures are identical
5. **Status codes:** Similar semantics but different enum names (e.g., `Solve_Succeeded` → `RIPOPT_SOLVE_SUCCEEDED`)
6. **Infinity:** Ipopt uses ±2e19 by default; ripopt uses ±1e30 in bounds and ±1e19 for `nlp_*_bound_inf`
7. **Linking:** `-lipopt` → `-lripopt`

### Compile and run the examples

```bash
cargo build --release

# HS071 — constrained NLP with inequality + equality constraints
cc examples/c_api_test.c -I. -Ltarget/release -lripopt \
   -Wl,-rpath,$(pwd)/target/release -o c_api_test -lm
./c_api_test

# Rosenbrock — unconstrained optimization
cc examples/c_rosenbrock.c -I. -Ltarget/release -lripopt \
   -Wl,-rpath,$(pwd)/target/release -o c_rosenbrock -lm
./c_rosenbrock

# HS035 — bound-constrained QP with inequality
cc examples/c_hs035.c -I. -Ltarget/release -lripopt \
   -Wl,-rpath,$(pwd)/target/release -o c_hs035 -lm
./c_hs035

# Full multiplier extraction and options demonstration
cc examples/c_example_with_options.c -I. -Ltarget/release -lripopt \
   -Wl,-rpath,$(pwd)/target/release -o c_example_with_options -lm
./c_example_with_options
```

## Examples

### Rust

```bash
# Rosenbrock function (unconstrained with bounds)
cargo run --example rosenbrock

# HS071 (constrained NLP with inequalities)
cargo run --example hs071

# Benchmark timing across 5 problems
cargo run --release --example benchmark
```

### C

See [Compile and run the examples](#compile-and-run-the-examples) above for build instructions. The C examples are:

| Example | Problem type | Demonstrates |
|---------|-------------|-------------|
| `c_api_test.c` | HS071 (constrained) | Basic usage, all 5 callbacks |
| `c_rosenbrock.c` | Rosenbrock (unconstrained) | No constraints, no bounds |
| `c_hs035.c` | HS035 (bounds + inequality) | Bound multipliers, constraint multiplier |
| `c_example_with_options.c` | HS071 (multiple solves) | Options tuning, multiplier extraction, status interpretation |

## Tests

```bash
cargo test
```

162 tests total:
- **114 unit tests**: Dense LDL factorization, convergence checking, filter line search, fraction-to-boundary, KKT assembly, restoration, preprocessing, linearity detection, SQP, linear solver
- **12 C API tests**: FFI integration tests
- **21 integration tests**: Rosenbrock, SimpleQP, HS071, HS035, PureBoundConstrained, MultipleEqualityConstraints, NE-to-LS reformulation, and more
- **15 HS regression tests**: Selected Hock-Schittkowski problems for regression checking

## Architecture

```
src/
  lib.rs              Public API (solve function, re-exports)
  c_api.rs            C FFI layer (extern "C" functions, ripopt.h)
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
  sqp.rs              SQP fallback for constrained problems
  slack_formulation.rs     Explicit slack reformulation fallback
  preprocessing.rs    Fixed variable elimination, redundant constraint removal, bound tightening
  linearity.rs        Near-linear constraint detection
  warmstart.rs        Warm-start initialization
  linear_solver/
    mod.rs            LinearSolver trait, SymmetricMatrix
    dense.rs          Dense LDL^T (Bunch-Kaufman) factorization

tests/
  correctness.rs      Integration tests (21 NLP problems)
  hs_regression.rs    HS suite regression tests (15 problems)
  c_api.rs            C API integration tests (11 tests via FFI)

examples/
  rosenbrock.rs       Unconstrained optimization
  hs071.rs            Constrained NLP
  benchmark.rs        Timing benchmark
  c_api_test.c        HS071 via the C API
  c_rosenbrock.c      Unconstrained Rosenbrock via C API
  c_hs035.c           Bound-constrained QP via C API
  c_example_with_options.c  Options and multiplier extraction demo
```

## Algorithm Details

### Preprocessing

Before solving, ripopt automatically analyzes the problem to reduce its size:

1. **Fixed variable elimination**: Variables with `x_l == x_u` are removed and set to their fixed values in all evaluations
2. **Redundant constraint removal**: Duplicate constraints (same Jacobian structure, values, and bounds) are eliminated
3. **Bound tightening**: Single-variable linear constraints are used to tighten variable bounds

The reduced problem is solved and the solution is mapped back to the original dimensions. Disable with `enable_preprocessing: false`.

### Near-Linear Constraint Detection

The Jacobian is evaluated at two points to identify linear constraints (where all Jacobian entries remain constant). For linear constraints, the Hessian contribution `lambda[i] * nabla^2 g_i` is exactly zero and is skipped, reducing computation in the Hessian evaluation.

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
3. **SQP**: Equality-constrained Sequential Quadratic Programming with l1 merit function line search
4. **Explicit slack reformulation**: Converts g(x) to g(x)-s=0 with bounds on s, stabilizing multiplier oscillation at degenerate points
5. **Best-du tracking**: Throughout the solve, tracks the iterate with lowest dual infeasibility and recovers it at max iterations

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
