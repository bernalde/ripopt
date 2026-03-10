# ripopt Solver Diagnostics Guide

## Overview

`SolverDiagnostics` captures structured data about solver behavior during a
solve. It is available both programmatically via `result.diagnostics` and as a
stderr summary block (printed when `print_level >= 5`).

The diagnostic block looks like:

```
--- ripopt diagnostics ---
status: Optimal
iterations: 8
wall_time: 0.001s
final_mu: 4.14e-9
final_primal_inf: 5.55e-17
final_dual_inf: 1.04e-12
final_compl: 7.75e-9
restoration_count: 0
nlp_restoration_count: 0
mu_mode_switches: 2
filter_rejects: 0
watchdog_activations: 0
soc_corrections: 0
--- end diagnostics ---
```

## Fields

| Field | Meaning |
|---|---|
| `restoration_count` | GN (Gauss-Newton) restoration entries |
| `nlp_restoration_count` | Full NLP restoration entries (heavier, Ipopt-style) |
| `mu_mode_switches` | Barrier mode transitions (Free <-> Fixed) |
| `filter_rejects` | Line search failures (backtracking exhausted) |
| `watchdog_activations` | Watchdog triggered by consecutive short steps |
| `soc_corrections` | Second-order corrections accepted |
| `final_mu` | Barrier parameter at termination |
| `final_primal_inf` | Constraint violation at termination |
| `final_dual_inf` | Dual infeasibility (stationarity error) |
| `final_compl` | Complementarity error at termination |
| `wall_time_secs` | Total wall-clock time |
| `fallback_used` | Which fallback succeeded, if any (`lbfgs_hessian`, `augmented_lagrangian`, `sqp`, `slack`) |

## Reading the diagnostics

**Healthy solve** (HS071-like): 0 restorations, 0 filter rejects, 2-4 mu mode
switches, `final_mu` near `1e-9`, `final_primal_inf` and `final_dual_inf` both
below `tol`.

**Struggling solve**: Many filter rejects, multiple restorations, `final_mu`
stuck above `1e-4`, or a fallback was used.

**Key patterns and what to try:**

| Pattern | Likely cause | Options to adjust |
|---|---|---|
| `filter_rejects` > 5 | Line search fighting constraints | Increase `mu_init`, reduce `kappa` |
| `restoration_count` > 3 | Repeated feasibility recovery | Try `enable_slack_fallback`, or increase `mu_init` |
| `mu_mode_switches` > 10 | Free/Fixed cycling | Set `mu_strategy_adaptive: false` for monotone mode |
| `final_mu` stuck > 1e-4 | Barrier parameter not decreasing | Increase `max_iter`, reduce `mu_linear_decrease_factor` |
| `fallback_used: Some(...)` | Primary IPM failed | Examine which fallback; consider changing Hessian strategy |
| `soc_corrections` > 0 | Nonlinear constraints causing step rejection | Normal; increase `max_soc` if filter rejects are also high |
| `watchdog_activations` > 0 | Tiny steps detected | Try `hessian_approximation_lbfgs: true` |

---

## Example 1: Easy problem (HS071)

HS071 is a 4-variable, 2-constraint nonlinear program. ripopt solves it in ~8
iterations with default options.

```rust
use ripopt::{NlpProblem, SolveStatus, SolverOptions};

// min  x1*x4*(x1+x2+x3) + x3
// s.t. x1*x2*x3*x4 >= 25
//      x1^2+x2^2+x3^2+x4^2 = 40
//      1 <= xi <= 5
// x0 = (1, 5, 5, 1),  f* ~ 17.014

struct Hs071;

impl NlpProblem for Hs071 {
    fn num_variables(&self) -> usize { 4 }
    fn num_constraints(&self) -> usize { 2 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..4 { x_l[i] = 1.0; x_u[i] = 5.0; }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l[0] = 25.0; g_u[0] = f64::INFINITY;  // product >= 25
        g_l[1] = 40.0; g_u[1] = 40.0;            // sum of squares = 40
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 1.0; x0[1] = 5.0; x0[2] = 5.0; x0[3] = 1.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = x[3] * (2.0 * x[0] + x[1] + x[2]);
        grad[1] = x[0] * x[3];
        grad[2] = x[0] * x[3] + 1.0;
        grad[3] = x[0] * (x[0] + x[1] + x[2]);
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] * x[1] * x[2] * x[3];
        g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0,0,0,0, 1,1,1,1], vec![0,1,2,3, 0,1,2,3])
    }

    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        vals[0] = x[1]*x[2]*x[3]; vals[1] = x[0]*x[2]*x[3];
        vals[2] = x[0]*x[1]*x[3]; vals[3] = x[0]*x[1]*x[2];
        vals[4] = 2.0*x[0]; vals[5] = 2.0*x[1];
        vals[6] = 2.0*x[2]; vals[7] = 2.0*x[3];
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0,1,1,2,2,2,3,3,3,3], vec![0,0,1,0,1,2,0,1,2,3])
    }

    fn hessian_values(&self, x: &[f64], s: f64, l: &[f64], v: &mut [f64]) {
        v[0] = s*2.0*x[3] + l[1]*2.0;
        v[1] = s*x[3] + l[0]*x[2]*x[3];
        v[2] = l[1]*2.0;
        v[3] = s*x[3] + l[0]*x[1]*x[3];
        v[4] = l[0]*x[0]*x[3];
        v[5] = l[1]*2.0;
        v[6] = s*(2.0*x[0]+x[1]+x[2]) + l[0]*x[1]*x[2];
        v[7] = s*x[0] + l[0]*x[0]*x[2];
        v[8] = s*x[0] + l[0]*x[0]*x[1];
        v[9] = l[1]*2.0;
    }
}

fn main() {
    let result = ripopt::solve(&Hs071, &SolverOptions::default());

    assert_eq!(result.status, SolveStatus::Optimal);
    assert!((result.objective - 17.014).abs() < 0.01);

    // Diagnostics: expect clean solve
    let d = &result.diagnostics;
    println!("iterations: {}", result.iterations);   // ~8
    println!("filter_rejects: {}", d.filter_rejects); // 0
    println!("restorations: {}", d.restoration_count); // 0
}
```

**Expected diagnostics:**
```
status: Optimal
iterations: 8
filter_rejects: 0
restoration_count: 0
mu_mode_switches: 2
final_mu: ~4e-9
```

No drama. The solver converges in a straight line.

---

## Example 2: Hard problem (TP374 — unsolved)

TP374 has 10 variables, 35 nonlinear inequality constraints involving
trigonometric sums. ripopt hits `MaxIterations` at 2999 iterations.
Known optimal: f* = 0.233264.

```rust
use ripopt::{NlpProblem, SolverOptions};
use std::f64::consts::PI;

struct TP374;

fn tp374_a(z: f64, x: &[f64]) -> f64 {
    (1..=9).map(|k| x[k-1] * (k as f64 * z).cos()).sum()
}
fn tp374_b(z: f64, x: &[f64]) -> f64 {
    (1..=9).map(|k| x[k-1] * (k as f64 * z).sin()).sum()
}
fn tp374_g(z: f64, x: &[f64]) -> f64 {
    let (a, b) = (tp374_a(z, x), tp374_b(z, x));
    a*a + b*b
}

impl NlpProblem for TP374 {
    fn num_variables(&self) -> usize { 10 }
    fn num_constraints(&self) -> usize { 35 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..10 { x_l[i] = f64::NEG_INFINITY; x_u[i] = f64::INFINITY; }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        for i in 0..35 { g_l[i] = 0.0; g_u[i] = f64::INFINITY; }
    }

    fn initial_point(&self, x0: &mut [f64]) {
        for i in 0..10 { x0[i] = 0.1; }
    }

    fn objective(&self, x: &[f64]) -> f64 { x[9] }

    fn gradient(&self, _x: &[f64], grad: &mut [f64]) {
        for i in 0..9 { grad[i] = 0.0; }
        grad[9] = 1.0;
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        for i in 0..10 {
            let z = PI/4.0 * (i as f64 * 0.1);
            g[i] = tp374_g(z, x) - (1.0 - x[9]).powi(2);
        }
        for i in 10..20 {
            let z = PI/4.0 * ((i-10) as f64 * 0.1);
            g[i] = (1.0 + x[9]).powi(2) - tp374_g(z, x);
        }
        for i in 20..35 {
            let z = PI/4.0 * (1.2 + (i-20) as f64 * 0.2);
            g[i] = x[9].powi(2) - tp374_g(z, x);
        }
    }

    // ... jacobian and hessian implementations omitted for brevity
    // (see examples/debug_tp374.rs for the full implementation)
#   fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
#       let mut rows = Vec::new(); let mut cols = Vec::new();
#       for i in 0..35 { for j in 0..10 { rows.push(i); cols.push(j); } }
#       (rows, cols)
#   }
#   fn jacobian_values(&self, _x: &[f64], _v: &mut [f64]) { /* ... */ }
#   fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
#       let mut rows = Vec::new(); let mut cols = Vec::new();
#       for i in 0..10 { for j in 0..=i { rows.push(i); cols.push(j); } }
#       (rows, cols)
#   }
#   fn hessian_values(&self, _x: &[f64], _s: f64, _l: &[f64], _v: &mut [f64]) { /* ... */ }
}

fn main() {
    // Attempt 1: default options
    let opts1 = SolverOptions {
        print_level: 5,
        max_iter: 3000,
        ..SolverOptions::default()
    };
    let r1 = ripopt::solve(&TP374, &opts1);
    println!("Attempt 1: {:?}, obj={:.6}, iters={}", r1.status, r1.objective, r1.iterations);
    println!("  filter_rejects={}, restorations={}, mu_switches={}",
        r1.diagnostics.filter_rejects,
        r1.diagnostics.restoration_count,
        r1.diagnostics.mu_mode_switches);

    // Attempt 2: adjust based on diagnostics
    // If filter_rejects is high -> raise mu_init, lower kappa
    // If mu stuck -> try monotone mode
    let opts2 = SolverOptions {
        print_level: 5,
        max_iter: 3000,
        mu_init: 1.0,
        kappa: 3.0,
        hessian_approximation_lbfgs: true,
        ..SolverOptions::default()
    };
    let r2 = ripopt::solve(&TP374, &opts2);
    println!("Attempt 2: {:?}, obj={:.6}, iters={}", r2.status, r2.objective, r2.iterations);
}
```

**Why it's hard:** The 35 trigonometric inequality constraints create a
narrow feasible region with many near-degenerate active sets. The solver
accumulates filter rejects and restorations, and the barrier parameter
fails to decrease toward zero. Constraint multipliers grow to ~1e11,
signaling numerical instability.

**What diagnostics tell you:**
```
status: MaxIterations
iterations: 2999
filter_rejects: ~50+
restoration_count: ~10+
mu_mode_switches: ~20+
final_mu: ~1e-2 (stuck, should be ~1e-9)
final_primal_inf: ~1e-1 (not feasible)
```

**Strategies to try (informed by diagnostics):**

1. **High filter_rejects** -> Increase `mu_init` to 1.0 or 10.0, giving the
   solver more room for infeasible steps early on
2. **Restorations dominating** -> Enable `enable_slack_fallback: true` to
   reformulate inequalities with explicit slacks
3. **mu stuck** -> Try `mu_strategy_adaptive: false` for monotone decrease,
   or reduce `kappa` to slow mu reduction
4. **Large multipliers** -> Try `hessian_approximation_lbfgs: true` to avoid
   ill-conditioned exact Hessians
5. **Still stuck** -> Try a different starting point; the initial `x0 = 0.1`
   may be in a bad basin

---

## LLM-Guided Steering Workflow

The intended workflow for Claude Code or similar agents:

```
1. Call ripopt::solve(&problem, &default_options)
2. Read result.diagnostics
3. IF status == Optimal: done
4. ELSE: reason about diagnostics pattern
   - "15 filter_rejects + 3 restorations + mu stuck at 1e-3"
   - -> "Line search is fighting the filter. Try higher mu_init
        and slack reformulation."
5. Adjust SolverOptions accordingly
6. Re-solve and compare
7. Repeat (with memory of what was tried)
```

The key SolverOptions knobs for steering:

| Option | Default | Effect |
|---|---|---|
| `mu_init` | 0.1 | Higher = more room for infeasible exploration |
| `kappa` | 10.0 | Lower = slower mu decrease (more conservative) |
| `mu_strategy_adaptive` | true | false = monotone decrease (simpler, sometimes better) |
| `mu_linear_decrease_factor` | 0.2 | Controls monotone mu reduction speed |
| `max_iter` | 3000 | Budget |
| `max_soc` | 4 | More SOC = better for nonlinear constraints |
| `hessian_approximation_lbfgs` | false | true = skip exact Hessian (robustness) |
| `enable_slack_fallback` | true | Explicit slacks for inequalities |
| `enable_al_fallback` | true | Augmented Lagrangian for equalities |
| `enable_sqp_fallback` | true | SQP for small constrained problems |
| `mehrotra_pc` | false | Predictor-corrector (fewer iterations) |
| `gondzio_mcc_max` | 0 | Centrality corrections (better centering) |
| `warm_start` | false | Reuse previous solution as starting point |
