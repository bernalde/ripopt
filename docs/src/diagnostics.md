# Diagnostics

`SolverDiagnostics` captures structured data about solver behavior. It is available at `result.diagnostics` and printed to stderr when `print_level >= 5`.

## Diagnostic block format

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
| `status` | Final solve status |
| `iterations` | Total IPM iterations |
| `wall_time_secs` | Total wall-clock time |
| `final_mu` | Barrier parameter at termination |
| `final_primal_inf` | Constraint violation at termination |
| `final_dual_inf` | Dual infeasibility (stationarity error) at termination |
| `final_compl` | Complementarity error at termination |
| `restoration_count` | Gauss-Newton restoration entries |
| `nlp_restoration_count` | Full NLP restoration entries (heavier) |
| `mu_mode_switches` | Barrier mode transitions (Free ↔ Fixed) |
| `filter_rejects` | Line search failures (backtracking exhausted) |
| `watchdog_activations` | Watchdog triggered by consecutive short steps |
| `soc_corrections` | Second-order corrections accepted |
| `fallback_used` | Which fallback succeeded, if any (`lbfgs_hessian`, `augmented_lagrangian`, `sqp`, `slack`) |

## Interpreting the diagnostics

**Healthy solve** (HS071-like): 0 restorations, 0 filter rejects, 2–4 mu mode switches, `final_mu` near `1e-9`, `final_primal_inf` and `final_dual_inf` both below `tol`.

**Struggling solve**: Many filter rejects, multiple restorations, `final_mu` stuck above `1e-4`, or a fallback was used.

### Pattern → cause → fix

| Pattern | Likely cause | Options to try |
|---|---|---|
| `filter_rejects` > 5 | Line search fighting constraints | Increase `mu_init`, reduce `kappa` |
| `restoration_count` > 3 | Repeated feasibility recovery | Set `enable_slack_fallback: true`, increase `mu_init` |
| `mu_mode_switches` > 10 | Free/Fixed cycling | Set `mu_strategy_adaptive: false` |
| `final_mu` stuck > 1e-4 | Barrier parameter not decreasing | Increase `max_iter`, reduce `mu_linear_decrease_factor` |
| `fallback_used: Some(...)` | Primary IPM failed | Check which fallback; consider changing Hessian strategy |
| `soc_corrections` > 0 | Nonlinear constraints causing step rejection | Normal; increase `max_soc` if filter rejects are also high |
| `watchdog_activations` > 0 | Tiny steps detected | Try `hessian_approximation_lbfgs: true` |

## Example: healthy solve (HS071)

```rust
let result = ripopt::solve(&Hs071, &SolverOptions::default());
let d = &result.diagnostics;
assert_eq!(d.filter_rejects, 0);
assert_eq!(d.restoration_count, 0);
assert!(d.final_mu < 1e-8);
// iterations ≈ 8
```

## Example: struggling solve — reading and reacting

```rust
let r1 = ripopt::solve(&problem, &opts);

if r1.status != SolveStatus::Optimal {
    let d = &r1.diagnostics;
    let mut opts2 = opts.clone();

    if d.filter_rejects > 5 {
        opts2.mu_init = 1.0;
        opts2.kappa = 3.0;
    }
    if d.restoration_count > 3 {
        opts2.enable_slack_fallback = true;
    }
    if d.final_mu > 1e-4 {
        opts2.mu_strategy_adaptive = false;
    }
    if d.mu_mode_switches > 10 {
        opts2.hessian_approximation_lbfgs = true;
    }

    let r2 = ripopt::solve(&problem, &opts2);
}
```

## Using diagnostics with Claude Code

Claude Code can read the `--- ripopt diagnostics ---` block from stderr and automatically adjust solver options:

```bash
claude -p "
  Run: cargo run --example debug_tp374 2>&1
  Parse the diagnostics block.
  If not Optimal:
    - High filter_rejects → increase mu_init, decrease kappa
    - High restoration_count → try enable_slack_fallback
    - mu stuck high → try mu_strategy_adaptive: false
    - Large multipliers → try hessian_approximation_lbfgs: true
  Adjust options, re-run, compare. Up to 3 attempts.
"
```

The Rust code is a pure reporter. All intelligence — pattern matching, strategy selection — lives in Claude's reasoning.
