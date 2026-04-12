# Hock-Schittkowski Suite

The Hock-Schittkowski (HS) test set is a collection of 120 small, classic
nonlinear programming problems from Hock, W. and Schittkowski, K., *Test
Examples for Nonlinear Programming Codes*, Lecture Notes in Economics and
Mathematical Systems, Springer, 1981. These problems are the standard
regression suite for NLP solvers: they exercise every common structure
(bounds, equalities, inequalities, degenerate multipliers, non-convex
objectives, nonlinear constraints) while staying small enough to debug.

## Contents

- `generated/hs_problems.rs` — ripopt-compatible problem definitions (generated
  by `generate.py` from the original Fortran sources in `PROB.FOR`)
- `run_ripopt.rs` — the `hs_suite` binary (ripopt runner)
- `run_ipopt_native.rs` — the `ipopt_native` binary (Ipopt C API runner)
- `run_cyipopt.py` / `generated/hs_cyipopt.py` — Python cyipopt runner for
  cross-validation
- `compare.py` — side-by-side report generator
- `PERFORMANCE_COMPARISON.md` — full HS-specific report with per-problem
  timing and iteration counts
- `HS_VALIDATION_REPORT.md` — correctness check against known optima
- `hs_ripopt_results.json`, `hs_ipopt_native_results.json`,
  `hs_cyipopt_results.json` — per-solver results for the latest run
- `hs_ripopt_results_vX.Y.Z.json` — tagged historical results for regression
  tracking

## Prerequisites

None beyond ripopt itself (plus the `ipopt-native` feature if you want the
side-by-side comparison).

## How to run

From the repo root:

```bash
make hs-run
```

or equivalently the underlying cargo commands:

```bash
cargo run --release --bin hs_suite --features hs,ipopt-native \
    > benchmarks/hs/hs_ripopt_results.json
cargo run --release --bin ipopt_native --features ipopt-native \
    > benchmarks/hs/hs_ipopt_native_results.json
```

The runs take roughly two minutes each.

## Output

- `hs_ripopt_results.json` — ripopt results
- `hs_ipopt_native_results.json` — Ipopt C API results
- `hs_ripopt_stderr.txt`, `hs_ipopt_native_stderr.txt` — solver chatter

The HS suite feeds the composite `benchmarks/BENCHMARK_REPORT.md` and is
also summarised in detail in `PERFORMANCE_COMPARISON.md`.
