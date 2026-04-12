# Electrical Grid Suite (AC Optimal Power Flow)

AC Optimal Power Flow (AC OPF) is the canonical nonconvex NLP of the power
systems community: minimise generation cost subject to nonlinear AC power
flow equations, voltage limits, line thermal limits, and generator output
bounds. Problems come from the MATPOWER test case library (Zimmerman,
Murillo-Sanchez, Thomas, "MATPOWER: Steady-State Operations, Planning, and
Analysis Tools for Power Systems Research and Education", IEEE Trans. Power
Syst. 2011). ripopt uses the polar form with explicit bus voltage magnitude
and angle variables.

Problems currently range from case3_lmbd (3 buses, 9 variables) to
case30_ieee (30 buses, 60+ variables with up to 40 equality and 30
inequality constraints). These are good stress tests for dense KKT paths
because the constraint Jacobian is moderately dense and the Hessian has
many non-convex blocks from the voltage–angle coupling.

The suite name is "grid" rather than "opf" so that the directory layout
makes its purpose obvious at a glance — every problem here is an electrical
grid optimization.

## Contents

- `problems.rs` — MATPOWER test cases wrapped as ripopt `NlpProblem`
  instances; shared by `tests/grid.rs` and `examples/grid_benchmark.rs`
- `grid_results.json` — latest comparative results (ripopt vs Ipopt)
- `grid_benchmark_report.md` — per-problem analysis

## Prerequisites

None beyond ripopt itself (`ipopt-native` feature for the comparison).

## How to run

From the repo root:

```bash
make grid-run
```

or directly:

```bash
RESULTS_FILE=benchmarks/grid/grid_results.json \
    cargo run --release --features ipopt-native --example grid_benchmark
```

As tests:

```bash
cargo test --release --test grid
```

## Output

- `grid_results.json` — ripopt and Ipopt per-problem results
- `grid_stderr.txt` — solver chatter (gitignored)

This suite feeds the composite `benchmarks/BENCHMARK_REPORT.md`.
