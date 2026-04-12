# Electrolyte Thermodynamics Suite

Gibbs free energy minimisation for aqueous electrolyte systems. Given a set
of species and their reference-state thermodynamic data, find the
equilibrium composition subject to element-balance (mass conservation) and
charge-neutrality constraints. These problems are notoriously ill-conditioned
because species concentrations span many orders of magnitude and the
logarithmic activity terms blow up near zero concentration.

The 13 problems in `problems.rs` cover mixed acid-base, hydrolysis, and
redox systems drawn from classic chemical engineering thermodynamics
textbook examples and ripopt's own regression set.

## Contents

- `problems.rs` — NlpProblem definitions for all 13 systems; shared by
  `tests/electrolyte_thermo.rs` and `examples/electrolyte_benchmark.rs`
- `electrolyte_results.json` — latest comparative results
- `electrolyte_benchmark_report.md` — per-problem analysis

## Prerequisites

None beyond ripopt itself (`ipopt-native` feature required for the
comparative benchmark).

## How to run

From the repo root:

```bash
make electrolyte-run
```

or directly:

```bash
RESULTS_FILE=benchmarks/electrolyte/electrolyte_results.json \
    cargo run --release --features ipopt-native --example electrolyte_benchmark
```

The test suite form:

```bash
cargo test --release --test electrolyte_thermo
```

## Output

- `electrolyte_results.json` — ripopt and Ipopt per-problem results
- `electrolyte_stderr.txt` — solver chatter (gitignored)

This suite feeds the composite `benchmarks/BENCHMARK_REPORT.md`.
