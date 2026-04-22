# Benchmarks

ripopt is benchmarked against Ipopt (native C++ with MUMPS linear solver) on three standard test suites and several domain-specific problem sets. All benchmarks run on the same hardware; both solvers receive identical problem data through the same Rust trait interface.

## Hock-Schittkowski Suite (120 problems)

The HS suite is the classic test set for NLP solvers, covering small problems (n ≤ 15) with mixed equality/inequality constraints.

| Metric | ripopt | Ipopt (MUMPS) |
|---|---|---|
| Problems solved | **118/120 (98.3%)** | 116/120 (96.7%) |
| Solved by ripopt only | **2** | — |
| Solved by Ipopt only | — | 0 |

On 116 commonly solved problems (strict `Optimal` status required):

| Metric | Value |
|---|---|
| Geometric mean speedup | **15.0x** |
| Median speedup | **14.2x** |
| ripopt faster | 113/116 (97%) |
| ripopt 10x+ faster | 84/116 (72%) |
| Matching objectives (rel diff < 1e-4) | 110/116 (95%) |

Run: `make hs-run`

## CUTEst Benchmark Suite (727 problems)

CUTEst covers a wide range of problem types, sizes, and structures. Problems range from n=2 to n=10,000+.

| Metric | ripopt | Ipopt (MUMPS) |
|---|---|---|
| Total solved | **562/727 (77.3%)** | 561/727 (77.2%) |
| Solved by ripopt only | **37** | — |
| Solved by Ipopt only | — | 36 |

On 525 commonly solved problems:

| Metric | Value |
|---|---|
| Geometric mean speedup | **9.9x** |
| Median speedup | **18.9x** |
| ripopt faster | 440/525 (84%) |
| ripopt 10x+ faster | 332/525 (63%) |
| Matching objectives (rel diff < 1e-4) | 436/525 (83.0%) |

Run: `make benchmark` (full suite, ~2 hours) or individual problems:
```bash
cargo run --bin cutest_suite --features cutest,ipopt-native --release -- PROBLEM1 PROBLEM2
```

### Where ripopt is faster

1. **Small problems (n < 50).** Stack allocation and cache-efficient dense BK factorization avoid sparse overhead. 2–5x per-iteration speedup over Ipopt.
2. **Tall-narrow problems (m >> n, n ≤ 100).** Dense condensed KKT via Schur complement reduces an (n+m)×(n+m) sparse factorization to n×n dense. 100–800x speedup on EXPFITC (n=5, m=502), OET3 (n=4, m=1002).
3. **Better iteration counts.** Mehrotra predictor-corrector with Gondzio centrality corrections (enabled by default) cuts iterations by 20–40% on many problems.
4. **Fallback cascade.** NE-to-LS reformulation, two-phase restoration, and multi-solver fallback recover problems Ipopt cannot solve.

### Where Ipopt is faster

1. **Large sparse problems (n+m > 5,000).** Ipopt's Fortran MUMPS is 10–15x faster per factorization than rmumps on large systems.
2. **Some medium constrained problems.** A handful of problems (CORE1, HAIFAM, NET1) have higher per-iteration cost in ripopt's line search or fallback cascade.

## Large-Scale Benchmarks

Both solvers use the same Rust trait interface; ripopt uses rmumps, Ipopt uses Fortran MUMPS.

| Problem | n | m | ripopt | time | Ipopt | time | speedup |
|---|---|---|---|---|---|---|---|
| Rosenbrock 500 | 500 | 0 | Optimal | 0.003s | Optimal | 0.199s | **76.2x** |
| Bratu 1K | 1,000 | 998 | Optimal | 0.002s | Optimal | 0.002s | 1.1x |
| SparseQP 1K | 500 | 500 | Optimal | 0.176s | Optimal | 0.004s | 0.02x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 0.006s | Optimal | 0.002s | 0.4x |

Numbers above are fresh (v0.7.0). The larger problems (Poisson 2.5K,
Rosenbrock 5K, Bratu 10K, OptControl 20K, Poisson 50K, SparseQP 100K)
are not re-run for v0.7.0 — the stricter KKT backward-error probe
causes Poisson 2.5K to exhaust `max_iter` and the full sweep is gated
behind a separate investigation. Historical timings from v0.6.2 remain
in `benchmarks/large_scale/large_scale_results.txt` snapshots.

Run: `make benchmark`

## Domain-Specific Benchmarks

| Suite | Problems | ripopt | Ipopt | Notes |
|---|---|---|---|---|
| Electrolyte thermodynamics | 13 | **13/13 (100%)** | 12/13 (92.3%) | 17.5x geo mean; ripopt uniquely solves seawater speciation |
| Grid (AC Optimal Power Flow) | 4 | 3/4 (75%) | **4/4 (100%)** | 2.8x geo mean on 3 commonly-solved; case30_ieee regression |
| CHO parameter estimation | 1 | 0/1 | 0/1 | n=21,672, m=21,660; both hit iteration limit |
| Gas pipeline NLPs | 4 | see suite README | see suite README | PDE-discretized Euler equations on pipe networks (gaslib11/40). Standalone — does not feed `BENCHMARK_REPORT.md` |
| Water distribution NLPs | 6 | see suite README | see suite README | MINLPLib water network design instances. Standalone — does not feed `BENCHMARK_REPORT.md` |
