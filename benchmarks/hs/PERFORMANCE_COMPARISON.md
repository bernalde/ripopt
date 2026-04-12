# Performance Comparison: ripopt vs Ipopt (Native C++)

Benchmark of ripopt (pure Rust interior point optimizer) against Ipopt 3.14.19
(compiled C++, MUMPS linear solver) on the Hock-Schittkowski test suite (120 problems).

Both solvers use the same settings: adaptive mu strategy, tol=1e-8, max_iter=3000.
Strict `Optimal` status is required for a problem to count as solved. Platform:
Apple Mac Mini (aarch64-apple-darwin).

## Solve Success Rates

| Metric           | ripopt              | Ipopt               |
|------------------|---------------------|---------------------|
| Optimal          | 115                 | 116                 |
| **Total solved** | **115/120** (95.8%) | **116/120** (96.7%) |
| Failed           | 5                   | 4                   |
| Both solved      | 113                 | 113                 |

**ripopt-only solves (2):** HS214, HS223 — ripopt succeeds where Ipopt reports
`InvalidNumberDetected` or `Infeasible`.

**Ipopt-only solves (3):** HS004 (ripopt: `MaxIterations`), HS013
(ripopt: `NumericalError`), HS108 (ripopt: `NumericalError`).

**ripopt failure modes:** 3 `NumericalError`, 1 `MaxIterations`, 1 `RestorationFailed`.

## Solution Accuracy (where both solve)

On the 113 problems where both solvers reach `Optimal`:

| Metric                     | Value           |
|----------------------------|-----------------|
| Matching (rel diff < 1e-4) | 107/113 (94.7%) |

Relative objective difference is `|f_ripopt - f_ipopt| / max(|f_ripopt|, |f_ipopt|, 1)`.
The 6 non-matching cases reach valid KKT points but at different local optima.

## Iteration Counts (where both solve)

| Statistic | ripopt | Ipopt |
|-----------|--------|-------|
| Mean      | 28.1   | 13.0  |
| Median    | 13     | 10    |

Ipopt uses fewer iterations on average, reflecting its more mature barrier parameter
strategy and line search tuning. ripopt's lower per-iteration cost more than
compensates in wall-clock time.

## Solve Time (where both solve)

| Statistic | ripopt | Ipopt   | Speedup |
|-----------|--------|---------|---------|
| Median    | 112us  | 2.0ms   | 15.1x   |
| Total     | 42.2ms | 264.2ms | 6.3x    |

- **Geometric mean speedup**: **14.0x**
- **Median speedup**: **15.1x**
- ripopt faster: **111/113** (98%)
- ripopt 10x+ faster: **80/113** (71%)
- Ipopt faster: 2/113

The speed advantage on HS comes from:

1. **Lower per-iteration overhead.** Stack allocation, no C/Fortran interop, and
   cache-efficient dense Bunch-Kaufman factorization dominate on small problems
   (HS problems are all n ≤ 15).
2. **Mehrotra predictor-corrector** with Gondzio centrality corrections (enabled
   by default), which reduces iteration counts on many well-conditioned problems.
3. **Dense linear algebra throughout.** All HS problems fit below the sparse
   threshold (n+m < 110), so ripopt never pays the sparse symbolic-analysis cost.
