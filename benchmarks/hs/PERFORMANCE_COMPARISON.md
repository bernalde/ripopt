# Performance Comparison: ripopt vs Ipopt (Native C++)

Benchmark of ripopt (pure Rust interior point optimizer) against Ipopt 3.14.19
(compiled C++, MUMPS linear solver) on the Hock-Schittkowski test suite (120 problems).

Both solvers use the same settings: adaptive mu strategy, tol=1e-8, max_iter=3000.
Strict `Optimal` status is required for a problem to count as solved. Platform:
Apple Mac Mini (aarch64-apple-darwin).

## Solve Success Rates

| Metric           | ripopt              | Ipopt               |
|------------------|---------------------|---------------------|
| Optimal          | 116                 | 116                 |
| **Total solved** | **116/120** (96.7%) | **116/120** (96.7%) |
| Failed           | 4                   | 4                   |
| Both solved      | 114                 | 114                 |

**ripopt-only solves (2):** HS214, HS223 — ripopt succeeds where Ipopt reports
`InvalidNumberDetected` or `Infeasible`.

**Ipopt-only solves (2):** HS013 and HS225 (ripopt: `NumericalError` on both).

**ripopt failure modes:** 3 `NumericalError`, 1 `RestorationFailed`.

## Solution Accuracy (where both solve)

On the 114 problems where both solvers reach `Optimal`:

| Metric                     | Value           |
|----------------------------|-----------------|
| Matching (rel diff < 1e-4) | 106/114 (93.0%) |

Relative objective difference is `|f_ripopt - f_ipopt| / max(|f_ripopt|, |f_ipopt|, 1)`.
The 8 non-matching cases reach valid KKT points but at different local optima.

## Iteration Counts (where both solve)

| Statistic | ripopt | Ipopt |
|-----------|--------|-------|
| Mean      | 15.1   | 13.0  |
| Median    | 13     | 10    |

Ipopt uses fewer iterations on average, reflecting its more mature barrier parameter
strategy and line search tuning. ripopt's lower per-iteration cost more than
compensates in wall-clock time.

## Solve Time (where both solve)

| Statistic | ripopt | Ipopt   | Speedup |
|-----------|--------|---------|---------|
| Median    | 94us   | 1.8ms   | 17.8x   |
| Total     | 23.1ms | 269.4ms | 11.7x   |

- **Geometric mean speedup**: **17.7x**
- **Median speedup**: **17.8x**
- ripopt faster: **113/114** (99%)
- ripopt 10x+ faster: **95/114** (83%)
- Ipopt faster: 1/114

The speed advantage on HS comes from:

1. **Lower per-iteration overhead.** Stack allocation, no C/Fortran interop, and
   cache-efficient dense Bunch-Kaufman factorization dominate on small problems
   (HS problems are all n ≤ 15).
2. **Mehrotra predictor-corrector** with Gondzio centrality corrections (enabled
   by default), which reduces iteration counts on many well-conditioned problems.
3. **Dense linear algebra throughout.** All HS problems fit below the sparse
   threshold (n+m < 110), so ripopt never pays the sparse symbolic-analysis cost.
