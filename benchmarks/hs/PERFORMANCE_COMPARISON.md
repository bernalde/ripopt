# Performance Comparison: ripopt vs Ipopt (Native C++)

Benchmark of ripopt (pure Rust interior point optimizer) against Ipopt 3.14.19
(compiled C++, MUMPS linear solver) on the Hock-Schittkowski test suite (120 problems).

Both solvers use the same settings: adaptive mu strategy, tol=1e-8, max_iter=3000.
Strict `Optimal` status is required for a problem to count as solved. Platform:
Apple Mac Mini (aarch64-apple-darwin).

## Solve Success Rates

| Metric           | ripopt              | Ipopt               |
|------------------|---------------------|---------------------|
| Optimal          | 118                 | 116                 |
| **Total solved** | **118/120** (98.3%) | **116/120** (96.7%) |
| Failed           | 2                   | 4                   |
| Both solved      | 116                 | 116                 |

**ripopt-only solves (2):** HS214, HS223 — ripopt succeeds where Ipopt reports
`InvalidNumberDetected` or `Infeasible`.

**Ipopt-only solves:** 0.

## Solution Accuracy (where both solve)

On the 116 problems where both solvers reach `Optimal`:

| Metric                     | Value            |
|----------------------------|------------------|
| Matching (rel diff < 1e-4) | 110/116 (94.8%)  |

Relative objective difference is `|f_ripopt - f_ipopt| / max(|f_ripopt|, |f_ipopt|, 1)`.
The 6 non-matching cases reach valid KKT points but at different local optima.

## Iteration Counts (where both solve)

| Statistic | ripopt | Ipopt |
|-----------|--------|-------|
| Mean      | 70.8   | 13.3  |
| Median    | 13     | 10    |

Ipopt uses fewer iterations on average, reflecting its more mature barrier parameter
strategy and line search tuning. ripopt's lower per-iteration cost more than
compensates in wall-clock time.

## Solve Time (where both solve)

| Statistic | ripopt | Ipopt   | Speedup |
|-----------|--------|---------|---------|
| Median    | 118us  | 1.9ms   | 16.2x   |
| Total     | 53.6ms | 281.9ms | 5.3x    |

- **Geometric mean speedup**: **15.0x**
- **Median speedup**: **14.2x**
- ripopt faster: **113/116** (97%)
- ripopt 10x+ faster: **84/116** (72%)
- Ipopt faster: 3/116

The speed advantage on HS comes from:

1. **Lower per-iteration overhead.** Stack allocation, no C/Fortran interop, and
   cache-efficient dense Bunch-Kaufman factorization dominate on small problems
   (HS problems are all n ≤ 15).
2. **Mehrotra predictor-corrector** with Gondzio centrality corrections (enabled
   by default), which reduces iteration counts on many well-conditioned problems.
3. **Dense linear algebra throughout.** All HS problems fit below the sparse
   threshold (n+m < 110), so ripopt never pays the sparse symbolic-analysis cost.
