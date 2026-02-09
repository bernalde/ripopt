# Performance Comparison: ripopt vs Ipopt (Native C++)

Benchmark of ripopt (pure Rust interior point optimizer) against Ipopt 3.14.19
(compiled C++, MUMPS linear solver) on the Hock-Schittkowski test suite (120 problems).

Both solvers use the same settings: adaptive mu strategy, tol=1e-8, max_iter=3000.
Timing: minimum of 5 runs per problem. Platform: Apple M3 (aarch64-apple-darwin).

## Solve Success Rates

| Metric           | ripopt              | Ipopt               |
|------------------|---------------------|---------------------|
| Optimal          | 117                 | 118                 |
| Acceptable       | 2                   | 0                   |
| **Total solved** | **119/120** (99.2%) | **118/120** (98.3%) |
| Failed           | 1                   | 2                   |
| Both solved      | 117                 | 117                 |

## Solution Accuracy (where both solve)

Relative objective difference: |f_ripopt - f_ipopt| / max(|f_ripopt|, |f_ipopt|, 1)

| Metric                     | Value           |
|----------------------------|-----------------|
| Matching (rel diff < 1e-4) | 107/117 (91.5%) |
| Mean rel diff              | 5.00e-02        |
| Median rel diff            | 2.88e-09        |
| 95th percentile            | 3.72e-01        |
| Max rel diff               | 1.00e+00        |

## Iteration Counts (where both solve)

| Statistic | ripopt | Ipopt |
|-----------|--------|-------|
| Mean      | 14.5   | 13.3  |
| Median    | 11     | 10    |
| 95th pctl | 36     | 39    |
| Max       | 73     | 68    |
| Total     | 1691   | 1561  |

- ripopt uses fewer iterations: 38/117 problems
- Ipopt uses fewer iterations: 57/117 problems
- Same iteration count: 22/117 problems

## Solve Time (where both solve)

| Statistic | ripopt | Ipopt   | Speedup |
|-----------|--------|---------|---------|
| Mean      | 31us   | 2.2ms   | 70x     |
| Median    | 13us   | 1.8ms   | 138x    |
| Min       | 3us    | 370us   |         |
| Max       | 1.0ms  | 11.4ms  |         |
| Total     | 3.7ms  | 255.9ms | 70x     |

- **Geometric mean speedup: 124.4x** (ripopt faster)
- ripopt faster on 117/117 problems
- Ipopt faster on 0/117 problems

## Per-Iteration Time

| Statistic | ripopt | Ipopt |
|-----------|--------|-------|
| Mean      | 2us    | 187us |
| Median    | 1us    | 183us |

ripopt's per-iteration advantage comes from its lean dense implementation:
direct array operations and dense LDL^T factorization with minimal overhead.
Ipopt's higher per-iteration cost reflects its large-scale infrastructure:
sparse matrix management (MUMPS), memory allocation, and option handling.

## Breakdown by Problem Type

| Category      | Count | ripopt mean | Ipopt mean | Geo mean speedup |
|---------------|-------|-------------|------------|------------------|
| Constrained   | 98    | 34us        | 2.2ms      | 116.5x           |
| Unconstrained | 19    | 16us        | 2.3ms      | 175.2x           |

## Speedup Distribution

| Speedup range      | Count | Percentage |
|--------------------|-------|------------|
| 1-10x              | 2     | 2%         |
| 10-50x             | 6     | 5%         |
| 50-100x            | 26    | 22%        |
| 100-200x           | 61    | 52%        |
| 200-500x           | 20    | 17%        |
| 500-1000x          | 2     | 2%         |
| <1x (Ipopt faster) | 0     | 0%         |

## Detailed Per-Problem Results

| TP# | n  | m  | ripopt        | Ipopt                 | r_iter | i_iter | r_time  | i_time  | Speedup | Obj Diff |
|-----|----|----|---------------|-----------------------|--------|--------|---------|---------|---------|----------|
| 001 | 2  | 0  | Optimal       | Optimal               | 26     | 28     | 15us    | 4.6ms   | 314x    | 2.3e-14  |
| 002 | 2  | 0  | Optimal       | Optimal               | 8      | 10     | 6us     | 1.8ms   | 314x    | 5.5e-09  |
| 003 | 2  | 0  | Optimal       | Optimal               | 5      | 4      | 4us     | 841us   | 189x    | 1.0e-08  |
| 004 | 2  | 0  | Optimal       | Optimal               | 5      | 4      | 4us     | 846us   | 191x    | 2.8e-08  |
| 005 | 2  | 0  | Optimal       | Optimal               | 10     | 7      | 9us     | 1.2ms   | 144x    | 2.3e-16  |
| 006 | 2  | 1  | Optimal       | Optimal               | 2      | 5      | 5us     | 957us   | 205x    | 0.0e+00  |
| 007 | 2  | 1  | Optimal       | Optimal               | 13     | 27     | 13us    | 3.8ms   | 288x    | 8.6e-10  |
| 009 | 2  | 1  | Optimal       | Optimal               | 3      | 3      | 6us     | 642us   | 109x    | 7.1e-15  |
| 010 | 2  | 1  | Optimal       | Optimal               | 19     | 12     | 14us    | 2.0ms   | 143x    | 5.0e-09  |
| 011 | 2  | 1  | Optimal       | Optimal               | 6      | 6      | 6us     | 1.2ms   | 202x    | 3.6e-09  |
| 012 | 2  | 1  | Optimal       | Optimal               | 7      | 6      | 8us     | 1.2ms   | 146x    | 1.7e-10  |
| 013 | 2  | 1  | Acceptable    | Optimal               | 35     | 47     | 25us    | 7.2ms   | 288x    | 9.4e-04  |
| 014 | 2  | 2  | Optimal       | Optimal               | 5      | 5      | 6us     | 1.1ms   | 174x    | 1.3e-08  |
| 015 | 2  | 2  | Optimal       | Optimal               | 9      | 11     | 15us    | 1.9ms   | 123x    | 8.0e-08  |
| 016 | 2  | 2  | Optimal       | Optimal               | 12     | 10     | 11us    | 1.9ms   | 164x    | 9.9e-01  |
| 017 | 2  | 2  | Optimal       | Optimal               | 16     | 18     | 13us    | 3.0ms   | 228x    | 2.6e-08  |
| 018 | 2  | 2  | Optimal       | Optimal               | 12     | 10     | 12us    | 1.7ms   | 145x    | 7.1e-10  |
| 019 | 2  | 2  | Optimal       | Optimal               | 15     | 12     | 14us    | 2.2ms   | 160x    | 3.2e-09  |
| 020 | 2  | 3  | Optimal       | Optimal               | 11     | 11     | 21us    | 2.2ms   | 107x    | 7.0e-08  |
| 021 | 2  | 1  | Optimal       | Optimal               | 9      | 6      | 8us     | 1.2ms   | 153x    | 1.4e-10  |
| 022 | 2  | 2  | Optimal       | Optimal               | 6      | 5      | 6us     | 1.0ms   | 166x    | 3.2e-05  |
| 023 | 2  | 5  | Optimal       | Optimal               | 25     | 9      | 42us    | 1.6ms   | 37x     | 8.7e-01  |
| 024 | 2  | 3  | Optimal       | Optimal               | 12     | 14     | 24us    | 2.7ms   | 115x    | 1.5e-08  |
| 026 | 3  | 1  | Optimal       | Optimal               | 20     | 25     | 17us    | 2.9ms   | 178x    | 6.7e-16  |
| 027 | 3  | 1  | Optimal       | Optimal               | 27     | 21     | 33us    | 2.8ms   | 84x     | 3.9e-11  |
| 028 | 3  | 1  | Optimal       | Optimal               | 1      | 1      | 4us     | 432us   | 123x    | 0.0e+00  |
| 029 | 3  | 1  | Optimal       | Optimal               | 13     | 7      | 17us    | 1.4ms   | 79x     | 3.1e-10  |
| 030 | 3  | 1  | Optimal       | Optimal               | 10     | 7      | 9us     | 1.3ms   | 144x    | 1.4e-06  |
| 031 | 3  | 1  | Optimal       | Optimal               | 11     | 6      | 11us    | 1.1ms   | 104x    | 1.0e-08  |
| 032 | 3  | 2  | Optimal       | Optimal               | 15     | 15     | 14us    | 2.4ms   | 175x    | 4.6e-08  |
| 033 | 3  | 2  | Optimal       | Optimal               | 11     | 9      | 17us    | 1.6ms   | 95x     | 2.8e-08  |
| 034 | 3  | 2  | Optimal       | Optimal               | 11     | 7      | 15us    | 1.3ms   | 83x     | 1.1e-08  |
| 035 | 3  | 1  | Optimal       | Optimal               | 9      | 7      | 9us     | 1.3ms   | 146x    | 3.8e-09  |
| 036 | 3  | 1  | Optimal       | Optimal               | 13     | 11     | 16us    | 2.0ms   | 126x    | 6.3e-09  |
| 037 | 3  | 2  | Optimal       | Optimal               | 9      | 11     | 14us    | 2.0ms   | 147x    | 4.2e-10  |
| 038 | 4  | 0  | Optimal       | Optimal               | 9      | 39     | 9us     | 6.3ms   | 674x    | 1.0e+00  |
| 039 | 4  | 2  | Optimal       | Optimal               | 23     | 13     | 33us    | 1.8ms   | 54x     | 6.6e-09  |
| 040 | 4  | 3  | Optimal       | Optimal               | 4      | 3      | 8us     | 619us   | 83x     | 7.9e-11  |
| 041 | 4  | 1  | Optimal       | Optimal               | 11     | 7      | 12us    | 1.3ms   | 111x    | 1.6e-09  |
| 042 | 4  | 2  | Optimal       | Optimal               | 5      | 4      | 7us     | 721us   | 105x    | 4.7e-12  |
| 043 | 4  | 3  | Optimal       | Optimal               | 8      | 8      | 13us    | 1.5ms   | 120x    | 6.8e-10  |
| 044 | 4  | 6  | Optimal       | Optimal               | 15     | 24     | 35us    | 4.3ms   | 123x    | 1.0e-08  |
| 045 | 5  | 0  | Optimal       | Optimal               | 73     | 11     | 99us    | 2.1ms   | 21x     | 6.5e-08  |
| 046 | 5  | 2  | Optimal       | Optimal               | 25     | 19     | 26us    | 2.5ms   | 96x     | 2.1e-02  |
| 047 | 5  | 3  | Optimal       | Optimal               | 13     | 15     | 26us    | 2.1ms   | 79x     | 4.4e-16  |
| 048 | 5  | 2  | Optimal       | Optimal               | 1      | 1      | 4us     | 408us   | 97x     | 0.0e+00  |
| 049 | 5  | 2  | Optimal       | Optimal               | 20     | 19     | 19us    | 2.2ms   | 119x    | 8.5e-12  |
| 050 | 5  | 3  | Optimal       | Optimal               | 8      | 8      | 11us    | 1.1ms   | 103x    | 0.0e+00  |
| 051 | 5  | 3  | Optimal       | Optimal               | 1      | 1      | 4us     | 411us   | 92x     | 0.0e+00  |
| 052 | 5  | 3  | Optimal       | Optimal               | 1      | 1      | 4us     | 389us   | 90x     | 5.0e-16  |
| 053 | 5  | 3  | Optimal       | Optimal               | 8      | 6      | 15us    | 1.1ms   | 75x     | 6.3e-13  |
| 056 | 7  | 4  | Optimal       | Optimal               | 2      | 3      | 10us    | 829us   | 85x     | 9.4e-49  |
| 058 | 2  | 3  | Optimal       | Optimal               | 11     | 11     | 13us    | 2.0ms   | 152x    | 1.3e-08  |
| 060 | 3  | 1  | Optimal       | Optimal               | 8      | 6      | 8us     | 1.1ms   | 135x    | 1.2e-13  |
| 061 | 3  | 2  | Optimal       | Optimal               | 9      | 10     | 11us    | 1.3ms   | 115x    | 1.2e-15  |
| 063 | 3  | 2  | Optimal       | Optimal               | 9      | 5      | 12us    | 1.0ms   | 84x     | 1.2e-16  |
| 064 | 3  | 1  | Optimal       | Optimal               | 16     | 16     | 14us    | 2.6ms   | 187x    | 1.6e-09  |
| 065 | 3  | 1  | Optimal       | Optimal               | 11     | 16     | 13us    | 3.0ms   | 226x    | 6.5e-09  |
| 066 | 3  | 2  | Optimal       | Optimal               | 10     | 10     | 12us    | 1.6ms   | 135x    | 1.1e-08  |
| 071 | 4  | 2  | Optimal       | Optimal               | 10     | 12     | 12us    | 2.0ms   | 167x    | 4.5e-09  |
| 072 | 4  | 2  | Optimal       | Optimal               | 28     | 16     | 34us    | 2.6ms   | 77x     | 6.7e-07  |
| 076 | 4  | 3  | Optimal       | Optimal               | 11     | 7      | 13us    | 1.3ms   | 99x     | 4.7e-09  |
| 077 | 5  | 2  | Optimal       | Optimal               | 9      | 11     | 13us    | 1.5ms   | 115x    | 2.6e-11  |
| 078 | 5  | 3  | Optimal       | Optimal               | 4      | 4      | 8us     | 710us   | 91x     | 3.7e-10  |
| 079 | 5  | 3  | Optimal       | Optimal               | 4      | 4      | 8us     | 728us   | 97x     | 2.6e-10  |
| 080 | 5  | 3  | Optimal       | Optimal               | 9      | 5      | 15us    | 1.1ms   | 72x     | 5.4e-13  |
| 081 | 5  | 3  | Optimal       | Optimal               | 73     | 68     | 174us   | 11.4ms  | 65x     | 9.5e-01  |
| 106 | 8  | 6  | Optimal       | Optimal               | 17     | 18     | 51us    | 2.8ms   | 55x     | 1.3e-08  |
| 108 | 9  | 13 | Optimal       | Optimal               | 19     | 23     | 124us   | 4.2ms   | 34x     | 5.7e-09  |
| 113 | 10 | 8  | Optimal       | Optimal               | 13     | 9      | 57us    | 1.7ms   | 30x     | 1.6e-09  |
| 114 | 10 | 11 | Optimal       | Optimal               | 30     | 13     | 147us   | 2.2ms   | 15x     | 1.1e-07  |
| 116 | 13 | 15 | Optimal       | Optimal               | 56     | 17     | 497us   | 3.2ms   | 6x      | 3.7e-07  |
| 201 | 2  | 0  | Optimal       | Optimal               | 1      | 1      | 3us     | 396us   | 127x    | 0.0e+00  |
| 206 | 2  | 0  | Optimal       | Optimal               | 4      | 4      | 4us     | 663us   | 156x    | 8.9e-16  |
| 211 | 2  | 0  | Optimal       | Optimal               | 27     | 27     | 13us    | 3.6ms   | 273x    | 1.2e-14  |
| 212 | 2  | 0  | Optimal       | Optimal               | 8      | 11     | 6us     | 1.6ms   | 249x    | 5.2e-23  |
| 213 | 2  | 0  | Acceptable    | Optimal               | 44     | 22     | 38us    | 2.4ms   | 63x     | 6.2e-06  |
| 214 | 2  | 0  | Optimal       | InvalidNumberDetected | 42     | 32     | 50us    | 6.7ms   | N/A     | N/A      |
| 215 | 2  | 1  | Optimal       | Optimal               | 9      | 13     | 8us     | 2.3ms   | 279x    | 1.4e-08  |
| 216 | 2  | 1  | Optimal       | Optimal               | 9      | 7      | 10us    | 1.4ms   | 144x    | 4.6e-13  |
| 217 | 2  | 2  | Optimal       | Optimal               | 10     | 8      | 9us     | 1.6ms   | 178x    | 1.2e-08  |
| 218 | 2  | 1  | Optimal       | Optimal               | 19     | 6      | 14us    | 1.3ms   | 91x     | 1.0e-08  |
| 219 | 4  | 2  | Optimal       | Optimal               | 26     | 50     | 27us    | 5.9ms   | 220x    | 1.5e-13  |
| 220 | 2  | 1  | Optimal       | Optimal               | 4      | 3      | 5us     | 781us   | 146x    | 1.0e-08  |
| 221 | 2  | 1  | Optimal       | Optimal               | 8      | 37     | 7us     | 5.2ms   | 731x    | 1.0e-08  |
| 223 | 2  | 2  | Optimal       | Infeasible            | 4      | 1567   | 10us    | 593.0ms | N/A     | N/A      |
| 224 | 2  | 4  | Optimal       | Optimal               | 15     | 9      | 15us    | 1.8ms   | 118x    | 8.3e-10  |
| 225 | 2  | 5  | Optimal       | Optimal               | 13     | 10     | 19us    | 1.8ms   | 94x     | 7.9e-01  |
| 226 | 2  | 2  | Optimal       | Optimal               | 13     | 11     | 14us    | 2.0ms   | 151x    | 5.0e-09  |
| 227 | 2  | 2  | Optimal       | Optimal               | 6      | 10     | 9us     | 1.9ms   | 215x    | 2.0e-08  |
| 228 | 2  | 2  | Optimal       | Optimal               | 5      | 7      | 8us     | 1.3ms   | 178x    | 2.2e-09  |
| 229 | 2  | 0  | Optimal       | Optimal               | 16     | 21     | 11us    | 3.7ms   | 354x    | 4.2e-03  |
| 230 | 2  | 2  | Optimal       | Optimal               | 2      | 7      | 4us     | 1.3ms   | 352x    | 9.7e-09  |
| 232 | 2  | 3  | Optimal       | Optimal               | 9      | 8      | 14us    | 1.5ms   | 111x    | 1.8e-08  |
| 234 | 2  | 1  | Optimal       | Optimal               | 9      | 18     | 9us     | 2.9ms   | 326x    | 2.7e-01  |
| 235 | 3  | 1  | Optimal       | Optimal               | 29     | 13     | 32us    | 1.9ms   | 58x     | 3.6e-14  |
| 240 | 3  | 0  | Optimal       | Optimal               | 1      | 1      | 3us     | 370us   | 111x    | 3.2e-29  |
| 248 | 3  | 2  | Optimal       | Optimal               | 10     | 14     | 14us    | 2.4ms   | 177x    | 5.5e-08  |
| 249 | 3  | 1  | Optimal       | Optimal               | 22     | 7      | 15us    | 1.3ms   | 84x     | 1.0e-08  |
| 250 | 3  | 2  | Optimal       | Optimal               | 13     | 14     | 18us    | 2.5ms   | 140x    | 6.3e-09  |
| 251 | 3  | 1  | Optimal       | Optimal               | 9      | 10     | 12us    | 1.9ms   | 159x    | 4.2e-10  |
| 252 | 3  | 1  | Optimal       | Optimal               | 21     | 15     | 15us    | 2.4ms   | 165x    | 4.0e-10  |
| 254 | 3  | 2  | Optimal       | Optimal               | 17     | 29     | 14us    | 4.3ms   | 307x    | 6.8e-12  |
| 255 | 4  | 0  | Optimal       | Optimal               | 13     | 13     | 13us    | 2.5ms   | 191x    | 2.0e-08  |
| 256 | 4  | 0  | Optimal       | Optimal               | 20     | 19     | 12us    | 2.2ms   | 174x    | 5.3e-12  |
| 257 | 4  | 0  | Optimal       | Optimal               | 12     | 7      | 12us    | 1.4ms   | 120x    | 4.9e-15  |
| 258 | 4  | 0  | Optimal       | Optimal               | 42     | 40     | 25us    | 5.5ms   | 217x    | 7.8e-14  |
| 259 | 4  | 0  | Optimal       | Optimal               | 11     | 9      | 9us     | 1.9ms   | 207x    | 2.3e-14  |
| 262 | 4  | 4  | Optimal       | Optimal               | 16     | 7      | 29us    | 1.3ms   | 47x     | 2.4e-09  |
| 263 | 4  | 4  | Optimal       | Optimal               | 21     | 58     | 23us    | 8.5ms   | 360x    | 2.9e-09  |
| 264 | 4  | 3  | Optimal       | Optimal               | 8      | 13     | 14us    | 2.4ms   | 173x    | 6.7e-10  |
| 270 | 5  | 1  | Optimal       | Optimal               | 15     | 11     | 19us    | 2.4ms   | 125x    | 2.1e-12  |
| 325 | 2  | 3  | Optimal       | Optimal               | 12     | 14     | 22us    | 3.2ms   | 146x    | 1.9e-09  |
| 335 | 3  | 2  | Optimal       | Optimal               | 26     | 25     | 22us    | 2.9ms   | 129x    | 4.3e-13  |
| 338 | 3  | 2  | Optimal       | Optimal               | 26     | 45     | 34us    | 5.7ms   | 171x    | 4.8e-16  |
| 339 | 3  | 1  | Optimal       | Optimal               | 10     | 7      | 12us    | 1.3ms   | 111x    | 3.3e-09  |
| 344 | 3  | 1  | Optimal       | Optimal               | 7      | 7      | 7us     | 973us   | 134x    | 5.9e-14  |
| 354 | 4  | 1  | Optimal       | Optimal               | 9      | 10     | 10us    | 1.7ms   | 170x    | 4.9e-09  |
| 374 | 10 | 35 | MaxIterations | Optimal               | 2999   | 92     | 224.9ms | 25.7ms  | N/A     | N/A      |
| 376 | 10 | 15 | Optimal       | Optimal               | 59     | 23     | 1.0ms   | 4.7ms   | 5x      | 9.6e-01  |

## Failure Analysis

### ripopt-only failures (1)

| TP# | n  | m  | ripopt status | Ipopt objective |
|-----|----|----|---------------|-----------------|
| 374 | 10 | 35 | MaxIterations | 2.332635e-01    |

### Ipopt-only failures (2)

| TP# | n | m | Ipopt status          | ripopt objective |
|-----|---|---|-----------------------|------------------|
| 214 | 2 | 0 | InvalidNumberDetected | 0.000000e+00     |
| 223 | 2 | 2 | Infeasible            | N/A              |

## Notes

- These are small dense problems (n=2-13, m=0-35). For large sparse problems,
  Ipopt's MUMPS sparse solver would have an advantage over ripopt's O(n^3) dense LDL^T.
- ripopt's speed advantage is most relevant for embedded optimization, real-time
  control, or applications solving many small NLPs where per-solve overhead dominates.
- Ipopt version: 3.14.19 (Homebrew, aarch64-apple-darwin, MUMPS linear solver)
- ripopt version: 0.1.0 (pure Rust, no external linear algebra dependencies)
