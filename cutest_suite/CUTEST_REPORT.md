# CUTEst Benchmark Report

Comparison of ripopt vs Ipopt (C++) on the CUTEst test set.

## Executive Summary

- **Total problems**: 35
- **ripopt solved**: 32/35 (91.4%)
- **Ipopt solved**: 35/35 (100.0%)
- **Both solved**: 32/35
- **Matching solutions** (rel obj diff < 1e-4): 25/32

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (25 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 9.88e-07 |
| Median | 4.03e-15 |
| Max    | 1.74e-05 |

**All both-solved** (32 problems, including 7 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 3.97e-02 |
| Median | 7.86e-09 |
| Max    | 9.92e-01 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 17 | 16 | 17 | 16 | 12 |
| unconstrained | 18 | 16 | 18 | 16 | 13 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 44us | 796us | 18.0x | PASS |
| AIRCRFTB | 8 | 0 | Optimal | Optimal | 2.63e-21 | 49 | 15 | 66us | 3.9ms | 58.9x | PASS |
| ALLINIT | 4 | 0 | Optimal | Optimal | 0.00e+00 | 15 | 20 | 24us | 4.9ms | 199.5x | PASS |
| BENNETT5LS | 3 | 0 | Acceptable | Optimal | 1.74e-05 | 40 | 21 | 3.6ms | 5.1ms | 1.4x | PASS |
| BIGGS3 | 6 | 0 | Optimal | Optimal | 5.29e-24 | 19 | 9 | 36us | 2.5ms | 69.2x | PASS |
| BIGGS5 | 6 | 0 | Optimal | Optimal | 5.66e-03 | 50 | 20 | 91us | 4.5ms | 50.0x | MISMATCH |
| BOX2 | 3 | 0 | Acceptable | Optimal | 6.33e-14 | 15 | 8 | 23us | 1.7ms | 72.7x | PASS |
| BQPGABIM | 50 | 0 | Acceptable | Optimal | 4.34e-06 | 11 | 12 | 10.0ms | 3.5ms | 0.3x | PASS |
| CHWIRUT2LS | 3 | 0 | Acceptable | Optimal | 8.86e-16 | 21 | 6 | 178us | 1.5ms | 8.5x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 11.2ms | 1.8ms | 0.2x | PASS |
| DECONVC | 63 | 1 | Optimal | Optimal | 1.16e-03 | 35 | 31 | 6.0ms | 11.0ms | 1.8x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 2 | 26 | 899us | 31.3ms | 34.8x | PASS |
| DECONVU | 63 | 0 | Acceptable | Optimal | 2.37e-07 | 100 | 279 | 729us | 95.6ms | 131.1x | PASS |
| DJTL | 2 | 0 | Acceptable | Optimal | 1.42e-15 | 2044 | 1524 | 2.6ms | 207.7ms | 78.4x | PASS |
| DNIEPER | 61 | 24 | Acceptable | Optimal | 2.40e-03 | 6993 | 23 | 1.84s | 5.9ms | 0.0x | MISMATCH |
| HAHN1LS | 7 | 0 | Acceptable | Optimal | 1.52e-01 | 195 | 78 | 3.0ms | 22.9ms | 7.7x | MISMATCH |
| HS109 | 9 | 10 | NumericalErr | Optimal | N/A | 81 | 14 | 102.6ms | 3.3ms | 0.0x | ripopt_FAIL |
| HS67 | 3 | 14 | Optimal | Optimal | 2.99e-09 | 12 | 9 | 172us | 2.7ms | 15.5x | PASS |
| HS83 | 5 | 3 | Acceptable | Optimal | 2.47e-06 | 8991 | 9 | 305.4ms | 2.0ms | 0.0x | PASS |
| HS84 | 5 | 3 | Optimal | Optimal | 7.79e-03 | 2 | 9 | 2.5ms | 2.9ms | 1.2x | MISMATCH |
| HS85 | 5 | 21 | Optimal | Optimal | 2.07e-08 | 58 | 13 | 30.5ms | 4.5ms | 0.1x | PASS |
| HS99EXP | 31 | 21 | Optimal | Optimal | 0.00e+00 | 9 | 17 | 823us | 5.6ms | 6.8x | PASS |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 13.3ms | 2.1ms | 0.2x | PASS |
| MEYER3 | 3 | 0 | MaxIteration | Optimal | N/A | 2999 | 193 | 72.2ms | 37.5ms | 0.5x | ripopt_FAIL |
| MGH10LS | 3 | 0 | NumericalErr | Optimal | N/A | 251 | 1598 | 57.1ms | 314.8ms | 5.5x | ripopt_FAIL |
| MINSURF | 64 | 0 | Optimal | Optimal | 0.00e+00 | 17 | 4 | 88us | 2.1ms | 24.2x | PASS |
| OPTCNTRL | 32 | 20 | Optimal | Optimal | 2.03e-09 | 36 | 9 | 2.8ms | 2.9ms | 1.0x | PASS |
| OSBORNEA | 5 | 0 | Acceptable | Optimal | 1.88e-18 | 93 | 64 | 477us | 15.0ms | 31.5x | PASS |
| QC | 9 | 4 | Optimal | Optimal | 1.08e-01 | 0 | 44 | 1.2ms | 13.4ms | 10.8x | MISMATCH |
| QCNEW | 9 | 3 | Optimal | Optimal | 1.06e-07 | 999 | 6 | 68.4ms | 1.5ms | 0.0x | PASS |
| RAT43LS | 4 | 0 | Optimal | Optimal | 9.92e-01 | 3 | 34 | 15us | 9.2ms | 621.7x | MISMATCH |
| SIM2BQP | 2 | 0 | Optimal | Optimal | 7.86e-09 | 6 | 5 | 48us | 1.5ms | 31.9x | PASS |
| SSINE | 3 | 2 | Acceptable | Optimal | 0.00e+00 | 2029 | 224 | 33.1ms | 42.9ms | 1.3x | PASS |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 4.03e-15 | 20 | 19 | 171.1ms | 4.2ms | 0.0x | PASS |
| TRIGGER | 7 | 6 | Acceptable | Optimal | 0.00e+00 | 15 | 15 | 283us | 3.9ms | 13.7x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 685.8 | 80.3 |
| Median | 31 | 15 |
| Total  | 21944 | 2570 |

- ripopt fewer iterations: 8/32
- Ipopt fewer iterations: 22/32
- Tied: 2/32

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 78.4ms | 16.3ms |
| Median | 1.2ms | 4.2ms |
| Total  | 2.51s | 521.3ms |

- Geometric mean speedup (Ipopt_time/ripopt_time): **4.00x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 24/32 problems
- Ipopt faster: 8/32 problems
- Overall speedup (total time): 0.21x

## Failure Analysis

### Problems where only ripopt fails (3)

| Problem | n | m | ripopt status | Ipopt obj |
|---------|---|---|---------------|-----------|
| HS109 | 9 | 10 | NumericalError | 5.362069e+03 |
| MEYER3 | 3 | 0 | MaxIterations | 8.794586e+01 |
| MGH10LS | 3 | 0 | NumericalError | 8.794586e+01 |

### Objective mismatches (7)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 5
- **Convergence gap** (one Acceptable): 2
- **Better objective found by**: ripopt 1, Ipopt 6

| Problem | ripopt obj | Ipopt obj | Rel Diff | r_status | i_status | Better |
|---------|-----------|-----------|----------|----------|----------|--------|
| RAT43LS | 1.076462e+06 | 8.786405e+03 | 9.92e-01 | Optimal | Optimal | ipopt |
| HAHN1LS | 3.937328e+01 | 3.338424e+01 | 1.52e-01 | Acceptable | Optimal | ipopt |
| QC | -8.533467e+02 | -9.565379e+02 | 1.08e-01 | Optimal | Optimal | ipopt |
| HS84 | -5.239215e+06 | -5.280335e+06 | 7.79e-03 | Optimal | Optimal | ipopt |
| BIGGS5 | 5.655650e-03 | 1.088199e-19 | 5.66e-03 | Optimal | Optimal | ipopt |
| DNIEPER | 1.869902e+04 | 1.874401e+04 | 2.40e-03 | Acceptable | Optimal | ripopt |
| DECONVC | 3.734409e-03 | 2.569475e-03 | 1.16e-03 | Optimal | Optimal | ipopt |

---
*Generated by cutest_suite/compare.py*