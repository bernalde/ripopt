# HS Test Suite Validation Report

Comparison of ripopt vs cyipopt (C++ Ipopt) on the Hock-Schittkowski test suite.

## Executive Summary

- **Total problems**: 120
- **ripopt solved**: 119/120 (99.2%)
- **cyipopt solved**: 118/120 (98.3%)
- **Both solved**: 117/120 (97.5%)
- **Matching solutions** (rel obj diff < 1e-4): 107/117 (91.5%)

## Accuracy Statistics (where both solve)

| Metric | Objective Rel Diff | Solution Max Diff |
|--------|-------------------|-------------------|
| Mean   | 5.00e-02 | 1.82e+00 |
| Median | 3.18e-09 | 1.13e-08 |
| Max    | 1.00e+00 | 2.02e+02 |

## Category Breakdown

| Category | Total | ripopt | cyipopt | Both | Match |
|----------|-------|--------|---------|------|-------|
| constrained | 100 | 99 | 99 | 98 | 90 |
| unconstrained | 20 | 20 | 19 | 19 | 17 |

## Detailed Results

| TP# | n | m | ripopt | cyipopt | Obj Diff | r_iter | c_iter | r_time | c_time | Speedup | Status |
|-----|---|---|--------|---------|----------|--------|--------|--------|--------|---------|--------|
| 001 | 2 | 0 | Optimal | Optimal | 1.42e-14 | 26 | 29 | 15us | 4.8ms | 328.6x | PASS |
| 002 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 8 | 11 | 6us | 1.7ms | 300.0x | PASS |
| 003 | 2 | 0 | Optimal | Optimal | 9.99e-09 | 5 | 5 | 4us | 712us | 159.6x | PASS |
| 004 | 2 | 0 | Optimal | Optimal | 2.81e-08 | 5 | 5 | 4us | 689us | 155.9x | PASS |
| 005 | 2 | 0 | Optimal | Optimal | 4.64e-16 | 10 | 8 | 9us | 1.1ms | 133.1x | PASS |
| 006 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 2 | 6 | 5us | 841us | 180.2x | PASS |
| 007 | 2 | 1 | Optimal | Optimal | 8.63e-10 | 13 | 28 | 13us | 4.0ms | 304.9x | PASS |
| 009 | 2 | 1 | Optimal | Optimal | 7.11e-15 | 3 | 4 | 6us | 562us | 95.7x | PASS |
| 010 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 13 | 14us | 1.9ms | 133.8x | PASS |
| 011 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 7 | 6us | 1.0ms | 180.7x | PASS |
| 012 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 7 | 7 | 8us | 1.1ms | 132.8x | PASS |
| 013 | 2 | 1 | Acceptable | Optimal | 9.43e-04 | 35 | 49 | 25us | 7.2ms | 288.2x | MISMATCH |
| 014 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 6 | 6us | 935us | 152.6x | PASS |
| 015 | 2 | 2 | Optimal | Optimal | 7.99e-08 | 9 | 12 | 15us | 1.8ms | 116.7x | PASS |
| 016 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 12 | 11 | 11us | 1.8ms | 161.3x | MISMATCH |
| 017 | 2 | 2 | Optimal | Optimal | 2.58e-08 | 16 | 19 | 13us | 3.0ms | 232.9x | PASS |
| 018 | 2 | 2 | Optimal | Optimal | 7.10e-10 | 12 | 11 | 12us | 1.6ms | 137.9x | PASS |
| 019 | 2 | 2 | Optimal | Optimal | 3.18e-09 | 15 | 13 | 14us | 2.1ms | 148.1x | PASS |
| 020 | 2 | 3 | Optimal | Optimal | 6.97e-08 | 11 | 12 | 21us | 2.1ms | 98.4x | PASS |
| 021 | 2 | 1 | Optimal | Optimal | 1.37e-10 | 9 | 7 | 8us | 1.0ms | 135.3x | PASS |
| 022 | 2 | 2 | Optimal | Optimal | 3.21e-05 | 6 | 6 | 6us | 897us | 147.5x | PASS |
| 023 | 2 | 5 | Optimal | Optimal | 8.69e-01 | 25 | 10 | 42us | 1.5ms | 35.1x | MISMATCH |
| 024 | 2 | 3 | Optimal | Optimal | 1.51e-08 | 12 | 15 | 24us | 2.6ms | 108.6x | PASS |
| 026 | 3 | 1 | Optimal | Optimal | 8.88e-16 | 20 | 26 | 17us | 2.9ms | 174.7x | PASS |
| 027 | 3 | 1 | Optimal | Optimal | 3.88e-11 | 27 | 22 | 33us | 2.7ms | 80.8x | PASS |
| 028 | 3 | 1 | Optimal | Optimal | 0.00e+00 | 1 | 2 | 4us | 279us | 79.8x | PASS |
| 029 | 3 | 1 | Optimal | Optimal | 3.12e-10 | 13 | 8 | 17us | 1.2ms | 71.8x | PASS |
| 030 | 3 | 1 | Optimal | Optimal | 1.44e-06 | 10 | 8 | 9us | 1.2ms | 131.9x | PASS |
| 031 | 3 | 1 | Optimal | Optimal | 1.02e-08 | 11 | 7 | 11us | 1.1ms | 95.5x | PASS |
| 032 | 3 | 2 | Optimal | Optimal | 4.56e-08 | 15 | 16 | 14us | 2.4ms | 174.8x | PASS |
| 033 | 3 | 2 | Optimal | Optimal | 2.84e-08 | 11 | 10 | 17us | 1.5ms | 91.1x | PASS |
| 034 | 3 | 2 | Optimal | Optimal | 1.14e-08 | 11 | 8 | 15us | 1.2ms | 80.2x | PASS |
| 035 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 8 | 9us | 1.2ms | 139.9x | PASS |
| 036 | 3 | 1 | Optimal | Optimal | 6.31e-09 | 13 | 12 | 16us | 2.0ms | 126.0x | PASS |
| 037 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 9 | 12 | 14us | 2.0ms | 144.8x | PASS |
| 038 | 4 | 0 | Optimal | Optimal | 1.00e+00 | 9 | 40 | 9us | 6.4ms | 685.2x | MISMATCH |
| 039 | 4 | 2 | Optimal | Optimal | 6.57e-09 | 23 | 14 | 33us | 1.6ms | 49.0x | PASS |
| 040 | 4 | 3 | Optimal | Optimal | 7.91e-11 | 4 | 4 | 8us | 495us | 66.1x | PASS |
| 041 | 4 | 1 | Optimal | Optimal | 1.61e-09 | 11 | 8 | 12us | 1.2ms | 100.5x | PASS |
| 042 | 4 | 2 | Optimal | Optimal | 4.73e-12 | 5 | 5 | 7us | 593us | 86.3x | PASS |
| 043 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 8 | 9 | 13us | 1.4ms | 110.7x | PASS |
| 044 | 4 | 6 | Optimal | Optimal | 1.01e-08 | 15 | 25 | 35us | 4.1ms | 117.6x | PASS |
| 045 | 5 | 0 | Optimal | Optimal | 6.52e-08 | 73 | 12 | 99us | 1.9ms | 19.5x | PASS |
| 046 | 5 | 2 | Optimal | Optimal | 2.11e-02 | 25 | 20 | 26us | 2.7ms | 102.9x | MISMATCH |
| 047 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 13 | 16 | 26us | 2.1ms | 81.2x | PASS |
| 048 | 5 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 2 | 4us | 283us | 67.3x | PASS |
| 049 | 5 | 2 | Optimal | Optimal | 8.51e-12 | 20 | 20 | 19us | 2.6ms | 140.0x | PASS |
| 050 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 9 | 11us | 1.1ms | 96.0x | PASS |
| 051 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 2 | 4us | 278us | 62.4x | PASS |
| 052 | 5 | 3 | Optimal | Optimal | 5.00e-16 | 1 | 2 | 4us | 278us | 64.1x | PASS |
| 053 | 5 | 3 | Optimal | Optimal | 6.30e-13 | 8 | 7 | 15us | 990us | 66.9x | PASS |
| 056 | 7 | 4 | Optimal | Optimal | 9.39e-49 | 2 | 4 | 10us | 735us | 75.1x | PASS |
| 058 | 2 | 3 | Optimal | Optimal | 1.25e-08 | 11 | 12 | 13us | 1.9ms | 146.8x | PASS |
| 060 | 3 | 1 | Optimal | Optimal | 1.20e-13 | 8 | 7 | 8us | 1.1ms | 124.3x | PASS |
| 061 | 3 | 2 | Optimal | Optimal | 1.19e-15 | 9 | 11 | 11us | 1.2ms | 108.3x | PASS |
| 063 | 3 | 2 | Optimal | Optimal | 1.18e-16 | 9 | 6 | 12us | 894us | 75.3x | PASS |
| 064 | 3 | 1 | Optimal | Optimal | 1.63e-09 | 16 | 17 | 14us | 2.7ms | 192.6x | PASS |
| 065 | 3 | 1 | Optimal | Optimal | 6.48e-09 | 11 | 17 | 13us | 2.9ms | 220.4x | PASS |
| 066 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 10 | 11 | 12us | 1.6ms | 131.7x | PASS |
| 071 | 4 | 2 | Optimal | Optimal | 4.54e-09 | 10 | 13 | 12us | 1.9ms | 161.2x | PASS |
| 072 | 4 | 2 | Optimal | Optimal | 6.70e-07 | 28 | 17 | 34us | 2.5ms | 75.6x | PASS |
| 076 | 4 | 3 | Optimal | Optimal | 4.73e-09 | 11 | 8 | 13us | 1.2ms | 91.4x | PASS |
| 077 | 5 | 2 | Optimal | Optimal | 2.61e-11 | 9 | 12 | 13us | 1.5ms | 117.0x | PASS |
| 078 | 5 | 3 | Optimal | Optimal | 3.72e-10 | 4 | 5 | 8us | 630us | 80.4x | PASS |
| 079 | 5 | 3 | Optimal | Optimal | 2.60e-10 | 4 | 5 | 8us | 629us | 83.8x | PASS |
| 080 | 5 | 3 | Optimal | Optimal | 5.45e-13 | 9 | 6 | 15us | 1.0ms | 70.8x | PASS |
| 081 | 5 | 3 | Optimal | Optimal | 9.46e-01 | 73 | 70 | 174us | 13.2ms | 75.9x | MISMATCH |
| 106 | 8 | 6 | Optimal | Optimal | 1.27e-08 | 17 | 19 | 51us | 2.9ms | 56.1x | PASS |
| 108 | 9 | 13 | Optimal | Optimal | 5.70e-09 | 19 | 24 | 124us | 4.5ms | 36.5x | PASS |
| 113 | 10 | 8 | Optimal | Optimal | 1.63e-09 | 13 | 10 | 57us | 1.7ms | 30.5x | PASS |
| 114 | 10 | 11 | Optimal | Optimal | 1.06e-07 | 30 | 14 | 147us | 2.4ms | 16.4x | PASS |
| 116 | 13 | 15 | Optimal | Optimal | 3.69e-07 | 56 | 18 | 497us | 3.1ms | 6.2x | PASS |
| 201 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1 | 2 | 3us | 215us | 68.9x | PASS |
| 206 | 2 | 0 | Optimal | Optimal | 3.77e-15 | 4 | 5 | 4us | 536us | 126.2x | PASS |
| 211 | 2 | 0 | Optimal | Optimal | 9.77e-15 | 27 | 28 | 13us | 3.7ms | 281.0x | PASS |
| 212 | 2 | 0 | Optimal | Optimal | 5.24e-23 | 8 | 12 | 6us | 1.7ms | 261.5x | PASS |
| 213 | 2 | 0 | Acceptable | Optimal | 6.18e-06 | 44 | 23 | 38us | 2.9ms | 75.5x | PASS |
| 214 | 2 | 0 | Optimal | IpoptStatus( | N/A | 42 | 31 | 50us | 7.3ms | 147.9x | MISMATCH |
| 215 | 2 | 1 | Optimal | Optimal | 1.36e-08 | 9 | 14 | 8us | 2.0ms | 253.4x | PASS |
| 216 | 2 | 1 | Optimal | Optimal | 1.42e-13 | 9 | 8 | 10us | 1.2ms | 125.5x | PASS |
| 217 | 2 | 2 | Optimal | Optimal | 1.22e-08 | 10 | 9 | 9us | 1.3ms | 149.9x | PASS |
| 218 | 2 | 1 | Optimal | Optimal | 1.00e-08 | 19 | 7 | 14us | 1.0ms | 75.0x | PASS |
| 219 | 4 | 2 | Optimal | Optimal | 1.50e-13 | 26 | 51 | 27us | 6.1ms | 226.1x | PASS |
| 220 | 2 | 1 | Optimal | Optimal | 1.00e-08 | 4 | 4 | 5us | 664us | 124.5x | PASS |
| 221 | 2 | 1 | Optimal | Optimal | 9.99e-09 | 8 | 38 | 7us | 5.2ms | 720.9x | PASS |
| 223 | 2 | 2 | Optimal | Infeasible | N/A | 4 | 1582 | 10us | 639.5ms | 64759.5x | MISMATCH |
| 224 | 2 | 4 | Optimal | Optimal | 8.28e-10 | 15 | 10 | 15us | 1.5ms | 97.5x | PASS |
| 225 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 13 | 11 | 19us | 1.6ms | 84.2x | MISMATCH |
| 226 | 2 | 2 | Optimal | Optimal | 4.99e-09 | 13 | 12 | 14us | 1.8ms | 134.7x | PASS |
| 227 | 2 | 2 | Optimal | Optimal | 2.00e-08 | 6 | 11 | 9us | 1.7ms | 192.2x | PASS |
| 228 | 2 | 2 | Optimal | Optimal | 2.16e-09 | 5 | 8 | 8us | 1.2ms | 159.3x | PASS |
| 229 | 2 | 0 | Optimal | Optimal | 4.19e-03 | 16 | 22 | 11us | 3.7ms | 352.6x | MISMATCH |
| 230 | 2 | 2 | Optimal | Optimal | 9.67e-09 | 2 | 8 | 4us | 1.2ms | 335.5x | PASS |
| 232 | 2 | 3 | Optimal | Optimal | 1.82e-08 | 9 | 9 | 14us | 1.5ms | 109.9x | PASS |
| 234 | 2 | 1 | Optimal | Optimal | 2.68e-01 | 9 | 19 | 9us | 2.9ms | 327.4x | MISMATCH |
| 235 | 3 | 1 | Optimal | Optimal | 3.62e-14 | 29 | 14 | 32us | 1.8ms | 56.1x | PASS |
| 240 | 3 | 0 | Optimal | Optimal | 3.23e-29 | 1 | 2 | 3us | 221us | 66.4x | PASS |
| 248 | 3 | 2 | Optimal | Optimal | 5.46e-08 | 10 | 15 | 14us | 2.4ms | 172.7x | PASS |
| 249 | 3 | 1 | Optimal | Optimal | 1.00e-08 | 22 | 8 | 15us | 1.2ms | 77.3x | PASS |
| 250 | 3 | 2 | Optimal | Optimal | 6.32e-09 | 13 | 15 | 18us | 2.5ms | 137.4x | PASS |
| 251 | 3 | 1 | Optimal | Optimal | 4.17e-10 | 9 | 11 | 12us | 1.8ms | 151.7x | PASS |
| 252 | 3 | 1 | Optimal | Optimal | 4.00e-10 | 21 | 16 | 15us | 2.3ms | 157.5x | PASS |
| 254 | 3 | 2 | Optimal | Optimal | 1.37e-11 | 17 | 26 | 14us | 3.7ms | 264.0x | PASS |
| 255 | 4 | 0 | Optimal | Optimal | 1.99e-08 | 13 | 14 | 13us | 2.3ms | 172.7x | PASS |
| 256 | 4 | 0 | Optimal | Optimal | 5.35e-12 | 20 | 20 | 12us | 2.2ms | 176.4x | PASS |
| 257 | 4 | 0 | Optimal | Optimal | 2.31e-14 | 12 | 8 | 12us | 1.4ms | 124.0x | PASS |
| 258 | 4 | 0 | Optimal | Optimal | 5.35e-14 | 42 | 41 | 25us | 5.6ms | 221.9x | PASS |
| 259 | 4 | 0 | Optimal | Optimal | 2.02e-14 | 11 | 10 | 9us | 1.8ms | 196.9x | PASS |
| 262 | 4 | 4 | Optimal | Optimal | 2.39e-09 | 16 | 8 | 29us | 1.2ms | 42.5x | PASS |
| 263 | 4 | 4 | Optimal | Acceptable | 1.29e-08 | 21 | 75 | 23us | 11.5ms | 490.4x | PASS |
| 264 | 4 | 3 | Optimal | Optimal | 6.66e-10 | 8 | 14 | 14us | 2.2ms | 159.7x | PASS |
| 270 | 5 | 1 | Optimal | Optimal | 2.12e-12 | 15 | 12 | 19us | 2.3ms | 119.5x | PASS |
| 325 | 2 | 3 | Optimal | Optimal | 1.91e-09 | 12 | 16 | 22us | 2.8ms | 127.9x | PASS |
| 335 | 3 | 2 | Optimal | Optimal | 4.32e-13 | 26 | 26 | 22us | 2.9ms | 128.5x | PASS |
| 338 | 3 | 2 | Optimal | Optimal | 4.85e-16 | 26 | 46 | 34us | 6.0ms | 178.8x | PASS |
| 339 | 3 | 1 | Optimal | Optimal | 3.35e-09 | 10 | 8 | 12us | 1.2ms | 104.3x | PASS |
| 344 | 3 | 1 | Optimal | Optimal | 6.22e-14 | 7 | 8 | 7us | 962us | 132.6x | PASS |
| 354 | 4 | 1 | Optimal | Optimal | 4.87e-09 | 9 | 11 | 10us | 1.7ms | 172.4x | PASS |
| 374 | 10 | 35 | MaxIteration | Optimal | N/A | 2999 | 103 | 224.9ms | 171.0ms | 0.8x | ripopt_FAIL |
| 376 | 10 | 15 | Optimal | Optimal | 9.62e-01 | 59 | 24 | 1.0ms | 5.2ms | 5.2x | MISMATCH |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | cyipopt |
|--------|--------|---------|
| Mean   | 14.5 | 14.5 |
| Median | 11 | 11 |
| Max    | 73 | 75 |
| Total  | 1691 | 1693 |

- ripopt uses fewer iterations: 60/117 problems
- cyipopt uses fewer iterations: 46/117 problems
- Same iteration count: 11/117 problems

### Timing Comparison

| Metric | ripopt | cyipopt |
|--------|--------|---------|
| Mean   | 31us | 2.2ms |
| Median | 13us | 1.7ms |
| Max    | 1.0ms | 13.2ms |
| Total  | 3.7ms | 254.0ms |

- Geometric mean speedup (cyipopt_time/ripopt_time): **116.38x**
  - \>1 means ripopt is faster, <1 means cyipopt is faster
- ripopt faster: 117/117 problems
- cyipopt faster: 0/117 problems
- Overall speedup (total time): 69.10x


## Failure Analysis

### Problems where only ripopt fails (1)

| TP# | n | m | ripopt status | cyipopt obj |
|-----|---|---|---------------|-------------|
| 374 | 10 | 35 | MaxIterations | 2.332635e-01 |

### Problems where only cyipopt fails (2)

| TP# | n | m | cyipopt status | ripopt obj |
|-----|---|---|----------------|------------|
| 214 | 2 | 0 | IpoptStatus(-13) | 0.000000e+00 |
| 223 | 2 | 2 | Infeasible | N/A |

### Objective mismatches (both solve but differ > 1e-4) (10)

| TP# | ripopt obj | cyipopt obj | Rel Diff |
|-----|-----------|-------------|----------|
| 013 | 9.936354e-01 | 9.945785e-01 | 9.43e-04 |
| 016 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 |
| 023 | 1.522111e+01 | 2.000000e+00 | 8.69e-01 |
| 038 | 7.875146e+00 | -7.105427e-15 | 1.00e+00 |
| 046 | 2.220446e-16 | 2.108893e-02 | 2.11e-02 |
| 081 | 9.999998e-01 | 5.394985e-02 | 9.46e-01 |
| 225 | 9.472136e+00 | 2.000000e+00 | 7.89e-01 |
| 229 | 4.190708e-03 | 1.421085e-14 | 4.19e-03 |
| 234 | -5.323739e-01 | -8.000000e-01 | 2.68e-01 |
| 376 | -6.060761e+01 | -1.614959e+03 | 9.62e-01 |

---
*Generated by hs_suite/compare.py*