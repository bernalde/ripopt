# HS Test Suite Validation Report

Comparison of ripopt vs cyipopt (C++ Ipopt) on the Hock-Schittkowski test suite.

## Executive Summary

- **Total problems**: 120
- **ripopt solved**: 116/120 (96.7%)
- **cyipopt solved**: 116/120 (96.7%)
- **Both solved**: 114/120 (95.0%)
- **Matching solutions** (rel obj diff < 1e-4): 106/114 (93.0%)

## Accuracy Statistics (where both solve)

| Metric | Objective Rel Diff | Solution Max Diff |
|--------|-------------------|-------------------|
| Mean   | 3.46e-02 | 8.98e-02 |
| Median | 3.18e-09 | 1.05e-08 |
| Max    | 1.00e+00 | 2.62e+00 |

## Category Breakdown

| Category | Total | ripopt | cyipopt | Both | Match |
|----------|-------|--------|---------|------|-------|
| constrained | 100 | 96 | 97 | 95 | 89 |
| unconstrained | 20 | 20 | 19 | 19 | 17 |

## Detailed Results

| TP# | n | m | ripopt | cyipopt | Obj Diff | x Diff | Status |
|-----|---|---|--------|---------|----------|--------|--------|
| 001 | 2 | 0 | Optimal | Optimal | 1.42e-14 | 1.77e-10 | PASS |
| 002 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 1.49e-08 | PASS |
| 003 | 2 | 0 | Optimal | Optimal | 9.99e-09 | 9.99e-09 | PASS |
| 004 | 2 | 0 | Optimal | Optimal | 2.81e-08 | 3.51e-08 | PASS |
| 005 | 2 | 0 | Optimal | Optimal | 4.64e-16 | 6.14e-09 | PASS |
| 006 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 2.22e-16 | PASS |
| 007 | 2 | 1 | Optimal | Optimal | 8.63e-10 | 1.49e-09 | PASS |
| 009 | 2 | 1 | Optimal | Optimal | 7.11e-15 | 5.68e-14 | PASS |
| 010 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 4.99e-09 | PASS |
| 011 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 3.99e-09 | PASS |
| 012 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 4.73e-10 | PASS |
| 013 | 2 | 1 | Acceptable | Optimal | 9.43e-04 | 4.73e-04 | MISMATCH |
| 014 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 7.56e-09 | PASS |
| 015 | 2 | 2 | Optimal | Optimal | 7.99e-08 | 5.99e-08 | PASS |
| 016 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 1.00e+00 | MISMATCH |
| 017 | 2 | 2 | Optimal | Optimal | 2.58e-08 | 2.13e-05 | PASS |
| 018 | 2 | 2 | Optimal | Optimal | 7.10e-10 | 4.14e-09 | PASS |
| 019 | 2 | 2 | Optimal | Optimal | 3.18e-09 | 1.96e-08 | PASS |
| 020 | 2 | 3 | Optimal | Optimal | 6.97e-08 | 1.15e-08 | PASS |
| 021 | 2 | 1 | Optimal | Optimal | 1.37e-10 | 3.43e-07 | PASS |
| 022 | 2 | 2 | Optimal | Optimal | 3.21e-05 | 3.20e-05 | PASS |
| 023 | 2 | 5 | Optimal | Optimal | 8.69e-01 | 1.76e+00 | MISMATCH |
| 024 | 2 | 3 | Optimal | Optimal | 1.51e-08 | 8.74e-09 | PASS |
| 026 | 3 | 1 | Optimal | Optimal | 8.88e-16 | 9.67e-06 | PASS |
| 027 | 3 | 1 | Optimal | Optimal | 3.88e-11 | 7.80e-11 | PASS |
| 028 | 3 | 1 | Optimal | Optimal | 0.00e+00 | 1.78e-15 | PASS |
| 029 | 3 | 1 | Optimal | Optimal | 3.12e-10 | 4.15e-10 | PASS |
| 030 | 3 | 1 | Optimal | Optimal | 1.44e-06 | 1.08e-03 | PASS |
| 031 | 3 | 1 | Optimal | Optimal | 1.02e-08 | 1.13e-08 | PASS |
| 032 | 3 | 2 | Optimal | Optimal | 4.56e-08 | 1.94e-05 | PASS |
| 033 | 3 | 2 | Optimal | Optimal | 2.84e-08 | 1.19e-08 | PASS |
| 034 | 3 | 2 | Optimal | Optimal | 1.14e-08 | 1.23e-08 | PASS |
| 035 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 8.59e-09 | PASS |
| 036 | 3 | 1 | Optimal | Optimal | 6.31e-09 | 2.05e-07 | PASS |
| 037 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 3.48e-09 | PASS |
| 038 | 4 | 0 | Optimal | Optimal | 1.00e+00 | 2.02e+00 | MISMATCH |
| 039 | 4 | 2 | Optimal | Optimal | 6.57e-09 | 1.75e-08 | PASS |
| 040 | 4 | 3 | Optimal | Optimal | 7.91e-11 | 3.98e-10 | PASS |
| 041 | 4 | 1 | Optimal | Optimal | 1.61e-09 | 3.72e-06 | PASS |
| 042 | 4 | 2 | Optimal | Optimal | 4.73e-12 | 3.54e-11 | PASS |
| 043 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 5.26e-07 | PASS |
| 044 | 4 | 6 | Optimal | Optimal | 1.01e-08 | 1.01e-08 | PASS |
| 045 | 5 | 0 | Optimal | Optimal | 6.52e-08 | 6.54e-08 | PASS |
| 046 | 5 | 2 | Optimal | Optimal | 2.11e-02 | 1.58e+00 | MISMATCH |
| 047 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4.92e-10 | PASS |
| 048 | 5 | 2 | Optimal | Optimal | 0.00e+00 | 8.88e-16 | PASS |
| 049 | 5 | 2 | Optimal | Optimal | 8.51e-12 | 1.20e-03 | PASS |
| 050 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 2.22e-16 | PASS |
| 051 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 0.00e+00 | PASS |
| 052 | 5 | 3 | Optimal | Optimal | 5.00e-16 | 4.44e-16 | PASS |
| 053 | 5 | 3 | Optimal | Optimal | 6.30e-13 | 8.57e-07 | PASS |
| 056 | 7 | 4 | Optimal | Optimal | 9.39e-49 | 1.13e-16 | PASS |
| 058 | 2 | 3 | Optimal | Optimal | 1.25e-08 | 2.33e-08 | PASS |
| 060 | 3 | 1 | Optimal | Optimal | 1.20e-13 | 9.01e-10 | PASS |
| 061 | 3 | 2 | Optimal | Optimal | 1.19e-15 | 2.40e-14 | PASS |
| 063 | 3 | 2 | Optimal | Optimal | 1.18e-16 | 2.46e-09 | PASS |
| 064 | 3 | 1 | Optimal | Optimal | 1.63e-09 | 1.04e-06 | PASS |
| 065 | 3 | 1 | Optimal | Optimal | 6.48e-09 | 2.97e-08 | PASS |
| 066 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 2.81e-07 | PASS |
| 071 | 4 | 2 | Optimal | Optimal | 4.54e-09 | 7.98e-07 | PASS |
| 072 | 4 | 2 | Optimal | Optimal | 6.70e-07 | 2.34e-04 | PASS |
| 076 | 4 | 3 | Optimal | Optimal | 4.73e-09 | 1.05e-08 | PASS |
| 077 | 5 | 2 | Optimal | Optimal | 2.61e-11 | 4.58e-09 | PASS |
| 078 | 5 | 3 | Optimal | Optimal | 3.72e-10 | 2.87e-10 | PASS |
| 079 | 5 | 3 | Optimal | Optimal | 2.60e-10 | 4.65e-09 | PASS |
| 080 | 5 | 3 | Optimal | Optimal | 5.45e-13 | 2.74e-07 | PASS |
| 081 | 5 | 3 | RestorationF | Optimal | N/A | N/A | ripopt_FAIL |
| 106 | 8 | 6 | Optimal | Optimal | 1.27e-08 | 9.30e-05 | PASS |
| 108 | 9 | 13 | Optimal | Optimal | 5.70e-09 | 6.15e-02 | PASS |
| 113 | 10 | 8 | Optimal | Optimal | 1.63e-09 | 1.39e-09 | PASS |
| 114 | 10 | 11 | Optimal | Optimal | 1.06e-07 | 1.10e-03 | PASS |
| 116 | 13 | 15 | RestorationF | IpoptStatus( | N/A | N/A | BOTH_FAIL |
| 201 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 0.00e+00 | PASS |
| 206 | 2 | 0 | Optimal | Optimal | 3.77e-15 | 2.22e-16 | PASS |
| 211 | 2 | 0 | Optimal | Optimal | 9.77e-15 | 1.93e-14 | PASS |
| 212 | 2 | 0 | Optimal | Optimal | 5.24e-23 | 1.85e-12 | PASS |
| 213 | 2 | 0 | Acceptable | Optimal | 6.18e-06 | 5.53e-02 | PASS |
| 214 | 2 | 0 | Optimal | IpoptStatus( | N/A | N/A | MISMATCH |
| 215 | 2 | 1 | Optimal | Optimal | 1.36e-08 | 5.43e-05 | PASS |
| 216 | 2 | 1 | Optimal | Optimal | 1.42e-13 | 8.42e-12 | PASS |
| 217 | 2 | 2 | Optimal | Optimal | 1.22e-08 | 1.62e-08 | PASS |
| 218 | 2 | 1 | Optimal | Optimal | 1.00e-08 | 9.31e-05 | PASS |
| 219 | 4 | 2 | Optimal | Optimal | 1.50e-13 | 4.01e-13 | PASS |
| 220 | 2 | 1 | Optimal | Optimal | 1.00e-08 | 1.00e-08 | PASS |
| 221 | 2 | 1 | Optimal | Optimal | 9.99e-09 | 9.99e-09 | PASS |
| 223 | 2 | 2 | Optimal | Infeasible | N/A | N/A | MISMATCH |
| 224 | 2 | 4 | Optimal | Optimal | 8.28e-10 | 5.56e-09 | PASS |
| 225 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 2.62e+00 | MISMATCH |
| 226 | 2 | 2 | Optimal | Optimal | 4.99e-09 | 3.53e-09 | PASS |
| 227 | 2 | 2 | Optimal | Optimal | 2.00e-08 | 9.99e-09 | PASS |
| 228 | 2 | 2 | Optimal | Optimal | 2.16e-09 | 1.71e-05 | PASS |
| 229 | 2 | 0 | Optimal | Optimal | 4.19e-03 | 1.25e-01 | MISMATCH |
| 230 | 2 | 2 | Optimal | Optimal | 9.67e-09 | 9.67e-09 | PASS |
| 232 | 2 | 3 | Optimal | Optimal | 1.82e-08 | 1.05e-08 | PASS |
| 234 | 2 | 1 | Optimal | Optimal | 2.68e-01 | 7.29e-01 | MISMATCH |
| 235 | 3 | 1 | Optimal | Optimal | 3.62e-14 | 3.88e-11 | PASS |
| 240 | 3 | 0 | Optimal | Optimal | 3.23e-29 | 3.00e-15 | PASS |
| 248 | 3 | 2 | Optimal | Optimal | 5.46e-08 | 7.14e-08 | PASS |
| 249 | 3 | 1 | Optimal | Optimal | 1.00e-08 | 3.06e-08 | PASS |
| 250 | 3 | 2 | Optimal | Optimal | 6.32e-09 | 2.05e-07 | PASS |
| 251 | 3 | 1 | Optimal | Optimal | 4.17e-10 | 3.48e-09 | PASS |
| 252 | 3 | 1 | Optimal | Optimal | 4.00e-10 | 5.36e-05 | PASS |
| 254 | 3 | 2 | Optimal | Optimal | 1.37e-11 | 9.68e-05 | PASS |
| 255 | 4 | 0 | Optimal | Optimal | 1.99e-08 | 1.02e-07 | PASS |
| 256 | 4 | 0 | Optimal | Optimal | 5.35e-12 | 3.58e-04 | PASS |
| 257 | 4 | 0 | Optimal | Optimal | 2.31e-14 | 8.93e-09 | PASS |
| 258 | 4 | 0 | Optimal | Optimal | 5.35e-14 | 5.91e-14 | PASS |
| 259 | 4 | 0 | Optimal | Optimal | 2.02e-14 | 2.93e-11 | PASS |
| 262 | 4 | 4 | Optimal | Optimal | 4.79e-09 | 5.05e-08 | PASS |
| 263 | 4 | 4 | Optimal | Acceptable | 1.29e-08 | 5.02e-05 | PASS |
| 264 | 4 | 3 | Optimal | Optimal | 6.66e-10 | 2.44e-05 | PASS |
| 270 | 5 | 1 | Optimal | Optimal | 2.12e-12 | 2.82e-01 | PASS |
| 325 | 2 | 3 | Optimal | Optimal | 7.22e-10 | 5.86e-10 | PASS |
| 335 | 3 | 2 | Optimal | Optimal | 4.32e-13 | 1.09e-12 | PASS |
| 338 | 3 | 2 | Optimal | Optimal | 4.85e-16 | 1.11e-15 | PASS |
| 339 | 3 | 1 | Optimal | Optimal | 3.35e-09 | 1.24e-07 | PASS |
| 344 | 3 | 1 | Optimal | Optimal | 6.22e-14 | 4.11e-12 | PASS |
| 354 | 4 | 1 | Optimal | Optimal | 4.87e-09 | 6.58e-09 | PASS |
| 374 | 10 | 35 | MaxIteration | Error: Faile | N/A | N/A | BOTH_FAIL |
| 376 | 10 | 15 | NumericalErr | Optimal | N/A | N/A | ripopt_FAIL |

## Performance Comparison (where both solve)

| Metric | ripopt iters | cyipopt iters |
|--------|-------------|---------------|
| Mean   | 13.4 | 13.9 |
| Median | 11 | 11 |
| Max    | 73 | 75 |

## Failure Analysis

### Problems where only ripopt fails (2)

| TP# | n | m | ripopt status | cyipopt obj |
|-----|---|---|---------------|-------------|
| 081 | 5 | 3 | RestorationFailed | 5.394985e-02 |
| 376 | 10 | 15 | NumericalError | -1.614959e+03 |

### Problems where only cyipopt fails (2)

| TP# | n | m | cyipopt status | ripopt obj |
|-----|---|---|----------------|------------|
| 214 | 2 | 0 | IpoptStatus(-13) | 0.000000e+00 |
| 223 | 2 | 2 | Infeasible | N/A |

### Problems where both fail (2)

| TP# | n | m | ripopt status | cyipopt status |
|-----|---|---|---------------|----------------|
| 116 | 13 | 15 | RestorationFailed | IpoptStatus(-11) |
| 374 | 10 | 35 | MaxIterations | Error: Failed to create NLP problem. Make sure inputs are okay! |

### Objective mismatches (both solve but differ > 1e-4) (8)

| TP# | ripopt obj | cyipopt obj | Rel Diff |
|-----|-----------|-------------|----------|
| 013 | 9.936354e-01 | 9.945785e-01 | 9.43e-04 |
| 016 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 |
| 023 | 1.522111e+01 | 2.000000e+00 | 8.69e-01 |
| 038 | 7.875146e+00 | -7.105427e-15 | 1.00e+00 |
| 046 | 2.220446e-16 | 2.108893e-02 | 2.11e-02 |
| 225 | 9.472136e+00 | 2.000000e+00 | 7.89e-01 |
| 229 | 4.190708e-03 | 1.421085e-14 | 4.19e-03 |
| 234 | -5.323739e-01 | -8.000000e-01 | 2.68e-01 |

---
*Generated by hs_suite/compare.py*