# ripopt Benchmark Report

Generated: 2026-04-18 21:33:08

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **672/865** (77.7%) | **689/865** (79.7%) |
| Acceptable | 0 | 3 |
| Total solved (Optimal + Acceptable) | 672 (77.7%) | 692 (80.0%) |
| Solved exclusively | 40 | 60 |
| Both solved | 632 | |
| Matching objectives (< 0.01%) | 540/632 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 115 (95.8%) | 116 (96.7%) | 2 | 3 | 113 | 107/113 |
| CUTEst | 727 | 540 (74.3%) | 560 (77.0%) | 37 | 57 | 503 | 419/503 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 11/12 |
| Grid | 4 | 4 (100.0%) | 4 (100.0%) | 0 | 0 | 4 | 3/4 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 113 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 150us | 2.8ms |
| Total time | 399.3ms | 381.0ms |
| Mean iterations | 22.8 | 13.0 |
| Median iterations | 13 | 10 |

- **Geometric mean speedup**: 13.9x
- **Median speedup**: 16.0x
- ripopt faster: 110/113 (97%)
- ripopt 10x+ faster: 85/113
- Ipopt faster: 3/113

## CUTEst Suite — Performance

On 503 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 194us | 3.2ms |
| Total time | 89.48s | 11.52s |
| Mean iterations | 91.4 | 28.7 |
| Median iterations | 17 | 12 |

- **Geometric mean speedup**: 7.9x
- **Median speedup**: 18.3x
- ripopt faster: 407/503 (81%)
- ripopt 10x+ faster: 301/503
- Ipopt faster: 96/503

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 105us | 1.3ms |
| Total time | 2.1ms | 46.6ms |
| Mean iterations | 23.1 | 19.2 |
| Median iterations | 9 | 7 |

- **Geometric mean speedup**: 18.8x
- **Median speedup**: 21.2x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 7/12
- Ipopt faster: 0/12

## Grid Suite — Performance

On 4 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 102.0ms | 8.0ms |
| Total time | 460.4ms | 69.5ms |
| Mean iterations | 35.0 | 12.5 |
| Median iterations | 39 | 14 |

- **Geometric mean speedup**: 0.2x
- **Median speedup**: 0.5x
- ripopt faster: 1/4 (25%)
- ripopt 10x+ faster: 0/4
- Ipopt faster: 3/4

## Failure Analysis

### HS Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-11) | 0 | 1 |
| IpoptStatus(-199) | 0 | 1 |
| MaxIterations | 1 | 0 |
| NumericalError | 3 | 0 |
| RestorationFailed | 1 | 0 |

### CUTEst Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| ErrorInStepComputation | 0 | 2 |
| Infeasible | 0 | 10 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-10) | 0 | 123 |
| IpoptStatus(3) | 0 | 1 |
| IpoptStatus(4) | 0 | 2 |
| LocalInfeasibility | 56 | 0 |
| MaxIterations | 22 | 13 |
| NumericalError | 88 | 0 |
| RestorationFailed | 9 | 4 |
| Timeout | 12 | 11 |

### Electrolyte Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |

### CHO Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| MaxIter | 0 | 1 |
| MaxIterations | 1 | 0 |

## Regressions (Ipopt solves, ripopt fails)

| Problem | Suite | n | m | ripopt status | Ipopt obj |
|---------|-------|---|---|--------------|-----------|
| ACOPP30 | CUTEst | 72 | 142 | Timeout | 5.768924e+02 |
| ACOPR30 | CUTEst | 72 | 172 | NumericalError | 5.768924e+02 |
| ALLINITA | CUTEst | 4 | 4 | NumericalError | 3.329611e+01 |
| ALLINITC | CUTEst | 4 | 1 | MaxIterations | 3.049261e+01 |
| BATCH | CUTEst | 48 | 73 | NumericalError | 2.591803e+05 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CHWIRUT1LS | CUTEst | 3 | 0 | MaxIterations | 2.384477e+03 |
| CONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | RestorationFailed | 7.862467e-01 |
| DECONVBNE | CUTEst | 63 | 40 | NumericalError | 0.000000e+00 |
| DEGENLPA | CUTEst | 20 | 15 | NumericalError | 3.054881e+00 |
| DEMBO7 | CUTEst | 16 | 20 | NumericalError | 1.747870e+02 |
| DISCS | CUTEst | 36 | 66 | NumericalError | 1.200007e+01 |
| DUALC8 | CUTEst | 8 | 503 | NumericalError | 1.830936e+04 |
| ELATTAR | CUTEst | 7 | 102 | NumericalError | 7.420618e+01 |
| ERRINBAR | CUTEst | 18 | 9 | NumericalError | 2.804526e+01 |
| EXPFITC | CUTEst | 5 | 502 | NumericalError | 2.330263e-02 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HATFLDFL | CUTEst | 3 | 0 | NumericalError | 6.016849e-05 |
| HET-Z | CUTEst | 2 | 1002 | NumericalError | 1.000000e+00 |
| HIMMELBK | CUTEst | 24 | 14 | NumericalError | 5.181432e-02 |
| HS004 | HS | 2 | 0 | MaxIterations | 2.666667e+00 |
| HS013 | HS | 2 | 1 | NumericalError | 9.945785e-01 |
| HS108 | HS | 9 | 13 | NumericalError | -8.660254e-01 |
| HS13 | CUTEst | 2 | 1 | NumericalError | 9.945785e-01 |
| HS4 | CUTEst | 2 | 0 | MaxIterations | 2.666667e+00 |
| HS83 | CUTEst | 5 | 3 | NumericalError | -3.066554e+04 |
| HS84 | CUTEst | 5 | 3 | NumericalError | -5.280335e+06 |
| HS89 | CUTEst | 3 | 1 | MaxIterations | 1.362646e+00 |
| HS92 | CUTEst | 6 | 1 | MaxIterations | 1.362646e+00 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.486788e-17 |
| HYDCAR6LS | CUTEst | 29 | 0 | NumericalError | 3.040548e-20 |
| LAKES | CUTEst | 90 | 78 | NumericalError | 3.505251e+05 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| MAKELA3 | CUTEst | 21 | 20 | NumericalError | -9.810659e-09 |
| MCONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| MGH10LS | CUTEst | 3 | 0 | MaxIterations | 8.794586e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| NET1 | CUTEst | 48 | 57 | NumericalError | 9.411943e+05 |
| OET2 | CUTEst | 3 | 1002 | NumericalError | 8.715962e-02 |
| OET4 | CUTEst | 4 | 1002 | Timeout | 4.295421e-03 |
| OET5 | CUTEst | 5 | 1002 | MaxIterations | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | NumericalError | 4.465915e-05 |
| OPTCNTRL | CUTEst | 32 | 20 | NumericalError | 5.500000e+02 |
| PALMER5B | CUTEst | 9 | 0 | NumericalError | 9.752496e-03 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | MaxIterations | -9.136404e-03 |
| ROSZMAN1LS | CUTEst | 4 | 0 | NumericalError | 4.948485e-04 |
| SNAKE | CUTEst | 2 | 2 | NumericalError | -1.999999e-04 |
| SPANHYD | CUTEst | 97 | 33 | NumericalError | 2.397380e+02 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| SYNTHES2 | CUTEst | 11 | 14 | NumericalError | -5.544063e-01 |
| TRO3X3 | CUTEst | 30 | 13 | NumericalError | 8.967478e+00 |
| VESUVIALS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| VESUVIOLS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| WATER | CUTEst | 31 | 10 | NumericalError | 1.054938e+04 |

## Wins (ripopt solves, Ipopt fails) — 40 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| AVION2 | CUTEst | 49 | 15 | MaxIterations | 2.943684e+08 |
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651C | CUTEst | 7 | 56 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651D | CUTEst | 7 | 67 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651E | CUTEst | 7 | 64 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 2.569511e-03 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA2NE | CUTEst | 5 | 16 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | CUTEst | 9 | 3 | ErrorInStepComputation | -1.030758e+03 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996944e+03 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE5 | CUTEst | 2 | 4 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.126754e+00 |
| LOGHAIRY | CUTEst | 2 | 0 | MaxIterations | 1.823216e-01 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -2.163128e+12 |
| MGH17S | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 2.128087e+00 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 6.547916e+00 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| RAT43 | CUTEst | 4 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 5.980467e+00 |
| SPIRAL | CUTEst | 3 | 2 | MaxIterations | -1.132223e-09 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.281501e+00 |
| TAXR13322 | CUTEst | 72 | 1261 | Timeout | -7.462711e+03 |
| TRO6X2 | CUTEst | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000002e+00 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.003s | Optimal | 749 | 0.206s | 61.4x |
| SparseQP 1K | 500 | 500 | Optimal | 7 | 0.046s | Optimal | 6 | 0.004s | 0.1x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.003s | 1.4x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.003s | 0.4x |

ripopt: **4/4 solved** in 0.1s total
Ipopt: **4/4 solved** in 0.2s total

---
*Generated by benchmark_report.py*