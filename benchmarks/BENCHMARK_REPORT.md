# ripopt Benchmark Report

Generated: 2026-04-20 06:40:00

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **701/865** (81.0%) | **690/865** (79.8%) |
| Acceptable | 4 | 3 |
| Total solved (Optimal + Acceptable) | 705 (81.5%) | 693 (80.1%) |
| Solved exclusively | 41 | 29 |
| Both solved | 664 | |
| Matching objectives (< 0.01%) | 539/664 | |
| Acceptable at worse local min | 2 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 118 (98.3%) | 116 (96.7%) | 2 | 0 | 116 | 104/116 |
| CUTEst | 727 | 570 (78.4%) | 561 (77.2%) | 38 | 29 | 532 | 420/532 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| Grid | 4 | 4 (100.0%) | 4 (100.0%) | 0 | 0 | 4 | 3/4 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 116 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 139us | 2.4ms |
| Total time | 108.2ms | 346.0ms |
| Mean iterations | 78.3 | 13.3 |
| Median iterations | 14 | 10 |

- **Geometric mean speedup**: 13.7x
- **Median speedup**: 14.7x
- ripopt faster: 112/116 (97%)
- ripopt 10x+ faster: 82/116
- Ipopt faster: 4/116

## CUTEst Suite — Performance

On 532 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 183us | 2.7ms |
| Total time | 313.61s | 24.33s |
| Mean iterations | 130.7 | 35.1 |
| Median iterations | 18 | 12 |

- **Geometric mean speedup**: 6.5x
- **Median speedup**: 15.4x
- ripopt faster: 418/532 (79%)
- ripopt 10x+ faster: 301/532
- Ipopt faster: 114/532

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 183us | 2.4ms |
| Total time | 3.0ms | 74.6ms |
| Mean iterations | 22.6 | 19.2 |
| Median iterations | 8 | 7 |

- **Geometric mean speedup**: 20.3x
- **Median speedup**: 17.4x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 8/12
- Ipopt faster: 0/12

## Grid Suite — Performance

On 4 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 245.2ms | 10.2ms |
| Total time | 31.06s | 89.8ms |
| Mean iterations | 70.2 | 12.5 |
| Median iterations | 82 | 14 |

- **Geometric mean speedup**: 0.1x
- **Median speedup**: 0.2x
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
| NumericalError | 1 | 0 |
| RestorationFailed | 1 | 0 |

### CUTEst Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Crash(exit code 101) | 1 | 0 |
| ErrorInStepComputation | 0 | 2 |
| Infeasible | 0 | 10 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-10) | 0 | 123 |
| IpoptStatus(3) | 0 | 1 |
| IpoptStatus(4) | 0 | 2 |
| LocalInfeasibility | 59 | 0 |
| MaxIterations | 12 | 13 |
| NumericalError | 60 | 0 |
| RestorationFailed | 12 | 4 |
| Timeout | 13 | 10 |

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
| ACOPP30 | CUTEst | 72 | 142 | NumericalError | 5.768924e+02 |
| AIRPORT | CUTEst | 84 | 42 | Timeout | 4.795270e+04 |
| BATCH | CUTEst | 48 | 73 | NumericalError | 2.591803e+05 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | RestorationFailed | 7.862467e-01 |
| DISCS | CUTEst | 36 | 66 | NumericalError | 1.200007e+01 |
| EXPFITC | CUTEst | 5 | 502 | NumericalError | 2.330263e-02 |
| FLETCHER | CUTEst | 4 | 4 | NumericalError | 1.165685e+01 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HATFLDFL | CUTEst | 3 | 0 | NumericalError | 6.016849e-05 |
| HATFLDH | CUTEst | 4 | 7 | NumericalError | -2.450000e+01 |
| HET-Z | CUTEst | 2 | 1002 | Timeout | 1.000000e+00 |
| HIMMELBI | CUTEst | 100 | 12 | NumericalError | -1.735570e+03 |
| HS118 | CUTEst | 15 | 17 | NumericalError | 6.648204e+02 |
| HS83 | CUTEst | 5 | 3 | NumericalError | -3.066554e+04 |
| LAKES | CUTEst | 90 | 78 | RestorationFailed | 3.505251e+05 |
| MAKELA3 | CUTEst | 21 | 20 | NumericalError | -9.810659e-09 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| OET4 | CUTEst | 4 | 1002 | NumericalError | 4.295421e-03 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | NumericalError | -9.136404e-03 |
| SPANHYD | CUTEst | 97 | 33 | NumericalError | 2.397380e+02 |
| TAXR13322 | CUTEst | 72 | 1261 | NumericalError | -6.449419e+04 |
| TENBARS2 | CUTEst | 18 | 8 | NumericalError | 2.302548e+03 |
| TRO3X3 | CUTEst | 30 | 13 | NumericalError | 8.967478e+00 |
| VESUVIALS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |

## Wins (ripopt solves, Ipopt fails) — 41 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651B | CUTEst | 7 | 66 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651C | CUTEst | 7 | 56 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651D | CUTEst | 7 | 67 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 1.209250e-08 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA2NE | CUTEst | 5 | 16 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | CUTEst | 9 | 3 | ErrorInStepComputation | -8.274800e+02 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996921e+03 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.167666e+00 |
| LOGHAIRY | CUTEst | 2 | 0 | MaxIterations | 1.823216e-01 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -3.438801e+07 |
| MGH17 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| MGH17S | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7E | CUTEst | 8 | 0 | MaxIterations | 1.015395e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 6.547916e+00 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 5.616003e+00 |
| SPIRAL | CUTEst | 3 | 2 | MaxIterations | -1.132223e-09 |
| SSI | CUTEst | 3 | 0 | MaxIterations | 3.363463e-04 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.344196e+00 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 9.999767e-01 |

## Different Local Minima — 2 problems

ripopt converged (Acceptable) but to a different — usually worse — local
minimum than Ipopt found. Both solvers satisfied first-order KKT conditions
at their respective solutions. For nonconvex problems this is expected;
for convex problems it indicates the solver trajectory went astray.

| Problem | Suite | n | m | ripopt obj | Ipopt obj | Rel. error |
|---------|-------|---|---|------------|-----------|------------|
| HS98 | CUTEst | 6 | 4 | 4.368247e+00 | 3.135806e+00 | 28.2% |
| DEMBO7 | CUTEst | 16 | 20 | 1.784522e+02 | 1.747870e+02 | 2.1% |

## Acceptable (not Optimal) — 4 problems

These problems converged within relaxed tolerances but not strict tolerances.

| Problem | Suite | n | m | Ipopt status | ripopt obj | Ipopt obj |
|---------|-------|---|---|-------------|------------|-----------|
| ACOPR30 | CUTEst | 72 | 172 | Optimal | 5.815650e+02 | 5.768924e+02 |
| ALLINITC | CUTEst | 4 | 1 | Optimal | 3.048094e+01 | 3.049261e+01 |
| DEMBO7 | CUTEst | 16 | 20 | Optimal | 1.784522e+02 | 1.747870e+02 |
| HS98 | CUTEst | 6 | 4 | Optimal | 4.368247e+00 | 3.135806e+00 |

---
*Generated by benchmark_report.py*