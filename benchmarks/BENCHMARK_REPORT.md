# ripopt Benchmark Report

Generated: 2026-04-10 16:56:46

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **701/865** (81.0%) | **688/865** (79.5%) |
| Acceptable | 0 | 5 |
| Total solved (Optimal + Acceptable) | 701 (81.0%) | 693 (80.1%) |
| Solved exclusively | 46 | 38 |
| Both solved | 655 | |
| Matching objectives (< 0.01%) | 550/655 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 116 (96.7%) | 116 (96.7%) | 2 | 2 | 114 | 106/114 |
| CUTEst | 727 | 569 (78.3%) | 561 (77.2%) | 43 | 35 | 526 | 430/526 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| Grid | 4 | 3 (75.0%) | 4 (100.0%) | 0 | 1 | 3 | 2/3 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 114 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 94us | 1.8ms |
| Total time | 23.1ms | 269.4ms |
| Mean iterations | 15.1 | 13.0 |
| Median iterations | 13 | 10 |

- **Geometric mean speedup**: 17.7x
- **Median speedup**: 17.8x
- ripopt faster: 113/114 (99%)
- ripopt 10x+ faster: 95/114
- Ipopt faster: 1/114

## CUTEst Suite — Performance

On 526 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 101us | 2.3ms |
| Total time | 114.15s | 8.58s |
| Mean iterations | 57.6 | 38.1 |
| Median iterations | 16 | 12 |

- **Geometric mean speedup**: 11.1x
- **Median speedup**: 23.1x
- ripopt faster: 455/526 (87%)
- ripopt 10x+ faster: 341/526
- Ipopt faster: 71/526

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 92us | 1.3ms |
| Total time | 2.0ms | 46.0ms |
| Mean iterations | 23.5 | 19.2 |
| Median iterations | 11 | 7 |

- **Geometric mean speedup**: 20.3x
- **Median speedup**: 19.6x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 7/12
- Ipopt faster: 0/12

## Grid Suite — Performance

On 3 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 6.8ms | 3.6ms |
| Total time | 16.7ms | 15.2ms |
| Mean iterations | 191.3 | 12.0 |
| Median iterations | 28 | 11 |

- **Geometric mean speedup**: 1.0x
- **Median speedup**: 1.3x
- ripopt faster: 2/3 (67%)
- ripopt 10x+ faster: 0/3
- Ipopt faster: 1/3

## Failure Analysis

### HS Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-11) | 0 | 1 |
| IpoptStatus(-199) | 0 | 1 |
| NumericalError | 3 | 0 |
| RestorationFailed | 1 | 0 |

### CUTEst Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| ErrorInStepComputation | 0 | 2 |
| Infeasible | 0 | 11 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-10) | 0 | 123 |
| IpoptStatus(3) | 0 | 1 |
| IpoptStatus(4) | 0 | 2 |
| LocalInfeasibility | 66 | 0 |
| MaxIterations | 6 | 11 |
| NumericalError | 63 | 0 |
| RestorationFailed | 9 | 4 |
| Timeout | 14 | 11 |

### Electrolyte Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |

### Grid Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| MaxIterations | 1 | 0 |

### CHO Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| MaxIter | 0 | 1 |
| MaxIterations | 1 | 0 |

## Regressions (Ipopt solves, ripopt fails)

| Problem | Suite | n | m | ripopt status | Ipopt obj |
|---------|-------|---|---|--------------|-----------|
| ACOPP30 | CUTEst | 72 | 142 | NumericalError | 5.768924e+02 |
| ACOPR30 | CUTEst | 72 | 172 | NumericalError | 5.768924e+02 |
| ALLINITC | CUTEst | 4 | 1 | NumericalError | 3.049261e+01 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | RestorationFailed | 7.862467e-01 |
| DECONVBNE | CUTEst | 63 | 40 | NumericalError | 0.000000e+00 |
| DISCS | CUTEst | 36 | 66 | NumericalError | 1.200007e+01 |
| DUALC8 | CUTEst | 8 | 503 | Timeout | 1.830936e+04 |
| FLETCHER | CUTEst | 4 | 4 | NumericalError | 1.165685e+01 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HIMMELBI | CUTEst | 100 | 12 | NumericalError | -1.735570e+03 |
| HS013 | HS | 2 | 1 | NumericalError | 9.945785e-01 |
| HS13 | CUTEst | 2 | 1 | NumericalError | 9.945785e-01 |
| HS225 | HS | 2 | 5 | NumericalError | 2.000000e+00 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.967522e-15 |
| KIRBY2LS | CUTEst | 5 | 0 | MaxIterations | 3.905074e+00 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| MAKELA3 | CUTEst | 21 | 20 | NumericalError | -9.810659e-09 |
| MGH10LS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| MUONSINELS | CUTEst | 1 | 0 | NumericalError | 4.387412e+04 |
| OET2 | CUTEst | 3 | 1002 | NumericalError | 8.715962e-02 |
| OET4 | CUTEst | 4 | 1002 | Timeout | 4.295421e-03 |
| OET5 | CUTEst | 5 | 1002 | MaxIterations | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | NumericalError | 4.465915e-05 |
| PALMER3 | CUTEst | 4 | 0 | NumericalError | 2.416980e+03 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | NumericalError | -9.136404e-03 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| SWOPF | CUTEst | 83 | 92 | NumericalError | 6.786018e-02 |
| TAXR13322 | CUTEst | 72 | 1261 | NumericalError | -3.429089e+02 |
| THURBERLS | CUTEst | 7 | 0 | NumericalError | 5.642708e+03 |
| VESUVIOLS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| case30_ieee | Grid | 72 | 142 | MaxIterations | 8.208515e+03 |

## Wins (ripopt solves, Ipopt fails) — 46 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| AVION2 | CUTEst | 49 | 15 | MaxIterations | 2.146364e+08 |
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BIGGS6NE | CUTEst | 6 | 13 | IpoptStatus(-10) | 0.000000e+00 |
| BLEACHNG | CUTEst | 17 | 0 | Timeout | 1.823872e+04 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651A | CUTEst | 7 | 61 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651D | CUTEst | 7 | 67 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 3.233018e-03 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DIAMON3DLS | CUTEst | 99 | 0 | Timeout | 5.754759e+02 |
| DMN15102LS | CUTEst | 66 | 0 | Timeout | 8.950802e+03 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | CUTEst | 9 | 3 | ErrorInStepComputation | -8.295477e+02 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996943e+03 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.213249e+00 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -1.850243e+09 |
| MGH17 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 7.795015e-02 |
| PALMER7A | CUTEst | 6 | 0 | MaxIterations | 1.089545e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -4.400000e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| RAT43 | CUTEst | 4 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | -7.564135e-09 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.294720e+00 |
| TAX13322 | CUTEst | 72 | 1261 | Timeout | -3.027143e+03 |
| TRO4X4 | CUTEst | 63 | 25 | IpoptStatus(4) | 8.999997e+00 |
| TRO6X2 | CUTEst | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

---
*Generated by benchmark_report.py*