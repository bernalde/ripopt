# ripopt Benchmark Report

Generated: 2026-04-21 20:12:50

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **694/865** (80.2%) | **688/865** (79.5%) |
| Acceptable | 2 | 5 |
| Total solved (Optimal + Acceptable) | 696 (80.5%) | 693 (80.1%) |
| Solved exclusively | 40 | 37 |
| Both solved | 656 | |
| Matching objectives (< 0.01%) | 561/656 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 118 (98.3%) | 116 (96.7%) | 2 | 0 | 116 | 110/116 |
| CUTEst | 727 | 562 (77.3%) | 561 (77.2%) | 37 | 36 | 525 | 436/525 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| Grid | 4 | 3 (75.0%) | 4 (100.0%) | 0 | 1 | 3 | 3/3 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 116 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 119us | 1.9ms |
| Total time | 53.6ms | 281.9ms |
| Mean iterations | 70.8 | 13.3 |
| Median iterations | 13 | 10 |

- **Geometric mean speedup**: 15.0x
- **Median speedup**: 14.2x
- ripopt faster: 113/116 (97%)
- ripopt 10x+ faster: 84/116
- Ipopt faster: 3/116

## CUTEst Suite — Performance

On 525 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 97us | 2.2ms |
| Total time | 147.02s | 12.38s |
| Mean iterations | 84.1 | 43.0 |
| Median iterations | 15 | 12 |

- **Geometric mean speedup**: 9.9x
- **Median speedup**: 18.9x
- ripopt faster: 440/525 (84%)
- ripopt 10x+ faster: 332/525
- Ipopt faster: 85/525

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 135us | 1.3ms |
| Total time | 3.7ms | 49.6ms |
| Mean iterations | 29.3 | 19.2 |
| Median iterations | 8 | 7 |

- **Geometric mean speedup**: 17.5x
- **Median speedup**: 25.2x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 6/12
- Ipopt faster: 0/12

## Grid Suite — Performance

On 3 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 1.6ms | 5.6ms |
| Total time | 7.7ms | 17.5ms |
| Mean iterations | 15.7 | 12.0 |
| Median iterations | 14 | 11 |

- **Geometric mean speedup**: 2.8x
- **Median speedup**: 2.4x
- ripopt faster: 3/3 (100%)
- ripopt 10x+ faster: 0/3
- Ipopt faster: 0/3

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
| ErrorInStepComputation | 0 | 2 |
| Infeasible | 0 | 11 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-10) | 0 | 123 |
| IpoptStatus(3) | 0 | 1 |
| IpoptStatus(4) | 0 | 2 |
| LocalInfeasibility | 58 | 0 |
| MaxIterations | 12 | 12 |
| NumericalError | 72 | 0 |
| RestorationFailed | 9 | 4 |
| Timeout | 14 | 10 |

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
| ACOPR14 | CUTEst | 38 | 82 | NumericalError | 8.081526e+03 |
| ACOPR30 | CUTEst | 72 | 172 | NumericalError | 5.768924e+02 |
| AIRPORT | CUTEst | 84 | 42 | NumericalError | 4.795270e+04 |
| ANTWERP | CUTEst | 27 | 10 | NumericalError | 3.245241e+03 |
| BATCH | CUTEst | 48 | 73 | Timeout | 2.591803e+05 |
| BIGGSC4 | CUTEst | 4 | 7 | NumericalError | -2.450000e+01 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| FEEDLOC | CUTEst | 90 | 259 | NumericalError | -9.539854e-10 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HATFLDFL | CUTEst | 3 | 0 | NumericalError | 6.016804e-05 |
| HATFLDH | CUTEst | 4 | 7 | NumericalError | -2.450000e+01 |
| HIMMELBI | CUTEst | 100 | 12 | NumericalError | -1.735570e+03 |
| HS116 | CUTEst | 13 | 14 | NumericalError | 9.758747e+01 |
| HS118 | CUTEst | 15 | 17 | NumericalError | 6.648204e+02 |
| HS54 | CUTEst | 6 | 1 | NumericalError | -9.080748e-01 |
| HS83 | CUTEst | 5 | 3 | NumericalError | -3.066554e+04 |
| HS84 | CUTEst | 5 | 3 | NumericalError | -5.280335e+06 |
| HS85 | CUTEst | 5 | 21 | NumericalError | -2.215605e+00 |
| HS99EXP | CUTEst | 31 | 21 | NumericalError | -1.260006e+12 |
| LAKES | CUTEst | 90 | 78 | NumericalError | 3.505251e+05 |
| LINSPANH | CUTEst | 97 | 33 | NumericalError | -7.700005e+01 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| MISTAKE | CUTEst | 9 | 13 | NumericalError | -1.000000e+00 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| OET2 | CUTEst | 3 | 1002 | NumericalError | 8.715962e-02 |
| OET5 | CUTEst | 5 | 1002 | NumericalError | 2.650077e-03 |
| OET7 | CUTEst | 7 | 1002 | Timeout | 4.465915e-05 |
| QCNEW | CUTEst | 9 | 3 | NumericalError | -8.065219e+02 |
| QPNBLEND | CUTEst | 83 | 74 | MaxIterations | -9.136404e-03 |
| SPANHYD | CUTEst | 97 | 33 | MaxIterations | 2.397380e+02 |
| SWOPF | CUTEst | 83 | 92 | NumericalError | 6.786018e-02 |
| TAXR13322 | CUTEst | 72 | 1261 | NumericalError | -3.429089e+02 |
| VESUVIALS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| case30_ieee | Grid | 72 | 142 | MaxIterations | 8.208515e+03 |

## Wins (ripopt solves, Ipopt fails) — 40 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651B | CUTEst | 7 | 66 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651C | CUTEst | 7 | 56 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651E | CUTEst | 7 | 64 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 3.233123e-03 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA2NE | CUTEst | 5 | 16 | IpoptStatus(-10) | 0.000000e+00 |
| DMN37142LS | CUTEst | 66 | 0 | Timeout | 2.937451e+02 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.213785e+00 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -5.148561e+08 |
| MGH17 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| MGH17S | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7E | CUTEst | 8 | 0 | MaxIterations | 1.015395e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 6.104491e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -4.400000e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | 1.259274e-10 |
| SSI | CUTEst | 3 | 0 | MaxIterations | 3.363463e-04 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.340731e+00 |
| TAX13322 | CUTEst | 72 | 1261 | MaxIterations | -3.633286e+04 |
| TRO6X2 | CUTEst | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

## Acceptable (not Optimal) — 2 problems

These problems converged within relaxed tolerances but not strict tolerances.

| Problem | Suite | n | m | Ipopt status | ripopt obj | Ipopt obj |
|---------|-------|---|---|-------------|------------|-----------|
| ALLINITA | CUTEst | 4 | 4 | Optimal | 3.328898e+01 | 3.329611e+01 |
| ALLINITC | CUTEst | 4 | 1 | Optimal | 3.047967e+01 | 3.049261e+01 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.003s | Optimal | 749 | 0.199s | 76.2x |
| SparseQP 1K | 500 | 500 | Optimal | 13 | 0.176s | Optimal | 6 | 0.004s | 0.0x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.1x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.002s | 0.4x |

ripopt: **4/4 solved** in 0.2s total
Ipopt: **4/4 solved** in 0.2s total

---
*Generated by benchmark_report.py*