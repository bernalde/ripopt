# ripopt Benchmark Report

Generated: 2026-04-25 11:42:51

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **692/865** (80.0%) | **688/865** (79.5%) |
| Acceptable | 2 | 5 |
| Total solved (Optimal + Acceptable) | 694 (80.2%) | 693 (80.1%) |
| Solved exclusively | 42 | 41 |
| Both solved | 652 | |
| Matching objectives (< 0.01%) | 573/652 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 118 (98.3%) | 116 (96.7%) | 2 | 0 | 116 | 111/116 |
| CUTEst | 727 | 559 (76.9%) | 561 (77.2%) | 38 | 40 | 521 | 447/521 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| Grid | 4 | 3 (75.0%) | 4 (100.0%) | 0 | 1 | 3 | 3/3 |
| CHO | 1 | 1 (100.0%) | 0 (0.0%) | 1 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 116 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 83us | 1.9ms |
| Total time | 46.4ms | 285.7ms |
| Mean iterations | 67.0 | 13.3 |
| Median iterations | 13 | 10 |

- **Geometric mean speedup**: 20.8x
- **Median speedup**: 20.9x
- ripopt faster: 114/116 (98%)
- ripopt 10x+ faster: 103/116
- Ipopt faster: 2/116

## CUTEst Suite — Performance

On 521 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 98us | 2.3ms |
| Total time | 67.94s | 20.45s |
| Mean iterations | 62.3 | 39.5 |
| Median iterations | 16 | 12 |

- **Geometric mean speedup**: 12.0x
- **Median speedup**: 21.0x
- ripopt faster: 454/521 (87%)
- ripopt 10x+ faster: 340/521
- Ipopt faster: 67/521

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 92us | 1.1ms |
| Total time | 1.1ms | 46.0ms |
| Mean iterations | 22.0 | 19.2 |
| Median iterations | 8 | 7 |

- **Geometric mean speedup**: 25.4x
- **Median speedup**: 22.0x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 11/12
- Ipopt faster: 0/12

## Grid Suite — Performance

On 3 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 1.2ms | 3.3ms |
| Total time | 6.7ms | 13.6ms |
| Mean iterations | 15.7 | 12.0 |
| Median iterations | 14 | 11 |

- **Geometric mean speedup**: 2.8x
- **Median speedup**: 2.8x
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
| LocalInfeasibility | 76 | 0 |
| MaxIterations | 5 | 11 |
| NumericalError | 77 | 0 |
| RestorationFailed | 1 | 4 |
| Timeout | 9 | 11 |

### Electrolyte Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |

### Grid Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| NumericalError | 1 | 0 |

### CHO Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| MaxIter | 0 | 1 |

## Regressions (Ipopt solves, ripopt fails)

| Problem | Suite | n | m | ripopt status | Ipopt obj |
|---------|-------|---|---|--------------|-----------|
| ACOPP30 | CUTEst | 72 | 142 | NumericalError | 5.768924e+02 |
| ACOPR14 | CUTEst | 38 | 82 | NumericalError | 8.081526e+03 |
| ACOPR30 | CUTEst | 72 | 172 | NumericalError | 5.768924e+02 |
| AIRPORT | CUTEst | 84 | 42 | NumericalError | 4.795270e+04 |
| BATCH | CUTEst | 48 | 73 | NumericalError | 2.591803e+05 |
| BIGGSC4 | CUTEst | 4 | 7 | NumericalError | -2.450000e+01 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CERI651BLS | CUTEst | 7 | 0 | NumericalError | 9.589532e+01 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | NumericalError | 7.862467e-01 |
| DECONVBNE | CUTEst | 63 | 40 | NumericalError | 0.000000e+00 |
| DUALC8 | CUTEst | 8 | 503 | NumericalError | 1.830936e+04 |
| FEEDLOC | CUTEst | 90 | 259 | NumericalError | -9.539854e-10 |
| FLETCHER | CUTEst | 4 | 4 | NumericalError | 1.165685e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HATFLDFL | CUTEst | 3 | 0 | NumericalError | 6.016804e-05 |
| HATFLDH | CUTEst | 4 | 7 | NumericalError | -2.450000e+01 |
| HIMMELBI | CUTEst | 100 | 12 | NumericalError | -1.735570e+03 |
| HS118 | CUTEst | 15 | 17 | NumericalError | 6.648204e+02 |
| HS83 | CUTEst | 5 | 3 | NumericalError | -3.066554e+04 |
| HS84 | CUTEst | 5 | 3 | NumericalError | -5.280335e+06 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.967522e-15 |
| LAKES | CUTEst | 90 | 78 | NumericalError | 3.505251e+05 |
| LINSPANH | CUTEst | 97 | 33 | NumericalError | -7.700005e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| NET1 | CUTEst | 48 | 57 | NumericalError | 9.411943e+05 |
| OET2 | CUTEst | 3 | 1002 | NumericalError | 8.715962e-02 |
| OET4 | CUTEst | 4 | 1002 | NumericalError | 4.295421e-03 |
| OET5 | CUTEst | 5 | 1002 | MaxIterations | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | NumericalError | 4.465915e-05 |
| QCNEW | CUTEst | 9 | 3 | NumericalError | -8.065219e+02 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | NumericalError | -9.136404e-03 |
| SNAKE | CUTEst | 2 | 2 | NumericalError | -1.999999e-04 |
| SPANHYD | CUTEst | 97 | 33 | NumericalError | 2.397380e+02 |
| SSINE | CUTEst | 3 | 2 | MaxIterations | 0.000000e+00 |
| SWOPF | CUTEst | 83 | 92 | NumericalError | 6.786018e-02 |
| TFI1 | CUTEst | 3 | 101 | NumericalError | 5.334687e+00 |
| case30_ieee | Grid | 72 | 142 | NumericalError | 8.208515e+03 |

## Wins (ripopt solves, Ipopt fails) — 42 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BENNETT5 | CUTEst | 3 | 154 | IpoptStatus(-10) | 0.000000e+00 |
| BIGGS6NE | CUTEst | 6 | 13 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651C | CUTEst | 7 | 56 | IpoptStatus(-10) | 0.000000e+00 |
| CHO parmest | CHO | 21672 | 21660 | MaxIter | 6.767313e+04 |
| CRESC132 | CUTEst | 6 | 2654 | Timeout | 7.068752e-01 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 5.414582e-12 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA2NE | CUTEst | 5 | 16 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 3.000000e+00 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -1.071076e+07 |
| MGH17 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | CUTEst | 6 | 0 | MaxIterations | 1.089545e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 7.252303e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -4.371195e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| RAT43 | CUTEst | 4 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | -5.423447e-09 |
| SSI | CUTEst | 3 | 0 | MaxIterations | 3.363463e-04 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.347328e+00 |
| TRO4X4 | CUTEst | 63 | 25 | IpoptStatus(4) | 8.999997e+00 |
| TRO6X2 | CUTEst | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

## Acceptable (not Optimal) — 2 problems

These problems converged within relaxed tolerances but not strict tolerances.

| Problem | Suite | n | m | Ipopt status | ripopt obj | Ipopt obj |
|---------|-------|---|---|-------------|------------|-----------|
| ALLINITA | CUTEst | 4 | 4 | Optimal | 3.329676e+01 | 3.329611e+01 |
| ALLINITC | CUTEst | 4 | 1 | Optimal | 3.049299e+01 | 3.049261e+01 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.001s | Optimal | 749 | 0.201s | 169.1x |
| SparseQP 1K | 500 | 500 | Optimal | 13 | 0.183s | Optimal | 6 | 0.004s | 0.0x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.1x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.003s | 0.5x |

ripopt: **4/4 solved** in 0.2s total
Ipopt: **4/4 solved** in 0.2s total

---
*Generated by benchmark_report.py*