# ripopt Benchmark Report

Generated: 2026-03-28 20:35:12

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **697/865** (80.6%) | **688/865** (79.5%) |
| Acceptable | 0 | 5 |
| Total solved (Optimal + Acceptable) | 697 (80.6%) | 693 (80.1%) |
| Solved exclusively | 41 | 37 |
| Both solved | 656 | |
| Matching objectives (< 0.01%) | 541/656 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 118 (98.3%) | 116 (96.7%) | 2 | 0 | 116 | 108/116 |
| CUTEst | 727 | 562 (77.3%) | 561 (77.2%) | 38 | 37 | 524 | 419/524 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| OPF | 4 | 4 (100.0%) | 4 (100.0%) | 0 | 0 | 4 | 2/4 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 116 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 100us | 1.8ms |
| Total time | 66.6ms | 280.5ms |
| Mean iterations | 19.4 | 13.3 |
| Median iterations | 11 | 10 |

- **Geometric mean speedup**: 15.7x
- **Median speedup**: 16.8x
- ripopt faster: 113/116 (97%)
- ripopt 10x+ faster: 90/116
- Ipopt faster: 3/116

## CUTEst Suite — Performance

On 524 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 92us | 2.3ms |
| Total time | 103.21s | 9.27s |
| Mean iterations | 61.2 | 39.3 |
| Median iterations | 15 | 12 |

- **Geometric mean speedup**: 10.4x
- **Median speedup**: 24.3x
- ripopt faster: 442/524 (84%)
- ripopt 10x+ faster: 337/524
- Ipopt faster: 82/524

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 121us | 1.8ms |
| Total time | 3.3ms | 58.7ms |
| Mean iterations | 22.4 | 19.2 |
| Median iterations | 8 | 7 |

- **Geometric mean speedup**: 20.1x
- **Median speedup**: 15.7x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 9/12
- Ipopt faster: 0/12

## OPF Suite — Performance

On 4 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 149.3ms | 8.4ms |
| Total time | 30.30s | 67.3ms |
| Mean iterations | 37.2 | 12.5 |
| Median iterations | 31 | 14 |

- **Geometric mean speedup**: 0.0x
- **Median speedup**: 0.1x
- ripopt faster: 0/4 (0%)
- ripopt 10x+ faster: 0/4
- Ipopt faster: 4/4

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
| LocalInfeasibility | 69 | 0 |
| MaxIterations | 6 | 12 |
| NumericalError | 63 | 0 |
| RestorationFailed | 15 | 4 |
| Timeout | 12 | 10 |

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
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | RestorationFailed | 7.862467e-01 |
| DISCS | CUTEst | 36 | 66 | NumericalError | 1.200007e+01 |
| ELATTAR | CUTEst | 7 | 102 | NumericalError | 7.420618e+01 |
| FEEDLOC | CUTEst | 90 | 259 | NumericalError | -9.539854e-10 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HS84 | CUTEst | 5 | 3 | NumericalError | -5.280335e+06 |
| HS99EXP | CUTEst | 31 | 21 | NumericalError | -1.260006e+12 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.967522e-15 |
| HYDCAR20 | CUTEst | 99 | 99 | NumericalError | 0.000000e+00 |
| KIRBY2LS | CUTEst | 5 | 0 | MaxIterations | 3.905074e+00 |
| LAKES | CUTEst | 90 | 78 | RestorationFailed | 3.505251e+05 |
| LINSPANH | CUTEst | 97 | 33 | NumericalError | -7.700005e+01 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| MCONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| MGH10LS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| MUONSINELS | CUTEst | 1 | 0 | NumericalError | 4.387412e+04 |
| NET1 | CUTEst | 48 | 57 | NumericalError | 9.411943e+05 |
| OET5 | CUTEst | 5 | 1002 | MaxIterations | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | NumericalError | 4.465915e-05 |
| PFIT4 | CUTEst | 3 | 3 | NumericalError | 0.000000e+00 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | NumericalError | -9.136404e-03 |
| SSINE | CUTEst | 3 | 2 | NumericalError | 0.000000e+00 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| SWOPF | CUTEst | 83 | 92 | RestorationFailed | 6.786018e-02 |
| TAXR13322 | CUTEst | 72 | 1261 | NumericalError | -3.429089e+02 |
| THURBERLS | CUTEst | 7 | 0 | NumericalError | 5.642708e+03 |
| TRO3X3 | CUTEst | 30 | 13 | NumericalError | 8.967478e+00 |
| VESUVIOLS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |

## Wins (ripopt solves, Ipopt fails) — 41 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| AVION2 | CUTEst | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BIGGS6NE | CUTEst | 6 | 13 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651A | CUTEst | 7 | 61 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651E | CUTEst | 7 | 64 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 7.989309e-06 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DIAMON3DLS | CUTEst | 99 | 0 | Timeout | 5.754759e+02 |
| DMN15102LS | CUTEst | 66 | 0 | Timeout | 8.950802e+03 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | CUTEst | 9 | 3 | ErrorInStepComputation | -8.274326e+02 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996944e+03 |
| KOEBHELBNE | CUTEst | 3 | 156 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.126755e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 7.795015e-02 |
| PALMER7A | CUTEst | 6 | 0 | MaxIterations | 1.089545e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | -1.108995e-10 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.320459e+00 |
| TAX13322 | CUTEst | 72 | 1261 | MaxIterations | -1.018578e+04 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.003s | Optimal | 749 | 0.202s | 70.1x |
| SparseQP 1K | 500 | 500 | Optimal | 8 | 0.003s | Optimal | 6 | 0.004s | 1.4x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.1x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 5 | 0.007s | Optimal | 1 | 0.002s | 0.3x |
| Rosenbrock 5K | 5,000 | 0 | Optimal | 82 | 256.060s | Failed | 3000 | 3.699s | 0.0x |
| Poisson 2.5K | 5,000 | 2,500 | NumericalError | 81 | 855.404s | Optimal | 1 | 0.009s | N/A |
| Bratu 10K | 10,000 | 9,998 | Optimal | 10 | 0.132s | Optimal | 2 | 0.013s | 0.1x |
| OptControl 20K | 19,999 | 10,000 | Optimal | 6 | 0.249s | Optimal | 1 | 0.022s | 0.1x |
| Poisson 50K | 49,928 | 24,964 | Optimal | 1 | 1.787s | Optimal | 1 | 0.130s | 0.1x |
| SparseQP 100K | 50,000 | 50,000 | Optimal | 8 | 4.636s | Optimal | 6 | 0.310s | 0.1x |

ripopt: **9/10 solved** in 1118.3s total
Ipopt: **9/10 solved** in 4.4s total

---
*Generated by benchmark_report.py*