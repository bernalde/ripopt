# ripopt Benchmark Report

Generated: 2026-04-08 18:28:09

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **703/865** (81.3%) | **688/865** (79.5%) |
| Acceptable | 0 | 5 |
| Total solved (Optimal + Acceptable) | 703 (81.3%) | 693 (80.1%) |
| Solved exclusively | 42 | 32 |
| Both solved | 661 | |
| Matching objectives (< 0.01%) | 555/661 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 117 (97.5%) | 116 (96.7%) | 2 | 1 | 115 | 109/115 |
| CUTEst | 727 | 571 (78.5%) | 561 (77.2%) | 40 | 30 | 531 | 433/531 |
| Electrolyte | 13 | 12 (92.3%) | 12 (92.3%) | 0 | 0 | 12 | 12/12 |
| OPF | 4 | 3 (75.0%) | 4 (100.0%) | 0 | 1 | 3 | 1/3 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 115 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 100us | 1.8ms |
| Total time | 38.1ms | 272.4ms |
| Mean iterations | 18.8 | 13.3 |
| Median iterations | 13 | 10 |

- **Geometric mean speedup**: 16.3x
- **Median speedup**: 14.7x
- ripopt faster: 113/115 (98%)
- ripopt 10x+ faster: 92/115
- Ipopt faster: 2/115

## CUTEst Suite — Performance

On 531 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 120us | 2.4ms |
| Total time | 69.11s | 11.61s |
| Mean iterations | 63.3 | 39.1 |
| Median iterations | 15 | 12 |

- **Geometric mean speedup**: 9.5x
- **Median speedup**: 19.8x
- ripopt faster: 450/531 (85%)
- ripopt 10x+ faster: 334/531
- Ipopt faster: 81/531

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 107us | 1.7ms |
| Total time | 2.2ms | 56.8ms |
| Mean iterations | 23.3 | 19.2 |
| Median iterations | 9 | 7 |

- **Geometric mean speedup**: 22.6x
- **Median speedup**: 22.1x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 10/12
- Ipopt faster: 0/12

## OPF Suite — Performance

On 3 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 9.8ms | 4.9ms |
| Total time | 51.2ms | 18.0ms |
| Mean iterations | 146.0 | 12.0 |
| Median iterations | 70 | 11 |

- **Geometric mean speedup**: 0.5x
- **Median speedup**: 0.9x
- ripopt faster: 1/3 (33%)
- ripopt 10x+ faster: 0/3
- Ipopt faster: 2/3

## Failure Analysis

### HS Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-11) | 0 | 1 |
| IpoptStatus(-199) | 0 | 1 |
| NumericalError | 2 | 0 |
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
| LocalInfeasibility | 68 | 0 |
| MaxIterations | 5 | 12 |
| NumericalError | 59 | 0 |
| RestorationFailed | 11 | 4 |
| Timeout | 13 | 10 |

### Electrolyte Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |
| NumericalError | 1 | 0 |

### OPF Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| NumericalError | 1 | 0 |

### CHO Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| MaxIter | 0 | 1 |
| MaxIterations | 1 | 0 |

## Regressions (Ipopt solves, ripopt fails)

| Problem | Suite | n | m | ripopt status | Ipopt obj |
|---------|-------|---|---|--------------|-----------|
| ACOPP30 | CUTEst | 72 | 142 | NumericalError | 5.768924e+02 |
| ACOPR30 | CUTEst | 72 | 172 | Timeout | 5.768924e+02 |
| BATCH | CUTEst | 48 | 73 | NumericalError | 2.591803e+05 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | RestorationFailed | 7.862467e-01 |
| DISCS | CUTEst | 36 | 66 | NumericalError | 1.200007e+01 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HS225 | HS | 2 | 5 | NumericalError | 2.000000e+00 |
| HS99EXP | CUTEst | 31 | 21 | NumericalError | -1.260006e+12 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.967522e-15 |
| KIRBY2LS | CUTEst | 5 | 0 | MaxIterations | 3.905074e+00 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| MCONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| MGH10LS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| MUONSINELS | CUTEst | 1 | 0 | NumericalError | 4.387412e+04 |
| OET2 | CUTEst | 3 | 1002 | NumericalError | 8.715962e-02 |
| OET5 | CUTEst | 5 | 1002 | Timeout | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| PALMER3 | CUTEst | 4 | 0 | NumericalError | 2.416980e+03 |
| QPCBLEND | CUTEst | 83 | 74 | NumericalError | -7.842801e-03 |
| QPNBLEND | CUTEst | 83 | 74 | NumericalError | -9.136404e-03 |
| SPANHYD | CUTEst | 97 | 33 | NumericalError | 2.397380e+02 |
| SSINE | CUTEst | 3 | 2 | NumericalError | 0.000000e+00 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| SWOPF | CUTEst | 83 | 92 | NumericalError | 6.786018e-02 |
| THURBERLS | CUTEst | 7 | 0 | NumericalError | 5.642708e+03 |
| VESUVIOLS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| case30_ieee | OPF | 72 | 142 | NumericalError | 8.208515e+03 |

## Wins (ripopt solves, Ipopt fails) — 42 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BENNETT5 | CUTEst | 3 | 154 | IpoptStatus(-10) | 0.000000e+00 |
| BIGGS6NE | CUTEst | 6 | 13 | IpoptStatus(-10) | 0.000000e+00 |
| BLEACHNG | CUTEst | 17 | 0 | Timeout | 1.823872e+04 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651A | CUTEst | 7 | 61 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651E | CUTEst | 7 | 64 | IpoptStatus(-10) | 0.000000e+00 |
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
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996927e+03 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.126755e+00 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -4.511177e+10 |
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
| SPIRAL | CUTEst | 3 | 2 | Infeasible | -1.109001e-10 |
| TAX13322 | CUTEst | 72 | 1261 | MaxIterations | -3.662555e+03 |
| TRO4X4 | CUTEst | 63 | 25 | IpoptStatus(4) | 8.999997e+00 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.002s | Optimal | 749 | 0.192s | 81.0x |
| SparseQP 1K | 500 | 500 | Optimal | 10 | 0.004s | Optimal | 6 | 0.004s | 0.8x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.0x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.003s | 0.5x |

ripopt: **4/4 solved** in 0.0s total
Ipopt: **4/4 solved** in 0.2s total

---
*Generated by benchmark_report.py*