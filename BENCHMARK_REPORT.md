# ripopt Benchmark Report

Generated: 2026-03-16 19:41:55

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Optimal | **646/865** (74.7%) | **688/865** (79.5%) |
| Acceptable | 0 | 5 |
| Total solved (Optimal + Acceptable) | 646 (74.7%) | 693 (80.1%) |
| Solved exclusively | 32 | 79 |
| Both solved | 614 | |
| Matching objectives (< 0.01%) | 510/614 | |

> **Note:** ripopt uses fallback strategies (L-BFGS Hessian, AL, SQP, slack
> reformulation) that Ipopt does not have, which accounts for much of the
> Acceptable count difference. The "Different Local Minima" section below
> lists Acceptable solutions where ripopt converged to a worse local minimum.

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 113 (94.2%) | 116 (96.7%) | 2 | 5 | 111 | 104/111 |
| CUTEst | 727 | 516 (71.0%) | 561 (77.2%) | 29 | 74 | 487 | 394/487 |
| Electrolyte | 13 | 13 (100.0%) | 12 (92.3%) | 1 | 0 | 12 | 12/12 |
| OPF | 4 | 4 (100.0%) | 4 (100.0%) | 0 | 0 | 4 | 0/4 |
| CHO | 1 | 0 (0.0%) | 0 (0.0%) | 0 | 0 | 0 | 0/1 |

## HS Suite — Performance

On 111 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 103us | 1.8ms |
| Total time | 24.1ms | 261.6ms |
| Mean iterations | 13.9 | 12.9 |
| Median iterations | 11 | 10 |

- **Geometric mean speedup**: 16.8x
- **Median speedup**: 15.7x
- ripopt faster: 110/111 (99%)
- ripopt 10x+ faster: 83/111
- Ipopt faster: 1/111

## CUTEst Suite — Performance

On 487 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 130us | 2.8ms |
| Total time | 56.50s | 10.65s |
| Mean iterations | 45.1 | 34.8 |
| Median iterations | 14 | 11 |

- **Geometric mean speedup**: 11.2x
- **Median speedup**: 21.0x
- ripopt faster: 425/487 (87%)
- ripopt 10x+ faster: 321/487
- Ipopt faster: 62/487

## Electrolyte Suite — Performance

On 12 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 141us | 1.7ms |
| Total time | 1.6ms | 56.2ms |
| Mean iterations | 22.0 | 19.2 |
| Median iterations | 8 | 7 |

- **Geometric mean speedup**: 23.7x
- **Median speedup**: 21.7x
- ripopt faster: 12/12 (100%)
- ripopt 10x+ faster: 9/12
- Ipopt faster: 0/12

## OPF Suite — Performance

On 4 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 336.2ms | 7.8ms |
| Total time | 20.71s | 68.5ms |
| Mean iterations | 111.8 | 12.5 |
| Median iterations | 46 | 14 |

- **Geometric mean speedup**: 0.1x
- **Median speedup**: 0.3x
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
| NumericalError | 6 | 0 |
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
| MaxIterations | 4 | 12 |
| NumericalError | 112 | 0 |
| RestorationFailed | 10 | 4 |
| Timeout | 16 | 10 |

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
| AIRPORT | CUTEst | 84 | 42 | NumericalError | 4.795270e+04 |
| ALLINITA | CUTEst | 4 | 4 | NumericalError | 3.329611e+01 |
| ALLINITC | CUTEst | 4 | 1 | NumericalError | 3.049261e+01 |
| ANTWERP | CUTEst | 27 | 10 | NumericalError | 3.245241e+03 |
| BATCH | CUTEst | 48 | 73 | NumericalError | 2.591803e+05 |
| BT8 | CUTEst | 5 | 2 | NumericalError | 1.000000e+00 |
| CERI651ALS | CUTEst | 7 | 0 | NumericalError | 3.348152e+02 |
| CERI651DLS | CUTEst | 7 | 0 | NumericalError | 1.401753e+01 |
| CERI651ELS | CUTEst | 7 | 0 | NumericalError | 2.368372e+01 |
| CONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| CORE1 | CUTEst | 65 | 59 | NumericalError | 9.105624e+01 |
| DEGENLPB | CUTEst | 20 | 15 | NumericalError | -3.076401e+01 |
| DEMBO7 | CUTEst | 16 | 20 | NumericalError | 1.747870e+02 |
| DENSCHND | CUTEst | 3 | 0 | NumericalError | 2.221899e-04 |
| DISCS | CUTEst | 36 | 66 | Timeout | 1.200007e+01 |
| DJTL | CUTEst | 2 | 0 | NumericalError | -8.951545e+03 |
| DNIEPER | CUTEst | 61 | 24 | NumericalError | 1.874401e+04 |
| DUALC1 | CUTEst | 9 | 215 | NumericalError | 6.155210e+03 |
| DUALC8 | CUTEst | 8 | 503 | Timeout | 1.830936e+04 |
| ENSOLS | CUTEst | 9 | 0 | NumericalError | 7.885398e+02 |
| FEEDLOC | CUTEst | 90 | 259 | NumericalError | -9.539854e-10 |
| GAUSS2LS | CUTEst | 8 | 0 | NumericalError | 1.247528e+03 |
| HAHN1LS | CUTEst | 7 | 0 | NumericalError | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | NumericalError | -4.500036e+01 |
| HALDMADS | CUTEst | 6 | 42 | NumericalError | 2.218282e+00 |
| HEART6LS | CUTEst | 6 | 0 | NumericalError | 9.463828e-23 |
| HS013 | HS | 2 | 1 | NumericalError | 9.945785e-01 |
| HS026 | HS | 3 | 1 | NumericalError | 6.661338e-16 |
| HS045 | HS | 5 | 0 | NumericalError | 1.000000e+00 |
| HS108 | CUTEst | 9 | 13 | NumericalError | -5.000000e-01 |
| HS116 | CUTEst | 13 | 14 | NumericalError | 9.758747e+01 |
| HS13 | CUTEst | 2 | 1 | NumericalError | 9.945785e-01 |
| HS213 | HS | 2 | 0 | NumericalError | 6.182847e-06 |
| HS225 | HS | 2 | 5 | NumericalError | 2.000000e+00 |
| HS26 | CUTEst | 3 | 1 | NumericalError | 1.291384e-16 |
| HS44NEW | CUTEst | 4 | 6 | NumericalError | -1.500000e+01 |
| HS45 | CUTEst | 5 | 0 | NumericalError | 1.000000e+00 |
| HS46 | CUTEst | 5 | 2 | NumericalError | 8.553352e-16 |
| HS84 | CUTEst | 5 | 3 | NumericalError | -5.280335e+06 |
| HS97 | CUTEst | 6 | 4 | NumericalError | 3.135806e+00 |
| HS99EXP | CUTEst | 31 | 21 | NumericalError | -1.260006e+12 |
| HYDC20LS | CUTEst | 99 | 0 | NumericalError | 2.967522e-15 |
| HYDCAR6LS | CUTEst | 29 | 0 | NumericalError | 2.863560e-18 |
| KIRBY2LS | CUTEst | 5 | 0 | MaxIterations | 3.905074e+00 |
| LINSPANH | CUTEst | 97 | 33 | NumericalError | -7.700005e+01 |
| LRCOVTYPE | CUTEst | 54 | 0 | Timeout | 5.723072e-01 |
| LSC2LS | CUTEst | 3 | 0 | NumericalError | 1.333439e+01 |
| LSNNODOC | CUTEst | 5 | 4 | NumericalError | 1.231124e+02 |
| MATRIX2 | CUTEst | 6 | 2 | NumericalError | 6.962253e-12 |
| MCONCON | CUTEst | 15 | 11 | NumericalError | -6.230796e+03 |
| MEYER3 | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MGH10LS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MGH10SLS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| MISTAKE | CUTEst | 9 | 13 | NumericalError | -1.000000e+00 |
| MSS1 | CUTEst | 90 | 73 | NumericalError | -1.400000e+01 |
| MUONSINELS | CUTEst | 1 | 0 | NumericalError | 4.387412e+04 |
| NET1 | CUTEst | 48 | 57 | NumericalError | 9.411943e+05 |
| OET2 | CUTEst | 3 | 1002 | Timeout | 8.715962e-02 |
| OET4 | CUTEst | 4 | 1002 | NumericalError | 4.295421e-03 |
| OET5 | CUTEst | 5 | 1002 | NumericalError | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | NumericalError | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | NumericalError | 4.465915e-05 |
| OPTCNTRL | CUTEst | 32 | 20 | NumericalError | 5.500000e+02 |
| PFIT4 | CUTEst | 3 | 3 | NumericalError | 0.000000e+00 |
| POLAK4 | CUTEst | 3 | 3 | NumericalError | -9.965951e-09 |
| POLAK5 | CUTEst | 3 | 2 | NumericalError | 5.000000e+01 |
| QCNEW | CUTEst | 9 | 3 | NumericalError | -8.065219e+02 |
| SPANHYD | CUTEst | 97 | 33 | NumericalError | 2.397380e+02 |
| SSINE | CUTEst | 3 | 2 | NumericalError | 0.000000e+00 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| TENBARS1 | CUTEst | 18 | 9 | NumericalError | 2.295373e+03 |
| TFI1 | CUTEst | 3 | 101 | NumericalError | 5.334687e+00 |
| THURBERLS | CUTEst | 7 | 0 | NumericalError | 5.642708e+03 |
| TRO3X3 | CUTEst | 30 | 13 | NumericalError | 8.967478e+00 |
| VESUVIALS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| VESUVIOLS | CUTEst | 8 | 0 | NumericalError | 9.914100e+02 |
| VESUVIOULS | CUTEst | 8 | 0 | NumericalError | 4.771138e-01 |
| WATER | CUTEst | 31 | 10 | NumericalError | 1.054938e+04 |

## Wins (ripopt solves, Ipopt fails) — 32 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BIGGS6NE | CUTEst | 6 | 13 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651D | CUTEst | 7 | 67 | IpoptStatus(-10) | 0.000000e+00 |
| CERI651E | CUTEst | 7 | 64 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 3.233018e-03 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| KOEBHELBNE | CUTEst | 3 | 156 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -7.941372e+07 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | -1.109070e-10 |
| Seawater speciation | Electrolyte | 15 | 8 | Infeasible | -1.348257e+00 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Optimal | 81 | 0.002s | Optimal | 749 | 0.196s | 85.7x |
| SparseQP 1K | 500 | 500 | Optimal | 8 | 0.008s | Optimal | 6 | 0.004s | 0.4x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.1x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.002s | 0.4x |
| Rosenbrock 5K | 5,000 | 0 | NumericalError | 82 | 17.234s | Failed | 3000 | 3.624s | 0.2x |
| Poisson 2.5K | 5,000 | 2,500 | Optimal | 1 | 0.026s | Optimal | 1 | 0.009s | 0.4x |
| Bratu 10K | 10,000 | 9,998 | Optimal | 11 | 0.125s | Optimal | 2 | 0.012s | 0.1x |
| OptControl 20K | 19,999 | 10,000 | Optimal | 1 | 0.192s | Optimal | 1 | 0.019s | 0.1x |
| Poisson 50K | 49,928 | 24,964 | Optimal | 1 | 1.724s | Optimal | 1 | 0.121s | 0.1x |
| SparseQP 100K | 50,000 | 50,000 | Optimal | 8 | 4.736s | Optimal | 6 | 0.309s | 0.1x |

ripopt: **9/10 solved** in 24.1s total
Ipopt: **9/10 solved** in 4.3s total

---
*Generated by benchmark_report.py*