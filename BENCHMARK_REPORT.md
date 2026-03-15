# ripopt Benchmark Report

Generated: 2026-03-15 10:43:06

## Executive Summary

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Total solved | **715/847** (84.4%) | 677/847 (79.9%) |
| Optimal | 475 | 672 |
| Acceptable | 240 | 5 |
| Solved exclusively | 48 | 10 |
| Both solved | 667 | |
| Matching objectives | 529/667 | |

## Per-Suite Summary

| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |
|-------|----------|--------------|-------------|-------------|------------|------------|-------|
| HS | 120 | 119 (99.2%) | 116 (96.7%) | 3 | 0 | 116 | 108/116 |
| CUTEst | 727 | 596 (82.0%) | 561 (77.2%) | 45 | 10 | 551 | 421/551 |

## HS Suite — Performance

On 116 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 100us | 1.9ms |
| Total time | 46.0ms | 281.9ms |
| Mean iterations | 16.4 | 13.3 |
| Median iterations | 11 | 10 |

- **Geometric mean speedup**: 16.5x
- **Median speedup**: 17.8x
- ripopt faster: 113/116 (97%)
- ripopt 10x+ faster: 92/116
- Ipopt faster: 3/116

## CUTEst Suite — Performance

On 551 commonly-solved problems:

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Median time | 97us | 2.3ms |
| Total time | 106.32s | 19.06s |
| Mean iterations | 130.8 | 41.1 |
| Median iterations | 18 | 12 |

- **Geometric mean speedup**: 9.1x
- **Median speedup**: 26.2x
- ripopt faster: 446/551 (81%)
- ripopt 10x+ faster: 350/551
- Ipopt faster: 105/551

## Failure Analysis

### HS Suite

| Failure Mode | ripopt | Ipopt |
|-------------|--------|-------|
| Infeasible | 0 | 1 |
| InvalidNumberDetected | 0 | 1 |
| IpoptStatus(-11) | 0 | 1 |
| IpoptStatus(-199) | 0 | 1 |
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
| LocalInfeasibility | 89 | 0 |
| MaxIterations | 4 | 12 |
| NumericalError | 17 | 0 |
| RestorationFailed | 3 | 4 |
| Timeout | 18 | 10 |

## Regressions (Ipopt solves, ripopt fails)

| Problem | Suite | n | m | ripopt status | Ipopt obj |
|---------|-------|---|---|--------------|-----------|
| ACOPR30 | CUTEst | 72 | 172 | Timeout | 5.768924e+02 |
| CORE1 | CUTEst | 65 | 59 | Timeout | 9.105624e+01 |
| CRESC50 | CUTEst | 6 | 100 | Timeout | 7.862467e-01 |
| DISCS | CUTEst | 36 | 66 | Timeout | 1.200007e+01 |
| FEEDLOC | CUTEst | 90 | 259 | Timeout | -9.539854e-10 |
| HS116 | CUTEst | 13 | 14 | NumericalError | 9.758747e+01 |
| MEYER3 | CUTEst | 3 | 0 | MaxIterations | 8.794586e+01 |
| MGH10LS | CUTEst | 3 | 0 | NumericalError | 8.794586e+01 |
| STRATEC | CUTEst | 10 | 0 | NumericalError | 2.212262e+03 |
| TRO3X3 | CUTEst | 30 | 13 | NumericalError | 8.967478e+00 |

## Wins (ripopt solves, Ipopt fails) — 48 problems

| Problem | Suite | n | m | Ipopt status | ripopt obj |
|---------|-------|---|---|-------------|------------|
| ARGAUSS | CUTEst | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| AVION2 | CUTEst | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | CUTEst | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 3.894509e-02 |
| DENSCHNBNE | CUTEst | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DENSCHNENE | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| DEVGLA1NE | CUTEst | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | CUTEst | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | CUTEst | 9 | 3 | ErrorInStepComputation | -8.295477e+02 |
| EXP2NE | CUTEst | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | CUTEst | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3LS | CUTEst | 6 | 0 | MaxIterations | 3.462949e-01 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS214 | HS | 2 | 0 | InvalidNumberDetected | 0.000000e+00 |
| HS223 | HS | 2 | 2 | Infeasible | -0.000000e+00 |
| HS25NE | CUTEst | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HS374 | HS | 10 | 35 | IpoptStatus(-199) | -5.497558e+14 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996858e+03 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS2 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS3 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | CUTEst | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 |
| LHAIFAM | CUTEst | 99 | 150 | InvalidNumberDetected | 6.931472e-01 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -3.245744e+00 |
| NYSTROM5 | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | CUTEst | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| OSBORNE1 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5A | CUTEst | 8 | 0 | MaxIterations | 3.047214e-01 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 7.795015e-02 |
| PALMER7A | CUTEst | 6 | 0 | MaxIterations | 1.089545e+01 |
| PALMER7E | CUTEst | 8 | 0 | MaxIterations | 1.016079e+01 |
| PFIT1 | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | CUTEst | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | CUTEst | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | CUTEst | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | CUTEst | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | CUTEst | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | CUTEst | 3 | 2 | Infeasible | 1.259927e-01 |
| SSI | CUTEst | 3 | 0 | MaxIterations | 3.363463e-04 |
| TAX13322 | CUTEst | 72 | 1261 | MaxIterations | -3.726455e+03 |
| TRO4X4 | CUTEst | 63 | 25 | IpoptStatus(4) | 9.284887e+00 |
| WACHBIEG | CUTEst | 3 | 2 | Infeasible | 1.000000e+00 |
| YFITNE | CUTEst | 3 | 17 | IpoptStatus(-10) | 0.000000e+00 |

## Acceptable (not Optimal) — 240 problems

These problems converged within relaxed tolerances but not strict tolerances.

| Problem | Suite | n | m | Ipopt status | ripopt obj | Ipopt obj |
|---------|-------|---|---|-------------|------------|-----------|
| 3PK | CUTEst | 30 | 0 | Optimal | 5.682734e+00 | 1.720119e+00 |
| ACOPR14 | CUTEst | 38 | 82 | Optimal | 8.081526e+03 | 8.081526e+03 |
| AIRPORT | CUTEst | 84 | 42 | Optimal | 4.791640e+04 | 4.795270e+04 |
| ALLINITA | CUTEst | 4 | 4 | Optimal | 3.329576e+01 | 3.329611e+01 |
| ALLINITC | CUTEst | 4 | 1 | Optimal | 3.049348e+01 | 3.049261e+01 |
| ANTWERP | CUTEst | 27 | 10 | Optimal | 3.245242e+03 | 3.245241e+03 |
| ARGAUSS | CUTEst | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| AVION2 | CUTEst | 49 | 15 | MaxIterations | 9.468013e+07 | 9.468013e+07 |
| BATCH | CUTEst | 48 | 73 | Optimal | 2.591866e+05 | 2.591803e+05 |
| BENNETT5LS | CUTEst | 3 | 0 | Optimal | 5.389456e-04 | 5.563289e-04 |
| BIGGSC4 | CUTEst | 4 | 7 | Optimal | -2.449980e+01 | -2.450000e+01 |
| BOX2 | CUTEst | 3 | 0 | Optimal | 6.334773e-14 | 6.251257e-23 |
| BOXBODLS | CUTEst | 2 | 0 | Optimal | 1.168009e+03 | 9.771500e+03 |
| BQPGABIM | CUTEst | 50 | 0 | Optimal | -3.774370e-05 | -3.790637e-05 |
| BQPGASIM | CUTEst | 50 | 0 | Optimal | -4.767945e-05 | -5.519663e-05 |
| BRKMCC | CUTEst | 2 | 0 | Optimal | 1.690427e-01 | 1.690427e-01 |
| BROWNBS | CUTEst | 2 | 0 | Optimal | 0.000000e+00 | 0.000000e+00 |
| BT13 | CUTEst | 5 | 1 | Optimal | 9.320027e-12 | -9.990000e-09 |
| BT8 | CUTEst | 5 | 2 | Optimal | 1.000000e+00 | 1.000000e+00 |
| CERI651ALS | CUTEst | 7 | 0 | Optimal | 3.348152e+02 | 3.348152e+02 |
| CERI651BLS | CUTEst | 7 | 0 | Optimal | 9.589532e+01 | 9.589532e+01 |
| CERI651CLS | CUTEst | 7 | 0 | Optimal | 5.657972e+01 | 5.657971e+01 |
| CERI651DLS | CUTEst | 7 | 0 | Optimal | 1.401753e+01 | 1.401753e+01 |
| CERI651ELS | CUTEst | 7 | 0 | Optimal | 2.368373e+01 | 2.368372e+01 |
| CHWIRUT2LS | CUTEst | 3 | 0 | Optimal | 5.130480e+02 | 5.130480e+02 |
| CLUSTERLS | CUTEst | 2 | 0 | Optimal | 1.958029e-14 | 2.724155e-18 |
| CONCON | CUTEst | 15 | 11 | Optimal | -6.230796e+03 | -6.230796e+03 |
| COOLHANSLS | CUTEst | 9 | 0 | Optimal | 2.587513e-05 | 1.208512e-18 |
| DECONVB | CUTEst | 63 | 0 | MaxIterations | 3.894509e-02 | 2.569475e-03 |
| DECONVBNE | CUTEst | 63 | 40 | Optimal | 0.000000e+00 | 0.000000e+00 |
| DECONVU | CUTEst | 63 | 0 | Optimal | 2.369401e-07 | 4.146188e-13 |
| DEGENLPA | CUTEst | 20 | 15 | Optimal | 3.060930e+00 | 3.054881e+00 |
| DEGENLPB | CUTEst | 20 | 15 | Optimal | -3.072266e+01 | -3.076401e+01 |
| DEMBO7 | CUTEst | 16 | 20 | Optimal | 1.747891e+02 | 1.747870e+02 |
| DENSCHND | CUTEst | 3 | 0 | Optimal | 3.820110e-09 | 2.221899e-04 |
| DENSCHNDNE | CUTEst | 3 | 3 | Acceptable | 0.000000e+00 | 0.000000e+00 |
| DENSCHNENE | CUTEst | 3 | 3 | Infeasible | 0.000000e+00 | 0.000000e+00 |
| DEVGLA2 | CUTEst | 5 | 0 | Optimal | 2.501146e-02 | 6.672171e-19 |
| DEVGLA2B | CUTEst | 5 | 0 | Optimal | 1.064046e+04 | 1.064046e+04 |
| DJTL | CUTEst | 2 | 0 | Acceptable | -8.951545e+03 | -8.951545e+03 |
| DNIEPER | CUTEst | 61 | 24 | Optimal | 1.870876e+04 | 1.874401e+04 |
| DUAL1 | CUTEst | 85 | 1 | Optimal | 3.501517e-02 | 3.501296e-02 |
| DUAL2 | CUTEst | 96 | 1 | Optimal | 3.373375e-02 | 3.373367e-02 |
| DUAL4 | CUTEst | 75 | 1 | Optimal | 7.460949e-01 | 7.460906e-01 |
| DUALC1 | CUTEst | 9 | 215 | Optimal | 6.155278e+03 | 6.155210e+03 |
| DUALC2 | CUTEst | 7 | 229 | Optimal | 3.551309e+03 | 3.551303e+03 |
| DUALC5 | CUTEst | 8 | 278 | Optimal | 4.272326e+02 | 4.272325e+02 |
| DUALC8 | CUTEst | 8 | 503 | Optimal | 2.125352e+04 | 1.830936e+04 |
| ECKERLE4LS | CUTEst | 3 | 0 | Optimal | 6.996961e-01 | 1.463589e-03 |
| EG1 | CUTEst | 3 | 0 | Optimal | -1.132800e+00 | -1.429307e+00 |
| ELATVIDUB | CUTEst | 2 | 0 | Optimal | 1.712780e+00 | 5.475112e+01 |
| ENSOLS | CUTEst | 9 | 0 | Optimal | 7.885398e+02 | 7.885398e+02 |
| ERRINBAR | CUTEst | 18 | 9 | Optimal | 2.804526e+01 | 2.804526e+01 |
| FBRAIN3LS | CUTEst | 6 | 0 | MaxIterations | 3.462949e-01 | 2.420258e-01 |
| GAUSS1LS | CUTEst | 8 | 0 | Optimal | 1.315823e+03 | 1.315822e+03 |
| GAUSS2LS | CUTEst | 8 | 0 | Optimal | 1.247528e+03 | 1.247528e+03 |
| GAUSS3LS | CUTEst | 8 | 0 | Optimal | 1.244485e+03 | 1.244485e+03 |
| GOULDQP1 | CUTEst | 32 | 17 | Optimal | -3.485319e+03 | -3.485333e+03 |
| GROUPING | CUTEst | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 | 0.000000e+00 |
| HAHN1LS | CUTEst | 7 | 0 | Optimal | 3.355356e+01 | 3.338424e+01 |
| HAIFAM | CUTEst | 99 | 150 | Optimal | 8.034756e-01 | -4.500036e+01 |
| HALDMADS | CUTEst | 6 | 42 | Optimal | 1.807809e-01 | 2.218282e+00 |
| HATFLDB | CUTEst | 4 | 0 | Optimal | 5.574709e-03 | 5.572808e-03 |
| HATFLDC | CUTEst | 25 | 0 | Optimal | 3.346038e-14 | 4.609311e-22 |
| HATFLDFL | CUTEst | 3 | 0 | Optimal | 6.180997e-05 | 6.016804e-05 |
| HATFLDGLS | CUTEst | 25 | 0 | Optimal | 2.193720e-13 | 5.431689e-27 |
| HATFLDH | CUTEst | 4 | 7 | Optimal | -2.450304e+01 | -2.450000e+01 |
| HEART6LS | CUTEst | 6 | 0 | Optimal | 3.304817e-02 | 9.463828e-23 |
| HEART8LS | CUTEst | 8 | 0 | Optimal | 1.140151e+00 | 5.961138e-29 |
| HIMMELBF | CUTEst | 4 | 0 | Optimal | 3.185717e+02 | 3.185717e+02 |
| HIMMELBJ | CUTEst | 45 | 14 | ErrorInStepComputation | -1.910345e+03 | -1.910345e+03 |
| HIMMELP1 | CUTEst | 2 | 0 | Optimal | -6.205394e+01 | -6.205394e+01 |
| HS002 | HS | 2 | 0 | Optimal | 5.044083e-02 | 4.941229e+00 |
| HS013 | HS | 2 | 1 | Optimal | 9.936482e-01 | 9.945785e-01 |
| HS026 | HS | 3 | 1 | Optimal | 2.944756e-12 | 6.661338e-16 |
| HS033 | HS | 3 | 2 | Optimal | -4.585786e+00 | -4.585787e+00 |
| HS045 | HS | 5 | 0 | Optimal | 1.000000e+00 | 1.000000e+00 |
| HS049 | HS | 5 | 2 | Optimal | 2.716664e-10 | 1.060052e-11 |
| HS107 | CUTEst | 9 | 6 | Optimal | 5.055012e+03 | 5.055012e+03 |
| HS108 | HS | 9 | 13 | Optimal | -8.660254e-01 | -8.660254e-01 |
| HS108 | CUTEst | 9 | 13 | Optimal | -8.660286e-01 | -5.000000e-01 |
| HS117 | CUTEst | 15 | 5 | Optimal | 3.234868e+01 | 3.234868e+01 |
| HS118 | CUTEst | 15 | 17 | Optimal | 6.648214e+02 | 6.648204e+02 |
| HS119 | CUTEst | 16 | 8 | Optimal | 2.449128e+02 | 2.448997e+02 |
| HS13 | CUTEst | 2 | 1 | Optimal | 9.936482e-01 | 9.945785e-01 |
| HS2 | CUTEst | 2 | 0 | Optimal | 5.044083e-02 | 4.941229e+00 |
| HS206 | HS | 2 | 0 | Optimal | -1.953993e-14 | -4.662937e-15 |
| HS213 | HS | 2 | 0 | Optimal | 1.338776e-09 | 6.182847e-06 |
| HS21MOD | CUTEst | 7 | 1 | Optimal | -9.595986e+01 | -9.596000e+01 |
| HS221 | HS | 2 | 1 | Optimal | -9.999997e-01 | -1.000000e+00 |
| HS225 | HS | 2 | 5 | Optimal | 2.000000e+00 | 2.000000e+00 |
| HS25 | CUTEst | 3 | 0 | Optimal | 6.143011e-05 | 3.082490e-20 |
| HS252 | HS | 3 | 1 | Optimal | 4.000000e-02 | 4.000000e-02 |
| HS254 | HS | 3 | 2 | Optimal | -1.732053e+00 | -1.732051e+00 |
| HS256 | HS | 4 | 0 | Optimal | 5.853569e-12 | 6.666382e-12 |
| HS257 | HS | 4 | 0 | Optimal | -5.884182e-14 | -7.061018e-14 |
| HS26 | CUTEst | 3 | 1 | Optimal | 2.944558e-12 | 1.291384e-16 |
| HS262 | HS | 4 | 4 | Optimal | -9.999994e+00 | -1.000000e+01 |
| HS263 | HS | 4 | 4 | Optimal | -1.000000e+00 | -1.000000e+00 |
| HS33 | CUTEst | 3 | 2 | Optimal | -4.585786e+00 | -4.585787e+00 |
| HS376 | HS | 10 | 15 | Optimal | -1.614947e+03 | -1.614959e+03 |
| HS44NEW | CUTEst | 4 | 6 | Optimal | -1.499987e+01 | -1.500000e+01 |
| HS45 | CUTEst | 5 | 0 | Optimal | 1.000000e+00 | 1.000000e+00 |
| HS46 | CUTEst | 5 | 2 | Optimal | 3.434445e-11 | 8.553352e-16 |
| HS47 | CUTEst | 5 | 3 | Optimal | 3.286145e-12 | 6.575160e-14 |
| HS49 | CUTEst | 5 | 2 | Optimal | 2.716654e-10 | 1.059996e-11 |
| HS76I | CUTEst | 4 | 3 | Optimal | -4.681818e+00 | -4.681818e+00 |
| HS83 | CUTEst | 5 | 3 | Optimal | -3.066561e+04 | -3.066554e+04 |
| HS87 | CUTEst | 6 | 4 | MaxIterations | 8.996858e+03 | 8.913970e+03 |
| HS97 | CUTEst | 6 | 4 | Optimal | 3.135809e+00 | 3.135806e+00 |
| HS98 | CUTEst | 6 | 4 | Optimal | 3.135810e+00 | 3.135806e+00 |
| HUMPS | CUTEst | 2 | 0 | Optimal | 2.780287e-14 | 2.761039e-17 |
| HYDC20LS | CUTEst | 99 | 0 | Optimal | 5.531032e-01 | 2.967522e-15 |
| HYDCAR6LS | CUTEst | 29 | 0 | Optimal | 3.053228e-01 | 2.863560e-18 |
| KIRBY2LS | CUTEst | 5 | 0 | Optimal | 3.905074e+00 | 3.905074e+00 |
| KOEBHELB | CUTEst | 3 | 0 | Optimal | 1.122151e+02 | 7.751635e+01 |
| LANCZOS1 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| LANCZOS1LS | CUTEst | 6 | 0 | Optimal | 4.345769e-06 | 4.715816e-21 |
| LANCZOS2 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| LANCZOS2LS | CUTEst | 6 | 0 | Optimal | 4.365561e-06 | 2.229944e-11 |
| LANCZOS3 | CUTEst | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| LAUNCH | CUTEst | 25 | 28 | Optimal | 9.004976e+00 | 9.004902e+00 |
| LEVYMONT5 | CUTEst | 2 | 0 | Optimal | 7.773811e+00 | 1.239502e-25 |
| LEWISPOL | CUTEst | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 | 0.000000e+00 |
| LINSPANH | CUTEst | 97 | 33 | Optimal | -7.700000e+01 | -7.700005e+01 |
| LOADBAL | CUTEst | 31 | 31 | Optimal | 4.528529e-01 | 4.528510e-01 |
| LOGHAIRY | CUTEst | 2 | 0 | Optimal | 5.350676e+00 | 1.823216e-01 |
| LOOTSMA | CUTEst | 3 | 2 | Optimal | 1.414214e+00 | 1.414213e+00 |
| LRCOVTYPE | CUTEst | 54 | 0 | Optimal | 8.693444e-01 | 5.723072e-01 |
| LRIJCNN1 | CUTEst | 22 | 0 | Optimal | 2.730781e-01 | 2.671576e-01 |
| LSC2LS | CUTEst | 3 | 0 | Optimal | 1.333398e+01 | 1.333439e+01 |
| LSNNODOC | CUTEst | 5 | 4 | Optimal | 1.231129e+02 | 1.231124e+02 |
| MATRIX2 | CUTEst | 6 | 2 | Optimal | 7.910335e-10 | 6.962253e-12 |
| MAXLIKA | CUTEst | 8 | 0 | Optimal | 1.136307e+03 | 1.136307e+03 |
| MCONCON | CUTEst | 15 | 11 | Optimal | -6.230796e+03 | -6.230796e+03 |
| MESH | CUTEst | 41 | 48 | IpoptStatus(4) | -3.245744e+00 | -2.201155e+37 |
| METHANB8LS | CUTEst | 31 | 0 | Optimal | 1.485006e-03 | 3.659179e-26 |
| METHANL8LS | CUTEst | 31 | 0 | Optimal | 1.892675e-01 | 5.699326e-17 |
| MGH09LS | CUTEst | 4 | 0 | Optimal | 1.019673e-03 | 3.075056e-04 |
| MGH10SLS | CUTEst | 3 | 0 | Optimal | 1.417868e+09 | 8.794586e+01 |
| MGH17LS | CUTEst | 5 | 0 | Optimal | 5.872152e-05 | 7.898394e-05 |
| MISRA1ALS | CUTEst | 2 | 0 | Optimal | 1.245514e-01 | 1.245514e-01 |
| MISRA1CLS | CUTEst | 2 | 0 | Optimal | 4.096684e-02 | 4.096684e-02 |
| MISRA1DLS | CUTEst | 2 | 0 | Optimal | 5.641930e-02 | 5.641930e-02 |
| MISTAKE | CUTEst | 9 | 13 | Optimal | -1.000005e+00 | -1.000000e+00 |
| MRIBASIS | CUTEst | 36 | 55 | Optimal | 1.821790e+01 | 1.821790e+01 |
| MSS1 | CUTEst | 90 | 73 | Optimal | -1.600000e+01 | -1.400000e+01 |
| MUONSINELS | CUTEst | 1 | 0 | Optimal | 5.114204e+04 | 4.387412e+04 |
| NET1 | CUTEst | 48 | 57 | Optimal | 9.411943e+05 | 9.411943e+05 |
| OET2 | CUTEst | 3 | 1002 | Optimal | 8.709576e-02 | 8.715962e-02 |
| OET4 | CUTEst | 4 | 1002 | Optimal | 5.011192e-03 | 4.295421e-03 |
| OET5 | CUTEst | 5 | 1002 | Optimal | 5.972873e-03 | 2.650077e-03 |
| OET6 | CUTEst | 5 | 1002 | Optimal | 8.704970e-02 | 2.069727e-03 |
| OET7 | CUTEst | 7 | 1002 | Optimal | 8.710336e-02 | 4.465915e-05 |
| OPTCNTRL | CUTEst | 32 | 20 | Optimal | 5.500000e+02 | 5.500000e+02 |
| OSBORNE1 | CUTEst | 5 | 33 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| OSBORNEA | CUTEst | 5 | 0 | Optimal | 5.464895e-05 | 5.464895e-05 |
| OSBORNEB | CUTEst | 11 | 0 | Optimal | 4.013774e-02 | 4.013774e-02 |
| PALMER1 | CUTEst | 4 | 0 | Optimal | 1.175460e+04 | 1.175460e+04 |
| PALMER1A | CUTEst | 6 | 0 | Optimal | 9.013922e-02 | 8.988363e-02 |
| PALMER1B | CUTEst | 4 | 0 | Optimal | 3.447355e+00 | 3.447355e+00 |
| PALMER1C | CUTEst | 8 | 0 | Optimal | 9.767045e-02 | 9.759806e-02 |
| PALMER1D | CUTEst | 7 | 0 | Optimal | 6.526826e-01 | 6.526826e-01 |
| PALMER1E | CUTEst | 8 | 0 | Optimal | 1.145320e-01 | 8.352685e-04 |
| PALMER2A | CUTEst | 6 | 0 | Optimal | 1.741402e-02 | 1.710972e-02 |
| PALMER2B | CUTEst | 4 | 0 | Optimal | 6.232670e-01 | 6.232670e-01 |
| PALMER2C | CUTEst | 8 | 0 | Optimal | 1.406377e-01 | 1.436889e-02 |
| PALMER2E | CUTEst | 8 | 0 | Optimal | 2.681950e-02 | 2.065001e-04 |
| PALMER3 | CUTEst | 4 | 0 | Optimal | 2.265958e+03 | 2.416980e+03 |
| PALMER3A | CUTEst | 6 | 0 | Optimal | 2.177909e-02 | 2.043142e-02 |
| PALMER3C | CUTEst | 8 | 0 | Optimal | 1.576191e-01 | 1.953764e-02 |
| PALMER3E | CUTEst | 8 | 0 | Optimal | 3.022987e-02 | 5.074083e-05 |
| PALMER4 | CUTEst | 4 | 0 | Optimal | 2.424016e+03 | 2.285383e+03 |
| PALMER4A | CUTEst | 6 | 0 | Optimal | 4.075083e-02 | 4.060614e-02 |
| PALMER4B | CUTEst | 4 | 0 | Optimal | 6.835139e+00 | 6.835139e+00 |
| PALMER4C | CUTEst | 8 | 0 | Optimal | 3.729344e-01 | 5.031069e-02 |
| PALMER4E | CUTEst | 8 | 0 | Optimal | 6.354217e-02 | 1.480042e-04 |
| PALMER5A | CUTEst | 8 | 0 | MaxIterations | 3.047214e-01 | 3.872650e-02 |
| PALMER5B | CUTEst | 9 | 0 | Optimal | 1.066612e-01 | 9.752496e-03 |
| PALMER5D | CUTEst | 4 | 0 | Optimal | 8.733940e+01 | 8.733940e+01 |
| PALMER5E | CUTEst | 8 | 0 | MaxIterations | 7.795015e-02 | 2.086295e-02 |
| PALMER6A | CUTEst | 6 | 0 | Optimal | 6.210127e-02 | 5.594884e-02 |
| PALMER6C | CUTEst | 8 | 0 | Optimal | 1.209246e-01 | 1.638742e-02 |
| PALMER6E | CUTEst | 8 | 0 | Optimal | 4.369251e-02 | 2.239551e-04 |
| PALMER7A | CUTEst | 6 | 0 | MaxIterations | 1.089545e+01 | 1.033491e+01 |
| PALMER7C | CUTEst | 8 | 0 | Optimal | 4.320687e+00 | 6.019856e-01 |
| PALMER7E | CUTEst | 8 | 0 | MaxIterations | 1.016079e+01 | 6.572691e+00 |
| PALMER8A | CUTEst | 6 | 0 | Optimal | 7.400970e-02 | 7.400970e-02 |
| PALMER8C | CUTEst | 8 | 0 | Optimal | 5.905795e-01 | 1.597681e-01 |
| PALMER8E | CUTEst | 8 | 0 | Optimal | 1.089516e-01 | 6.339308e-03 |
| PARKCH | CUTEst | 15 | 0 | Optimal | 1.623746e+03 | 1.623743e+03 |
| PFIT1LS | CUTEst | 3 | 0 | Optimal | 5.913197e-06 | 1.376372e-20 |
| PFIT2LS | CUTEst | 3 | 0 | Optimal | 1.844522e-03 | 1.896973e-20 |
| PFIT3LS | CUTEst | 3 | 0 | Optimal | 1.081693e-02 | 3.063567e-20 |
| PFIT4LS | CUTEst | 3 | 0 | Optimal | 6.275102e-04 | 6.664019e-20 |
| POLAK4 | CUTEst | 3 | 3 | Optimal | -5.434980e-09 | -9.965951e-09 |
| POLAK5 | CUTEst | 3 | 2 | Optimal | 5.000000e+01 | 5.000000e+01 |
| PORTFL1 | CUTEst | 12 | 1 | Optimal | 2.048692e-02 | 2.048627e-02 |
| PORTFL2 | CUTEst | 12 | 1 | Optimal | 2.968973e-02 | 2.968924e-02 |
| PORTFL3 | CUTEst | 12 | 1 | Optimal | 3.275144e-02 | 3.274971e-02 |
| PORTFL4 | CUTEst | 12 | 1 | Optimal | 2.631113e-02 | 2.630695e-02 |
| PORTFL6 | CUTEst | 12 | 1 | Optimal | 2.579262e-02 | 2.579180e-02 |
| PRICE4 | CUTEst | 2 | 0 | Optimal | 1.847439e-11 | 4.597166e-24 |
| PRICE4B | CUTEst | 2 | 0 | Optimal | 2.188710e-11 | 7.746134e-24 |
| PRICE4NE | CUTEst | 2 | 2 | Acceptable | 0.000000e+00 | 0.000000e+00 |
| PRODPL0 | CUTEst | 60 | 29 | Optimal | 5.879094e+01 | 5.879010e+01 |
| PRODPL1 | CUTEst | 60 | 29 | Optimal | 3.573928e+01 | 3.573897e+01 |
| PSPDOC | CUTEst | 4 | 0 | Optimal | 2.414235e+00 | 2.414214e+00 |
| RECIPE | CUTEst | 3 | 3 | Optimal | 0.000000e+00 | 0.000000e+00 |
| RECIPELS | CUTEst | 3 | 0 | Optimal | 4.147763e-12 | 1.585696e-11 |
| RES | CUTEst | 20 | 14 | Optimal | 0.000000e+00 | 0.000000e+00 |
| RK23 | CUTEst | 17 | 11 | Optimal | 8.333684e-02 | 8.333327e-02 |
| ROSZMAN1LS | CUTEst | 4 | 0 | Optimal | 3.951903e-02 | 4.948485e-04 |
| SANTALS | CUTEst | 21 | 0 | Optimal | 1.224627e-05 | 1.224358e-05 |
| SISSER | CUTEst | 2 | 0 | Optimal | 1.229019e-10 | 6.331104e-13 |
| SISSER2 | CUTEst | 2 | 0 | Optimal | 2.057380e-11 | 8.121606e-13 |
| SPANHYD | CUTEst | 97 | 33 | Optimal | 2.397380e+02 | 2.397380e+02 |
| SSI | CUTEst | 3 | 0 | MaxIterations | 3.363463e-04 | 1.381400e-09 |
| SSINE | CUTEst | 3 | 2 | Optimal | 0.000000e+00 | 0.000000e+00 |
| SYNTHES1 | CUTEst | 6 | 6 | Optimal | 7.592849e-01 | 7.592842e-01 |
| SYNTHES2 | CUTEst | 11 | 14 | Optimal | -5.544049e-01 | -5.544063e-01 |
| SYNTHES3 | CUTEst | 17 | 23 | Optimal | 1.508396e+01 | 1.508219e+01 |
| TAXR13322 | CUTEst | 72 | 1261 | Acceptable | -3.611652e+02 | -3.429089e+02 |
| TFI1 | CUTEst | 3 | 101 | Optimal | 2.876995e+02 | 5.334687e+00 |
| THURBERLS | CUTEst | 7 | 0 | Optimal | 5.642708e+03 | 5.642708e+03 |
| TOINTGOR | CUTEst | 50 | 0 | Optimal | 1.373905e+03 | 1.373905e+03 |
| TOINTPSP | CUTEst | 50 | 0 | Optimal | 2.255604e+02 | 2.255604e+02 |
| TOINTQOR | CUTEst | 50 | 0 | Optimal | 1.175472e+03 | 1.175472e+03 |
| TRIGGER | CUTEst | 7 | 6 | Optimal | 0.000000e+00 | 0.000000e+00 |
| TRO4X4 | CUTEst | 63 | 25 | IpoptStatus(4) | 9.284887e+00 | -1.957476e+20 |
| TRUSPYR1 | CUTEst | 11 | 4 | Optimal | 1.122874e+01 | 1.122874e+01 |
| VESUVIALS | CUTEst | 8 | 0 | Optimal | 1.500440e+03 | 9.914100e+02 |
| VESUVIOLS | CUTEst | 8 | 0 | Optimal | 9.914100e+02 | 9.914100e+02 |
| VESUVIOULS | CUTEst | 8 | 0 | Optimal | 4.771138e-01 | 4.771138e-01 |
| VIBRBEAM | CUTEst | 8 | 0 | Optimal | 9.119267e+00 | 3.322376e-01 |
| WATER | CUTEst | 31 | 10 | Optimal | 1.054938e+04 | 1.054938e+04 |
| WEEDS | CUTEst | 3 | 0 | Optimal | 2.587277e+00 | 2.587277e+00 |
| WOMFLET | CUTEst | 3 | 3 | Optimal | -1.723896e-19 | 6.050000e+00 |
| YFITNE | CUTEst | 3 | 17 | IpoptStatus(-10) | 0.000000e+00 | 0.000000e+00 |
| ZY2 | CUTEst | 3 | 2 | Optimal | 2.000001e+00 | 2.000000e+00 |

## Large-Scale Synthetic Problems — ripopt vs Ipopt

Synthetic problems with known structure, up to 100K variables.
Both solvers receive the exact same NlpProblem struct via the Rust trait interface.

| Problem | n | m | ripopt | iters | time | Ipopt | iters | time | speedup |
|---------|---|---|--------|-------|------|-------|-------|------|---------|
| Rosenbrock 500 | 500 | 0 | Acceptable | 70 | 0.002s | Optimal | 749 | 0.210s | 88.3x |
| SparseQP 1K | 500 | 500 | Optimal | 8 | 0.009s | Optimal | 6 | 0.004s | 0.4x |
| Bratu 1K | 1,000 | 998 | Optimal | 3 | 0.002s | Optimal | 2 | 0.002s | 1.1x |
| OptControl 2.5K | 2,499 | 1,250 | Optimal | 1 | 0.006s | Optimal | 1 | 0.002s | 0.4x |
| Rosenbrock 5K | 5,000 | 0 | Acceptable | 74 | 0.010s | Failed | 3000 | 3.730s | 383.7x |
| Poisson 2.5K | 5,000 | 2,500 | Optimal | 1 | 0.026s | Optimal | 1 | 0.010s | 0.4x |
| Bratu 10K | 10,000 | 9,998 | Optimal | 11 | 0.127s | Optimal | 2 | 0.013s | 0.1x |
| OptControl 20K | 19,999 | 10,000 | Optimal | 1 | 0.198s | Optimal | 1 | 0.020s | 0.1x |
| Poisson 50K | 49,928 | 24,964 | Optimal | 1 | 1.714s | Optimal | 1 | 0.130s | 0.1x |
| SparseQP 100K | 50,000 | 50,000 | Optimal | 8 | 4.870s | Optimal | 6 | 0.312s | 0.1x |

ripopt: **10/10 solved** in 7.0s total
Ipopt: **9/10 solved** in 4.4s total

---
*Generated by benchmark_report.py*