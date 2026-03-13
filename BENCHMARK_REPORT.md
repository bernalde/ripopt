# ripopt Benchmark Report

Generated: 2026-03-12 18:22:06

## Executive Summary

| Metric              | ripopt              | Ipopt           |
|---------------------|---------------------|-----------------|
| Total solved        | **703/847** (83.0%) | 679/847 (80.2%) |
| Optimal             | 510                 | 676             |
| Acceptable          | 193                 | 3               |
| Solved exclusively  | 43                  | 19              |
| Both solved         | 660                 |                 |
| Matching objectives | 549/660             |                 |

## Per-Suite Summary

| Suite  | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match   |
|--------|----------|---------------|--------------|-------------|------------|-------------|---------|
| HS     | 120      | 120 (100.0%)  | 118 (98.3%)  | 2           | 0          | 118         | 111/118 |
| CUTEst | 727      | 583 (80.2%)   | 561 (77.2%)  | 41          | 19         | 542         | 438/542 |

## HS Suite — Performance

On 118 commonly-solved problems:

| Metric            | ripopt  | Ipopt   |
|-------------------|---------|---------|
| Median time       | 29us    | 2.4ms   |
| Total time        | 183.0ms | 373.6ms |
| Mean iterations   | 16.3    | 14.0    |
| Median iterations | 11      | 10      |

- **Geometric mean speedup**: 54.2x
- **Median speedup**: 75.7x
- ripopt faster: 111/118 (94%)
- ripopt 10x+ faster: 107/118
- Ipopt faster: 7/118

## CUTEst Suite — Performance

On 542 commonly-solved problems:

| Metric            | ripopt  | Ipopt  |
|-------------------|---------|--------|
| Median time       | 146us   | 3.0ms  |
| Total time        | 125.02s | 31.59s |
| Mean iterations   | 169.6   | 37.2   |
| Median iterations | 17      | 12     |


- **Geometric mean speedup**: 7.5x
- **Median speedup**: 22.8x
- ripopt faster: 420/542 (77%)
- ripopt 10x+ faster: 340/542
- Ipopt faster: 122/542

## Failure Analysis

### HS Suite

| Failure Mode          | ripopt | Ipopt |
|-----------------------|--------|-------|
| Infeasible            | 0      | 1     |
| InvalidNumberDetected | 0      | 1     |

### CUTEst Suite

| Failure Mode           | ripopt | Ipopt |
|------------------------|--------|-------|
| ErrorInStepComputation | 0      | 2     |
| Infeasible             | 0      | 10    |
| InvalidNumberDetected  | 0      | 1     |
| IpoptStatus(-10)       | 0      | 123   |
| IpoptStatus(3)         | 0      | 1     |
| IpoptStatus(4)         | 0      | 2     |
| LocalInfeasibility     | 84     | 0     |
| MaxIterations          | 7      | 13    |
| NumericalError         | 22     | 0     |
| RestorationFailed      | 3      | 4     |
| Timeout                | 27     | 10    |
| Unbounded              | 1      | 0     |

## Regressions (Ipopt solves, ripopt fails)

| Problem    | Suite  | n  | m    | ripopt status  | Ipopt obj     |
|------------|--------|----|------|----------------|---------------|
| ACOPP30    | CUTEst | 72 | 142  | Timeout        | 5.768924e+02  |
| ACOPR30    | CUTEst | 72 | 172  | NumericalError | 5.768924e+02  |
| CERI651ALS | CUTEst | 7  | 0    | MaxIterations  | 3.348152e+02  |
| CRESC50    | CUTEst | 6  | 100  | Timeout        | 7.862467e-01  |
| DEMBO7     | CUTEst | 16 | 20   | NumericalError | 1.747870e+02  |
| DISCS      | CUTEst | 36 | 66   | Timeout        | 1.200007e+01  |
| DUALC8     | CUTEst | 8  | 503  | Timeout        | 1.830936e+04  |
| FEEDLOC    | CUTEst | 90 | 259  | Timeout        | -9.550425e-10 |
| HS109      | CUTEst | 9  | 10   | NumericalError | 5.362069e+03  |
| HS116      | CUTEst | 13 | 14   | NumericalError | 9.758747e+01  |
| HYDC20LS   | CUTEst | 99 | 0    | MaxIterations  | 2.486788e-17  |
| HYDCAR6LS  | CUTEst | 29 | 0    | MaxIterations  | 3.040548e-20  |
| LRCOVTYPE  | CUTEst | 54 | 0    | Timeout        | 5.723072e-01  |
| MEYER3     | CUTEst | 3  | 0    | MaxIterations  | 8.794586e+01  |
| MGH10LS    | CUTEst | 3  | 0    | NumericalError | 8.794586e+01  |
| OET5       | CUTEst | 5  | 1002 | Timeout        | 2.650077e-03  |
| PENTAGON   | CUTEst | 6  | 15   | Unbounded      | 1.365217e-04  |
| STRATEC    | CUTEst | 10 | 0    | NumericalError | 2.212262e+03  |
| VESUVIOLS  | CUTEst | 8  | 0    | NumericalError | 9.914100e+02  |

## Wins (ripopt solves, Ipopt fails) — 43 problems

| Problem    | Suite  | n   | m    | Ipopt status           | ripopt obj    |
|------------|--------|-----|------|------------------------|---------------|
| ARGAUSS    | CUTEst | 3   | 15   | IpoptStatus(-10)       | 0.000000e+00  |
| AVION2     | CUTEst | 49  | 15   | MaxIterations          | 9.468013e+07  |
| BEALENE    | CUTEst | 2   | 3    | IpoptStatus(-10)       | 0.000000e+00  |
| BIGGS6NE   | CUTEst | 6   | 13   | IpoptStatus(-10)       | 0.000000e+00  |
| BOX3NE     | CUTEst | 3   | 10   | IpoptStatus(-10)       | 0.000000e+00  |
| BROWNBSNE  | CUTEst | 2   | 3    | IpoptStatus(-10)       | 0.000000e+00  |
| DECONVB    | CUTEst | 63  | 0    | MaxIterations          | 1.291462e-08  |
| DENSCHNBNE | CUTEst | 2   | 3    | IpoptStatus(-10)       | 0.000000e+00  |
| DENSCHNENE | CUTEst | 3   | 3    | Infeasible             | 0.000000e+00  |
| DEVGLA1NE  | CUTEst | 4   | 24   | IpoptStatus(-10)       | 0.000000e+00  |
| ENGVAL2NE  | CUTEst | 3   | 5    | IpoptStatus(-10)       | 0.000000e+00  |
| EQC        | CUTEst | 9   | 3    | ErrorInStepComputation | -8.295477e+02 |
| EXP2NE     | CUTEst | 2   | 10   | IpoptStatus(-10)       | 0.000000e+00  |
| FBRAIN3    | CUTEst | 6   | 2211 | IpoptStatus(-10)       | 0.000000e+00  |
| GROUPING   | CUTEst | 100 | 125  | IpoptStatus(-10)       | 1.385040e+01  |
| GULFNE     | CUTEst | 3   | 99   | IpoptStatus(-10)       | 0.000000e+00  |
| HIMMELBJ   | CUTEst | 45  | 14   | ErrorInStepComputation | -1.910345e+03 |
| HS214      | HS     | 2   | 0    | InvalidNumberDetected  | 0.000000e+00  |
| HS223      | HS     | 2   | 2    | Infeasible             | -0.000000e+00 |
| HS87       | CUTEst | 6   | 4    | MaxIterations          | 8.996943e+03  |
| LANCZOS1   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  |
| LANCZOS2   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  |
| LANCZOS3   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  |
| LEVYMONE6  | CUTEst | 3   | 6    | IpoptStatus(-10)       | 0.000000e+00  |
| LEWISPOL   | CUTEst | 6   | 9    | IpoptStatus(-10)       | 1.212776e+00  |
| LHAIFAM    | CUTEst | 99  | 150  | InvalidNumberDetected  | 6.931472e-01  |
| LOGHAIRY   | CUTEst | 2   | 0    | MaxIterations          | 1.823216e-01  |
| MESH       | CUTEst | 41  | 48   | IpoptStatus(4)         | -1.129753e+06 |
| NYSTROM5   | CUTEst | 18  | 20   | IpoptStatus(-10)       | 0.000000e+00  |
| NYSTROM5C  | CUTEst | 18  | 20   | IpoptStatus(-10)       | 0.000000e+00  |
| OSBORNE1   | CUTEst | 5   | 33   | IpoptStatus(-10)       | 0.000000e+00  |
| PALMER2ENE | CUTEst | 8   | 23   | IpoptStatus(-10)       | 0.000000e+00  |
| PALMER5E   | CUTEst | 8   | 0    | MaxIterations          | 2.128087e+00  |
| PFIT1      | CUTEst | 3   | 3    | Infeasible             | 0.000000e+00  |
| PFIT2      | CUTEst | 3   | 3    | RestorationFailed      | 0.000000e+00  |
| POLAK3     | CUTEst | 12  | 10   | MaxIterations          | 7.084571e+00  |
| POWELLSQ   | CUTEst | 2   | 2    | Infeasible             | 0.000000e+00  |
| ROBOT      | CUTEst | 14  | 2    | IpoptStatus(3)         | 6.593299e+00  |
| SPIRAL     | CUTEst | 3   | 2    | MaxIterations          | 6.131462e-09  |
| TAX13322   | CUTEst | 72  | 1261 | MaxIterations          | -3.130697e+02 |
| TRO6X2     | CUTEst | 45  | 21   | RestorationFailed      | 1.225000e+03  |
| WACHBIEG   | CUTEst | 3   | 2    | Infeasible             | 1.000000e+00  |
| YFITNE     | CUTEst | 3   | 17   | IpoptStatus(-10)       | 0.000000e+00  |

## Acceptable (not Optimal) — 193 problems

These problems converged within relaxed tolerances but not strict tolerances.

| Problem    | Suite  | n   | m    | Ipopt status           | ripopt obj    | Ipopt obj     |
|------------|--------|-----|------|------------------------|---------------|---------------|
| ACOPR14    | CUTEst | 38  | 82   | Optimal                | 8.082583e+03  | 8.081526e+03  |
| ALLINITA   | CUTEst | 4   | 4    | Optimal                | 3.329605e+01  | 3.329611e+01  |
| ALLINITC   | CUTEst | 4   | 1    | Optimal                | 3.049283e+01  | 3.049261e+01  |
| ANTWERP    | CUTEst | 27  | 10   | Optimal                | 3.246598e+03  | 3.245241e+03  |
| ARGAUSS    | CUTEst | 3   | 15   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| AVION2     | CUTEst | 49  | 15   | MaxIterations          | 9.468013e+07  | 9.468013e+07  |
| BATCH      | CUTEst | 48  | 73   | Optimal                | 2.587698e+05  | 2.591803e+05  |
| BENNETT5LS | CUTEst | 3   | 0    | Optimal                | 5.389456e-04  | 5.563289e-04  |
| BOX2       | CUTEst | 3   | 0    | Optimal                | 6.334773e-14  | 6.251257e-23  |
| BOXBODLS   | CUTEst | 2   | 0    | Optimal                | 1.168009e+03  | 9.771500e+03  |
| BQPGABIM   | CUTEst | 50  | 0    | Optimal                | -3.356783e-05 | -3.790637e-05 |
| BQPGASIM   | CUTEst | 50  | 0    | Optimal                | -4.907072e-05 | -5.519663e-05 |
| BRKMCC     | CUTEst | 2   | 0    | Optimal                | 1.690427e-01  | 1.690427e-01  |
| BROWNBS    | CUTEst | 2   | 0    | Optimal                | 0.000000e+00  | 0.000000e+00  |
| BT13       | CUTEst | 5   | 1    | Optimal                | 9.320027e-12  | -9.990000e-09 |
| BT8        | CUTEst | 5   | 2    | Optimal                | 1.000000e+00  | 1.000000e+00  |
| CERI651DLS | CUTEst | 7   | 0    | Optimal                | 1.401753e+01  | 1.401753e+01  |
| CERI651ELS | CUTEst | 7   | 0    | Optimal                | 2.368372e+01  | 2.368372e+01  |
| CHWIRUT2LS | CUTEst | 3   | 0    | Optimal                | 5.130480e+02  | 5.130480e+02  |
| CLUSTERLS  | CUTEst | 2   | 0    | Optimal                | 1.958029e-14  | 2.724155e-18  |
| CONCON     | CUTEst | 15  | 11   | Optimal                | -6.230795e+03 | -6.230796e+03 |
| COOLHANSLS | CUTEst | 9   | 0    | Optimal                | 2.246968e-08  | 1.208368e-18  |
| DECONVB    | CUTEst | 63  | 0    | MaxIterations          | 1.291462e-08  | 2.569475e-03  |
| DECONVBNE  | CUTEst | 63  | 40   | Optimal                | 0.000000e+00  | 0.000000e+00  |
| DECONVU    | CUTEst | 63  | 0    | Optimal                | 2.369401e-07  | 8.777426e-12  |
| DEGENLPA   | CUTEst | 20  | 15   | Optimal                | 3.060434e+00  | 3.054881e+00  |
| DEGENLPB   | CUTEst | 20  | 15   | Optimal                | -3.069162e+01 | -3.076401e+01 |
| DENSCHND   | CUTEst | 3   | 0    | Optimal                | 3.820110e-09  | 2.221899e-04  |
| DENSCHNDNE | CUTEst | 3   | 3    | Acceptable             | 0.000000e+00  | 0.000000e+00  |
| DENSCHNENE | CUTEst | 3   | 3    | Infeasible             | 0.000000e+00  | 0.000000e+00  |
| DEVGLA2    | CUTEst | 5   | 0    | Optimal                | 2.501146e-02  | 6.672171e-19  |
| DEVGLA2B   | CUTEst | 5   | 0    | Optimal                | 1.064046e+04  | 1.064046e+04  |
| DGOSPEC    | CUTEst | 3   | 0    | Optimal                | -9.887540e+02 | -9.933506e+02 |
| DJTL       | CUTEst | 2   | 0    | Optimal                | -8.951545e+03 | -8.951545e+03 |
| DNIEPER    | CUTEst | 61  | 24   | Optimal                | 1.869902e+04  | 1.874401e+04  |
| DUAL1      | CUTEst | 85  | 1    | Optimal                | 3.501684e-02  | 3.501296e-02  |
| DUAL2      | CUTEst | 96  | 1    | Optimal                | 3.373368e-02  | 3.373367e-02  |
| DUAL4      | CUTEst | 75  | 1    | Optimal                | 7.460908e-01  | 7.460906e-01  |
| DUALC1     | CUTEst | 9   | 215  | Optimal                | 6.734849e+03  | 6.155210e+03  |
| DUALC2     | CUTEst | 7   | 229  | Optimal                | 4.036311e+03  | 3.551303e+03  |
| DUALC5     | CUTEst | 8   | 278  | Optimal                | 4.272326e+02  | 4.272325e+02  |
| ECKERLE4LS | CUTEst | 3   | 0    | Optimal                | 6.996961e-01  | 1.463589e-03  |
| ELATVIDUB  | CUTEst | 2   | 0    | Optimal                | 1.712780e+00  | 5.475112e+01  |
| ENSOLS     | CUTEst | 9   | 0    | Optimal                | 7.885398e+02  | 7.885398e+02  |
| ERRINBAR   | CUTEst | 18  | 9    | Optimal                | 2.804561e+01  | 2.804526e+01  |
| EXPFITC    | CUTEst | 5   | 502  | Optimal                | 1.283923e+01  | 2.330263e-02  |
| GAUSS2LS   | CUTEst | 8   | 0    | Optimal                | 1.247528e+03  | 1.247528e+03  |
| GOULDQP1   | CUTEst | 32  | 17   | Optimal                | -3.485332e+03 | -3.485333e+03 |
| GROUPING   | CUTEst | 100 | 125  | IpoptStatus(-10)       | 1.385040e+01  | 0.000000e+00  |
| HAHN1LS    | CUTEst | 7   | 0    | Optimal                | 3.937328e+01  | 3.338424e+01  |
| HAIFAM     | CUTEst | 99  | 150  | Optimal                | -4.500032e+01 | -4.500036e+01 |
| HALDMADS   | CUTEst | 6   | 42   | Optimal                | 9.942038e-02  | 2.218282e+00  |
| HATFLDC    | CUTEst | 25  | 0    | Optimal                | 3.346038e-14  | 4.609311e-22  |
| HATFLDFL   | CUTEst | 3   | 0    | Optimal                | 6.657448e-05  | 6.016849e-05  |
| HATFLDGLS  | CUTEst | 25  | 0    | Optimal                | 2.193720e-13  | 5.450251e-27  |
| HATFLDH    | CUTEst | 4   | 7    | Optimal                | -2.450304e+01 | -2.450000e+01 |
| HEART8LS   | CUTEst | 8   | 0    | Optimal                | 5.223015e-15  | 2.199258e-29  |
| HIMMELBF   | CUTEst | 4   | 0    | Optimal                | 3.185717e+02  | 3.185717e+02  |
| HIMMELBJ   | CUTEst | 45  | 14   | ErrorInStepComputation | -1.910345e+03 | -1.910345e+03 |
| HS013      | HS     | 2   | 1    | Optimal                | 9.933710e-01  | 9.945785e-01  |
| HS026      | HS     | 3   | 1    | Optimal                | 2.944534e-12  | 6.661338e-16  |
| HS030      | HS     | 3   | 1    | Optimal                | 1.000002e+00  | 1.000000e+00  |
| HS032      | HS     | 3   | 2    | Optimal                | 1.000001e+00  | 1.000000e+00  |
| HS033      | HS     | 3   | 2    | Optimal                | -4.585786e+00 | -4.585787e+00 |
| HS044      | HS     | 4   | 6    | Optimal                | -1.499996e+01 | -1.500000e+01 |
| HS045      | HS     | 5   | 0    | Optimal                | 1.000314e+00  | 1.000000e+00  |
| HS049      | HS     | 5   | 2    | Optimal                | 2.716664e-10  | 1.060052e-11  |
| HS107      | CUTEst | 9   | 6    | Optimal                | 5.055012e+03  | 5.055012e+03  |
| HS108      | HS     | 9   | 13   | Optimal                | -8.660254e-01 | -8.660254e-01 |
| HS108      | CUTEst | 9   | 13   | Optimal                | -8.660306e-01 | -5.000000e-01 |
| HS116      | HS     | 13  | 15   | Optimal                | 9.758677e+01  | 9.758747e+01  |
| HS117      | CUTEst | 15  | 5    | Optimal                | 3.234869e+01  | 3.234868e+01  |
| HS118      | CUTEst | 15  | 17   | Optimal                | 6.648214e+02  | 6.648204e+02  |
| HS119      | CUTEst | 16  | 8    | Optimal                | 2.449212e+02  | 2.448997e+02  |
| HS13       | CUTEst | 2   | 1    | Optimal                | 9.933710e-01  | 9.945785e-01  |
| HS206      | HS     | 2   | 0    | Optimal                | -1.953993e-14 | -4.662937e-15 |
| HS213      | HS     | 2   | 0    | Optimal                | 1.338776e-09  | 6.182847e-06  |
| HS21MOD    | CUTEst | 7   | 1    | Optimal                | -9.596000e+01 | -9.596000e+01 |
| HS25       | CUTEst | 3   | 0    | Optimal                | 6.082309e-05  | 3.082490e-20  |
| HS252      | HS     | 3   | 1    | Optimal                | 4.000000e-02  | 4.000000e-02  |
| HS254      | HS     | 3   | 2    | Optimal                | -1.732053e+00 | -1.732051e+00 |
| HS256      | HS     | 4   | 0    | Optimal                | 5.853569e-12  | 6.666382e-12  |
| HS257      | HS     | 4   | 0    | Optimal                | -5.884182e-14 | -7.061018e-14 |
| HS26       | CUTEst | 3   | 1    | Optimal                | 2.944558e-12  | 1.291384e-16  |
| HS263      | HS     | 4   | 4    | Optimal                | -1.000000e+00 | -1.000000e+00 |
| HS30       | CUTEst | 3   | 1    | Optimal                | 1.000002e+00  | 1.000000e+00  |
| HS32       | CUTEst | 3   | 2    | Optimal                | 1.000001e+00  | 1.000000e+00  |
| HS33       | CUTEst | 3   | 2    | Optimal                | -4.585786e+00 | -4.585787e+00 |
| HS374      | HS     | 10  | 35   | Optimal                | 2.332776e-01  | 2.332635e-01  |
| HS376      | HS     | 10  | 15   | Optimal                | -1.614957e+03 | -1.614959e+03 |
| HS44       | CUTEst | 4   | 6    | Optimal                | -1.499996e+01 | -1.500000e+01 |
| HS44NEW    | CUTEst | 4   | 6    | Optimal                | -1.499983e+01 | -1.500000e+01 |
| HS45       | CUTEst | 5   | 0    | Optimal                | 1.000314e+00  | 1.000000e+00  |
| HS46       | CUTEst | 5   | 2    | Optimal                | 3.434445e-11  | 8.553352e-16  |
| HS47       | CUTEst | 5   | 3    | Optimal                | 3.286145e-12  | 6.575160e-14  |
| HS49       | CUTEst | 5   | 2    | Optimal                | 2.716654e-10  | 1.059996e-11  |
| HS55       | CUTEst | 6   | 6    | Optimal                | 6.666667e+00  | 6.333333e+00  |
| HS83       | CUTEst | 5   | 3    | Optimal                | -3.066561e+04 | -3.066554e+04 |
| HS87       | CUTEst | 6   | 4    | MaxIterations          | 8.996943e+03  | 8.913970e+03  |
| HS96       | CUTEst | 6   | 4    | Optimal                | 1.573793e-02  | 1.561775e-02  |
| HS97       | CUTEst | 6   | 4    | Optimal                | 3.135809e+00  | 3.135806e+00  |
| KIRBY2LS   | CUTEst | 5   | 0    | Optimal                | 3.905074e+00  | 3.905074e+00  |
| LANCZOS1   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| LANCZOS1LS | CUTEst | 6   | 0    | Optimal                | 4.345769e-06  | 2.220745e-16  |
| LANCZOS2   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| LANCZOS2LS | CUTEst | 6   | 0    | Optimal                | 4.365561e-06  | 2.229943e-11  |
| LANCZOS3   | CUTEst | 6   | 24   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| LAUNCH     | CUTEst | 25  | 28   | Optimal                | 1.111913e+01  | 9.004902e+00  |
| LEVYMONT5  | CUTEst | 2   | 0    | Optimal                | 7.773811e+00  | 1.239502e-25  |
| LEWISPOL   | CUTEst | 6   | 9    | IpoptStatus(-10)       | 1.212776e+00  | 0.000000e+00  |
| LINSPANH   | CUTEst | 97  | 33   | Optimal                | -7.700002e+01 | -7.700005e+01 |
| LRIJCNN1   | CUTEst | 22  | 0    | Optimal                | 2.730781e-01  | 2.671576e-01  |
| LSC2LS     | CUTEst | 3   | 0    | Optimal                | 1.333399e+01  | 1.333438e+01  |
| LSNNODOC   | CUTEst | 5   | 4    | Optimal                | 1.231129e+02  | 1.231124e+02  |
| MATRIX2    | CUTEst | 6   | 2    | Optimal                | 7.898155e-10  | 6.962253e-12  |
| MAXLIKA    | CUTEst | 8   | 0    | Optimal                | 1.149346e+03  | 1.136307e+03  |
| MCONCON    | CUTEst | 15  | 11   | Optimal                | -6.230795e+03 | -6.230796e+03 |
| MGH09LS    | CUTEst | 4   | 0    | Optimal                | 1.019673e-03  | 3.075056e-04  |
| MGH10SLS   | CUTEst | 3   | 0    | Optimal                | 1.417866e+09  | 8.794586e+01  |
| MGH17LS    | CUTEst | 5   | 0    | Optimal                | 2.451828e-02  | 7.928588e-05  |
| MISRA1ALS  | CUTEst | 2   | 0    | Optimal                | 1.245514e-01  | 1.245514e-01  |
| MISRA1BLS  | CUTEst | 2   | 0    | Optimal                | 7.546468e-02  | 7.546468e-02  |
| MISRA1CLS  | CUTEst | 2   | 0    | Optimal                | 4.096684e-02  | 4.096684e-02  |
| MISRA1DLS  | CUTEst | 2   | 0    | Optimal                | 5.641930e-02  | 5.641930e-02  |
| MISTAKE    | CUTEst | 9   | 13   | Optimal                | -1.000001e+00 | -1.000000e+00 |
| NET1       | CUTEst | 48  | 57   | Optimal                | 9.411943e+05  | 9.411943e+05  |
| OET2       | CUTEst | 3   | 1002 | Optimal                | 8.709576e-02  | 8.715962e-02  |
| OET3       | CUTEst | 4   | 1002 | Optimal                | 4.585075e-03  | 4.505043e-03  |
| OET4       | CUTEst | 4   | 1002 | Optimal                | 5.011192e-03  | 4.295421e-03  |
| OET6       | CUTEst | 5   | 1002 | Optimal                | 8.704970e-02  | 2.069727e-03  |
| OET7       | CUTEst | 7   | 1002 | Optimal                | 8.710336e-02  | 4.465915e-05  |
| OSBORNE1   | CUTEst | 5   | 33   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| OSBORNEA   | CUTEst | 5   | 0    | Optimal                | 5.464895e-05  | 5.464895e-05  |
| OSBORNEB   | CUTEst | 11  | 0    | Optimal                | 4.013774e-02  | 4.013774e-02  |
| PALMER1A   | CUTEst | 6   | 0    | Optimal                | 8.988363e-02  | 8.988363e-02  |
| PALMER1E   | CUTEst | 8   | 0    | Optimal                | 8.352685e-04  | 8.352685e-04  |
| PALMER2    | CUTEst | 4   | 0    | Optimal                | 3.651098e+03  | 3.651098e+03  |
| PALMER2A   | CUTEst | 6   | 0    | Optimal                | 1.710972e-02  | 1.710972e-02  |
| PALMER2B   | CUTEst | 4   | 0    | Optimal                | 6.232670e-01  | 6.232670e-01  |
| PALMER2C   | CUTEst | 8   | 0    | Optimal                | 1.436889e-02  | 1.436889e-02  |
| PALMER2E   | CUTEst | 8   | 0    | Optimal                | 2.065001e-04  | 1.163043e-01  |
| PALMER2ENE | CUTEst | 8   | 23   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| PALMER3    | CUTEst | 4   | 0    | Optimal                | 2.265958e+03  | 2.416980e+03  |
| PALMER3A   | CUTEst | 6   | 0    | Optimal                | 2.043142e-02  | 2.043142e-02  |
| PALMER3E   | CUTEst | 8   | 0    | Optimal                | 5.074083e-05  | 5.074083e-05  |
| PALMER4A   | CUTEst | 6   | 0    | Optimal                | 4.060614e-02  | 4.060614e-02  |
| PALMER4B   | CUTEst | 4   | 0    | Optimal                | 6.835139e+00  | 6.835139e+00  |
| PALMER4E   | CUTEst | 8   | 0    | Optimal                | 1.480042e-04  | 1.480042e-04  |
| PALMER5B   | CUTEst | 9   | 0    | Optimal                | 1.523681e-02  | 9.752496e-03  |
| PALMER5E   | CUTEst | 8   | 0    | MaxIterations          | 2.128087e+00  | 2.086281e-02  |
| PALMER6C   | CUTEst | 8   | 0    | Optimal                | 1.638742e-02  | 1.638742e-02  |
| PALMER6E   | CUTEst | 8   | 0    | Optimal                | 2.239551e-04  | 2.239551e-04  |
| PALMER7C   | CUTEst | 8   | 0    | Optimal                | 6.019856e-01  | 6.019856e-01  |
| PALMER8A   | CUTEst | 6   | 0    | Optimal                | 7.400970e-02  | 7.400970e-02  |
| PALMER8C   | CUTEst | 8   | 0    | Optimal                | 1.597681e-01  | 1.597681e-01  |
| PALMER8E   | CUTEst | 8   | 0    | Optimal                | 6.339308e-03  | 6.339308e-03  |
| PARKCH     | CUTEst | 15  | 0    | Optimal                | 1.623743e+03  | 1.623743e+03  |
| POLAK4     | CUTEst | 3   | 3    | Optimal                | -5.434980e-09 | -9.965951e-09 |
| POLAK5     | CUTEst | 3   | 2    | Optimal                | 5.000000e+01  | 5.000000e+01  |
| PORTFL1    | CUTEst | 12  | 1    | Optimal                | 2.048756e-02  | 2.048627e-02  |
| PORTFL2    | CUTEst | 12  | 1    | Optimal                | 2.968968e-02  | 2.968924e-02  |
| PORTFL3    | CUTEst | 12  | 1    | Optimal                | 3.275121e-02  | 3.274971e-02  |
| PORTFL4    | CUTEst | 12  | 1    | Optimal                | 2.631078e-02  | 2.630695e-02  |
| PORTFL6    | CUTEst | 12  | 1    | Optimal                | 2.579230e-02  | 2.579180e-02  |
| PRICE4     | CUTEst | 2   | 0    | Optimal                | 1.847440e-11  | 4.597166e-24  |
| PRICE4B    | CUTEst | 2   | 0    | Optimal                | 2.188723e-11  | 7.746134e-24  |
| PRICE4NE   | CUTEst | 2   | 2    | Acceptable             | 0.000000e+00  | 0.000000e+00  |
| PRODPL0    | CUTEst | 60  | 29   | Optimal                | 5.879018e+01  | 5.879010e+01  |
| PRODPL1    | CUTEst | 60  | 29   | Optimal                | 3.574225e+01  | 3.573897e+01  |
| RECIPE     | CUTEst | 3   | 3    | Optimal                | 0.000000e+00  | 0.000000e+00  |
| RECIPELS   | CUTEst | 3   | 0    | Optimal                | 4.147763e-12  | 1.585696e-11  |
| RES        | CUTEst | 20  | 14   | Optimal                | 0.000000e+00  | 0.000000e+00  |
| RK23       | CUTEst | 17  | 11   | Optimal                | 8.339525e-02  | 8.333327e-02  |
| ROSZMAN1LS | CUTEst | 4   | 0    | Optimal                | 3.951903e-02  | 4.948485e-04  |
| SANTALS    | CUTEst | 21  | 0    | Optimal                | 1.224608e-05  | 1.224358e-05  |
| SIPOW3     | CUTEst | 4   | 2000 | Optimal                | 5.346586e-01  | 5.346586e-01  |
| SISSER     | CUTEst | 2   | 0    | Optimal                | 1.229019e-10  | 6.331104e-13  |
| SISSER2    | CUTEst | 2   | 0    | Optimal                | 2.057380e-11  | 8.121606e-13  |
| SPANHYD    | CUTEst | 97  | 33   | Optimal                | 2.397380e+02  | 2.397380e+02  |
| SSINE      | CUTEst | 3   | 2    | Optimal                | 0.000000e+00  | 0.000000e+00  |
| SYNTHES2   | CUTEst | 11  | 14   | Optimal                | -5.544022e-01 | -5.544063e-01 |
| TENBARS1   | CUTEst | 18  | 9    | Optimal                | 2.302549e+03  | 2.295373e+03  |
| THURBERLS  | CUTEst | 7   | 0    | Optimal                | 5.642708e+03  | 5.642708e+03  |
| TOINTGOR   | CUTEst | 50  | 0    | Optimal                | 1.373905e+03  | 1.373905e+03  |
| TOINTPSP   | CUTEst | 50  | 0    | Optimal                | 2.255604e+02  | 2.255604e+02  |
| TOINTQOR   | CUTEst | 50  | 0    | Optimal                | 1.175472e+03  | 1.175472e+03  |
| TRIGGER    | CUTEst | 7   | 6    | Optimal                | 0.000000e+00  | 0.000000e+00  |
| TRO6X2     | CUTEst | 45  | 21   | RestorationFailed      | 1.225000e+03  | -3.177966e+14 |
| VESUVIALS  | CUTEst | 8   | 0    | Optimal                | 1.500440e+03  | 9.914100e+02  |
| VESUVIOULS | CUTEst | 8   | 0    | Optimal                | 4.771193e-01  | 4.771138e-01  |
| WATER      | CUTEst | 31  | 10   | Optimal                | 1.054938e+04  | 1.054938e+04  |
| YFITNE     | CUTEst | 3   | 17   | IpoptStatus(-10)       | 0.000000e+00  | 0.000000e+00  |
| ZY2        | CUTEst | 3   | 2    | Optimal                | 2.000084e+00  | 2.000000e+00  |

---
*Generated by benchmark_report.py*