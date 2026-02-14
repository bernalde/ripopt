# CUTEst Benchmark Report

Comparison of ripopt vs Ipopt (C++) on the CUTEst test set.

## Executive Summary

- **Total problems**: 727
- **ripopt solved**: 586/727 (80.6%)
- **Ipopt solved**: 560/727 (77.0%)
- **Both solved**: 548/727
- **Matching solutions** (rel obj diff < 1e-4): 455/548

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (455 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 1.17e-06 |
| Median | 1.53e-13 |
| Max    | 8.77e-05 |

**All both-solved** (548 problems, including 93 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 8.23e-02 |
| Median | 3.10e-10 |
| Max    | 2.00e+00 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 493 | 367 | 343 | 333 | 268 |
| unconstrained | 234 | 219 | 217 | 215 | 187 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| 3PK | 30 | 0 | Optimal | Optimal | 1.42e-15 | 9 | 9 | 160us | 2.6ms | 16.0x | PASS |
| ACOPP14 | 38 | 68 | Optimal | Optimal | 9.79e-10 | 16 | 9 | 1.7ms | 4.0ms | 2.3x | PASS |
| ACOPP30 | 72 | 142 | Optimal | Optimal | 6.40e-09 | 47 | 13 | 10.0ms | 6.5ms | 0.7x | PASS |
| ACOPR14 | 38 | 82 | MaxIteration | Optimal | N/A | 2999 | 13 | 1.51s | 5.4ms | 0.0x | ripopt_FAIL |
| ACOPR30 | 72 | 172 | RestorationF | Optimal | N/A | 258 | 221 | 1.45s | 109.4ms | 0.1x | ripopt_FAIL |
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 52us | 1.3ms | 24.8x | PASS |
| AIRCRFTB | 8 | 0 | Optimal | Optimal | 4.50e-21 | 11 | 15 | 76us | 3.1ms | 40.9x | PASS |
| AIRPORT | 84 | 42 | Optimal | Optimal | 4.99e-09 | 12 | 13 | 3.7ms | 6.0ms | 1.6x | PASS |
| AKIVA | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 76us | 1.6ms | 20.6x | PASS |
| ALLINIT | 4 | 0 | Acceptable | Optimal | 3.69e-09 | 10 | 20 | 60us | 4.0ms | 66.1x | PASS |
| ALLINITA | 4 | 4 | Acceptable | Optimal | 9.53e-06 | 27 | 12 | 156us | 2.8ms | 17.8x | PASS |
| ALLINITC | 4 | 1 | Acceptable | Optimal | 4.03e-07 | 26 | 17 | 137us | 3.5ms | 25.5x | PASS |
| ALLINITU | 4 | 0 | Optimal | Optimal | 3.09e-16 | 8 | 14 | 52us | 2.9ms | 55.2x | PASS |
| ALSOTAME | 2 | 1 | Optimal | Optimal | 1.47e-08 | 10 | 8 | 57us | 2.2ms | 37.7x | PASS |
| ANTWERP | 27 | 10 | Optimal | Optimal | 7.28e-08 | 44 | 108 | 1.6ms | 22.8ms | 14.2x | PASS |
| ARGAUSS | 3 | 15 | Acceptable | IpoptStatus( | N/A | 2 | 0 | 74us | 984us | 13.3x | ipopt_FAIL |
| AVGASA | 8 | 10 | Acceptable | Optimal | 2.70e-01 | 11 | 9 | 109us | 2.5ms | 22.7x | MISMATCH |
| AVGASB | 8 | 10 | Acceptable | Optimal | 1.71e-01 | 19 | 11 | 148us | 2.7ms | 18.2x | MISMATCH |
| AVION2 | 49 | 15 | Acceptable | MaxIteration | N/A | 26 | 3000 | 3.2ms | 643.9ms | 202.0x | ipopt_FAIL |
| BA-L1 | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 6 | 624us | 2.1ms | 3.4x | PASS |
| BA-L1LS | 57 | 0 | Optimal | Optimal | 7.65e-21 | 7 | 10 | 663us | 2.9ms | 4.3x | PASS |
| BA-L1SP | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.4ms | 2.9ms | 2.1x | PASS |
| BA-L1SPLS | 57 | 0 | Optimal | Optimal | 6.48e-17 | 23 | 9 | 7.0ms | 4.9ms | 0.7x | PASS |
| BARD | 3 | 0 | Optimal | Optimal | 1.73e-18 | 8 | 8 | 50us | 1.9ms | 37.5x | PASS |
| BARDNE | 3 | 15 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 77us | 593us | 7.7x | BOTH_FAIL |
| BATCH | 48 | 73 | Acceptable | Optimal | 2.29e-05 | 201 | 29 | 10.2ms | 8.3ms | 0.8x | PASS |
| BEALE | 2 | 0 | Optimal | Optimal | 4.34e-18 | 7 | 8 | 68us | 2.2ms | 32.3x | PASS |
| BEALENE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 56us | 590us | 10.6x | ipopt_FAIL |
| BENNETT5 | 3 | 154 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 251.7ms | 663us | 0.0x | BOTH_FAIL |
| BENNETT5LS | 3 | 0 | Optimal | Optimal | 1.00e+00 | 14 | 21 | 584us | 4.5ms | 7.7x | MISMATCH |
| BIGGS3 | 6 | 0 | Optimal | Optimal | 1.30e-18 | 8 | 9 | 76us | 2.5ms | 32.5x | PASS |
| BIGGS5 | 6 | 0 | Optimal | Optimal | 4.74e-20 | 25 | 20 | 156us | 4.0ms | 25.9x | PASS |
| BIGGS6 | 6 | 0 | Optimal | Optimal | 4.32e-20 | 89 | 79 | 388us | 12.4ms | 31.9x | PASS |
| BIGGS6NE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 591 | 0 | 3.8ms | 586us | 0.2x | BOTH_FAIL |
| BIGGSC4 | 4 | 7 | Acceptable | Optimal | 8.72e-01 | 2999 | 17 | 10.3ms | 3.7ms | 0.4x | MISMATCH |
| BLEACHNG | 17 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| BOOTH | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 27us | 966us | 36.2x | PASS |
| BOX2 | 3 | 0 | Optimal | Optimal | 7.83e-18 | 9 | 8 | 64us | 1.9ms | 30.1x | PASS |
| BOX3 | 3 | 0 | Optimal | Optimal | 1.69e-25 | 11 | 9 | 62us | 2.2ms | 34.7x | PASS |
| BOX3NE | 3 | 10 | Optimal | IpoptStatus( | N/A | 11 | 0 | 90us | 573us | 6.4x | ipopt_FAIL |
| BOXBOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 3 | 0 | 43us | 594us | 13.8x | BOTH_FAIL |
| BOXBODLS | 2 | 0 | Optimal | Optimal | 4.40e-13 | 3 | 13 | 33us | 3.1ms | 96.1x | PASS |
| BQP1VAR | 1 | 0 | Optimal | Optimal | 9.98e-09 | 7 | 5 | 37us | 1.7ms | 45.1x | PASS |
| BQPGABIM | 50 | 0 | Acceptable | Optimal | 4.40e-06 | 11 | 12 | 371us | 3.1ms | 8.4x | PASS |
| BQPGASIM | 50 | 0 | Acceptable | Optimal | 6.13e-06 | 11 | 12 | 371us | 3.1ms | 8.5x | PASS |
| BRANIN | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 51us | 2.1ms | 41.2x | PASS |
| BRKMCC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 29us | 1.2ms | 42.0x | PASS |
| BROWNBS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 47us | 1.7ms | 37.4x | PASS |
| BROWNBSNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 8 | 0 | 57us | 594us | 10.4x | ipopt_FAIL |
| BROWNDEN | 4 | 0 | Optimal | Optimal | 1.70e-16 | 8 | 8 | 65us | 1.8ms | 27.5x | PASS |
| BROWNDENE | 4 | 20 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 94us | 581us | 6.2x | BOTH_FAIL |
| BT1 | 2 | 1 | Optimal | Optimal | 2.40e-09 | 17 | 7 | 104us | 1.9ms | 18.8x | PASS |
| BT10 | 2 | 2 | Optimal | Optimal | 2.79e-09 | 7 | 6 | 46us | 1.7ms | 37.0x | PASS |
| BT11 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 8 | 78us | 1.9ms | 24.1x | PASS |
| BT12 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 48us | 1.4ms | 30.1x | PASS |
| BT13 | 5 | 1 | Acceptable | Optimal | 1.00e-08 | 27 | 24 | 103us | 4.7ms | 46.2x | PASS |
| BT2 | 3 | 1 | Optimal | Optimal | 0.00e+00 | 12 | 12 | 62us | 2.3ms | 36.8x | PASS |
| BT3 | 5 | 3 | Optimal | Optimal | 1.22e-14 | 1 | 1 | 35us | 1.0ms | 29.4x | PASS |
| BT4 | 3 | 2 | Optimal | Optimal | 9.19e-01 | 6 | 9 | 62us | 2.2ms | 34.9x | MISMATCH |
| BT5 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 56us | 1.7ms | 30.7x | PASS |
| BT6 | 5 | 2 | Optimal | Optimal | 3.46e-12 | 9 | 13 | 71us | 2.5ms | 35.8x | PASS |
| BT7 | 5 | 3 | Optimal | Optimal | 1.50e-01 | 25 | 16 | 160us | 3.4ms | 21.3x | MISMATCH |
| BT8 | 5 | 2 | Acceptable | Optimal | 3.73e-09 | 32 | 14 | 133us | 2.5ms | 18.9x | PASS |
| BT9 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 80us | 2.5ms | 31.9x | PASS |
| BURKEHAN | 1 | 1 | RestorationF | Infeasible | N/A | 229 | 11 | 32.5ms | 3.3ms | 0.1x | BOTH_FAIL |
| BYRDSPHR | 3 | 2 | Optimal | Optimal | 3.95e-10 | 23 | 12 | 150us | 2.8ms | 18.6x | PASS |
| CAMEL6 | 2 | 0 | Optimal | Optimal | 7.91e-01 | 9 | 8 | 48us | 2.2ms | 47.3x | MISMATCH |
| CANTILVR | 5 | 1 | Optimal | Optimal | 3.23e-09 | 31 | 11 | 115us | 2.8ms | 24.1x | PASS |
| CB2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 8 | 8 | 57us | 2.2ms | 39.4x | MISMATCH |
| CB3 | 3 | 3 | Optimal | Optimal | 4.30e-09 | 8 | 8 | 55us | 2.3ms | 41.3x | PASS |
| CERI651A | 7 | 61 | LocalInfeasi | IpoptStatus( | N/A | 264 | 0 | 6.8ms | 610us | 0.1x | BOTH_FAIL |
| CERI651ALS | 7 | 0 | Acceptable | Optimal | 8.56e-08 | 283 | 95 | 4.9ms | 15.8ms | 3.2x | PASS |
| CERI651B | 7 | 66 | LocalInfeasi | IpoptStatus( | N/A | 83 | 0 | 2.3ms | 630us | 0.3x | BOTH_FAIL |
| CERI651BLS | 7 | 0 | Optimal | Optimal | 1.23e-08 | 89 | 56 | 1.6ms | 9.3ms | 6.0x | PASS |
| CERI651C | 7 | 56 | LocalInfeasi | IpoptStatus( | N/A | 171 | 0 | 4.0ms | 647us | 0.2x | BOTH_FAIL |
| CERI651CLS | 7 | 0 | Optimal | Optimal | 3.90e-09 | 178 | 53 | 2.6ms | 8.0ms | 3.0x | PASS |
| CERI651D | 7 | 67 | LocalInfeasi | IpoptStatus( | N/A | 91 | 0 | 2.5ms | 613us | 0.2x | BOTH_FAIL |
| CERI651DLS | 7 | 0 | Acceptable | Optimal | 3.08e-10 | 98 | 60 | 1.7ms | 10.7ms | 6.1x | PASS |
| CERI651E | 7 | 64 | LocalInfeasi | IpoptStatus( | N/A | 52 | 0 | 1.5ms | 641us | 0.4x | BOTH_FAIL |
| CERI651ELS | 7 | 0 | Acceptable | Optimal | 1.03e-09 | 83 | 45 | 1.4ms | 6.9ms | 5.0x | PASS |
| CHACONN1 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 5 | 6 | 50us | 1.9ms | 38.3x | MISMATCH |
| CHACONN2 | 3 | 3 | Optimal | Optimal | 4.49e-09 | 6 | 6 | 52us | 1.9ms | 36.9x | PASS |
| CHWIRUT1 | 3 | 214 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 1.1ms | 590us | 0.5x | BOTH_FAIL |
| CHWIRUT1LS | 3 | 0 | Acceptable | Optimal | 9.64e-01 | 21 | 6 | 563us | 2.0ms | 3.5x | MISMATCH |
| CHWIRUT2 | 3 | 54 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 341us | 596us | 1.7x | BOTH_FAIL |
| CHWIRUT2LS | 3 | 0 | Acceptable | Optimal | 9.62e-01 | 2999 | 6 | 35.0ms | 2.0ms | 0.1x | MISMATCH |
| CLIFF | 2 | 0 | Optimal | Optimal | 7.45e-03 | 27 | 23 | 66us | 3.4ms | 52.2x | MISMATCH |
| CLUSTER | 2 | 2 | Acceptable | Optimal | 0.00e+00 | 12 | 9 | 68us | 2.3ms | 33.7x | PASS |
| CLUSTERLS | 2 | 0 | Optimal | Optimal | 2.72e-18 | 13 | 17 | 52us | 2.9ms | 56.2x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 297us | 2.2ms | 7.5x | PASS |
| CONGIGMZ | 3 | 5 | Optimal | Optimal | 3.19e-09 | 7 | 20 | 135us | 4.1ms | 30.6x | PASS |
| COOLHANS | 9 | 9 | Optimal | Optimal | 0.00e+00 | 22 | 9 | 221us | 2.0ms | 9.2x | PASS |
| COOLHANSLS | 9 | 0 | Optimal | Optimal | 1.21e-18 | 23 | 25 | 170us | 4.3ms | 25.4x | PASS |
| CORE1 | 65 | 59 | Optimal | Optimal | 4.04e-09 | 24 | 33 | 51.2ms | 8.1ms | 0.2x | PASS |
| CRESC100 | 6 | 200 | Optimal | Infeasible | N/A | 199 | 155 | 96.0ms | 116.5ms | 1.2x | ipopt_FAIL |
| CRESC132 | 6 | 2654 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC4 | 6 | 8 | Optimal | Optimal | 1.35e-08 | 21 | 64 | 360us | 12.9ms | 35.9x | PASS |
| CRESC50 | 6 | 100 | RestorationF | Optimal | N/A | 91 | 194 | 269.0ms | 82.9ms | 0.3x | ripopt_FAIL |
| CSFI1 | 5 | 4 | Optimal | Optimal | 1.46e-08 | 26 | 11 | 144us | 2.9ms | 20.4x | PASS |
| CSFI2 | 5 | 4 | Optimal | Optimal | 1.48e-08 | 18 | 14 | 271us | 3.4ms | 12.7x | PASS |
| CUBE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 27 | 27 | 64us | 4.5ms | 70.7x | PASS |
| CUBENE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 27 | 1 | 96us | 999us | 10.5x | PASS |
| DALLASS | 46 | 31 | Optimal | Optimal | 1.41e-07 | 29 | 22 | 3.6ms | 5.7ms | 1.6x | PASS |
| DANIWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 73us | 607us | 8.3x | BOTH_FAIL |
| DANIWOODLS | 2 | 0 | Optimal | Optimal | 2.60e-18 | 10 | 10 | 49us | 2.3ms | 45.9x | PASS |
| DANWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 1 | 0 | 38us | 649us | 17.0x | BOTH_FAIL |
| DANWOODLS | 2 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 11 | 30us | 2.5ms | 83.9x | MISMATCH |
| DECONVB | 63 | 0 | Optimal | MaxIteration | N/A | 74 | 3000 | 9.4ms | 713.1ms | 75.5x | ipopt_FAIL |
| DECONVBNE | 63 | 40 | Acceptable | Optimal | 0.00e+00 | 253 | 505 | 40.1ms | 160.9ms | 4.0x | PASS |
| DECONVC | 63 | 1 | Acceptable | Optimal | 1.36e-03 | 37 | 31 | 5.9ms | 10.2ms | 1.7x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 2 | 26 | 532us | 26.8ms | 50.4x | PASS |
| DECONVU | 63 | 0 | Acceptable | Optimal | 7.56e-09 | 35 | 333 | 5.8ms | 85.7ms | 14.8x | PASS |
| DEGENLPA | 20 | 15 | Acceptable | Optimal | 1.81e-03 | 40 | 18 | 747us | 4.1ms | 5.5x | MISMATCH |
| DEGENLPB | 20 | 15 | Acceptable | Optimal | 2.35e-03 | 33 | 19 | 550us | 4.3ms | 7.9x | MISMATCH |
| DEMBO7 | 16 | 20 | Acceptable | Optimal | 1.03e-07 | 225 | 45 | 8.0ms | 9.3ms | 1.2x | PASS |
| DEMYMALO | 3 | 3 | Optimal | Optimal | 2.96e-09 | 12 | 9 | 69us | 2.6ms | 38.4x | PASS |
| DENSCHNA | 2 | 0 | Optimal | Optimal | 5.88e-39 | 6 | 6 | 37us | 1.5ms | 40.7x | PASS |
| DENSCHNB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 40us | 1.9ms | 47.1x | PASS |
| DENSCHNBNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 50us | 602us | 12.1x | ipopt_FAIL |
| DENSCHNC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 43us | 2.0ms | 47.5x | PASS |
| DENSCHNCNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 10 | 7 | 56us | 1.8ms | 31.7x | PASS |
| DENSCHND | 3 | 0 | Acceptable | Optimal | 2.22e-04 | 30 | 26 | 81us | 4.2ms | 52.5x | MISMATCH |
| DENSCHNDNE | 3 | 3 | Acceptable | Acceptable | 0.00e+00 | 19 | 22 | 166us | 3.7ms | 22.5x | PASS |
| DENSCHNE | 3 | 0 | Optimal | Optimal | 1.86e-17 | 10 | 14 | 46us | 3.1ms | 68.3x | PASS |
| DENSCHNENE | 3 | 3 | Optimal | Infeasible | N/A | 10 | 10 | 60us | 2.8ms | 47.1x | ipopt_FAIL |
| DENSCHNF | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 37us | 1.5ms | 42.2x | PASS |
| DENSCHNFNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 47us | 1.5ms | 32.1x | PASS |
| DEVGLA1 | 4 | 0 | Optimal | Optimal | 3.42e-23 | 21 | 23 | 205us | 4.3ms | 21.2x | PASS |
| DEVGLA1B | 4 | 0 | Acceptable | Optimal | 1.00e+00 | 28 | 20 | 236us | 5.0ms | 21.4x | MISMATCH |
| DEVGLA1NE | 4 | 24 | Optimal | IpoptStatus( | N/A | 16 | 0 | 243us | 590us | 2.4x | ipopt_FAIL |
| DEVGLA2 | 5 | 0 | Optimal | Optimal | 1.00e+00 | 16 | 13 | 185us | 2.7ms | 14.9x | MISMATCH |
| DEVGLA2B | 5 | 0 | Acceptable | Optimal | 2.58e-07 | 14 | 24 | 158us | 5.3ms | 33.2x | PASS |
| DEVGLA2NE | 5 | 16 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 286us | 597us | 2.1x | BOTH_FAIL |
| DGOSPEC | 3 | 0 | Acceptable | Optimal | 4.63e-03 | 10 | 27 | 56us | 5.7ms | 101.8x | MISMATCH |
| DIAMON2D | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3D | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIPIGRI | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 118us | 2.5ms | 21.3x | PASS |
| DISC2 | 29 | 23 | Optimal | Optimal | 9.11e-10 | 30 | 24 | 3.2ms | 6.2ms | 1.9x | PASS |
| DISCS | 36 | 66 | Optimal | Optimal | 1.11e-08 | 79 | 184 | 25.1ms | 70.6ms | 2.8x | PASS |
| DIXCHLNG | 10 | 5 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 188us | 2.3ms | 12.0x | PASS |
| DJTL | 2 | 0 | Acceptable | Acceptable | 0.00e+00 | 2999 | 1538 | 7.5ms | 159.3ms | 21.3x | PASS |
| DMN15102 | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15102LS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15103 | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15103LS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15332 | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15332LS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15333 | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN15333LS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN37142 | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN37142LS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN37143 | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DMN37143LS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DNIEPER | 61 | 24 | Acceptable | Optimal | 4.28e-08 | 338 | 23 | 45.8ms | 5.9ms | 0.1x | PASS |
| DUAL1 | 85 | 1 | Acceptable | Optimal | 3.88e-06 | 12 | 15 | 3.1ms | 6.7ms | 2.1x | PASS |
| DUAL2 | 96 | 1 | Acceptable | Optimal | 8.07e-09 | 12 | 12 | 4.2ms | 6.4ms | 1.5x | PASS |
| DUAL4 | 75 | 1 | Acceptable | Optimal | 1.40e-07 | 11 | 12 | 2.1ms | 5.0ms | 2.3x | PASS |
| DUALC1 | 9 | 215 | Acceptable | Optimal | 6.76e-06 | 47 | 18 | 10.6ms | 11.6ms | 1.1x | PASS |
| DUALC2 | 7 | 229 | Acceptable | Optimal | 1.62e-06 | 17 | 12 | 4.3ms | 8.2ms | 1.9x | PASS |
| DUALC5 | 8 | 278 | Acceptable | Optimal | 1.83e-07 | 12 | 11 | 4.2ms | 9.1ms | 2.2x | PASS |
| DUALC8 | 8 | 503 | Acceptable | Optimal | 1.27e-02 | 20 | 13 | 16.4ms | 15.4ms | 0.9x | MISMATCH |
| ECKERLE4 | 3 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 214us | 604us | 2.8x | BOTH_FAIL |
| ECKERLE4LS | 3 | 0 | Acceptable | Optimal | 4.97e-01 | 25 | 36 | 163us | 6.8ms | 41.5x | MISMATCH |
| EG1 | 3 | 0 | Optimal | Optimal | 2.07e-01 | 9 | 8 | 51us | 2.4ms | 48.0x | MISMATCH |
| EGGCRATE | 2 | 0 | Optimal | Optimal | 5.00e-01 | 7 | 5 | 48us | 1.5ms | 32.0x | MISMATCH |
| EGGCRATEB | 2 | 0 | Optimal | Optimal | 5.62e-16 | 10 | 6 | 46us | 1.9ms | 42.2x | PASS |
| EGGCRATENE | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 64us | 580us | 9.0x | BOTH_FAIL |
| ELATTAR | 7 | 102 | Acceptable | Optimal | 1.00e+00 | 65 | 81 | 7.5ms | 37.3ms | 5.0x | MISMATCH |
| ELATVIDU | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 11 | 45us | 2.3ms | 50.0x | PASS |
| ELATVIDUB | 2 | 0 | Optimal | Optimal | 1.03e-11 | 10 | 11 | 44us | 2.6ms | 58.3x | PASS |
| ELATVIDUNE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 61us | 598us | 9.8x | BOTH_FAIL |
| ENGVAL2 | 3 | 0 | Optimal | Optimal | 1.70e-20 | 20 | 21 | 76us | 3.8ms | 50.5x | PASS |
| ENGVAL2NE | 3 | 5 | Optimal | IpoptStatus( | N/A | 17 | 0 | 97us | 583us | 6.0x | ipopt_FAIL |
| ENSO | 9 | 168 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 1.2ms | 659us | 0.5x | BOTH_FAIL |
| ENSOLS | 9 | 0 | Acceptable | Optimal | 5.77e-16 | 13 | 7 | 1.7ms | 2.9ms | 1.7x | PASS |
| EQC | 9 | 3 | Acceptable | ErrorInStepC | N/A | 7 | 15 | 163us | 5.0ms | 30.3x | ipopt_FAIL |
| ERRINBAR | 18 | 9 | Acceptable | Optimal | 4.70e-07 | 68 | 37 | 782us | 7.7ms | 9.9x | PASS |
| EXP2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 44us | 1.8ms | 41.5x | PASS |
| EXP2B | 2 | 0 | Optimal | Optimal | 2.25e-15 | 9 | 7 | 50us | 2.0ms | 39.2x | PASS |
| EXP2NE | 2 | 10 | Optimal | IpoptStatus( | N/A | 7 | 0 | 70us | 603us | 8.7x | ipopt_FAIL |
| EXPFIT | 2 | 0 | Optimal | Optimal | 8.33e-17 | 5 | 8 | 45us | 2.2ms | 49.4x | PASS |
| EXPFITA | 5 | 22 | Optimal | Optimal | 3.60e-03 | 18 | 13 | 294us | 3.4ms | 11.6x | MISMATCH |
| EXPFITB | 5 | 102 | Optimal | Optimal | 1.92e-10 | 74 | 16 | 4.3ms | 5.9ms | 1.4x | PASS |
| EXPFITC | 5 | 502 | Optimal | Optimal | 5.30e-05 | 221 | 18 | 67.1ms | 17.0ms | 0.3x | PASS |
| EXPFITNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 68us | 607us | 9.0x | BOTH_FAIL |
| EXTRASIM | 2 | 1 | Optimal | Optimal | 1.82e-08 | 4 | 3 | 39us | 1.5ms | 37.5x | PASS |
| FBRAIN | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 5.2ms | 864us | 0.2x | BOTH_FAIL |
| FBRAIN2 | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 40.8ms | 1.1ms | 0.0x | BOTH_FAIL |
| FBRAIN2LS | 4 | 0 | Optimal | Optimal | 5.11e-10 | 10 | 10 | 7.5ms | 9.3ms | 1.2x | PASS |
| FBRAIN2NE | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 12.7ms | 1.2ms | 0.1x | BOTH_FAIL |
| FBRAIN3 | 6 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 5.29s | 1.3ms | 0.0x | BOTH_FAIL |
| FBRAIN3LS | 6 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 3.58s | 3.82s | 1.1x | BOTH_FAIL |
| FBRAINLS | 2 | 0 | Optimal | Optimal | 6.11e-16 | 9 | 7 | 3.7ms | 4.6ms | 1.2x | PASS |
| FBRAINNE | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 6.4ms | 919us | 0.1x | BOTH_FAIL |
| FCCU | 19 | 8 | Optimal | Optimal | 1.75e-15 | 10 | 9 | 146us | 2.5ms | 17.3x | PASS |
| FEEDLOC | 90 | 259 | Optimal | Optimal | 1.10e-08 | 38 | 23 | 88.0ms | 15.9ms | 0.2x | PASS |
| FLETCHER | 4 | 4 | Optimal | Optimal | 1.23e-08 | 31 | 28 | 170us | 5.9ms | 34.4x | PASS |
| FLT | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 56us | 1.9ms | 33.2x | PASS |
| GAUSS1 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 708us | 710us | 1.0x | BOTH_FAIL |
| GAUSS1LS | 8 | 0 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 443us | 1.8ms | 4.1x | PASS |
| GAUSS2 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 1.4ms | 640us | 0.5x | BOTH_FAIL |
| GAUSS2LS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 10 | 5 | 1.0ms | 1.9ms | 1.9x | PASS |
| GAUSS3 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 917us | 633us | 0.7x | BOTH_FAIL |
| GAUSS3LS | 8 | 0 | Optimal | Optimal | 3.65e-16 | 7 | 11 | 594us | 3.1ms | 5.3x | PASS |
| GAUSSIAN | 3 | 0 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 37us | 1.2ms | 31.5x | PASS |
| GBRAIN | 2 | 2200 | LocalInfeasi | IpoptStatus( | N/A | 6 | 0 | 4.7ms | 923us | 0.2x | BOTH_FAIL |
| GBRAINLS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 2.5ms | 3.8ms | 1.5x | PASS |
| GENHS28 | 10 | 8 | Optimal | Optimal | 1.22e-15 | 1 | 1 | 48us | 1.1ms | 22.6x | PASS |
| GIGOMEZ1 | 3 | 3 | Optimal | Optimal | 2.87e-09 | 10 | 13 | 69us | 3.2ms | 45.9x | PASS |
| GIGOMEZ2 | 3 | 3 | Optimal | Optimal | 5.23e-09 | 25 | 7 | 256us | 2.2ms | 8.6x | PASS |
| GIGOMEZ3 | 3 | 3 | Optimal | Optimal | 4.08e-09 | 10 | 8 | 70us | 2.3ms | 32.7x | PASS |
| GOFFIN | 51 | 50 | Acceptable | Optimal | 1.00e+00 | 2999 | 7 | 6.38s | 4.5ms | 0.0x | MISMATCH |
| GOTTFR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 5 | 65us | 1.7ms | 25.8x | PASS |
| GOULDQP1 | 32 | 17 | Acceptable | Optimal | 4.28e-07 | 34 | 15 | 1.2ms | 3.9ms | 3.3x | PASS |
| GROUPING | 100 | 125 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 2.1ms | 608us | 0.3x | ipopt_FAIL |
| GROWTH | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 262us | 580us | 2.2x | BOTH_FAIL |
| GROWTHLS | 3 | 0 | Optimal | Optimal | 1.77e-15 | 43 | 71 | 163us | 12.1ms | 74.3x | PASS |
| GULF | 3 | 0 | Optimal | Optimal | 3.04e-22 | 22 | 28 | 692us | 5.9ms | 8.5x | PASS |
| GULFNE | 3 | 99 | Optimal | IpoptStatus( | N/A | 22 | 0 | 1.1ms | 612us | 0.5x | ipopt_FAIL |
| HAHN1 | 7 | 236 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 284.5ms | 709us | 0.0x | BOTH_FAIL |
| HAHN1LS | 7 | 0 | Acceptable | Optimal | 7.55e-02 | 2999 | 78 | 224.7ms | 17.7ms | 0.1x | MISMATCH |
| HAIFAM | 99 | 150 | Acceptable | Optimal | 8.28e-07 | 27 | 40 | 6.0ms | 16.3ms | 2.7x | PASS |
| HAIFAS | 13 | 9 | Optimal | Optimal | 9.97e-09 | 43 | 16 | 503us | 4.4ms | 8.8x | PASS |
| HAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 24 | 62 | 82us | 10.8ms | 131.6x | PASS |
| HALDMADS | 6 | 42 | Optimal | Optimal | 1.00e+00 | 34 | 8 | 3.0ms | 3.3ms | 1.1x | MISMATCH |
| HART6 | 6 | 0 | Optimal | Optimal | 1.34e-16 | 10 | 7 | 68us | 2.4ms | 35.8x | PASS |
| HATFLDA | 4 | 0 | Optimal | Optimal | 1.58e-13 | 9 | 13 | 48us | 3.0ms | 62.1x | PASS |
| HATFLDANE | 4 | 4 | Acceptable | Optimal | 0.00e+00 | 9 | 6 | 58us | 2.0ms | 34.8x | PASS |
| HATFLDB | 4 | 0 | Optimal | Optimal | 3.90e-09 | 9 | 8 | 48us | 2.1ms | 44.7x | PASS |
| HATFLDBNE | 4 | 4 | MaxIteration | Infeasible | N/A | 2999 | 13 | 66.3ms | 3.8ms | 0.1x | BOTH_FAIL |
| HATFLDC | 25 | 0 | Optimal | Optimal | 6.53e-16 | 9 | 5 | 110us | 1.8ms | 16.2x | PASS |
| HATFLDCNE | 25 | 25 | Acceptable | Optimal | 0.00e+00 | 9 | 4 | 199us | 1.8ms | 9.0x | PASS |
| HATFLDD | 3 | 0 | Optimal | Optimal | 1.30e-19 | 21 | 21 | 89us | 3.6ms | 40.6x | PASS |
| HATFLDDNE | 3 | 10 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 129us | 602us | 4.7x | BOTH_FAIL |
| HATFLDE | 3 | 0 | Optimal | Optimal | 2.72e-20 | 19 | 20 | 105us | 3.2ms | 30.7x | PASS |
| HATFLDENE | 3 | 21 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 188us | 650us | 3.5x | BOTH_FAIL |
| HATFLDF | 3 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 135 | 92us | 24.1ms | 262.7x | PASS |
| HATFLDFL | 3 | 0 | Optimal | Optimal | 1.03e-08 | 494 | 1281 | 744us | 197.7ms | 265.8x | PASS |
| HATFLDFLNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 15 | 15 | 1.2ms | 3.7ms | 3.0x | PASS |
| HATFLDFLS | 3 | 0 | Optimal | Optimal | 3.79e-18 | 18 | 36 | 66us | 6.6ms | 99.4x | PASS |
| HATFLDG | 25 | 25 | Optimal | Optimal | 0.00e+00 | 13 | 7 | 290us | 2.1ms | 7.3x | PASS |
| HATFLDGLS | 25 | 0 | Optimal | Optimal | 1.37e-16 | 13 | 14 | 212us | 3.0ms | 14.0x | PASS |
| HATFLDH | 4 | 7 | MaxIteration | Optimal | N/A | 2999 | 17 | 9.9ms | 3.8ms | 0.4x | ripopt_FAIL |
| HEART6 | 6 | 6 | Optimal | Optimal | 0.00e+00 | 1014 | 22 | 5.3ms | 5.7ms | 1.1x | PASS |
| HEART6LS | 6 | 0 | Optimal | Optimal | 9.46e-23 | 1616 | 875 | 5.8ms | 140.7ms | 24.1x | PASS |
| HEART8 | 8 | 8 | Optimal | Optimal | 0.00e+00 | 66 | 12 | 391us | 3.1ms | 7.8x | PASS |
| HEART8LS | 8 | 0 | Optimal | Optimal | 3.37e-25 | 67 | 106 | 308us | 18.1ms | 58.7x | PASS |
| HELIX | 3 | 0 | Optimal | Optimal | 6.06e-25 | 8 | 13 | 58us | 2.8ms | 47.5x | PASS |
| HELIXNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 12 | 7 | 78us | 1.8ms | 23.4x | PASS |
| HET-Z | 2 | 1002 | Optimal | Optimal | 2.70e-02 | 65 | 11 | 42.7ms | 20.9ms | 0.5x | MISMATCH |
| HIELOW | 3 | 0 | Optimal | Optimal | 5.46e-15 | 7 | 8 | 10.0ms | 11.6ms | 1.2x | PASS |
| HIMMELBA | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 29us | 1.0ms | 35.4x | PASS |
| HIMMELBB | 2 | 0 | Optimal | Optimal | 1.40e-17 | 10 | 18 | 54us | 3.4ms | 62.5x | PASS |
| HIMMELBC | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 53us | 1.8ms | 33.6x | PASS |
| HIMMELBCLS | 2 | 0 | Optimal | Optimal | 5.80e-25 | 8 | 6 | 43us | 1.6ms | 37.0x | PASS |
| HIMMELBD | 2 | 2 | RestorationF | Infeasible | N/A | 14 | 22 | 883us | 5.2ms | 5.9x | BOTH_FAIL |
| HIMMELBE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 5 | 2 | 49us | 1.1ms | 23.6x | PASS |
| HIMMELBF | 4 | 0 | Acceptable | Optimal | 3.57e-15 | 46 | 75 | 162us | 11.1ms | 68.9x | PASS |
| HIMMELBFNE | 4 | 7 | LocalInfeasi | IpoptStatus( | N/A | 39 | 0 | 214us | 574us | 2.7x | BOTH_FAIL |
| HIMMELBG | 2 | 0 | Optimal | Optimal | 3.63e-22 | 6 | 6 | 40us | 1.8ms | 45.2x | PASS |
| HIMMELBH | 2 | 0 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 34us | 1.5ms | 43.2x | PASS |
| HIMMELBI | 100 | 12 | Optimal | Optimal | 5.80e-10 | 24 | 13 | 1.1ms | 4.0ms | 3.5x | PASS |
| HIMMELBJ | 45 | 14 | Acceptable | ErrorInStepC | N/A | 34 | 580 | 2.4ms | 138.3ms | 57.0x | ipopt_FAIL |
| HIMMELBK | 24 | 14 | Optimal | Optimal | 4.89e-08 | 19 | 18 | 624us | 4.3ms | 6.8x | PASS |
| HIMMELP1 | 2 | 0 | Optimal | Optimal | 1.83e-15 | 11 | 10 | 48us | 2.8ms | 59.7x | PASS |
| HIMMELP2 | 2 | 1 | Optimal | Optimal | 8.68e-01 | 12 | 17 | 70us | 4.1ms | 57.8x | MISMATCH |
| HIMMELP3 | 2 | 2 | Optimal | Optimal | 8.66e-01 | 0 | 11 | 24us | 2.8ms | 118.5x | MISMATCH |
| HIMMELP4 | 2 | 3 | Optimal | Optimal | 1.23e-08 | 9 | 23 | 71us | 5.1ms | 70.9x | PASS |
| HIMMELP5 | 2 | 3 | Optimal | Optimal | 7.50e-01 | 9 | 46 | 70us | 8.7ms | 123.0x | MISMATCH |
| HIMMELP6 | 2 | 5 | Optimal | Optimal | 7.50e-01 | 11 | 31 | 98us | 6.6ms | 67.7x | MISMATCH |
| HONG | 4 | 1 | Optimal | Optimal | 4.72e-16 | 9 | 7 | 61us | 2.0ms | 33.7x | PASS |
| HS1 | 2 | 0 | Optimal | Optimal | 9.22e-21 | 25 | 28 | 67us | 5.8ms | 86.2x | PASS |
| HS10 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 12 | 76us | 2.8ms | 37.2x | PASS |
| HS100 | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 122us | 2.7ms | 22.6x | PASS |
| HS100LNP | 7 | 2 | Optimal | Optimal | 1.67e-16 | 6 | 20 | 79us | 3.5ms | 44.6x | PASS |
| HS100MOD | 7 | 4 | Optimal | Optimal | 1.73e-11 | 9 | 14 | 114us | 3.6ms | 31.9x | PASS |
| HS101 | 7 | 5 | Optimal | Optimal | 4.60e-08 | 49 | 39 | 8.5ms | 9.7ms | 1.1x | PASS |
| HS102 | 7 | 5 | Optimal | Optimal | 4.27e-08 | 24 | 52 | 313us | 10.6ms | 33.9x | PASS |
| HS103 | 7 | 5 | Optimal | Optimal | 3.36e-08 | 21 | 21 | 263us | 4.9ms | 18.6x | PASS |
| HS104 | 8 | 5 | Optimal | Optimal | 2.69e-08 | 14 | 8 | 126us | 2.4ms | 18.7x | PASS |
| HS105 | 8 | 1 | Optimal | Optimal | 9.48e-12 | 30 | 23 | 3.5ms | 7.7ms | 2.2x | PASS |
| HS106 | 8 | 6 | Optimal | Optimal | 1.74e-08 | 18 | 18 | 126us | 4.1ms | 32.1x | PASS |
| HS107 | 9 | 6 | Acceptable | Optimal | 1.29e-07 | 52 | 7 | 359us | 2.3ms | 6.3x | PASS |
| HS108 | 9 | 13 | Optimal | Optimal | 8.77e-09 | 2596 | 11 | 15.2ms | 3.2ms | 0.2x | PASS |
| HS109 | 9 | 10 | MaxIteration | Optimal | N/A | 2999 | 14 | 22.8ms | 3.4ms | 0.1x | ripopt_FAIL |
| HS11 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 6 | 43us | 1.9ms | 44.3x | PASS |
| HS111 | 10 | 3 | Optimal | Optimal | 1.05e-11 | 11 | 15 | 155us | 3.6ms | 23.4x | PASS |
| HS111LNP | 10 | 3 | Optimal | Optimal | 1.03e-10 | 11 | 15 | 153us | 2.9ms | 18.7x | PASS |
| HS112 | 10 | 3 | Optimal | Optimal | 6.49e-14 | 10 | 10 | 145us | 2.6ms | 17.9x | PASS |
| HS113 | 10 | 8 | Optimal | Optimal | 1.68e-09 | 14 | 9 | 179us | 2.7ms | 15.3x | PASS |
| HS114 | 10 | 11 | Optimal | Optimal | 1.04e-07 | 19 | 13 | 206us | 3.5ms | 17.2x | PASS |
| HS116 | 13 | 14 | Acceptable | Optimal | 2.03e-04 | 150 | 19 | 1.9ms | 4.4ms | 2.3x | MISMATCH |
| HS117 | 15 | 5 | Acceptable | Optimal | 3.56e-07 | 38 | 19 | 380us | 4.5ms | 11.8x | PASS |
| HS118 | 15 | 17 | Optimal | Optimal | 1.29e-09 | 624 | 10 | 11.5ms | 2.9ms | 0.2x | PASS |
| HS119 | 16 | 8 | Acceptable | Optimal | 8.77e-05 | 21 | 17 | 399us | 4.0ms | 9.9x | PASS |
| HS12 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 7 | 6 | 57us | 1.9ms | 33.7x | PASS |
| HS13 | 2 | 1 | Acceptable | Optimal | 7.19e-04 | 24 | 47 | 84us | 8.6ms | 102.1x | MISMATCH |
| HS14 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 5 | 46us | 1.8ms | 40.3x | PASS |
| HS15 | 2 | 2 | Optimal | Optimal | 8.08e-08 | 11 | 13 | 85us | 3.0ms | 36.0x | PASS |
| HS16 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 10 | 10 | 69us | 2.7ms | 39.6x | MISMATCH |
| HS17 | 2 | 2 | Optimal | Optimal | 2.01e-04 | 12 | 22 | 71us | 4.5ms | 62.9x | MISMATCH |
| HS18 | 2 | 2 | Optimal | Optimal | 1.16e-09 | 15 | 10 | 81us | 2.6ms | 32.0x | PASS |
| HS19 | 2 | 2 | Optimal | Optimal | 3.09e-09 | 20 | 12 | 89us | 3.0ms | 34.1x | PASS |
| HS1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 25 | 30 | 86us | 6.9ms | 80.5x | PASS |
| HS2 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 10 | 10 | 43us | 2.6ms | 60.7x | PASS |
| HS20 | 2 | 3 | Optimal | Optimal | 6.17e-08 | 18 | 5 | 105us | 1.8ms | 17.3x | PASS |
| HS21 | 2 | 1 | Optimal | Optimal | 1.38e-10 | 9 | 6 | 52us | 1.9ms | 36.0x | PASS |
| HS21MOD | 7 | 1 | Acceptable | Optimal | 1.70e-08 | 11 | 13 | 73us | 3.1ms | 42.2x | PASS |
| HS22 | 2 | 2 | Optimal | Optimal | 1.16e-08 | 5 | 5 | 45us | 1.8ms | 40.5x | PASS |
| HS23 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 61 | 9 | 854us | 2.5ms | 2.9x | MISMATCH |
| HS24 | 2 | 3 | Optimal | Optimal | 1.80e-08 | 7 | 14 | 67us | 3.6ms | 53.2x | PASS |
| HS25 | 3 | 0 | Acceptable | Optimal | 6.08e-05 | 16 | 27 | 340us | 6.5ms | 19.1x | PASS |
| HS25NE | 3 | 99 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 539us | 612us | 1.1x | BOTH_FAIL |
| HS26 | 3 | 1 | Acceptable | Optimal | 2.94e-12 | 18 | 25 | 74us | 3.8ms | 51.9x | PASS |
| HS268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 46us | 3.4ms | 74.3x | MISMATCH |
| HS27 | 3 | 1 | Optimal | Optimal | 1.65e-14 | 12 | 57 | 85us | 9.4ms | 109.8x | PASS |
| HS28 | 3 | 1 | Optimal | Optimal | 9.24e-31 | 1 | 1 | 35us | 1.1ms | 30.7x | PASS |
| HS29 | 3 | 1 | Optimal | Optimal | 3.10e-10 | 18 | 7 | 93us | 2.1ms | 22.6x | PASS |
| HS2NE | 2 | 2 | RestorationF | Infeasible | N/A | 200 | 12 | 8.8ms | 3.6ms | 0.4x | BOTH_FAIL |
| HS3 | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 36us | 1.6ms | 44.9x | PASS |
| HS30 | 3 | 1 | Acceptable | Optimal | 1.45e-05 | 9 | 7 | 56us | 2.1ms | 38.1x | PASS |
| HS31 | 3 | 1 | Optimal | Optimal | 8.67e-09 | 11 | 6 | 63us | 2.1ms | 32.8x | PASS |
| HS32 | 3 | 2 | Acceptable | Optimal | 5.45e-05 | 12 | 15 | 81us | 3.6ms | 45.0x | PASS |
| HS33 | 3 | 2 | Acceptable | Optimal | 8.21e-08 | 14 | 9 | 95us | 2.5ms | 25.9x | PASS |
| HS34 | 3 | 2 | Optimal | Optimal | 9.46e-09 | 13 | 7 | 74us | 2.1ms | 28.1x | PASS |
| HS35 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 7 | 58us | 2.1ms | 35.8x | PASS |
| HS35I | 3 | 1 | Optimal | Optimal | 3.44e-09 | 9 | 7 | 58us | 2.1ms | 36.5x | PASS |
| HS35MOD | 3 | 1 | Acceptable | Optimal | 4.67e-07 | 8 | 14 | 56us | 3.2ms | 56.2x | PASS |
| HS36 | 3 | 1 | Optimal | Optimal | 6.33e-09 | 16 | 11 | 76us | 3.1ms | 41.1x | PASS |
| HS37 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 10 | 11 | 83us | 3.2ms | 38.6x | PASS |
| HS38 | 4 | 0 | Optimal | Optimal | 4.52e-11 | 37 | 39 | 105us | 7.8ms | 74.5x | PASS |
| HS39 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 79us | 2.4ms | 31.1x | PASS |
| HS3MOD | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 35us | 1.6ms | 44.9x | PASS |
| HS4 | 2 | 0 | Optimal | Optimal | 1.97e-08 | 6 | 4 | 37us | 1.6ms | 41.8x | PASS |
| HS40 | 4 | 3 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 46us | 1.3ms | 28.4x | PASS |
| HS41 | 4 | 1 | Optimal | Optimal | 1.15e-09 | 11 | 7 | 77us | 2.2ms | 28.6x | PASS |
| HS42 | 4 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 47us | 1.4ms | 30.6x | PASS |
| HS43 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 8 | 8 | 80us | 2.4ms | 29.7x | PASS |
| HS44 | 4 | 6 | Acceptable | Optimal | 2.97e-06 | 28 | 24 | 156us | 5.4ms | 34.9x | PASS |
| HS44NEW | 4 | 6 | Acceptable | Optimal | 1.15e-05 | 19 | 18 | 139us | 4.3ms | 31.2x | PASS |
| HS45 | 5 | 0 | Acceptable | Optimal | 3.14e-04 | 31 | 11 | 89us | 2.9ms | 32.5x | MISMATCH |
| HS46 | 5 | 2 | Acceptable | Optimal | 3.43e-11 | 17 | 19 | 96us | 3.1ms | 32.1x | PASS |
| HS47 | 5 | 3 | Acceptable | Optimal | 3.22e-12 | 17 | 19 | 98us | 3.2ms | 32.6x | PASS |
| HS48 | 5 | 2 | Optimal | Optimal | 4.44e-31 | 1 | 1 | 38us | 1.0ms | 27.4x | PASS |
| HS49 | 5 | 2 | Acceptable | Optimal | 2.61e-10 | 17 | 19 | 84us | 3.2ms | 37.7x | PASS |
| HS5 | 2 | 0 | Optimal | Optimal | 2.32e-16 | 9 | 7 | 47us | 2.1ms | 44.1x | PASS |
| HS50 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 75us | 2.1ms | 27.6x | PASS |
| HS51 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 38us | 1.0ms | 26.7x | PASS |
| HS52 | 5 | 3 | Optimal | Optimal | 1.17e-15 | 1 | 1 | 39us | 1.0ms | 26.6x | PASS |
| HS53 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 6 | 71us | 1.9ms | 26.6x | PASS |
| HS54 | 6 | 1 | Optimal | Optimal | 7.51e-01 | 10 | 15 | 81us | 3.5ms | 43.7x | MISMATCH |
| HS55 | 6 | 6 | Optimal | Optimal | 2.04e-02 | 11 | 18 | 92us | 5.1ms | 55.5x | MISMATCH |
| HS56 | 7 | 4 | Optimal | Optimal | 1.99e-14 | 5 | 10 | 70us | 2.3ms | 32.8x | PASS |
| HS57 | 2 | 1 | Optimal | Optimal | 6.11e-15 | 10 | 10 | 77us | 2.3ms | 29.8x | PASS |
| HS59 | 2 | 3 | Acceptable | Optimal | 9.28e-04 | 10 | 17 | 98us | 4.5ms | 45.5x | MISMATCH |
| HS6 | 2 | 1 | Optimal | Optimal | 4.93e-32 | 7 | 5 | 61us | 1.8ms | 28.8x | PASS |
| HS60 | 3 | 1 | Optimal | Optimal | 1.19e-13 | 8 | 6 | 57us | 2.0ms | 34.6x | PASS |
| HS61 | 3 | 2 | Optimal | Optimal | 7.91e-16 | 9 | 10 | 65us | 2.1ms | 31.8x | PASS |
| HS62 | 3 | 1 | Optimal | Optimal | 1.38e-16 | 9 | 6 | 67us | 2.2ms | 33.2x | PASS |
| HS63 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 5 | 65us | 1.8ms | 27.8x | PASS |
| HS64 | 3 | 1 | Optimal | Optimal | 3.62e-09 | 18 | 16 | 81us | 3.5ms | 43.9x | PASS |
| HS65 | 3 | 1 | Optimal | Optimal | 6.05e-09 | 15 | 16 | 85us | 3.9ms | 46.2x | PASS |
| HS66 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 10 | 10 | 66us | 2.5ms | 37.7x | PASS |
| HS67 | 3 | 14 | Optimal | Optimal | 8.73e-02 | 17 | 9 | 252us | 2.7ms | 10.5x | MISMATCH |
| HS68 | 4 | 2 | Optimal | Optimal | 5.21e-11 | 18 | 16 | 113us | 3.8ms | 33.7x | PASS |
| HS69 | 4 | 2 | Optimal | Optimal | 2.85e-15 | 11 | 10 | 84us | 2.7ms | 32.3x | PASS |
| HS7 | 2 | 1 | Optimal | Optimal | 5.30e-12 | 10 | 27 | 68us | 4.8ms | 71.4x | PASS |
| HS70 | 4 | 1 | Optimal | Optimal | 1.80e-01 | 10 | 46 | 176us | 8.7ms | 49.3x | MISMATCH |
| HS71 | 4 | 2 | Optimal | Optimal | 1.01e-09 | 10 | 8 | 72us | 2.4ms | 33.4x | PASS |
| HS72 | 4 | 2 | Optimal | Optimal | 6.77e-07 | 18 | 16 | 89us | 3.4ms | 38.8x | PASS |
| HS73 | 4 | 3 | Optimal | Optimal | 1.30e-09 | 11 | 8 | 76us | 2.4ms | 31.6x | PASS |
| HS74 | 4 | 5 | Optimal | Optimal | 0.00e+00 | 11 | 8 | 89us | 2.5ms | 28.2x | PASS |
| HS75 | 4 | 5 | Optimal | Optimal | 2.54e-09 | 12 | 8 | 92us | 2.3ms | 25.6x | PASS |
| HS76 | 4 | 3 | Optimal | Optimal | 5.89e-09 | 9 | 7 | 68us | 2.2ms | 32.3x | PASS |
| HS76I | 4 | 3 | Optimal | Optimal | 3.59e-09 | 10 | 6 | 74us | 2.0ms | 26.7x | PASS |
| HS77 | 5 | 2 | Optimal | Optimal | 1.99e-11 | 9 | 11 | 72us | 2.3ms | 31.9x | PASS |
| HS78 | 5 | 3 | Optimal | Optimal | 1.52e-16 | 4 | 4 | 49us | 1.5ms | 30.0x | PASS |
| HS79 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 54us | 1.4ms | 26.1x | PASS |
| HS8 | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 47us | 1.5ms | 33.0x | PASS |
| HS80 | 5 | 3 | Optimal | Optimal | 6.75e-13 | 9 | 5 | 81us | 1.8ms | 22.4x | PASS |
| HS81 | 5 | 3 | Optimal | Optimal | 3.46e-14 | 31 | 68 | 196us | 13.3ms | 67.5x | PASS |
| HS83 | 5 | 3 | MaxIteration | Optimal | N/A | 2999 | 9 | 11.3ms | 2.4ms | 0.2x | ripopt_FAIL |
| HS84 | 5 | 3 | Acceptable | Optimal | 1.04e-05 | 2999 | 9 | 183.8ms | 2.7ms | 0.0x | PASS |
| HS85 | 5 | 21 | Acceptable | Optimal | 1.22e-03 | 2999 | 13 | 1.94s | 4.6ms | 0.0x | MISMATCH |
| HS86 | 5 | 10 | Optimal | Optimal | 2.58e-09 | 11 | 10 | 131us | 2.6ms | 19.9x | PASS |
| HS87 | 6 | 4 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 12.8ms | 474.8ms | 37.2x | BOTH_FAIL |
| HS88 | 2 | 1 | Optimal | Optimal | 6.34e-06 | 24 | 18 | 1.0ms | 4.4ms | 4.4x | PASS |
| HS89 | 3 | 1 | Optimal | Optimal | 3.42e-06 | 45 | 15 | 4.6ms | 4.3ms | 0.9x | PASS |
| HS9 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 41us | 1.3ms | 32.6x | PASS |
| HS90 | 4 | 1 | Optimal | Optimal | 5.20e-06 | 18 | 16 | 1.4ms | 4.8ms | 3.5x | PASS |
| HS91 | 5 | 1 | Optimal | Optimal | 6.44e-06 | 18 | 16 | 1.8ms | 5.1ms | 2.8x | PASS |
| HS92 | 6 | 1 | Optimal | Optimal | 6.80e-06 | 16 | 35 | 1.8ms | 10.7ms | 6.0x | PASS |
| HS93 | 6 | 2 | Optimal | Optimal | 9.89e-09 | 14 | 7 | 113us | 2.5ms | 22.0x | PASS |
| HS95 | 6 | 4 | Acceptable | Optimal | 6.92e-04 | 46 | 9 | 199us | 2.5ms | 12.6x | MISMATCH |
| HS96 | 6 | 4 | Acceptable | Optimal | 6.88e-04 | 46 | 8 | 198us | 2.3ms | 11.8x | MISMATCH |
| HS97 | 6 | 4 | Acceptable | Optimal | 1.05e-06 | 48 | 24 | 282us | 5.2ms | 18.5x | PASS |
| HS98 | 6 | 4 | Acceptable | Optimal | 1.25e-04 | 156 | 13 | 552us | 3.1ms | 5.6x | MISMATCH |
| HS99 | 7 | 2 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 75us | 1.8ms | 24.2x | PASS |
| HS99EXP | 31 | 21 | Acceptable | Optimal | 1.78e-14 | 21 | 17 | 901us | 3.9ms | 4.3x | PASS |
| HUBFIT | 2 | 1 | Optimal | Optimal | 2.24e-09 | 7 | 7 | 47us | 2.0ms | 43.0x | PASS |
| HUMPS | 2 | 0 | Optimal | Optimal | 2.76e-17 | 40 | 1533 | 138us | 215.4ms | 1560.7x | PASS |
| HYDC20LS | 99 | 0 | Acceptable | Optimal | 2.98e-01 | 2999 | 639 | 998.0ms | 172.1ms | 0.2x | MISMATCH |
| HYDCAR20 | 99 | 99 | Optimal | Optimal | 0.00e+00 | 11 | 9 | 1.70s | 4.0ms | 0.0x | PASS |
| HYDCAR6 | 29 | 29 | Optimal | Optimal | 0.00e+00 | 1298 | 5 | 39.6ms | 2.0ms | 0.1x | PASS |
| HYDCAR6LS | 29 | 0 | Optimal | Optimal | 2.86e-18 | 1322 | 149 | 28.1ms | 30.0ms | 1.1x | PASS |
| HYPCIR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 53us | 1.7ms | 32.9x | PASS |
| JENSMP | 2 | 0 | Optimal | Optimal | 3.43e-16 | 10 | 9 | 52us | 2.0ms | 38.2x | PASS |
| JENSMPNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 80us | 616us | 7.7x | BOTH_FAIL |
| JUDGE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 53us | 2.0ms | 38.0x | PASS |
| JUDGEB | 2 | 0 | Optimal | Optimal | 2.21e-16 | 9 | 9 | 55us | 2.3ms | 42.8x | PASS |
| JUDGENE | 2 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 81us | 641us | 7.9x | BOTH_FAIL |
| KIRBY2 | 5 | 151 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 993us | 612us | 0.6x | BOTH_FAIL |
| KIRBY2LS | 5 | 0 | Acceptable | Optimal | 1.21e-14 | 26 | 11 | 675us | 2.7ms | 3.9x | PASS |
| KIWCRESC | 3 | 2 | Optimal | Optimal | 1.07e-08 | 13 | 8 | 78us | 2.5ms | 32.6x | PASS |
| KOEBHELB | 3 | 0 | Optimal | Optimal | 7.33e-16 | 843 | 71 | 12.4ms | 15.7ms | 1.3x | PASS |
| KOEBHELBNE | 3 | 156 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 2.3ms | 638us | 0.3x | BOTH_FAIL |
| KOWOSB | 4 | 0 | Optimal | Optimal | 5.96e-19 | 7 | 8 | 53us | 2.2ms | 42.1x | PASS |
| KOWOSBNE | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 77us | 579us | 7.5x | BOTH_FAIL |
| KSIP | 20 | 1001 | Optimal | Optimal | 9.79e-01 | 1 | 22 | 13.8ms | 71.3ms | 5.2x | MISMATCH |
| LAKES | 90 | 78 | Optimal | Optimal | 1.53e-13 | 13 | 11 | 1.2ms | 4.0ms | 3.2x | PASS |
| LANCZOS1 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 41 | 0 | 518us | 629us | 1.2x | ipopt_FAIL |
| LANCZOS1LS | 6 | 0 | Acceptable | Optimal | 5.33e-09 | 105 | 115 | 680us | 20.1ms | 29.5x | PASS |
| LANCZOS2 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 70 | 0 | 713us | 585us | 0.8x | ipopt_FAIL |
| LANCZOS2LS | 6 | 0 | Acceptable | Optimal | 7.41e-09 | 97 | 101 | 621us | 17.4ms | 27.9x | PASS |
| LANCZOS3 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 30 | 0 | 357us | 613us | 1.7x | ipopt_FAIL |
| LANCZOS3LS | 6 | 0 | Acceptable | Optimal | 1.08e-08 | 78 | 174 | 540us | 30.8ms | 57.0x | PASS |
| LAUNCH | 25 | 28 | Acceptable | Optimal | 3.45e-04 | 2999 | 12 | 179.7ms | 3.5ms | 0.0x | MISMATCH |
| LEVYMONE10 | 10 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 111us | 569us | 5.1x | BOTH_FAIL |
| LEVYMONE5 | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 72us | 618us | 8.5x | BOTH_FAIL |
| LEVYMONE6 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 79us | 576us | 7.3x | BOTH_FAIL |
| LEVYMONE7 | 4 | 8 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 84us | 589us | 7.0x | BOTH_FAIL |
| LEVYMONE8 | 5 | 10 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 77us | 601us | 7.8x | BOTH_FAIL |
| LEVYMONE9 | 8 | 16 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 103us | 599us | 5.8x | BOTH_FAIL |
| LEVYMONT10 | 10 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 4 | 74us | 1.6ms | 21.7x | PASS |
| LEVYMONT5 | 2 | 0 | Acceptable | Optimal | 1.00e+00 | 10 | 10 | 56us | 2.7ms | 47.5x | MISMATCH |
| LEVYMONT6 | 3 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 8 | 61us | 2.6ms | 42.9x | PASS |
| LEVYMONT7 | 4 | 0 | Optimal | Optimal | 1.42e-16 | 10 | 7 | 58us | 2.4ms | 40.7x | PASS |
| LEVYMONT8 | 5 | 0 | Optimal | Optimal | 1.64e-16 | 8 | 4 | 55us | 1.6ms | 29.0x | PASS |
| LEVYMONT9 | 8 | 0 | Optimal | Optimal | 1.71e-16 | 8 | 4 | 57us | 1.6ms | 28.4x | PASS |
| LEWISPOL | 6 | 9 | Acceptable | IpoptStatus( | N/A | 10 | 0 | 986us | 583us | 0.6x | ipopt_FAIL |
| LHAIFAM | 99 | 150 | MaxIteration | InvalidNumbe | N/A | 2999 | 0 | 908.4ms | 740us | 0.0x | BOTH_FAIL |
| LIN | 4 | 2 | Optimal | Optimal | 2.03e-03 | 9 | 7 | 79us | 2.0ms | 25.4x | MISMATCH |
| LINSPANH | 97 | 33 | Acceptable | Optimal | 5.76e-07 | 2999 | 24 | 115.7ms | 6.3ms | 0.1x | PASS |
| LOADBAL | 31 | 31 | Optimal | Optimal | 6.55e-09 | 15 | 13 | 883us | 3.6ms | 4.1x | PASS |
| LOGHAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 55 | 2747 | 230us | 402.3ms | 1752.6x | PASS |
| LOGROS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 50 | 49 | 106us | 10.5ms | 98.7x | PASS |
| LOOTSMA | 3 | 2 | Optimal | Optimal | 8.76e-08 | 26 | 13 | 127us | 3.4ms | 26.7x | PASS |
| LOTSCHD | 12 | 7 | Optimal | Optimal | 4.47e-10 | 15 | 9 | 155us | 2.6ms | 17.0x | PASS |
| LRCOVTYPE | 54 | 0 | Optimal | Optimal | 1.78e-02 | 65 | 33 | 14.73s | 6.45s | 0.4x | MISMATCH |
| LRIJCNN1 | 22 | 0 | Optimal | Optimal | 1.16e-14 | 18 | 11 | 355.3ms | 187.7ms | 0.5x | PASS |
| LSC1 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 97us | 586us | 6.1x | BOTH_FAIL |
| LSC1LS | 3 | 0 | Acceptable | Optimal | 2.30e-15 | 17 | 16 | 88us | 3.6ms | 40.6x | PASS |
| LSC2 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 142us | 573us | 4.0x | BOTH_FAIL |
| LSC2LS | 3 | 0 | Acceptable | Optimal | 2.32e-04 | 32 | 38 | 101us | 5.7ms | 56.8x | MISMATCH |
| LSNNODOC | 5 | 4 | Acceptable | Optimal | 3.42e-06 | 13 | 10 | 95us | 2.7ms | 28.0x | PASS |
| LSQFIT | 2 | 1 | Optimal | Optimal | 4.24e-09 | 7 | 7 | 54us | 2.4ms | 44.0x | PASS |
| MADSEN | 3 | 6 | Optimal | Optimal | 2.67e-01 | 20 | 18 | 192us | 4.5ms | 23.3x | MISMATCH |
| MAKELA1 | 3 | 2 | Optimal | Optimal | 2.00e+00 | 6 | 12 | 59us | 3.1ms | 51.5x | MISMATCH |
| MAKELA2 | 3 | 3 | Optimal | Optimal | 1.27e-01 | 3 | 6 | 42us | 1.9ms | 44.6x | MISMATCH |
| MAKELA3 | 21 | 20 | Acceptable | Optimal | 8.53e-09 | 25 | 11 | 491us | 3.1ms | 6.2x | PASS |
| MAKELA4 | 21 | 40 | Optimal | Optimal | 6.88e-01 | 1 | 5 | 103us | 2.0ms | 19.4x | MISMATCH |
| MARATOS | 2 | 1 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 36us | 1.4ms | 37.9x | PASS |
| MARATOSB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 670 | 672 | 808us | 100.3ms | 124.1x | PASS |
| MATRIX2 | 6 | 2 | Acceptable | Optimal | 7.84e-10 | 12 | 42 | 83us | 7.5ms | 91.1x | PASS |
| MAXLIKA | 8 | 0 | Acceptable | Optimal | 1.13e-02 | 20 | 23 | 1.9ms | 7.1ms | 3.6x | MISMATCH |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 289us | 2.1ms | 7.3x | PASS |
| MDHOLE | 2 | 0 | Optimal | Optimal | 9.98e-09 | 35 | 42 | 80us | 9.0ms | 112.2x | PASS |
| MESH | 41 | 48 | Acceptable | IpoptStatus( | N/A | 71 | 79 | 17.8ms | 22.0ms | 1.2x | ipopt_FAIL |
| METHANB8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 9 | 3 | 365us | 1.5ms | 4.2x | PASS |
| METHANB8LS | 31 | 0 | Optimal | Optimal | 5.35e-26 | 9 | 8 | 252us | 2.3ms | 9.2x | PASS |
| METHANL8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 608 | 4 | 19.4ms | 1.7ms | 0.1x | PASS |
| METHANL8LS | 31 | 0 | Optimal | Optimal | 5.70e-17 | 977 | 40 | 21.6ms | 8.6ms | 0.4x | PASS |
| MEXHAT | 2 | 0 | Optimal | Optimal | 6.77e-10 | 28 | 26 | 72us | 4.2ms | 58.3x | PASS |
| MEYER3 | 3 | 0 | Acceptable | Optimal | 2.23e-12 | 2999 | 194 | 10.8ms | 30.2ms | 2.8x | PASS |
| MEYER3NE | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 16.6ms | 590us | 0.0x | BOTH_FAIL |
| MGH09 | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 301us | 575us | 1.9x | BOTH_FAIL |
| MGH09LS | 4 | 0 | Optimal | Optimal | 5.96e-19 | 58 | 72 | 192us | 11.9ms | 62.1x | PASS |
| MGH10 | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 617 | 0 | 3.0ms | 597us | 0.2x | BOTH_FAIL |
| MGH10LS | 3 | 0 | Acceptable | Optimal | 4.94e-12 | 2999 | 1828 | 10.4ms | 278.4ms | 26.8x | PASS |
| MGH10S | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 13.4ms | 654us | 0.0x | BOTH_FAIL |
| MGH10SLS | 3 | 0 | MaxIteration | Optimal | N/A | 2999 | 354 | 8.6ms | 53.3ms | 6.2x | ripopt_FAIL |
| MGH17 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 381us | 581us | 1.5x | BOTH_FAIL |
| MGH17LS | 5 | 0 | Acceptable | Optimal | 6.06e-07 | 37 | 47 | 295us | 9.5ms | 32.3x | PASS |
| MGH17S | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 755us | 595us | 0.8x | BOTH_FAIL |
| MGH17SLS | 5 | 0 | Optimal | Optimal | 2.45e-02 | 54 | 41 | 387us | 8.1ms | 20.9x | MISMATCH |
| MIFFLIN1 | 3 | 2 | Optimal | Optimal | 9.24e-09 | 6 | 5 | 57us | 1.7ms | 30.3x | PASS |
| MIFFLIN2 | 3 | 2 | Optimal | Optimal | 9.97e-09 | 15 | 11 | 76us | 2.8ms | 37.6x | PASS |
| MINMAXBD | 5 | 20 | Optimal | Optimal | 8.55e-11 | 186 | 25 | 3.6ms | 6.3ms | 1.8x | PASS |
| MINMAXRB | 3 | 4 | Optimal | Optimal | 9.85e-09 | 3 | 8 | 47us | 2.2ms | 47.6x | PASS |
| MINSURF | 64 | 0 | Acceptable | Optimal | 0.00e+00 | 5 | 4 | 387us | 1.8ms | 4.7x | PASS |
| MISRA1A | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 174us | 588us | 3.4x | BOTH_FAIL |
| MISRA1ALS | 2 | 0 | Acceptable | Optimal | 1.91e-14 | 35 | 40 | 126us | 7.0ms | 55.8x | PASS |
| MISRA1B | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 178us | 632us | 3.5x | BOTH_FAIL |
| MISRA1BLS | 2 | 0 | Optimal | Optimal | 4.76e-14 | 25 | 34 | 94us | 5.8ms | 62.1x | PASS |
| MISRA1C | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 132us | 579us | 4.4x | BOTH_FAIL |
| MISRA1CLS | 2 | 0 | Acceptable | Optimal | 4.63e-14 | 18 | 14 | 92us | 3.2ms | 34.2x | PASS |
| MISRA1D | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 24 | 0 | 146us | 609us | 4.2x | BOTH_FAIL |
| MISRA1DLS | 2 | 0 | Optimal | Optimal | 9.02e-17 | 20 | 30 | 92us | 5.5ms | 59.2x | PASS |
| MISTAKE | 9 | 13 | Acceptable | Optimal | 5.00e-01 | 81 | 16 | 1.0ms | 4.2ms | 4.2x | MISMATCH |
| MRIBASIS | 36 | 55 | Acceptable | Optimal | 1.00e-08 | 36 | 15 | 6.8ms | 5.0ms | 0.7x | PASS |
| MSS1 | 90 | 73 | Acceptable | Optimal | 1.25e-01 | 72 | 95 | 38.0ms | 51.2ms | 1.3x | MISMATCH |
| MUONSINE | 1 | 512 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 2.2ms | 679us | 0.3x | BOTH_FAIL |
| MUONSINELS | 1 | 0 | Acceptable | Optimal | 1.36e-01 | 11 | 8 | 1.2ms | 2.2ms | 1.8x | MISMATCH |
| MWRIGHT | 5 | 3 | Optimal | Optimal | 9.48e-01 | 13 | 10 | 82us | 2.2ms | 26.3x | MISMATCH |
| NASH | 72 | 24 | RestorationF | Infeasible | N/A | 47 | 45 | 23.8ms | 12.5ms | 0.5x | BOTH_FAIL |
| NELSON | 3 | 128 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 102.6ms | 660us | 0.0x | BOTH_FAIL |
| NET1 | 48 | 57 | Acceptable | Optimal | 1.19e-07 | 2999 | 26 | 159.9ms | 6.9ms | 0.0x | PASS |
| NYSTROM5 | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 367us | 585us | 1.6x | ipopt_FAIL |
| NYSTROM5C | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 334us | 516us | 1.5x | ipopt_FAIL |
| ODFITS | 10 | 6 | Optimal | Optimal | 1.91e-16 | 11 | 8 | 107us | 2.3ms | 21.9x | PASS |
| OET1 | 3 | 1002 | Optimal | Optimal | 2.80e-01 | 32 | 33 | 26.5ms | 48.5ms | 1.8x | MISMATCH |
| OET2 | 3 | 1002 | MaxIteration | Optimal | N/A | 2999 | 181 | 1.55s | 275.8ms | 0.2x | ripopt_FAIL |
| OET3 | 4 | 1002 | Optimal | Optimal | 1.00e+00 | 5 | 13 | 10.4ms | 21.9ms | 2.1x | MISMATCH |
| OET4 | 4 | 1002 | Optimal | Optimal | 9.15e-01 | 25 | 165 | 28.0ms | 251.3ms | 9.0x | MISMATCH |
| OET5 | 5 | 1002 | Optimal | Optimal | 1.00e+00 | 211 | 64 | 33.80s | 119.1ms | 0.0x | MISMATCH |
| OET6 | 5 | 1002 | MaxIteration | Optimal | N/A | 2999 | 126 | 2.76s | 378.4ms | 0.1x | ripopt_FAIL |
| OET7 | 7 | 1002 | MaxIteration | Optimal | N/A | 2999 | 193 | 3.56s | 545.5ms | 0.2x | ripopt_FAIL |
| OPTCNTRL | 32 | 20 | Acceptable | Optimal | 2.02e-08 | 40 | 9 | 1.3ms | 2.5ms | 2.0x | PASS |
| OPTPRLOC | 30 | 30 | Optimal | Optimal | 1.22e-08 | 57 | 13 | 2.8ms | 4.0ms | 1.4x | PASS |
| ORTHREGB | 27 | 6 | Optimal | Optimal | 4.20e-19 | 2 | 2 | 99us | 1.3ms | 13.4x | PASS |
| OSBORNE1 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 2982 | 0 | 30.7ms | 682us | 0.0x | BOTH_FAIL |
| OSBORNE2 | 11 | 65 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 694us | 590us | 0.8x | BOTH_FAIL |
| OSBORNEA | 5 | 0 | MaxIteration | Optimal | N/A | 2999 | 64 | 18.5ms | 11.5ms | 0.6x | ripopt_FAIL |
| OSBORNEB | 11 | 0 | Optimal | Optimal | 4.86e-17 | 16 | 19 | 599us | 4.0ms | 6.7x | PASS |
| OSLBQP | 8 | 0 | Acceptable | Optimal | 7.24e-07 | 13 | 15 | 58us | 3.2ms | 54.5x | PASS |
| PALMER1 | 4 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 13 | 89us | 3.0ms | 33.2x | PASS |
| PALMER1A | 6 | 0 | Optimal | Optimal | 1.98e-14 | 47 | 48 | 353us | 10.2ms | 28.8x | PASS |
| PALMER1ANE | 6 | 35 | LocalInfeasi | IpoptStatus( | N/A | 45 | 0 | 491us | 571us | 1.2x | BOTH_FAIL |
| PALMER1B | 4 | 0 | Optimal | Optimal | 4.12e-15 | 18 | 17 | 138us | 4.0ms | 28.6x | PASS |
| PALMER1BNE | 4 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 212us | 1.1ms | 5.0x | BOTH_FAIL |
| PALMER1C | 8 | 0 | Optimal | Optimal | 1.79e-13 | 4 | 1 | 67us | 1.0ms | 15.0x | PASS |
| PALMER1D | 7 | 0 | Optimal | Optimal | 3.71e-14 | 2 | 1 | 61us | 1.0ms | 16.4x | PASS |
| PALMER1E | 8 | 0 | Optimal | Optimal | 4.05e-13 | 119 | 55 | 989us | 11.7ms | 11.8x | PASS |
| PALMER1ENE | 8 | 35 | LocalInfeasi | IpoptStatus( | N/A | 121 | 0 | 1.5ms | 584us | 0.4x | BOTH_FAIL |
| PALMER1NE | 4 | 31 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 176us | 633us | 3.6x | BOTH_FAIL |
| PALMER2 | 4 | 0 | Optimal | Optimal | 4.98e-16 | 17 | 28 | 118us | 7.2ms | 61.1x | PASS |
| PALMER2A | 6 | 0 | Optimal | Optimal | 7.77e-16 | 71 | 91 | 384us | 20.5ms | 53.2x | PASS |
| PALMER2ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 500us | 607us | 1.2x | BOTH_FAIL |
| PALMER2B | 4 | 0 | Optimal | Optimal | 5.44e-15 | 17 | 15 | 115us | 4.1ms | 35.5x | PASS |
| PALMER2BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 195us | 572us | 2.9x | BOTH_FAIL |
| PALMER2C | 8 | 0 | Optimal | Optimal | 5.44e-15 | 1 | 1 | 38us | 927us | 24.2x | PASS |
| PALMER2E | 8 | 0 | Optimal | Optimal | 2.09e-12 | 113 | 114 | 765us | 25.4ms | 33.2x | PASS |
| PALMER2ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 109 | 0 | 998us | 630us | 0.6x | BOTH_FAIL |
| PALMER2NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 181us | 592us | 3.3x | BOTH_FAIL |
| PALMER3 | 4 | 0 | Acceptable | Optimal | 6.25e-02 | 20 | 44 | 143us | 8.6ms | 60.0x | MISMATCH |
| PALMER3A | 6 | 0 | Optimal | Optimal | 3.89e-16 | 66 | 73 | 338us | 15.8ms | 46.8x | PASS |
| PALMER3ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 74 | 0 | 660us | 712us | 1.1x | BOTH_FAIL |
| PALMER3B | 4 | 0 | Optimal | Optimal | 1.47e-15 | 15 | 15 | 111us | 3.9ms | 34.8x | PASS |
| PALMER3BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 142us | 568us | 4.0x | BOTH_FAIL |
| PALMER3C | 8 | 0 | Optimal | Optimal | 3.09e-15 | 1 | 1 | 38us | 912us | 23.9x | PASS |
| PALMER3E | 8 | 0 | Optimal | Optimal | 1.94e-13 | 28 | 32 | 221us | 6.7ms | 30.3x | PASS |
| PALMER3ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 328us | 533us | 1.6x | BOTH_FAIL |
| PALMER3NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 270us | 577us | 2.1x | BOTH_FAIL |
| PALMER4 | 4 | 0 | Optimal | Optimal | 5.72e-02 | 29 | 16 | 162us | 4.2ms | 26.0x | MISMATCH |
| PALMER4A | 6 | 0 | Optimal | Optimal | 1.79e-15 | 48 | 53 | 266us | 10.7ms | 40.4x | PASS |
| PALMER4ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 57 | 0 | 454us | 560us | 1.2x | BOTH_FAIL |
| PALMER4B | 4 | 0 | Optimal | Optimal | 1.82e-15 | 14 | 16 | 110us | 4.3ms | 39.1x | PASS |
| PALMER4BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 124us | 574us | 4.6x | BOTH_FAIL |
| PALMER4C | 8 | 0 | Optimal | Optimal | 3.29e-15 | 1 | 1 | 39us | 1.0ms | 25.8x | PASS |
| PALMER4E | 8 | 0 | Optimal | Optimal | 8.84e-16 | 25 | 25 | 201us | 5.8ms | 29.0x | PASS |
| PALMER4ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 330us | 610us | 1.9x | BOTH_FAIL |
| PALMER4NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 227us | 578us | 2.5x | BOTH_FAIL |
| PALMER5A | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 13.3ms | 646.2ms | 48.4x | BOTH_FAIL |
| PALMER5ANE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 17.3ms | 586us | 0.0x | BOTH_FAIL |
| PALMER5B | 9 | 0 | Acceptable | Optimal | 1.57e-13 | 57 | 113 | 300us | 23.4ms | 78.0x | PASS |
| PALMER5BNE | 9 | 12 | LocalInfeasi | IpoptStatus( | N/A | 65 | 0 | 451us | 715us | 1.6x | BOTH_FAIL |
| PALMER5C | 6 | 0 | Optimal | Optimal | 2.50e-15 | 1 | 1 | 33us | 983us | 30.0x | PASS |
| PALMER5D | 4 | 0 | Optimal | Optimal | 3.25e-16 | 1 | 1 | 30us | 982us | 32.9x | PASS |
| PALMER5E | 8 | 0 | Acceptable | MaxIteration | N/A | 14 | 3000 | 110us | 485.0ms | 4424.4x | ipopt_FAIL |
| PALMER5ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 132us | 556us | 4.2x | BOTH_FAIL |
| PALMER6A | 6 | 0 | Optimal | Optimal | 2.42e-15 | 115 | 105 | 440us | 20.4ms | 46.3x | PASS |
| PALMER6ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 104 | 0 | 573us | 624us | 1.1x | BOTH_FAIL |
| PALMER6C | 8 | 0 | Optimal | Optimal | 2.30e-15 | 1 | 1 | 37us | 950us | 26.0x | PASS |
| PALMER6E | 8 | 0 | Optimal | Optimal | 2.54e-11 | 37 | 30 | 198us | 6.9ms | 34.7x | PASS |
| PALMER6ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 285us | 618us | 2.2x | BOTH_FAIL |
| PALMER7A | 6 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 10.0ms | 501.1ms | 49.9x | ipopt_FAIL |
| PALMER7ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 14.2ms | 647us | 0.0x | BOTH_FAIL |
| PALMER7C | 8 | 0 | Optimal | Optimal | 2.59e-13 | 2 | 1 | 42us | 991us | 23.5x | PASS |
| PALMER7E | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 13.4ms | 640.0ms | 47.6x | BOTH_FAIL |
| PALMER7ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 18.3ms | 669us | 0.0x | BOTH_FAIL |
| PALMER8A | 6 | 0 | Optimal | Optimal | 1.90e-15 | 28 | 36 | 138us | 8.5ms | 61.5x | PASS |
| PALMER8ANE | 6 | 12 | LocalInfeasi | IpoptStatus( | N/A | 36 | 0 | 237us | 585us | 2.5x | BOTH_FAIL |
| PALMER8C | 8 | 0 | Optimal | Optimal | 8.33e-15 | 1 | 1 | 37us | 992us | 26.9x | PASS |
| PALMER8E | 8 | 0 | Optimal | Optimal | 1.65e-17 | 30 | 23 | 175us | 5.5ms | 31.3x | PASS |
| PALMER8ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 215us | 589us | 2.7x | BOTH_FAIL |
| PARKCH | 15 | 0 | Acceptable | Optimal | 1.27e-14 | 19 | 17 | 4.57s | 3.82s | 0.8x | PASS |
| PENTAGON | 6 | 15 | Optimal | Optimal | 2.18e-11 | 13 | 19 | 195us | 5.2ms | 26.7x | PASS |
| PFIT1 | 3 | 3 | Optimal | Infeasible | N/A | 257 | 266 | 680us | 47.0ms | 69.1x | ipopt_FAIL |
| PFIT1LS | 3 | 0 | Optimal | Optimal | 1.64e-20 | 226 | 263 | 447us | 47.2ms | 105.5x | PASS |
| PFIT2 | 3 | 3 | Optimal | RestorationF | N/A | 111 | 247 | 326us | 49.7ms | 152.2x | ipopt_FAIL |
| PFIT2LS | 3 | 0 | Optimal | Optimal | 1.47e-20 | 106 | 82 | 220us | 15.7ms | 71.4x | PASS |
| PFIT3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 115 | 133 | 341us | 28.3ms | 83.0x | PASS |
| PFIT3LS | 3 | 0 | Optimal | Optimal | 2.90e-20 | 119 | 132 | 245us | 24.8ms | 101.4x | PASS |
| PFIT4 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 209 | 190 | 562us | 37.8ms | 67.3x | PASS |
| PFIT4LS | 3 | 0 | Optimal | Optimal | 6.57e-20 | 206 | 215 | 390us | 38.7ms | 99.3x | PASS |
| POLAK1 | 3 | 2 | Optimal | Optimal | 3.67e-09 | 9 | 5 | 66us | 1.8ms | 27.8x | PASS |
| POLAK2 | 11 | 2 | Optimal | Optimal | 3.59e-11 | 30 | 10 | 194us | 2.8ms | 14.5x | PASS |
| POLAK3 | 12 | 10 | Optimal | MaxIteration | N/A | 15 | 3000 | 382us | 677.6ms | 1774.8x | ipopt_FAIL |
| POLAK4 | 3 | 3 | Acceptable | Optimal | 4.53e-09 | 14 | 4 | 75us | 1.7ms | 22.2x | PASS |
| POLAK5 | 3 | 2 | Acceptable | Optimal | 1.73e-10 | 30 | 31 | 152us | 6.1ms | 40.5x | PASS |
| POLAK6 | 5 | 4 | Optimal | MaxIteration | N/A | 13 | 3000 | 116us | 851.1ms | 7319.0x | ipopt_FAIL |
| PORTFL1 | 12 | 1 | Acceptable | Optimal | 1.28e-06 | 10 | 9 | 323us | 3.1ms | 9.5x | PASS |
| PORTFL2 | 12 | 1 | Acceptable | Optimal | 4.37e-07 | 10 | 8 | 271us | 2.4ms | 8.9x | PASS |
| PORTFL3 | 12 | 1 | Acceptable | Optimal | 1.50e-06 | 10 | 9 | 278us | 2.9ms | 10.3x | PASS |
| PORTFL4 | 12 | 1 | Acceptable | Optimal | 3.83e-06 | 9 | 8 | 249us | 2.8ms | 11.1x | PASS |
| PORTFL6 | 12 | 1 | Acceptable | Optimal | 4.96e-07 | 10 | 8 | 281us | 2.5ms | 9.0x | PASS |
| POWELLBS | 2 | 2 | Optimal | Optimal | 0.00e+00 | 89 | 11 | 225us | 2.3ms | 10.0x | PASS |
| POWELLBSLS | 2 | 0 | Optimal | Optimal | 6.26e-26 | 90 | 91 | 145us | 13.8ms | 95.0x | PASS |
| POWELLSQ | 2 | 2 | Acceptable | Infeasible | N/A | 9 | 29 | 64us | 6.1ms | 95.4x | ipopt_FAIL |
| POWELLSQLS | 2 | 0 | Acceptable | Optimal | 6.89e-11 | 525 | 10 | 622us | 2.4ms | 3.9x | PASS |
| PRICE3NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 60us | 1.7ms | 28.2x | PASS |
| PRICE4 | 2 | 0 | Optimal | Optimal | 1.01e-22 | 13 | 8 | 59us | 2.0ms | 33.0x | PASS |
| PRICE4B | 2 | 0 | Optimal | Optimal | 3.11e-12 | 10 | 8 | 47us | 2.1ms | 44.5x | PASS |
| PRICE4NE | 2 | 2 | Optimal | Acceptable | 0.00e+00 | 10 | 23 | 64us | 4.0ms | 63.0x | PASS |
| PRODPL0 | 60 | 29 | Acceptable | Optimal | 1.33e-06 | 20 | 15 | 2.9ms | 4.4ms | 1.5x | PASS |
| PRODPL1 | 60 | 29 | Acceptable | Optimal | 9.38e-06 | 44 | 28 | 6.5ms | 6.9ms | 1.1x | PASS |
| PSPDOC | 4 | 0 | Optimal | Optimal | 3.50e-09 | 7 | 5 | 43us | 1.9ms | 43.8x | PASS |
| PT | 2 | 501 | Optimal | Optimal | 8.94e-03 | 4 | 106 | 2.0ms | 84.0ms | 42.6x | MISMATCH |
| QC | 9 | 4 | Acceptable | Optimal | 2.42e-02 | 12 | 44 | 159us | 10.1ms | 63.1x | MISMATCH |
| QCNEW | 9 | 3 | MaxIteration | Optimal | N/A | 2999 | 6 | 55.7ms | 2.0ms | 0.0x | ripopt_FAIL |
| QPCBLEND | 83 | 74 | Optimal | Optimal | 5.56e-07 | 37 | 19 | 3.0ms | 6.7ms | 2.2x | PASS |
| QPNBLEND | 83 | 74 | Optimal | Optimal | 4.77e-07 | 50 | 18 | 4.2ms | 6.0ms | 1.4x | PASS |
| RAT42 | 3 | 9 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 142us | 561us | 4.0x | BOTH_FAIL |
| RAT42LS | 3 | 0 | Optimal | Optimal | 8.82e-16 | 22 | 28 | 100us | 4.7ms | 46.6x | PASS |
| RAT43 | 4 | 15 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 158us | 613us | 3.9x | BOTH_FAIL |
| RAT43LS | 4 | 0 | Acceptable | Optimal | 9.65e-01 | 2999 | 34 | 27.0ms | 6.0ms | 0.2x | MISMATCH |
| RECIPE | 3 | 3 | Acceptable | Optimal | 0.00e+00 | 19 | 16 | 90us | 3.2ms | 35.8x | PASS |
| RECIPELS | 3 | 0 | Acceptable | Optimal | 1.17e-10 | 19 | 29 | 63us | 5.6ms | 89.0x | PASS |
| RES | 20 | 14 | Acceptable | Optimal | 0.00e+00 | 12 | 10 | 233us | 2.5ms | 10.6x | PASS |
| RK23 | 17 | 11 | Acceptable | Optimal | 6.20e-05 | 17 | 10 | 260us | 2.9ms | 11.0x | PASS |
| ROBOT | 14 | 2 | Optimal | IpoptStatus( | N/A | 13 | 18 | 225us | 4.9ms | 21.8x | ipopt_FAIL |
| ROSENBR | 2 | 0 | Optimal | Optimal | 0.00e+00 | 21 | 21 | 52us | 3.7ms | 71.0x | PASS |
| ROSENBRTU | 2 | 0 | Optimal | Optimal | 1.52e-24 | 45 | 87 | 100us | 14.3ms | 142.8x | PASS |
| ROSENMMX | 5 | 4 | Optimal | Optimal | 2.27e-10 | 10 | 13 | 111us | 3.6ms | 32.3x | PASS |
| ROSZMAN1 | 4 | 25 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 134us | 596us | 4.4x | BOTH_FAIL |
| ROSZMAN1LS | 4 | 0 | Optimal | Optimal | 7.59e-19 | 52 | 28 | 271us | 5.4ms | 20.0x | PASS |
| RSNBRNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 82us | 1.1ms | 13.5x | PASS |
| S268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 43us | 3.2ms | 74.4x | MISMATCH |
| S308 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 49us | 2.2ms | 45.2x | PASS |
| S308NE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 51us | 522us | 10.2x | BOTH_FAIL |
| S316-322 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 60us | 2.0ms | 33.4x | PASS |
| S365 | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 28.5ms | 1.4ms | 0.0x | BOTH_FAIL |
| S365MOD | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 29.3ms | 1.5ms | 0.1x | BOTH_FAIL |
| SANTA | 21 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 703us | 649us | 0.9x | BOTH_FAIL |
| SANTALS | 21 | 0 | Optimal | Optimal | 2.36e-09 | 35 | 31 | 458us | 8.0ms | 17.4x | PASS |
| SIM2BQP | 2 | 0 | Acceptable | Optimal | 2.25e-06 | 8 | 5 | 48us | 1.8ms | 36.8x | PASS |
| SIMBQP | 2 | 0 | Optimal | Optimal | 8.29e-09 | 6 | 5 | 38us | 1.7ms | 45.4x | PASS |
| SIMPLLPA | 2 | 2 | Optimal | Optimal | 1.18e-08 | 6 | 8 | 46us | 2.3ms | 50.6x | PASS |
| SIMPLLPB | 2 | 3 | Optimal | Optimal | 9.07e-09 | 3 | 7 | 42us | 2.3ms | 54.5x | PASS |
| SINEVAL | 2 | 0 | Optimal | Optimal | 1.17e-40 | 42 | 42 | 92us | 7.4ms | 80.2x | PASS |
| SINVALNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 79us | 982us | 12.4x | PASS |
| SIPOW1 | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 81 | 33.8ms | 222.2ms | 6.6x | MISMATCH |
| SIPOW1M | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 88 | 39.2ms | 234.9ms | 6.0x | MISMATCH |
| SIPOW2 | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 69 | 32.9ms | 177.1ms | 5.4x | MISMATCH |
| SIPOW2M | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 73 | 37.1ms | 184.8ms | 5.0x | MISMATCH |
| SIPOW3 | 4 | 2000 | Optimal | Optimal | 7.10e-01 | 1 | 12 | 36.3ms | 36.9ms | 1.0x | MISMATCH |
| SIPOW4 | 4 | 2000 | Optimal | Optimal | 8.69e-01 | 1 | 11 | 48.9ms | 34.0ms | 0.7x | MISMATCH |
| SISSER | 2 | 0 | Acceptable | Optimal | 8.15e-11 | 15 | 18 | 52us | 2.9ms | 55.4x | PASS |
| SISSER2 | 2 | 0 | Acceptable | Optimal | 7.49e-11 | 16 | 20 | 55us | 3.3ms | 60.6x | PASS |
| SNAIL | 2 | 0 | Optimal | Optimal | 1.63e-26 | 63 | 63 | 124us | 9.9ms | 79.8x | PASS |
| SNAKE | 2 | 2 | Optimal | Optimal | 2.00e-04 | 5 | 8 | 63us | 2.6ms | 41.2x | MISMATCH |
| SPANHYD | 97 | 33 | Optimal | Optimal | 2.37e-16 | 67 | 20 | 7.8ms | 6.5ms | 0.8x | PASS |
| SPIRAL | 3 | 2 | Acceptable | Infeasible | N/A | 116 | 370 | 365us | 63.5ms | 174.0x | ipopt_FAIL |
| SSI | 3 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 3.8ms | 437.9ms | 115.2x | ipopt_FAIL |
| SSINE | 3 | 2 | Acceptable | Optimal | 0.00e+00 | 77 | 224 | 341us | 34.6ms | 101.4x | PASS |
| STANCMIN | 3 | 2 | Optimal | Optimal | 9.32e-09 | 10 | 9 | 79us | 2.6ms | 32.4x | PASS |
| STRATEC | 10 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| STREG | 4 | 0 | Optimal | Optimal | 8.90e-02 | 21 | 13 | 68us | 2.8ms | 41.9x | MISMATCH |
| STREGNE | 4 | 2 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 37us | 1.2ms | 31.7x | PASS |
| SUPERSIM | 2 | 2 | Optimal | Optimal | 2.22e-16 | 7 | 1 | 61us | 1.1ms | 18.5x | PASS |
| SWOPF | 83 | 92 | Optimal | Optimal | 9.19e-10 | 15 | 13 | 1.4ms | 4.3ms | 3.2x | PASS |
| SYNTHES1 | 6 | 6 | Acceptable | Optimal | 1.08e-05 | 13 | 8 | 99us | 2.4ms | 24.1x | PASS |
| SYNTHES2 | 11 | 14 | Optimal | Optimal | 8.24e-07 | 21 | 14 | 221us | 3.4ms | 15.2x | PASS |
| SYNTHES3 | 17 | 23 | Optimal | Optimal | 3.99e-08 | 36 | 13 | 706us | 3.4ms | 4.8x | PASS |
| TAME | 2 | 1 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 43us | 1.6ms | 37.2x | PASS |
| TAX13322 | 72 | 1261 | Acceptable | MaxIteration | N/A | 112 | 3000 | 264.3ms | 19.67s | 74.4x | ipopt_FAIL |
| TAXR13322 | 72 | 1261 | Acceptable | Acceptable | 9.29e-01 | 189 | 56 | 427.0ms | 2.70s | 6.3x | MISMATCH |
| TENBARS1 | 18 | 9 | Acceptable | Optimal | 1.19e-08 | 181 | 39 | 2.7ms | 7.8ms | 2.9x | PASS |
| TENBARS2 | 18 | 8 | Acceptable | Optimal | 1.16e-08 | 125 | 33 | 4.0ms | 6.9ms | 1.7x | PASS |
| TENBARS3 | 18 | 8 | Optimal | Optimal | 1.01e-08 | 25 | 34 | 329us | 7.2ms | 21.9x | PASS |
| TENBARS4 | 18 | 9 | Optimal | Optimal | 1.01e-10 | 14 | 14 | 233us | 3.9ms | 16.7x | PASS |
| TFI1 | 3 | 101 | Optimal | Optimal | 1.00e+00 | 14 | 19 | 1.4ms | 7.2ms | 5.1x | MISMATCH |
| TFI2 | 3 | 101 | Optimal | Optimal | 2.79e-03 | 48 | 8 | 2.5ms | 3.3ms | 1.3x | MISMATCH |
| TFI3 | 3 | 101 | Optimal | Optimal | 6.17e-09 | 79 | 13 | 3.5ms | 4.7ms | 1.4x | PASS |
| THURBER | 7 | 37 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 300us | 583us | 1.9x | BOTH_FAIL |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 1.24e-14 | 30 | 19 | 357us | 3.9ms | 11.0x | PASS |
| TOINTGOR | 50 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 252us | 1.9ms | 7.5x | PASS |
| TOINTPSP | 50 | 0 | Optimal | Optimal | 0.00e+00 | 21 | 20 | 616us | 5.1ms | 8.3x | PASS |
| TOINTQOR | 50 | 0 | Optimal | Optimal | 1.93e-16 | 1 | 1 | 72us | 1.1ms | 14.7x | PASS |
| TRIGGER | 7 | 6 | Acceptable | Optimal | 0.00e+00 | 15 | 15 | 109us | 3.2ms | 29.1x | PASS |
| TRO3X3 | 30 | 13 | Acceptable | Optimal | 3.58e-03 | 2999 | 47 | 64.2ms | 10.5ms | 0.2x | MISMATCH |
| TRO4X4 | 63 | 25 | Acceptable | IpoptStatus( | N/A | 2999 | 157 | 318.6ms | 42.5ms | 0.1x | ipopt_FAIL |
| TRO6X2 | 45 | 21 | Acceptable | RestorationF | N/A | 25 | 353 | 2.0ms | 91.7ms | 45.7x | ipopt_FAIL |
| TRUSPYR1 | 11 | 4 | Optimal | Optimal | 2.03e-08 | 30 | 10 | 288us | 2.7ms | 9.2x | PASS |
| TRUSPYR2 | 11 | 11 | Acceptable | Optimal | 1.72e-08 | 27 | 13 | 402us | 3.4ms | 8.4x | PASS |
| TRY-B | 2 | 1 | Optimal | Optimal | 5.01e-19 | 13 | 23 | 66us | 4.9ms | 74.9x | PASS |
| TWOBARS | 2 | 2 | Optimal | Optimal | 7.13e-01 | 19 | 8 | 112us | 2.2ms | 20.2x | MISMATCH |
| VESUVIA | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 462 | 0 | 170.2ms | 935us | 0.0x | BOTH_FAIL |
| VESUVIALS | 8 | 0 | Acceptable | Optimal | 3.39e-01 | 330 | 48 | 72.1ms | 18.5ms | 0.3x | MISMATCH |
| VESUVIO | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 1.55s | 841us | 0.0x | BOTH_FAIL |
| VESUVIOLS | 8 | 0 | Acceptable | Optimal | 2.41e-15 | 2999 | 10 | 1.04s | 5.2ms | 0.0x | PASS |
| VESUVIOU | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 9.1ms | 821us | 0.1x | BOTH_FAIL |
| VESUVIOULS | 8 | 0 | Acceptable | Optimal | 5.55e-16 | 22 | 8 | 5.6ms | 3.9ms | 0.7x | PASS |
| VIBRBEAM | 8 | 0 | Optimal | Optimal | 9.48e-01 | 115 | 58 | 2.4ms | 9.8ms | 4.1x | MISMATCH |
| VIBRBEAMNE | 8 | 30 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 1.3ms | 576us | 0.4x | BOTH_FAIL |
| WACHBIEG | 3 | 2 | Optimal | Infeasible | N/A | 12 | 15 | 189us | 4.0ms | 21.1x | ipopt_FAIL |
| WATER | 31 | 10 | Acceptable | Optimal | 1.99e-07 | 37 | 17 | 1.0ms | 3.9ms | 3.8x | PASS |
| WAYSEA1 | 2 | 0 | Optimal | Optimal | 2.69e-15 | 15 | 14 | 50us | 2.5ms | 49.2x | PASS |
| WAYSEA1B | 2 | 0 | Optimal | Optimal | 7.97e-09 | 13 | 14 | 46us | 3.3ms | 70.4x | PASS |
| WAYSEA1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 15 | 7 | 67us | 1.7ms | 25.5x | PASS |
| WAYSEA2 | 2 | 0 | Optimal | Optimal | 9.85e-18 | 23 | 22 | 69us | 3.3ms | 47.8x | PASS |
| WAYSEA2B | 2 | 0 | Optimal | Optimal | 1.61e-11 | 21 | 22 | 62us | 4.3ms | 69.5x | PASS |
| WAYSEA2NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 22 | 11 | 79us | 2.2ms | 28.0x | PASS |
| WEEDS | 3 | 0 | Optimal | Optimal | 3.43e-16 | 23 | 28 | 110us | 6.5ms | 58.6x | PASS |
| WEEDSNE | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 192us | 566us | 3.0x | BOTH_FAIL |
| WOMFLET | 3 | 3 | Acceptable | Optimal | 1.00e+00 | 14 | 8 | 88us | 2.2ms | 25.2x | MISMATCH |
| YFIT | 3 | 0 | Optimal | Optimal | 1.37e-19 | 36 | 36 | 146us | 7.5ms | 50.9x | PASS |
| YFITNE | 3 | 17 | Acceptable | IpoptStatus( | N/A | 36 | 0 | 243us | 606us | 2.5x | ipopt_FAIL |
| YFITU | 3 | 0 | Optimal | Optimal | 6.69e-21 | 36 | 36 | 147us | 5.9ms | 40.0x | PASS |
| ZANGWIL2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 25us | 940us | 37.4x | PASS |
| ZANGWIL3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 28us | 953us | 34.4x | PASS |
| ZECEVIC2 | 2 | 2 | Optimal | Optimal | 5.57e-10 | 11 | 8 | 60us | 2.2ms | 36.9x | PASS |
| ZECEVIC3 | 2 | 2 | Optimal | Optimal | 8.24e-10 | 12 | 17 | 91us | 3.6ms | 40.1x | PASS |
| ZECEVIC4 | 2 | 2 | Optimal | Optimal | 2.62e-09 | 11 | 10 | 68us | 2.6ms | 38.3x | PASS |
| ZY2 | 3 | 2 | Acceptable | Optimal | 4.67e-05 | 13 | 14 | 82us | 3.5ms | 42.9x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 136.2 | 43.8 |
| Median | 14 | 13 |
| Total  | 74643 | 23994 |

- ripopt fewer iterations: 168/548
- Ipopt fewer iterations: 292/548
- Tied: 88/548

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 124.5ms | 35.1ms |
| Median | 112us | 3.2ms |
| Total  | 68.23s | 19.26s |

- Geometric mean speedup (Ipopt_time/ripopt_time): **14.86x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 510/548 problems
- Ipopt faster: 38/548 problems
- Overall speedup (total time): 0.28x

## Failure Analysis

### Problems where only ripopt fails (12)

| Problem | n | m | ripopt status | Ipopt obj |
|---------|---|---|---------------|-----------|
| ACOPR14 | 38 | 82 | MaxIterations | 8.081526e+03 |
| ACOPR30 | 72 | 172 | RestorationFailed | 5.768924e+02 |
| CRESC50 | 6 | 100 | RestorationFailed | 7.862467e-01 |
| HATFLDH | 4 | 7 | MaxIterations | -2.450000e+01 |
| HS109 | 9 | 10 | MaxIterations | 5.362069e+03 |
| HS83 | 5 | 3 | MaxIterations | -3.066554e+04 |
| MGH10SLS | 3 | 0 | MaxIterations | 8.794586e+01 |
| OET2 | 3 | 1002 | MaxIterations | 8.715962e-02 |
| OET6 | 5 | 1002 | MaxIterations | 2.069727e-03 |
| OET7 | 7 | 1002 | MaxIterations | 4.465915e-05 |
| OSBORNEA | 5 | 0 | MaxIterations | 5.464895e-05 |
| QCNEW | 9 | 3 | MaxIterations | -8.065219e+02 |

### Problems where only Ipopt fails (38)

| Problem | n | m | Ipopt status | ripopt obj |
|---------|---|---|--------------|------------|
| ARGAUSS | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| AVION2 | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CRESC100 | 6 | 200 | Infeasible | 5.676027e-01 |
| DECONVB | 63 | 0 | MaxIterations | 3.233019e-03 |
| DENSCHNBNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DENSCHNENE | 3 | 3 | Infeasible | 0.000000e+00 |
| DEVGLA1NE | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | 9 | 3 | ErrorInStepComputation | -8.274326e+02 |
| EXP2NE | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| LANCZOS1 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS2 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS3 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 |
| MESH | 41 | 48 | IpoptStatus(4) | -1.798364e+08 |
| NYSTROM5 | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | 6 | 0 | MaxIterations | 1.033491e+01 |
| PFIT1 | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | 3 | 2 | Infeasible | 2.322254e-12 |
| SSI | 3 | 0 | MaxIterations | 1.376194e-09 |
| TAX13322 | 72 | 1261 | MaxIterations | -7.250498e+03 |
| TRO4X4 | 63 | 25 | IpoptStatus(4) | 8.999898e+00 |
| TRO6X2 | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | 3 | 2 | Infeasible | 1.000000e+00 |
| YFITNE | 3 | 17 | IpoptStatus(-10) | 0.000000e+00 |

### Problems where both fail (129)

| Problem | n | m | ripopt status | Ipopt status |
|---------|---|---|---------------|--------------|
| BARDNE | 3 | 15 | LocalInfeasibility | IpoptStatus(-10) |
| BENNETT5 | 3 | 154 | LocalInfeasibility | IpoptStatus(-10) |
| BIGGS6NE | 6 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| BLEACHNG | 17 | 0 | Timeout | Timeout |
| BOXBOD | 2 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| BROWNDENE | 4 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| BURKEHAN | 1 | 1 | RestorationFailed | Infeasible |
| CERI651A | 7 | 61 | LocalInfeasibility | IpoptStatus(-10) |
| CERI651B | 7 | 66 | LocalInfeasibility | IpoptStatus(-10) |
| CERI651C | 7 | 56 | LocalInfeasibility | IpoptStatus(-10) |
| CERI651D | 7 | 67 | LocalInfeasibility | IpoptStatus(-10) |
| CERI651E | 7 | 64 | LocalInfeasibility | IpoptStatus(-10) |
| CHWIRUT1 | 3 | 214 | LocalInfeasibility | IpoptStatus(-10) |
| CHWIRUT2 | 3 | 54 | LocalInfeasibility | IpoptStatus(-10) |
| CRESC132 | 6 | 2654 | Timeout | Timeout |
| DANIWOOD | 2 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| DANWOOD | 2 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| DEVGLA2NE | 5 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| DIAMON2D | 66 | 4643 | Timeout | Timeout |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout |
| DIAMON3D | 99 | 4643 | Timeout | Timeout |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout |
| DMN15102 | 66 | 4643 | Timeout | Timeout |
| DMN15102LS | 66 | 0 | Timeout | Timeout |
| DMN15103 | 99 | 4643 | Timeout | Timeout |
| DMN15103LS | 99 | 0 | Timeout | Timeout |
| DMN15332 | 66 | 4643 | Timeout | Timeout |
| DMN15332LS | 66 | 0 | Timeout | Timeout |
| DMN15333 | 99 | 4643 | Timeout | Timeout |
| DMN15333LS | 99 | 0 | Timeout | Timeout |
| DMN37142 | 66 | 4643 | Timeout | Timeout |
| DMN37142LS | 66 | 0 | Timeout | Timeout |
| DMN37143 | 99 | 4643 | Timeout | Timeout |
| DMN37143LS | 99 | 0 | Timeout | Timeout |
| ECKERLE4 | 3 | 35 | LocalInfeasibility | IpoptStatus(-10) |
| EGGCRATENE | 2 | 4 | LocalInfeasibility | IpoptStatus(-10) |
| ELATVIDUNE | 2 | 3 | LocalInfeasibility | IpoptStatus(-10) |
| ENSO | 9 | 168 | LocalInfeasibility | IpoptStatus(-10) |
| EXPFITNE | 2 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| FBRAIN | 2 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| FBRAIN2 | 4 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| FBRAIN2NE | 4 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| FBRAIN3 | 6 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| FBRAIN3LS | 6 | 0 | MaxIterations | MaxIterations |
| FBRAINNE | 2 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS1 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS2 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS3 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GBRAIN | 2 | 2200 | LocalInfeasibility | IpoptStatus(-10) |
| GROWTH | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| HAHN1 | 7 | 236 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDBNE | 4 | 4 | MaxIterations | Infeasible |
| HATFLDDNE | 3 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDENE | 3 | 21 | LocalInfeasibility | IpoptStatus(-10) |
| HIMMELBD | 2 | 2 | RestorationFailed | Infeasible |
| HIMMELBFNE | 4 | 7 | LocalInfeasibility | IpoptStatus(-10) |
| HS25NE | 3 | 99 | LocalInfeasibility | IpoptStatus(-10) |
| HS2NE | 2 | 2 | RestorationFailed | Infeasible |
| HS87 | 6 | 4 | MaxIterations | MaxIterations |
| JENSMPNE | 2 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| JUDGENE | 2 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| KIRBY2 | 5 | 151 | LocalInfeasibility | IpoptStatus(-10) |
| KOEBHELBNE | 3 | 156 | LocalInfeasibility | IpoptStatus(-10) |
| KOWOSBNE | 4 | 11 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE10 | 10 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE5 | 2 | 4 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE6 | 3 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE7 | 4 | 8 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE8 | 5 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE9 | 8 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| LHAIFAM | 99 | 150 | MaxIterations | InvalidNumberDetected |
| LSC1 | 3 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| LSC2 | 3 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| MEYER3NE | 3 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| MGH09 | 4 | 11 | LocalInfeasibility | IpoptStatus(-10) |
| MGH10 | 3 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| MGH10S | 3 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| MGH17 | 5 | 33 | LocalInfeasibility | IpoptStatus(-10) |
| MGH17S | 5 | 33 | LocalInfeasibility | IpoptStatus(-10) |
| MISRA1A | 2 | 14 | LocalInfeasibility | IpoptStatus(-10) |
| MISRA1B | 2 | 14 | LocalInfeasibility | IpoptStatus(-10) |
| MISRA1C | 2 | 14 | LocalInfeasibility | IpoptStatus(-10) |
| MISRA1D | 2 | 14 | LocalInfeasibility | IpoptStatus(-10) |
| MUONSINE | 1 | 512 | LocalInfeasibility | IpoptStatus(-10) |
| NASH | 72 | 24 | RestorationFailed | Infeasible |
| NELSON | 3 | 128 | LocalInfeasibility | IpoptStatus(-10) |
| OSBORNE1 | 5 | 33 | LocalInfeasibility | IpoptStatus(-10) |
| OSBORNE2 | 11 | 65 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER1ANE | 6 | 35 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER1BNE | 4 | 35 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER1ENE | 8 | 35 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER1NE | 4 | 31 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER2ANE | 6 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER2BNE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER2ENE | 8 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER2NE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER3ANE | 6 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER3BNE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER3ENE | 8 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER3NE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER4ANE | 6 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER4BNE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER4ENE | 8 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER4NE | 4 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER5A | 8 | 0 | MaxIterations | MaxIterations |
| PALMER5ANE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER5BNE | 9 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER5ENE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER6ANE | 6 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER6ENE | 8 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER7ANE | 6 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER7E | 8 | 0 | MaxIterations | MaxIterations |
| PALMER7ENE | 8 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER8ANE | 6 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER8ENE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| RAT42 | 3 | 9 | LocalInfeasibility | IpoptStatus(-10) |
| RAT43 | 4 | 15 | LocalInfeasibility | IpoptStatus(-10) |
| ROSZMAN1 | 4 | 25 | LocalInfeasibility | IpoptStatus(-10) |
| S308NE | 2 | 3 | LocalInfeasibility | IpoptStatus(-10) |
| S365 | 7 | 5 | MaxIterations | RestorationFailed |
| S365MOD | 7 | 5 | MaxIterations | RestorationFailed |
| SANTA | 21 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| STRATEC | 10 | 0 | Timeout | Timeout |
| THURBER | 7 | 37 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIA | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIO | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIOU | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VIBRBEAMNE | 8 | 30 | LocalInfeasibility | IpoptStatus(-10) |
| WEEDSNE | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |

### Objective mismatches (93)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 54
- **Convergence gap** (one Acceptable): 39
- **Better objective found by**: ripopt 17, Ipopt 76

| Problem | ripopt obj | Ipopt obj | Rel Diff | r_status | i_status | Better |
|---------|-----------|-----------|----------|----------|----------|--------|
| MAKELA1 | 1.414214e+00 | -1.414214e+00 | 2.00e+00 | Optimal | Optimal | ipopt |
| SIPOW1 | 4.338956e-01 | -1.000000e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW1M | 4.338960e-01 | -1.000001e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW2M | 3.310254e-01 | -1.000005e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| SIPOW2 | 3.310227e-01 | -1.000000e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| GOFFIN | 1.918601e+01 | -9.500002e-09 | 1.00e+00 | Acceptable | Optimal | ipopt |
| DEVGLA1B | 1.057042e+05 | 8.216237e-18 | 1.00e+00 | Acceptable | Optimal | ipopt |
| DEVGLA2 | 5.932686e+01 | 6.672171e-19 | 1.00e+00 | Optimal | Optimal | ipopt |
| LEVYMONT5 | 1.248706e+01 | 1.239502e-25 | 1.00e+00 | Acceptable | Optimal | ipopt |
| OET5 | 2.909604e+11 | 2.650077e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| WOMFLET | 1.010458e-11 | 6.050000e+00 | 1.00e+00 | Acceptable | Optimal | ripopt |
| OET3 | 1.528803e+07 | 4.505043e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| BENNETT5LS | 1.613684e+05 | 5.563289e-04 | 1.00e+00 | Optimal | Optimal | ipopt |
| ELATTAR | 3.978474e+09 | 7.420618e+01 | 1.00e+00 | Acceptable | Optimal | ipopt |
| HS268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| S268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| DANWOODLS | 1.039178e+02 | 4.317308e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| HALDMADS | 1.223712e-04 | 2.218282e+00 | 1.00e+00 | Optimal | Optimal | ripopt |
| TFI1 | 2.459480e+04 | 5.334687e+00 | 1.00e+00 | Optimal | Optimal | ipopt |
| HS16 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 | Optimal | Optimal | ipopt |
| KSIP | 2.768715e+01 | 5.757979e-01 | 9.79e-01 | Optimal | Optimal | ipopt |
| RAT43LS | 2.525083e+05 | 8.786405e+03 | 9.65e-01 | Acceptable | Optimal | ipopt |
| CHWIRUT1LS | 6.551998e+04 | 2.384477e+03 | 9.64e-01 | Acceptable | Optimal | ipopt |
| CHWIRUT2LS | 1.355314e+04 | 5.130480e+02 | 9.62e-01 | Acceptable | Optimal | ipopt |
| MWRIGHT | 1.288383e+00 | 2.497881e+01 | 9.48e-01 | Optimal | Optimal | ripopt |
| VIBRBEAM | 6.346254e+00 | 3.322376e-01 | 9.48e-01 | Optimal | Optimal | ipopt |
| TAXR13322 | -4.807354e+03 | -3.429089e+02 | 9.29e-01 | Acceptable | Acceptable | ripopt |
| BT4 | -4.551055e+01 | -3.704768e+00 | 9.19e-01 | Optimal | Optimal | ripopt |
| OET4 | 9.191546e-01 | 4.295421e-03 | 9.15e-01 | Optimal | Optimal | ipopt |
| BIGGSC4 | -3.128034e+00 | -2.450000e+01 | 8.72e-01 | Acceptable | Optimal | ipopt |
| SIPOW4 | 2.080345e+00 | 2.723620e-01 | 8.69e-01 | Optimal | Optimal | ipopt |
| HIMMELP2 | -6.205394e+01 | -8.198044e+00 | 8.68e-01 | Optimal | Optimal | ripopt |
| HIMMELP3 | -7.913699e+00 | -5.901318e+01 | 8.66e-01 | Optimal | Optimal | ipopt |
| CAMEL6 | -2.154638e-01 | -1.031628e+00 | 7.91e-01 | Optimal | Optimal | ipopt |
| HS23 | 9.472140e+00 | 2.000000e+00 | 7.89e-01 | Optimal | Optimal | ipopt |
| HS54 | -1.566691e-01 | -9.080748e-01 | 7.51e-01 | Optimal | Optimal | ipopt |
| HIMMELP6 | -1.475339e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| HIMMELP5 | -1.475901e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| TWOBARS | 5.257563e+00 | 1.508652e+00 | 7.13e-01 | Optimal | Optimal | ipopt |
| SIPOW3 | 1.846738e+00 | 5.346586e-01 | 7.10e-01 | Optimal | Optimal | ipopt |
| MAKELA4 | 6.877365e-01 | -9.600000e-09 | 6.88e-01 | Optimal | Optimal | ipopt |
| MISTAKE | -5.000000e-01 | -1.000000e+00 | 5.00e-01 | Acceptable | Optimal | ipopt |
| EGGCRATE | 1.897639e+01 | 9.488197e+00 | 5.00e-01 | Optimal | Optimal | ipopt |
| ECKERLE4LS | 4.988568e-01 | 1.463589e-03 | 4.97e-01 | Acceptable | Optimal | ipopt |
| VESUVIALS | 1.500440e+03 | 9.914100e+02 | 3.39e-01 | Acceptable | Optimal | ipopt |
| HYDC20LS | 2.976453e-01 | 2.967522e-15 | 2.98e-01 | Acceptable | Optimal | ipopt |
| OET1 | 8.183445e-01 | 5.382431e-01 | 2.80e-01 | Optimal | Optimal | ipopt |
| AVGASA | -3.383299e+00 | -4.631926e+00 | 2.70e-01 | Acceptable | Optimal | ipopt |
| MADSEN | 8.836579e-01 | 6.164324e-01 | 2.67e-01 | Optimal | Optimal | ipopt |
| EG1 | -1.132801e+00 | -1.429307e+00 | 2.07e-01 | Optimal | Optimal | ipopt |
| HS70 | 1.877383e-01 | 7.498464e-03 | 1.80e-01 | Optimal | Optimal | ipopt |
| AVGASB | -3.717055e+00 | -4.483219e+00 | 1.71e-01 | Acceptable | Optimal | ipopt |
| BT7 | 3.603798e+02 | 3.065000e+02 | 1.50e-01 | Optimal | Optimal | ipopt |
| MUONSINELS | 5.080769e+04 | 4.387412e+04 | 1.36e-01 | Acceptable | Optimal | ipopt |
| MAKELA2 | 8.244898e+00 | 7.200000e+00 | 1.27e-01 | Optimal | Optimal | ipopt |
| MSS1 | -1.600000e+01 | -1.400000e+01 | 1.25e-01 | Acceptable | Optimal | ripopt |
| STREG | 3.743976e-21 | 8.901950e-02 | 8.90e-02 | Optimal | Optimal | ripopt |
| HS67 | -1.060612e+03 | -1.162119e+03 | 8.73e-02 | Optimal | Optimal | ipopt |
| HAHN1LS | 3.086398e+01 | 3.338424e+01 | 7.55e-02 | Acceptable | Optimal | ripopt |
| PALMER3 | 2.265958e+03 | 2.416980e+03 | 6.25e-02 | Acceptable | Optimal | ripopt |
| PALMER4 | 2.424016e+03 | 2.285383e+03 | 5.72e-02 | Optimal | Optimal | ipopt |
| HET-Z | 1.027774e+00 | 1.000000e+00 | 2.70e-02 | Optimal | Optimal | ipopt |
| MGH17SLS | 5.464895e-05 | 2.451788e-02 | 2.45e-02 | Optimal | Optimal | ripopt |
| QC | -9.333976e+02 | -9.565379e+02 | 2.42e-02 | Acceptable | Optimal | ipopt |
| CB2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CHACONN1 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| HS55 | 6.666667e+00 | 6.805833e+00 | 2.04e-02 | Optimal | Optimal | ripopt |
| LRCOVTYPE | 5.901541e-01 | 5.723072e-01 | 1.78e-02 | Optimal | Optimal | ipopt |
| DUALC8 | 1.854450e+04 | 1.830936e+04 | 1.27e-02 | Acceptable | Optimal | ipopt |
| MAXLIKA | 1.149346e+03 | 1.136307e+03 | 1.13e-02 | Acceptable | Optimal | ipopt |
| PT | 1.873320e-01 | 1.783942e-01 | 8.94e-03 | Optimal | Optimal | ipopt |
| CLIFF | 1.997866e-01 | 2.072380e-01 | 7.45e-03 | Optimal | Optimal | ripopt |
| DGOSPEC | -9.887540e+02 | -9.933506e+02 | 4.63e-03 | Acceptable | Optimal | ipopt |
| EXPFITA | 4.737412e-03 | 1.136646e-03 | 3.60e-03 | Optimal | Optimal | ipopt |
| TRO3X3 | 8.999704e+00 | 8.967478e+00 | 3.58e-03 | Acceptable | Optimal | ipopt |
| TFI2 | 6.518174e-01 | 6.490311e-01 | 2.79e-03 | Optimal | Optimal | ipopt |
| DEGENLPB | -3.069162e+01 | -3.076401e+01 | 2.35e-03 | Acceptable | Optimal | ipopt |
| LIN | -1.960628e-02 | -1.757754e-02 | 2.03e-03 | Optimal | Optimal | ripopt |
| DEGENLPA | 3.060434e+00 | 3.054881e+00 | 1.81e-03 | Acceptable | Optimal | ipopt |
| DECONVC | 3.927874e-03 | 2.569475e-03 | 1.36e-03 | Acceptable | Optimal | ipopt |
| HS85 | -2.212898e+00 | -2.215605e+00 | 1.22e-03 | Acceptable | Optimal | ipopt |
| HS59 | -6.743243e+00 | -6.749505e+00 | 9.28e-04 | Acceptable | Optimal | ipopt |
| HS13 | 9.938594e-01 | 9.945785e-01 | 7.19e-04 | Acceptable | Optimal | ripopt |
| HS95 | 1.630984e-02 | 1.561772e-02 | 6.92e-04 | Acceptable | Optimal | ipopt |
| HS96 | 1.630574e-02 | 1.561775e-02 | 6.88e-04 | Acceptable | Optimal | ipopt |
| LAUNCH | 9.008011e+00 | 9.004902e+00 | 3.45e-04 | Acceptable | Optimal | ipopt |
| HS45 | 1.000314e+00 | 1.000000e+00 | 3.14e-04 | Acceptable | Optimal | ipopt |
| LSC2LS | 1.333749e+01 | 1.333439e+01 | 2.32e-04 | Acceptable | Optimal | ipopt |
| DENSCHND | 6.940852e-09 | 2.221899e-04 | 2.22e-04 | Acceptable | Optimal | ripopt |
| HS116 | 9.756763e+01 | 9.758747e+01 | 2.03e-04 | Acceptable | Optimal | ripopt |
| HS17 | 1.000201e+00 | 1.000000e+00 | 2.01e-04 | Optimal | Optimal | ipopt |
| SNAKE | -2.601155e-10 | -1.999999e-04 | 2.00e-04 | Optimal | Optimal | ipopt |
| HS98 | 3.136197e+00 | 3.135806e+00 | 1.25e-04 | Acceptable | Optimal | ipopt |

---
*Generated by cutest_suite/compare.py*