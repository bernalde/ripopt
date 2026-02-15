# CUTEst Benchmark Report

Comparison of ripopt vs Ipopt (C++) on the CUTEst test set.

## Executive Summary

- **Total problems**: 727
- **ripopt solved**: 583/727 (80.2%)
- **Ipopt solved**: 558/727 (76.8%)
- **Both solved**: 544/727
- **Matching solutions** (rel obj diff < 1e-4): 452/544

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (452 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 1.23e-06 |
| Median | 1.58e-13 |
| Max    | 8.50e-05 |

**All both-solved** (544 problems, including 92 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 7.51e-02 |
| Median | 3.08e-10 |
| Max    | 2.00e+00 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 493 | 365 | 341 | 330 | 262 |
| unconstrained | 234 | 218 | 217 | 214 | 190 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| 3PK | 30 | 0 | Optimal | Optimal | 1.42e-15 | 9 | 9 | 205us | 2.5ms | 12.2x | PASS |
| ACOPP14 | 38 | 68 | Optimal | Optimal | 1.05e-09 | 152 | 9 | 16.0ms | 4.1ms | 0.3x | PASS |
| ACOPP30 | 72 | 142 | Optimal | Optimal | 2.40e-09 | 59 | 13 | 13.7ms | 6.5ms | 0.5x | PASS |
| ACOPR14 | 38 | 82 | Acceptable | Optimal | 1.20e-03 | 520 | 13 | 76.6ms | 5.3ms | 0.1x | MISMATCH |
| ACOPR30 | 72 | 172 | Optimal | Optimal | 1.40e-02 | 2443 | 221 | 14.36s | 110.6ms | 0.0x | MISMATCH |
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 57us | 1.3ms | 23.3x | PASS |
| AIRCRFTB | 8 | 0 | Optimal | Optimal | 4.50e-21 | 11 | 15 | 79us | 3.1ms | 39.5x | PASS |
| AIRPORT | 84 | 42 | Optimal | Optimal | 5.34e-09 | 17 | 13 | 4.9ms | 6.0ms | 1.2x | PASS |
| AKIVA | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 82us | 1.6ms | 19.3x | PASS |
| ALLINIT | 4 | 0 | Acceptable | Optimal | 3.69e-09 | 10 | 20 | 63us | 4.0ms | 64.4x | PASS |
| ALLINITA | 4 | 4 | Acceptable | Optimal | 1.42e-05 | 28 | 12 | 174us | 2.9ms | 16.7x | PASS |
| ALLINITC | 4 | 1 | Acceptable | Optimal | 3.25e-05 | 27 | 17 | 140us | 3.4ms | 24.4x | PASS |
| ALLINITU | 4 | 0 | Optimal | Optimal | 3.09e-16 | 8 | 14 | 56us | 2.7ms | 49.4x | PASS |
| ALSOTAME | 2 | 1 | Optimal | Optimal | 1.47e-08 | 10 | 8 | 59us | 2.2ms | 37.0x | PASS |
| ANTWERP | 27 | 10 | Acceptable | Optimal | 4.01e-03 | 148 | 108 | 5.2ms | 23.0ms | 4.4x | MISMATCH |
| ARGAUSS | 3 | 15 | Acceptable | IpoptStatus( | N/A | 2 | 0 | 52us | 564us | 10.8x | ipopt_FAIL |
| AVGASA | 8 | 10 | Optimal | Optimal | 1.25e-08 | 29 | 9 | 202us | 2.5ms | 12.2x | PASS |
| AVGASB | 8 | 10 | Optimal | Optimal | 3.43e-07 | 30 | 11 | 220us | 2.9ms | 13.2x | PASS |
| AVION2 | 49 | 15 | Acceptable | MaxIteration | N/A | 21 | 3000 | 2.7ms | 646.9ms | 240.8x | ipopt_FAIL |
| BA-L1 | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 6 | 632us | 2.3ms | 3.6x | PASS |
| BA-L1LS | 57 | 0 | Optimal | Optimal | 7.65e-21 | 7 | 10 | 679us | 2.9ms | 4.2x | PASS |
| BA-L1SP | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.4ms | 3.1ms | 2.2x | PASS |
| BA-L1SPLS | 57 | 0 | Optimal | Optimal | 6.48e-17 | 23 | 9 | 7.1ms | 4.8ms | 0.7x | PASS |
| BARD | 3 | 0 | Optimal | Optimal | 1.73e-18 | 8 | 8 | 56us | 2.0ms | 34.8x | PASS |
| BARDNE | 3 | 15 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 85us | 586us | 6.9x | BOTH_FAIL |
| BATCH | 48 | 73 | Acceptable | Optimal | 3.60e-08 | 161 | 29 | 7.4ms | 7.8ms | 1.0x | PASS |
| BEALE | 2 | 0 | Optimal | Optimal | 4.34e-18 | 7 | 8 | 44us | 2.2ms | 50.8x | PASS |
| BEALENE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 53us | 578us | 10.9x | ipopt_FAIL |
| BENNETT5 | 3 | 154 | LocalInfeasi | IpoptStatus( | N/A | 581 | 0 | 23.1ms | 634us | 0.0x | BOTH_FAIL |
| BENNETT5LS | 3 | 0 | Optimal | Optimal | 3.23e-05 | 425 | 21 | 9.2ms | 4.4ms | 0.5x | PASS |
| BIGGS3 | 6 | 0 | Optimal | Optimal | 1.30e-18 | 8 | 9 | 75us | 2.5ms | 32.9x | PASS |
| BIGGS5 | 6 | 0 | Optimal | Optimal | 4.74e-20 | 25 | 20 | 159us | 4.1ms | 25.6x | PASS |
| BIGGS6 | 6 | 0 | Optimal | Optimal | 4.32e-20 | 89 | 79 | 404us | 12.2ms | 30.2x | PASS |
| BIGGS6NE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 591 | 0 | 3.8ms | 571us | 0.1x | BOTH_FAIL |
| BIGGSC4 | 4 | 7 | Acceptable | Optimal | 8.72e-01 | 2999 | 17 | 10.1ms | 3.9ms | 0.4x | MISMATCH |
| BLEACHNG | 17 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| BOOTH | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 30us | 976us | 32.6x | PASS |
| BOX2 | 3 | 0 | Optimal | Optimal | 7.83e-18 | 9 | 8 | 68us | 2.0ms | 29.4x | PASS |
| BOX3 | 3 | 0 | Optimal | Optimal | 1.69e-25 | 11 | 9 | 63us | 2.2ms | 34.2x | PASS |
| BOX3NE | 3 | 10 | Optimal | IpoptStatus( | N/A | 11 | 0 | 97us | 601us | 6.2x | ipopt_FAIL |
| BOXBOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 3 | 0 | 50us | 599us | 11.9x | BOTH_FAIL |
| BOXBODLS | 2 | 0 | Optimal | Optimal | 4.40e-13 | 3 | 13 | 46us | 3.0ms | 64.3x | PASS |
| BQP1VAR | 1 | 0 | Optimal | Optimal | 9.98e-09 | 7 | 5 | 36us | 1.7ms | 47.0x | PASS |
| BQPGABIM | 50 | 0 | Acceptable | Optimal | 4.40e-06 | 11 | 12 | 400us | 3.3ms | 8.3x | PASS |
| BQPGASIM | 50 | 0 | Acceptable | Optimal | 6.13e-06 | 11 | 12 | 395us | 3.1ms | 7.9x | PASS |
| BRANIN | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 55us | 2.2ms | 39.7x | PASS |
| BRKMCC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 36us | 1.3ms | 35.7x | PASS |
| BROWNBS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 49us | 1.7ms | 35.3x | PASS |
| BROWNBSNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 8 | 0 | 64us | 594us | 9.3x | ipopt_FAIL |
| BROWNDEN | 4 | 0 | Optimal | Optimal | 1.70e-16 | 8 | 8 | 66us | 1.8ms | 26.8x | PASS |
| BROWNDENE | 4 | 20 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 103us | 583us | 5.7x | BOTH_FAIL |
| BT1 | 2 | 1 | Optimal | Optimal | 2.47e-12 | 12 | 7 | 71us | 1.9ms | 27.2x | PASS |
| BT10 | 2 | 2 | Optimal | Optimal | 2.79e-09 | 7 | 6 | 44us | 1.7ms | 38.5x | PASS |
| BT11 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 8 | 77us | 1.9ms | 25.3x | PASS |
| BT12 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 51us | 1.4ms | 28.2x | PASS |
| BT13 | 5 | 1 | Acceptable | Optimal | 1.00e-08 | 27 | 24 | 108us | 4.7ms | 43.2x | PASS |
| BT2 | 3 | 1 | Optimal | Optimal | 3.80e-12 | 11 | 12 | 58us | 2.3ms | 39.5x | PASS |
| BT3 | 5 | 3 | Optimal | Optimal | 1.22e-14 | 1 | 1 | 35us | 1.0ms | 29.2x | PASS |
| BT4 | 3 | 2 | Optimal | Optimal | 9.19e-01 | 7 | 9 | 63us | 2.3ms | 36.0x | MISMATCH |
| BT5 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 59us | 1.8ms | 29.9x | PASS |
| BT6 | 5 | 2 | Optimal | Optimal | 6.16e-12 | 9 | 13 | 71us | 2.6ms | 37.0x | PASS |
| BT7 | 5 | 3 | Optimal | Optimal | 1.15e-14 | 25 | 16 | 198us | 3.2ms | 16.4x | PASS |
| BT8 | 5 | 2 | Acceptable | Optimal | 3.73e-09 | 32 | 14 | 143us | 2.6ms | 18.0x | PASS |
| BT9 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 77us | 2.4ms | 31.7x | PASS |
| BURKEHAN | 1 | 1 | RestorationF | Infeasible | N/A | 232 | 11 | 15.9ms | 3.2ms | 0.2x | BOTH_FAIL |
| BYRDSPHR | 3 | 2 | Optimal | Optimal | 4.31e-13 | 37 | 12 | 181us | 2.7ms | 14.9x | PASS |
| CAMEL6 | 2 | 0 | Optimal | Optimal | 7.91e-01 | 9 | 8 | 51us | 2.3ms | 45.0x | MISMATCH |
| CANTILVR | 5 | 1 | Optimal | Optimal | 3.30e-09 | 14 | 11 | 70us | 2.8ms | 39.5x | PASS |
| CB2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 8 | 8 | 59us | 2.3ms | 38.9x | MISMATCH |
| CB3 | 3 | 3 | Optimal | Optimal | 4.30e-09 | 8 | 8 | 60us | 2.2ms | 37.3x | PASS |
| CERI651A | 7 | 61 | LocalInfeasi | IpoptStatus( | N/A | 264 | 0 | 6.8ms | 597us | 0.1x | BOTH_FAIL |
| CERI651ALS | 7 | 0 | Acceptable | Optimal | 8.56e-08 | 283 | 95 | 4.9ms | 16.1ms | 3.3x | PASS |
| CERI651B | 7 | 66 | LocalInfeasi | IpoptStatus( | N/A | 83 | 0 | 2.3ms | 611us | 0.3x | BOTH_FAIL |
| CERI651BLS | 7 | 0 | Optimal | Optimal | 1.23e-08 | 89 | 56 | 1.5ms | 9.2ms | 6.0x | PASS |
| CERI651C | 7 | 56 | LocalInfeasi | IpoptStatus( | N/A | 171 | 0 | 3.9ms | 603us | 0.2x | BOTH_FAIL |
| CERI651CLS | 7 | 0 | Optimal | Optimal | 3.90e-09 | 178 | 53 | 2.6ms | 8.3ms | 3.2x | PASS |
| CERI651D | 7 | 67 | LocalInfeasi | IpoptStatus( | N/A | 91 | 0 | 2.7ms | 656us | 0.2x | BOTH_FAIL |
| CERI651DLS | 7 | 0 | Acceptable | Optimal | 3.08e-10 | 98 | 60 | 1.8ms | 10.5ms | 5.9x | PASS |
| CERI651E | 7 | 64 | LocalInfeasi | IpoptStatus( | N/A | 52 | 0 | 1.5ms | 617us | 0.4x | BOTH_FAIL |
| CERI651ELS | 7 | 0 | Acceptable | Optimal | 1.03e-09 | 83 | 45 | 1.4ms | 7.2ms | 5.1x | PASS |
| CHACONN1 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 5 | 6 | 53us | 2.0ms | 38.3x | MISMATCH |
| CHACONN2 | 3 | 3 | Optimal | Optimal | 4.49e-09 | 6 | 6 | 54us | 1.9ms | 35.4x | PASS |
| CHWIRUT1 | 3 | 214 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 674us | 636us | 0.9x | BOTH_FAIL |
| CHWIRUT1LS | 3 | 0 | Optimal | Optimal | 1.91e-16 | 12 | 6 | 399us | 2.1ms | 5.2x | PASS |
| CHWIRUT2 | 3 | 54 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 185us | 596us | 3.2x | BOTH_FAIL |
| CHWIRUT2LS | 3 | 0 | Optimal | Optimal | 2.22e-16 | 12 | 6 | 142us | 2.0ms | 13.7x | PASS |
| CLIFF | 2 | 0 | Optimal | Optimal | 7.45e-03 | 27 | 23 | 74us | 3.5ms | 46.9x | MISMATCH |
| CLUSTER | 2 | 2 | Acceptable | Optimal | 0.00e+00 | 12 | 9 | 74us | 2.2ms | 29.6x | PASS |
| CLUSTERLS | 2 | 0 | Optimal | Optimal | 2.72e-18 | 13 | 17 | 60us | 3.0ms | 49.3x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.03e-08 | 31 | 7 | 293us | 2.2ms | 7.6x | PASS |
| CONGIGMZ | 3 | 5 | Optimal | Optimal | 3.19e-09 | 7 | 20 | 122us | 4.1ms | 33.8x | PASS |
| COOLHANS | 9 | 9 | Optimal | Optimal | 0.00e+00 | 22 | 9 | 215us | 2.1ms | 9.9x | PASS |
| COOLHANSLS | 9 | 0 | Optimal | Optimal | 1.21e-18 | 23 | 25 | 179us | 4.6ms | 25.7x | PASS |
| CORE1 | 65 | 59 | Optimal | Optimal | 4.43e-09 | 25 | 33 | 49.7ms | 8.1ms | 0.2x | PASS |
| CRESC100 | 6 | 200 | Optimal | Infeasible | N/A | 158 | 155 | 805.9ms | 116.0ms | 0.1x | ipopt_FAIL |
| CRESC132 | 6 | 2654 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC4 | 6 | 8 | Optimal | Optimal | 2.06e-08 | 23 | 64 | 383us | 13.4ms | 35.0x | PASS |
| CRESC50 | 6 | 100 | Optimal | Optimal | 2.48e-08 | 175 | 194 | 184.8ms | 81.1ms | 0.4x | PASS |
| CSFI1 | 5 | 4 | Optimal | Optimal | 1.46e-08 | 26 | 11 | 163us | 3.0ms | 18.5x | PASS |
| CSFI2 | 5 | 4 | Optimal | Optimal | 1.49e-08 | 17 | 14 | 217us | 3.4ms | 15.9x | PASS |
| CUBE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 27 | 27 | 67us | 4.7ms | 69.9x | PASS |
| CUBENE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 27 | 1 | 90us | 1.0ms | 11.3x | PASS |
| DALLASS | 46 | 31 | Optimal | Optimal | 4.38e-08 | 29 | 22 | 3.6ms | 5.3ms | 1.5x | PASS |
| DANIWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 77us | 604us | 7.9x | BOTH_FAIL |
| DANIWOODLS | 2 | 0 | Optimal | Optimal | 2.60e-18 | 10 | 10 | 54us | 2.3ms | 42.7x | PASS |
| DANWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 1 | 0 | 39us | 585us | 14.9x | BOTH_FAIL |
| DANWOODLS | 2 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 11 | 31us | 2.5ms | 79.0x | MISMATCH |
| DECONVB | 63 | 0 | Optimal | MaxIteration | N/A | 74 | 3000 | 9.4ms | 707.0ms | 74.9x | ipopt_FAIL |
| DECONVBNE | 63 | 40 | Acceptable | Optimal | 0.00e+00 | 177 | 505 | 24.9ms | 161.4ms | 6.5x | PASS |
| DECONVC | 63 | 1 | Acceptable | Optimal | 1.23e-03 | 50 | 31 | 7.5ms | 9.2ms | 1.2x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 2 | 26 | 513us | 25.0ms | 48.7x | PASS |
| DECONVU | 63 | 0 | Acceptable | Optimal | 7.56e-09 | 35 | 333 | 5.7ms | 85.3ms | 14.9x | PASS |
| DEGENLPA | 20 | 15 | Acceptable | Optimal | 1.41e-01 | 26 | 18 | 11.7ms | 4.1ms | 0.4x | MISMATCH |
| DEGENLPB | 20 | 15 | Acceptable | Optimal | 1.08e-03 | 37 | 19 | 643us | 4.2ms | 6.5x | MISMATCH |
| DEMBO7 | 16 | 20 | Acceptable | Optimal | 5.43e-02 | 26 | 45 | 566us | 9.1ms | 16.0x | MISMATCH |
| DEMYMALO | 3 | 3 | Optimal | Optimal | 2.96e-09 | 12 | 9 | 69us | 2.5ms | 36.7x | PASS |
| DENSCHNA | 2 | 0 | Optimal | Optimal | 5.88e-39 | 6 | 6 | 42us | 1.6ms | 38.2x | PASS |
| DENSCHNB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 47us | 1.9ms | 40.4x | PASS |
| DENSCHNBNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 54us | 585us | 10.8x | ipopt_FAIL |
| DENSCHNC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 49us | 2.1ms | 42.7x | PASS |
| DENSCHNCNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 10 | 7 | 66us | 1.9ms | 28.2x | PASS |
| DENSCHND | 3 | 0 | Acceptable | Optimal | 2.22e-04 | 30 | 26 | 86us | 4.3ms | 50.0x | MISMATCH |
| DENSCHNDNE | 3 | 3 | LocalInfeasi | Acceptable | N/A | 29 | 22 | 118us | 3.9ms | 32.8x | ripopt_FAIL |
| DENSCHNE | 3 | 0 | Optimal | Optimal | 1.86e-17 | 10 | 14 | 48us | 3.2ms | 66.0x | PASS |
| DENSCHNENE | 3 | 3 | Optimal | Infeasible | N/A | 10 | 10 | 63us | 2.8ms | 44.8x | ipopt_FAIL |
| DENSCHNF | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 38us | 1.6ms | 42.0x | PASS |
| DENSCHNFNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 54us | 1.5ms | 28.5x | PASS |
| DEVGLA1 | 4 | 0 | Optimal | Optimal | 3.42e-23 | 21 | 23 | 199us | 4.4ms | 22.0x | PASS |
| DEVGLA1B | 4 | 0 | Acceptable | Optimal | 1.00e+00 | 28 | 20 | 253us | 5.1ms | 20.1x | MISMATCH |
| DEVGLA1NE | 4 | 24 | Optimal | IpoptStatus( | N/A | 16 | 0 | 251us | 590us | 2.3x | ipopt_FAIL |
| DEVGLA2 | 5 | 0 | Optimal | Optimal | 1.00e+00 | 16 | 13 | 190us | 2.8ms | 14.7x | MISMATCH |
| DEVGLA2B | 5 | 0 | Acceptable | Optimal | 2.58e-07 | 14 | 24 | 160us | 5.4ms | 33.9x | PASS |
| DEVGLA2NE | 5 | 16 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 274us | 598us | 2.2x | BOTH_FAIL |
| DGOSPEC | 3 | 0 | Acceptable | Optimal | 4.63e-03 | 10 | 27 | 56us | 5.6ms | 99.3x | MISMATCH |
| DIAMON2D | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3D | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIPIGRI | 7 | 4 | Optimal | Optimal | 1.77e-11 | 14 | 9 | 129us | 2.7ms | 20.9x | PASS |
| DISC2 | 29 | 23 | Optimal | Optimal | 2.87e-09 | 21 | 24 | 1.6ms | 6.4ms | 4.1x | PASS |
| DISCS | 36 | 66 | Acceptable | Optimal | 1.74e-05 | 150 | 184 | 17.7ms | 72.7ms | 4.1x | PASS |
| DIXCHLNG | 10 | 5 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 148us | 2.3ms | 15.3x | PASS |
| DJTL | 2 | 0 | Acceptable | Acceptable | 1.63e-15 | 2999 | 1538 | 8.2ms | 159.9ms | 19.6x | PASS |
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
| DNIEPER | 61 | 24 | Acceptable | Optimal | 4.28e-08 | 200 | 23 | 23.4ms | 5.4ms | 0.2x | PASS |
| DUAL1 | 85 | 1 | Acceptable | Optimal | 3.88e-06 | 12 | 15 | 3.0ms | 6.7ms | 2.2x | PASS |
| DUAL2 | 96 | 1 | Acceptable | Optimal | 8.07e-09 | 12 | 12 | 4.0ms | 6.4ms | 1.6x | PASS |
| DUAL4 | 75 | 1 | Acceptable | Optimal | 1.40e-07 | 11 | 12 | 2.1ms | 5.0ms | 2.3x | PASS |
| DUALC1 | 9 | 215 | Acceptable | Optimal | 1.06e-05 | 59 | 18 | 13.1ms | 11.6ms | 0.9x | PASS |
| DUALC2 | 7 | 229 | Acceptable | Optimal | 1.07e-05 | 25 | 12 | 5.6ms | 8.2ms | 1.5x | PASS |
| DUALC5 | 8 | 278 | Acceptable | Optimal | 1.83e-07 | 12 | 11 | 4.1ms | 8.8ms | 2.1x | PASS |
| DUALC8 | 8 | 503 | Acceptable | Optimal | 6.45e-08 | 50 | 13 | 29.3ms | 15.2ms | 0.5x | PASS |
| ECKERLE4 | 3 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 218us | 576us | 2.6x | BOTH_FAIL |
| ECKERLE4LS | 3 | 0 | Acceptable | Optimal | 4.97e-01 | 25 | 36 | 162us | 6.4ms | 39.7x | MISMATCH |
| EG1 | 3 | 0 | Optimal | Optimal | 2.07e-01 | 9 | 8 | 48us | 2.3ms | 47.8x | MISMATCH |
| EGGCRATE | 2 | 0 | Optimal | Optimal | 5.00e-01 | 7 | 5 | 49us | 1.6ms | 32.2x | MISMATCH |
| EGGCRATEB | 2 | 0 | Optimal | Optimal | 5.62e-16 | 10 | 6 | 48us | 2.0ms | 41.1x | PASS |
| EGGCRATENE | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 62us | 549us | 8.8x | BOTH_FAIL |
| ELATTAR | 7 | 102 | Optimal | Optimal | 9.86e-01 | 98 | 81 | 12.4ms | 35.0ms | 2.8x | MISMATCH |
| ELATVIDU | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 11 | 50us | 2.1ms | 42.5x | PASS |
| ELATVIDUB | 2 | 0 | Optimal | Optimal | 1.03e-11 | 10 | 11 | 40us | 2.6ms | 63.2x | PASS |
| ELATVIDUNE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 62us | 558us | 9.0x | BOTH_FAIL |
| ENGVAL2 | 3 | 0 | Optimal | Optimal | 1.70e-20 | 20 | 21 | 75us | 3.7ms | 49.3x | PASS |
| ENGVAL2NE | 3 | 5 | Optimal | IpoptStatus( | N/A | 17 | 0 | 96us | 527us | 5.5x | ipopt_FAIL |
| ENSO | 9 | 168 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 1.2ms | 681us | 0.5x | BOTH_FAIL |
| ENSOLS | 9 | 0 | Acceptable | Optimal | 5.77e-16 | 13 | 7 | 1.6ms | 2.8ms | 1.7x | PASS |
| EQC | 9 | 3 | Acceptable | ErrorInStepC | N/A | 7 | 15 | 154us | 4.8ms | 31.3x | ipopt_FAIL |
| ERRINBAR | 18 | 9 | Optimal | Optimal | 1.55e-08 | 35 | 37 | 469us | 7.7ms | 16.3x | PASS |
| EXP2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 40us | 1.7ms | 42.8x | PASS |
| EXP2B | 2 | 0 | Optimal | Optimal | 2.25e-15 | 9 | 7 | 45us | 1.8ms | 40.8x | PASS |
| EXP2NE | 2 | 10 | Optimal | IpoptStatus( | N/A | 7 | 0 | 70us | 560us | 8.0x | ipopt_FAIL |
| EXPFIT | 2 | 0 | Optimal | Optimal | 8.33e-17 | 5 | 8 | 40us | 1.9ms | 48.2x | PASS |
| EXPFITA | 5 | 22 | Optimal | Optimal | 1.40e-03 | 28 | 13 | 452us | 3.4ms | 7.5x | MISMATCH |
| EXPFITB | 5 | 102 | Optimal | Optimal | 3.69e-08 | 92 | 16 | 5.0ms | 6.0ms | 1.2x | PASS |
| EXPFITC | 5 | 502 | MaxIteration | Optimal | N/A | 2999 | 18 | 1.34s | 16.6ms | 0.0x | ripopt_FAIL |
| EXPFITNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 61us | 575us | 9.4x | BOTH_FAIL |
| EXTRASIM | 2 | 1 | Optimal | Optimal | 1.07e-08 | 4 | 3 | 34us | 1.4ms | 41.2x | PASS |
| FBRAIN | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 5.4ms | 886us | 0.2x | BOTH_FAIL |
| FBRAIN2 | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 41.1ms | 1.1ms | 0.0x | BOTH_FAIL |
| FBRAIN2LS | 4 | 0 | Optimal | Optimal | 5.11e-10 | 10 | 10 | 7.4ms | 9.4ms | 1.3x | PASS |
| FBRAIN2NE | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 12.7ms | 1.1ms | 0.1x | BOTH_FAIL |
| FBRAIN3 | 6 | 2211 | Optimal | IpoptStatus( | N/A | 1 | 0 | 5.45s | 1.3ms | 0.0x | ipopt_FAIL |
| FBRAIN3LS | 6 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 3.65s | 3.76s | 1.0x | BOTH_FAIL |
| FBRAINLS | 2 | 0 | Optimal | Optimal | 6.11e-16 | 9 | 7 | 3.6ms | 4.4ms | 1.2x | PASS |
| FBRAINNE | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 6.5ms | 1.0ms | 0.2x | BOTH_FAIL |
| FCCU | 19 | 8 | Optimal | Optimal | 1.75e-15 | 10 | 9 | 133us | 2.5ms | 18.6x | PASS |
| FEEDLOC | 90 | 259 | Optimal | Optimal | 1.10e-08 | 124 | 23 | 69.6ms | 15.1ms | 0.2x | PASS |
| FLETCHER | 4 | 4 | Optimal | Optimal | 1.11e-08 | 61 | 28 | 310us | 5.8ms | 18.8x | PASS |
| FLT | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 42us | 1.6ms | 36.8x | PASS |
| GAUSS1 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 706us | 655us | 0.9x | BOTH_FAIL |
| GAUSS1LS | 8 | 0 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 459us | 1.9ms | 4.1x | PASS |
| GAUSS2 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 1.3ms | 601us | 0.4x | BOTH_FAIL |
| GAUSS2LS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 10 | 5 | 1.0ms | 1.8ms | 1.8x | PASS |
| GAUSS3 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 926us | 647us | 0.7x | BOTH_FAIL |
| GAUSS3LS | 8 | 0 | Optimal | Optimal | 3.65e-16 | 7 | 11 | 586us | 3.2ms | 5.4x | PASS |
| GAUSSIAN | 3 | 0 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 32us | 1.0ms | 31.9x | PASS |
| GBRAIN | 2 | 2200 | LocalInfeasi | IpoptStatus( | N/A | 6 | 0 | 4.7ms | 944us | 0.2x | BOTH_FAIL |
| GBRAINLS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 2.5ms | 3.7ms | 1.5x | PASS |
| GENHS28 | 10 | 8 | Optimal | Optimal | 1.22e-15 | 1 | 1 | 44us | 1.1ms | 24.3x | PASS |
| GIGOMEZ1 | 3 | 3 | Optimal | Optimal | 1.00e+00 | 11 | 13 | 66us | 3.1ms | 46.5x | MISMATCH |
| GIGOMEZ2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 21 | 7 | 149us | 2.1ms | 14.0x | MISMATCH |
| GIGOMEZ3 | 3 | 3 | Optimal | Optimal | 4.08e-09 | 10 | 8 | 63us | 2.3ms | 37.0x | PASS |
| GOFFIN | 51 | 50 | Acceptable | Optimal | 1.00e+00 | 2999 | 7 | 6.25s | 4.4ms | 0.0x | MISMATCH |
| GOTTFR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 5 | 60us | 1.6ms | 26.0x | PASS |
| GOULDQP1 | 32 | 17 | Acceptable | Optimal | 1.52e-07 | 55 | 15 | 1.7ms | 3.6ms | 2.1x | PASS |
| GROUPING | 100 | 125 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 2.1ms | 623us | 0.3x | ipopt_FAIL |
| GROWTH | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 249us | 512us | 2.1x | BOTH_FAIL |
| GROWTHLS | 3 | 0 | Optimal | Optimal | 1.77e-15 | 43 | 71 | 163us | 12.2ms | 74.9x | PASS |
| GULF | 3 | 0 | Optimal | Optimal | 3.04e-22 | 22 | 28 | 672us | 5.8ms | 8.6x | PASS |
| GULFNE | 3 | 99 | Optimal | IpoptStatus( | N/A | 22 | 0 | 1.2ms | 588us | 0.5x | ipopt_FAIL |
| HAHN1 | 7 | 236 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 940.4ms | 684us | 0.0x | BOTH_FAIL |
| HAHN1LS | 7 | 0 | MaxIteration | Optimal | N/A | 2999 | 78 | 235.4ms | 17.4ms | 0.1x | ripopt_FAIL |
| HAIFAM | 99 | 150 | Acceptable | Optimal | 1.47e-06 | 2726 | 40 | 559.6ms | 15.6ms | 0.0x | PASS |
| HAIFAS | 13 | 9 | Optimal | Optimal | 9.97e-09 | 28 | 16 | 325us | 4.0ms | 12.4x | PASS |
| HAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 24 | 62 | 82us | 10.8ms | 131.1x | PASS |
| HALDMADS | 6 | 42 | Optimal | Optimal | 6.31e-01 | 157 | 8 | 12.2ms | 3.3ms | 0.3x | MISMATCH |
| HART6 | 6 | 0 | Optimal | Optimal | 1.34e-16 | 10 | 7 | 64us | 2.3ms | 35.5x | PASS |
| HATFLDA | 4 | 0 | Optimal | Optimal | 1.58e-13 | 9 | 13 | 44us | 2.9ms | 66.6x | PASS |
| HATFLDANE | 4 | 4 | Acceptable | Optimal | 0.00e+00 | 9 | 6 | 66us | 2.0ms | 30.1x | PASS |
| HATFLDB | 4 | 0 | Optimal | Optimal | 3.90e-09 | 9 | 8 | 43us | 2.2ms | 49.8x | PASS |
| HATFLDBNE | 4 | 4 | LocalInfeasi | Infeasible | N/A | 12 | 13 | 77us | 3.5ms | 45.9x | BOTH_FAIL |
| HATFLDC | 25 | 0 | Optimal | Optimal | 6.53e-16 | 9 | 5 | 102us | 1.8ms | 17.4x | PASS |
| HATFLDCNE | 25 | 25 | Acceptable | Optimal | 0.00e+00 | 9 | 4 | 177us | 1.7ms | 9.5x | PASS |
| HATFLDD | 3 | 0 | Optimal | Optimal | 1.30e-19 | 21 | 21 | 86us | 3.6ms | 41.9x | PASS |
| HATFLDDNE | 3 | 10 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 131us | 586us | 4.5x | BOTH_FAIL |
| HATFLDE | 3 | 0 | Optimal | Optimal | 2.72e-20 | 19 | 20 | 98us | 3.3ms | 33.7x | PASS |
| HATFLDENE | 3 | 21 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 181us | 607us | 3.3x | BOTH_FAIL |
| HATFLDF | 3 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 135 | 83us | 23.8ms | 286.5x | PASS |
| HATFLDFL | 3 | 0 | Optimal | Optimal | 1.03e-08 | 494 | 1281 | 744us | 189.3ms | 254.4x | PASS |
| HATFLDFLNE | 3 | 3 | LocalInfeasi | Optimal | N/A | 486 | 15 | 1.0ms | 3.5ms | 3.4x | ripopt_FAIL |
| HATFLDFLS | 3 | 0 | Optimal | Optimal | 3.79e-18 | 18 | 36 | 66us | 6.1ms | 93.0x | PASS |
| HATFLDG | 25 | 25 | Optimal | Optimal | 0.00e+00 | 13 | 7 | 285us | 2.1ms | 7.3x | PASS |
| HATFLDGLS | 25 | 0 | Optimal | Optimal | 1.37e-16 | 13 | 14 | 187us | 3.0ms | 16.2x | PASS |
| HATFLDH | 4 | 7 | MaxIteration | Optimal | N/A | 2999 | 17 | 9.7ms | 3.8ms | 0.4x | ripopt_FAIL |
| HEART6 | 6 | 6 | Optimal | Optimal | 0.00e+00 | 956 | 22 | 5.1ms | 5.9ms | 1.2x | PASS |
| HEART6LS | 6 | 0 | Optimal | Optimal | 9.46e-23 | 1544 | 875 | 5.6ms | 140.5ms | 25.0x | PASS |
| HEART8 | 8 | 8 | Optimal | Optimal | 0.00e+00 | 66 | 12 | 388us | 2.9ms | 7.5x | PASS |
| HEART8LS | 8 | 0 | Optimal | Optimal | 3.37e-25 | 67 | 106 | 287us | 17.9ms | 62.4x | PASS |
| HELIX | 3 | 0 | Optimal | Optimal | 6.06e-25 | 8 | 13 | 51us | 2.8ms | 54.4x | PASS |
| HELIXNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 12 | 7 | 72us | 1.8ms | 24.7x | PASS |
| HET-Z | 2 | 1002 | Optimal | Optimal | 4.60e-05 | 55 | 11 | 30.7ms | 19.5ms | 0.6x | PASS |
| HIELOW | 3 | 0 | Optimal | Optimal | 5.46e-15 | 7 | 8 | 9.5ms | 11.7ms | 1.2x | PASS |
| HIMMELBA | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 34us | 975us | 28.6x | PASS |
| HIMMELBB | 2 | 0 | Optimal | Optimal | 1.40e-17 | 10 | 18 | 50us | 3.3ms | 66.2x | PASS |
| HIMMELBC | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 55us | 1.7ms | 30.8x | PASS |
| HIMMELBCLS | 2 | 0 | Optimal | Optimal | 5.80e-25 | 8 | 6 | 45us | 1.7ms | 37.2x | PASS |
| HIMMELBD | 2 | 2 | LocalInfeasi | Infeasible | N/A | 20 | 22 | 78us | 5.2ms | 67.1x | BOTH_FAIL |
| HIMMELBE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 5 | 2 | 48us | 1.2ms | 24.4x | PASS |
| HIMMELBF | 4 | 0 | Acceptable | Optimal | 3.57e-15 | 44 | 75 | 160us | 11.6ms | 72.4x | PASS |
| HIMMELBFNE | 4 | 7 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 208us | 596us | 2.9x | BOTH_FAIL |
| HIMMELBG | 2 | 0 | Optimal | Optimal | 3.63e-22 | 6 | 6 | 40us | 1.8ms | 45.7x | PASS |
| HIMMELBH | 2 | 0 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 34us | 1.6ms | 46.1x | PASS |
| HIMMELBI | 100 | 12 | Optimal | Optimal | 1.13e-10 | 22 | 13 | 913us | 4.0ms | 4.4x | PASS |
| HIMMELBJ | 45 | 14 | Acceptable | ErrorInStepC | N/A | 31 | 580 | 2.2ms | 138.4ms | 63.2x | ipopt_FAIL |
| HIMMELBK | 24 | 14 | Optimal | Optimal | 5.04e-08 | 17 | 18 | 519us | 4.3ms | 8.3x | PASS |
| HIMMELP1 | 2 | 0 | Optimal | Optimal | 1.83e-15 | 11 | 10 | 47us | 2.6ms | 55.4x | PASS |
| HIMMELP2 | 2 | 1 | Optimal | Optimal | 8.68e-01 | 12 | 17 | 67us | 4.1ms | 61.4x | MISMATCH |
| HIMMELP3 | 2 | 2 | Optimal | Optimal | 8.66e-01 | 0 | 11 | 22us | 2.8ms | 125.1x | MISMATCH |
| HIMMELP4 | 2 | 3 | Optimal | Optimal | 1.23e-08 | 9 | 23 | 75us | 4.9ms | 65.9x | PASS |
| HIMMELP5 | 2 | 3 | Acceptable | Optimal | 6.39e-01 | 10 | 46 | 75us | 8.6ms | 115.5x | MISMATCH |
| HIMMELP6 | 2 | 5 | Optimal | Optimal | 1.22e-08 | 13 | 31 | 113us | 6.5ms | 57.7x | PASS |
| HONG | 4 | 1 | Optimal | Optimal | 7.87e-16 | 9 | 7 | 59us | 2.1ms | 34.6x | PASS |
| HS1 | 2 | 0 | Optimal | Optimal | 9.22e-21 | 25 | 28 | 58us | 5.6ms | 96.2x | PASS |
| HS10 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 12 | 72us | 2.8ms | 38.7x | PASS |
| HS100 | 7 | 4 | Optimal | Optimal | 1.77e-11 | 14 | 9 | 122us | 2.5ms | 20.7x | PASS |
| HS100LNP | 7 | 2 | Optimal | Optimal | 1.67e-16 | 6 | 20 | 69us | 3.3ms | 47.3x | PASS |
| HS100MOD | 7 | 4 | Optimal | Optimal | 6.73e-13 | 15 | 14 | 131us | 3.3ms | 25.1x | PASS |
| HS101 | 7 | 5 | Optimal | Optimal | 4.60e-08 | 29 | 39 | 1.7ms | 9.7ms | 5.7x | PASS |
| HS102 | 7 | 5 | Optimal | Optimal | 3.85e-08 | 33 | 52 | 812us | 10.0ms | 12.3x | PASS |
| HS103 | 7 | 5 | Optimal | Optimal | 2.63e-08 | 26 | 21 | 326us | 4.8ms | 14.6x | PASS |
| HS104 | 8 | 5 | Optimal | Optimal | 2.45e-08 | 32 | 8 | 319us | 2.3ms | 7.1x | PASS |
| HS105 | 8 | 1 | Acceptable | Optimal | 5.34e-11 | 27 | 23 | 2.9ms | 7.2ms | 2.5x | PASS |
| HS106 | 8 | 6 | Optimal | Optimal | 1.74e-08 | 17 | 18 | 125us | 3.8ms | 30.5x | PASS |
| HS107 | 9 | 6 | Acceptable | Optimal | 3.11e-06 | 26 | 7 | 200us | 2.1ms | 10.6x | PASS |
| HS108 | 9 | 13 | Optimal | Optimal | 8.77e-09 | 2596 | 11 | 15.0ms | 3.2ms | 0.2x | PASS |
| HS109 | 9 | 10 | Acceptable | Optimal | 1.30e-08 | 29 | 14 | 276us | 3.3ms | 11.9x | PASS |
| HS11 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 6 | 43us | 1.8ms | 43.2x | PASS |
| HS111 | 10 | 3 | Optimal | Optimal | 1.80e-11 | 13 | 15 | 158us | 3.6ms | 22.5x | PASS |
| HS111LNP | 10 | 3 | Optimal | Optimal | 1.03e-10 | 13 | 15 | 161us | 2.9ms | 18.1x | PASS |
| HS112 | 10 | 3 | Optimal | Optimal | 6.49e-14 | 10 | 10 | 128us | 2.5ms | 19.6x | PASS |
| HS113 | 10 | 8 | Optimal | Optimal | 1.69e-09 | 35 | 9 | 345us | 2.6ms | 7.6x | PASS |
| HS114 | 10 | 11 | Optimal | Optimal | 1.07e-07 | 122 | 13 | 1.0ms | 3.2ms | 3.1x | PASS |
| HS116 | 13 | 14 | Acceptable | Optimal | 6.97e-04 | 297 | 19 | 4.5ms | 4.5ms | 1.0x | MISMATCH |
| HS117 | 15 | 5 | Acceptable | Optimal | 2.55e-06 | 19 | 19 | 209us | 4.5ms | 21.7x | PASS |
| HS118 | 15 | 17 | Acceptable | Optimal | 1.65e-02 | 2999 | 10 | 105.2ms | 2.8ms | 0.0x | MISMATCH |
| HS119 | 16 | 8 | Acceptable | Optimal | 8.50e-05 | 21 | 17 | 400us | 4.0ms | 10.1x | PASS |
| HS12 | 2 | 1 | Optimal | Optimal | 1.57e-10 | 10 | 6 | 63us | 1.9ms | 30.4x | PASS |
| HS13 | 2 | 1 | Acceptable | Optimal | 1.82e-04 | 32 | 47 | 96us | 8.3ms | 86.5x | MISMATCH |
| HS14 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 5 | 42us | 1.8ms | 41.8x | PASS |
| HS15 | 2 | 2 | Optimal | Optimal | 6.59e-08 | 12 | 13 | 67us | 2.9ms | 44.0x | PASS |
| HS16 | 2 | 2 | Acceptable | Optimal | 9.89e-01 | 11 | 10 | 68us | 2.7ms | 39.9x | MISMATCH |
| HS17 | 2 | 2 | Optimal | Optimal | 9.30e-06 | 15 | 22 | 73us | 4.5ms | 61.5x | PASS |
| HS18 | 2 | 2 | Optimal | Optimal | 7.79e-10 | 14 | 10 | 68us | 2.4ms | 35.4x | PASS |
| HS19 | 2 | 2 | Optimal | Optimal | 3.34e-09 | 18 | 12 | 88us | 3.0ms | 33.8x | PASS |
| HS1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 25 | 30 | 94us | 6.9ms | 73.6x | PASS |
| HS2 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 10 | 10 | 43us | 2.5ms | 58.7x | PASS |
| HS20 | 2 | 3 | Optimal | Optimal | 6.53e-08 | 11 | 5 | 67us | 1.7ms | 25.3x | PASS |
| HS21 | 2 | 1 | Optimal | Optimal | 8.01e-12 | 7 | 6 | 47us | 2.0ms | 41.6x | PASS |
| HS21MOD | 7 | 1 | Acceptable | Optimal | 1.74e-08 | 11 | 13 | 61us | 3.0ms | 48.9x | PASS |
| HS22 | 2 | 2 | Optimal | Optimal | 1.16e-08 | 5 | 5 | 45us | 1.8ms | 40.9x | PASS |
| HS23 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 271 | 9 | 876us | 2.3ms | 2.7x | MISMATCH |
| HS24 | 2 | 3 | Optimal | Optimal | 1.66e-08 | 6 | 14 | 60us | 3.4ms | 56.9x | PASS |
| HS25 | 3 | 0 | Acceptable | Optimal | 6.08e-05 | 16 | 27 | 337us | 6.2ms | 18.4x | PASS |
| HS25NE | 3 | 99 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 522us | 581us | 1.1x | BOTH_FAIL |
| HS26 | 3 | 1 | Acceptable | Optimal | 5.46e-12 | 17 | 25 | 67us | 3.7ms | 55.2x | PASS |
| HS268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 42us | 3.2ms | 76.8x | MISMATCH |
| HS27 | 3 | 1 | Optimal | Optimal | 1.73e-14 | 16 | 57 | 76us | 8.9ms | 117.3x | PASS |
| HS28 | 3 | 1 | Optimal | Optimal | 9.24e-31 | 1 | 1 | 32us | 1.0ms | 31.4x | PASS |
| HS29 | 3 | 1 | Optimal | Optimal | 3.12e-10 | 15 | 7 | 81us | 2.1ms | 26.5x | PASS |
| HS2NE | 2 | 2 | LocalInfeasi | Infeasible | N/A | 10 | 12 | 54us | 3.3ms | 61.4x | BOTH_FAIL |
| HS3 | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 37us | 1.5ms | 41.0x | PASS |
| HS30 | 3 | 1 | Acceptable | Optimal | 6.94e-09 | 9 | 7 | 55us | 2.0ms | 36.8x | PASS |
| HS31 | 3 | 1 | Optimal | Optimal | 8.94e-09 | 11 | 6 | 55us | 1.8ms | 33.5x | PASS |
| HS32 | 3 | 2 | Acceptable | Optimal | 6.20e-05 | 13 | 15 | 70us | 3.4ms | 48.2x | PASS |
| HS33 | 3 | 2 | Optimal | Optimal | 2.57e-08 | 12 | 9 | 77us | 2.4ms | 31.6x | PASS |
| HS34 | 3 | 2 | Optimal | Optimal | 3.19e-09 | 14 | 7 | 71us | 2.1ms | 28.9x | PASS |
| HS35 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 7 | 58us | 2.1ms | 36.5x | PASS |
| HS35I | 3 | 1 | Optimal | Optimal | 3.44e-09 | 9 | 7 | 56us | 2.1ms | 37.9x | PASS |
| HS35MOD | 3 | 1 | Acceptable | Optimal | 4.67e-07 | 8 | 14 | 48us | 3.0ms | 63.2x | PASS |
| HS36 | 3 | 1 | Optimal | Optimal | 6.33e-09 | 16 | 11 | 75us | 2.9ms | 38.8x | PASS |
| HS37 | 3 | 2 | Optimal | Optimal | 3.06e-10 | 9 | 11 | 65us | 3.0ms | 46.3x | PASS |
| HS38 | 4 | 0 | Optimal | Optimal | 4.52e-11 | 37 | 39 | 99us | 7.7ms | 77.5x | PASS |
| HS39 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 79us | 2.5ms | 31.7x | PASS |
| HS3MOD | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 35us | 1.5ms | 42.3x | PASS |
| HS4 | 2 | 0 | Optimal | Optimal | 1.97e-08 | 6 | 4 | 37us | 1.5ms | 40.9x | PASS |
| HS40 | 4 | 3 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 37us | 1.2ms | 31.1x | PASS |
| HS41 | 4 | 1 | Acceptable | Optimal | 1.41e-08 | 12 | 7 | 66us | 2.0ms | 31.0x | PASS |
| HS42 | 4 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 50us | 1.5ms | 29.2x | PASS |
| HS43 | 4 | 3 | Optimal | Optimal | 6.58e-03 | 8 | 8 | 73us | 2.3ms | 31.4x | MISMATCH |
| HS44 | 4 | 6 | Acceptable | Optimal | 7.91e-06 | 266 | 24 | 822us | 5.3ms | 6.4x | PASS |
| HS44NEW | 4 | 6 | Acceptable | Optimal | 1.55e-05 | 21 | 18 | 135us | 4.4ms | 32.9x | PASS |
| HS45 | 5 | 0 | Acceptable | Optimal | 3.14e-04 | 31 | 11 | 85us | 2.8ms | 33.0x | MISMATCH |
| HS46 | 5 | 2 | Acceptable | Optimal | 3.87e-11 | 17 | 19 | 92us | 3.1ms | 33.9x | PASS |
| HS47 | 5 | 3 | Acceptable | Optimal | 3.22e-12 | 17 | 19 | 92us | 3.0ms | 32.1x | PASS |
| HS48 | 5 | 2 | Optimal | Optimal | 4.44e-31 | 1 | 1 | 37us | 1.1ms | 29.1x | PASS |
| HS49 | 5 | 2 | Acceptable | Optimal | 2.61e-10 | 17 | 19 | 81us | 3.0ms | 36.8x | PASS |
| HS5 | 2 | 0 | Optimal | Optimal | 2.32e-16 | 9 | 7 | 42us | 2.1ms | 48.7x | PASS |
| HS50 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 64us | 1.9ms | 29.9x | PASS |
| HS51 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 42us | 1.0ms | 24.0x | PASS |
| HS52 | 5 | 3 | Optimal | Optimal | 1.17e-15 | 1 | 1 | 35us | 1.0ms | 28.6x | PASS |
| HS53 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 6 | 59us | 1.8ms | 30.2x | PASS |
| HS54 | 6 | 1 | Optimal | Optimal | 7.51e-01 | 10 | 15 | 75us | 3.4ms | 46.0x | MISMATCH |
| HS55 | 6 | 6 | Acceptable | Optimal | 2.04e-02 | 11 | 18 | 90us | 5.1ms | 56.3x | MISMATCH |
| HS56 | 7 | 4 | Optimal | Optimal | 1.21e-13 | 6 | 10 | 72us | 2.3ms | 31.2x | PASS |
| HS57 | 2 | 1 | Optimal | Optimal | 1.06e-13 | 9 | 10 | 69us | 2.3ms | 34.0x | PASS |
| HS59 | 2 | 3 | Optimal | Optimal | 1.34e+00 | 9 | 17 | 77us | 4.0ms | 52.3x | MISMATCH |
| HS6 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 62us | 1.8ms | 28.1x | PASS |
| HS60 | 3 | 1 | Optimal | Optimal | 1.19e-13 | 8 | 6 | 53us | 1.9ms | 35.6x | PASS |
| HS61 | 3 | 2 | Optimal | Optimal | 7.91e-16 | 9 | 10 | 59us | 2.0ms | 33.4x | PASS |
| HS62 | 3 | 1 | Optimal | Optimal | 1.38e-16 | 9 | 6 | 58us | 2.0ms | 34.7x | PASS |
| HS63 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 5 | 60us | 1.7ms | 28.3x | PASS |
| HS64 | 3 | 1 | Optimal | Optimal | 3.58e-09 | 19 | 16 | 79us | 3.6ms | 45.3x | PASS |
| HS65 | 3 | 1 | Optimal | Optimal | 3.63e-09 | 12 | 16 | 60us | 3.8ms | 63.3x | PASS |
| HS66 | 3 | 2 | Optimal | Optimal | 1.12e-08 | 10 | 10 | 61us | 2.4ms | 38.5x | PASS |
| HS67 | 3 | 14 | Optimal | Optimal | 2.29e-01 | 24 | 9 | 320us | 2.5ms | 7.8x | MISMATCH |
| HS68 | 4 | 2 | Optimal | Optimal | 4.76e-12 | 18 | 16 | 95us | 3.6ms | 38.1x | PASS |
| HS69 | 4 | 2 | Optimal | Optimal | 2.61e-15 | 11 | 10 | 66us | 2.6ms | 39.4x | PASS |
| HS7 | 2 | 1 | Optimal | Optimal | 5.30e-12 | 10 | 27 | 56us | 5.0ms | 88.7x | PASS |
| HS70 | 4 | 1 | Optimal | Optimal | 1.72e-01 | 11 | 46 | 170us | 8.8ms | 51.6x | MISMATCH |
| HS71 | 4 | 2 | Optimal | Optimal | 8.97e-10 | 11 | 8 | 70us | 2.3ms | 33.4x | PASS |
| HS72 | 4 | 2 | Optimal | Optimal | 3.14e-07 | 13 | 16 | 69us | 3.5ms | 49.9x | PASS |
| HS73 | 4 | 3 | Optimal | Optimal | 1.24e-09 | 14 | 8 | 79us | 2.3ms | 28.8x | PASS |
| HS74 | 4 | 5 | Optimal | Optimal | 1.77e-16 | 16 | 8 | 114us | 2.2ms | 19.6x | PASS |
| HS75 | 4 | 5 | Optimal | Optimal | 5.37e-09 | 16 | 8 | 100us | 2.2ms | 22.0x | PASS |
| HS76 | 4 | 3 | Optimal | Optimal | 6.09e-09 | 9 | 7 | 62us | 2.1ms | 34.4x | PASS |
| HS76I | 4 | 3 | Optimal | Optimal | 3.36e-08 | 9 | 6 | 63us | 1.9ms | 30.2x | PASS |
| HS77 | 5 | 2 | Optimal | Optimal | 1.99e-11 | 9 | 11 | 68us | 2.3ms | 33.7x | PASS |
| HS78 | 5 | 3 | Optimal | Optimal | 1.52e-16 | 4 | 4 | 48us | 1.4ms | 28.5x | PASS |
| HS79 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 50us | 1.4ms | 28.0x | PASS |
| HS8 | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 45us | 1.5ms | 33.9x | PASS |
| HS80 | 5 | 3 | Optimal | Optimal | 5.69e-14 | 10 | 5 | 73us | 1.7ms | 23.6x | PASS |
| HS81 | 5 | 3 | Optimal | Optimal | 3.85e-01 | 111 | 68 | 509us | 12.3ms | 24.2x | MISMATCH |
| HS83 | 5 | 3 | MaxIteration | Optimal | N/A | 2999 | 9 | 8.8ms | 2.5ms | 0.3x | ripopt_FAIL |
| HS84 | 5 | 3 | MaxIteration | Optimal | N/A | 2999 | 9 | 118.7ms | 2.6ms | 0.0x | ripopt_FAIL |
| HS85 | 5 | 21 | Acceptable | Optimal | 9.51e-02 | 55 | 13 | 15.4ms | 4.3ms | 0.3x | MISMATCH |
| HS86 | 5 | 10 | Acceptable | Optimal | 5.19e-09 | 13 | 10 | 136us | 2.8ms | 20.4x | PASS |
| HS87 | 6 | 4 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 10.9ms | 466.2ms | 42.7x | BOTH_FAIL |
| HS88 | 2 | 1 | Optimal | Optimal | 7.70e-06 | 16 | 18 | 713us | 4.4ms | 6.2x | PASS |
| HS89 | 3 | 1 | Optimal | Optimal | 7.69e-06 | 100 | 15 | 24.1ms | 4.4ms | 0.2x | PASS |
| HS9 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 41us | 1.4ms | 34.0x | PASS |
| HS90 | 4 | 1 | Optimal | Optimal | 5.48e-06 | 17 | 16 | 1.2ms | 4.5ms | 3.8x | PASS |
| HS91 | 5 | 1 | Optimal | Optimal | 7.70e-06 | 18 | 16 | 1.5ms | 4.9ms | 3.2x | PASS |
| HS92 | 6 | 1 | Optimal | Optimal | 6.87e-06 | 16 | 35 | 1.8ms | 10.4ms | 5.8x | PASS |
| HS93 | 6 | 2 | Optimal | Optimal | 9.89e-09 | 80 | 7 | 643us | 2.2ms | 3.4x | PASS |
| HS95 | 6 | 4 | Acceptable | Optimal | 5.90e-04 | 16 | 9 | 116us | 2.5ms | 21.8x | MISMATCH |
| HS96 | 6 | 4 | Acceptable | Optimal | 6.28e-04 | 16 | 8 | 110us | 2.3ms | 20.9x | MISMATCH |
| HS97 | 6 | 4 | Acceptable | Optimal | 2.30e-01 | 20 | 24 | 143us | 5.0ms | 34.9x | MISMATCH |
| HS98 | 6 | 4 | Acceptable | Optimal | 2.30e-01 | 17 | 13 | 132us | 3.0ms | 22.6x | MISMATCH |
| HS99 | 7 | 2 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 74us | 1.8ms | 24.0x | PASS |
| HS99EXP | 31 | 21 | Acceptable | Optimal | 1.86e-14 | 17 | 17 | 913us | 4.0ms | 4.3x | PASS |
| HUBFIT | 2 | 1 | Optimal | Optimal | 2.61e-09 | 7 | 7 | 46us | 2.0ms | 43.5x | PASS |
| HUMPS | 2 | 0 | Optimal | Optimal | 2.76e-17 | 38 | 1533 | 118us | 214.0ms | 1815.1x | PASS |
| HYDC20LS | 99 | 0 | Acceptable | Optimal | 2.98e-01 | 2999 | 639 | 1.01s | 165.4ms | 0.2x | MISMATCH |
| HYDCAR20 | 99 | 99 | LocalInfeasi | Optimal | N/A | 2999 | 9 | 1.73s | 3.3ms | 0.0x | ripopt_FAIL |
| HYDCAR6 | 29 | 29 | Optimal | Optimal | 0.00e+00 | 1298 | 5 | 38.5ms | 1.9ms | 0.1x | PASS |
| HYDCAR6LS | 29 | 0 | Optimal | Optimal | 2.86e-18 | 1322 | 149 | 28.2ms | 30.0ms | 1.1x | PASS |
| HYPCIR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 53us | 1.7ms | 31.8x | PASS |
| JENSMP | 2 | 0 | Optimal | Optimal | 3.43e-16 | 10 | 9 | 50us | 1.9ms | 37.9x | PASS |
| JENSMPNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 78us | 568us | 7.3x | BOTH_FAIL |
| JUDGE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 57us | 1.9ms | 33.4x | PASS |
| JUDGEB | 2 | 0 | Optimal | Optimal | 2.21e-16 | 9 | 9 | 50us | 2.2ms | 44.6x | PASS |
| JUDGENE | 2 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 80us | 576us | 7.2x | BOTH_FAIL |
| KIRBY2 | 5 | 151 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 979us | 547us | 0.6x | BOTH_FAIL |
| KIRBY2LS | 5 | 0 | Acceptable | Optimal | 1.21e-14 | 26 | 11 | 650us | 2.7ms | 4.2x | PASS |
| KIWCRESC | 3 | 2 | Optimal | Optimal | 1.00e-08 | 21 | 8 | 90us | 2.3ms | 25.4x | PASS |
| KOEBHELB | 3 | 0 | Optimal | Optimal | 7.33e-16 | 843 | 71 | 12.6ms | 15.6ms | 1.2x | PASS |
| KOEBHELBNE | 3 | 156 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 2.2ms | 606us | 0.3x | BOTH_FAIL |
| KOWOSB | 4 | 0 | Optimal | Optimal | 5.96e-19 | 7 | 8 | 55us | 2.3ms | 41.0x | PASS |
| KOWOSBNE | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 80us | 585us | 7.3x | BOTH_FAIL |
| KSIP | 20 | 1001 | Optimal | Optimal | 9.79e-01 | 1 | 22 | 13.3ms | 69.3ms | 5.2x | MISMATCH |
| LAKES | 90 | 78 | Optimal | Optimal | 1.36e-13 | 16 | 11 | 1.4ms | 3.9ms | 2.9x | PASS |
| LANCZOS1 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 41 | 0 | 473us | 606us | 1.3x | ipopt_FAIL |
| LANCZOS1LS | 6 | 0 | Acceptable | Optimal | 5.33e-09 | 105 | 115 | 670us | 19.3ms | 28.8x | PASS |
| LANCZOS2 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 70 | 0 | 704us | 572us | 0.8x | ipopt_FAIL |
| LANCZOS2LS | 6 | 0 | Acceptable | Optimal | 7.41e-09 | 97 | 101 | 596us | 17.6ms | 29.6x | PASS |
| LANCZOS3 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 30 | 0 | 343us | 570us | 1.7x | ipopt_FAIL |
| LANCZOS3LS | 6 | 0 | Acceptable | Optimal | 1.08e-08 | 78 | 174 | 494us | 29.5ms | 59.8x | PASS |
| LAUNCH | 25 | 28 | Acceptable | Optimal | 1.13e-05 | 13 | 12 | 540us | 3.6ms | 6.6x | PASS |
| LEVYMONE10 | 10 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 110us | 580us | 5.3x | BOTH_FAIL |
| LEVYMONE5 | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 66us | 537us | 8.1x | BOTH_FAIL |
| LEVYMONE6 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 71us | 537us | 7.6x | BOTH_FAIL |
| LEVYMONE7 | 4 | 8 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 80us | 584us | 7.3x | BOTH_FAIL |
| LEVYMONE8 | 5 | 10 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 73us | 549us | 7.5x | BOTH_FAIL |
| LEVYMONE9 | 8 | 16 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 95us | 579us | 6.1x | BOTH_FAIL |
| LEVYMONT10 | 10 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 4 | 60us | 1.5ms | 25.3x | PASS |
| LEVYMONT5 | 2 | 0 | Acceptable | Optimal | 1.00e+00 | 10 | 10 | 53us | 2.7ms | 50.7x | MISMATCH |
| LEVYMONT6 | 3 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 8 | 55us | 2.3ms | 41.3x | PASS |
| LEVYMONT7 | 4 | 0 | Optimal | Optimal | 1.42e-16 | 10 | 7 | 54us | 2.3ms | 43.6x | PASS |
| LEVYMONT8 | 5 | 0 | Optimal | Optimal | 1.64e-16 | 8 | 4 | 47us | 1.5ms | 31.5x | PASS |
| LEVYMONT9 | 8 | 0 | Optimal | Optimal | 1.71e-16 | 8 | 4 | 55us | 1.5ms | 27.5x | PASS |
| LEWISPOL | 6 | 9 | Acceptable | IpoptStatus( | N/A | 10 | 0 | 982us | 548us | 0.6x | ipopt_FAIL |
| LHAIFAM | 99 | 150 | MaxIteration | InvalidNumbe | N/A | 2999 | 0 | 905.7ms | 740us | 0.0x | BOTH_FAIL |
| LIN | 4 | 2 | Optimal | Optimal | 2.03e-03 | 9 | 7 | 76us | 2.1ms | 27.2x | MISMATCH |
| LINSPANH | 97 | 33 | Acceptable | Optimal | 5.76e-07 | 2999 | 24 | 110.2ms | 6.0ms | 0.1x | PASS |
| LOADBAL | 31 | 31 | Optimal | Optimal | 1.06e-08 | 15 | 13 | 848us | 3.6ms | 4.3x | PASS |
| LOGHAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 55 | 2747 | 152us | 389.1ms | 2553.4x | PASS |
| LOGROS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 50 | 49 | 101us | 10.1ms | 100.0x | PASS |
| LOOTSMA | 3 | 2 | Optimal | Optimal | 7.42e-04 | 11 | 13 | 76us | 3.3ms | 43.3x | MISMATCH |
| LOTSCHD | 12 | 7 | Optimal | Optimal | 4.44e-10 | 14 | 9 | 139us | 2.5ms | 17.6x | PASS |
| LRCOVTYPE | 54 | 0 | Optimal | Optimal | 1.78e-02 | 61 | 33 | 13.74s | 6.38s | 0.5x | MISMATCH |
| LRIJCNN1 | 22 | 0 | Optimal | Optimal | 1.16e-14 | 18 | 11 | 340.8ms | 184.0ms | 0.5x | PASS |
| LSC1 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 91us | 580us | 6.3x | BOTH_FAIL |
| LSC1LS | 3 | 0 | Acceptable | Optimal | 2.30e-15 | 17 | 16 | 72us | 3.3ms | 46.1x | PASS |
| LSC2 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 136us | 605us | 4.5x | BOTH_FAIL |
| LSC2LS | 3 | 0 | Acceptable | Optimal | 2.32e-04 | 32 | 38 | 85us | 5.5ms | 64.5x | MISMATCH |
| LSNNODOC | 5 | 4 | Acceptable | Optimal | 3.42e-06 | 13 | 10 | 95us | 2.7ms | 28.4x | PASS |
| LSQFIT | 2 | 1 | Optimal | Optimal | 4.47e-09 | 7 | 7 | 48us | 2.0ms | 41.8x | PASS |
| MADSEN | 3 | 6 | Optimal | Optimal | 8.18e-09 | 18 | 18 | 138us | 4.0ms | 28.7x | PASS |
| MAKELA1 | 3 | 2 | Optimal | Optimal | 2.00e+00 | 6 | 12 | 59us | 3.1ms | 51.8x | MISMATCH |
| MAKELA2 | 3 | 3 | Optimal | Optimal | 1.27e-01 | 3 | 6 | 39us | 1.9ms | 47.3x | MISMATCH |
| MAKELA3 | 21 | 20 | Acceptable | Optimal | 3.19e-09 | 137 | 11 | 2.3ms | 3.0ms | 1.3x | PASS |
| MAKELA4 | 21 | 40 | Optimal | Optimal | 6.88e-01 | 1 | 5 | 104us | 2.0ms | 19.0x | MISMATCH |
| MARATOS | 2 | 1 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 39us | 1.4ms | 36.6x | PASS |
| MARATOSB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 670 | 672 | 803us | 99.5ms | 124.0x | PASS |
| MATRIX2 | 6 | 2 | Acceptable | Optimal | 3.77e-09 | 17 | 42 | 91us | 7.8ms | 85.8x | PASS |
| MAXLIKA | 8 | 0 | Acceptable | Optimal | 1.13e-02 | 20 | 23 | 2.0ms | 7.2ms | 3.6x | MISMATCH |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.03e-08 | 31 | 7 | 318us | 2.2ms | 7.0x | PASS |
| MDHOLE | 2 | 0 | Optimal | Optimal | 9.98e-09 | 35 | 42 | 84us | 8.9ms | 105.7x | PASS |
| MESH | 41 | 48 | Acceptable | IpoptStatus( | N/A | 49 | 79 | 12.2ms | 20.9ms | 1.7x | ipopt_FAIL |
| METHANB8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 9 | 3 | 357us | 1.5ms | 4.3x | PASS |
| METHANB8LS | 31 | 0 | Optimal | Optimal | 5.35e-26 | 9 | 8 | 246us | 2.1ms | 8.4x | PASS |
| METHANL8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 608 | 4 | 19.9ms | 1.7ms | 0.1x | PASS |
| METHANL8LS | 31 | 0 | Optimal | Optimal | 5.70e-17 | 977 | 40 | 21.7ms | 9.1ms | 0.4x | PASS |
| MEXHAT | 2 | 0 | Optimal | Optimal | 6.77e-10 | 28 | 26 | 73us | 4.3ms | 58.6x | PASS |
| MEYER3 | 3 | 0 | Acceptable | Optimal | 2.23e-12 | 2999 | 194 | 11.2ms | 30.5ms | 2.7x | PASS |
| MEYER3NE | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 20.7ms | 651us | 0.0x | BOTH_FAIL |
| MGH09 | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 302us | 565us | 1.9x | BOTH_FAIL |
| MGH09LS | 4 | 0 | Optimal | Optimal | 5.96e-19 | 58 | 72 | 194us | 12.0ms | 61.9x | PASS |
| MGH10 | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 18.1ms | 599us | 0.0x | BOTH_FAIL |
| MGH10LS | 3 | 0 | Acceptable | Optimal | 4.94e-12 | 2999 | 1828 | 10.9ms | 276.6ms | 25.4x | PASS |
| MGH10S | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 14.5ms | 581us | 0.0x | BOTH_FAIL |
| MGH10SLS | 3 | 0 | MaxIteration | Optimal | N/A | 2999 | 354 | 8.6ms | 53.9ms | 6.3x | ripopt_FAIL |
| MGH17 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 378us | 585us | 1.5x | BOTH_FAIL |
| MGH17LS | 5 | 0 | Acceptable | Optimal | 6.06e-07 | 37 | 47 | 271us | 9.3ms | 34.2x | PASS |
| MGH17S | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 750us | 597us | 0.8x | BOTH_FAIL |
| MGH17SLS | 5 | 0 | Optimal | Optimal | 2.45e-02 | 54 | 41 | 386us | 8.3ms | 21.5x | MISMATCH |
| MIFFLIN1 | 3 | 2 | Optimal | Optimal | 6.15e-09 | 13 | 5 | 83us | 1.9ms | 22.5x | PASS |
| MIFFLIN2 | 3 | 2 | Optimal | Optimal | 9.97e-09 | 15 | 11 | 83us | 3.0ms | 35.7x | PASS |
| MINMAXBD | 5 | 20 | MaxIteration | Optimal | N/A | 2999 | 25 | 60.5ms | 6.4ms | 0.1x | ripopt_FAIL |
| MINMAXRB | 3 | 4 | Optimal | Optimal | 9.85e-09 | 3 | 8 | 43us | 2.3ms | 52.8x | PASS |
| MINSURF | 64 | 0 | Acceptable | Optimal | 0.00e+00 | 5 | 4 | 384us | 1.8ms | 4.8x | PASS |
| MISRA1A | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 171us | 579us | 3.4x | BOTH_FAIL |
| MISRA1ALS | 2 | 0 | Acceptable | Optimal | 1.91e-14 | 35 | 40 | 128us | 6.7ms | 52.2x | PASS |
| MISRA1B | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 171us | 556us | 3.3x | BOTH_FAIL |
| MISRA1BLS | 2 | 0 | Optimal | Optimal | 4.76e-14 | 25 | 34 | 92us | 5.6ms | 61.6x | PASS |
| MISRA1C | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 131us | 560us | 4.3x | BOTH_FAIL |
| MISRA1CLS | 2 | 0 | Acceptable | Optimal | 4.63e-14 | 18 | 14 | 89us | 3.1ms | 34.9x | PASS |
| MISRA1D | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 24 | 0 | 140us | 568us | 4.1x | BOTH_FAIL |
| MISRA1DLS | 2 | 0 | Optimal | Optimal | 9.02e-17 | 20 | 30 | 80us | 5.1ms | 63.8x | PASS |
| MISTAKE | 9 | 13 | Acceptable | Optimal | 5.00e-01 | 81 | 16 | 1.0ms | 4.0ms | 3.9x | MISMATCH |
| MRIBASIS | 36 | 55 | Acceptable | Optimal | 1.00e-08 | 29 | 15 | 4.3ms | 4.9ms | 1.1x | PASS |
| MSS1 | 90 | 73 | Acceptable | Optimal | 1.25e-01 | 81 | 95 | 39.7ms | 50.2ms | 1.3x | MISMATCH |
| MUONSINE | 1 | 512 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 1.5ms | 613us | 0.4x | BOTH_FAIL |
| MUONSINELS | 1 | 0 | Acceptable | Optimal | 1.31e-01 | 10 | 8 | 922us | 2.3ms | 2.5x | MISMATCH |
| MWRIGHT | 5 | 3 | Optimal | Optimal | 9.48e-01 | 13 | 10 | 80us | 2.2ms | 27.1x | MISMATCH |
| NASH | 72 | 24 | RestorationF | Infeasible | N/A | 180 | 45 | 83.8ms | 12.8ms | 0.2x | BOTH_FAIL |
| NELSON | 3 | 128 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 110.0ms | 683us | 0.0x | BOTH_FAIL |
| NET1 | 48 | 57 | Acceptable | Optimal | 7.74e-08 | 30 | 26 | 70.3ms | 6.6ms | 0.1x | PASS |
| NYSTROM5 | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 353us | 574us | 1.6x | ipopt_FAIL |
| NYSTROM5C | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 373us | 622us | 1.7x | ipopt_FAIL |
| ODFITS | 10 | 6 | Optimal | Optimal | 1.91e-16 | 11 | 8 | 99us | 2.3ms | 22.9x | PASS |
| OET1 | 3 | 1002 | Optimal | Optimal | 8.40e-02 | 113 | 33 | 86.1ms | 47.8ms | 0.6x | MISMATCH |
| OET2 | 3 | 1002 | Acceptable | Optimal | 8.54e-09 | 1681 | 181 | 1.75s | 278.0ms | 0.2x | PASS |
| OET3 | 4 | 1002 | Optimal | Optimal | 9.88e-09 | 53 | 13 | 35.8ms | 21.7ms | 0.6x | PASS |
| OET4 | 4 | 1002 | Optimal | Optimal | 7.66e-09 | 409 | 165 | 391.4ms | 252.9ms | 0.6x | PASS |
| OET5 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET6 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET7 | 7 | 1002 | MaxIteration | Optimal | N/A | 2999 | 193 | 54.92s | 530.7ms | 0.0x | ripopt_FAIL |
| OPTCNTRL | 32 | 20 | Acceptable | Optimal | 1.64e-08 | 40 | 9 | 1.4ms | 2.6ms | 1.9x | PASS |
| OPTPRLOC | 30 | 30 | Acceptable | Optimal | 2.27e-07 | 29 | 13 | 1.3ms | 3.7ms | 2.9x | PASS |
| ORTHREGB | 27 | 6 | Optimal | Optimal | 4.20e-19 | 2 | 2 | 88us | 1.2ms | 14.0x | PASS |
| OSBORNE1 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 2982 | 0 | 30.5ms | 616us | 0.0x | BOTH_FAIL |
| OSBORNE2 | 11 | 65 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 718us | 570us | 0.8x | BOTH_FAIL |
| OSBORNEA | 5 | 0 | MaxIteration | Optimal | N/A | 2999 | 64 | 18.4ms | 10.8ms | 0.6x | ripopt_FAIL |
| OSBORNEB | 11 | 0 | Optimal | Optimal | 4.86e-17 | 16 | 19 | 493us | 3.9ms | 7.8x | PASS |
| OSLBQP | 8 | 0 | Acceptable | Optimal | 7.24e-07 | 13 | 15 | 58us | 3.2ms | 55.8x | PASS |
| PALMER1 | 4 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 13 | 96us | 3.2ms | 32.9x | PASS |
| PALMER1A | 6 | 0 | Optimal | Optimal | 1.98e-14 | 47 | 48 | 326us | 10.1ms | 31.0x | PASS |
| PALMER1ANE | 6 | 35 | LocalInfeasi | IpoptStatus( | N/A | 45 | 0 | 508us | 565us | 1.1x | BOTH_FAIL |
| PALMER1B | 4 | 0 | Optimal | Optimal | 4.12e-15 | 18 | 17 | 128us | 3.6ms | 28.3x | PASS |
| PALMER1BNE | 4 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 202us | 558us | 2.8x | BOTH_FAIL |
| PALMER1C | 8 | 0 | Optimal | Optimal | 1.79e-13 | 4 | 1 | 63us | 898us | 14.2x | PASS |
| PALMER1D | 7 | 0 | Optimal | Optimal | 3.71e-14 | 2 | 1 | 49us | 942us | 19.1x | PASS |
| PALMER1E | 8 | 0 | Optimal | Optimal | 4.05e-13 | 119 | 55 | 1.0ms | 11.5ms | 11.3x | PASS |
| PALMER1ENE | 8 | 35 | LocalInfeasi | IpoptStatus( | N/A | 121 | 0 | 1.5ms | 566us | 0.4x | BOTH_FAIL |
| PALMER1NE | 4 | 31 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 168us | 542us | 3.2x | BOTH_FAIL |
| PALMER2 | 4 | 0 | Optimal | Optimal | 4.98e-16 | 17 | 28 | 112us | 6.9ms | 61.9x | PASS |
| PALMER2A | 6 | 0 | Optimal | Optimal | 7.77e-16 | 71 | 91 | 353us | 19.5ms | 55.3x | PASS |
| PALMER2ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 469us | 535us | 1.1x | BOTH_FAIL |
| PALMER2B | 4 | 0 | Optimal | Optimal | 5.44e-15 | 17 | 15 | 100us | 3.7ms | 36.8x | PASS |
| PALMER2BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 183us | 507us | 2.8x | BOTH_FAIL |
| PALMER2C | 8 | 0 | Optimal | Optimal | 5.44e-15 | 1 | 1 | 38us | 935us | 24.6x | PASS |
| PALMER2E | 8 | 0 | Optimal | Optimal | 2.09e-12 | 113 | 114 | 707us | 24.2ms | 34.2x | PASS |
| PALMER2ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 109 | 0 | 1.0ms | 584us | 0.6x | BOTH_FAIL |
| PALMER2NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 173us | 499us | 2.9x | BOTH_FAIL |
| PALMER3 | 4 | 0 | Acceptable | Optimal | 6.25e-02 | 20 | 44 | 140us | 8.2ms | 58.4x | MISMATCH |
| PALMER3A | 6 | 0 | Optimal | Optimal | 3.89e-16 | 66 | 73 | 345us | 15.0ms | 43.4x | PASS |
| PALMER3ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 74 | 0 | 593us | 535us | 0.9x | BOTH_FAIL |
| PALMER3B | 4 | 0 | Optimal | Optimal | 1.47e-15 | 15 | 15 | 106us | 3.7ms | 35.1x | PASS |
| PALMER3BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 124us | 524us | 4.2x | BOTH_FAIL |
| PALMER3C | 8 | 0 | Optimal | Optimal | 3.09e-15 | 1 | 1 | 38us | 952us | 25.2x | PASS |
| PALMER3E | 8 | 0 | Optimal | Optimal | 1.94e-13 | 28 | 32 | 208us | 6.3ms | 30.5x | PASS |
| PALMER3ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 339us | 559us | 1.6x | BOTH_FAIL |
| PALMER3NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 261us | 526us | 2.0x | BOTH_FAIL |
| PALMER4 | 4 | 0 | Optimal | Optimal | 5.72e-02 | 29 | 16 | 167us | 4.2ms | 25.0x | MISMATCH |
| PALMER4A | 6 | 0 | Optimal | Optimal | 1.79e-15 | 48 | 53 | 271us | 10.9ms | 40.3x | PASS |
| PALMER4ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 57 | 0 | 466us | 553us | 1.2x | BOTH_FAIL |
| PALMER4B | 4 | 0 | Optimal | Optimal | 1.82e-15 | 14 | 16 | 87us | 4.1ms | 47.4x | PASS |
| PALMER4BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 135us | 553us | 4.1x | BOTH_FAIL |
| PALMER4C | 8 | 0 | Optimal | Optimal | 3.29e-15 | 1 | 1 | 36us | 918us | 25.7x | PASS |
| PALMER4E | 8 | 0 | Optimal | Optimal | 8.84e-16 | 25 | 25 | 199us | 5.3ms | 26.6x | PASS |
| PALMER4ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 329us | 536us | 1.6x | BOTH_FAIL |
| PALMER4NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 213us | 553us | 2.6x | BOTH_FAIL |
| PALMER5A | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 13.1ms | 632.1ms | 48.3x | BOTH_FAIL |
| PALMER5ANE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 913.4ms | 606us | 0.0x | BOTH_FAIL |
| PALMER5B | 9 | 0 | Acceptable | Optimal | 1.57e-13 | 57 | 113 | 329us | 21.7ms | 65.8x | PASS |
| PALMER5BNE | 9 | 12 | LocalInfeasi | IpoptStatus( | N/A | 65 | 0 | 453us | 574us | 1.3x | BOTH_FAIL |
| PALMER5C | 6 | 0 | Optimal | Optimal | 2.50e-15 | 1 | 1 | 30us | 911us | 30.2x | PASS |
| PALMER5D | 4 | 0 | Optimal | Optimal | 3.25e-16 | 1 | 1 | 29us | 911us | 31.5x | PASS |
| PALMER5E | 8 | 0 | Acceptable | MaxIteration | N/A | 14 | 3000 | 100us | 477.6ms | 4770.5x | ipopt_FAIL |
| PALMER5ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 132us | 575us | 4.4x | BOTH_FAIL |
| PALMER6A | 6 | 0 | Optimal | Optimal | 2.42e-15 | 115 | 105 | 427us | 20.3ms | 47.5x | PASS |
| PALMER6ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 104 | 0 | 544us | 529us | 1.0x | BOTH_FAIL |
| PALMER6C | 8 | 0 | Optimal | Optimal | 2.30e-15 | 1 | 1 | 35us | 914us | 25.9x | PASS |
| PALMER6E | 8 | 0 | Optimal | Optimal | 2.54e-11 | 37 | 30 | 193us | 6.3ms | 32.6x | PASS |
| PALMER6ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 266us | 550us | 2.1x | BOTH_FAIL |
| PALMER7A | 6 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 9.9ms | 482.8ms | 48.8x | ipopt_FAIL |
| PALMER7ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 14.1ms | 594us | 0.0x | BOTH_FAIL |
| PALMER7C | 8 | 0 | Optimal | Optimal | 2.59e-13 | 2 | 1 | 39us | 896us | 22.9x | PASS |
| PALMER7E | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 13.4ms | 619.4ms | 46.3x | BOTH_FAIL |
| PALMER7ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 56.9ms | 612us | 0.0x | BOTH_FAIL |
| PALMER8A | 6 | 0 | Optimal | Optimal | 1.90e-15 | 28 | 36 | 148us | 7.9ms | 53.7x | PASS |
| PALMER8ANE | 6 | 12 | LocalInfeasi | IpoptStatus( | N/A | 36 | 0 | 237us | 582us | 2.5x | BOTH_FAIL |
| PALMER8C | 8 | 0 | Optimal | Optimal | 8.33e-15 | 1 | 1 | 33us | 877us | 26.2x | PASS |
| PALMER8E | 8 | 0 | Optimal | Optimal | 1.65e-17 | 30 | 23 | 160us | 4.9ms | 30.9x | PASS |
| PALMER8ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 208us | 502us | 2.4x | BOTH_FAIL |
| PARKCH | 15 | 0 | Acceptable | Optimal | 1.27e-14 | 19 | 17 | 4.38s | 3.74s | 0.9x | PASS |
| PENTAGON | 6 | 15 | Optimal | Optimal | 1.47e-05 | 13 | 19 | 212us | 4.6ms | 21.6x | PASS |
| PFIT1 | 3 | 3 | Optimal | Infeasible | N/A | 257 | 266 | 675us | 45.0ms | 66.6x | ipopt_FAIL |
| PFIT1LS | 3 | 0 | Optimal | Optimal | 1.64e-20 | 226 | 263 | 417us | 46.9ms | 112.3x | PASS |
| PFIT2 | 3 | 3 | Optimal | RestorationF | N/A | 111 | 247 | 318us | 47.0ms | 147.7x | ipopt_FAIL |
| PFIT2LS | 3 | 0 | Optimal | Optimal | 1.47e-20 | 106 | 82 | 214us | 15.1ms | 70.6x | PASS |
| PFIT3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 115 | 133 | 337us | 26.3ms | 78.3x | PASS |
| PFIT3LS | 3 | 0 | Optimal | Optimal | 2.90e-20 | 119 | 132 | 237us | 23.2ms | 98.1x | PASS |
| PFIT4 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 209 | 190 | 542us | 35.9ms | 66.3x | PASS |
| PFIT4LS | 3 | 0 | Optimal | Optimal | 6.57e-20 | 206 | 215 | 389us | 38.0ms | 97.6x | PASS |
| POLAK1 | 3 | 2 | Optimal | Optimal | 3.67e-09 | 15 | 5 | 140us | 1.8ms | 13.0x | PASS |
| POLAK2 | 11 | 2 | Optimal | Optimal | 3.59e-11 | 12 | 10 | 111us | 2.7ms | 23.9x | PASS |
| POLAK3 | 12 | 10 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 107.1ms | 656.6ms | 6.1x | ipopt_FAIL |
| POLAK4 | 3 | 3 | Acceptable | Optimal | 4.53e-09 | 14 | 4 | 84us | 1.6ms | 19.4x | PASS |
| POLAK5 | 3 | 2 | Optimal | Optimal | 2.00e-10 | 5 | 31 | 51us | 5.7ms | 111.5x | PASS |
| POLAK6 | 5 | 4 | Optimal | MaxIteration | N/A | 13 | 3000 | 95us | 835.0ms | 8751.0x | ipopt_FAIL |
| PORTFL1 | 12 | 1 | Acceptable | Optimal | 1.28e-06 | 10 | 9 | 284us | 2.5ms | 8.9x | PASS |
| PORTFL2 | 12 | 1 | Acceptable | Optimal | 4.37e-07 | 10 | 8 | 261us | 2.3ms | 8.9x | PASS |
| PORTFL3 | 12 | 1 | Acceptable | Optimal | 1.50e-06 | 10 | 9 | 258us | 2.5ms | 9.5x | PASS |
| PORTFL4 | 12 | 1 | Acceptable | Optimal | 3.83e-06 | 9 | 8 | 243us | 2.3ms | 9.6x | PASS |
| PORTFL6 | 12 | 1 | Acceptable | Optimal | 4.91e-07 | 10 | 8 | 261us | 2.2ms | 8.4x | PASS |
| POWELLBS | 2 | 2 | Optimal | Optimal | 0.00e+00 | 89 | 11 | 207us | 2.1ms | 9.9x | PASS |
| POWELLBSLS | 2 | 0 | Optimal | Optimal | 6.26e-26 | 90 | 91 | 147us | 13.5ms | 92.3x | PASS |
| POWELLSQ | 2 | 2 | Acceptable | Infeasible | N/A | 9 | 29 | 63us | 5.7ms | 89.7x | ipopt_FAIL |
| POWELLSQLS | 2 | 0 | Acceptable | Optimal | 6.89e-11 | 525 | 10 | 662us | 2.3ms | 3.4x | PASS |
| PRICE3NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 59us | 1.7ms | 28.2x | PASS |
| PRICE4 | 2 | 0 | Optimal | Optimal | 1.01e-22 | 13 | 8 | 56us | 1.8ms | 31.7x | PASS |
| PRICE4B | 2 | 0 | Optimal | Optimal | 3.11e-12 | 10 | 8 | 51us | 2.2ms | 43.6x | PASS |
| PRICE4NE | 2 | 2 | Optimal | Acceptable | 0.00e+00 | 10 | 23 | 61us | 3.9ms | 64.1x | PASS |
| PRODPL0 | 60 | 29 | Acceptable | Optimal | 1.12e-07 | 540 | 15 | 61.4ms | 4.2ms | 0.1x | PASS |
| PRODPL1 | 60 | 29 | Optimal | Optimal | 1.95e-01 | 28 | 28 | 3.0ms | 7.1ms | 2.4x | MISMATCH |
| PSPDOC | 4 | 0 | Optimal | Optimal | 3.50e-09 | 7 | 5 | 46us | 1.7ms | 38.2x | PASS |
| PT | 2 | 501 | Optimal | Optimal | 8.94e-03 | 4 | 106 | 1.8ms | 78.5ms | 42.6x | MISMATCH |
| QC | 9 | 4 | Acceptable | Optimal | 2.42e-02 | 12 | 44 | 146us | 9.4ms | 64.1x | MISMATCH |
| QCNEW | 9 | 3 | MaxIteration | Optimal | N/A | 2999 | 6 | 54.6ms | 1.9ms | 0.0x | ripopt_FAIL |
| QPCBLEND | 83 | 74 | Acceptable | Optimal | 1.83e-03 | 26 | 19 | 2.2ms | 5.6ms | 2.6x | MISMATCH |
| QPNBLEND | 83 | 74 | Acceptable | Optimal | 3.81e-05 | 31 | 18 | 2.6ms | 5.4ms | 2.1x | PASS |
| RAT42 | 3 | 9 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 136us | 519us | 3.8x | BOTH_FAIL |
| RAT42LS | 3 | 0 | Optimal | Optimal | 8.82e-16 | 22 | 28 | 101us | 4.6ms | 45.9x | PASS |
| RAT43 | 4 | 15 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 152us | 542us | 3.6x | BOTH_FAIL |
| RAT43LS | 4 | 0 | Acceptable | Optimal | 9.65e-01 | 2999 | 34 | 26.9ms | 5.8ms | 0.2x | MISMATCH |
| RECIPE | 3 | 3 | Acceptable | Optimal | 0.00e+00 | 19 | 16 | 80us | 3.1ms | 39.0x | PASS |
| RECIPELS | 3 | 0 | Acceptable | Optimal | 1.17e-10 | 19 | 29 | 59us | 5.3ms | 90.3x | PASS |
| RES | 20 | 14 | Acceptable | Optimal | 0.00e+00 | 11 | 10 | 222us | 2.4ms | 10.9x | PASS |
| RK23 | 17 | 11 | Acceptable | Optimal | 2.29e-06 | 27 | 10 | 450us | 2.9ms | 6.6x | PASS |
| ROBOT | 14 | 2 | Acceptable | IpoptStatus( | N/A | 9 | 18 | 171us | 4.6ms | 27.1x | ipopt_FAIL |
| ROSENBR | 2 | 0 | Optimal | Optimal | 0.00e+00 | 21 | 21 | 56us | 3.8ms | 67.7x | PASS |
| ROSENBRTU | 2 | 0 | Optimal | Optimal | 1.52e-24 | 45 | 87 | 93us | 13.9ms | 150.5x | PASS |
| ROSENMMX | 5 | 4 | Optimal | Optimal | 2.27e-10 | 148 | 13 | 699us | 3.4ms | 4.8x | PASS |
| ROSZMAN1 | 4 | 25 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 136us | 559us | 4.1x | BOTH_FAIL |
| ROSZMAN1LS | 4 | 0 | Optimal | Optimal | 7.59e-19 | 52 | 28 | 260us | 5.0ms | 19.4x | PASS |
| RSNBRNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 75us | 948us | 12.6x | PASS |
| S268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 40us | 3.1ms | 77.3x | MISMATCH |
| S308 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 40us | 2.0ms | 49.5x | PASS |
| S308NE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 67us | 508us | 7.6x | BOTH_FAIL |
| S316-322 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 10 | 7 | 61us | 1.9ms | 30.7x | PASS |
| S365 | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 28.2ms | 1.5ms | 0.1x | BOTH_FAIL |
| S365MOD | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 27.0ms | 1.4ms | 0.1x | BOTH_FAIL |
| SANTA | 21 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 680us | 546us | 0.8x | BOTH_FAIL |
| SANTALS | 21 | 0 | Optimal | Optimal | 2.36e-09 | 35 | 31 | 471us | 7.4ms | 15.8x | PASS |
| SIM2BQP | 2 | 0 | Acceptable | Optimal | 2.25e-06 | 8 | 5 | 46us | 1.7ms | 36.1x | PASS |
| SIMBQP | 2 | 0 | Optimal | Optimal | 8.29e-09 | 6 | 5 | 33us | 1.6ms | 48.5x | PASS |
| SIMPLLPA | 2 | 2 | Optimal | Optimal | 1.18e-08 | 6 | 8 | 42us | 2.3ms | 54.4x | PASS |
| SIMPLLPB | 2 | 3 | Optimal | Optimal | 1.90e-08 | 7 | 7 | 59us | 2.1ms | 34.7x | PASS |
| SINEVAL | 2 | 0 | Optimal | Optimal | 1.17e-40 | 42 | 42 | 82us | 6.8ms | 83.0x | PASS |
| SINVALNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 85us | 984us | 11.6x | PASS |
| SIPOW1 | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 81 | 30.2ms | 211.6ms | 7.0x | MISMATCH |
| SIPOW1M | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 88 | 33.0ms | 225.1ms | 6.8x | MISMATCH |
| SIPOW2 | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 69 | 29.7ms | 171.7ms | 5.8x | MISMATCH |
| SIPOW2M | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 73 | 33.3ms | 177.5ms | 5.3x | MISMATCH |
| SIPOW3 | 4 | 2000 | Optimal | Optimal | 7.10e-01 | 1 | 12 | 31.3ms | 35.4ms | 1.1x | MISMATCH |
| SIPOW4 | 4 | 2000 | Optimal | Optimal | 8.69e-01 | 1 | 11 | 44.2ms | 33.6ms | 0.8x | MISMATCH |
| SISSER | 2 | 0 | Acceptable | Optimal | 8.15e-11 | 15 | 18 | 54us | 2.8ms | 51.2x | PASS |
| SISSER2 | 2 | 0 | Acceptable | Optimal | 7.49e-11 | 16 | 20 | 56us | 3.3ms | 58.7x | PASS |
| SNAIL | 2 | 0 | Optimal | Optimal | 1.63e-26 | 63 | 63 | 121us | 10.1ms | 83.2x | PASS |
| SNAKE | 2 | 2 | Optimal | Optimal | 1.38e-02 | 717 | 8 | 16.9ms | 2.4ms | 0.1x | MISMATCH |
| SPANHYD | 97 | 33 | Acceptable | Optimal | 1.38e-10 | 61 | 20 | 7.0ms | 6.0ms | 0.9x | PASS |
| SPIRAL | 3 | 2 | Acceptable | Infeasible | N/A | 124 | 370 | 319us | 59.9ms | 187.4x | ipopt_FAIL |
| SSI | 3 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 3.6ms | 435.1ms | 122.5x | ipopt_FAIL |
| SSINE | 3 | 2 | MaxIteration | Optimal | N/A | 2999 | 224 | 16.5ms | 33.6ms | 2.0x | ripopt_FAIL |
| STANCMIN | 3 | 2 | Optimal | Optimal | 9.32e-09 | 10 | 9 | 61us | 2.4ms | 39.5x | PASS |
| STRATEC | 10 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| STREG | 4 | 0 | Optimal | Optimal | 8.90e-02 | 21 | 13 | 62us | 2.7ms | 44.3x | MISMATCH |
| STREGNE | 4 | 2 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 35us | 1.1ms | 31.8x | PASS |
| SUPERSIM | 2 | 2 | Optimal | Optimal | 2.22e-16 | 7 | 1 | 53us | 1.0ms | 19.7x | PASS |
| SWOPF | 83 | 92 | Acceptable | Optimal | 5.18e-07 | 14 | 13 | 1.3ms | 4.2ms | 3.3x | PASS |
| SYNTHES1 | 6 | 6 | Acceptable | Optimal | 3.50e-07 | 13 | 8 | 89us | 2.2ms | 24.7x | PASS |
| SYNTHES2 | 11 | 14 | Acceptable | Optimal | 9.03e-07 | 22 | 14 | 257us | 3.3ms | 13.0x | PASS |
| SYNTHES3 | 17 | 23 | Acceptable | Optimal | 7.93e-01 | 95 | 13 | 2.0ms | 3.3ms | 1.7x | MISMATCH |
| TAME | 2 | 1 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 46us | 1.6ms | 35.5x | PASS |
| TAX13322 | 72 | 1261 | Acceptable | MaxIteration | N/A | 243 | 3000 | 415.7ms | 19.23s | 46.2x | ipopt_FAIL |
| TAXR13322 | 72 | 1261 | Optimal | Acceptable | 9.56e-01 | 138 | 56 | 526.2ms | 2.66s | 5.1x | MISMATCH |
| TENBARS1 | 18 | 9 | Acceptable | Optimal | 3.12e-03 | 36 | 39 | 424us | 7.8ms | 18.3x | MISMATCH |
| TENBARS2 | 18 | 8 | Acceptable | Optimal | 6.50e-08 | 48 | 33 | 576us | 6.8ms | 11.8x | PASS |
| TENBARS3 | 18 | 8 | Optimal | Optimal | 1.01e-08 | 19 | 34 | 236us | 7.0ms | 29.6x | PASS |
| TENBARS4 | 18 | 9 | Optimal | Optimal | 1.79e-10 | 40 | 14 | 648us | 3.7ms | 5.7x | PASS |
| TFI1 | 3 | 101 | Optimal | Optimal | 9.76e-01 | 72 | 19 | 15.3ms | 7.1ms | 0.5x | MISMATCH |
| TFI2 | 3 | 101 | Optimal | Optimal | 9.85e-03 | 255 | 8 | 11.5ms | 3.2ms | 0.3x | MISMATCH |
| TFI3 | 3 | 101 | Optimal | Optimal | 2.81e-02 | 33 | 13 | 1.7ms | 4.7ms | 2.7x | MISMATCH |
| THURBER | 7 | 37 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 301us | 566us | 1.9x | BOTH_FAIL |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 1.02e-14 | 2999 | 19 | 36.9ms | 3.9ms | 0.1x | PASS |
| TOINTGOR | 50 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 253us | 1.8ms | 7.3x | PASS |
| TOINTPSP | 50 | 0 | Optimal | Optimal | 0.00e+00 | 19 | 20 | 572us | 5.2ms | 9.2x | PASS |
| TOINTQOR | 50 | 0 | Optimal | Optimal | 1.93e-16 | 1 | 1 | 72us | 965us | 13.4x | PASS |
| TRIGGER | 7 | 6 | Acceptable | Optimal | 0.00e+00 | 15 | 15 | 105us | 3.1ms | 29.1x | PASS |
| TRO3X3 | 30 | 13 | Acceptable | Optimal | 3.58e-03 | 2999 | 47 | 59.9ms | 10.4ms | 0.2x | MISMATCH |
| TRO4X4 | 63 | 25 | Acceptable | IpoptStatus( | N/A | 108 | 157 | 34.9ms | 42.8ms | 1.2x | ipopt_FAIL |
| TRO6X2 | 45 | 21 | Acceptable | RestorationF | N/A | 2999 | 353 | 146.7ms | 88.7ms | 0.6x | ipopt_FAIL |
| TRUSPYR1 | 11 | 4 | Optimal | Optimal | 1.25e-08 | 35 | 10 | 280us | 2.6ms | 9.2x | PASS |
| TRUSPYR2 | 11 | 11 | Optimal | Optimal | 1.46e-08 | 123 | 13 | 3.9ms | 3.3ms | 0.9x | PASS |
| TRY-B | 2 | 1 | Optimal | Optimal | 6.52e-17 | 12 | 23 | 65us | 4.9ms | 74.3x | PASS |
| TWOBARS | 2 | 2 | Optimal | Optimal | 9.02e-09 | 19 | 8 | 84us | 2.2ms | 26.4x | PASS |
| VESUVIA | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 462 | 0 | 161.5ms | 863us | 0.0x | BOTH_FAIL |
| VESUVIALS | 8 | 0 | Acceptable | Optimal | 3.39e-01 | 330 | 48 | 71.2ms | 18.2ms | 0.3x | MISMATCH |
| VESUVIO | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 1.60s | 805us | 0.0x | BOTH_FAIL |
| VESUVIOLS | 8 | 0 | Acceptable | Optimal | 2.41e-15 | 2999 | 10 | 1.08s | 5.0ms | 0.0x | PASS |
| VESUVIOU | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 8.9ms | 768us | 0.1x | BOTH_FAIL |
| VESUVIOULS | 8 | 0 | Acceptable | Optimal | 5.55e-16 | 22 | 8 | 5.5ms | 3.8ms | 0.7x | PASS |
| VIBRBEAM | 8 | 0 | Optimal | Optimal | 9.48e-01 | 115 | 58 | 2.4ms | 10.1ms | 4.2x | MISMATCH |
| VIBRBEAMNE | 8 | 30 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 1.3ms | 554us | 0.4x | BOTH_FAIL |
| WACHBIEG | 3 | 2 | Acceptable | Infeasible | N/A | 27 | 15 | 233us | 4.0ms | 17.0x | ipopt_FAIL |
| WATER | 31 | 10 | Acceptable | Optimal | 7.93e-09 | 94 | 17 | 2.7ms | 3.8ms | 1.4x | PASS |
| WAYSEA1 | 2 | 0 | Optimal | Optimal | 2.69e-15 | 15 | 14 | 56us | 2.4ms | 42.4x | PASS |
| WAYSEA1B | 2 | 0 | Optimal | Optimal | 7.97e-09 | 13 | 14 | 43us | 2.9ms | 68.7x | PASS |
| WAYSEA1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 15 | 7 | 68us | 1.7ms | 24.6x | PASS |
| WAYSEA2 | 2 | 0 | Optimal | Optimal | 9.85e-18 | 23 | 22 | 68us | 3.3ms | 48.5x | PASS |
| WAYSEA2B | 2 | 0 | Optimal | Optimal | 1.61e-11 | 21 | 22 | 62us | 4.2ms | 67.0x | PASS |
| WAYSEA2NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 22 | 11 | 80us | 2.2ms | 27.2x | PASS |
| WEEDS | 3 | 0 | Optimal | Optimal | 3.43e-16 | 23 | 28 | 108us | 6.4ms | 59.2x | PASS |
| WEEDSNE | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 190us | 542us | 2.8x | BOTH_FAIL |
| WOMFLET | 3 | 3 | Acceptable | Optimal | 1.00e+00 | 14 | 8 | 82us | 2.3ms | 27.4x | MISMATCH |
| YFIT | 3 | 0 | Optimal | Optimal | 1.37e-19 | 36 | 36 | 149us | 7.3ms | 49.0x | PASS |
| YFITNE | 3 | 17 | Acceptable | IpoptStatus( | N/A | 36 | 0 | 233us | 532us | 2.3x | ipopt_FAIL |
| YFITU | 3 | 0 | Optimal | Optimal | 6.69e-21 | 36 | 36 | 140us | 5.7ms | 40.5x | PASS |
| ZANGWIL2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 30us | 954us | 31.5x | PASS |
| ZANGWIL3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 32us | 969us | 30.4x | PASS |
| ZECEVIC2 | 2 | 2 | Optimal | Optimal | 6.97e-10 | 11 | 8 | 64us | 2.1ms | 32.1x | PASS |
| ZECEVIC3 | 2 | 2 | Optimal | Optimal | 6.87e-10 | 12 | 17 | 82us | 3.6ms | 43.9x | PASS |
| ZECEVIC4 | 2 | 2 | Optimal | Optimal | 2.68e-09 | 10 | 10 | 61us | 2.6ms | 42.9x | PASS |
| ZY2 | 3 | 2 | Acceptable | Optimal | 5.33e-06 | 14 | 14 | 82us | 3.6ms | 43.6x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 132.6 | 44.4 |
| Median | 15 | 13 |
| Total  | 72137 | 24153 |

- ripopt fewer iterations: 166/544
- Ipopt fewer iterations: 291/544
- Tied: 87/544

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 85.3ms | 35.3ms |
| Median | 108us | 3.1ms |
| Total  | 46.42s | 19.20s |

- Geometric mean speedup (Ipopt_time/ripopt_time): **14.56x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 496/544 problems
- Ipopt faster: 48/544 problems
- Overall speedup (total time): 0.41x

## Failure Analysis

### Problems where only ripopt fails (14)

| Problem | n | m | ripopt status | Ipopt obj |
|---------|---|---|---------------|-----------|
| DENSCHNDNE | 3 | 3 | LocalInfeasibility | 0.000000e+00 |
| EXPFITC | 5 | 502 | MaxIterations | 2.330262e-02 |
| HAHN1LS | 7 | 0 | MaxIterations | 3.338424e+01 |
| HATFLDFLNE | 3 | 3 | LocalInfeasibility | 0.000000e+00 |
| HATFLDH | 4 | 7 | MaxIterations | -2.450000e+01 |
| HS83 | 5 | 3 | MaxIterations | -3.066554e+04 |
| HS84 | 5 | 3 | MaxIterations | -5.280335e+06 |
| HYDCAR20 | 99 | 99 | LocalInfeasibility | 0.000000e+00 |
| MGH10SLS | 3 | 0 | MaxIterations | 8.794586e+01 |
| MINMAXBD | 5 | 20 | MaxIterations | 1.157064e+02 |
| OET7 | 7 | 1002 | MaxIterations | 4.465915e-05 |
| OSBORNEA | 5 | 0 | MaxIterations | 5.464895e-05 |
| QCNEW | 9 | 3 | MaxIterations | -8.065219e+02 |
| SSINE | 3 | 2 | MaxIterations | 0.000000e+00 |

### Problems where only Ipopt fails (39)

| Problem | n | m | Ipopt status | ripopt obj |
|---------|---|---|--------------|------------|
| ARGAUSS | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| AVION2 | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CRESC100 | 6 | 200 | Infeasible | 7.445687e-01 |
| DECONVB | 63 | 0 | MaxIterations | 3.233019e-03 |
| DENSCHNBNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DENSCHNENE | 3 | 3 | Infeasible | 0.000000e+00 |
| DEVGLA1NE | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | 9 | 3 | ErrorInStepComputation | -8.274326e+02 |
| EXP2NE | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| FBRAIN3 | 6 | 2211 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| LANCZOS1 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS2 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS3 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 |
| MESH | 41 | 48 | IpoptStatus(4) | -3.499682e+06 |
| NYSTROM5 | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5E | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | 6 | 0 | MaxIterations | 1.033491e+01 |
| PFIT1 | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | 12 | 10 | MaxIterations | 5.020111e+01 |
| POLAK6 | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | 3 | 2 | Infeasible | 2.563322e-12 |
| SSI | 3 | 0 | MaxIterations | 1.376194e-09 |
| TAX13322 | 72 | 1261 | MaxIterations | -2.688286e+04 |
| TRO4X4 | 63 | 25 | IpoptStatus(4) | 8.999999e+00 |
| TRO6X2 | 45 | 21 | RestorationFailed | 1.224983e+03 |
| WACHBIEG | 3 | 2 | Infeasible | 1.000006e+00 |
| YFITNE | 3 | 17 | IpoptStatus(-10) | 0.000000e+00 |

### Problems where both fail (130)

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
| FBRAIN3LS | 6 | 0 | MaxIterations | MaxIterations |
| FBRAINNE | 2 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS1 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS2 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS3 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GBRAIN | 2 | 2200 | LocalInfeasibility | IpoptStatus(-10) |
| GROWTH | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| HAHN1 | 7 | 236 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDBNE | 4 | 4 | LocalInfeasibility | Infeasible |
| HATFLDDNE | 3 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDENE | 3 | 21 | LocalInfeasibility | IpoptStatus(-10) |
| HIMMELBD | 2 | 2 | LocalInfeasibility | Infeasible |
| HIMMELBFNE | 4 | 7 | LocalInfeasibility | IpoptStatus(-10) |
| HS25NE | 3 | 99 | LocalInfeasibility | IpoptStatus(-10) |
| HS2NE | 2 | 2 | LocalInfeasibility | Infeasible |
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
| OET5 | 5 | 1002 | Timeout | Timeout |
| OET6 | 5 | 1002 | Timeout | Timeout |
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

### Objective mismatches (92)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 51
- **Convergence gap** (one Acceptable): 41
- **Better objective found by**: ripopt 17, Ipopt 75

| Problem | ripopt obj | Ipopt obj | Rel Diff | r_status | i_status | Better |
|---------|-----------|-----------|----------|----------|----------|--------|
| MAKELA1 | 1.414214e+00 | -1.414214e+00 | 2.00e+00 | Optimal | Optimal | ipopt |
| SIPOW1 | 4.338956e-01 | -1.000000e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW1M | 4.338960e-01 | -1.000001e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| HS59 | 1.982326e+01 | -6.749505e+00 | 1.34e+00 | Optimal | Optimal | ipopt |
| SIPOW2M | 3.310254e-01 | -1.000005e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| SIPOW2 | 3.310227e-01 | -1.000000e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| GOFFIN | 1.918594e+01 | -9.500002e-09 | 1.00e+00 | Acceptable | Optimal | ipopt |
| DEVGLA1B | 1.057042e+05 | 8.216237e-18 | 1.00e+00 | Acceptable | Optimal | ipopt |
| DEVGLA2 | 5.932686e+01 | 6.672171e-19 | 1.00e+00 | Optimal | Optimal | ipopt |
| LEVYMONT5 | 1.248706e+01 | 1.239502e-25 | 1.00e+00 | Acceptable | Optimal | ipopt |
| WOMFLET | 9.966429e-12 | 6.050000e+00 | 1.00e+00 | Acceptable | Optimal | ripopt |
| GIGOMEZ1 | -5.836123e-10 | -3.000000e+00 | 1.00e+00 | Optimal | Optimal | ipopt |
| HS268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| S268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| DANWOODLS | 1.039178e+02 | 4.317308e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| HS16 | 2.316106e+01 | 2.500000e-01 | 9.89e-01 | Acceptable | Optimal | ipopt |
| ELATTAR | 1.054115e+00 | 7.420618e+01 | 9.86e-01 | Optimal | Optimal | ripopt |
| KSIP | 2.768715e+01 | 5.757979e-01 | 9.79e-01 | Optimal | Optimal | ipopt |
| TFI1 | 2.269758e+02 | 5.334687e+00 | 9.76e-01 | Optimal | Optimal | ipopt |
| RAT43LS | 2.525083e+05 | 8.786405e+03 | 9.65e-01 | Acceptable | Optimal | ipopt |
| TAXR13322 | -7.850173e+03 | -3.429089e+02 | 9.56e-01 | Optimal | Acceptable | ripopt |
| MWRIGHT | 1.288383e+00 | 2.497881e+01 | 9.48e-01 | Optimal | Optimal | ripopt |
| VIBRBEAM | 6.346254e+00 | 3.322376e-01 | 9.48e-01 | Optimal | Optimal | ipopt |
| BT4 | -4.551055e+01 | -3.704768e+00 | 9.19e-01 | Optimal | Optimal | ripopt |
| BIGGSC4 | -3.125001e+00 | -2.450000e+01 | 8.72e-01 | Acceptable | Optimal | ipopt |
| SIPOW4 | 2.080345e+00 | 2.723620e-01 | 8.69e-01 | Optimal | Optimal | ipopt |
| HIMMELP2 | -6.205394e+01 | -8.198044e+00 | 8.68e-01 | Optimal | Optimal | ripopt |
| HIMMELP3 | -7.913699e+00 | -5.901318e+01 | 8.66e-01 | Optimal | Optimal | ipopt |
| SYNTHES3 | 7.275706e+01 | 1.508219e+01 | 7.93e-01 | Acceptable | Optimal | ipopt |
| CAMEL6 | -2.154638e-01 | -1.031628e+00 | 7.91e-01 | Optimal | Optimal | ipopt |
| HS23 | 9.472223e+00 | 2.000000e+00 | 7.89e-01 | Optimal | Optimal | ipopt |
| HS54 | -1.566691e-01 | -9.080748e-01 | 7.51e-01 | Optimal | Optimal | ipopt |
| SIPOW3 | 1.846738e+00 | 5.346586e-01 | 7.10e-01 | Optimal | Optimal | ipopt |
| MAKELA4 | 6.877365e-01 | -9.600000e-09 | 6.88e-01 | Optimal | Optimal | ipopt |
| HIMMELP5 | -2.131055e+01 | -5.901318e+01 | 6.39e-01 | Acceptable | Optimal | ipopt |
| HALDMADS | 8.180683e-01 | 2.218282e+00 | 6.31e-01 | Optimal | Optimal | ripopt |
| MISTAKE | -5.000000e-01 | -1.000000e+00 | 5.00e-01 | Acceptable | Optimal | ipopt |
| EGGCRATE | 1.897639e+01 | 9.488197e+00 | 5.00e-01 | Optimal | Optimal | ipopt |
| ECKERLE4LS | 4.988568e-01 | 1.463589e-03 | 4.97e-01 | Acceptable | Optimal | ipopt |
| HS81 | 4.388512e-01 | 5.394985e-02 | 3.85e-01 | Optimal | Optimal | ipopt |
| VESUVIALS | 1.500440e+03 | 9.914100e+02 | 3.39e-01 | Acceptable | Optimal | ipopt |
| HYDC20LS | 2.976453e-01 | 2.967522e-15 | 2.98e-01 | Acceptable | Optimal | ipopt |
| HS98 | 4.071290e+00 | 3.135806e+00 | 2.30e-01 | Acceptable | Optimal | ipopt |
| HS97 | 4.071255e+00 | 3.135806e+00 | 2.30e-01 | Acceptable | Optimal | ipopt |
| HS67 | -8.955413e+02 | -1.162119e+03 | 2.29e-01 | Optimal | Optimal | ipopt |
| EG1 | -1.132801e+00 | -1.429307e+00 | 2.07e-01 | Optimal | Optimal | ipopt |
| PRODPL1 | 4.439207e+01 | 3.573897e+01 | 1.95e-01 | Optimal | Optimal | ipopt |
| HS70 | 1.798079e-01 | 7.498464e-03 | 1.72e-01 | Optimal | Optimal | ipopt |
| DEGENLPA | 3.557317e+00 | 3.054881e+00 | 1.41e-01 | Acceptable | Optimal | ipopt |
| MUONSINELS | 5.050172e+04 | 4.387412e+04 | 1.31e-01 | Acceptable | Optimal | ipopt |
| MAKELA2 | 8.244898e+00 | 7.200000e+00 | 1.27e-01 | Optimal | Optimal | ipopt |
| MSS1 | -1.600000e+01 | -1.400000e+01 | 1.25e-01 | Acceptable | Optimal | ripopt |
| HS85 | -2.004792e+00 | -2.215605e+00 | 9.51e-02 | Acceptable | Optimal | ipopt |
| STREG | 3.743976e-21 | 8.901950e-02 | 8.90e-02 | Optimal | Optimal | ripopt |
| OET1 | 6.222150e-01 | 5.382431e-01 | 8.40e-02 | Optimal | Optimal | ipopt |
| PALMER3 | 2.265958e+03 | 2.416980e+03 | 6.25e-02 | Acceptable | Optimal | ripopt |
| PALMER4 | 2.424016e+03 | 2.285383e+03 | 5.72e-02 | Optimal | Optimal | ipopt |
| DEMBO7 | 1.848169e+02 | 1.747870e+02 | 5.43e-02 | Acceptable | Optimal | ipopt |
| TFI3 | 4.425675e+00 | 4.301158e+00 | 2.81e-02 | Optimal | Optimal | ipopt |
| MGH17SLS | 5.464895e-05 | 2.451788e-02 | 2.45e-02 | Optimal | Optimal | ripopt |
| QC | -9.333976e+02 | -9.565379e+02 | 2.42e-02 | Acceptable | Optimal | ipopt |
| CB2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CHACONN1 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| GIGOMEZ2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| HS55 | 6.666671e+00 | 6.805833e+00 | 2.04e-02 | Acceptable | Optimal | ripopt |
| LRCOVTYPE | 5.901541e-01 | 5.723072e-01 | 1.78e-02 | Optimal | Optimal | ipopt |
| HS118 | 6.759936e+02 | 6.648204e+02 | 1.65e-02 | Acceptable | Optimal | ipopt |
| ACOPR30 | 5.850632e+02 | 5.768924e+02 | 1.40e-02 | Optimal | Optimal | ipopt |
| SNAKE | 1.355757e-02 | -1.999999e-04 | 1.38e-02 | Optimal | Optimal | ipopt |
| MAXLIKA | 1.149346e+03 | 1.136307e+03 | 1.13e-02 | Acceptable | Optimal | ipopt |
| TFI2 | 6.588814e-01 | 6.490311e-01 | 9.85e-03 | Optimal | Optimal | ipopt |
| PT | 1.873320e-01 | 1.783942e-01 | 8.94e-03 | Optimal | Optimal | ipopt |
| CLIFF | 1.997866e-01 | 2.072380e-01 | 7.45e-03 | Optimal | Optimal | ripopt |
| HS43 | -4.371038e+01 | -4.400000e+01 | 6.58e-03 | Optimal | Optimal | ipopt |
| DGOSPEC | -9.887540e+02 | -9.933506e+02 | 4.63e-03 | Acceptable | Optimal | ipopt |
| ANTWERP | 3.258300e+03 | 3.245241e+03 | 4.01e-03 | Acceptable | Optimal | ipopt |
| TRO3X3 | 8.999741e+00 | 8.967478e+00 | 3.58e-03 | Acceptable | Optimal | ipopt |
| TENBARS1 | 2.302549e+03 | 2.295373e+03 | 3.12e-03 | Acceptable | Optimal | ipopt |
| LIN | -1.960628e-02 | -1.757754e-02 | 2.03e-03 | Optimal | Optimal | ripopt |
| QPCBLEND | -6.015132e-03 | -7.842801e-03 | 1.83e-03 | Acceptable | Optimal | ipopt |
| EXPFITA | 2.536784e-03 | 1.136646e-03 | 1.40e-03 | Optimal | Optimal | ipopt |
| DECONVC | 3.797755e-03 | 2.569475e-03 | 1.23e-03 | Acceptable | Optimal | ipopt |
| ACOPR14 | 8.091232e+03 | 8.081526e+03 | 1.20e-03 | Acceptable | Optimal | ipopt |
| DEGENLPB | -3.073079e+01 | -3.076401e+01 | 1.08e-03 | Acceptable | Optimal | ipopt |
| LOOTSMA | 1.415263e+00 | 1.414213e+00 | 7.42e-04 | Optimal | Optimal | ipopt |
| HS116 | 9.751941e+01 | 9.758747e+01 | 6.97e-04 | Acceptable | Optimal | ripopt |
| HS96 | 1.624621e-02 | 1.561775e-02 | 6.28e-04 | Acceptable | Optimal | ipopt |
| HS95 | 1.620728e-02 | 1.561772e-02 | 5.90e-04 | Acceptable | Optimal | ipopt |
| HS45 | 1.000314e+00 | 1.000000e+00 | 3.14e-04 | Acceptable | Optimal | ipopt |
| LSC2LS | 1.333749e+01 | 1.333439e+01 | 2.32e-04 | Acceptable | Optimal | ipopt |
| DENSCHND | 6.940852e-09 | 2.221899e-04 | 2.22e-04 | Acceptable | Optimal | ripopt |
| HS13 | 9.943965e-01 | 9.945785e-01 | 1.82e-04 | Acceptable | Optimal | ripopt |

---
*Generated by cutest_suite/compare.py*