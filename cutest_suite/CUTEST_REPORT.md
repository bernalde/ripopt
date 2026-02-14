# CUTEst Benchmark Report

Comparison of ripopt vs Ipopt (C++) on the CUTEst test set.

## Executive Summary

- **Total problems**: 727
- **ripopt solved**: 571/727 (78.5%)
- **Ipopt solved**: 556/727 (76.5%)
- **Both solved**: 533/727
- **Matching solutions** (rel obj diff < 1e-4): 438/533

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (438 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 1.37e-06 |
| Median | 1.57e-13 |
| Max    | 7.63e-05 |

**All both-solved** (533 problems, including 95 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 8.57e-02 |
| Median | 5.11e-10 |
| Max    | 2.00e+00 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 493 | 359 | 340 | 326 | 258 |
| unconstrained | 234 | 212 | 216 | 207 | 180 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| 3PK | 30 | 0 | Optimal | Optimal | 1.42e-15 | 9 | 9 | 109us | 2.2ms | 20.4x | PASS |
| ACOPP14 | 38 | 68 | Optimal | Optimal | 9.79e-10 | 16 | 9 | 1.5ms | 3.6ms | 2.4x | PASS |
| ACOPP30 | 72 | 142 | Optimal | Optimal | 6.40e-09 | 47 | 13 | 10.2ms | 6.2ms | 0.6x | PASS |
| ACOPR14 | 38 | 82 | Acceptable | Optimal | 3.22e-03 | 2999 | 13 | 582.7ms | 4.9ms | 0.0x | MISMATCH |
| ACOPR30 | 72 | 172 | Optimal | Optimal | 6.70e-03 | 428 | 221 | 3.66s | 115.0ms | 0.0x | MISMATCH |
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 9 | 3 | 1.5ms | 776us | 0.5x | PASS |
| AIRCRFTB | 8 | 0 | MaxIteration | Optimal | N/A | 2999 | 15 | 19.3ms | 2.7ms | 0.1x | ripopt_FAIL |
| AIRPORT | 84 | 42 | Optimal | Optimal | 4.99e-09 | 12 | 13 | 3.4ms | 5.8ms | 1.7x | PASS |
| AKIVA | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 58us | 1.0ms | 17.3x | PASS |
| ALLINIT | 4 | 0 | MaxIteration | Optimal | N/A | 2999 | 20 | 18.2ms | 3.4ms | 0.2x | ripopt_FAIL |
| ALLINITA | 4 | 4 | MaxIteration | Optimal | N/A | 2999 | 12 | 22.9ms | 2.3ms | 0.1x | ripopt_FAIL |
| ALLINITC | 4 | 1 | Acceptable | Optimal | 8.85e-06 | 25 | 17 | 235us | 3.0ms | 12.8x | PASS |
| ALLINITU | 4 | 0 | Optimal | Optimal | 3.09e-16 | 8 | 14 | 21us | 2.3ms | 109.2x | PASS |
| ALSOTAME | 2 | 1 | Optimal | Optimal | 1.47e-08 | 10 | 8 | 24us | 1.7ms | 71.7x | PASS |
| ANTWERP | 27 | 10 | Optimal | Optimal | 7.48e-08 | 47 | 108 | 1.7ms | 22.9ms | 13.8x | PASS |
| ARGAUSS | 3 | 15 | Acceptable | IpoptStatus( | N/A | 2 | 0 | 15us | 300us | 20.6x | ipopt_FAIL |
| AVGASA | 8 | 10 | Acceptable | Optimal | 2.70e-01 | 11 | 9 | 59us | 2.0ms | 33.6x | MISMATCH |
| AVGASB | 8 | 10 | Acceptable | Optimal | 1.71e-01 | 19 | 11 | 101us | 2.5ms | 24.7x | MISMATCH |
| AVION2 | 49 | 15 | Acceptable | MaxIteration | N/A | 26 | 3000 | 3.2ms | 660.9ms | 208.9x | ipopt_FAIL |
| BA-L1 | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 6 | 546us | 1.8ms | 3.3x | PASS |
| BA-L1LS | 57 | 0 | Optimal | Optimal | 7.65e-21 | 7 | 10 | 615us | 2.5ms | 4.0x | PASS |
| BA-L1SP | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.4ms | 2.6ms | 2.0x | PASS |
| BA-L1SPLS | 57 | 0 | Optimal | Optimal | 6.48e-17 | 23 | 9 | 7.0ms | 4.2ms | 0.6x | PASS |
| BARD | 3 | 0 | Optimal | Optimal | 1.73e-18 | 8 | 8 | 23us | 1.4ms | 59.9x | PASS |
| BARDNE | 3 | 15 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 24us | 279us | 11.6x | BOTH_FAIL |
| BATCH | 48 | 73 | Acceptable | Optimal | 4.26e-05 | 2999 | 29 | 130.5ms | 7.2ms | 0.1x | PASS |
| BEALE | 2 | 0 | Optimal | Optimal | 4.34e-18 | 7 | 8 | 15us | 1.7ms | 109.6x | PASS |
| BEALENE | 2 | 3 | Optimal | IpoptStatus( | N/A | 6 | 0 | 17us | 282us | 16.8x | ipopt_FAIL |
| BENNETT5 | 3 | 154 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 297us | 303us | 1.0x | BOTH_FAIL |
| BENNETT5LS | 3 | 0 | Optimal | Optimal | 1.00e+00 | 14 | 21 | 511us | 3.9ms | 7.5x | MISMATCH |
| BIGGS3 | 6 | 0 | MaxIteration | Optimal | N/A | 2999 | 9 | 32.4ms | 2.0ms | 0.1x | ripopt_FAIL |
| BIGGS5 | 6 | 0 | MaxIteration | Optimal | N/A | 2999 | 20 | 30.3ms | 3.8ms | 0.1x | ripopt_FAIL |
| BIGGS6 | 6 | 0 | Optimal | Optimal | 4.32e-20 | 89 | 79 | 358us | 12.1ms | 33.7x | PASS |
| BIGGS6NE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 391us | 271us | 0.7x | BOTH_FAIL |
| BIGGSC4 | 4 | 7 | Acceptable | Optimal | 8.72e-01 | 2999 | 17 | 9.4ms | 3.4ms | 0.4x | MISMATCH |
| BLEACHNG | 17 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| BOOTH | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 6us | 475us | 86.3x | PASS |
| BOX2 | 3 | 0 | MaxIteration | Optimal | N/A | 2999 | 8 | 23.1ms | 1.4ms | 0.1x | ripopt_FAIL |
| BOX3 | 3 | 0 | Optimal | Optimal | 1.69e-25 | 11 | 9 | 32us | 1.6ms | 52.3x | PASS |
| BOX3NE | 3 | 10 | Optimal | IpoptStatus( | N/A | 5 | 0 | 24us | 303us | 12.8x | ipopt_FAIL |
| BOXBOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 55us | 307us | 5.5x | BOTH_FAIL |
| BOXBODLS | 2 | 0 | Optimal | Optimal | 4.40e-13 | 3 | 13 | 10us | 2.7ms | 265.0x | PASS |
| BQP1VAR | 1 | 0 | Optimal | Optimal | 9.98e-09 | 7 | 5 | 9us | 1.2ms | 124.7x | PASS |
| BQPGABIM | 50 | 0 | Optimal | Optimal | 1.46e-05 | 42 | 12 | 1.5ms | 2.7ms | 1.8x | PASS |
| BQPGASIM | 50 | 0 | Acceptable | Optimal | 6.13e-06 | 11 | 12 | 336us | 2.8ms | 8.3x | PASS |
| BRANIN | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 18us | 1.6ms | 88.7x | PASS |
| BRKMCC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 8us | 819us | 105.1x | PASS |
| BROWNBS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 17us | 1.1ms | 64.8x | PASS |
| BROWNBSNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 13 | 0 | 29us | 268us | 9.4x | ipopt_FAIL |
| BROWNDEN | 4 | 0 | Optimal | Optimal | 1.70e-16 | 8 | 8 | 32us | 1.2ms | 37.5x | PASS |
| BROWNDENE | 4 | 20 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 25.4ms | 312us | 0.0x | BOTH_FAIL |
| BT1 | 2 | 1 | Optimal | Optimal | 2.40e-09 | 17 | 7 | 44us | 1.7ms | 39.0x | PASS |
| BT10 | 2 | 2 | Optimal | Optimal | 2.79e-09 | 7 | 6 | 20us | 1.1ms | 52.8x | PASS |
| BT11 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 8 | 24us | 1.3ms | 53.3x | PASS |
| BT12 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 18us | 814us | 44.9x | PASS |
| BT13 | 5 | 1 | Acceptable | Optimal | 1.00e-08 | 27 | 24 | 61us | 4.2ms | 68.8x | PASS |
| BT2 | 3 | 1 | Optimal | Optimal | 0.00e+00 | 12 | 12 | 25us | 1.8ms | 71.7x | PASS |
| BT3 | 5 | 3 | Optimal | Optimal | 1.22e-14 | 1 | 1 | 10us | 491us | 48.7x | PASS |
| BT4 | 3 | 2 | Optimal | Optimal | 9.19e-01 | 6 | 9 | 18us | 1.6ms | 89.7x | MISMATCH |
| BT5 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 20us | 1.2ms | 61.4x | PASS |
| BT6 | 5 | 2 | Optimal | Optimal | 3.46e-12 | 9 | 13 | 35us | 2.3ms | 65.9x | PASS |
| BT7 | 5 | 3 | Optimal | Optimal | 1.50e-01 | 23 | 16 | 90us | 3.0ms | 33.2x | MISMATCH |
| BT8 | 5 | 2 | Acceptable | Optimal | 3.73e-09 | 32 | 14 | 105us | 2.1ms | 20.4x | PASS |
| BT9 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 38us | 1.9ms | 49.6x | PASS |
| BURKEHAN | 1 | 1 | RestorationF | Infeasible | N/A | 229 | 11 | 3.8ms | 2.8ms | 0.7x | BOTH_FAIL |
| BYRDSPHR | 3 | 2 | Optimal | Optimal | 4.40e-10 | 23 | 12 | 78us | 2.3ms | 29.9x | PASS |
| CAMEL6 | 2 | 0 | Optimal | Optimal | 7.91e-01 | 9 | 8 | 17us | 1.7ms | 97.7x | MISMATCH |
| CANTILVR | 5 | 1 | Optimal | Optimal | 3.23e-09 | 31 | 11 | 73us | 2.5ms | 34.2x | PASS |
| CB2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 8 | 8 | 21us | 1.9ms | 90.1x | MISMATCH |
| CB3 | 3 | 3 | Optimal | Optimal | 4.30e-09 | 8 | 8 | 23us | 1.9ms | 84.1x | PASS |
| CERI651A | 7 | 61 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 102.6ms | 310us | 0.0x | BOTH_FAIL |
| CERI651ALS | 7 | 0 | Acceptable | Optimal | 8.56e-08 | 283 | 95 | 5.4ms | 16.0ms | 3.0x | PASS |
| CERI651B | 7 | 66 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 218.9ms | 294us | 0.0x | BOTH_FAIL |
| CERI651BLS | 7 | 0 | Optimal | Optimal | 1.23e-08 | 89 | 56 | 1.6ms | 9.0ms | 5.5x | PASS |
| CERI651C | 7 | 56 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 156.2ms | 316us | 0.0x | BOTH_FAIL |
| CERI651CLS | 7 | 0 | Optimal | Optimal | 3.90e-09 | 178 | 53 | 2.7ms | 7.8ms | 2.9x | PASS |
| CERI651D | 7 | 67 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 65.3ms | 298us | 0.0x | BOTH_FAIL |
| CERI651DLS | 7 | 0 | Acceptable | Optimal | 3.08e-10 | 98 | 60 | 1.7ms | 10.7ms | 6.3x | PASS |
| CERI651E | 7 | 64 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 486us | 290us | 0.6x | BOTH_FAIL |
| CERI651ELS | 7 | 0 | Acceptable | Optimal | 1.03e-09 | 83 | 45 | 1.3ms | 6.4ms | 4.8x | PASS |
| CHACONN1 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 5 | 6 | 17us | 1.3ms | 76.9x | MISMATCH |
| CHACONN2 | 3 | 3 | Optimal | Optimal | 4.49e-09 | 6 | 6 | 18us | 1.3ms | 68.7x | PASS |
| CHWIRUT1 | 3 | 214 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 745us | 286us | 0.4x | BOTH_FAIL |
| CHWIRUT1LS | 3 | 0 | Acceptable | Optimal | 9.64e-01 | 21 | 6 | 477us | 1.6ms | 3.3x | MISMATCH |
| CHWIRUT2 | 3 | 54 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 120us | 284us | 2.4x | BOTH_FAIL |
| CHWIRUT2LS | 3 | 0 | Acceptable | Optimal | 9.62e-01 | 2999 | 6 | 35.3ms | 1.3ms | 0.0x | MISMATCH |
| CLIFF | 2 | 0 | Optimal | Optimal | 7.45e-03 | 27 | 23 | 32us | 3.0ms | 91.7x | MISMATCH |
| CLUSTER | 2 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 23us | 1.6ms | 68.2x | PASS |
| CLUSTERLS | 2 | 0 | Optimal | Optimal | 2.72e-18 | 13 | 17 | 23us | 2.3ms | 103.9x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 220us | 1.6ms | 7.4x | PASS |
| CONGIGMZ | 3 | 5 | Optimal | Optimal | 3.19e-09 | 7 | 20 | 55us | 3.9ms | 70.9x | PASS |
| COOLHANS | 9 | 9 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 66us | 1.5ms | 22.7x | PASS |
| COOLHANSLS | 9 | 0 | Optimal | Optimal | 1.21e-18 | 23 | 25 | 127us | 3.9ms | 30.7x | PASS |
| CORE1 | 65 | 59 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC100 | 6 | 200 | Acceptable | Infeasible | N/A | 691 | 155 | 1.19s | 116.4ms | 0.1x | ipopt_FAIL |
| CRESC132 | 6 | 2654 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC4 | 6 | 8 | Optimal | Optimal | 1.35e-08 | 21 | 64 | 280us | 11.8ms | 42.1x | PASS |
| CRESC50 | 6 | 100 | Acceptable | Optimal | 9.65e-02 | 2999 | 194 | 2.21s | 79.8ms | 0.0x | MISMATCH |
| CSFI1 | 5 | 4 | Optimal | Optimal | 1.47e-08 | 32 | 11 | 105us | 2.2ms | 21.3x | PASS |
| CSFI2 | 5 | 4 | Optimal | Optimal | 1.48e-08 | 18 | 14 | 130us | 2.9ms | 22.6x | PASS |
| CUBE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 27 | 27 | 31us | 3.9ms | 123.8x | PASS |
| CUBENE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 14 | 1 | 28us | 470us | 17.0x | PASS |
| DALLASS | 46 | 31 | Optimal | Optimal | 1.41e-07 | 29 | 22 | 3.6ms | 4.8ms | 1.3x | PASS |
| DANIWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 19us | 277us | 14.5x | BOTH_FAIL |
| DANIWOODLS | 2 | 0 | Optimal | Optimal | 2.60e-18 | 10 | 10 | 20us | 1.7ms | 85.7x | PASS |
| DANWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 6 | 0 | 20us | 283us | 14.1x | BOTH_FAIL |
| DANWOODLS | 2 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 11 | 5us | 1.8ms | 352.1x | MISMATCH |
| DECONVB | 63 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 357.9ms | 705.7ms | 2.0x | ipopt_FAIL |
| DECONVBNE | 63 | 40 | Optimal | Optimal | 0.00e+00 | 68 | 505 | 113.5ms | 161.5ms | 1.4x | PASS |
| DECONVC | 63 | 1 | Acceptable | Optimal | 2.57e-03 | 145 | 31 | 20.5ms | 8.2ms | 0.4x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 1 | 26 | 380us | 23.7ms | 62.3x | PASS |
| DECONVU | 63 | 0 | MaxIteration | Optimal | N/A | 2999 | 333 | 447.6ms | 85.3ms | 0.2x | ripopt_FAIL |
| DEGENLPA | 20 | 15 | Acceptable | Optimal | 1.81e-03 | 44 | 18 | 611us | 3.5ms | 5.7x | MISMATCH |
| DEGENLPB | 20 | 15 | Acceptable | Optimal | 2.35e-03 | 33 | 19 | 478us | 3.5ms | 7.3x | MISMATCH |
| DEMBO7 | 16 | 20 | Acceptable | Optimal | 7.63e-05 | 259 | 45 | 6.4ms | 8.2ms | 1.3x | PASS |
| DEMYMALO | 3 | 3 | Optimal | Optimal | 2.96e-09 | 12 | 9 | 26us | 2.0ms | 74.5x | PASS |
| DENSCHNA | 2 | 0 | Optimal | Optimal | 5.88e-39 | 6 | 6 | 11us | 951us | 83.6x | PASS |
| DENSCHNB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 13us | 1.3ms | 99.8x | PASS |
| DENSCHNBNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 5 | 0 | 14us | 285us | 20.8x | ipopt_FAIL |
| DENSCHNC | 2 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 15us | 1.4ms | 92.5x | PASS |
| DENSCHNCNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 18us | 1.2ms | 64.7x | PASS |
| DENSCHND | 3 | 0 | Acceptable | Optimal | 2.22e-04 | 30 | 26 | 46us | 3.5ms | 75.2x | MISMATCH |
| DENSCHNDNE | 3 | 3 | Acceptable | Acceptable | 0.00e+00 | 19 | 22 | 39us | 3.2ms | 81.7x | PASS |
| DENSCHNE | 3 | 0 | Optimal | Optimal | 1.86e-17 | 10 | 14 | 17us | 2.4ms | 144.8x | PASS |
| DENSCHNENE | 3 | 3 | RestorationF | Infeasible | N/A | 10 | 10 | 191us | 2.0ms | 10.6x | BOTH_FAIL |
| DENSCHNF | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 14us | 1.0ms | 76.8x | PASS |
| DENSCHNFNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 16us | 935us | 58.6x | PASS |
| DEVGLA1 | 4 | 0 | Optimal | Optimal | 3.42e-23 | 21 | 23 | 149us | 3.9ms | 25.9x | PASS |
| DEVGLA1B | 4 | 0 | Acceptable | Optimal | 1.00e+00 | 28 | 20 | 187us | 4.3ms | 23.0x | MISMATCH |
| DEVGLA1NE | 4 | 24 | Optimal | IpoptStatus( | N/A | 13 | 0 | 124us | 271us | 2.2x | ipopt_FAIL |
| DEVGLA2 | 5 | 0 | Optimal | Optimal | 1.00e+00 | 16 | 13 | 135us | 2.2ms | 15.9x | MISMATCH |
| DEVGLA2B | 5 | 0 | Acceptable | Optimal | 2.58e-07 | 14 | 24 | 109us | 4.4ms | 40.1x | PASS |
| DEVGLA2NE | 5 | 16 | Optimal | IpoptStatus( | N/A | 11 | 0 | 101us | 301us | 3.0x | ipopt_FAIL |
| DGOSPEC | 3 | 0 | Acceptable | Optimal | 4.63e-03 | 10 | 27 | 21us | 4.9ms | 237.5x | MISMATCH |
| DIAMON2D | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3D | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIPIGRI | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 51us | 2.0ms | 39.0x | PASS |
| DISC2 | 29 | 23 | Optimal | Optimal | 9.11e-10 | 30 | 24 | 3.0ms | 5.3ms | 1.8x | PASS |
| DISCS | 36 | 66 | RestorationF | Optimal | N/A | 7 | 184 | 192.0ms | 72.8ms | 0.4x | ripopt_FAIL |
| DIXCHLNG | 10 | 5 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 78us | 1.6ms | 20.1x | PASS |
| DJTL | 2 | 0 | Acceptable | Acceptable | 0.00e+00 | 2999 | 1538 | 7.2ms | 161.1ms | 22.5x | PASS |
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
| DNIEPER | 61 | 24 | MaxIteration | Optimal | N/A | 2999 | 23 | 359.7ms | 4.9ms | 0.0x | ripopt_FAIL |
| DUAL1 | 85 | 1 | Acceptable | Optimal | 3.88e-06 | 12 | 15 | 3.0ms | 6.2ms | 2.1x | PASS |
| DUAL2 | 96 | 1 | Acceptable | Optimal | 8.07e-09 | 12 | 12 | 3.9ms | 5.7ms | 1.5x | PASS |
| DUAL4 | 75 | 1 | Acceptable | Optimal | 1.40e-07 | 11 | 12 | 2.0ms | 4.3ms | 2.2x | PASS |
| DUALC1 | 9 | 215 | Acceptable | Optimal | 6.76e-06 | 47 | 18 | 10.3ms | 10.5ms | 1.0x | PASS |
| DUALC2 | 7 | 229 | Acceptable | Optimal | 1.62e-06 | 17 | 12 | 4.0ms | 7.3ms | 1.8x | PASS |
| DUALC5 | 8 | 278 | Acceptable | Optimal | 1.83e-07 | 12 | 11 | 3.9ms | 7.9ms | 2.0x | PASS |
| DUALC8 | 8 | 503 | Acceptable | Optimal | 1.27e-02 | 20 | 13 | 16.2ms | 15.2ms | 0.9x | MISMATCH |
| ECKERLE4 | 3 | 35 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 60us | 282us | 4.7x | BOTH_FAIL |
| ECKERLE4LS | 3 | 0 | Acceptable | Optimal | 4.97e-01 | 25 | 36 | 126us | 5.7ms | 45.5x | MISMATCH |
| EG1 | 3 | 0 | Optimal | Optimal | 2.07e-01 | 9 | 8 | 17us | 1.8ms | 104.9x | MISMATCH |
| EGGCRATE | 2 | 0 | Optimal | Optimal | 5.00e-01 | 7 | 5 | 20us | 1.2ms | 59.6x | MISMATCH |
| EGGCRATEB | 2 | 0 | Optimal | Optimal | 5.62e-16 | 10 | 6 | 17us | 1.4ms | 79.6x | PASS |
| EGGCRATENE | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 19us | 271us | 14.3x | BOTH_FAIL |
| ELATTAR | 7 | 102 | Optimal | Optimal | 9.86e-01 | 90 | 81 | 11.5ms | 35.7ms | 3.1x | MISMATCH |
| ELATVIDU | 2 | 0 | Optimal | Optimal | 0.00e+00 | 11 | 11 | 15us | 1.6ms | 105.1x | PASS |
| ELATVIDUB | 2 | 0 | Optimal | Optimal | 1.03e-11 | 10 | 11 | 16us | 2.0ms | 126.0x | PASS |
| ELATVIDUNE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 152 | 0 | 228us | 280us | 1.2x | BOTH_FAIL |
| ENGVAL2 | 3 | 0 | Optimal | Optimal | 1.70e-20 | 20 | 21 | 34us | 3.0ms | 87.5x | PASS |
| ENGVAL2NE | 3 | 5 | Optimal | IpoptStatus( | N/A | 10 | 0 | 24us | 271us | 11.0x | ipopt_FAIL |
| ENSO | 9 | 168 | LocalInfeasi | IpoptStatus( | N/A | 39 | 0 | 3.6ms | 319us | 0.1x | BOTH_FAIL |
| ENSOLS | 9 | 0 | Acceptable | Optimal | 5.77e-16 | 13 | 7 | 1.6ms | 2.1ms | 1.3x | PASS |
| EQC | 9 | 3 | Optimal | ErrorInStepC | N/A | 2 | 15 | 23us | 4.1ms | 174.1x | ipopt_FAIL |
| ERRINBAR | 18 | 9 | Acceptable | Optimal | 2.38e-07 | 68 | 37 | 729us | 6.8ms | 9.4x | PASS |
| EXP2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 20us | 1.2ms | 62.7x | PASS |
| EXP2B | 2 | 0 | Optimal | Optimal | 2.25e-15 | 9 | 7 | 20us | 1.4ms | 68.4x | PASS |
| EXP2NE | 2 | 10 | Optimal | IpoptStatus( | N/A | 5 | 0 | 21us | 278us | 13.3x | ipopt_FAIL |
| EXPFIT | 2 | 0 | Optimal | Optimal | 8.33e-17 | 5 | 8 | 20us | 1.3ms | 65.9x | PASS |
| EXPFITA | 5 | 22 | Optimal | Optimal | 3.92e-01 | 56 | 13 | 613us | 2.6ms | 4.3x | MISMATCH |
| EXPFITB | 5 | 102 | Optimal | Optimal | 4.46e-05 | 224 | 16 | 13.6ms | 5.3ms | 0.4x | PASS |
| EXPFITC | 5 | 502 | Optimal | Optimal | 7.81e-03 | 30 | 18 | 12.7ms | 16.2ms | 1.3x | MISMATCH |
| EXPFITNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 34us | 280us | 8.3x | BOTH_FAIL |
| EXTRASIM | 2 | 1 | Optimal | Optimal | 1.82e-08 | 4 | 3 | 14us | 861us | 61.7x | PASS |
| FBRAIN | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 2.7ms | 494us | 0.2x | BOTH_FAIL |
| FBRAIN2 | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 58 | 0 | 61.9ms | 708us | 0.0x | BOTH_FAIL |
| FBRAIN2LS | 4 | 0 | Optimal | Optimal | 5.11e-10 | 10 | 10 | 7.2ms | 8.7ms | 1.2x | PASS |
| FBRAIN2NE | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 8.0ms | 678us | 0.1x | BOTH_FAIL |
| FBRAIN3 | 6 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 10.42s | 900us | 0.0x | BOTH_FAIL |
| FBRAIN3LS | 6 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 3.65s | 3.76s | 1.0x | BOTH_FAIL |
| FBRAINLS | 2 | 0 | Optimal | Optimal | 6.11e-16 | 9 | 7 | 3.4ms | 3.7ms | 1.1x | PASS |
| FBRAINNE | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 4.1ms | 508us | 0.1x | BOTH_FAIL |
| FCCU | 19 | 8 | Optimal | Optimal | 1.75e-15 | 10 | 9 | 81us | 1.8ms | 22.4x | PASS |
| FEEDLOC | 90 | 259 | RestorationF | Optimal | N/A | 8 | 23 | 1.07s | 14.0ms | 0.0x | ripopt_FAIL |
| FLETCHER | 4 | 4 | Optimal | Optimal | 4.03e-01 | 50 | 28 | 204us | 5.0ms | 24.6x | MISMATCH |
| FLT | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 17us | 1.1ms | 63.2x | PASS |
| GAUSS1 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 1.0ms | 322us | 0.3x | BOTH_FAIL |
| GAUSS1LS | 8 | 0 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 374us | 1.3ms | 3.4x | PASS |
| GAUSS2 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 1.1ms | 337us | 0.3x | BOTH_FAIL |
| GAUSS2LS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 10 | 5 | 859us | 1.2ms | 1.4x | PASS |
| GAUSS3 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 1.5ms | 314us | 0.2x | BOTH_FAIL |
| GAUSS3LS | 8 | 0 | Optimal | Optimal | 3.65e-16 | 7 | 11 | 507us | 2.5ms | 5.0x | PASS |
| GAUSSIAN | 3 | 0 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 10us | 577us | 58.2x | PASS |
| GBRAIN | 2 | 2200 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 4.0ms | 515us | 0.1x | BOTH_FAIL |
| GBRAINLS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 2.4ms | 3.0ms | 1.2x | PASS |
| GENHS28 | 10 | 8 | Optimal | Optimal | 1.22e-15 | 1 | 1 | 15us | 500us | 32.9x | PASS |
| GIGOMEZ1 | 3 | 3 | Optimal | Optimal | 1.00e+00 | 11 | 13 | 30us | 2.6ms | 84.8x | MISMATCH |
| GIGOMEZ2 | 3 | 3 | Optimal | Optimal | 5.23e-09 | 25 | 7 | 179us | 1.4ms | 8.1x | PASS |
| GIGOMEZ3 | 3 | 3 | Optimal | Optimal | 4.08e-09 | 10 | 8 | 27us | 1.6ms | 61.6x | PASS |
| GOFFIN | 51 | 50 | Acceptable | Optimal | 1.00e+00 | 2999 | 7 | 6.36s | 3.5ms | 0.0x | MISMATCH |
| GOTTFR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 15us | 971us | 64.5x | PASS |
| GOULDQP1 | 32 | 17 | Acceptable | Optimal | 4.28e-07 | 34 | 15 | 1.1ms | 3.1ms | 2.9x | PASS |
| GROUPING | 100 | 125 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 1.8ms | 292us | 0.2x | ipopt_FAIL |
| GROWTH | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 82us | 271us | 3.3x | BOTH_FAIL |
| GROWTHLS | 3 | 0 | Optimal | Optimal | 1.77e-15 | 43 | 71 | 122us | 10.5ms | 86.4x | PASS |
| GULF | 3 | 0 | Optimal | Optimal | 3.04e-22 | 22 | 28 | 622us | 5.0ms | 8.1x | PASS |
| GULFNE | 3 | 99 | LocalInfeasi | IpoptStatus( | N/A | 25 | 0 | 1.1ms | 291us | 0.3x | BOTH_FAIL |
| HAHN1 | 7 | 236 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 1.0ms | 300us | 0.3x | BOTH_FAIL |
| HAHN1LS | 7 | 0 | Acceptable | Optimal | 7.55e-02 | 2999 | 78 | 213.2ms | 16.8ms | 0.1x | MISMATCH |
| HAIFAM | 99 | 150 | Acceptable | Optimal | 8.45e-07 | 29 | 40 | 5.4ms | 15.2ms | 2.8x | PASS |
| HAIFAS | 13 | 9 | Optimal | Optimal | 9.97e-09 | 43 | 16 | 390us | 3.4ms | 8.7x | PASS |
| HAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 24 | 62 | 43us | 9.6ms | 224.5x | PASS |
| HALDMADS | 6 | 42 | Optimal | Optimal | 1.00e+00 | 34 | 8 | 2.8ms | 2.6ms | 0.9x | MISMATCH |
| HART6 | 6 | 0 | Optimal | Optimal | 1.34e-16 | 10 | 7 | 28us | 1.6ms | 59.0x | PASS |
| HATFLDA | 4 | 0 | Optimal | Optimal | 1.58e-13 | 9 | 13 | 18us | 2.2ms | 126.2x | PASS |
| HATFLDANE | 4 | 4 | Optimal | Optimal | 0.00e+00 | 10 | 6 | 25us | 1.4ms | 54.9x | PASS |
| HATFLDB | 4 | 0 | Optimal | Optimal | 3.90e-09 | 9 | 8 | 16us | 1.6ms | 94.5x | PASS |
| HATFLDBNE | 4 | 4 | MaxIteration | Infeasible | N/A | 2999 | 13 | 47.9ms | 2.8ms | 0.1x | BOTH_FAIL |
| HATFLDC | 25 | 0 | Optimal | Optimal | 6.53e-16 | 9 | 5 | 66us | 1.2ms | 18.5x | PASS |
| HATFLDCNE | 25 | 25 | Optimal | Optimal | 0.00e+00 | 7 | 4 | 160us | 1.2ms | 7.6x | PASS |
| HATFLDD | 3 | 0 | Optimal | Optimal | 1.30e-19 | 21 | 21 | 50us | 2.8ms | 55.5x | PASS |
| HATFLDDNE | 3 | 10 | LocalInfeasi | IpoptStatus( | N/A | 6 | 0 | 22us | 266us | 12.0x | BOTH_FAIL |
| HATFLDE | 3 | 0 | Optimal | Optimal | 2.72e-20 | 19 | 20 | 69us | 2.7ms | 38.9x | PASS |
| HATFLDENE | 3 | 21 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 31us | 279us | 9.1x | BOTH_FAIL |
| HATFLDF | 3 | 3 | Optimal | Optimal | 0.00e+00 | 133 | 135 | 575us | 21.9ms | 38.2x | PASS |
| HATFLDFL | 3 | 0 | Optimal | Optimal | 1.03e-08 | 494 | 1281 | 678us | 184.7ms | 272.3x | PASS |
| HATFLDFLNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 15 | 15 | 50us | 2.7ms | 54.2x | PASS |
| HATFLDFLS | 3 | 0 | Optimal | Optimal | 3.79e-18 | 18 | 36 | 30us | 5.3ms | 175.3x | PASS |
| HATFLDG | 25 | 25 | Optimal | Optimal | 0.00e+00 | 9 | 7 | 265us | 1.5ms | 5.7x | PASS |
| HATFLDGLS | 25 | 0 | Optimal | Optimal | 1.37e-16 | 13 | 14 | 135us | 2.3ms | 17.2x | PASS |
| HATFLDH | 4 | 7 | MaxIteration | Optimal | N/A | 2999 | 17 | 9.3ms | 3.1ms | 0.3x | ripopt_FAIL |
| HEART6 | 6 | 6 | NumericalErr | Optimal | N/A | 3 | 22 | 226us | 5.0ms | 22.2x | ripopt_FAIL |
| HEART6LS | 6 | 0 | Optimal | Optimal | 9.46e-23 | 1616 | 875 | 5.6ms | 132.1ms | 23.6x | PASS |
| HEART8 | 8 | 8 | Optimal | Optimal | 0.00e+00 | 4 | 12 | 89us | 2.2ms | 25.2x | PASS |
| HEART8LS | 8 | 0 | Optimal | Optimal | 3.37e-25 | 67 | 106 | 236us | 16.3ms | 68.9x | PASS |
| HELIX | 3 | 0 | Optimal | Optimal | 6.06e-25 | 8 | 13 | 14us | 2.1ms | 149.7x | PASS |
| HELIXNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 23us | 1.2ms | 54.1x | PASS |
| HET-Z | 2 | 1002 | Optimal | Optimal | 1.17e-02 | 74 | 11 | 39.7ms | 18.7ms | 0.5x | MISMATCH |
| HIELOW | 3 | 0 | Optimal | Optimal | 5.46e-15 | 7 | 8 | 9.2ms | 10.7ms | 1.2x | PASS |
| HIMMELBA | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 7us | 441us | 63.8x | PASS |
| HIMMELBB | 2 | 0 | Optimal | Optimal | 1.40e-17 | 10 | 18 | 21us | 2.7ms | 132.3x | PASS |
| HIMMELBC | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 19us | 1.1ms | 59.5x | PASS |
| HIMMELBCLS | 2 | 0 | Optimal | Optimal | 5.80e-25 | 8 | 6 | 14us | 1.1ms | 77.4x | PASS |
| HIMMELBD | 2 | 2 | RestorationF | Infeasible | N/A | 14 | 22 | 726us | 4.5ms | 6.2x | BOTH_FAIL |
| HIMMELBE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 9us | 633us | 69.3x | PASS |
| HIMMELBF | 4 | 0 | Acceptable | Optimal | 3.57e-15 | 46 | 75 | 120us | 10.2ms | 85.1x | PASS |
| HIMMELBFNE | 4 | 7 | LocalInfeasi | IpoptStatus( | N/A | 26 | 0 | 84us | 282us | 3.4x | BOTH_FAIL |
| HIMMELBG | 2 | 0 | Optimal | Optimal | 3.63e-22 | 6 | 6 | 12us | 1.2ms | 98.0x | PASS |
| HIMMELBH | 2 | 0 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 10us | 1.0ms | 101.0x | PASS |
| HIMMELBI | 100 | 12 | Optimal | Optimal | 5.80e-10 | 24 | 13 | 855us | 3.3ms | 3.8x | PASS |
| HIMMELBJ | 45 | 14 | RestorationF | ErrorInStepC | N/A | 14 | 580 | 63.0ms | 132.7ms | 2.1x | BOTH_FAIL |
| HIMMELBK | 24 | 14 | Optimal | Optimal | 4.89e-08 | 19 | 18 | 533us | 3.5ms | 6.6x | PASS |
| HIMMELP1 | 2 | 0 | Optimal | Optimal | 1.83e-15 | 11 | 10 | 16us | 1.9ms | 118.1x | PASS |
| HIMMELP2 | 2 | 1 | Optimal | Optimal | 8.68e-01 | 12 | 17 | 27us | 3.3ms | 123.3x | MISMATCH |
| HIMMELP3 | 2 | 2 | Optimal | Optimal | 8.66e-01 | 0 | 11 | 4us | 2.1ms | 546.7x | MISMATCH |
| HIMMELP4 | 2 | 3 | Optimal | Optimal | 1.23e-08 | 9 | 23 | 29us | 4.1ms | 142.5x | PASS |
| HIMMELP5 | 2 | 3 | Optimal | Optimal | 7.50e-01 | 9 | 46 | 25us | 7.5ms | 294.3x | MISMATCH |
| HIMMELP6 | 2 | 5 | Optimal | Optimal | 7.50e-01 | 11 | 31 | 37us | 5.8ms | 156.1x | MISMATCH |
| HONG | 4 | 1 | Optimal | Optimal | 4.72e-16 | 9 | 7 | 24us | 1.5ms | 63.4x | PASS |
| HS1 | 2 | 0 | Optimal | Optimal | 9.22e-21 | 25 | 28 | 30us | 4.8ms | 161.0x | PASS |
| HS10 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 12 | 30us | 2.2ms | 72.8x | PASS |
| HS100 | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 52us | 1.9ms | 37.2x | PASS |
| HS100LNP | 7 | 2 | Optimal | Optimal | 1.67e-16 | 6 | 20 | 24us | 2.7ms | 112.2x | PASS |
| HS100MOD | 7 | 4 | Optimal | Optimal | 1.73e-11 | 9 | 14 | 45us | 2.6ms | 58.7x | PASS |
| HS101 | 7 | 5 | Optimal | Optimal | 4.60e-08 | 101 | 39 | 8.7ms | 8.6ms | 1.0x | PASS |
| HS102 | 7 | 5 | Optimal | Optimal | 4.27e-08 | 24 | 52 | 231us | 9.2ms | 39.8x | PASS |
| HS103 | 7 | 5 | Optimal | Optimal | 3.36e-08 | 21 | 21 | 186us | 4.0ms | 21.7x | PASS |
| HS104 | 8 | 5 | Optimal | Optimal | 2.69e-08 | 14 | 8 | 67us | 1.7ms | 25.6x | PASS |
| HS105 | 8 | 1 | Optimal | Optimal | 9.48e-12 | 30 | 23 | 3.1ms | 6.6ms | 2.1x | PASS |
| HS106 | 8 | 6 | Optimal | Optimal | 1.74e-08 | 18 | 18 | 66us | 3.1ms | 46.5x | PASS |
| HS107 | 9 | 6 | Acceptable | Optimal | 1.76e-07 | 2506 | 7 | 13.6ms | 1.5ms | 0.1x | PASS |
| HS108 | 9 | 13 | Optimal | Optimal | 8.77e-09 | 2596 | 11 | 14.7ms | 2.4ms | 0.2x | PASS |
| HS109 | 9 | 10 | Optimal | Optimal | 1.62e-05 | 31 | 14 | 539us | 2.6ms | 4.9x | PASS |
| HS11 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 6 | 17us | 1.3ms | 76.9x | PASS |
| HS111 | 10 | 3 | Optimal | Optimal | 1.05e-11 | 11 | 15 | 85us | 2.8ms | 33.0x | PASS |
| HS111LNP | 10 | 3 | Optimal | Optimal | 1.03e-10 | 11 | 15 | 85us | 2.2ms | 26.2x | PASS |
| HS112 | 10 | 3 | Optimal | Optimal | 6.49e-14 | 10 | 10 | 79us | 1.9ms | 24.5x | PASS |
| HS113 | 10 | 8 | Optimal | Optimal | 1.68e-09 | 14 | 9 | 97us | 2.0ms | 20.4x | PASS |
| HS114 | 10 | 11 | Optimal | Optimal | 1.04e-07 | 19 | 13 | 123us | 2.4ms | 19.5x | PASS |
| HS116 | 13 | 14 | Acceptable | Optimal | 7.65e-04 | 2999 | 19 | 164.7ms | 3.5ms | 0.0x | MISMATCH |
| HS117 | 15 | 5 | Acceptable | Optimal | 3.56e-07 | 38 | 19 | 304us | 3.6ms | 11.8x | PASS |
| HS118 | 15 | 17 | Optimal | Optimal | 1.29e-09 | 624 | 10 | 9.9ms | 2.1ms | 0.2x | PASS |
| HS119 | 16 | 8 | Acceptable | Optimal | 1.05e-04 | 21 | 17 | 328us | 3.2ms | 9.7x | MISMATCH |
| HS12 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 7 | 6 | 19us | 1.3ms | 70.3x | PASS |
| HS13 | 2 | 1 | Acceptable | Optimal | 7.19e-04 | 24 | 47 | 40us | 7.5ms | 188.5x | MISMATCH |
| HS14 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 5 | 14us | 1.2ms | 85.0x | PASS |
| HS15 | 2 | 2 | Optimal | Optimal | 8.08e-08 | 11 | 13 | 36us | 2.3ms | 65.0x | PASS |
| HS16 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 10 | 10 | 21us | 2.2ms | 103.8x | MISMATCH |
| HS17 | 2 | 2 | Optimal | Optimal | 2.01e-04 | 12 | 22 | 25us | 3.7ms | 150.8x | MISMATCH |
| HS18 | 2 | 2 | Optimal | Optimal | 1.16e-09 | 15 | 10 | 34us | 1.9ms | 57.2x | PASS |
| HS19 | 2 | 2 | Optimal | Optimal | 3.09e-09 | 20 | 12 | 40us | 2.4ms | 59.1x | PASS |
| HS1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 30 | 19us | 5.9ms | 313.7x | PASS |
| HS2 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 10 | 10 | 15us | 2.0ms | 134.8x | PASS |
| HS20 | 2 | 3 | Optimal | Optimal | 6.17e-08 | 18 | 5 | 51us | 1.1ms | 22.7x | PASS |
| HS21 | 2 | 1 | Optimal | Optimal | 1.38e-10 | 9 | 6 | 20us | 1.4ms | 71.2x | PASS |
| HS21MOD | 7 | 1 | Acceptable | Optimal | 1.70e-08 | 11 | 13 | 29us | 2.4ms | 82.3x | PASS |
| HS22 | 2 | 2 | Optimal | Optimal | 1.16e-08 | 5 | 5 | 18us | 1.1ms | 62.9x | PASS |
| HS23 | 2 | 5 | MaxIteration | Optimal | N/A | 2999 | 9 | 23.9ms | 1.6ms | 0.1x | ripopt_FAIL |
| HS24 | 2 | 3 | Optimal | Optimal | 1.80e-08 | 7 | 14 | 24us | 2.8ms | 116.2x | PASS |
| HS25 | 3 | 0 | Acceptable | Optimal | 6.08e-05 | 16 | 27 | 295us | 5.4ms | 18.2x | PASS |
| HS25NE | 3 | 99 | Acceptable | IpoptStatus( | N/A | 15 | 0 | 368us | 299us | 0.8x | ipopt_FAIL |
| HS26 | 3 | 1 | Acceptable | Optimal | 2.94e-12 | 18 | 25 | 34us | 3.0ms | 89.3x | PASS |
| HS268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 14us | 2.5ms | 181.7x | MISMATCH |
| HS27 | 3 | 1 | Optimal | Optimal | 1.65e-14 | 12 | 57 | 31us | 7.7ms | 245.6x | PASS |
| HS28 | 3 | 1 | Optimal | Optimal | 9.24e-31 | 1 | 1 | 8us | 510us | 65.1x | PASS |
| HS29 | 3 | 1 | Optimal | Optimal | 3.12e-10 | 18 | 7 | 42us | 1.5ms | 34.9x | PASS |
| HS2NE | 2 | 2 | RestorationF | Infeasible | N/A | 200 | 12 | 3.7ms | 2.6ms | 0.7x | BOTH_FAIL |
| HS3 | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 12us | 928us | 80.4x | PASS |
| HS30 | 3 | 1 | Acceptable | Optimal | 1.45e-05 | 9 | 7 | 21us | 1.4ms | 69.3x | PASS |
| HS31 | 3 | 1 | Optimal | Optimal | 8.67e-09 | 11 | 6 | 23us | 1.3ms | 58.0x | PASS |
| HS32 | 3 | 2 | Acceptable | Optimal | 5.45e-05 | 12 | 15 | 29us | 2.7ms | 93.0x | PASS |
| HS33 | 3 | 2 | Acceptable | Optimal | 8.21e-08 | 14 | 9 | 37us | 1.8ms | 48.2x | PASS |
| HS34 | 3 | 2 | Optimal | Optimal | 9.46e-09 | 13 | 7 | 28us | 1.4ms | 51.6x | PASS |
| HS35 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 7 | 18us | 1.5ms | 78.9x | PASS |
| HS35I | 3 | 1 | Optimal | Optimal | 3.44e-09 | 9 | 7 | 20us | 1.4ms | 71.5x | PASS |
| HS35MOD | 3 | 1 | Optimal | Optimal | 9.62e-01 | 1 | 14 | 8us | 2.4ms | 314.0x | MISMATCH |
| HS36 | 3 | 1 | Optimal | Optimal | 6.33e-09 | 16 | 11 | 32us | 2.3ms | 71.7x | PASS |
| HS37 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 10 | 11 | 27us | 2.3ms | 83.6x | PASS |
| HS38 | 4 | 0 | Optimal | Optimal | 4.52e-11 | 37 | 39 | 59us | 6.6ms | 112.1x | PASS |
| HS39 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 32us | 1.9ms | 58.2x | PASS |
| HS3MOD | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 9us | 925us | 103.3x | PASS |
| HS4 | 2 | 0 | Optimal | Optimal | 1.97e-08 | 6 | 4 | 11us | 925us | 83.7x | PASS |
| HS40 | 4 | 3 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 14us | 701us | 49.9x | PASS |
| HS41 | 4 | 1 | Optimal | Optimal | 1.15e-09 | 11 | 7 | 23us | 1.5ms | 63.2x | PASS |
| HS42 | 4 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 18us | 812us | 44.9x | PASS |
| HS43 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 8 | 8 | 27us | 1.6ms | 58.9x | PASS |
| HS44 | 4 | 6 | Acceptable | Optimal | 2.97e-06 | 28 | 24 | 90us | 4.4ms | 49.3x | PASS |
| HS44NEW | 4 | 6 | Acceptable | Optimal | 1.48e-05 | 19 | 18 | 71us | 3.6ms | 50.6x | PASS |
| HS45 | 5 | 0 | Acceptable | Optimal | 3.14e-04 | 31 | 11 | 48us | 2.2ms | 44.5x | MISMATCH |
| HS46 | 5 | 2 | Acceptable | Optimal | 3.43e-11 | 17 | 19 | 45us | 2.4ms | 53.3x | PASS |
| HS47 | 5 | 3 | Acceptable | Optimal | 3.22e-12 | 17 | 19 | 43us | 2.5ms | 58.2x | PASS |
| HS48 | 5 | 2 | Optimal | Optimal | 4.44e-31 | 1 | 1 | 9us | 513us | 54.5x | PASS |
| HS49 | 5 | 2 | Acceptable | Optimal | 2.61e-10 | 17 | 19 | 36us | 2.6ms | 72.1x | PASS |
| HS5 | 2 | 0 | Optimal | Optimal | 2.32e-16 | 9 | 7 | 13us | 1.4ms | 106.4x | PASS |
| HS50 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 26us | 1.3ms | 52.5x | PASS |
| HS51 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 11us | 492us | 45.1x | PASS |
| HS52 | 5 | 3 | Optimal | Optimal | 1.17e-15 | 1 | 1 | 12us | 493us | 42.5x | PASS |
| HS53 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 6 | 25us | 1.2ms | 49.6x | PASS |
| HS54 | 6 | 1 | Optimal | Optimal | 7.51e-01 | 10 | 15 | 33us | 2.8ms | 83.8x | MISMATCH |
| HS55 | 6 | 6 | Optimal | Optimal | 2.04e-02 | 11 | 18 | 41us | 4.3ms | 104.7x | MISMATCH |
| HS56 | 7 | 4 | Optimal | Optimal | 2.45e-14 | 6 | 10 | 31us | 1.7ms | 56.1x | PASS |
| HS57 | 2 | 1 | Optimal | Optimal | 6.11e-15 | 10 | 10 | 40us | 1.7ms | 43.0x | PASS |
| HS59 | 2 | 3 | Acceptable | Optimal | 9.28e-04 | 10 | 17 | 42us | 3.3ms | 77.8x | MISMATCH |
| HS6 | 2 | 1 | Optimal | Optimal | 4.93e-32 | 7 | 5 | 20us | 1.1ms | 52.4x | PASS |
| HS60 | 3 | 1 | Optimal | Optimal | 1.19e-13 | 8 | 6 | 21us | 1.3ms | 61.2x | PASS |
| HS61 | 3 | 2 | Optimal | Optimal | 7.91e-16 | 9 | 10 | 21us | 1.4ms | 65.6x | PASS |
| HS62 | 3 | 1 | Optimal | Optimal | 1.38e-16 | 9 | 6 | 23us | 1.4ms | 60.2x | PASS |
| HS63 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 5 | 24us | 1.2ms | 48.9x | PASS |
| HS64 | 3 | 1 | Optimal | Optimal | 3.62e-09 | 18 | 16 | 36us | 2.9ms | 80.3x | PASS |
| HS65 | 3 | 1 | Optimal | Optimal | 6.05e-09 | 15 | 16 | 38us | 3.1ms | 82.6x | PASS |
| HS66 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 10 | 10 | 22us | 1.8ms | 80.9x | PASS |
| HS67 | 3 | 14 | Acceptable | Optimal | 3.01e-09 | 20 | 9 | 179us | 1.9ms | 10.6x | PASS |
| HS68 | 4 | 2 | Optimal | Optimal | 5.21e-11 | 18 | 16 | 50us | 2.8ms | 56.5x | PASS |
| HS69 | 4 | 2 | Optimal | Optimal | 2.85e-15 | 11 | 10 | 32us | 2.0ms | 64.1x | PASS |
| HS7 | 2 | 1 | Optimal | Optimal | 5.30e-12 | 10 | 27 | 19us | 4.0ms | 205.7x | PASS |
| HS70 | 4 | 1 | Optimal | Optimal | 1.80e-01 | 10 | 46 | 117us | 7.9ms | 68.0x | MISMATCH |
| HS71 | 4 | 2 | Optimal | Optimal | 1.01e-09 | 10 | 8 | 28us | 1.7ms | 60.0x | PASS |
| HS72 | 4 | 2 | Optimal | Optimal | 6.77e-07 | 18 | 16 | 40us | 2.8ms | 70.1x | PASS |
| HS73 | 4 | 3 | Optimal | Optimal | 1.30e-09 | 11 | 8 | 28us | 1.6ms | 56.1x | PASS |
| HS74 | 4 | 5 | Optimal | Optimal | 0.00e+00 | 11 | 8 | 37us | 1.7ms | 46.1x | PASS |
| HS75 | 4 | 5 | Optimal | Optimal | 2.54e-09 | 12 | 8 | 47us | 1.7ms | 35.2x | PASS |
| HS76 | 4 | 3 | Optimal | Optimal | 5.89e-09 | 9 | 7 | 24us | 1.5ms | 61.0x | PASS |
| HS76I | 4 | 3 | Optimal | Optimal | 3.59e-09 | 10 | 6 | 27us | 1.3ms | 48.0x | PASS |
| HS77 | 5 | 2 | Optimal | Optimal | 1.99e-11 | 9 | 11 | 27us | 1.7ms | 61.3x | PASS |
| HS78 | 5 | 3 | Optimal | Optimal | 1.52e-16 | 4 | 4 | 15us | 811us | 54.2x | PASS |
| HS79 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 18us | 839us | 46.1x | PASS |
| HS8 | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 20us | 950us | 47.1x | PASS |
| HS80 | 5 | 3 | Optimal | Optimal | 6.75e-13 | 9 | 5 | 27us | 1.1ms | 42.9x | PASS |
| HS81 | 5 | 3 | Optimal | Optimal | 3.46e-14 | 31 | 68 | 104us | 12.1ms | 115.9x | PASS |
| HS83 | 5 | 3 | MaxIteration | Optimal | N/A | 2999 | 9 | 9.2ms | 2.0ms | 0.2x | ripopt_FAIL |
| HS84 | 5 | 3 | Acceptable | Optimal | 2.27e-06 | 2999 | 9 | 164.1ms | 2.0ms | 0.0x | PASS |
| HS85 | 5 | 21 | Optimal | Optimal | 2.07e-03 | 36 | 13 | 10.0ms | 3.7ms | 0.4x | MISMATCH |
| HS86 | 5 | 10 | Optimal | Optimal | 2.58e-09 | 11 | 10 | 64us | 2.0ms | 31.3x | PASS |
| HS87 | 6 | 4 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 1.05s | 468.5ms | 0.4x | BOTH_FAIL |
| HS88 | 2 | 1 | Optimal | Optimal | 6.34e-06 | 24 | 18 | 924us | 3.5ms | 3.8x | PASS |
| HS89 | 3 | 1 | Optimal | Optimal | 3.42e-06 | 45 | 15 | 4.0ms | 3.4ms | 0.8x | PASS |
| HS9 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 13us | 742us | 57.6x | PASS |
| HS90 | 4 | 1 | Optimal | Optimal | 5.20e-06 | 18 | 16 | 1.3ms | 4.3ms | 3.4x | PASS |
| HS91 | 5 | 1 | Optimal | Optimal | 6.44e-06 | 18 | 16 | 1.6ms | 4.1ms | 2.5x | PASS |
| HS92 | 6 | 1 | Optimal | Optimal | 6.80e-06 | 16 | 35 | 1.7ms | 10.0ms | 5.8x | PASS |
| HS93 | 6 | 2 | Optimal | Optimal | 9.89e-09 | 14 | 7 | 72us | 1.9ms | 25.8x | PASS |
| HS95 | 6 | 4 | Acceptable | Optimal | 6.92e-04 | 46 | 9 | 140us | 2.0ms | 14.3x | MISMATCH |
| HS96 | 6 | 4 | Acceptable | Optimal | 6.88e-04 | 46 | 8 | 144us | 1.8ms | 12.3x | MISMATCH |
| HS97 | 6 | 4 | Acceptable | Optimal | 1.05e-06 | 48 | 24 | 221us | 4.5ms | 20.6x | PASS |
| HS98 | 6 | 4 | Acceptable | Optimal | 1.34e-04 | 155 | 13 | 488us | 2.4ms | 4.9x | MISMATCH |
| HS99 | 7 | 2 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 32us | 1.2ms | 36.1x | PASS |
| HS99EXP | 31 | 21 | RestorationF | Optimal | N/A | 6 | 17 | 61.7ms | 3.4ms | 0.1x | ripopt_FAIL |
| HUBFIT | 2 | 1 | Optimal | Optimal | 2.24e-09 | 7 | 7 | 19us | 1.8ms | 94.9x | PASS |
| HUMPS | 2 | 0 | Optimal | Optimal | 2.76e-17 | 40 | 1533 | 88us | 205.2ms | 2319.3x | PASS |
| HYDC20LS | 99 | 0 | Acceptable | Optimal | 2.98e-01 | 2999 | 639 | 972.0ms | 161.5ms | 0.2x | MISMATCH |
| HYDCAR20 | 99 | 99 | Optimal | Optimal | 0.00e+00 | 11 | 9 | 1.8ms | 2.8ms | 1.5x | PASS |
| HYDCAR6 | 29 | 29 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 227us | 1.3ms | 5.6x | PASS |
| HYDCAR6LS | 29 | 0 | Optimal | Optimal | 2.86e-18 | 1322 | 149 | 27.0ms | 28.1ms | 1.0x | PASS |
| HYPCIR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 18us | 984us | 56.0x | PASS |
| JENSMP | 2 | 0 | Optimal | Optimal | 3.43e-16 | 10 | 9 | 20us | 1.4ms | 68.8x | PASS |
| JENSMPNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 23.7ms | 270us | 0.0x | BOTH_FAIL |
| JUDGE | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 23us | 1.3ms | 56.1x | PASS |
| JUDGEB | 2 | 0 | Optimal | Optimal | 2.21e-16 | 9 | 9 | 29us | 1.7ms | 58.1x | PASS |
| JUDGENE | 2 | 20 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 44us | 267us | 6.1x | BOTH_FAIL |
| KIRBY2 | 5 | 151 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 121.3ms | 280us | 0.0x | BOTH_FAIL |
| KIRBY2LS | 5 | 0 | Acceptable | Optimal | 1.21e-14 | 26 | 11 | 580us | 2.1ms | 3.6x | PASS |
| KIWCRESC | 3 | 2 | Optimal | Optimal | 1.07e-08 | 13 | 8 | 30us | 1.7ms | 57.9x | PASS |
| KOEBHELB | 3 | 0 | Optimal | Optimal | 7.33e-16 | 843 | 71 | 12.3ms | 14.0ms | 1.1x | PASS |
| KOEBHELBNE | 3 | 156 | LocalInfeasi | IpoptStatus( | N/A | 45 | 0 | 968us | 281us | 0.3x | BOTH_FAIL |
| KOWOSB | 4 | 0 | Optimal | Optimal | 5.96e-19 | 7 | 8 | 24us | 1.7ms | 72.0x | PASS |
| KOWOSBNE | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 44us | 274us | 6.1x | BOTH_FAIL |
| KSIP | 20 | 1001 | Optimal | Optimal | 9.79e-01 | 1 | 22 | 11.7ms | 67.9ms | 5.8x | MISMATCH |
| LAKES | 90 | 78 | Optimal | Optimal | 1.53e-13 | 13 | 11 | 957us | 3.2ms | 3.3x | PASS |
| LANCZOS1 | 6 | 24 | Optimal | IpoptStatus( | N/A | 13 | 0 | 108us | 283us | 2.6x | ipopt_FAIL |
| LANCZOS1LS | 6 | 0 | Acceptable | Optimal | 5.33e-09 | 105 | 115 | 596us | 18.4ms | 31.0x | PASS |
| LANCZOS2 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 12 | 0 | 99us | 279us | 2.8x | ipopt_FAIL |
| LANCZOS2LS | 6 | 0 | Acceptable | Optimal | 7.41e-09 | 97 | 101 | 568us | 16.3ms | 28.8x | PASS |
| LANCZOS3 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 12 | 0 | 98us | 285us | 2.9x | ipopt_FAIL |
| LANCZOS3LS | 6 | 0 | Acceptable | Optimal | 1.08e-08 | 78 | 174 | 460us | 28.0ms | 60.9x | PASS |
| LAUNCH | 25 | 28 | Acceptable | Optimal | 1.26e-06 | 72 | 12 | 2.9ms | 2.7ms | 0.9x | PASS |
| LEVYMONE10 | 10 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 48us | 278us | 5.8x | BOTH_FAIL |
| LEVYMONE5 | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 17us | 282us | 16.1x | BOTH_FAIL |
| LEVYMONE6 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 25us | 286us | 11.2x | BOTH_FAIL |
| LEVYMONE7 | 4 | 8 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 25us | 287us | 11.5x | BOTH_FAIL |
| LEVYMONE8 | 5 | 10 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 26us | 273us | 10.4x | BOTH_FAIL |
| LEVYMONE9 | 8 | 16 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 40us | 301us | 7.5x | BOTH_FAIL |
| LEVYMONT10 | 10 | 0 | Optimal | Optimal | 0.00e+00 | 8 | 4 | 28us | 945us | 33.3x | PASS |
| LEVYMONT5 | 2 | 0 | Acceptable | Optimal | 1.00e+00 | 10 | 10 | 20us | 2.0ms | 101.1x | MISMATCH |
| LEVYMONT6 | 3 | 0 | Optimal | Optimal | 0.00e+00 | 10 | 8 | 20us | 1.6ms | 83.0x | PASS |
| LEVYMONT7 | 4 | 0 | Optimal | Optimal | 1.42e-16 | 10 | 7 | 23us | 1.6ms | 71.3x | PASS |
| LEVYMONT8 | 5 | 0 | Optimal | Optimal | 1.64e-16 | 8 | 4 | 18us | 915us | 49.8x | PASS |
| LEVYMONT9 | 8 | 0 | Optimal | Optimal | 1.71e-16 | 8 | 4 | 24us | 944us | 40.0x | PASS |
| LEWISPOL | 6 | 9 | Acceptable | IpoptStatus( | N/A | 10 | 0 | 924us | 279us | 0.3x | ipopt_FAIL |
| LHAIFAM | 99 | 150 | MaxIteration | InvalidNumbe | N/A | 2999 | 0 | 911.0ms | 304us | 0.0x | BOTH_FAIL |
| LIN | 4 | 2 | Optimal | Optimal | 2.03e-03 | 9 | 7 | 30us | 1.4ms | 45.2x | MISMATCH |
| LINSPANH | 97 | 33 | Acceptable | Optimal | 5.90e-07 | 2999 | 24 | 5.78s | 5.2ms | 0.0x | PASS |
| LOADBAL | 31 | 31 | Optimal | Optimal | 6.55e-09 | 15 | 13 | 848us | 3.0ms | 3.5x | PASS |
| LOGHAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 55 | 2747 | 110us | 382.3ms | 3467.1x | PASS |
| LOGROS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 50 | 49 | 59us | 9.5ms | 160.4x | PASS |
| LOOTSMA | 3 | 2 | Optimal | Optimal | 9.02e-08 | 26 | 13 | 62us | 2.6ms | 42.0x | PASS |
| LOTSCHD | 12 | 7 | Optimal | Optimal | 4.47e-10 | 15 | 9 | 80us | 1.8ms | 22.6x | PASS |
| LRCOVTYPE | 54 | 0 | Optimal | N/A | N/A | 65 | 0 | 14.39s | N/A | N/A | ipopt_FAIL |
| LRIJCNN1 | 22 | 0 | Optimal | Optimal | 1.16e-14 | 18 | 11 | 336.2ms | 180.6ms | 0.5x | PASS |
| LSC1 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 26us | 278us | 10.5x | BOTH_FAIL |
| LSC1LS | 3 | 0 | Acceptable | Optimal | 2.30e-15 | 17 | 16 | 40us | 2.6ms | 64.7x | PASS |
| LSC2 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 129 | 0 | 275us | 282us | 1.0x | BOTH_FAIL |
| LSC2LS | 3 | 0 | Acceptable | Optimal | 2.32e-04 | 32 | 38 | 55us | 4.9ms | 90.0x | MISMATCH |
| LSNNODOC | 5 | 4 | Acceptable | Optimal | 3.42e-06 | 13 | 10 | 41us | 2.0ms | 49.8x | PASS |
| LSQFIT | 2 | 1 | Optimal | Optimal | 4.24e-09 | 7 | 7 | 17us | 1.5ms | 84.7x | PASS |
| MADSEN | 3 | 6 | Optimal | Optimal | 9.76e-09 | 35 | 18 | 136us | 3.3ms | 24.5x | PASS |
| MAKELA1 | 3 | 2 | Optimal | Optimal | 2.00e+00 | 6 | 12 | 25us | 2.4ms | 96.7x | MISMATCH |
| MAKELA2 | 3 | 3 | Optimal | Optimal | 1.27e-01 | 3 | 6 | 15us | 1.3ms | 90.9x | MISMATCH |
| MAKELA3 | 21 | 20 | Acceptable | Optimal | 8.52e-09 | 25 | 11 | 410us | 2.5ms | 6.0x | PASS |
| MAKELA4 | 21 | 40 | Optimal | Optimal | 6.88e-01 | 1 | 5 | 65us | 1.5ms | 22.8x | MISMATCH |
| MARATOS | 2 | 1 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 12us | 822us | 66.7x | PASS |
| MARATOSB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 670 | 672 | 738us | 98.7ms | 133.8x | PASS |
| MATRIX2 | 6 | 2 | Acceptable | Optimal | 7.84e-10 | 12 | 42 | 38us | 6.7ms | 177.2x | PASS |
| MAXLIKA | 8 | 0 | Acceptable | Optimal | 1.13e-02 | 20 | 23 | 1.9ms | 6.4ms | 3.4x | MISMATCH |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 222us | 1.5ms | 6.9x | PASS |
| MDHOLE | 2 | 0 | Optimal | Optimal | 9.98e-09 | 35 | 42 | 45us | 8.7ms | 193.8x | PASS |
| MESH | 41 | 48 | Acceptable | IpoptStatus( | N/A | 19 | 79 | 4.1ms | 19.5ms | 4.8x | ipopt_FAIL |
| METHANB8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 168us | 872us | 5.2x | PASS |
| METHANB8LS | 31 | 0 | Optimal | Optimal | 5.35e-26 | 9 | 8 | 179us | 1.5ms | 8.3x | PASS |
| METHANL8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 197us | 1.0ms | 5.1x | PASS |
| METHANL8LS | 31 | 0 | Optimal | Optimal | 5.70e-17 | 977 | 40 | 21.5ms | 7.6ms | 0.4x | PASS |
| MEXHAT | 2 | 0 | Optimal | Optimal | 6.77e-10 | 28 | 26 | 34us | 3.7ms | 107.8x | PASS |
| MEYER3 | 3 | 0 | Acceptable | Optimal | 2.23e-12 | 2999 | 194 | 9.9ms | 28.7ms | 2.9x | PASS |
| MEYER3NE | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 12.4ms | 285us | 0.0x | BOTH_FAIL |
| MGH09 | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 130us | 277us | 2.1x | BOTH_FAIL |
| MGH09LS | 4 | 0 | Optimal | Optimal | 5.96e-19 | 58 | 72 | 147us | 11.0ms | 75.0x | PASS |
| MGH10 | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 1 | 0 | 18us | 277us | 15.6x | BOTH_FAIL |
| MGH10LS | 3 | 0 | Acceptable | Optimal | 4.94e-12 | 2999 | 1828 | 9.7ms | 270.3ms | 27.8x | PASS |
| MGH10S | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 1 | 0 | 13us | 289us | 22.8x | BOTH_FAIL |
| MGH10SLS | 3 | 0 | MaxIteration | Optimal | N/A | 2999 | 354 | 8.4ms | 51.4ms | 6.1x | ripopt_FAIL |
| MGH17 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 101us | 285us | 2.8x | BOTH_FAIL |
| MGH17LS | 5 | 0 | Acceptable | Optimal | 6.06e-07 | 37 | 47 | 231us | 8.3ms | 36.0x | PASS |
| MGH17S | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 101us | 281us | 2.8x | BOTH_FAIL |
| MGH17SLS | 5 | 0 | Optimal | Optimal | 2.45e-02 | 54 | 41 | 341us | 7.6ms | 22.2x | MISMATCH |
| MIFFLIN1 | 3 | 2 | Optimal | Optimal | 9.24e-09 | 6 | 5 | 20us | 1.1ms | 57.9x | PASS |
| MIFFLIN2 | 3 | 2 | Optimal | Optimal | 9.97e-09 | 15 | 11 | 30us | 2.3ms | 76.0x | PASS |
| MINMAXBD | 5 | 20 | Optimal | Optimal | 8.54e-11 | 762 | 25 | 10.4ms | 5.8ms | 0.6x | PASS |
| MINMAXRB | 3 | 4 | Optimal | Optimal | 9.85e-09 | 3 | 8 | 16us | 1.7ms | 105.3x | PASS |
| MINSURF | 64 | 0 | MaxIteration | Optimal | N/A | 2999 | 4 | 200.4ms | 1.2ms | 0.0x | ripopt_FAIL |
| MISRA1A | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 69us | 278us | 4.0x | BOTH_FAIL |
| MISRA1ALS | 2 | 0 | Acceptable | Optimal | 1.91e-14 | 35 | 40 | 91us | 6.1ms | 67.5x | PASS |
| MISRA1B | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 60us | 282us | 4.7x | BOTH_FAIL |
| MISRA1BLS | 2 | 0 | Optimal | Optimal | 4.76e-14 | 25 | 34 | 55us | 4.9ms | 88.6x | PASS |
| MISRA1C | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 48us | 279us | 5.8x | BOTH_FAIL |
| MISRA1CLS | 2 | 0 | Acceptable | Optimal | 4.63e-14 | 18 | 14 | 50us | 2.4ms | 49.0x | PASS |
| MISRA1D | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 51us | 290us | 5.7x | BOTH_FAIL |
| MISRA1DLS | 2 | 0 | Optimal | Optimal | 9.02e-17 | 20 | 30 | 45us | 4.4ms | 97.7x | PASS |
| MISTAKE | 9 | 13 | Acceptable | Optimal | 5.00e-01 | 81 | 16 | 949us | 3.3ms | 3.5x | MISMATCH |
| MRIBASIS | 36 | 55 | Acceptable | Optimal | 1.00e-08 | 66 | 15 | 24.6ms | 4.3ms | 0.2x | PASS |
| MSS1 | 90 | 73 | Acceptable | Optimal | 1.25e-01 | 211 | 95 | 82.6ms | 51.9ms | 0.6x | MISMATCH |
| MUONSINE | 1 | 512 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 528.2ms | 332us | 0.0x | BOTH_FAIL |
| MUONSINELS | 1 | 0 | Acceptable | Optimal | 1.36e-01 | 11 | 8 | 1.2ms | 1.6ms | 1.4x | MISMATCH |
| MWRIGHT | 5 | 3 | Optimal | Optimal | 9.48e-01 | 13 | 10 | 38us | 1.6ms | 41.7x | MISMATCH |
| NASH | 72 | 24 | RestorationF | Infeasible | N/A | 7 | 45 | 12.4ms | 11.5ms | 0.9x | BOTH_FAIL |
| NELSON | 3 | 128 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 322us | 293us | 0.9x | BOTH_FAIL |
| NET1 | 48 | 57 | RestorationF | Optimal | N/A | 14 | 26 | 1.06s | 5.5ms | 0.0x | ripopt_FAIL |
| NYSTROM5 | 18 | 20 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 59.0ms | 274us | 0.0x | BOTH_FAIL |
| NYSTROM5C | 18 | 20 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 61.2ms | 278us | 0.0x | BOTH_FAIL |
| ODFITS | 10 | 6 | Optimal | Optimal | 1.91e-16 | 11 | 8 | 49us | 1.7ms | 35.5x | PASS |
| OET1 | 3 | 1002 | Optimal | Optimal | 1.05e-01 | 13 | 33 | 10.1ms | 46.8ms | 4.6x | MISMATCH |
| OET2 | 3 | 1002 | MaxIteration | Optimal | N/A | 2999 | 181 | 1.57s | 288.1ms | 0.2x | ripopt_FAIL |
| OET3 | 4 | 1002 | Optimal | Optimal | 1.00e+00 | 5 | 13 | 8.6ms | 22.5ms | 2.6x | MISMATCH |
| OET4 | 4 | 1002 | Acceptable | Optimal | 9.98e-01 | 368 | 165 | 330.6ms | 258.7ms | 0.8x | MISMATCH |
| OET5 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET6 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET7 | 7 | 1002 | Acceptable | Optimal | 1.00e+00 | 166 | 193 | 145.2ms | 532.5ms | 3.7x | MISMATCH |
| OPTCNTRL | 32 | 20 | MaxIteration | Optimal | N/A | 2999 | 9 | 95.5ms | 1.9ms | 0.0x | ripopt_FAIL |
| OPTPRLOC | 30 | 30 | Optimal | Optimal | 1.22e-08 | 57 | 13 | 2.6ms | 3.2ms | 1.2x | PASS |
| ORTHREGB | 27 | 6 | Optimal | Optimal | 4.20e-19 | 2 | 2 | 47us | 739us | 15.7x | PASS |
| OSBORNE1 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 68us | 278us | 4.1x | BOTH_FAIL |
| OSBORNE2 | 11 | 65 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 405us | 274us | 0.7x | BOTH_FAIL |
| OSBORNEA | 5 | 0 | MaxIteration | Optimal | N/A | 2999 | 64 | 18.0ms | 10.0ms | 0.6x | ripopt_FAIL |
| OSBORNEB | 11 | 0 | Optimal | Optimal | 4.86e-17 | 16 | 19 | 446us | 3.2ms | 7.1x | PASS |
| OSLBQP | 8 | 0 | Acceptable | Optimal | 7.24e-07 | 13 | 15 | 22us | 2.6ms | 114.8x | PASS |
| PALMER1 | 4 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 13 | 61us | 2.5ms | 40.9x | PASS |
| PALMER1A | 6 | 0 | Optimal | Optimal | 1.98e-14 | 47 | 48 | 291us | 9.2ms | 31.5x | PASS |
| PALMER1ANE | 6 | 35 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 187us | 317us | 1.7x | BOTH_FAIL |
| PALMER1B | 4 | 0 | Optimal | Optimal | 4.12e-15 | 18 | 17 | 92us | 3.0ms | 32.9x | PASS |
| PALMER1BNE | 4 | 35 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 81us | 282us | 3.5x | BOTH_FAIL |
| PALMER1C | 8 | 0 | Optimal | Optimal | 1.79e-13 | 4 | 1 | 32us | 445us | 14.1x | PASS |
| PALMER1D | 7 | 0 | Optimal | Optimal | 3.71e-14 | 2 | 1 | 17us | 471us | 27.6x | PASS |
| PALMER1E | 8 | 0 | Optimal | Optimal | 4.05e-13 | 119 | 55 | 971us | 10.3ms | 10.6x | PASS |
| PALMER1ENE | 8 | 35 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 266us | 318us | 1.2x | BOTH_FAIL |
| PALMER1NE | 4 | 31 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 209us | 277us | 1.3x | BOTH_FAIL |
| PALMER2 | 4 | 0 | Optimal | Optimal | 4.98e-16 | 17 | 28 | 75us | 6.2ms | 83.1x | PASS |
| PALMER2A | 6 | 0 | Optimal | Optimal | 7.77e-16 | 71 | 91 | 323us | 18.6ms | 57.8x | PASS |
| PALMER2ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 61us | 278us | 4.5x | BOTH_FAIL |
| PALMER2B | 4 | 0 | Optimal | Optimal | 5.44e-15 | 17 | 15 | 64us | 3.2ms | 48.9x | PASS |
| PALMER2BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 48us | 284us | 5.9x | BOTH_FAIL |
| PALMER2C | 8 | 0 | Optimal | Optimal | 5.44e-15 | 1 | 1 | 12us | 457us | 39.1x | PASS |
| PALMER2E | 8 | 0 | Optimal | Optimal | 2.09e-12 | 113 | 114 | 672us | 23.7ms | 35.2x | PASS |
| PALMER2ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 66 | 0 | 588us | 283us | 0.5x | BOTH_FAIL |
| PALMER2NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 379 | 0 | 2.3ms | 293us | 0.1x | BOTH_FAIL |
| PALMER3 | 4 | 0 | Acceptable | Optimal | 6.25e-02 | 20 | 44 | 93us | 7.3ms | 78.5x | MISMATCH |
| PALMER3A | 6 | 0 | Optimal | Optimal | 3.89e-16 | 66 | 73 | 297us | 14.0ms | 47.2x | PASS |
| PALMER3ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 100us | 282us | 2.8x | BOTH_FAIL |
| PALMER3B | 4 | 0 | Optimal | Optimal | 1.47e-15 | 15 | 15 | 64us | 3.1ms | 47.9x | PASS |
| PALMER3BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 59us | 284us | 4.8x | BOTH_FAIL |
| PALMER3C | 8 | 0 | Optimal | Optimal | 3.09e-15 | 1 | 1 | 11us | 444us | 41.2x | PASS |
| PALMER3E | 8 | 0 | Optimal | Optimal | 1.94e-13 | 28 | 32 | 165us | 5.5ms | 33.4x | PASS |
| PALMER3ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 444 | 0 | 3.2ms | 288us | 0.1x | BOTH_FAIL |
| PALMER3NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 20 | 0 | 109us | 297us | 2.7x | BOTH_FAIL |
| PALMER4 | 4 | 0 | Optimal | Optimal | 5.72e-02 | 29 | 16 | 121us | 3.4ms | 28.0x | MISMATCH |
| PALMER4A | 6 | 0 | Optimal | Optimal | 1.79e-15 | 48 | 53 | 225us | 10.0ms | 44.4x | PASS |
| PALMER4ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 65us | 274us | 4.2x | BOTH_FAIL |
| PALMER4B | 4 | 0 | Optimal | Optimal | 1.82e-15 | 14 | 16 | 52us | 3.3ms | 63.9x | PASS |
| PALMER4BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 63us | 275us | 4.4x | BOTH_FAIL |
| PALMER4C | 8 | 0 | Optimal | Optimal | 3.29e-15 | 1 | 1 | 12us | 441us | 38.2x | PASS |
| PALMER4E | 8 | 0 | Optimal | Optimal | 8.84e-16 | 25 | 25 | 153us | 4.7ms | 30.8x | PASS |
| PALMER4ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 129us | 270us | 2.1x | BOTH_FAIL |
| PALMER4NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 26 | 0 | 146us | 287us | 2.0x | BOTH_FAIL |
| PALMER5A | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 12.8ms | 632.0ms | 49.5x | BOTH_FAIL |
| PALMER5ANE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 17.7ms | 291us | 0.0x | BOTH_FAIL |
| PALMER5B | 9 | 0 | Acceptable | Optimal | 1.57e-13 | 57 | 113 | 251us | 21.0ms | 83.6x | PASS |
| PALMER5BNE | 9 | 12 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 59us | 305us | 5.2x | BOTH_FAIL |
| PALMER5C | 6 | 0 | Optimal | Optimal | 2.50e-15 | 1 | 1 | 8us | 435us | 56.8x | PASS |
| PALMER5D | 4 | 0 | Optimal | Optimal | 3.25e-16 | 1 | 1 | 6us | 450us | 72.5x | PASS |
| PALMER5E | 8 | 0 | Acceptable | MaxIteration | N/A | 14 | 3000 | 60us | 475.8ms | 7930.4x | ipopt_FAIL |
| PALMER5ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 17.9ms | 270us | 0.0x | BOTH_FAIL |
| PALMER6A | 6 | 0 | Optimal | Optimal | 2.42e-15 | 115 | 105 | 360us | 19.2ms | 53.3x | PASS |
| PALMER6ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 156us | 278us | 1.8x | BOTH_FAIL |
| PALMER6C | 8 | 0 | Optimal | Optimal | 2.30e-15 | 1 | 1 | 10us | 450us | 45.0x | PASS |
| PALMER6E | 8 | 0 | Optimal | Optimal | 2.54e-11 | 37 | 30 | 152us | 5.7ms | 37.2x | PASS |
| PALMER6ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 26 | 0 | 159us | 286us | 1.8x | BOTH_FAIL |
| PALMER7A | 6 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 9.6ms | 485.1ms | 50.7x | ipopt_FAIL |
| PALMER7ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 262 | 0 | 1.3ms | 283us | 0.2x | BOTH_FAIL |
| PALMER7C | 8 | 0 | Optimal | Optimal | 2.59e-13 | 2 | 1 | 11us | 443us | 38.9x | PASS |
| PALMER7E | 8 | 0 | MaxIteration | MaxIteration | N/A | 2999 | 3000 | 13.0ms | 624.0ms | 47.9x | BOTH_FAIL |
| PALMER7ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 13.3ms | 277us | 0.0x | BOTH_FAIL |
| PALMER8A | 6 | 0 | Optimal | Optimal | 1.90e-15 | 28 | 36 | 92us | 7.5ms | 80.6x | PASS |
| PALMER8ANE | 6 | 12 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 78us | 283us | 3.6x | BOTH_FAIL |
| PALMER8C | 8 | 0 | Optimal | Optimal | 8.33e-15 | 1 | 1 | 9us | 437us | 51.1x | PASS |
| PALMER8E | 8 | 0 | Optimal | Optimal | 1.65e-17 | 30 | 23 | 121us | 4.2ms | 34.9x | PASS |
| PALMER8ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 20 | 0 | 118us | 272us | 2.3x | BOTH_FAIL |
| PARKCH | 15 | 0 | Acceptable | Optimal | 1.27e-14 | 19 | 17 | 4.59s | 3.84s | 0.8x | PASS |
| PENTAGON | 6 | 15 | Optimal | Optimal | 1.47e-05 | 13 | 19 | 132us | 3.9ms | 29.8x | PASS |
| PFIT1 | 3 | 3 | Optimal | Infeasible | N/A | 8 | 266 | 114us | 43.9ms | 383.9x | ipopt_FAIL |
| PFIT1LS | 3 | 0 | Optimal | Optimal | 1.64e-20 | 226 | 263 | 380us | 44.7ms | 117.8x | PASS |
| PFIT2 | 3 | 3 | Optimal | RestorationF | N/A | 14 | 247 | 45us | 45.1ms | 995.4x | ipopt_FAIL |
| PFIT2LS | 3 | 0 | Optimal | Optimal | 1.47e-20 | 106 | 82 | 179us | 13.8ms | 77.1x | PASS |
| PFIT3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 29 | 133 | 302us | 24.9ms | 82.5x | PASS |
| PFIT3LS | 3 | 0 | Optimal | Optimal | 2.90e-20 | 119 | 132 | 202us | 22.2ms | 110.1x | PASS |
| PFIT4 | 3 | 3 | NumericalErr | Optimal | N/A | 391 | 190 | 13.3ms | 35.0ms | 2.6x | ripopt_FAIL |
| PFIT4LS | 3 | 0 | Optimal | Optimal | 6.57e-20 | 206 | 215 | 358us | 37.2ms | 103.9x | PASS |
| POLAK1 | 3 | 2 | Optimal | Optimal | 3.67e-09 | 9 | 5 | 25us | 1.2ms | 46.3x | PASS |
| POLAK2 | 11 | 2 | Optimal | Optimal | 3.59e-11 | 30 | 10 | 124us | 2.2ms | 17.4x | PASS |
| POLAK3 | 12 | 10 | Optimal | MaxIteration | N/A | 39 | 3000 | 1.3ms | 656.3ms | 515.5x | ipopt_FAIL |
| POLAK4 | 3 | 3 | Acceptable | Optimal | 4.53e-09 | 14 | 4 | 29us | 1.0ms | 35.0x | PASS |
| POLAK5 | 3 | 2 | Acceptable | Optimal | 1.73e-10 | 30 | 31 | 95us | 5.1ms | 53.2x | PASS |
| POLAK6 | 5 | 4 | Optimal | MaxIteration | N/A | 13 | 3000 | 50us | 833.0ms | 16757.1x | ipopt_FAIL |
| PORTFL1 | 12 | 1 | Acceptable | Optimal | 1.28e-06 | 10 | 9 | 217us | 1.9ms | 8.7x | PASS |
| PORTFL2 | 12 | 1 | Acceptable | Optimal | 4.37e-07 | 10 | 8 | 196us | 1.8ms | 9.1x | PASS |
| PORTFL3 | 12 | 1 | Acceptable | Optimal | 1.50e-06 | 10 | 9 | 201us | 2.0ms | 9.8x | PASS |
| PORTFL4 | 12 | 1 | Acceptable | Optimal | 3.83e-06 | 9 | 8 | 180us | 1.8ms | 9.8x | PASS |
| PORTFL6 | 12 | 1 | Acceptable | Optimal | 4.91e-07 | 10 | 8 | 203us | 1.7ms | 8.5x | PASS |
| POWELLBS | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 11 | 24us | 1.5ms | 63.1x | PASS |
| POWELLBSLS | 2 | 0 | Optimal | Optimal | 6.26e-26 | 90 | 91 | 109us | 12.3ms | 112.9x | PASS |
| POWELLSQ | 2 | 2 | Optimal | Infeasible | N/A | 8 | 29 | 695us | 5.0ms | 7.2x | ipopt_FAIL |
| POWELLSQLS | 2 | 0 | Acceptable | Optimal | 6.89e-11 | 525 | 10 | 564us | 1.7ms | 3.1x | PASS |
| PRICE3NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 18us | 1.2ms | 67.4x | PASS |
| PRICE4 | 2 | 0 | Optimal | Optimal | 1.01e-22 | 13 | 8 | 21us | 1.3ms | 62.8x | PASS |
| PRICE4B | 2 | 0 | Optimal | Optimal | 3.11e-12 | 10 | 8 | 18us | 1.7ms | 97.4x | PASS |
| PRICE4NE | 2 | 2 | Acceptable | Acceptable | 0.00e+00 | 19 | 23 | 30us | 3.3ms | 111.3x | PASS |
| PRODPL0 | 60 | 29 | Acceptable | Optimal | 1.33e-06 | 20 | 15 | 2.9ms | 3.4ms | 1.2x | PASS |
| PRODPL1 | 60 | 29 | Acceptable | Optimal | 9.38e-06 | 44 | 28 | 6.3ms | 6.3ms | 1.0x | PASS |
| PSPDOC | 4 | 0 | Optimal | Optimal | 3.50e-09 | 7 | 5 | 14us | 1.2ms | 91.0x | PASS |
| PT | 2 | 501 | Optimal | Optimal | 8.94e-03 | 4 | 106 | 1.5ms | 78.7ms | 51.3x | MISMATCH |
| QC | 9 | 4 | Acceptable | Optimal | 1.61e-02 | 12 | 44 | 121us | 8.3ms | 68.5x | MISMATCH |
| QCNEW | 9 | 3 | Optimal | Optimal | 2.48e-03 | 8 | 6 | 135us | 1.3ms | 9.8x | MISMATCH |
| QPCBLEND | 83 | 74 | Optimal | Optimal | 4.11e-07 | 39 | 19 | 2.9ms | 5.0ms | 1.7x | PASS |
| QPNBLEND | 83 | 74 | Optimal | Optimal | 4.77e-07 | 50 | 18 | 4.0ms | 4.8ms | 1.2x | PASS |
| RAT42 | 3 | 9 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 80us | 300us | 3.7x | BOTH_FAIL |
| RAT42LS | 3 | 0 | Optimal | Optimal | 8.82e-16 | 22 | 28 | 66us | 4.1ms | 62.2x | PASS |
| RAT43 | 4 | 15 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 205us | 286us | 1.4x | BOTH_FAIL |
| RAT43LS | 4 | 0 | Acceptable | Optimal | 9.65e-01 | 2999 | 34 | 26.3ms | 5.0ms | 0.2x | MISMATCH |
| RECIPE | 3 | 3 | Acceptable | Optimal | 0.00e+00 | 14 | 16 | 28us | 2.5ms | 88.1x | PASS |
| RECIPELS | 3 | 0 | Acceptable | Optimal | 1.17e-10 | 19 | 29 | 28us | 4.6ms | 162.2x | PASS |
| RES | 20 | 14 | Acceptable | Optimal | 0.00e+00 | 12 | 10 | 164us | 1.8ms | 11.2x | PASS |
| RK23 | 17 | 11 | Acceptable | Optimal | 6.20e-05 | 17 | 10 | 196us | 2.4ms | 12.0x | PASS |
| ROBOT | 14 | 2 | Acceptable | IpoptStatus( | N/A | 2999 | 18 | 11.3ms | 3.9ms | 0.3x | ipopt_FAIL |
| ROSENBR | 2 | 0 | Optimal | Optimal | 0.00e+00 | 21 | 21 | 24us | 3.2ms | 133.3x | PASS |
| ROSENBRTU | 2 | 0 | Optimal | Optimal | 1.52e-24 | 45 | 87 | 54us | 13.0ms | 243.3x | PASS |
| ROSENMMX | 5 | 4 | Optimal | Optimal | 2.27e-10 | 11 | 13 | 42us | 2.7ms | 65.7x | PASS |
| ROSZMAN1 | 4 | 25 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 52us | 293us | 5.6x | BOTH_FAIL |
| ROSZMAN1LS | 4 | 0 | Optimal | Optimal | 7.59e-19 | 52 | 28 | 215us | 4.4ms | 20.5x | PASS |
| RSNBRNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 17 | 1 | 39us | 473us | 12.1x | PASS |
| S268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 14us | 2.6ms | 181.3x | MISMATCH |
| S308 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 14us | 1.5ms | 107.7x | PASS |
| S308NE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 27 | 0 | 42us | 281us | 6.7x | BOTH_FAIL |
| S316-322 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 8 | 7 | 19us | 1.3ms | 67.5x | PASS |
| S365 | 7 | 5 | Optimal | RestorationF | N/A | 1 | 1 | 26us | 802us | 30.6x | ipopt_FAIL |
| S365MOD | 7 | 5 | Optimal | RestorationF | N/A | 1 | 1 | 27us | 815us | 30.5x | ipopt_FAIL |
| SANTA | 21 | 23 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 223us | 272us | 1.2x | BOTH_FAIL |
| SANTALS | 21 | 0 | Optimal | Optimal | 2.36e-09 | 35 | 31 | 383us | 6.7ms | 17.4x | PASS |
| SIM2BQP | 2 | 0 | Optimal | Optimal | 7.75e-09 | 4 | 5 | 12us | 1.1ms | 88.6x | PASS |
| SIMBQP | 2 | 0 | Optimal | Optimal | 8.29e-09 | 6 | 5 | 11us | 1.1ms | 95.2x | PASS |
| SIMPLLPA | 2 | 2 | Optimal | Optimal | 1.18e-08 | 6 | 8 | 19us | 1.7ms | 90.0x | PASS |
| SIMPLLPB | 2 | 3 | Optimal | Optimal | 9.07e-09 | 3 | 7 | 14us | 1.5ms | 101.2x | PASS |
| SINEVAL | 2 | 0 | Optimal | Optimal | 1.17e-40 | 42 | 42 | 50us | 6.4ms | 128.8x | PASS |
| SINVALNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 12 | 1 | 29us | 484us | 16.9x | PASS |
| SIPOW1 | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 81 | 28.0ms | 211.0ms | 7.5x | MISMATCH |
| SIPOW1M | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 88 | 30.1ms | 228.8ms | 7.6x | MISMATCH |
| SIPOW2 | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 69 | 25.7ms | 171.9ms | 6.7x | MISMATCH |
| SIPOW2M | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 73 | 29.9ms | 178.5ms | 6.0x | MISMATCH |
| SIPOW3 | 4 | 2000 | Optimal | Optimal | 7.10e-01 | 1 | 12 | 26.7ms | 35.7ms | 1.3x | MISMATCH |
| SIPOW4 | 4 | 2000 | Optimal | Optimal | 8.69e-01 | 1 | 11 | 40.5ms | 32.8ms | 0.8x | MISMATCH |
| SISSER | 2 | 0 | Acceptable | Optimal | 8.15e-11 | 15 | 18 | 19us | 2.2ms | 118.7x | PASS |
| SISSER2 | 2 | 0 | Acceptable | Optimal | 7.49e-11 | 16 | 20 | 22us | 2.5ms | 116.1x | PASS |
| SNAIL | 2 | 0 | Optimal | Optimal | 1.63e-26 | 63 | 63 | 80us | 9.2ms | 115.4x | PASS |
| SNAKE | 2 | 2 | Optimal | Optimal | 2.00e-04 | 5 | 8 | 20us | 1.8ms | 88.8x | MISMATCH |
| SPANHYD | 97 | 33 | MaxIteration | Optimal | N/A | 2999 | 20 | 8.98s | 5.2ms | 0.0x | ripopt_FAIL |
| SPIRAL | 3 | 2 | Acceptable | Infeasible | N/A | 116 | 370 | 285us | 58.9ms | 206.4x | ipopt_FAIL |
| SSI | 3 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 3.2ms | 433.7ms | 133.9x | ipopt_FAIL |
| SSINE | 3 | 2 | Acceptable | Optimal | 0.00e+00 | 131 | 224 | 395us | 33.3ms | 84.4x | PASS |
| STANCMIN | 3 | 2 | Optimal | Optimal | 9.32e-09 | 10 | 9 | 24us | 1.7ms | 72.9x | PASS |
| STRATEC | 10 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| STREG | 4 | 0 | Optimal | Optimal | 8.90e-02 | 21 | 13 | 29us | 2.2ms | 76.8x | MISMATCH |
| STREGNE | 4 | 2 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 12us | 597us | 49.7x | PASS |
| SUPERSIM | 2 | 2 | Optimal | Optimal | 2.22e-16 | 7 | 1 | 23us | 596us | 25.5x | PASS |
| SWOPF | 83 | 92 | Optimal | Optimal | 9.19e-10 | 15 | 13 | 1.2ms | 3.5ms | 2.9x | PASS |
| SYNTHES1 | 6 | 6 | Acceptable | Optimal | 1.08e-05 | 13 | 8 | 43us | 1.7ms | 39.8x | PASS |
| SYNTHES2 | 11 | 14 | Optimal | Optimal | 1.38e-06 | 22 | 14 | 161us | 2.8ms | 17.3x | PASS |
| SYNTHES3 | 17 | 23 | Optimal | Optimal | 3.99e-08 | 36 | 13 | 610us | 2.7ms | 4.5x | PASS |
| TAME | 2 | 1 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 15us | 1.1ms | 72.8x | PASS |
| TAX13322 | 72 | 1261 | Optimal | MaxIteration | N/A | 90 | 3000 | 230.3ms | 19.52s | 84.8x | ipopt_FAIL |
| TAXR13322 | 72 | 1261 | Acceptable | Acceptable | 9.00e-01 | 229 | 56 | 470.1ms | 2.68s | 5.7x | MISMATCH |
| TENBARS1 | 18 | 9 | Acceptable | Optimal | 9.86e-09 | 181 | 39 | 2.3ms | 6.8ms | 2.9x | PASS |
| TENBARS2 | 18 | 8 | Optimal | Optimal | 1.00e-08 | 27 | 33 | 261us | 6.1ms | 23.2x | PASS |
| TENBARS3 | 18 | 8 | Optimal | Optimal | 1.01e-08 | 25 | 34 | 251us | 6.3ms | 25.1x | PASS |
| TENBARS4 | 18 | 9 | Acceptable | Optimal | 1.51e-08 | 161 | 14 | 152.2ms | 3.1ms | 0.0x | PASS |
| TFI1 | 3 | 101 | Optimal | Optimal | 7.91e-10 | 307 | 19 | 22.2ms | 6.5ms | 0.3x | PASS |
| TFI2 | 3 | 101 | Optimal | Optimal | 3.26e-05 | 64 | 8 | 2.6ms | 2.6ms | 1.0x | PASS |
| TFI3 | 3 | 101 | Optimal | Optimal | 6.17e-09 | 79 | 13 | 3.0ms | 4.0ms | 1.3x | PASS |
| THURBER | 7 | 37 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 42.2ms | 278us | 0.0x | BOTH_FAIL |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 1.24e-14 | 30 | 19 | 308us | 3.3ms | 10.8x | PASS |
| TOINTGOR | 50 | 0 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 196us | 1.3ms | 6.4x | PASS |
| TOINTPSP | 50 | 0 | Optimal | Optimal | 0.00e+00 | 21 | 20 | 545us | 4.5ms | 8.2x | PASS |
| TOINTQOR | 50 | 0 | Optimal | Optimal | 1.93e-16 | 1 | 1 | 34us | 510us | 15.0x | PASS |
| TRIGGER | 7 | 6 | Optimal | Optimal | 0.00e+00 | 1 | 15 | 26us | 2.5ms | 98.6x | PASS |
| TRO3X3 | 30 | 13 | Acceptable | Optimal | 3.88e-03 | 236 | 47 | 11.3ms | 9.9ms | 0.9x | MISMATCH |
| TRO4X4 | 63 | 25 | Acceptable | IpoptStatus( | N/A | 2999 | 157 | 310.6ms | 40.3ms | 0.1x | ipopt_FAIL |
| TRO6X2 | 45 | 21 | Acceptable | RestorationF | N/A | 25 | 353 | 1.9ms | 90.6ms | 48.4x | ipopt_FAIL |
| TRUSPYR1 | 11 | 4 | Acceptable | Optimal | 2.41e-08 | 55 | 10 | 334us | 1.9ms | 5.7x | PASS |
| TRUSPYR2 | 11 | 11 | Optimal | Optimal | 8.13e-08 | 24 | 13 | 256us | 2.7ms | 10.5x | PASS |
| TRY-B | 2 | 1 | Optimal | Optimal | 5.01e-19 | 13 | 23 | 26us | 4.2ms | 162.6x | PASS |
| TWOBARS | 2 | 2 | Optimal | Optimal | 7.13e-01 | 19 | 8 | 54us | 1.6ms | 30.1x | MISMATCH |
| VESUVIA | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 1.57s | 462us | 0.0x | BOTH_FAIL |
| VESUVIALS | 8 | 0 | Acceptable | Optimal | 3.39e-01 | 330 | 48 | 72.7ms | 17.4ms | 0.2x | MISMATCH |
| VESUVIO | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 1.64s | 465us | 0.0x | BOTH_FAIL |
| VESUVIOLS | 8 | 0 | Acceptable | Optimal | 2.41e-15 | 2999 | 10 | 1.02s | 4.4ms | 0.0x | PASS |
| VESUVIOU | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 7.5ms | 435us | 0.1x | BOTH_FAIL |
| VESUVIOULS | 8 | 0 | Acceptable | Optimal | 5.55e-16 | 22 | 8 | 5.2ms | 3.2ms | 0.6x | PASS |
| VIBRBEAM | 8 | 0 | Optimal | Optimal | 9.48e-01 | 115 | 58 | 2.4ms | 9.0ms | 3.8x | MISMATCH |
| VIBRBEAMNE | 8 | 30 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 105.1ms | 282us | 0.0x | BOTH_FAIL |
| WACHBIEG | 3 | 2 | Optimal | Infeasible | N/A | 87 | 15 | 704us | 3.2ms | 4.5x | ipopt_FAIL |
| WATER | 31 | 10 | Acceptable | Optimal | 1.99e-07 | 37 | 17 | 920us | 3.2ms | 3.5x | PASS |
| WAYSEA1 | 2 | 0 | Optimal | Optimal | 2.69e-15 | 15 | 14 | 19us | 1.8ms | 93.0x | PASS |
| WAYSEA1B | 2 | 0 | Optimal | Optimal | 7.97e-09 | 13 | 14 | 17us | 2.4ms | 138.4x | PASS |
| WAYSEA1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 17us | 1.1ms | 62.1x | PASS |
| WAYSEA2 | 2 | 0 | Optimal | Optimal | 9.85e-18 | 23 | 22 | 26us | 2.6ms | 96.5x | PASS |
| WAYSEA2B | 2 | 0 | Optimal | Optimal | 1.61e-11 | 21 | 22 | 27us | 3.5ms | 129.6x | PASS |
| WAYSEA2NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 11 | 20us | 1.7ms | 84.1x | PASS |
| WEEDS | 3 | 0 | Optimal | Optimal | 3.43e-16 | 23 | 28 | 62us | 5.7ms | 92.1x | PASS |
| WEEDSNE | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 64us | 283us | 4.4x | BOTH_FAIL |
| WOMFLET | 3 | 3 | Acceptable | Optimal | 1.00e+00 | 14 | 8 | 38us | 1.6ms | 42.9x | MISMATCH |
| YFIT | 3 | 0 | Optimal | Optimal | 1.37e-19 | 36 | 36 | 104us | 6.5ms | 62.3x | PASS |
| YFITNE | 3 | 17 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 33us | 276us | 8.5x | ipopt_FAIL |
| YFITU | 3 | 0 | Optimal | Optimal | 6.69e-21 | 36 | 36 | 100us | 5.2ms | 51.9x | PASS |
| ZANGWIL2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 3us | 461us | 172.9x | PASS |
| ZANGWIL3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 9us | 457us | 51.8x | PASS |
| ZECEVIC2 | 2 | 2 | Optimal | Optimal | 5.57e-10 | 11 | 8 | 23us | 1.5ms | 67.5x | PASS |
| ZECEVIC3 | 2 | 2 | Optimal | Optimal | 8.24e-10 | 12 | 17 | 37us | 2.8ms | 76.5x | PASS |
| ZECEVIC4 | 2 | 2 | Optimal | Optimal | 2.57e-09 | 11 | 10 | 26us | 2.0ms | 79.6x | PASS |
| ZY2 | 3 | 2 | Acceptable | Optimal | 4.67e-05 | 13 | 14 | 32us | 2.8ms | 86.3x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 140.6 | 44.2 |
| Median | 13 | 12 |
| Total  | 74955 | 23561 |

- ripopt fewer iterations: 163/533
- Ipopt fewer iterations: 271/533
- Tied: 99/533

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 53.2ms | 23.7ms |
| Median | 59us | 2.5ms |
| Total  | 28.38s | 12.66s |

- Geometric mean speedup (Ipopt_time/ripopt_time): **21.30x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 489/533 problems
- Ipopt faster: 44/533 problems
- Overall speedup (total time): 0.45x

## Failure Analysis

### Problems where only ripopt fails (23)

| Problem | n | m | ripopt status | Ipopt obj |
|---------|---|---|---------------|-----------|
| AIRCRFTB | 8 | 0 | MaxIterations | 4.790634e-25 |
| ALLINIT | 4 | 0 | MaxIterations | 1.670597e+01 |
| ALLINITA | 4 | 4 | MaxIterations | 3.329611e+01 |
| BIGGS3 | 6 | 0 | MaxIterations | 3.031135e-21 |
| BIGGS5 | 6 | 0 | MaxIterations | 1.088200e-19 |
| BOX2 | 3 | 0 | MaxIterations | 6.251257e-23 |
| DECONVU | 63 | 0 | MaxIterations | 4.146188e-13 |
| DISCS | 36 | 66 | RestorationFailed | 1.200007e+01 |
| DNIEPER | 61 | 24 | MaxIterations | 1.874401e+04 |
| FEEDLOC | 90 | 259 | RestorationFailed | -9.539854e-10 |
| HATFLDH | 4 | 7 | MaxIterations | -2.450000e+01 |
| HEART6 | 6 | 6 | NumericalError | 0.000000e+00 |
| HS23 | 2 | 5 | MaxIterations | 2.000000e+00 |
| HS83 | 5 | 3 | MaxIterations | -3.066554e+04 |
| HS99EXP | 31 | 21 | RestorationFailed | -1.260006e+12 |
| MGH10SLS | 3 | 0 | MaxIterations | 8.794586e+01 |
| MINSURF | 64 | 0 | MaxIterations | 1.000000e+00 |
| NET1 | 48 | 57 | RestorationFailed | 9.411943e+05 |
| OET2 | 3 | 1002 | MaxIterations | 8.715962e-02 |
| OPTCNTRL | 32 | 20 | MaxIterations | 5.500000e+02 |
| OSBORNEA | 5 | 0 | MaxIterations | 5.464895e-05 |
| PFIT4 | 3 | 3 | NumericalError | 0.000000e+00 |
| SPANHYD | 97 | 33 | MaxIterations | 2.397380e+02 |

### Problems where only Ipopt fails (38)

| Problem | n | m | Ipopt status | ripopt obj |
|---------|---|---|--------------|------------|
| ARGAUSS | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| AVION2 | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CRESC100 | 6 | 200 | Infeasible | 7.443592e-01 |
| DECONVB | 63 | 0 | MaxIterations | 2.336812e-09 |
| DENSCHNBNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA1NE | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| DEVGLA2NE | 5 | 16 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | 9 | 3 | ErrorInStepComputation | -8.278941e+02 |
| EXP2NE | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| HS25NE | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS1 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS2 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS3 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 |
| LRCOVTYPE | 54 | 0 | N/A | 5.901541e-01 |
| MESH | 41 | 48 | IpoptStatus(4) | -2.405958e-03 |
| PALMER5E | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | 6 | 0 | MaxIterations | 1.033491e+01 |
| PFIT1 | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | 12 | 10 | MaxIterations | 7.207531e+00 |
| POLAK6 | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| S365 | 7 | 5 | RestorationFailed | 1.000000e-16 |
| S365MOD | 7 | 5 | RestorationFailed | 2.500000e-01 |
| SPIRAL | 3 | 2 | Infeasible | 2.322254e-12 |
| SSI | 3 | 0 | MaxIterations | 1.376194e-09 |
| TAX13322 | 72 | 1261 | MaxIterations | -4.123356e+03 |
| TRO4X4 | 63 | 25 | IpoptStatus(4) | 8.999898e+00 |
| TRO6X2 | 45 | 21 | RestorationFailed | 1.225000e+03 |
| WACHBIEG | 3 | 2 | Infeasible | 1.000000e+00 |
| YFITNE | 3 | 17 | IpoptStatus(-10) | 0.000000e+00 |

### Problems where both fail (133)

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
| CORE1 | 65 | 59 | Timeout | Timeout |
| CRESC132 | 6 | 2654 | Timeout | Timeout |
| DANIWOOD | 2 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| DANWOOD | 2 | 6 | LocalInfeasibility | IpoptStatus(-10) |
| DENSCHNENE | 3 | 3 | RestorationFailed | Infeasible |
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
| GULFNE | 3 | 99 | LocalInfeasibility | IpoptStatus(-10) |
| HAHN1 | 7 | 236 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDBNE | 4 | 4 | MaxIterations | Infeasible |
| HATFLDDNE | 3 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDENE | 3 | 21 | LocalInfeasibility | IpoptStatus(-10) |
| HIMMELBD | 2 | 2 | RestorationFailed | Infeasible |
| HIMMELBFNE | 4 | 7 | LocalInfeasibility | IpoptStatus(-10) |
| HIMMELBJ | 45 | 14 | RestorationFailed | ErrorInStepComputation |
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
| NYSTROM5 | 18 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| NYSTROM5C | 18 | 20 | LocalInfeasibility | IpoptStatus(-10) |
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
| SANTA | 21 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| STRATEC | 10 | 0 | Timeout | Timeout |
| THURBER | 7 | 37 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIA | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIO | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIOU | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VIBRBEAMNE | 8 | 30 | LocalInfeasibility | IpoptStatus(-10) |
| WEEDSNE | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |

### Objective mismatches (95)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 54
- **Convergence gap** (one Acceptable): 41
- **Better objective found by**: ripopt 18, Ipopt 77

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
| WOMFLET | 1.010458e-11 | 6.050000e+00 | 1.00e+00 | Acceptable | Optimal | ripopt |
| OET7 | 1.141238e+06 | 4.465915e-05 | 1.00e+00 | Acceptable | Optimal | ipopt |
| GIGOMEZ1 | -5.720219e-10 | -3.000000e+00 | 1.00e+00 | Optimal | Optimal | ipopt |
| OET3 | 1.528803e+07 | 4.505043e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| BENNETT5LS | 1.613684e+05 | 5.563289e-04 | 1.00e+00 | Optimal | Optimal | ipopt |
| HS268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| S268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| DANWOODLS | 1.039178e+02 | 4.317308e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| HALDMADS | 1.223712e-04 | 2.218282e+00 | 1.00e+00 | Optimal | Optimal | ripopt |
| OET4 | 2.394093e+00 | 4.295421e-03 | 9.98e-01 | Acceptable | Optimal | ipopt |
| HS16 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 | Optimal | Optimal | ipopt |
| ELATTAR | 1.054115e+00 | 7.420618e+01 | 9.86e-01 | Optimal | Optimal | ripopt |
| KSIP | 2.768715e+01 | 5.757979e-01 | 9.79e-01 | Optimal | Optimal | ipopt |
| RAT43LS | 2.525083e+05 | 8.786405e+03 | 9.65e-01 | Acceptable | Optimal | ipopt |
| CHWIRUT1LS | 6.551998e+04 | 2.384477e+03 | 9.64e-01 | Acceptable | Optimal | ipopt |
| CHWIRUT2LS | 1.355314e+04 | 5.130480e+02 | 9.62e-01 | Acceptable | Optimal | ipopt |
| HS35MOD | 6.500000e+00 | 2.500000e-01 | 9.62e-01 | Optimal | Optimal | ipopt |
| MWRIGHT | 1.288383e+00 | 2.497881e+01 | 9.48e-01 | Optimal | Optimal | ripopt |
| VIBRBEAM | 6.346254e+00 | 3.322376e-01 | 9.48e-01 | Optimal | Optimal | ipopt |
| BT4 | -4.551055e+01 | -3.704768e+00 | 9.19e-01 | Optimal | Optimal | ripopt |
| TAXR13322 | -3.420730e+03 | -3.429089e+02 | 9.00e-01 | Acceptable | Acceptable | ripopt |
| BIGGSC4 | -3.128034e+00 | -2.450000e+01 | 8.72e-01 | Acceptable | Optimal | ipopt |
| SIPOW4 | 2.080345e+00 | 2.723620e-01 | 8.69e-01 | Optimal | Optimal | ipopt |
| HIMMELP2 | -6.205394e+01 | -8.198044e+00 | 8.68e-01 | Optimal | Optimal | ripopt |
| HIMMELP3 | -7.913699e+00 | -5.901318e+01 | 8.66e-01 | Optimal | Optimal | ipopt |
| CAMEL6 | -2.154638e-01 | -1.031628e+00 | 7.91e-01 | Optimal | Optimal | ipopt |
| HS54 | -1.566691e-01 | -9.080748e-01 | 7.51e-01 | Optimal | Optimal | ipopt |
| HIMMELP6 | -1.475339e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| HIMMELP5 | -1.475901e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| TWOBARS | 5.257563e+00 | 1.508652e+00 | 7.13e-01 | Optimal | Optimal | ipopt |
| SIPOW3 | 1.846738e+00 | 5.346586e-01 | 7.10e-01 | Optimal | Optimal | ipopt |
| MAKELA4 | 6.877365e-01 | -9.600000e-09 | 6.88e-01 | Optimal | Optimal | ipopt |
| MISTAKE | -5.000000e-01 | -1.000000e+00 | 5.00e-01 | Acceptable | Optimal | ipopt |
| EGGCRATE | 1.897639e+01 | 9.488197e+00 | 5.00e-01 | Optimal | Optimal | ipopt |
| ECKERLE4LS | 4.988568e-01 | 1.463589e-03 | 4.97e-01 | Acceptable | Optimal | ipopt |
| FLETCHER | 1.952537e+01 | 1.165685e+01 | 4.03e-01 | Optimal | Optimal | ipopt |
| EXPFITA | 3.932166e-01 | 1.136646e-03 | 3.92e-01 | Optimal | Optimal | ipopt |
| VESUVIALS | 1.500440e+03 | 9.914100e+02 | 3.39e-01 | Acceptable | Optimal | ipopt |
| HYDC20LS | 2.976453e-01 | 2.967522e-15 | 2.98e-01 | Acceptable | Optimal | ipopt |
| AVGASA | -3.383299e+00 | -4.631926e+00 | 2.70e-01 | Acceptable | Optimal | ipopt |
| EG1 | -1.132801e+00 | -1.429307e+00 | 2.07e-01 | Optimal | Optimal | ipopt |
| HS70 | 1.877383e-01 | 7.498464e-03 | 1.80e-01 | Optimal | Optimal | ipopt |
| AVGASB | -3.717055e+00 | -4.483219e+00 | 1.71e-01 | Acceptable | Optimal | ipopt |
| BT7 | 3.603798e+02 | 3.065000e+02 | 1.50e-01 | Optimal | Optimal | ipopt |
| MUONSINELS | 5.080769e+04 | 4.387412e+04 | 1.36e-01 | Acceptable | Optimal | ipopt |
| MAKELA2 | 8.244898e+00 | 7.200000e+00 | 1.27e-01 | Optimal | Optimal | ipopt |
| MSS1 | -1.600000e+01 | -1.400000e+01 | 1.25e-01 | Acceptable | Optimal | ripopt |
| OET1 | 6.431699e-01 | 5.382431e-01 | 1.05e-01 | Optimal | Optimal | ipopt |
| CRESC50 | 8.827223e-01 | 7.862467e-01 | 9.65e-02 | Acceptable | Optimal | ipopt |
| STREG | 3.743976e-21 | 8.901950e-02 | 8.90e-02 | Optimal | Optimal | ripopt |
| HAHN1LS | 3.086398e+01 | 3.338424e+01 | 7.55e-02 | Acceptable | Optimal | ripopt |
| PALMER3 | 2.265958e+03 | 2.416980e+03 | 6.25e-02 | Acceptable | Optimal | ripopt |
| PALMER4 | 2.424016e+03 | 2.285383e+03 | 5.72e-02 | Optimal | Optimal | ipopt |
| MGH17SLS | 5.464895e-05 | 2.451788e-02 | 2.45e-02 | Optimal | Optimal | ripopt |
| CB2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CHACONN1 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| HS55 | 6.666667e+00 | 6.805833e+00 | 2.04e-02 | Optimal | Optimal | ripopt |
| QC | -9.411446e+02 | -9.565379e+02 | 1.61e-02 | Acceptable | Optimal | ipopt |
| DUALC8 | 1.854450e+04 | 1.830936e+04 | 1.27e-02 | Acceptable | Optimal | ipopt |
| HET-Z | 1.011865e+00 | 1.000000e+00 | 1.17e-02 | Optimal | Optimal | ipopt |
| MAXLIKA | 1.149346e+03 | 1.136307e+03 | 1.13e-02 | Acceptable | Optimal | ipopt |
| PT | 1.873320e-01 | 1.783942e-01 | 8.94e-03 | Optimal | Optimal | ipopt |
| EXPFITC | 3.111184e-02 | 2.330262e-02 | 7.81e-03 | Optimal | Optimal | ipopt |
| CLIFF | 1.997866e-01 | 2.072380e-01 | 7.45e-03 | Optimal | Optimal | ripopt |
| ACOPR30 | 5.807813e+02 | 5.768924e+02 | 6.70e-03 | Optimal | Optimal | ipopt |
| DGOSPEC | -9.887540e+02 | -9.933506e+02 | 4.63e-03 | Acceptable | Optimal | ipopt |
| TRO3X3 | 9.002446e+00 | 8.967478e+00 | 3.88e-03 | Acceptable | Optimal | ipopt |
| ACOPR14 | 8.107633e+03 | 8.081526e+03 | 3.22e-03 | Acceptable | Optimal | ipopt |
| DECONVC | 3.880899e-10 | 2.569475e-03 | 2.57e-03 | Acceptable | Optimal | ripopt |
| QCNEW | -8.045191e+02 | -8.065219e+02 | 2.48e-03 | Optimal | Optimal | ipopt |
| DEGENLPB | -3.069162e+01 | -3.076401e+01 | 2.35e-03 | Acceptable | Optimal | ipopt |
| HS85 | -2.211026e+00 | -2.215605e+00 | 2.07e-03 | Optimal | Optimal | ipopt |
| LIN | -1.960628e-02 | -1.757754e-02 | 2.03e-03 | Optimal | Optimal | ripopt |
| DEGENLPA | 3.060434e+00 | 3.054881e+00 | 1.81e-03 | Acceptable | Optimal | ipopt |
| HS59 | -6.743243e+00 | -6.749505e+00 | 9.28e-04 | Acceptable | Optimal | ipopt |
| HS116 | 9.766220e+01 | 9.758747e+01 | 7.65e-04 | Acceptable | Optimal | ipopt |
| HS13 | 9.938594e-01 | 9.945785e-01 | 7.19e-04 | Acceptable | Optimal | ripopt |
| HS95 | 1.630984e-02 | 1.561772e-02 | 6.92e-04 | Acceptable | Optimal | ipopt |
| HS96 | 1.630574e-02 | 1.561775e-02 | 6.88e-04 | Acceptable | Optimal | ipopt |
| HS45 | 1.000314e+00 | 1.000000e+00 | 3.14e-04 | Acceptable | Optimal | ipopt |
| LSC2LS | 1.333749e+01 | 1.333439e+01 | 2.32e-04 | Acceptable | Optimal | ipopt |
| DENSCHND | 6.940852e-09 | 2.221899e-04 | 2.22e-04 | Acceptable | Optimal | ripopt |
| HS17 | 1.000201e+00 | 1.000000e+00 | 2.01e-04 | Optimal | Optimal | ipopt |
| SNAKE | -2.601155e-10 | -1.999999e-04 | 2.00e-04 | Optimal | Optimal | ipopt |
| HS98 | 3.136227e+00 | 3.135806e+00 | 1.34e-04 | Acceptable | Optimal | ipopt |
| HS119 | 2.449254e+02 | 2.448997e+02 | 1.05e-04 | Acceptable | Optimal | ipopt |

---
*Generated by cutest_suite/compare.py*