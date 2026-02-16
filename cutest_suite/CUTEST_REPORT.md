# CUTEst Benchmark Report

Comparison of ripopt vs Ipopt (C++) on the CUTEst test set.

## Executive Summary

- **Total problems**: 727
- **ripopt solved**: 597/727 (82.1%)
- **Ipopt solved**: 555/727 (76.3%)
- **Both solved**: 554/727
- **Matching solutions** (rel obj diff < 1e-4): 450/554

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (450 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 1.45e-06 |
| Median | 3.80e-12 |
| Max    | 9.46e-05 |

**All both-solved** (554 problems, including 104 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 7.92e-02 |
| Median | 2.21e-09 |
| Max    | 2.00e+00 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 493 | 376 | 340 | 339 | 274 |
| unconstrained | 234 | 221 | 215 | 215 | 176 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| 3PK | 30 | 0 | Optimal | Optimal | 1.42e-15 | 9 | 9 | 5.6ms | 2.0ms | 0.4x | PASS |
| ACOPP14 | 38 | 68 | Optimal | Optimal | 9.79e-10 | 16 | 9 | 5.8ms | 3.4ms | 0.6x | PASS |
| ACOPP30 | 72 | 142 | Optimal | Optimal | 6.40e-09 | 47 | 13 | 9.9ms | 6.1ms | 0.6x | PASS |
| ACOPR14 | 38 | 82 | Acceptable | Optimal | 8.78e-02 | 2320 | 13 | 1.14s | 4.8ms | 0.0x | MISMATCH |
| ACOPR30 | 72 | 172 | Acceptable | Optimal | 2.00e-07 | 557 | 221 | 1.78s | 116.9ms | 0.1x | PASS |
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 22us | 772us | 35.8x | PASS |
| AIRCRFTB | 8 | 0 | Acceptable | Optimal | 1.08e-16 | 687 | 15 | 8.3ms | 2.6ms | 0.3x | PASS |
| AIRPORT | 84 | 42 | Optimal | Optimal | 4.45e-09 | 14 | 13 | 4.1ms | 5.7ms | 1.4x | PASS |
| AKIVA | 2 | 0 | Optimal | Optimal | 2.88e-16 | 14 | 6 | 55us | 1.1ms | 19.4x | PASS |
| ALLINIT | 4 | 0 | Acceptable | Optimal | 3.69e-09 | 10 | 20 | 43.7ms | 3.5ms | 0.1x | PASS |
| ALLINITA | 4 | 4 | Acceptable | Optimal | 9.53e-06 | 27 | 12 | 101us | 2.4ms | 23.8x | PASS |
| ALLINITC | 4 | 1 | Acceptable | Optimal | 4.02e-07 | 26 | 17 | 88us | 2.9ms | 33.1x | PASS |
| ALLINITU | 4 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 14 | 11us | 2.3ms | 217.4x | PASS |
| ALSOTAME | 2 | 1 | Optimal | Optimal | 1.47e-08 | 10 | 8 | 23us | 1.6ms | 70.6x | PASS |
| ANTWERP | 27 | 10 | Acceptable | Optimal | 4.67e-03 | 131 | 108 | 4.5ms | 22.8ms | 5.1x | MISMATCH |
| ARGAUSS | 3 | 15 | Acceptable | IpoptStatus( | N/A | 2 | 0 | 21us | 292us | 13.7x | ipopt_FAIL |
| AVGASA | 8 | 10 | Acceptable | Optimal | 2.70e-01 | 11 | 9 | 59us | 2.2ms | 36.7x | MISMATCH |
| AVGASB | 8 | 10 | Acceptable | Optimal | 1.71e-01 | 19 | 11 | 105us | 2.5ms | 24.2x | MISMATCH |
| AVION2 | 49 | 15 | Acceptable | MaxIteration | N/A | 26 | 3000 | 3.2ms | 670.6ms | 207.3x | ipopt_FAIL |
| BA-L1 | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 6 | 558us | 1.7ms | 3.1x | PASS |
| BA-L1LS | 57 | 0 | Optimal | Optimal | 7.62e-21 | 28 | 10 | 229us | 2.3ms | 10.2x | PASS |
| BA-L1SP | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.3ms | 2.4ms | 1.9x | PASS |
| BA-L1SPLS | 57 | 0 | Optimal | Optimal | 6.48e-17 | 32 | 9 | 1.2ms | 4.2ms | 3.6x | PASS |
| BARD | 3 | 0 | Optimal | Optimal | 1.73e-17 | 24 | 8 | 20us | 1.3ms | 68.2x | PASS |
| BARDNE | 3 | 15 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 75us | 301us | 4.0x | BOTH_FAIL |
| BATCH | 48 | 73 | Acceptable | Optimal | 2.29e-05 | 201 | 29 | 9.5ms | 7.2ms | 0.8x | PASS |
| BEALE | 2 | 0 | Optimal | Optimal | 4.26e-18 | 19 | 8 | 13us | 1.7ms | 134.0x | PASS |
| BEALENE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 21us | 281us | 13.4x | ipopt_FAIL |
| BENNETT5 | 3 | 154 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 27.2ms | 329us | 0.0x | BOTH_FAIL |
| BENNETT5LS | 3 | 0 | Acceptable | Optimal | 1.74e-05 | 40 | 21 | 2.4ms | 3.8ms | 1.6x | PASS |
| BIGGS3 | 6 | 0 | Optimal | Optimal | 9.03e-21 | 41 | 9 | 399us | 1.9ms | 4.7x | PASS |
| BIGGS5 | 6 | 0 | Acceptable | Optimal | 5.66e-03 | 125 | 20 | 524us | 3.5ms | 6.6x | MISMATCH |
| BIGGS6 | 6 | 0 | Optimal | Optimal | 5.66e-03 | 41 | 79 | 46us | 12.2ms | 264.4x | MISMATCH |
| BIGGS6NE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 3.9ms | 279us | 0.1x | BOTH_FAIL |
| BIGGSC4 | 4 | 7 | Acceptable | Optimal | 8.72e-01 | 2999 | 17 | 9.1ms | 3.1ms | 0.3x | MISMATCH |
| BLEACHNG | 17 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| BOOTH | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 6us | 453us | 74.4x | PASS |
| BOX2 | 3 | 0 | Acceptable | Optimal | 2.86e-09 | 677 | 8 | 12.5ms | 1.3ms | 0.1x | PASS |
| BOX3 | 3 | 0 | Optimal | Optimal | 3.53e-19 | 23 | 9 | 20us | 1.6ms | 80.5x | PASS |
| BOX3NE | 3 | 10 | Optimal | IpoptStatus( | N/A | 11 | 0 | 50us | 288us | 5.8x | ipopt_FAIL |
| BOXBOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 49us | 292us | 5.9x | BOTH_FAIL |
| BOXBODLS | 2 | 0 | Acceptable | Optimal | 8.80e-01 | 27 | 13 | 103us | 2.3ms | 22.7x | MISMATCH |
| BQP1VAR | 1 | 0 | Optimal | Optimal | 9.98e-09 | 7 | 5 | 21us | 1.2ms | 56.1x | PASS |
| BQPGABIM | 50 | 0 | Acceptable | Optimal | 4.40e-06 | 11 | 12 | 4.8ms | 2.5ms | 0.5x | PASS |
| BQPGASIM | 50 | 0 | Acceptable | Optimal | 6.13e-06 | 11 | 12 | 5.1ms | 2.4ms | 0.5x | PASS |
| BRANIN | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 7 | 6us | 1.6ms | 286.9x | PASS |
| BRKMCC | 2 | 0 | Acceptable | Optimal | 8.33e-17 | 8 | 3 | 31us | 680us | 21.9x | PASS |
| BROWNBS | 2 | 0 | Acceptable | Optimal | 0.00e+00 | 12 | 7 | 7us | 1.6ms | 225.4x | PASS |
| BROWNBSNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 8 | 0 | 20us | 378us | 18.6x | ipopt_FAIL |
| BROWNDEN | 4 | 0 | Optimal | Optimal | 8.48e-16 | 23 | 8 | 30us | 1.3ms | 41.8x | PASS |
| BROWNDENE | 4 | 20 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 134us | 284us | 2.1x | BOTH_FAIL |
| BT1 | 2 | 1 | Optimal | Optimal | 2.40e-09 | 17 | 7 | 53us | 1.3ms | 25.1x | PASS |
| BT10 | 2 | 2 | Optimal | Optimal | 2.79e-09 | 7 | 6 | 17us | 1.1ms | 63.1x | PASS |
| BT11 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 8 | 29us | 1.2ms | 43.6x | PASS |
| BT12 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 16us | 1.0ms | 60.9x | PASS |
| BT13 | 5 | 1 | Acceptable | Optimal | 1.00e-08 | 27 | 24 | 59us | 4.2ms | 70.1x | PASS |
| BT2 | 3 | 1 | Optimal | Optimal | 3.80e-12 | 11 | 12 | 26us | 1.7ms | 66.7x | PASS |
| BT3 | 5 | 3 | Optimal | Optimal | 1.22e-14 | 1 | 1 | 10us | 544us | 56.0x | PASS |
| BT4 | 3 | 2 | Optimal | Optimal | 9.19e-01 | 6 | 9 | 24us | 1.7ms | 68.6x | MISMATCH |
| BT5 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 18us | 1.2ms | 64.3x | PASS |
| BT6 | 5 | 2 | Optimal | Optimal | 3.46e-12 | 9 | 13 | 33us | 1.9ms | 57.9x | PASS |
| BT7 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 16 | 68us | 2.9ms | 42.9x | PASS |
| BT8 | 5 | 2 | Acceptable | Optimal | 3.73e-09 | 32 | 14 | 86us | 2.3ms | 27.0x | PASS |
| BT9 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 35us | 1.9ms | 54.1x | PASS |
| BURKEHAN | 1 | 1 | RestorationF | Infeasible | N/A | 229 | 11 | 154.8ms | 2.6ms | 0.0x | BOTH_FAIL |
| BYRDSPHR | 3 | 2 | Optimal | Optimal | 3.95e-10 | 23 | 12 | 85us | 2.1ms | 24.4x | PASS |
| CAMEL6 | 2 | 0 | Optimal | Optimal | 4.30e-16 | 12 | 8 | 7us | 1.8ms | 253.5x | PASS |
| CANTILVR | 5 | 1 | Optimal | Optimal | 3.33e-09 | 31 | 11 | 67us | 2.1ms | 32.0x | PASS |
| CB2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 8 | 8 | 22us | 1.8ms | 80.2x | MISMATCH |
| CB3 | 3 | 3 | Optimal | Optimal | 4.30e-09 | 8 | 8 | 21us | 1.6ms | 78.6x | PASS |
| CERI651A | 7 | 61 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 660.5ms | 328us | 0.0x | BOTH_FAIL |
| CERI651ALS | 7 | 0 | Acceptable | Optimal | 8.56e-08 | 283 | 95 | 352.0ms | 15.6ms | 0.0x | PASS |
| CERI651B | 7 | 66 | LocalInfeasi | IpoptStatus( | N/A | 83 | 0 | 10.5ms | 302us | 0.0x | BOTH_FAIL |
| CERI651BLS | 7 | 0 | Optimal | Optimal | 1.23e-08 | 89 | 56 | 372.4ms | 8.6ms | 0.0x | PASS |
| CERI651C | 7 | 56 | LocalInfeasi | IpoptStatus( | N/A | 471 | 0 | 10.8ms | 279us | 0.0x | BOTH_FAIL |
| CERI651CLS | 7 | 0 | Acceptable | Optimal | 7.27e-08 | 482 | 53 | 2.6ms | 7.5ms | 2.9x | PASS |
| CERI651D | 7 | 67 | LocalInfeasi | IpoptStatus( | N/A | 91 | 0 | 4.3ms | 304us | 0.1x | BOTH_FAIL |
| CERI651DLS | 7 | 0 | Acceptable | Optimal | 3.08e-10 | 98 | 60 | 2.8ms | 10.2ms | 3.7x | PASS |
| CERI651E | 7 | 64 | LocalInfeasi | IpoptStatus( | N/A | 52 | 0 | 4.2ms | 301us | 0.1x | BOTH_FAIL |
| CERI651ELS | 7 | 0 | Acceptable | Optimal | 1.69e-07 | 166 | 45 | 1.9ms | 6.3ms | 3.4x | PASS |
| CHACONN1 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 5 | 6 | 17us | 1.4ms | 85.5x | MISMATCH |
| CHACONN2 | 3 | 3 | Optimal | Optimal | 4.49e-09 | 6 | 6 | 19us | 1.3ms | 71.3x | PASS |
| CHWIRUT1 | 3 | 214 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 3.0ms | 391us | 0.1x | BOTH_FAIL |
| CHWIRUT1LS | 3 | 0 | Optimal | Optimal | 7.63e-16 | 30 | 6 | 250us | 1.4ms | 5.7x | PASS |
| CHWIRUT2 | 3 | 54 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 500us | 284us | 0.6x | BOTH_FAIL |
| CHWIRUT2LS | 3 | 0 | Acceptable | Optimal | 8.86e-16 | 22 | 6 | 178us | 1.3ms | 7.5x | PASS |
| CLIFF | 2 | 0 | Optimal | Optimal | 7.45e-03 | 10 | 23 | 14us | 3.0ms | 209.0x | MISMATCH |
| CLUSTER | 2 | 2 | Acceptable | Optimal | 0.00e+00 | 12 | 9 | 29us | 1.7ms | 58.7x | PASS |
| CLUSTERLS | 2 | 0 | Acceptable | Optimal | 1.96e-14 | 19 | 17 | 10us | 2.5ms | 240.7x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 238us | 1.5ms | 6.5x | PASS |
| CONGIGMZ | 3 | 5 | Optimal | Optimal | 3.19e-09 | 7 | 20 | 68us | 3.4ms | 50.4x | PASS |
| COOLHANS | 9 | 9 | Optimal | Optimal | 0.00e+00 | 22 | 9 | 164us | 1.8ms | 10.7x | PASS |
| COOLHANSLS | 9 | 0 | Acceptable | Optimal | 2.25e-08 | 173 | 25 | 183us | 3.9ms | 21.2x | PASS |
| CORE1 | 65 | 59 | Optimal | Optimal | 4.04e-09 | 24 | 33 | 52.7ms | 8.5ms | 0.2x | PASS |
| CRESC100 | 6 | 200 | Acceptable | Infeasible | N/A | 169 | 155 | 655.4ms | 116.8ms | 0.2x | ipopt_FAIL |
| CRESC132 | 6 | 2654 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC4 | 6 | 8 | Optimal | Optimal | 1.35e-08 | 21 | 64 | 270us | 11.7ms | 43.2x | PASS |
| CRESC50 | 6 | 100 | Optimal | Optimal | 3.57e-08 | 153 | 194 | 200.5ms | 79.2ms | 0.4x | PASS |
| CSFI1 | 5 | 4 | Optimal | Optimal | 1.49e-08 | 18 | 11 | 54us | 2.2ms | 41.1x | PASS |
| CSFI2 | 5 | 4 | Optimal | Optimal | 1.50e-08 | 16 | 14 | 137us | 2.8ms | 20.3x | PASS |
| CUBE | 2 | 0 | Optimal | Optimal | 4.92e-25 | 34 | 27 | 18us | 4.0ms | 225.9x | PASS |
| CUBENE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 27 | 1 | 43us | 481us | 11.3x | PASS |
| DALLASS | 46 | 31 | Optimal | Optimal | 1.41e-07 | 29 | 22 | 3.5ms | 4.6ms | 1.3x | PASS |
| DANIWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 36us | 275us | 7.6x | BOTH_FAIL |
| DANIWOODLS | 2 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 10 | 1us | 1.7ms | 1387.6x | MISMATCH |
| DANWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 36us | 289us | 8.1x | BOTH_FAIL |
| DANWOODLS | 2 | 0 | Optimal | Optimal | 1.28e-16 | 19 | 11 | 19us | 1.8ms | 95.4x | PASS |
| DECONVB | 63 | 0 | Acceptable | MaxIteration | N/A | 306 | 3000 | 26.7ms | 699.2ms | 26.2x | ipopt_FAIL |
| DECONVBNE | 63 | 40 | Acceptable | Optimal | 0.00e+00 | 49 | 505 | 24.8ms | 157.6ms | 6.3x | PASS |
| DECONVC | 63 | 1 | Acceptable | Optimal | 1.39e-03 | 38 | 31 | 5.6ms | 8.3ms | 1.5x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 2 | 26 | 652us | 23.5ms | 36.0x | PASS |
| DECONVU | 63 | 0 | Acceptable | Optimal | 2.37e-07 | 100 | 333 | 502us | 83.7ms | 166.6x | PASS |
| DEGENLPA | 20 | 15 | Acceptable | Optimal | 1.81e-03 | 40 | 18 | 656us | 3.4ms | 5.3x | MISMATCH |
| DEGENLPB | 20 | 15 | Acceptable | Optimal | 2.35e-03 | 33 | 19 | 467us | 3.4ms | 7.4x | MISMATCH |
| DEMBO7 | 16 | 20 | Acceptable | Optimal | 7.72e-06 | 239 | 45 | 7.6ms | 8.3ms | 1.1x | PASS |
| DEMYMALO | 3 | 3 | Optimal | Optimal | 2.96e-09 | 12 | 9 | 26us | 1.9ms | 72.3x | PASS |
| DENSCHNA | 2 | 0 | Optimal | Optimal | 9.98e-19 | 10 | 6 | 6us | 962us | 166.1x | PASS |
| DENSCHNB | 2 | 0 | Optimal | Optimal | 1.88e-19 | 10 | 7 | 5us | 1.3ms | 272.7x | PASS |
| DENSCHNBNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 19us | 286us | 15.2x | ipopt_FAIL |
| DENSCHNC | 2 | 0 | Optimal | Optimal | 1.83e-01 | 17 | 10 | 12us | 1.4ms | 120.6x | MISMATCH |
| DENSCHNCNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 10 | 7 | 20us | 1.2ms | 60.4x | PASS |
| DENSCHND | 3 | 0 | Acceptable | Optimal | 2.22e-04 | 53 | 26 | 30us | 3.6ms | 119.5x | MISMATCH |
| DENSCHNDNE | 3 | 3 | Acceptable | Acceptable | 0.00e+00 | 19 | 22 | 103us | 3.3ms | 31.6x | PASS |
| DENSCHNE | 3 | 0 | Optimal | Optimal | 1.35e-17 | 28 | 14 | 15us | 2.5ms | 167.8x | PASS |
| DENSCHNENE | 3 | 3 | Optimal | Infeasible | N/A | 10 | 10 | 24us | 2.0ms | 84.2x | ipopt_FAIL |
| DENSCHNF | 2 | 0 | Optimal | Optimal | 4.55e-22 | 12 | 6 | 7us | 961us | 134.0x | PASS |
| DENSCHNFNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 16us | 957us | 60.3x | PASS |
| DEVGLA1 | 4 | 0 | Optimal | Optimal | 1.40e-22 | 46 | 23 | 121us | 3.6ms | 29.9x | PASS |
| DEVGLA1B | 4 | 0 | Optimal | Optimal | 1.11e-24 | 67 | 20 | 687us | 4.3ms | 6.3x | PASS |
| DEVGLA1NE | 4 | 24 | Optimal | IpoptStatus( | N/A | 16 | 0 | 191us | 279us | 1.5x | ipopt_FAIL |
| DEVGLA2 | 5 | 0 | Acceptable | Optimal | 2.50e-02 | 74 | 13 | 796us | 2.0ms | 2.6x | MISMATCH |
| DEVGLA2B | 5 | 0 | Acceptable | Optimal | 2.58e-07 | 14 | 24 | 246.3ms | 4.6ms | 0.0x | PASS |
| DEVGLA2NE | 5 | 16 | LocalInfeasi | IpoptStatus( | N/A | 101 | 0 | 718us | 289us | 0.4x | BOTH_FAIL |
| DGOSPEC | 3 | 0 | Acceptable | Optimal | 4.63e-03 | 10 | 27 | 22.9ms | 4.9ms | 0.2x | MISMATCH |
| DIAMON2D | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3D | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIPIGRI | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 66us | 2.0ms | 30.7x | PASS |
| DISC2 | 29 | 23 | Optimal | Optimal | 9.11e-10 | 30 | 24 | 3.2ms | 5.5ms | 1.7x | PASS |
| DISCS | 36 | 66 | Acceptable | Optimal | 2.30e-06 | 98 | 184 | 73.0ms | 70.4ms | 1.0x | PASS |
| DIXCHLNG | 10 | 5 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 83us | 1.6ms | 19.2x | PASS |
| DJTL | 2 | 0 | Acceptable | Acceptable | 1.22e-15 | 1959 | 1538 | 1.7ms | 159.0ms | 95.9x | PASS |
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
| DNIEPER | 61 | 24 | Acceptable | Optimal | 4.28e-08 | 1426 | 23 | 226.1ms | 5.0ms | 0.0x | PASS |
| DUAL1 | 85 | 1 | Acceptable | Optimal | 3.88e-06 | 12 | 15 | 2.9ms | 6.2ms | 2.1x | PASS |
| DUAL2 | 96 | 1 | Acceptable | Optimal | 8.07e-09 | 12 | 12 | 3.8ms | 5.8ms | 1.5x | PASS |
| DUAL4 | 75 | 1 | Acceptable | Optimal | 1.40e-07 | 11 | 12 | 2.0ms | 4.5ms | 2.2x | PASS |
| DUALC1 | 9 | 215 | Acceptable | Optimal | 6.76e-06 | 47 | 18 | 4.4ms | 10.9ms | 2.5x | PASS |
| DUALC2 | 7 | 229 | Acceptable | Optimal | 1.62e-06 | 17 | 12 | 2.2ms | 7.4ms | 3.4x | PASS |
| DUALC5 | 8 | 278 | Acceptable | Optimal | 1.83e-07 | 12 | 11 | 2.3ms | 8.4ms | 3.7x | PASS |
| DUALC8 | 8 | 503 | Acceptable | Optimal | 1.77e-07 | 24 | 13 | 14.6ms | 14.6ms | 1.0x | PASS |
| ECKERLE4 | 3 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 222us | 290us | 1.3x | BOTH_FAIL |
| ECKERLE4LS | 3 | 0 | Acceptable | Optimal | 6.98e-01 | 12 | 36 | 21us | 5.8ms | 274.4x | MISMATCH |
| EG1 | 3 | 0 | Optimal | Optimal | 2.07e-01 | 9 | 8 | 29.9ms | 1.8ms | 0.1x | MISMATCH |
| EGGCRATE | 2 | 0 | Optimal | Optimal | 5.62e-16 | 10 | 5 | 7us | 1.0ms | 153.7x | PASS |
| EGGCRATEB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 6 | 18us | 1.3ms | 72.1x | PASS |
| EGGCRATENE | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 36us | 287us | 7.9x | BOTH_FAIL |
| ELATTAR | 7 | 102 | Acceptable | Optimal | 1.00e+00 | 99 | 81 | 153.2ms | 34.3ms | 0.2x | MISMATCH |
| ELATVIDU | 2 | 0 | Optimal | Optimal | 9.69e-01 | 14 | 11 | 10us | 1.5ms | 158.9x | MISMATCH |
| ELATVIDUB | 2 | 0 | Acceptable | Optimal | 9.69e-01 | 18 | 11 | 34us | 2.0ms | 58.7x | MISMATCH |
| ELATVIDUNE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 43us | 283us | 6.7x | BOTH_FAIL |
| ENGVAL2 | 3 | 0 | Optimal | Optimal | 1.70e-20 | 31 | 21 | 18us | 3.0ms | 168.2x | PASS |
| ENGVAL2NE | 3 | 5 | Optimal | IpoptStatus( | N/A | 17 | 0 | 46us | 286us | 6.2x | ipopt_FAIL |
| ENSO | 9 | 168 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 8.2ms | 325us | 0.0x | BOTH_FAIL |
| ENSOLS | 9 | 0 | Acceptable | Optimal | 4.33e-16 | 27 | 7 | 1.8ms | 2.1ms | 1.2x | PASS |
| EQC | 9 | 3 | Acceptable | ErrorInStepC | N/A | 7 | 15 | 96us | 4.2ms | 43.4x | ipopt_FAIL |
| ERRINBAR | 18 | 9 | Acceptable | Optimal | 4.70e-07 | 68 | 37 | 728us | 6.7ms | 9.2x | PASS |
| EXP2 | 2 | 0 | Optimal | Optimal | 1.58e-20 | 13 | 7 | 11us | 1.3ms | 114.6x | PASS |
| EXP2B | 2 | 0 | Optimal | Optimal | 3.41e-21 | 13 | 7 | 9us | 1.4ms | 153.0x | PASS |
| EXP2NE | 2 | 10 | Optimal | IpoptStatus( | N/A | 7 | 0 | 30us | 283us | 9.3x | ipopt_FAIL |
| EXPFIT | 2 | 0 | Optimal | Optimal | 1.39e-16 | 14 | 8 | 15us | 1.3ms | 91.7x | PASS |
| EXPFITA | 5 | 22 | Optimal | Optimal | 3.40e-08 | 40 | 13 | 237us | 2.8ms | 11.9x | PASS |
| EXPFITB | 5 | 102 | Optimal | Optimal | 2.00e-03 | 189 | 16 | 7.5ms | 5.5ms | 0.7x | MISMATCH |
| EXPFITC | 5 | 502 | Optimal | Optimal | 9.42e-03 | 17 | 18 | 6.5ms | 16.6ms | 2.5x | MISMATCH |
| EXPFITNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 46us | 284us | 6.2x | BOTH_FAIL |
| EXTRASIM | 2 | 1 | Optimal | Optimal | 1.82e-08 | 4 | 3 | 11us | 878us | 77.2x | PASS |
| FBRAIN | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 11.0ms | 501us | 0.0x | BOTH_FAIL |
| FBRAIN2 | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 149.4ms | 678us | 0.0x | BOTH_FAIL |
| FBRAIN2LS | 4 | 0 | Optimal | Optimal | 5.11e-10 | 10 | 10 | 191.7ms | 8.7ms | 0.0x | PASS |
| FBRAIN2NE | 4 | 2211 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAIN3 | 6 | 2211 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAIN3LS | 6 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAINLS | 2 | 0 | Optimal | Optimal | 7.77e-16 | 14 | 7 | 3.9ms | 3.9ms | 1.0x | PASS |
| FBRAINNE | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 22.4ms | 507us | 0.0x | BOTH_FAIL |
| FCCU | 19 | 8 | Optimal | Optimal | 1.75e-15 | 10 | 9 | 83us | 1.9ms | 23.0x | PASS |
| FEEDLOC | 90 | 259 | Optimal | Optimal | 3.96e-08 | 72 | 23 | 165.1ms | 14.2ms | 0.1x | PASS |
| FLETCHER | 4 | 4 | Optimal | Optimal | 1.23e-08 | 31 | 28 | 113us | 4.9ms | 43.6x | PASS |
| FLT | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 18us | 1000us | 55.2x | PASS |
| GAUSS1 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 3.74s | 338us | 0.0x | BOTH_FAIL |
| GAUSS1LS | 8 | 0 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.41s | 1.5ms | 0.0x | PASS |
| GAUSS2 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 3.19s | 373us | 0.0x | BOTH_FAIL |
| GAUSS2LS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 10 | 5 | 1.31s | 1.3ms | 0.0x | PASS |
| GAUSS3 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 3.20s | 333us | 0.0x | BOTH_FAIL |
| GAUSS3LS | 8 | 0 | Optimal | Optimal | 3.65e-16 | 7 | 11 | 1.02s | 2.7ms | 0.0x | PASS |
| GAUSSIAN | 3 | 0 | Optimal | Optimal | 6.24e-18 | 3 | 2 | 4us | 532us | 118.2x | PASS |
| GBRAIN | 2 | 2200 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 31.1ms | 532us | 0.0x | BOTH_FAIL |
| GBRAINLS | 2 | 0 | Optimal | Optimal | 6.23e-16 | 12 | 6 | 4.7ms | 3.6ms | 0.8x | PASS |
| GENHS28 | 10 | 8 | Optimal | Optimal | 1.22e-15 | 1 | 1 | 18us | 581us | 32.9x | PASS |
| GIGOMEZ1 | 3 | 3 | Optimal | Optimal | 2.87e-09 | 10 | 13 | 27us | 2.5ms | 92.5x | PASS |
| GIGOMEZ2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 28 | 7 | 128us | 1.7ms | 13.4x | MISMATCH |
| GIGOMEZ3 | 3 | 3 | Optimal | Optimal | 4.08e-09 | 10 | 8 | 26us | 1.7ms | 66.0x | PASS |
| GOFFIN | 51 | 50 | Optimal | Optimal | 9.97e-09 | 4 | 7 | 1.3ms | 4.0ms | 3.0x | PASS |
| GOTTFR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 5 | 25us | 1.0ms | 41.8x | PASS |
| GOULDQP1 | 32 | 17 | Acceptable | Optimal | 4.28e-07 | 34 | 15 | 1.1ms | 3.4ms | 3.0x | PASS |
| GROUPING | 100 | 125 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 2.0ms | 312us | 0.2x | ipopt_FAIL |
| GROWTH | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 217us | 295us | 1.4x | BOTH_FAIL |
| GROWTHLS | 3 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 71 | 7us | 11.8ms | 1689.4x | MISMATCH |
| GULF | 3 | 0 | Optimal | Optimal | 1.82e-20 | 50 | 28 | 623us | 5.1ms | 8.2x | PASS |
| GULFNE | 3 | 99 | Optimal | IpoptStatus( | N/A | 22 | 0 | 1.1ms | 300us | 0.3x | ipopt_FAIL |
| HAHN1 | 7 | 236 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| HAHN1LS | 7 | 0 | Acceptable | Optimal | 7.55e-02 | 2999 | 78 | 235.1ms | 17.0ms | 0.1x | MISMATCH |
| HAIFAM | 99 | 150 | Acceptable | Optimal | 8.33e-07 | 27 | 40 | 5.9ms | 16.7ms | 2.8x | PASS |
| HAIFAS | 13 | 9 | Optimal | Optimal | 9.97e-09 | 52 | 16 | 583us | 3.4ms | 5.8x | PASS |
| HAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 53 | 62 | 33us | 10.5ms | 315.8x | PASS |
| HALDMADS | 6 | 42 | Acceptable | Optimal | 9.23e-01 | 86 | 8 | 2.2ms | 3.0ms | 1.3x | MISMATCH |
| HART6 | 6 | 0 | Optimal | Optimal | 0.00e+00 | 18 | 7 | 14us | 1.9ms | 134.8x | PASS |
| HATFLDA | 4 | 0 | Optimal | Optimal | 1.78e-19 | 49 | 13 | 26us | 2.3ms | 87.9x | PASS |
| HATFLDANE | 4 | 4 | Acceptable | Optimal | 0.00e+00 | 9 | 6 | 25us | 1.4ms | 56.3x | PASS |
| HATFLDB | 4 | 0 | Optimal | Optimal | 3.90e-09 | 9 | 8 | 23.3ms | 1.6ms | 0.1x | PASS |
| HATFLDBNE | 4 | 4 | MaxIteration | Infeasible | N/A | 2999 | 13 | 67.9ms | 2.8ms | 0.0x | BOTH_FAIL |
| HATFLDC | 25 | 0 | Acceptable | Optimal | 3.35e-14 | 25 | 5 | 30us | 1.3ms | 41.6x | PASS |
| HATFLDCNE | 25 | 25 | Acceptable | Optimal | 0.00e+00 | 9 | 4 | 121us | 1.2ms | 10.3x | PASS |
| HATFLDD | 3 | 0 | Optimal | Optimal | 1.89e-07 | 27 | 21 | 24us | 2.9ms | 123.0x | PASS |
| HATFLDDNE | 3 | 10 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 140us | 296us | 2.1x | BOTH_FAIL |
| HATFLDE | 3 | 0 | Optimal | Optimal | 2.22e-06 | 32 | 20 | 42us | 2.8ms | 66.0x | PASS |
| HATFLDENE | 3 | 21 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 226us | 335us | 1.5x | BOTH_FAIL |
| HATFLDF | 3 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 135 | 44us | 24.7ms | 559.8x | PASS |
| HATFLDFL | 3 | 0 | Optimal | Optimal | 1.03e-08 | 494 | 1281 | 10.5ms | 187.2ms | 17.8x | PASS |
| HATFLDFLNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 15 | 15 | 717us | 2.8ms | 3.9x | PASS |
| HATFLDFLS | 3 | 0 | Optimal | Optimal | 3.78e-18 | 73 | 36 | 40us | 5.6ms | 140.2x | PASS |
| HATFLDG | 25 | 25 | Optimal | Optimal | 0.00e+00 | 13 | 7 | 227us | 1.7ms | 7.4x | PASS |
| HATFLDGLS | 25 | 0 | Acceptable | Optimal | 2.19e-13 | 44 | 14 | 63us | 2.5ms | 39.5x | PASS |
| HATFLDH | 4 | 7 | Acceptable | Optimal | 1.24e-04 | 96 | 17 | 11.9ms | 3.5ms | 0.3x | MISMATCH |
| HEART6 | 6 | 6 | NumericalErr | Optimal | N/A | 3 | 22 | 1.6ms | 4.9ms | 3.1x | ripopt_FAIL |
| HEART6LS | 6 | 0 | Optimal | Optimal | 9.36e-23 | 907 | 875 | 707us | 137.1ms | 193.9x | PASS |
| HEART8 | 8 | 8 | Optimal | Optimal | 0.00e+00 | 66 | 12 | 321us | 2.2ms | 7.0x | PASS |
| HEART8LS | 8 | 0 | Acceptable | Optimal | 6.00e-18 | 551 | 106 | 443us | 17.2ms | 38.8x | PASS |
| HELIX | 3 | 0 | Optimal | Optimal | 9.65e-19 | 25 | 13 | 15us | 2.1ms | 140.2x | PASS |
| HELIXNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 12 | 7 | 34us | 1.3ms | 38.9x | PASS |
| HET-Z | 2 | 1002 | Optimal | Optimal | 1.85e-03 | 12 | 11 | 6.0ms | 19.1ms | 3.2x | MISMATCH |
| HIELOW | 3 | 0 | Optimal | Optimal | 5.46e-15 | 7 | 8 | 47.8ms | 11.0ms | 0.2x | PASS |
| HIMMELBA | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 6us | 450us | 70.6x | PASS |
| HIMMELBB | 2 | 0 | Optimal | Optimal | 1.39e-17 | 12 | 18 | 9us | 2.9ms | 327.2x | PASS |
| HIMMELBC | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 15us | 1.1ms | 71.7x | PASS |
| HIMMELBCLS | 2 | 0 | Optimal | Optimal | 1.82e-23 | 8 | 6 | 5us | 1.1ms | 221.5x | PASS |
| HIMMELBD | 2 | 2 | RestorationF | Infeasible | N/A | 14 | 22 | 809us | 4.5ms | 5.6x | BOTH_FAIL |
| HIMMELBE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 5 | 2 | 14us | 567us | 41.6x | PASS |
| HIMMELBF | 4 | 0 | Acceptable | Optimal | 4.10e-15 | 56 | 75 | 120us | 11.1ms | 92.0x | PASS |
| HIMMELBFNE | 4 | 7 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 466us | 298us | 0.6x | BOTH_FAIL |
| HIMMELBG | 2 | 0 | Optimal | Optimal | 3.53e-18 | 10 | 6 | 5us | 1.2ms | 274.1x | PASS |
| HIMMELBH | 2 | 0 | Optimal | Optimal | 4.44e-16 | 7 | 4 | 4us | 971us | 228.5x | PASS |
| HIMMELBI | 100 | 12 | Optimal | Optimal | 5.80e-10 | 24 | 13 | 899us | 3.4ms | 3.7x | PASS |
| HIMMELBJ | 45 | 14 | Acceptable | ErrorInStepC | N/A | 34 | 580 | 2.4ms | 138.0ms | 58.3x | ipopt_FAIL |
| HIMMELBK | 24 | 14 | Optimal | Optimal | 4.89e-08 | 19 | 18 | 558us | 3.8ms | 6.8x | PASS |
| HIMMELP1 | 2 | 0 | Acceptable | Optimal | 1.83e-15 | 21 | 10 | 55us | 2.2ms | 39.4x | PASS |
| HIMMELP2 | 2 | 1 | Optimal | Optimal | 8.68e-01 | 12 | 17 | 28us | 3.6ms | 125.7x | MISMATCH |
| HIMMELP3 | 2 | 2 | Optimal | Optimal | 8.66e-01 | 0 | 11 | 4us | 2.2ms | 502.2x | MISMATCH |
| HIMMELP4 | 2 | 3 | Optimal | Optimal | 1.23e-08 | 9 | 23 | 25us | 4.6ms | 184.2x | PASS |
| HIMMELP5 | 2 | 3 | Optimal | Optimal | 7.50e-01 | 9 | 46 | 31us | 8.3ms | 263.0x | MISMATCH |
| HIMMELP6 | 2 | 5 | Optimal | Optimal | 7.50e-01 | 11 | 31 | 37us | 6.0ms | 160.5x | MISMATCH |
| HONG | 4 | 1 | Optimal | Optimal | 4.72e-16 | 9 | 7 | 27us | 1.7ms | 62.5x | PASS |
| HS1 | 2 | 0 | Optimal | Optimal | 9.07e-21 | 25 | 28 | 14us | 5.2ms | 376.0x | PASS |
| HS10 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 12 | 33us | 2.5ms | 77.6x | PASS |
| HS100 | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 70us | 2.0ms | 29.1x | PASS |
| HS100LNP | 7 | 2 | Optimal | Optimal | 1.67e-16 | 6 | 20 | 26us | 2.7ms | 105.3x | PASS |
| HS100MOD | 7 | 4 | Optimal | Optimal | 1.73e-11 | 9 | 14 | 59us | 2.7ms | 46.5x | PASS |
| HS101 | 7 | 5 | Optimal | Optimal | 2.36e-08 | 30 | 39 | 285us | 9.3ms | 32.4x | PASS |
| HS102 | 7 | 5 | Optimal | Optimal | 4.27e-08 | 24 | 52 | 247us | 9.9ms | 40.1x | PASS |
| HS103 | 7 | 5 | Optimal | Optimal | 4.13e-08 | 26 | 21 | 349us | 4.6ms | 13.3x | PASS |
| HS104 | 8 | 5 | Optimal | Optimal | 1.50e-08 | 24 | 8 | 206us | 1.8ms | 8.9x | PASS |
| HS105 | 8 | 1 | Optimal | Optimal | 9.48e-12 | 30 | 23 | 3.5ms | 7.3ms | 2.1x | PASS |
| HS106 | 8 | 6 | Optimal | Optimal | 1.74e-08 | 18 | 18 | 66us | 3.2ms | 48.9x | PASS |
| HS107 | 9 | 6 | Acceptable | Optimal | 1.29e-07 | 52 | 7 | 298us | 1.8ms | 6.2x | PASS |
| HS108 | 9 | 13 | Optimal | Optimal | 8.77e-09 | 2596 | 11 | 15.5ms | 2.6ms | 0.2x | PASS |
| HS109 | 9 | 10 | Acceptable | Optimal | 2.16e-01 | 2999 | 14 | 22.5ms | 2.8ms | 0.1x | MISMATCH |
| HS11 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 6 | 14us | 1.3ms | 89.2x | PASS |
| HS111 | 10 | 3 | Optimal | Optimal | 1.05e-11 | 11 | 15 | 105us | 2.8ms | 27.0x | PASS |
| HS111LNP | 10 | 3 | Optimal | Optimal | 1.03e-10 | 11 | 15 | 90us | 2.5ms | 28.2x | PASS |
| HS112 | 10 | 3 | Optimal | Optimal | 6.49e-14 | 10 | 10 | 81us | 2.0ms | 24.8x | PASS |
| HS113 | 10 | 8 | Optimal | Optimal | 1.68e-09 | 14 | 9 | 109us | 1.9ms | 17.3x | PASS |
| HS114 | 10 | 11 | Optimal | Optimal | 1.04e-07 | 19 | 13 | 118us | 2.4ms | 20.5x | PASS |
| HS116 | 13 | 14 | Acceptable | Optimal | 4.43e-04 | 131 | 19 | 1.6ms | 3.6ms | 2.3x | MISMATCH |
| HS117 | 15 | 5 | Acceptable | Optimal | 3.56e-07 | 38 | 19 | 308us | 3.6ms | 11.8x | PASS |
| HS118 | 15 | 17 | Optimal | Optimal | 1.39e-09 | 236 | 10 | 4.0ms | 2.1ms | 0.5x | PASS |
| HS119 | 16 | 8 | Acceptable | Optimal | 8.77e-05 | 21 | 17 | 329us | 3.2ms | 9.7x | PASS |
| HS12 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 7 | 6 | 19us | 1.3ms | 65.2x | PASS |
| HS13 | 2 | 1 | Acceptable | Optimal | 7.19e-04 | 24 | 47 | 38us | 7.2ms | 189.1x | MISMATCH |
| HS14 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 5 | 16us | 1.2ms | 73.0x | PASS |
| HS15 | 2 | 2 | Optimal | Optimal | 8.08e-08 | 11 | 13 | 31us | 2.3ms | 75.0x | PASS |
| HS16 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 10 | 10 | 23us | 2.1ms | 89.8x | MISMATCH |
| HS17 | 2 | 2 | Optimal | Optimal | 2.01e-04 | 12 | 22 | 26us | 3.9ms | 150.4x | MISMATCH |
| HS18 | 2 | 2 | Optimal | Optimal | 9.03e-10 | 14 | 10 | 30us | 1.9ms | 63.1x | PASS |
| HS19 | 2 | 2 | Optimal | Optimal | 3.09e-09 | 20 | 12 | 39us | 2.3ms | 57.8x | PASS |
| HS1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 25 | 30 | 43us | 6.0ms | 139.1x | PASS |
| HS2 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 10 | 10 | 19.0ms | 2.1ms | 0.1x | PASS |
| HS20 | 2 | 3 | Optimal | Optimal | 6.17e-08 | 18 | 5 | 50us | 1.1ms | 22.3x | PASS |
| HS21 | 2 | 1 | Optimal | Optimal | 1.38e-10 | 9 | 6 | 19us | 1.3ms | 67.3x | PASS |
| HS21MOD | 7 | 1 | Acceptable | Optimal | 1.70e-08 | 11 | 13 | 27us | 2.5ms | 90.1x | PASS |
| HS22 | 2 | 2 | Optimal | Optimal | 1.16e-08 | 5 | 5 | 13us | 1.1ms | 84.2x | PASS |
| HS23 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 103 | 9 | 1.3ms | 1.7ms | 1.3x | MISMATCH |
| HS24 | 2 | 3 | Optimal | Optimal | 1.80e-08 | 7 | 14 | 24us | 2.8ms | 117.7x | PASS |
| HS25 | 3 | 0 | Acceptable | Optimal | 6.08e-05 | 16 | 27 | 624.0ms | 5.3ms | 0.0x | PASS |
| HS25NE | 3 | 99 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 490us | 307us | 0.6x | BOTH_FAIL |
| HS26 | 3 | 1 | Acceptable | Optimal | 2.94e-12 | 18 | 25 | 37us | 3.3ms | 89.0x | PASS |
| HS268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 18us | 2.5ms | 140.0x | MISMATCH |
| HS27 | 3 | 1 | Optimal | Optimal | 6.69e-13 | 17 | 57 | 42us | 7.8ms | 185.3x | PASS |
| HS28 | 3 | 1 | Optimal | Optimal | 9.24e-31 | 1 | 1 | 7us | 494us | 67.8x | PASS |
| HS29 | 3 | 1 | Optimal | Optimal | 3.10e-10 | 18 | 7 | 51us | 1.6ms | 30.7x | PASS |
| HS2NE | 2 | 2 | MaxIteration | Infeasible | N/A | 2999 | 12 | 69.5ms | 2.6ms | 0.0x | BOTH_FAIL |
| HS3 | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 44us | 901us | 20.4x | PASS |
| HS30 | 3 | 1 | Acceptable | Optimal | 1.45e-05 | 9 | 7 | 20us | 1.5ms | 72.8x | PASS |
| HS31 | 3 | 1 | Optimal | Optimal | 8.67e-09 | 11 | 6 | 23us | 1.4ms | 62.1x | PASS |
| HS32 | 3 | 2 | Acceptable | Optimal | 5.45e-05 | 12 | 15 | 32us | 2.7ms | 85.4x | PASS |
| HS33 | 3 | 2 | Acceptable | Optimal | 8.21e-08 | 14 | 9 | 42us | 1.8ms | 43.5x | PASS |
| HS34 | 3 | 2 | Optimal | Optimal | 9.46e-09 | 13 | 7 | 32us | 1.6ms | 49.6x | PASS |
| HS35 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 7 | 23us | 1.4ms | 61.9x | PASS |
| HS35I | 3 | 1 | Optimal | Optimal | 3.44e-09 | 9 | 7 | 22us | 1.6ms | 74.2x | PASS |
| HS35MOD | 3 | 1 | Acceptable | Optimal | 4.67e-07 | 8 | 14 | 19us | 2.4ms | 127.0x | PASS |
| HS36 | 3 | 1 | Optimal | Optimal | 6.33e-09 | 16 | 11 | 37us | 2.4ms | 66.3x | PASS |
| HS37 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 10 | 11 | 35us | 2.4ms | 67.7x | PASS |
| HS38 | 4 | 0 | Optimal | Optimal | 1.23e-21 | 24 | 39 | 18us | 6.6ms | 374.2x | PASS |
| HS39 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 34us | 1.9ms | 55.4x | PASS |
| HS3MOD | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 57us | 942us | 16.6x | PASS |
| HS4 | 2 | 0 | Optimal | Optimal | 1.97e-08 | 6 | 4 | 17.1ms | 983us | 0.1x | PASS |
| HS40 | 4 | 3 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 16us | 699us | 44.7x | PASS |
| HS41 | 4 | 1 | Optimal | Optimal | 1.15e-09 | 11 | 7 | 27us | 1.4ms | 52.2x | PASS |
| HS42 | 4 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 18us | 801us | 44.9x | PASS |
| HS43 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 8 | 8 | 32us | 1.8ms | 57.8x | PASS |
| HS44 | 4 | 6 | Acceptable | Optimal | 2.97e-06 | 28 | 24 | 100us | 4.5ms | 44.8x | PASS |
| HS44NEW | 4 | 6 | Acceptable | Optimal | 1.15e-05 | 19 | 18 | 75us | 3.9ms | 51.8x | PASS |
| HS45 | 5 | 0 | Acceptable | Optimal | 3.14e-04 | 31 | 11 | 15.8ms | 2.3ms | 0.1x | MISMATCH |
| HS46 | 5 | 2 | Acceptable | Optimal | 3.43e-11 | 17 | 19 | 50us | 2.5ms | 50.8x | PASS |
| HS47 | 5 | 3 | Acceptable | Optimal | 3.22e-12 | 17 | 19 | 47us | 2.6ms | 55.4x | PASS |
| HS48 | 5 | 2 | Optimal | Optimal | 4.44e-31 | 1 | 1 | 9us | 508us | 55.9x | PASS |
| HS49 | 5 | 2 | Acceptable | Optimal | 2.61e-10 | 17 | 19 | 38us | 2.6ms | 67.8x | PASS |
| HS5 | 2 | 0 | Optimal | Optimal | 2.32e-16 | 7 | 7 | 4us | 1.5ms | 402.1x | PASS |
| HS50 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 28us | 1.4ms | 50.6x | PASS |
| HS51 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 9us | 506us | 55.2x | PASS |
| HS52 | 5 | 3 | Optimal | Optimal | 1.17e-15 | 1 | 1 | 11us | 499us | 45.2x | PASS |
| HS53 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 6 | 22us | 1.2ms | 56.5x | PASS |
| HS54 | 6 | 1 | Optimal | Optimal | 7.51e-01 | 10 | 15 | 32us | 3.0ms | 93.6x | MISMATCH |
| HS55 | 6 | 6 | Optimal | Optimal | 2.04e-02 | 11 | 18 | 42us | 4.4ms | 104.4x | MISMATCH |
| HS56 | 7 | 4 | Optimal | Optimal | 1.99e-14 | 5 | 10 | 30us | 1.7ms | 58.1x | PASS |
| HS57 | 2 | 1 | Optimal | Optimal | 6.11e-15 | 10 | 10 | 40us | 1.7ms | 42.1x | PASS |
| HS59 | 2 | 3 | Acceptable | Optimal | 9.28e-04 | 10 | 17 | 44us | 3.4ms | 76.9x | MISMATCH |
| HS6 | 2 | 1 | Optimal | Optimal | 4.93e-32 | 7 | 5 | 20us | 1.1ms | 52.6x | PASS |
| HS60 | 3 | 1 | Optimal | Optimal | 1.19e-13 | 8 | 6 | 19us | 1.6ms | 84.7x | PASS |
| HS61 | 3 | 2 | Optimal | Optimal | 7.91e-16 | 9 | 10 | 23us | 1.5ms | 67.5x | PASS |
| HS62 | 3 | 1 | Optimal | Optimal | 1.38e-16 | 9 | 6 | 26us | 1.6ms | 60.5x | PASS |
| HS63 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 5 | 22us | 1.2ms | 52.5x | PASS |
| HS64 | 3 | 1 | Optimal | Optimal | 3.62e-09 | 18 | 16 | 36us | 3.0ms | 81.1x | PASS |
| HS65 | 3 | 1 | Optimal | Optimal | 7.39e-09 | 15 | 16 | 39us | 3.3ms | 84.6x | PASS |
| HS66 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 10 | 10 | 27us | 1.8ms | 65.6x | PASS |
| HS67 | 3 | 14 | Optimal | Optimal | 2.99e-09 | 12 | 9 | 102us | 1.9ms | 19.1x | PASS |
| HS68 | 4 | 2 | Optimal | Optimal | 5.21e-11 | 18 | 16 | 54us | 3.1ms | 57.7x | PASS |
| HS69 | 4 | 2 | Optimal | Optimal | 2.85e-15 | 11 | 10 | 35us | 2.0ms | 58.1x | PASS |
| HS7 | 2 | 1 | Optimal | Optimal | 5.30e-12 | 10 | 27 | 20us | 4.4ms | 223.7x | PASS |
| HS70 | 4 | 1 | Optimal | Optimal | 1.80e-01 | 10 | 46 | 119us | 8.5ms | 71.8x | MISMATCH |
| HS71 | 4 | 2 | Optimal | Optimal | 1.01e-09 | 10 | 8 | 31us | 1.9ms | 60.0x | PASS |
| HS72 | 4 | 2 | Optimal | Optimal | 6.74e-07 | 13 | 16 | 30us | 2.9ms | 96.2x | PASS |
| HS73 | 4 | 3 | Optimal | Optimal | 1.30e-09 | 11 | 8 | 33us | 2.0ms | 61.1x | PASS |
| HS74 | 4 | 5 | Optimal | Optimal | 0.00e+00 | 15 | 8 | 59us | 1.9ms | 32.7x | PASS |
| HS75 | 4 | 5 | Optimal | Optimal | 5.37e-09 | 15 | 8 | 52us | 1.7ms | 33.2x | PASS |
| HS76 | 4 | 3 | Optimal | Optimal | 5.89e-09 | 9 | 7 | 24us | 1.6ms | 64.0x | PASS |
| HS76I | 4 | 3 | Optimal | Optimal | 3.59e-09 | 10 | 6 | 26us | 1.3ms | 49.6x | PASS |
| HS77 | 5 | 2 | Optimal | Optimal | 1.99e-11 | 9 | 11 | 28us | 1.6ms | 59.2x | PASS |
| HS78 | 5 | 3 | Optimal | Optimal | 1.52e-16 | 4 | 4 | 18us | 820us | 45.4x | PASS |
| HS79 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 19us | 806us | 42.5x | PASS |
| HS8 | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 17us | 943us | 55.3x | PASS |
| HS80 | 5 | 3 | Optimal | Optimal | 6.75e-13 | 9 | 5 | 29us | 1.2ms | 42.2x | PASS |
| HS81 | 5 | 3 | Optimal | Optimal | 3.46e-14 | 31 | 68 | 108us | 12.0ms | 111.4x | PASS |
| HS83 | 5 | 3 | Acceptable | Optimal | 2.47e-06 | 8991 | 9 | 232.5ms | 1.8ms | 0.0x | PASS |
| HS84 | 5 | 3 | Acceptable | Optimal | 1.04e-05 | 2999 | 9 | 107.6ms | 2.0ms | 0.0x | PASS |
| HS85 | 5 | 21 | Acceptable | Optimal | 6.89e-03 | 2999 | 13 | 1.62s | 4.1ms | 0.0x | MISMATCH |
| HS86 | 5 | 10 | Optimal | Optimal | 2.21e-09 | 11 | 10 | 51us | 2.3ms | 43.9x | PASS |
| HS87 | 6 | 4 | Acceptable | MaxIteration | N/A | 5994 | 3000 | 205.3ms | 484.0ms | 2.4x | ipopt_FAIL |
| HS88 | 2 | 1 | Optimal | Optimal | 6.34e-06 | 24 | 18 | 966us | 3.7ms | 3.9x | PASS |
| HS89 | 3 | 1 | Optimal | Optimal | 3.42e-06 | 45 | 15 | 4.6ms | 3.5ms | 0.8x | PASS |
| HS9 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 13us | 778us | 59.5x | PASS |
| HS90 | 4 | 1 | Optimal | Optimal | 5.20e-06 | 18 | 16 | 1.3ms | 4.4ms | 3.5x | PASS |
| HS91 | 5 | 1 | Optimal | Optimal | 6.44e-06 | 18 | 16 | 1.8ms | 4.3ms | 2.4x | PASS |
| HS92 | 6 | 1 | Optimal | Optimal | 6.80e-06 | 16 | 35 | 1.8ms | 9.9ms | 5.6x | PASS |
| HS93 | 6 | 2 | Optimal | Optimal | 9.89e-09 | 14 | 7 | 65us | 1.7ms | 26.2x | PASS |
| HS95 | 6 | 4 | Acceptable | Optimal | 6.92e-04 | 46 | 9 | 147us | 2.0ms | 13.6x | MISMATCH |
| HS96 | 6 | 4 | Acceptable | Optimal | 6.88e-04 | 46 | 8 | 144us | 2.0ms | 14.0x | MISMATCH |
| HS97 | 6 | 4 | Acceptable | Optimal | 1.05e-06 | 48 | 24 | 225us | 4.4ms | 19.7x | PASS |
| HS98 | 6 | 4 | Acceptable | Optimal | 1.25e-04 | 156 | 13 | 506us | 2.6ms | 5.2x | MISMATCH |
| HS99 | 7 | 2 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 34us | 1.2ms | 34.5x | PASS |
| HS99EXP | 31 | 21 | Acceptable | Optimal | 1.78e-14 | 21 | 17 | 844us | 3.3ms | 4.0x | PASS |
| HUBFIT | 2 | 1 | Optimal | Optimal | 2.24e-09 | 7 | 7 | 20us | 1.6ms | 82.3x | PASS |
| HUMPS | 2 | 0 | Acceptable | Optimal | 2.78e-14 | 132 | 1533 | 92us | 213.6ms | 2317.6x | PASS |
| HYDC20LS | 99 | 0 | Acceptable | Optimal | 2.98e-01 | 2999 | 639 | 1.05s | 170.2ms | 0.2x | MISMATCH |
| HYDCAR20 | 99 | 99 | Optimal | Optimal | 0.00e+00 | 11 | 9 | 177.4ms | 3.0ms | 0.0x | PASS |
| HYDCAR6 | 29 | 29 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 9.8ms | 1.3ms | 0.1x | PASS |
| HYDCAR6LS | 29 | 0 | Optimal | Optimal | 2.86e-18 | 1322 | 149 | 35.4ms | 30.1ms | 0.8x | PASS |
| HYPCIR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 20us | 1.2ms | 59.0x | PASS |
| JENSMP | 2 | 0 | Optimal | Optimal | 9.38e-01 | 1 | 9 | 4us | 1.5ms | 377.6x | MISMATCH |
| JENSMPNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 46us | 293us | 6.4x | BOTH_FAIL |
| JUDGE | 2 | 0 | Optimal | Optimal | 2.15e-01 | 17 | 9 | 19us | 1.4ms | 73.0x | MISMATCH |
| JUDGEB | 2 | 0 | Optimal | Optimal | 2.15e-01 | 17 | 9 | 41us | 1.9ms | 45.7x | MISMATCH |
| JUDGENE | 2 | 20 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 79us | 293us | 3.7x | BOTH_FAIL |
| KIRBY2 | 5 | 151 | LocalInfeasi | IpoptStatus( | N/A | 118 | 0 | 4.9ms | 288us | 0.1x | BOTH_FAIL |
| KIRBY2LS | 5 | 0 | Acceptable | Optimal | 1.21e-14 | 26 | 11 | 1.9ms | 2.1ms | 1.1x | PASS |
| KIWCRESC | 3 | 2 | Optimal | Optimal | 1.07e-08 | 13 | 8 | 37us | 1.8ms | 48.0x | PASS |
| KOEBHELB | 3 | 0 | Optimal | Optimal | 1.83e-16 | 164 | 71 | 1.1ms | 14.7ms | 13.1x | PASS |
| KOEBHELBNE | 3 | 156 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 5.7ms | 281us | 0.0x | BOTH_FAIL |
| KOWOSB | 4 | 0 | Optimal | Optimal | 4.11e-17 | 21 | 8 | 19us | 1.8ms | 92.1x | PASS |
| KOWOSBNE | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 89us | 278us | 3.1x | BOTH_FAIL |
| KSIP | 20 | 1001 | Optimal | Optimal | 9.79e-01 | 1 | 22 | 12.2ms | 68.7ms | 5.6x | MISMATCH |
| LAKES | 90 | 78 | Optimal | Optimal | 1.53e-13 | 13 | 11 | 1.0ms | 3.2ms | 3.1x | PASS |
| LANCZOS1 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 41 | 0 | 422us | 321us | 0.8x | ipopt_FAIL |
| LANCZOS1LS | 6 | 0 | Acceptable | Optimal | 4.35e-06 | 55 | 115 | 97us | 19.6ms | 201.3x | PASS |
| LANCZOS2 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 70 | 0 | 665us | 357us | 0.5x | ipopt_FAIL |
| LANCZOS2LS | 6 | 0 | Acceptable | Optimal | 4.37e-06 | 54 | 101 | 92us | 16.3ms | 176.5x | PASS |
| LANCZOS3 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 30 | 0 | 292us | 286us | 1.0x | ipopt_FAIL |
| LANCZOS3LS | 6 | 0 | Optimal | Optimal | 4.33e-06 | 67 | 174 | 125us | 30.6ms | 244.2x | PASS |
| LAUNCH | 25 | 28 | Acceptable | Optimal | 3.45e-04 | 2999 | 12 | 189.9ms | 2.8ms | 0.0x | MISMATCH |
| LEVYMONE10 | 10 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 115us | 346us | 3.0x | BOTH_FAIL |
| LEVYMONE5 | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 40us | 310us | 7.8x | BOTH_FAIL |
| LEVYMONE6 | 3 | 6 | Optimal | IpoptStatus( | N/A | 29 | 0 | 59us | 281us | 4.8x | ipopt_FAIL |
| LEVYMONE7 | 4 | 8 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 49us | 284us | 5.8x | BOTH_FAIL |
| LEVYMONE8 | 5 | 10 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 53us | 282us | 5.4x | BOTH_FAIL |
| LEVYMONE9 | 8 | 16 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 78us | 336us | 4.3x | BOTH_FAIL |
| LEVYMONT10 | 10 | 0 | Optimal | Optimal | 1.74e-16 | 14 | 4 | 19us | 1.2ms | 61.0x | PASS |
| LEVYMONT5 | 2 | 0 | Acceptable | Optimal | 1.00e+00 | 10 | 10 | 19us | 2.1ms | 112.2x | MISMATCH |
| LEVYMONT6 | 3 | 0 | Optimal | Optimal | 1.00e+00 | 27 | 8 | 20us | 1.9ms | 96.2x | MISMATCH |
| LEVYMONT7 | 4 | 0 | Optimal | Optimal | 1.88e-01 | 14 | 7 | 10us | 1.7ms | 169.1x | MISMATCH |
| LEVYMONT8 | 5 | 0 | Optimal | Optimal | 9.71e-01 | 12 | 4 | 14us | 933us | 67.6x | MISMATCH |
| LEVYMONT9 | 8 | 0 | Optimal | Optimal | 6.77e-01 | 12 | 4 | 15us | 970us | 63.6x | MISMATCH |
| LEWISPOL | 6 | 9 | Acceptable | IpoptStatus( | N/A | 10 | 0 | 929us | 281us | 0.3x | ipopt_FAIL |
| LHAIFAM | 99 | 150 | Optimal | InvalidNumbe | N/A | 0 | 0 | 941.5ms | 350us | 0.0x | ipopt_FAIL |
| LIN | 4 | 2 | Optimal | Optimal | 2.03e-03 | 9 | 7 | 32us | 1.5ms | 47.9x | MISMATCH |
| LINSPANH | 97 | 33 | Acceptable | Optimal | 5.76e-07 | 2999 | 24 | 116.3ms | 5.2ms | 0.0x | PASS |
| LOADBAL | 31 | 31 | Optimal | Optimal | 6.55e-09 | 15 | 13 | 793us | 2.8ms | 3.6x | PASS |
| LOGHAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 408 | 2747 | 280us | 395.9ms | 1411.9x | PASS |
| LOGROS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 50 | 49 | 112us | 9.9ms | 88.1x | PASS |
| LOOTSMA | 3 | 2 | Optimal | Optimal | 8.76e-08 | 26 | 13 | 66us | 3.0ms | 45.0x | PASS |
| LOTSCHD | 12 | 7 | Optimal | Optimal | 4.44e-10 | 14 | 9 | 78us | 1.8ms | 23.7x | PASS |
| LRCOVTYPE | 54 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| LRIJCNN1 | 22 | 0 | Acceptable | Optimal | 5.92e-03 | 31 | 11 | 83.2ms | 190.1ms | 2.3x | MISMATCH |
| LSC1 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 55 | 0 | 102us | 285us | 2.8x | BOTH_FAIL |
| LSC1LS | 3 | 0 | Optimal | Optimal | 1.27e-15 | 49 | 16 | 33us | 3.2ms | 96.7x | PASS |
| LSC2 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 93 | 0 | 325us | 317us | 1.0x | BOTH_FAIL |
| LSC2LS | 3 | 0 | Acceptable | Optimal | 2.32e-04 | 32 | 38 | 29.5ms | 4.8ms | 0.2x | MISMATCH |
| LSNNODOC | 5 | 4 | Acceptable | Optimal | 3.42e-06 | 13 | 10 | 41us | 2.1ms | 50.8x | PASS |
| LSQFIT | 2 | 1 | Optimal | Optimal | 4.24e-09 | 7 | 7 | 19us | 1.4ms | 74.2x | PASS |
| MADSEN | 3 | 6 | Optimal | Optimal | 8.02e-09 | 29 | 18 | 82us | 3.4ms | 41.8x | PASS |
| MAKELA1 | 3 | 2 | Optimal | Optimal | 2.00e+00 | 6 | 12 | 23us | 2.4ms | 106.2x | MISMATCH |
| MAKELA2 | 3 | 3 | Optimal | Optimal | 1.27e-01 | 3 | 6 | 14us | 1.3ms | 93.3x | MISMATCH |
| MAKELA3 | 21 | 20 | Acceptable | Optimal | 7.98e-09 | 35 | 11 | 678us | 2.5ms | 3.8x | PASS |
| MAKELA4 | 21 | 40 | Optimal | Optimal | 6.88e-01 | 1 | 5 | 64us | 1.4ms | 21.0x | MISMATCH |
| MARATOS | 2 | 1 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 14us | 842us | 60.9x | PASS |
| MARATOSB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1107 | 672 | 524us | 101.6ms | 193.8x | PASS |
| MATRIX2 | 6 | 2 | Acceptable | Optimal | 7.84e-10 | 12 | 42 | 41us | 7.4ms | 178.9x | PASS |
| MAXLIKA | 8 | 0 | Acceptable | Optimal | 1.13e-02 | 20 | 23 | 2.06s | 6.7ms | 0.0x | MISMATCH |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 229us | 1.7ms | 7.2x | PASS |
| MDHOLE | 2 | 0 | Optimal | Optimal | 9.98e-09 | 35 | 42 | 172us | 8.8ms | 51.1x | PASS |
| MESH | 41 | 48 | Acceptable | IpoptStatus( | N/A | 71 | 79 | 17.4ms | 20.2ms | 1.2x | ipopt_FAIL |
| METHANB8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 9 | 3 | 286us | 885us | 3.1x | PASS |
| METHANB8LS | 31 | 0 | Optimal | Optimal | 5.35e-26 | 9 | 8 | 8.2ms | 1.5ms | 0.2x | PASS |
| METHANL8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 9.7ms | 1.0ms | 0.1x | PASS |
| METHANL8LS | 31 | 0 | Optimal | Optimal | 5.70e-17 | 977 | 40 | 30.0ms | 8.8ms | 0.3x | PASS |
| MEXHAT | 2 | 0 | Optimal | Optimal | 6.77e-10 | 45 | 26 | 28us | 4.0ms | 142.4x | PASS |
| MEYER3 | 3 | 0 | Acceptable | Optimal | 2.23e-12 | 2999 | 194 | 65.5ms | 29.6ms | 0.5x | PASS |
| MEYER3NE | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 150.0ms | 284us | 0.0x | BOTH_FAIL |
| MGH09 | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 424us | 277us | 0.7x | BOTH_FAIL |
| MGH09LS | 4 | 0 | Acceptable | Optimal | 7.12e-04 | 47 | 72 | 35us | 11.3ms | 321.3x | MISMATCH |
| MGH10 | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 112.4ms | 322us | 0.0x | BOTH_FAIL |
| MGH10LS | 3 | 0 | Acceptable | Optimal | 4.94e-12 | 2999 | 1828 | 49.6ms | 286.9ms | 5.8x | PASS |
| MGH10S | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 13.9ms | 275us | 0.0x | BOTH_FAIL |
| MGH10SLS | 3 | 0 | Acceptable | Optimal | 1.00e+00 | 1304 | 354 | 1.4ms | 56.8ms | 40.2x | MISMATCH |
| MGH17 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 875us | 291us | 0.3x | BOTH_FAIL |
| MGH17LS | 5 | 0 | Acceptable | Optimal | 2.43e-05 | 214 | 47 | 503us | 9.5ms | 18.9x | PASS |
| MGH17S | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 887us | 323us | 0.4x | BOTH_FAIL |
| MGH17SLS | 5 | 0 | Optimal | Optimal | 4.18e-07 | 71 | 41 | 204us | 7.9ms | 39.0x | PASS |
| MIFFLIN1 | 3 | 2 | Optimal | Optimal | 9.24e-09 | 6 | 5 | 26us | 1.3ms | 49.0x | PASS |
| MIFFLIN2 | 3 | 2 | Optimal | Optimal | 9.97e-09 | 15 | 11 | 38us | 2.5ms | 65.3x | PASS |
| MINMAXBD | 5 | 20 | Optimal | Optimal | 5.55e-11 | 469 | 25 | 4.4ms | 6.4ms | 1.5x | PASS |
| MINMAXRB | 3 | 4 | Optimal | Optimal | 9.85e-09 | 3 | 8 | 15us | 1.6ms | 109.1x | PASS |
| MINSURF | 64 | 0 | Optimal | Optimal | 0.00e+00 | 17 | 4 | 53us | 1.4ms | 26.1x | PASS |
| MISRA1A | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 473us | 280us | 0.6x | BOTH_FAIL |
| MISRA1ALS | 2 | 0 | Acceptable | Optimal | 1.91e-14 | 35 | 40 | 160us | 6.2ms | 38.8x | PASS |
| MISRA1B | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 284us | 315us | 1.1x | BOTH_FAIL |
| MISRA1BLS | 2 | 0 | Optimal | Optimal | 4.76e-14 | 25 | 34 | 113us | 5.2ms | 45.7x | PASS |
| MISRA1C | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 713us | 324us | 0.5x | BOTH_FAIL |
| MISRA1CLS | 2 | 0 | Acceptable | Optimal | 3.55e-14 | 59 | 14 | 162us | 2.5ms | 15.3x | PASS |
| MISRA1D | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 24 | 0 | 106.1ms | 274us | 0.0x | BOTH_FAIL |
| MISRA1DLS | 2 | 0 | Acceptable | Optimal | 4.95e-15 | 56 | 30 | 172us | 4.5ms | 26.5x | PASS |
| MISTAKE | 9 | 13 | Acceptable | Optimal | 5.00e-01 | 81 | 16 | 955us | 3.3ms | 3.4x | MISMATCH |
| MRIBASIS | 36 | 55 | Acceptable | Optimal | 1.00e-08 | 35 | 15 | 5.1ms | 4.6ms | 0.9x | PASS |
| MSS1 | 90 | 73 | Acceptable | Optimal | 1.25e-01 | 74 | 95 | 40.1ms | 51.7ms | 1.3x | MISMATCH |
| MUONSINE | 1 | 512 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 4.70s | 334us | 0.0x | BOTH_FAIL |
| MUONSINELS | 1 | 0 | Acceptable | Optimal | 1.42e-01 | 13 | 8 | 6.4ms | 2.0ms | 0.3x | MISMATCH |
| MWRIGHT | 5 | 3 | Optimal | Optimal | 9.48e-01 | 13 | 10 | 43us | 1.9ms | 43.6x | MISMATCH |
| NASH | 72 | 24 | RestorationF | Infeasible | N/A | 50 | 45 | 26.1ms | 12.0ms | 0.5x | BOTH_FAIL |
| NELSON | 3 | 128 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 152.7ms | 287us | 0.0x | BOTH_FAIL |
| NET1 | 48 | 57 | Acceptable | Optimal | 9.62e-08 | 45 | 26 | 78.1ms | 6.5ms | 0.1x | PASS |
| NYSTROM5 | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 283us | 311us | 1.1x | ipopt_FAIL |
| NYSTROM5C | 18 | 20 | Acceptable | IpoptStatus( | N/A | 18 | 0 | 312us | 277us | 0.9x | ipopt_FAIL |
| ODFITS | 10 | 6 | Optimal | Optimal | 1.91e-16 | 11 | 8 | 50us | 1.8ms | 36.5x | PASS |
| OET1 | 3 | 1002 | Optimal | Optimal | 4.53e-02 | 53 | 33 | 19.2ms | 48.9ms | 2.5x | MISMATCH |
| OET2 | 3 | 1002 | Acceptable | Optimal | 9.56e-01 | 1107 | 181 | 527.5ms | 292.2ms | 0.6x | MISMATCH |
| OET3 | 4 | 1002 | Optimal | Optimal | 5.47e-04 | 94 | 13 | 70.8ms | 23.3ms | 0.3x | MISMATCH |
| OET4 | 4 | 1002 | Optimal | Optimal | 7.66e-09 | 288 | 165 | 97.8ms | 262.8ms | 2.7x | PASS |
| OET5 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET6 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET7 | 7 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OPTCNTRL | 32 | 20 | Acceptable | Optimal | 2.02e-08 | 40 | 9 | 1.2ms | 2.0ms | 1.7x | PASS |
| OPTPRLOC | 30 | 30 | Optimal | Optimal | 1.22e-08 | 57 | 13 | 2.6ms | 3.2ms | 1.2x | PASS |
| ORTHREGB | 27 | 6 | Optimal | Optimal | 4.20e-19 | 2 | 2 | 53us | 864us | 16.3x | PASS |
| OSBORNE1 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 102 | 0 | 31.3ms | 363us | 0.0x | BOTH_FAIL |
| OSBORNE2 | 11 | 65 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 1.9ms | 295us | 0.2x | BOTH_FAIL |
| OSBORNEA | 5 | 0 | Acceptable | Optimal | 2.05e-18 | 99 | 64 | 323us | 10.6ms | 33.0x | PASS |
| OSBORNEB | 11 | 0 | Acceptable | Optimal | 2.34e-12 | 95 | 19 | 517us | 3.3ms | 6.5x | PASS |
| OSLBQP | 8 | 0 | Acceptable | Optimal | 7.24e-07 | 13 | 15 | 20.5ms | 2.7ms | 0.1x | PASS |
| PALMER1 | 4 | 0 | Acceptable | Optimal | 1.55e-16 | 34 | 13 | 96us | 2.7ms | 28.0x | PASS |
| PALMER1A | 6 | 0 | Acceptable | Optimal | 1.33e-14 | 173 | 48 | 677us | 9.9ms | 14.6x | PASS |
| PALMER1ANE | 6 | 35 | LocalInfeasi | IpoptStatus( | N/A | 187 | 0 | 1.4ms | 297us | 0.2x | BOTH_FAIL |
| PALMER1B | 4 | 0 | Acceptable | Optimal | 4.87e-14 | 51 | 17 | 314us | 3.2ms | 10.1x | PASS |
| PALMER1BNE | 4 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 1.2ms | 408us | 0.3x | BOTH_FAIL |
| PALMER1C | 8 | 0 | Optimal | Optimal | 1.79e-13 | 4 | 1 | 60.7ms | 440us | 0.0x | PASS |
| PALMER1D | 7 | 0 | Acceptable | Optimal | 1.65e-13 | 196 | 1 | 481us | 473us | 1.0x | PASS |
| PALMER1E | 8 | 0 | Acceptable | Optimal | 4.05e-13 | 675 | 55 | 1.7ms | 11.4ms | 6.6x | PASS |
| PALMER1ENE | 8 | 35 | LocalInfeasi | IpoptStatus( | N/A | 121 | 0 | 5.7ms | 353us | 0.1x | BOTH_FAIL |
| PALMER1NE | 4 | 31 | LocalInfeasi | IpoptStatus( | N/A | 25 | 0 | 293us | 292us | 1.0x | BOTH_FAIL |
| PALMER2 | 4 | 0 | Optimal | Optimal | 2.49e-16 | 24 | 28 | 43us | 6.2ms | 144.9x | PASS |
| PALMER2A | 6 | 0 | Optimal | Optimal | 1.78e-15 | 142 | 91 | 333us | 19.8ms | 59.5x | PASS |
| PALMER2ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 117 | 0 | 1.0ms | 322us | 0.3x | BOTH_FAIL |
| PALMER2B | 4 | 0 | Acceptable | Optimal | 3.44e-15 | 37 | 15 | 78us | 3.4ms | 43.1x | PASS |
| PALMER2BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 579us | 286us | 0.5x | BOTH_FAIL |
| PALMER2C | 8 | 0 | Acceptable | Optimal | 1.30e-14 | 965 | 1 | 1.2ms | 441us | 0.4x | PASS |
| PALMER2E | 8 | 0 | Acceptable | Optimal | 2.09e-12 | 633 | 114 | 1.0ms | 24.3ms | 23.9x | PASS |
| PALMER2ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 109 | 0 | 3.2ms | 350us | 0.1x | BOTH_FAIL |
| PALMER2NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 261us | 287us | 1.1x | BOTH_FAIL |
| PALMER3 | 4 | 0 | Acceptable | Optimal | 6.25e-02 | 20 | 44 | 75.8ms | 7.8ms | 0.1x | MISMATCH |
| PALMER3A | 6 | 0 | Acceptable | Optimal | 1.72e-15 | 139 | 73 | 252us | 15.6ms | 62.0x | PASS |
| PALMER3ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 74 | 0 | 1.0ms | 276us | 0.3x | BOTH_FAIL |
| PALMER3B | 4 | 0 | Optimal | Optimal | 8.40e-16 | 31 | 15 | 50us | 3.4ms | 66.5x | PASS |
| PALMER3BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 312us | 314us | 1.0x | BOTH_FAIL |
| PALMER3C | 8 | 0 | Acceptable | Optimal | 3.99e-12 | 680 | 1 | 925us | 514us | 0.6x | PASS |
| PALMER3E | 8 | 0 | Acceptable | Optimal | 1.13e-13 | 304 | 32 | 549us | 6.0ms | 10.8x | PASS |
| PALMER3ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 340 | 0 | 2.0ms | 303us | 0.2x | BOTH_FAIL |
| PALMER3NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 217.6ms | 295us | 0.0x | BOTH_FAIL |
| PALMER4 | 4 | 0 | Optimal | Optimal | 5.72e-02 | 29 | 16 | 76.9ms | 3.4ms | 0.0x | MISMATCH |
| PALMER4A | 6 | 0 | Acceptable | Optimal | 3.23e-15 | 115 | 53 | 164us | 10.5ms | 63.9x | PASS |
| PALMER4ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 57 | 0 | 1.0ms | 285us | 0.3x | BOTH_FAIL |
| PALMER4B | 4 | 0 | Acceptable | Optimal | 1.43e-15 | 31 | 16 | 72us | 3.4ms | 47.1x | PASS |
| PALMER4BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 228us | 323us | 1.4x | BOTH_FAIL |
| PALMER4C | 8 | 0 | Acceptable | Optimal | 3.12e-11 | 607 | 1 | 815us | 460us | 0.6x | PASS |
| PALMER4E | 8 | 0 | Acceptable | Optimal | 3.97e-13 | 290 | 25 | 455us | 5.5ms | 12.1x | PASS |
| PALMER4ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 513 | 0 | 2.3ms | 288us | 0.1x | BOTH_FAIL |
| PALMER4NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 217.0ms | 308us | 0.0x | BOTH_FAIL |
| PALMER5A | 8 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 16.1ms | 657.1ms | 40.8x | ipopt_FAIL |
| PALMER5ANE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 24.1ms | 352us | 0.0x | BOTH_FAIL |
| PALMER5B | 9 | 0 | Acceptable | Optimal | 5.49e-03 | 1132 | 113 | 1.3ms | 22.3ms | 17.4x | MISMATCH |
| PALMER5BNE | 9 | 12 | LocalInfeasi | IpoptStatus( | N/A | 65 | 0 | 2.8ms | 314us | 0.1x | BOTH_FAIL |
| PALMER5C | 6 | 0 | Optimal | Optimal | 8.35e-16 | 13 | 1 | 10us | 482us | 45.9x | PASS |
| PALMER5D | 4 | 0 | Acceptable | Optimal | 4.23e-15 | 21 | 1 | 66us | 586us | 8.9x | PASS |
| PALMER5E | 8 | 0 | Acceptable | MaxIteration | N/A | 14 | 3000 | 3.1ms | 496.2ms | 158.5x | ipopt_FAIL |
| PALMER5ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 6.8ms | 283us | 0.0x | BOTH_FAIL |
| PALMER6A | 6 | 0 | Acceptable | Optimal | 3.48e-15 | 229 | 105 | 323us | 20.5ms | 63.3x | PASS |
| PALMER6ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 104 | 0 | 1.2ms | 361us | 0.3x | BOTH_FAIL |
| PALMER6C | 8 | 0 | Acceptable | Optimal | 3.49e-12 | 499 | 1 | 418us | 538us | 1.3x | PASS |
| PALMER6E | 8 | 0 | Acceptable | Optimal | 2.53e-11 | 255 | 30 | 287us | 6.1ms | 21.2x | PASS |
| PALMER6ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 399 | 0 | 1.2ms | 287us | 0.2x | BOTH_FAIL |
| PALMER7A | 6 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 12.5ms | 495.2ms | 39.5x | ipopt_FAIL |
| PALMER7ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 19.9ms | 292us | 0.0x | BOTH_FAIL |
| PALMER7C | 8 | 0 | Optimal | Optimal | 2.59e-13 | 2 | 1 | 35.6ms | 459us | 0.0x | PASS |
| PALMER7E | 8 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 16.2ms | 642.6ms | 39.7x | ipopt_FAIL |
| PALMER7ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 25.0ms | 305us | 0.0x | BOTH_FAIL |
| PALMER8A | 6 | 0 | Acceptable | Optimal | 2.22e-16 | 84 | 36 | 124us | 7.8ms | 62.4x | PASS |
| PALMER8ANE | 6 | 12 | LocalInfeasi | IpoptStatus( | N/A | 36 | 0 | 462us | 296us | 0.6x | BOTH_FAIL |
| PALMER8C | 8 | 0 | Acceptable | Optimal | 1.03e-14 | 525 | 1 | 444us | 607us | 1.4x | PASS |
| PALMER8E | 8 | 0 | Acceptable | Optimal | 9.31e-14 | 257 | 23 | 274us | 4.4ms | 16.1x | PASS |
| PALMER8ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 1.1ms | 285us | 0.3x | BOTH_FAIL |
| PARKCH | 15 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| PENTAGON | 6 | 15 | Optimal | Optimal | 1.47e-05 | 13 | 19 | 70us | 5.3ms | 75.5x | PASS |
| PFIT1 | 3 | 3 | Optimal | Infeasible | N/A | 257 | 266 | 659us | 62.6ms | 95.0x | ipopt_FAIL |
| PFIT1LS | 3 | 0 | Optimal | Optimal | 1.37e-20 | 264 | 263 | 181us | 52.7ms | 290.7x | PASS |
| PFIT2 | 3 | 3 | Optimal | RestorationF | N/A | 111 | 247 | 297us | 54.2ms | 182.5x | ipopt_FAIL |
| PFIT2LS | 3 | 0 | Optimal | Optimal | 1.90e-20 | 529 | 82 | 2.3ms | 27.6ms | 12.1x | PASS |
| PFIT3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 115 | 133 | 317us | 30.4ms | 96.1x | PASS |
| PFIT3LS | 3 | 0 | Optimal | Optimal | 2.58e-20 | 772 | 132 | 532us | 24.7ms | 46.4x | PASS |
| PFIT4 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 209 | 190 | 534us | 42.0ms | 78.8x | PASS |
| PFIT4LS | 3 | 0 | Optimal | Optimal | 6.66e-20 | 950 | 215 | 2.5ms | 43.9ms | 17.9x | PASS |
| POLAK1 | 3 | 2 | Optimal | Optimal | 3.67e-09 | 9 | 5 | 28us | 1.4ms | 49.9x | PASS |
| POLAK2 | 11 | 2 | Optimal | Optimal | 3.59e-11 | 28 | 10 | 129us | 2.6ms | 20.1x | PASS |
| POLAK3 | 12 | 10 | Optimal | MaxIteration | N/A | 15 | 3000 | 325us | 1.10s | 3388.7x | ipopt_FAIL |
| POLAK4 | 3 | 3 | Acceptable | Optimal | 4.53e-09 | 14 | 4 | 35us | 1.3ms | 36.0x | PASS |
| POLAK5 | 3 | 2 | Acceptable | Optimal | 1.88e-10 | 52 | 31 | 163us | 5.6ms | 34.3x | PASS |
| POLAK6 | 5 | 4 | Optimal | MaxIteration | N/A | 13 | 3000 | 52us | 1.65s | 31746.2x | ipopt_FAIL |
| PORTFL1 | 12 | 1 | Acceptable | Optimal | 1.28e-06 | 10 | 9 | 264us | 2.2ms | 8.5x | PASS |
| PORTFL2 | 12 | 1 | Acceptable | Optimal | 4.37e-07 | 10 | 8 | 218us | 2.1ms | 9.6x | PASS |
| PORTFL3 | 12 | 1 | Acceptable | Optimal | 1.50e-06 | 10 | 9 | 223us | 2.3ms | 10.3x | PASS |
| PORTFL4 | 12 | 1 | Acceptable | Optimal | 3.83e-06 | 9 | 8 | 201us | 2.1ms | 10.6x | PASS |
| PORTFL6 | 12 | 1 | Acceptable | Optimal | 4.96e-07 | 10 | 8 | 217us | 2.0ms | 9.4x | PASS |
| POWELLBS | 2 | 2 | Optimal | Optimal | 0.00e+00 | 89 | 11 | 167us | 1.9ms | 11.3x | PASS |
| POWELLBSLS | 2 | 0 | Optimal | Optimal | 6.41e-26 | 159 | 91 | 87us | 14.4ms | 166.1x | PASS |
| POWELLSQ | 2 | 2 | Acceptable | Infeasible | N/A | 9 | 29 | 27us | 6.0ms | 222.4x | ipopt_FAIL |
| POWELLSQLS | 2 | 0 | Optimal | Optimal | 1.53e-14 | 13 | 10 | 9us | 2.2ms | 245.6x | PASS |
| PRICE3NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 29us | 1.4ms | 48.2x | PASS |
| PRICE4 | 2 | 0 | Acceptable | Optimal | 1.85e-11 | 56 | 8 | 28us | 1.7ms | 60.0x | PASS |
| PRICE4B | 2 | 0 | Acceptable | Optimal | 2.19e-11 | 59 | 8 | 32us | 2.0ms | 62.3x | PASS |
| PRICE4NE | 2 | 2 | Optimal | Acceptable | 0.00e+00 | 10 | 23 | 27us | 4.0ms | 145.2x | PASS |
| PRODPL0 | 60 | 29 | Acceptable | Optimal | 1.33e-06 | 20 | 15 | 3.3ms | 4.2ms | 1.3x | PASS |
| PRODPL1 | 60 | 29 | Acceptable | Optimal | 9.46e-05 | 231 | 28 | 32.2ms | 7.6ms | 0.2x | PASS |
| PSPDOC | 4 | 0 | Optimal | Optimal | 3.50e-09 | 7 | 5 | 3.0ms | 1.5ms | 0.5x | PASS |
| PT | 2 | 501 | Optimal | Optimal | 8.94e-03 | 4 | 106 | 1.7ms | 104.5ms | 62.0x | MISMATCH |
| QC | 9 | 4 | Acceptable | Optimal | 2.42e-02 | 12 | 44 | 239us | 9.9ms | 41.4x | MISMATCH |
| QCNEW | 9 | 3 | Acceptable | Optimal | 2.48e-03 | 2999 | 6 | 63.9ms | 1.6ms | 0.0x | MISMATCH |
| QPCBLEND | 83 | 74 | Optimal | Optimal | 5.56e-07 | 37 | 19 | 3.7ms | 6.0ms | 1.6x | PASS |
| QPNBLEND | 83 | 74 | Optimal | Optimal | 4.77e-07 | 50 | 18 | 4.8ms | 5.6ms | 1.2x | PASS |
| RAT42 | 3 | 9 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 112us | 364us | 3.2x | BOTH_FAIL |
| RAT42LS | 3 | 0 | Optimal | Optimal | 8.82e-16 | 22 | 28 | 80us | 4.7ms | 59.2x | PASS |
| RAT43 | 4 | 15 | LocalInfeasi | IpoptStatus( | N/A | 3 | 0 | 134us | 369us | 2.8x | BOTH_FAIL |
| RAT43LS | 4 | 0 | Optimal | Optimal | 9.92e-01 | 3 | 34 | 9us | 6.2ms | 712.7x | MISMATCH |
| RECIPE | 3 | 3 | Acceptable | Optimal | 0.00e+00 | 19 | 16 | 103us | 3.0ms | 29.0x | PASS |
| RECIPELS | 3 | 0 | Acceptable | Optimal | 1.17e-11 | 32 | 29 | 17us | 5.3ms | 315.1x | PASS |
| RES | 20 | 14 | Acceptable | Optimal | 0.00e+00 | 12 | 10 | 183us | 2.2ms | 11.9x | PASS |
| RK23 | 17 | 11 | Acceptable | Optimal | 6.20e-05 | 17 | 10 | 214us | 2.9ms | 13.7x | PASS |
| ROBOT | 14 | 2 | Optimal | IpoptStatus( | N/A | 13 | 18 | 176us | 4.9ms | 28.0x | ipopt_FAIL |
| ROSENBR | 2 | 0 | Optimal | Optimal | 1.65e-21 | 34 | 21 | 17us | 3.9ms | 225.3x | PASS |
| ROSENBRTU | 2 | 0 | Optimal | Optimal | 4.56e-26 | 90 | 87 | 45us | 14.9ms | 329.5x | PASS |
| ROSENMMX | 5 | 4 | Optimal | Optimal | 2.27e-10 | 10 | 13 | 48us | 3.2ms | 67.8x | PASS |
| ROSZMAN1 | 4 | 25 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 615us | 337us | 0.5x | BOTH_FAIL |
| ROSZMAN1LS | 4 | 0 | Acceptable | Optimal | 3.90e-02 | 118 | 28 | 189us | 5.0ms | 26.3x | MISMATCH |
| RSNBRNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 38us | 645us | 17.1x | PASS |
| S268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 16us | 3.0ms | 185.0x | MISMATCH |
| S308 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 13 | 9 | 8us | 1.7ms | 210.4x | PASS |
| S308NE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 34us | 352us | 10.3x | BOTH_FAIL |
| S316-322 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 9 | 7 | 23us | 1.5ms | 65.4x | PASS |
| S365 | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 43.6ms | 1000us | 0.0x | BOTH_FAIL |
| S365MOD | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 46.3ms | 1.3ms | 0.0x | BOTH_FAIL |
| SANTA | 21 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 1.3ms | 355us | 0.3x | BOTH_FAIL |
| SANTALS | 21 | 0 | Acceptable | Optimal | 2.69e-09 | 141 | 31 | 441us | 7.8ms | 17.7x | PASS |
| SIM2BQP | 2 | 0 | Acceptable | Optimal | 2.25e-06 | 8 | 5 | 27us | 1.3ms | 47.4x | PASS |
| SIMBQP | 2 | 0 | Optimal | Optimal | 8.29e-09 | 6 | 5 | 221us | 1.3ms | 5.8x | PASS |
| SIMPLLPA | 2 | 2 | Optimal | Optimal | 1.18e-08 | 6 | 8 | 15us | 2.1ms | 139.7x | PASS |
| SIMPLLPB | 2 | 3 | Optimal | Optimal | 9.07e-09 | 3 | 7 | 20us | 1.8ms | 90.8x | PASS |
| SINEVAL | 2 | 0 | Optimal | Optimal | 1.03e-19 | 75 | 42 | 38us | 7.1ms | 186.1x | PASS |
| SINVALNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 40us | 643us | 16.1x | PASS |
| SIPOW1 | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 81 | 35.2ms | 250.9ms | 7.1x | MISMATCH |
| SIPOW1M | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 88 | 33.0ms | 241.9ms | 7.3x | MISMATCH |
| SIPOW2 | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 69 | 27.9ms | 180.5ms | 6.5x | MISMATCH |
| SIPOW2M | 2 | 2000 | Optimal | Optimal | 1.33e+00 | 1 | 73 | 33.1ms | 188.7ms | 5.7x | MISMATCH |
| SIPOW3 | 4 | 2000 | Optimal | Optimal | 7.10e-01 | 1 | 12 | 30.2ms | 37.3ms | 1.2x | MISMATCH |
| SIPOW4 | 4 | 2000 | Optimal | Optimal | 8.69e-01 | 1 | 11 | 43.5ms | 34.6ms | 0.8x | MISMATCH |
| SISSER | 2 | 0 | Acceptable | Optimal | 1.22e-10 | 15 | 18 | 8us | 2.5ms | 307.6x | PASS |
| SISSER2 | 2 | 0 | Acceptable | Optimal | 1.98e-11 | 21 | 20 | 10us | 2.7ms | 271.1x | PASS |
| SNAIL | 2 | 0 | Optimal | Optimal | 3.13e-21 | 10 | 63 | 7us | 9.4ms | 1285.0x | PASS |
| SNAKE | 2 | 2 | Optimal | Optimal | 2.00e-04 | 5 | 8 | 23us | 1.9ms | 82.8x | MISMATCH |
| SPANHYD | 97 | 33 | Acceptable | Optimal | 5.73e-13 | 70 | 20 | 7.9ms | 5.4ms | 0.7x | PASS |
| SPIRAL | 3 | 2 | Acceptable | Infeasible | N/A | 116 | 370 | 288us | 60.7ms | 211.0x | ipopt_FAIL |
| SSI | 3 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 4.7ms | 438.4ms | 92.3x | ipopt_FAIL |
| SSINE | 3 | 2 | Acceptable | Optimal | 0.00e+00 | 2999 | 224 | 30.2ms | 32.9ms | 1.1x | PASS |
| STANCMIN | 3 | 2 | Optimal | Optimal | 9.32e-09 | 10 | 9 | 25us | 1.9ms | 76.8x | PASS |
| STRATEC | 10 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| STREG | 4 | 0 | Optimal | Optimal | 8.90e-02 | 71 | 13 | 39us | 2.2ms | 57.4x | MISMATCH |
| STREGNE | 4 | 2 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 12us | 623us | 51.2x | PASS |
| SUPERSIM | 2 | 2 | Optimal | Optimal | 2.22e-16 | 7 | 1 | 20us | 569us | 29.1x | PASS |
| SWOPF | 83 | 92 | Optimal | Optimal | 9.19e-10 | 15 | 13 | 1.2ms | 4.1ms | 3.4x | PASS |
| SYNTHES1 | 6 | 6 | Acceptable | Optimal | 1.08e-05 | 13 | 8 | 44us | 1.8ms | 40.1x | PASS |
| SYNTHES2 | 11 | 14 | Optimal | Optimal | 8.24e-07 | 21 | 14 | 151us | 2.9ms | 19.1x | PASS |
| SYNTHES3 | 17 | 23 | Optimal | Optimal | 3.99e-08 | 36 | 13 | 613us | 2.7ms | 4.3x | PASS |
| TAME | 2 | 1 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 13us | 1.1ms | 82.6x | PASS |
| TAX13322 | 72 | 1261 | Optimal | N/A | N/A | 103 | 0 | 153.3ms | N/A | N/A | ipopt_FAIL |
| TAXR13322 | 72 | 1261 | Acceptable | Acceptable | 9.29e-01 | 141 | 56 | 654.8ms | 2.77s | 4.2x | MISMATCH |
| TENBARS1 | 18 | 9 | Acceptable | Optimal | 1.11e-08 | 181 | 39 | 2.7ms | 7.4ms | 2.8x | PASS |
| TENBARS2 | 18 | 8 | Optimal | Optimal | 1.00e-08 | 27 | 33 | 256us | 6.3ms | 24.8x | PASS |
| TENBARS3 | 18 | 8 | Optimal | Optimal | 1.01e-08 | 25 | 34 | 261us | 6.6ms | 25.2x | PASS |
| TENBARS4 | 18 | 9 | Optimal | Optimal | 1.01e-10 | 14 | 14 | 158us | 3.1ms | 19.6x | PASS |
| TFI1 | 3 | 101 | Optimal | Optimal | 1.26e-10 | 112 | 19 | 3.2ms | 6.7ms | 2.1x | PASS |
| TFI2 | 3 | 101 | Optimal | Optimal | 5.14e-04 | 83 | 8 | 1.8ms | 2.8ms | 1.5x | MISMATCH |
| TFI3 | 3 | 101 | Optimal | Optimal | 6.13e-09 | 89 | 13 | 1.8ms | 4.1ms | 2.3x | PASS |
| THURBER | 7 | 37 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 414.0ms | 301us | 0.0x | BOTH_FAIL |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 1.02e-14 | 2999 | 19 | 164.2ms | 3.3ms | 0.0x | PASS |
| TOINTGOR | 50 | 0 | Acceptable | Optimal | 5.19e-13 | 105 | 7 | 267us | 1.3ms | 5.0x | PASS |
| TOINTPSP | 50 | 0 | Acceptable | Optimal | 1.90e-14 | 91 | 20 | 193us | 4.5ms | 23.5x | PASS |
| TOINTQOR | 50 | 0 | Acceptable | Optimal | 9.67e-16 | 37 | 1 | 73us | 497us | 6.8x | PASS |
| TRIGGER | 7 | 6 | Acceptable | Optimal | 0.00e+00 | 15 | 15 | 52us | 2.5ms | 47.4x | PASS |
| TRO3X3 | 30 | 13 | Acceptable | Optimal | 3.61e-03 | 228 | 47 | 12.3ms | 9.8ms | 0.8x | MISMATCH |
| TRO4X4 | 63 | 25 | Acceptable | IpoptStatus( | N/A | 2999 | 157 | 327.1ms | 43.3ms | 0.1x | ipopt_FAIL |
| TRO6X2 | 45 | 21 | Acceptable | RestorationF | N/A | 3048 | 353 | 363.1ms | 92.1ms | 0.3x | ipopt_FAIL |
| TRUSPYR1 | 11 | 4 | Optimal | Optimal | 2.03e-08 | 30 | 10 | 225us | 2.1ms | 9.2x | PASS |
| TRUSPYR2 | 11 | 11 | Acceptable | Optimal | 5.18e-08 | 141 | 13 | 2.8ms | 2.9ms | 1.1x | PASS |
| TRY-B | 2 | 1 | Optimal | Optimal | 5.01e-19 | 13 | 23 | 26us | 4.5ms | 173.5x | PASS |
| TWOBARS | 2 | 2 | Optimal | Optimal | 7.00e-01 | 19 | 8 | 57us | 1.7ms | 29.2x | MISMATCH |
| VESUVIA | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 15.29s | 487us | 0.0x | BOTH_FAIL |
| VESUVIALS | 8 | 0 | Acceptable | Optimal | 3.39e-01 | 330 | 48 | 99.3ms | 18.0ms | 0.2x | MISMATCH |
| VESUVIO | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 105 | 0 | 1.71s | 508us | 0.0x | BOTH_FAIL |
| VESUVIOLS | 8 | 0 | Acceptable | Optimal | 2.41e-15 | 2999 | 10 | 1.13s | 4.6ms | 0.0x | PASS |
| VESUVIOU | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 141 | 0 | 62.1ms | 468us | 0.0x | BOTH_FAIL |
| VESUVIOULS | 8 | 0 | Acceptable | Optimal | 5.55e-16 | 22 | 8 | 21.0ms | 3.5ms | 0.2x | PASS |
| VIBRBEAM | 8 | 0 | Optimal | Optimal | 9.48e-01 | 115 | 58 | 8.1ms | 9.7ms | 1.2x | MISMATCH |
| VIBRBEAMNE | 8 | 30 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 12.4ms | 412us | 0.0x | BOTH_FAIL |
| WACHBIEG | 3 | 2 | Optimal | Infeasible | N/A | 12 | 15 | 103us | 3.3ms | 32.1x | ipopt_FAIL |
| WATER | 31 | 10 | Acceptable | Optimal | 9.14e-08 | 34 | 17 | 863us | 3.6ms | 4.1x | PASS |
| WAYSEA1 | 2 | 0 | Optimal | Optimal | 2.69e-15 | 34 | 14 | 19us | 1.9ms | 103.4x | PASS |
| WAYSEA1B | 2 | 0 | Optimal | Optimal | 4.54e-16 | 31 | 14 | 18us | 2.6ms | 142.6x | PASS |
| WAYSEA1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 15 | 7 | 29us | 1.1ms | 38.4x | PASS |
| WAYSEA2 | 2 | 0 | Optimal | Optimal | 9.84e-18 | 19 | 22 | 12us | 2.9ms | 249.7x | PASS |
| WAYSEA2B | 2 | 0 | Optimal | Optimal | 9.15e-18 | 17 | 22 | 10us | 3.6ms | 347.2x | PASS |
| WAYSEA2NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 22 | 11 | 37us | 1.8ms | 46.9x | PASS |
| WEEDS | 3 | 0 | Acceptable | Optimal | 1.99e-14 | 225 | 28 | 412us | 6.0ms | 14.6x | PASS |
| WEEDSNE | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 70 | 0 | 565us | 345us | 0.6x | BOTH_FAIL |
| WOMFLET | 3 | 3 | Acceptable | Optimal | 1.00e+00 | 14 | 8 | 44us | 1.8ms | 40.9x | MISMATCH |
| YFIT | 3 | 0 | Optimal | Optimal | 1.33e-19 | 77 | 36 | 73us | 6.9ms | 95.0x | PASS |
| YFITNE | 3 | 17 | Acceptable | IpoptStatus( | N/A | 36 | 0 | 194us | 320us | 1.6x | ipopt_FAIL |
| YFITU | 3 | 0 | Optimal | Optimal | 5.31e-21 | 77 | 36 | 76us | 5.3ms | 69.6x | PASS |
| ZANGWIL2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 2 | 1 | 1us | 486us | 486.2x | PASS |
| ZANGWIL3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 7us | 511us | 68.8x | PASS |
| ZECEVIC2 | 2 | 2 | Optimal | Optimal | 5.57e-10 | 11 | 8 | 26us | 1.8ms | 67.0x | PASS |
| ZECEVIC3 | 2 | 2 | Optimal | Optimal | 8.24e-10 | 12 | 17 | 42us | 3.0ms | 72.5x | PASS |
| ZECEVIC4 | 2 | 2 | Optimal | Optimal | 2.62e-09 | 11 | 10 | 25us | 2.1ms | 85.1x | PASS |
| ZY2 | 3 | 2 | Acceptable | Optimal | 4.67e-05 | 13 | 14 | 35us | 2.9ms | 84.1x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 174.0 | 45.0 |
| Median | 18 | 13 |
| Total  | 96418 | 24931 |

- ripopt fewer iterations: 128/554
- Ipopt fewer iterations: 374/554
- Tied: 52/554

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 35.4ms | 16.8ms |
| Median | 86us | 2.6ms |
| Total  | 19.60s | 9.29s |

- Geometric mean speedup (Ipopt_time/ripopt_time): **14.38x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 472/554 problems
- Ipopt faster: 82/554 problems
- Overall speedup (total time): 0.47x

## Failure Analysis

### Problems where only ripopt fails (1)

| Problem | n | m | ripopt status | Ipopt obj |
|---------|---|---|---------------|-----------|
| HEART6 | 6 | 6 | NumericalError | 0.000000e+00 |

### Problems where only Ipopt fails (43)

| Problem | n | m | Ipopt status | ripopt obj |
|---------|---|---|--------------|------------|
| ARGAUSS | 3 | 15 | IpoptStatus(-10) | 0.000000e+00 |
| AVION2 | 49 | 15 | MaxIterations | 9.468013e+07 |
| BEALENE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| BOX3NE | 3 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| BROWNBSNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| CRESC100 | 6 | 200 | Infeasible | 7.443593e-01 |
| DECONVB | 63 | 0 | MaxIterations | 1.290455e-08 |
| DENSCHNBNE | 2 | 3 | IpoptStatus(-10) | 0.000000e+00 |
| DENSCHNENE | 3 | 3 | Infeasible | 0.000000e+00 |
| DEVGLA1NE | 4 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| ENGVAL2NE | 3 | 5 | IpoptStatus(-10) | 0.000000e+00 |
| EQC | 9 | 3 | ErrorInStepComputation | -8.274326e+02 |
| EXP2NE | 2 | 10 | IpoptStatus(-10) | 0.000000e+00 |
| GROUPING | 100 | 125 | IpoptStatus(-10) | 1.385040e+01 |
| GULFNE | 3 | 99 | IpoptStatus(-10) | 0.000000e+00 |
| HIMMELBJ | 45 | 14 | ErrorInStepComputation | -1.910345e+03 |
| HS87 | 6 | 4 | MaxIterations | 8.996858e+03 |
| LANCZOS1 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS2 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LANCZOS3 | 6 | 24 | IpoptStatus(-10) | 0.000000e+00 |
| LEVYMONE6 | 3 | 6 | IpoptStatus(-10) | 0.000000e+00 |
| LEWISPOL | 6 | 9 | IpoptStatus(-10) | 1.212776e+00 |
| LHAIFAM | 99 | 150 | InvalidNumberDetected | 6.931472e-01 |
| MESH | 41 | 48 | IpoptStatus(4) | -1.798364e+08 |
| NYSTROM5 | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5A | 8 | 0 | MaxIterations | 3.814524e-02 |
| PALMER5E | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | 6 | 0 | MaxIterations | 1.033491e+01 |
| PALMER7E | 8 | 0 | MaxIterations | 6.697960e+00 |
| PFIT1 | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | 3 | 2 | Infeasible | 2.322254e-12 |
| SSI | 3 | 0 | MaxIterations | 1.376194e-09 |
| TAX13322 | 72 | 1261 | N/A | -3.726455e+03 |
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
| FBRAIN2NE | 4 | 2211 | Timeout | Timeout |
| FBRAIN3 | 6 | 2211 | Timeout | Timeout |
| FBRAIN3LS | 6 | 0 | Timeout | Timeout |
| FBRAINNE | 2 | 2211 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS1 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS2 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GAUSS3 | 8 | 250 | LocalInfeasibility | IpoptStatus(-10) |
| GBRAIN | 2 | 2200 | LocalInfeasibility | IpoptStatus(-10) |
| GROWTH | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| HAHN1 | 7 | 236 | Timeout | Timeout |
| HATFLDBNE | 4 | 4 | MaxIterations | Infeasible |
| HATFLDDNE | 3 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| HATFLDENE | 3 | 21 | LocalInfeasibility | IpoptStatus(-10) |
| HIMMELBD | 2 | 2 | RestorationFailed | Infeasible |
| HIMMELBFNE | 4 | 7 | LocalInfeasibility | IpoptStatus(-10) |
| HS25NE | 3 | 99 | LocalInfeasibility | IpoptStatus(-10) |
| HS2NE | 2 | 2 | MaxIterations | Infeasible |
| JENSMPNE | 2 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| JUDGENE | 2 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| KIRBY2 | 5 | 151 | LocalInfeasibility | IpoptStatus(-10) |
| KOEBHELBNE | 3 | 156 | LocalInfeasibility | IpoptStatus(-10) |
| KOWOSBNE | 4 | 11 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE10 | 10 | 20 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE5 | 2 | 4 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE7 | 4 | 8 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE8 | 5 | 10 | LocalInfeasibility | IpoptStatus(-10) |
| LEVYMONE9 | 8 | 16 | LocalInfeasibility | IpoptStatus(-10) |
| LRCOVTYPE | 54 | 0 | Timeout | Timeout |
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
| OET7 | 7 | 1002 | Timeout | Timeout |
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
| PALMER5ANE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER5BNE | 9 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER5ENE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER6ANE | 6 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER6ENE | 8 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER7ANE | 6 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER7ENE | 8 | 13 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER8ANE | 6 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PALMER8ENE | 8 | 12 | LocalInfeasibility | IpoptStatus(-10) |
| PARKCH | 15 | 0 | Timeout | Timeout |
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

### Objective mismatches (104)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 55
- **Convergence gap** (one Acceptable): 49
- **Better objective found by**: ripopt 24, Ipopt 80

| Problem | ripopt obj | Ipopt obj | Rel Diff | r_status | i_status | Better |
|---------|-----------|-----------|----------|----------|----------|--------|
| MAKELA1 | 1.414214e+00 | -1.414214e+00 | 2.00e+00 | Optimal | Optimal | ipopt |
| SIPOW1 | 4.338956e-01 | -1.000000e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW1M | 4.338960e-01 | -1.000001e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW2M | 3.310254e-01 | -1.000005e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| SIPOW2 | 3.310227e-01 | -1.000000e+00 | 1.33e+00 | Optimal | Optimal | ipopt |
| LEVYMONT5 | 7.773811e+00 | 1.239502e-25 | 1.00e+00 | Acceptable | Optimal | ipopt |
| LEVYMONT6 | 1.073629e-22 | 1.250286e+01 | 1.00e+00 | Optimal | Optimal | ripopt |
| WOMFLET | 1.010458e-11 | 6.050000e+00 | 1.00e+00 | Acceptable | Optimal | ripopt |
| MGH10SLS | 1.417866e+09 | 8.794586e+01 | 1.00e+00 | Acceptable | Optimal | ipopt |
| HS268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| S268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| DANIWOODLS | 1.039178e+02 | 4.317308e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| ELATTAR | 5.266741e+05 | 7.420618e+01 | 1.00e+00 | Acceptable | Optimal | ipopt |
| GROWTHLS | 3.542149e+03 | 1.004041e+00 | 1.00e+00 | Optimal | Optimal | ipopt |
| RAT43LS | 1.076462e+06 | 8.786405e+03 | 9.92e-01 | Optimal | Optimal | ipopt |
| HS16 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 | Optimal | Optimal | ipopt |
| KSIP | 2.768715e+01 | 5.757979e-01 | 9.79e-01 | Optimal | Optimal | ipopt |
| LEVYMONT8 | 4.980481e+00 | 1.734956e+02 | 9.71e-01 | Optimal | Optimal | ripopt |
| ELATVIDU | 1.712780e+00 | 5.475112e+01 | 9.69e-01 | Optimal | Optimal | ripopt |
| ELATVIDUB | 1.712780e+00 | 5.475112e+01 | 9.69e-01 | Acceptable | Optimal | ripopt |
| OET2 | 2.000000e+00 | 8.715962e-02 | 9.56e-01 | Acceptable | Optimal | ipopt |
| MWRIGHT | 1.288383e+00 | 2.497881e+01 | 9.48e-01 | Optimal | Optimal | ripopt |
| VIBRBEAM | 6.346254e+00 | 3.322376e-01 | 9.48e-01 | Optimal | Optimal | ipopt |
| JENSMP | 2.020000e+03 | 1.243622e+02 | 9.38e-01 | Optimal | Optimal | ipopt |
| TAXR13322 | -4.830584e+03 | -3.429089e+02 | 9.29e-01 | Acceptable | Acceptable | ripopt |
| HALDMADS | 1.714635e-01 | 2.218282e+00 | 9.23e-01 | Acceptable | Optimal | ripopt |
| BT4 | -4.551055e+01 | -3.704768e+00 | 9.19e-01 | Optimal | Optimal | ripopt |
| BOXBODLS | 1.168009e+03 | 9.771500e+03 | 8.80e-01 | Acceptable | Optimal | ripopt |
| BIGGSC4 | -3.128034e+00 | -2.450000e+01 | 8.72e-01 | Acceptable | Optimal | ipopt |
| SIPOW4 | 2.080345e+00 | 2.723620e-01 | 8.69e-01 | Optimal | Optimal | ipopt |
| HIMMELP2 | -6.205394e+01 | -8.198044e+00 | 8.68e-01 | Optimal | Optimal | ripopt |
| HIMMELP3 | -7.913699e+00 | -5.901318e+01 | 8.66e-01 | Optimal | Optimal | ipopt |
| HS23 | 9.473230e+00 | 2.000000e+00 | 7.89e-01 | Optimal | Optimal | ipopt |
| HS54 | -1.566691e-01 | -9.080748e-01 | 7.51e-01 | Optimal | Optimal | ipopt |
| HIMMELP6 | -1.475339e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| HIMMELP5 | -1.475901e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| SIPOW3 | 1.846738e+00 | 5.346586e-01 | 7.10e-01 | Optimal | Optimal | ipopt |
| TWOBARS | 5.028823e+00 | 1.508652e+00 | 7.00e-01 | Optimal | Optimal | ipopt |
| ECKERLE4LS | 6.996961e-01 | 1.463589e-03 | 6.98e-01 | Acceptable | Optimal | ipopt |
| MAKELA4 | 6.877365e-01 | -9.600000e-09 | 6.88e-01 | Optimal | Optimal | ipopt |
| LEVYMONT9 | 5.367398e+01 | 1.661496e+02 | 6.77e-01 | Optimal | Optimal | ripopt |
| MISTAKE | -5.000000e-01 | -1.000000e+00 | 5.00e-01 | Acceptable | Optimal | ipopt |
| VESUVIALS | 1.500440e+03 | 9.914100e+02 | 3.39e-01 | Acceptable | Optimal | ipopt |
| HYDC20LS | 2.976453e-01 | 2.967522e-15 | 2.98e-01 | Acceptable | Optimal | ipopt |
| AVGASA | -3.383299e+00 | -4.631926e+00 | 2.70e-01 | Acceptable | Optimal | ipopt |
| HS109 | 6.841976e+03 | 5.362069e+03 | 2.16e-01 | Acceptable | Optimal | ipopt |
| JUDGEB | 2.048234e+01 | 1.608173e+01 | 2.15e-01 | Optimal | Optimal | ipopt |
| JUDGE | 2.048234e+01 | 1.608173e+01 | 2.15e-01 | Optimal | Optimal | ipopt |
| EG1 | -1.132801e+00 | -1.429307e+00 | 2.07e-01 | Optimal | Optimal | ipopt |
| LEVYMONT7 | 1.015401e+01 | 1.251076e+01 | 1.88e-01 | Optimal | Optimal | ripopt |
| DENSCHNC | 1.833617e-01 | 2.177679e-20 | 1.83e-01 | Optimal | Optimal | ipopt |
| HS70 | 1.877383e-01 | 7.498464e-03 | 1.80e-01 | Optimal | Optimal | ipopt |
| AVGASB | -3.717055e+00 | -4.483219e+00 | 1.71e-01 | Acceptable | Optimal | ipopt |
| MUONSINELS | 5.114204e+04 | 4.387412e+04 | 1.42e-01 | Acceptable | Optimal | ipopt |
| MAKELA2 | 8.244898e+00 | 7.200000e+00 | 1.27e-01 | Optimal | Optimal | ipopt |
| MSS1 | -1.600000e+01 | -1.400000e+01 | 1.25e-01 | Acceptable | Optimal | ripopt |
| STREG | 1.110999e-23 | 8.901950e-02 | 8.90e-02 | Optimal | Optimal | ripopt |
| ACOPR14 | 8.859569e+03 | 8.081526e+03 | 8.78e-02 | Acceptable | Optimal | ipopt |
| HAHN1LS | 3.086398e+01 | 3.338424e+01 | 7.55e-02 | Acceptable | Optimal | ripopt |
| PALMER3 | 2.265958e+03 | 2.416980e+03 | 6.25e-02 | Acceptable | Optimal | ripopt |
| PALMER4 | 2.424016e+03 | 2.285383e+03 | 5.72e-02 | Optimal | Optimal | ipopt |
| OET1 | 5.834999e-01 | 5.382431e-01 | 4.53e-02 | Optimal | Optimal | ipopt |
| ROSZMAN1LS | 3.951903e-02 | 4.948485e-04 | 3.90e-02 | Acceptable | Optimal | ipopt |
| DEVGLA2 | 2.501146e-02 | 6.672171e-19 | 2.50e-02 | Acceptable | Optimal | ipopt |
| QC | -9.333976e+02 | -9.565379e+02 | 2.42e-02 | Acceptable | Optimal | ipopt |
| GIGOMEZ2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CB2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CHACONN1 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| HS55 | 6.666667e+00 | 6.805833e+00 | 2.04e-02 | Optimal | Optimal | ripopt |
| MAXLIKA | 1.149346e+03 | 1.136307e+03 | 1.13e-02 | Acceptable | Optimal | ipopt |
| EXPFITC | 3.272073e-02 | 2.330262e-02 | 9.42e-03 | Optimal | Optimal | ipopt |
| PT | 1.873320e-01 | 1.783942e-01 | 8.94e-03 | Optimal | Optimal | ipopt |
| CLIFF | 1.997866e-01 | 2.072380e-01 | 7.45e-03 | Optimal | Optimal | ripopt |
| HS85 | -2.200340e+00 | -2.215605e+00 | 6.89e-03 | Acceptable | Optimal | ipopt |
| LRIJCNN1 | 2.730781e-01 | 2.671576e-01 | 5.92e-03 | Acceptable | Optimal | ipopt |
| BIGGS5 | 5.655657e-03 | 1.088200e-19 | 5.66e-03 | Acceptable | Optimal | ipopt |
| BIGGS6 | 5.655650e-03 | 2.316725e-22 | 5.66e-03 | Optimal | Optimal | ipopt |
| PALMER5B | 1.524120e-02 | 9.752496e-03 | 5.49e-03 | Acceptable | Optimal | ipopt |
| ANTWERP | 3.260473e+03 | 3.245241e+03 | 4.67e-03 | Acceptable | Optimal | ipopt |
| DGOSPEC | -9.887540e+02 | -9.933506e+02 | 4.63e-03 | Acceptable | Optimal | ipopt |
| TRO3X3 | 9.000000e+00 | 8.967478e+00 | 3.61e-03 | Acceptable | Optimal | ipopt |
| QCNEW | -8.045191e+02 | -8.065219e+02 | 2.48e-03 | Acceptable | Optimal | ipopt |
| DEGENLPB | -3.069162e+01 | -3.076401e+01 | 2.35e-03 | Acceptable | Optimal | ipopt |
| LIN | -1.960628e-02 | -1.757754e-02 | 2.03e-03 | Optimal | Optimal | ripopt |
| EXPFITB | 7.016641e-03 | 5.019367e-03 | 2.00e-03 | Optimal | Optimal | ipopt |
| HET-Z | 1.001858e+00 | 1.000000e+00 | 1.85e-03 | Optimal | Optimal | ipopt |
| DEGENLPA | 3.060434e+00 | 3.054881e+00 | 1.81e-03 | Acceptable | Optimal | ipopt |
| DECONVC | 3.957562e-03 | 2.569475e-03 | 1.39e-03 | Acceptable | Optimal | ipopt |
| HS59 | -6.743243e+00 | -6.749505e+00 | 9.28e-04 | Acceptable | Optimal | ipopt |
| HS13 | 9.938594e-01 | 9.945785e-01 | 7.19e-04 | Acceptable | Optimal | ripopt |
| MGH09LS | 1.019673e-03 | 3.075056e-04 | 7.12e-04 | Acceptable | Optimal | ipopt |
| HS95 | 1.630984e-02 | 1.561772e-02 | 6.92e-04 | Acceptable | Optimal | ipopt |
| HS96 | 1.630574e-02 | 1.561775e-02 | 6.88e-04 | Acceptable | Optimal | ipopt |
| OET3 | 5.052162e-03 | 4.505043e-03 | 5.47e-04 | Optimal | Optimal | ipopt |
| TFI2 | 6.495452e-01 | 6.490311e-01 | 5.14e-04 | Optimal | Optimal | ipopt |
| HS116 | 9.754427e+01 | 9.758747e+01 | 4.43e-04 | Acceptable | Optimal | ripopt |
| LAUNCH | 9.008011e+00 | 9.004902e+00 | 3.45e-04 | Acceptable | Optimal | ipopt |
| HS45 | 1.000314e+00 | 1.000000e+00 | 3.14e-04 | Acceptable | Optimal | ipopt |
| LSC2LS | 1.333749e+01 | 1.333439e+01 | 2.32e-04 | Acceptable | Optimal | ipopt |
| DENSCHND | 3.820110e-09 | 2.221899e-04 | 2.22e-04 | Acceptable | Optimal | ripopt |
| HS17 | 1.000201e+00 | 1.000000e+00 | 2.01e-04 | Optimal | Optimal | ipopt |
| SNAKE | -2.601155e-10 | -1.999999e-04 | 2.00e-04 | Optimal | Optimal | ipopt |
| HS98 | 3.136197e+00 | 3.135806e+00 | 1.25e-04 | Acceptable | Optimal | ipopt |
| HATFLDH | -2.450304e+01 | -2.450000e+01 | 1.24e-04 | Acceptable | Optimal | ripopt |

---
*Generated by cutest_suite/compare.py*