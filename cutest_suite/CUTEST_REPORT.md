Report written to cutest_suite/CUTEST_REPORT.md

Summary:
  Total: 727
  ripopt solved: 598/727
  Ipopt solved: 555/727
  Both solved: 555/727
  Matching (rel diff < 1e-4): 454/555
t solved**: 555/727 (76.3%)
- **Both solved**: 555/727
- **Matching solutions** (rel obj diff < 1e-4): 454/555

## Accuracy Statistics (where both solve)

Relative difference = |r_obj - i_obj| / max(|r_obj|, |i_obj|, 1.0).  
The 1.0 floor prevents near-zero objectives from inflating the metric.

**Matching solutions** (454 problems, rel diff < 1e-4):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 1.23e-06 |
| Median | 3.46e-12 |
| Max    | 9.46e-05 |

**All both-solved** (555 problems, including 101 mismatches):

| Metric | Rel Diff |
|--------|----------|
| Mean   | 7.09e-02 |
| Median | 1.01e-09 |
| Max    | 2.00e+00 |

## Category Breakdown

| Category | Total | ripopt | Ipopt | Both | Match |
|----------|-------|--------|-------|------|-------|
| constrained | 493 | 376 | 339 | 339 | 278 |
| unconstrained | 234 | 222 | 216 | 216 | 176 |

## Detailed Results

| Problem | n | m | ripopt | Ipopt | Obj Diff | r_iter | i_iter | r_time | i_time | Speedup | Status |
|---------|---|---|--------|-------|----------|--------|--------|--------|--------|---------|--------|
| 3PK | 30 | 0 | Optimal | Optimal | 1.42e-15 | 9 | 9 | 5.7ms | 1.7ms | 0.3x | PASS |
| ACOPP14 | 38 | 68 | Optimal | Optimal | 9.79e-10 | 16 | 9 | 5.9ms | 3.3ms | 0.6x | PASS |
| ACOPP30 | 72 | 142 | Optimal | Optimal | 3.46e-03 | 37 | 13 | 8.4ms | 6.1ms | 0.7x | MISMATCH |
| ACOPR14 | 38 | 82 | Acceptable | Optimal | 1.62e-04 | 2335 | 13 | 1.27s | 4.6ms | 0.0x | MISMATCH |
| ACOPR30 | 72 | 172 | Acceptable | Optimal | 1.46e-01 | 1047 | 221 | 3.33s | 120.6ms | 0.0x | MISMATCH |
| AIRCRFTA | 8 | 5 | Optimal | Optimal | 0.00e+00 | 5 | 3 | 72us | 1.8ms | 24.7x | PASS |
| AIRCRFTB | 8 | 0 | Optimal | Optimal | 2.63e-21 | 49 | 15 | 44us | 2.4ms | 55.1x | PASS |
| AIRPORT | 84 | 42 | Optimal | Optimal | 4.45e-09 | 14 | 13 | 4.3ms | 5.5ms | 1.3x | PASS |
| AKIVA | 2 | 0 | Optimal | Optimal | 2.88e-16 | 14 | 6 | 54us | 868us | 16.0x | PASS |
| ALLINIT | 4 | 0 | Optimal | Optimal | 0.00e+00 | 15 | 20 | 14us | 3.2ms | 228.8x | PASS |
| ALLINITA | 4 | 4 | Acceptable | Optimal | 1.79e-06 | 23 | 12 | 3.4ms | 2.4ms | 0.7x | PASS |
| ALLINITC | 4 | 1 | Acceptable | Optimal | 7.31e-06 | 25 | 17 | 4.5ms | 3.1ms | 0.7x | PASS |
| ALLINITU | 4 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 14 | 12us | 2.1ms | 167.4x | PASS |
| ALSOTAME | 2 | 1 | Optimal | Optimal | 1.47e-08 | 10 | 8 | 26us | 1.5ms | 59.5x | PASS |
| ANTWERP | 27 | 10 | Acceptable | Optimal | 4.18e-04 | 146 | 108 | 61.4ms | 23.1ms | 0.4x | MISMATCH |
| ARGAUSS | 3 | 15 | Acceptable | IpoptStatus( | N/A | 2 | 0 | 31us | 113us | 3.7x | ipopt_FAIL |
| AVGASA | 8 | 10 | Optimal | Optimal | 1.55e-08 | 23 | 9 | 4.6ms | 1.8ms | 0.4x | PASS |
| AVGASB | 8 | 10 | Optimal | Optimal | 1.38e-08 | 33 | 11 | 4.5ms | 2.1ms | 0.5x | PASS |
| AVION2 | 49 | 15 | Acceptable | MaxIteration | N/A | 26 | 3000 | 139.7ms | 691.5ms | 5.0x | ipopt_FAIL |
| BA-L1 | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 6 | 602us | 1.4ms | 2.3x | PASS |
| BA-L1LS | 57 | 0 | Optimal | Optimal | 7.62e-21 | 28 | 10 | 228us | 2.2ms | 9.6x | PASS |
| BA-L1SP | 57 | 12 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.5ms | 2.3ms | 1.5x | PASS |
| BA-L1SPLS | 57 | 0 | Optimal | Optimal | 6.48e-17 | 32 | 9 | 1.2ms | 4.0ms | 3.4x | PASS |
| BARD | 3 | 0 | Optimal | Optimal | 1.73e-17 | 24 | 8 | 19us | 1.2ms | 63.6x | PASS |
| BARDNE | 3 | 15 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 86us | 93us | 1.1x | BOTH_FAIL |
| BATCH | 48 | 73 | Acceptable | Optimal | 2.29e-05 | 201 | 29 | 212.0ms | 7.3ms | 0.0x | PASS |
| BEALE | 2 | 0 | Optimal | Optimal | 4.26e-18 | 19 | 8 | 14us | 1.5ms | 100.9x | PASS |
| BEALENE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 25us | 108us | 4.3x | ipopt_FAIL |
| BENNETT5 | 3 | 154 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 24.7ms | 104us | 0.0x | BOTH_FAIL |
| BENNETT5LS | 3 | 0 | Acceptable | Optimal | 1.74e-05 | 40 | 21 | 2.4ms | 4.2ms | 1.8x | PASS |
| BIGGS3 | 6 | 0 | Optimal | Optimal | 5.29e-24 | 19 | 9 | 27us | 1.8ms | 67.9x | PASS |
| BIGGS5 | 6 | 0 | Optimal | Optimal | 5.66e-03 | 50 | 20 | 67us | 3.5ms | 51.4x | MISMATCH |
| BIGGS6 | 6 | 0 | Optimal | Optimal | 5.66e-03 | 41 | 79 | 49us | 12.0ms | 246.2x | MISMATCH |
| BIGGS6NE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 3.8ms | 100us | 0.0x | BOTH_FAIL |
| BIGGSC4 | 4 | 7 | Acceptable | Optimal | 5.33e-08 | 3028 | 17 | 12.6ms | 3.0ms | 0.2x | PASS |
| BLEACHNG | 17 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| BOOTH | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 18us | 272us | 15.2x | PASS |
| BOX2 | 3 | 0 | Acceptable | Optimal | 6.33e-14 | 15 | 8 | 18us | 1.1ms | 63.4x | PASS |
| BOX3 | 3 | 0 | Optimal | Optimal | 3.53e-19 | 23 | 9 | 20us | 2.0ms | 99.4x | PASS |
| BOX3NE | 3 | 10 | Optimal | IpoptStatus( | N/A | 11 | 0 | 57us | 455us | 7.9x | ipopt_FAIL |
| BOXBOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 58us | 275us | 4.7x | BOTH_FAIL |
| BOXBODLS | 2 | 0 | Acceptable | Optimal | 8.80e-01 | 27 | 13 | 104us | 2.1ms | 20.6x | MISMATCH |
| BQP1VAR | 1 | 0 | Optimal | Optimal | 9.98e-09 | 7 | 5 | 27us | 1.2ms | 43.1x | PASS |
| BQPGABIM | 50 | 0 | Acceptable | Optimal | 4.34e-06 | 11 | 12 | 7.0ms | 2.8ms | 0.4x | PASS |
| BQPGASIM | 50 | 0 | Acceptable | Optimal | 6.13e-06 | 11 | 12 | 5.1ms | 2.3ms | 0.4x | PASS |
| BRANIN | 2 | 0 | Optimal | Optimal | 0.00e+00 | 9 | 7 | 7us | 1.3ms | 193.9x | PASS |
| BRKMCC | 2 | 0 | Acceptable | Optimal | 8.33e-17 | 8 | 3 | 33us | 520us | 15.7x | PASS |
| BROWNBS | 2 | 0 | Acceptable | Optimal | 0.00e+00 | 12 | 7 | 7us | 1.1ms | 143.4x | PASS |
| BROWNBSNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 8 | 0 | 24us | 95us | 3.9x | ipopt_FAIL |
| BROWNDEN | 4 | 0 | Optimal | Optimal | 8.48e-16 | 23 | 8 | 30us | 1.1ms | 35.6x | PASS |
| BROWNDENE | 4 | 20 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 148us | 108us | 0.7x | BOTH_FAIL |
| BT1 | 2 | 1 | Optimal | Optimal | 2.40e-09 | 17 | 7 | 60us | 1.3ms | 22.5x | PASS |
| BT10 | 2 | 2 | Optimal | Optimal | 2.79e-09 | 7 | 6 | 24us | 1.3ms | 51.6x | PASS |
| BT11 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 8 | 31us | 1.1ms | 35.2x | PASS |
| BT12 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 25us | 662us | 26.7x | PASS |
| BT13 | 5 | 1 | Acceptable | Optimal | 1.00e-08 | 27 | 24 | 2.4ms | 3.9ms | 1.6x | PASS |
| BT2 | 3 | 1 | Optimal | Optimal | 3.80e-12 | 11 | 12 | 31us | 1.6ms | 51.2x | PASS |
| BT3 | 5 | 3 | Optimal | Optimal | 1.22e-14 | 1 | 1 | 14us | 375us | 27.4x | PASS |
| BT4 | 3 | 2 | Optimal | Optimal | 9.19e-01 | 6 | 9 | 28us | 2.1ms | 75.7x | MISMATCH |
| BT5 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 7 | 7 | 40us | 2.9ms | 73.9x | PASS |
| BT6 | 5 | 2 | Optimal | Optimal | 3.46e-12 | 9 | 13 | 36us | 2.2ms | 60.2x | PASS |
| BT7 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 16 | 74us | 2.7ms | 36.8x | PASS |
| BT8 | 5 | 2 | Acceptable | Optimal | 3.73e-09 | 32 | 14 | 315us | 1.9ms | 6.1x | PASS |
| BT9 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 42us | 1.7ms | 41.5x | PASS |
| BURKEHAN | 1 | 1 | RestorationF | Infeasible | N/A | 229 | 11 | 159.9ms | 2.3ms | 0.0x | BOTH_FAIL |
| BYRDSPHR | 3 | 2 | Optimal | Optimal | 3.95e-10 | 23 | 12 | 90us | 1.9ms | 21.2x | PASS |
| CAMEL6 | 2 | 0 | Optimal | Optimal | 4.30e-16 | 12 | 8 | 8us | 1.6ms | 200.7x | PASS |
| CANTILVR | 5 | 1 | Optimal | Optimal | 3.33e-09 | 31 | 11 | 84us | 2.1ms | 25.4x | PASS |
| CB2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 9 | 8 | 28us | 1.5ms | 52.4x | MISMATCH |
| CB3 | 3 | 3 | Optimal | Optimal | 4.30e-09 | 8 | 8 | 25us | 1.5ms | 58.2x | PASS |
| CERI651A | 7 | 61 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 705.5ms | 108us | 0.0x | BOTH_FAIL |
| CERI651ALS | 7 | 0 | Acceptable | Optimal | 8.56e-08 | 283 | 95 | 364.3ms | 15.7ms | 0.0x | PASS |
| CERI651B | 7 | 66 | LocalInfeasi | IpoptStatus( | N/A | 83 | 0 | 10.9ms | 112us | 0.0x | BOTH_FAIL |
| CERI651BLS | 7 | 0 | Optimal | Optimal | 1.23e-08 | 89 | 56 | 390.3ms | 9.0ms | 0.0x | PASS |
| CERI651C | 7 | 56 | LocalInfeasi | IpoptStatus( | N/A | 471 | 0 | 11.2ms | 104us | 0.0x | BOTH_FAIL |
| CERI651CLS | 7 | 0 | Acceptable | Optimal | 7.27e-08 | 482 | 53 | 2.6ms | 7.5ms | 2.9x | PASS |
| CERI651D | 7 | 67 | LocalInfeasi | IpoptStatus( | N/A | 91 | 0 | 4.4ms | 152us | 0.0x | BOTH_FAIL |
| CERI651DLS | 7 | 0 | Acceptable | Optimal | 3.08e-10 | 98 | 60 | 2.7ms | 10.0ms | 3.7x | PASS |
| CERI651E | 7 | 64 | LocalInfeasi | IpoptStatus( | N/A | 52 | 0 | 4.3ms | 105us | 0.0x | BOTH_FAIL |
| CERI651ELS | 7 | 0 | Acceptable | Optimal | 1.69e-07 | 166 | 45 | 2.0ms | 6.2ms | 3.1x | PASS |
| CHACONN1 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 5 | 6 | 26us | 1.1ms | 44.2x | MISMATCH |
| CHACONN2 | 3 | 3 | Optimal | Optimal | 4.45e-09 | 7 | 6 | 25us | 1.2ms | 46.0x | PASS |
| CHWIRUT1 | 3 | 214 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 5.0ms | 113us | 0.0x | BOTH_FAIL |
| CHWIRUT1LS | 3 | 0 | Optimal | Optimal | 7.63e-16 | 30 | 6 | 258us | 1.2ms | 4.8x | PASS |
| CHWIRUT2 | 3 | 54 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 986us | 99us | 0.1x | BOTH_FAIL |
| CHWIRUT2LS | 3 | 0 | Acceptable | Optimal | 8.86e-16 | 22 | 6 | 171us | 1.2ms | 6.8x | PASS |
| CLIFF | 2 | 0 | Optimal | Optimal | 7.45e-03 | 10 | 23 | 23us | 3.2ms | 135.8x | MISMATCH |
| CLUSTER | 2 | 2 | Acceptable | Optimal | 0.00e+00 | 12 | 9 | 33us | 1.4ms | 41.7x | PASS |
| CLUSTERLS | 2 | 0 | Acceptable | Optimal | 1.96e-14 | 19 | 17 | 13us | 2.2ms | 165.1x | PASS |
| CONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 7.8ms | 1.4ms | 0.2x | PASS |
| CONGIGMZ | 3 | 5 | Optimal | Optimal | 3.19e-09 | 7 | 20 | 85us | 3.4ms | 40.2x | PASS |
| COOLHANS | 9 | 9 | Optimal | Optimal | 0.00e+00 | 22 | 9 | 184us | 1.3ms | 7.0x | PASS |
| COOLHANSLS | 9 | 0 | Acceptable | Optimal | 2.25e-08 | 173 | 25 | 190us | 3.7ms | 19.6x | PASS |
| CORE1 | 65 | 59 | Optimal | Optimal | 4.04e-09 | 24 | 33 | 51.7ms | 7.3ms | 0.1x | PASS |
| CRESC100 | 6 | 200 | Acceptable | Infeasible | N/A | 169 | 155 | 4.57s | 121.8ms | 0.0x | ipopt_FAIL |
| CRESC132 | 6 | 2654 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| CRESC4 | 6 | 8 | Optimal | Optimal | 1.35e-08 | 21 | 64 | 315us | 13.1ms | 41.6x | PASS |
| CRESC50 | 6 | 100 | Optimal | Optimal | 3.57e-08 | 153 | 194 | 212.7ms | 82.4ms | 0.4x | PASS |
| CSFI1 | 5 | 4 | Optimal | Optimal | 1.49e-08 | 18 | 11 | 63us | 2.4ms | 38.6x | PASS |
| CSFI2 | 5 | 4 | Optimal | Optimal | 1.50e-08 | 16 | 14 | 148us | 3.3ms | 22.3x | PASS |
| CUBE | 2 | 0 | Optimal | Optimal | 4.92e-25 | 34 | 27 | 17us | 3.9ms | 223.6x | PASS |
| CUBENE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 27 | 1 | 51us | 292us | 5.7x | PASS |
| DALLASS | 46 | 31 | Optimal | Optimal | 1.41e-07 | 29 | 22 | 3.6ms | 4.5ms | 1.3x | PASS |
| DANIWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 46us | 100us | 2.2x | BOTH_FAIL |
| DANIWOODLS | 2 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 10 | 2us | 1.7ms | 831.2x | MISMATCH |
| DANWOOD | 2 | 6 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 199us | 467us | 2.3x | BOTH_FAIL |
| DANWOODLS | 2 | 0 | Optimal | Optimal | 1.28e-16 | 19 | 11 | 21us | 1.8ms | 86.5x | PASS |
| DECONVB | 63 | 0 | Acceptable | MaxIteration | N/A | 306 | 3000 | 29.3ms | 757.8ms | 25.9x | ipopt_FAIL |
| DECONVBNE | 63 | 40 | Acceptable | Optimal | 0.00e+00 | 5397 | 505 | 3.38s | 169.5ms | 0.1x | PASS |
| DECONVC | 63 | 1 | Optimal | Optimal | 1.16e-03 | 35 | 31 | 3.7ms | 8.5ms | 2.3x | MISMATCH |
| DECONVNE | 63 | 40 | Optimal | Acceptable | 0.00e+00 | 2 | 26 | 567us | 27.0ms | 47.6x | PASS |
| DECONVU | 63 | 0 | Acceptable | Optimal | 2.37e-07 | 100 | 333 | 483us | 89.5ms | 185.3x | PASS |
| DEGENLPA | 20 | 15 | Acceptable | Optimal | 1.81e-03 | 40 | 18 | 10.4ms | 3.3ms | 0.3x | MISMATCH |
| DEGENLPB | 20 | 15 | Acceptable | Optimal | 2.35e-03 | 33 | 19 | 10.0ms | 3.5ms | 0.4x | MISMATCH |
| DEMBO7 | 16 | 20 | Acceptable | Optimal | 7.72e-06 | 239 | 45 | 37.7ms | 8.7ms | 0.2x | PASS |
| DEMYMALO | 3 | 3 | Optimal | Optimal | 2.96e-09 | 12 | 9 | 32us | 1.8ms | 56.5x | PASS |
| DENSCHNA | 2 | 0 | Optimal | Optimal | 9.98e-19 | 10 | 6 | 7us | 850us | 128.4x | PASS |
| DENSCHNB | 2 | 0 | Optimal | Optimal | 1.88e-19 | 10 | 7 | 6us | 1.1ms | 178.6x | PASS |
| DENSCHNBNE | 2 | 3 | Optimal | IpoptStatus( | N/A | 7 | 0 | 21us | 101us | 4.8x | ipopt_FAIL |
| DENSCHNC | 2 | 0 | Optimal | Optimal | 1.83e-01 | 17 | 10 | 13us | 1.2ms | 92.7x | MISMATCH |
| DENSCHNCNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 10 | 7 | 39us | 1.0ms | 26.3x | PASS |
| DENSCHND | 3 | 0 | Acceptable | Optimal | 2.22e-04 | 53 | 26 | 32us | 3.4ms | 106.1x | MISMATCH |
| DENSCHNDNE | 3 | 3 | Acceptable | Acceptable | 0.00e+00 | 19 | 22 | 119us | 3.0ms | 25.4x | PASS |
| DENSCHNE | 3 | 0 | Optimal | Optimal | 1.35e-17 | 28 | 14 | 15us | 2.4ms | 158.4x | PASS |
| DENSCHNENE | 3 | 3 | Optimal | Infeasible | N/A | 10 | 10 | 29us | 1.9ms | 65.6x | ipopt_FAIL |
| DENSCHNF | 2 | 0 | Optimal | Optimal | 4.55e-22 | 12 | 6 | 23us | 958us | 41.7x | PASS |
| DENSCHNFNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 19us | 866us | 45.4x | PASS |
| DEVGLA1 | 4 | 0 | Optimal | Optimal | 1.40e-22 | 46 | 23 | 135us | 3.5ms | 25.8x | PASS |
| DEVGLA1B | 4 | 0 | Optimal | Optimal | 1.11e-24 | 67 | 20 | 672us | 4.2ms | 6.3x | PASS |
| DEVGLA1NE | 4 | 24 | Optimal | IpoptStatus( | N/A | 16 | 0 | 211us | 119us | 0.6x | ipopt_FAIL |
| DEVGLA2 | 5 | 0 | Acceptable | Optimal | 2.50e-02 | 74 | 13 | 805us | 1.9ms | 2.4x | MISMATCH |
| DEVGLA2B | 5 | 0 | Acceptable | Optimal | 2.58e-07 | 14 | 24 | 254.2ms | 4.5ms | 0.0x | PASS |
| DEVGLA2NE | 5 | 16 | LocalInfeasi | IpoptStatus( | N/A | 101 | 0 | 770us | 97us | 0.1x | BOTH_FAIL |
| DGOSPEC | 3 | 0 | Acceptable | Optimal | 4.63e-03 | 10 | 27 | 23.6ms | 4.8ms | 0.2x | MISMATCH |
| DIAMON2D | 66 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON2DLS | 66 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3D | 99 | 4643 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIAMON3DLS | 99 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| DIPIGRI | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 102us | 1.9ms | 18.8x | PASS |
| DISC2 | 29 | 23 | Optimal | Optimal | 9.11e-10 | 30 | 24 | 3.6ms | 5.7ms | 1.6x | PASS |
| DISCS | 36 | 66 | Acceptable | Optimal | 1.13e-06 | 112 | 184 | 2.40s | 171.8ms | 0.1x | PASS |
| DIXCHLNG | 10 | 5 | Optimal | Optimal | 0.00e+00 | 10 | 10 | 97us | 1.8ms | 18.5x | PASS |
| DJTL | 2 | 0 | Acceptable | Acceptable | 1.22e-15 | 1959 | 1538 | 1.7ms | 369.9ms | 211.4x | PASS |
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
| DNIEPER | 61 | 24 | Acceptable | Optimal | 2.48e-07 | 1314 | 23 | 485.1ms | 4.7ms | 0.0x | PASS |
| DUAL1 | 85 | 1 | Acceptable | Optimal | 3.88e-06 | 12 | 15 | 258.5ms | 6.0ms | 0.0x | PASS |
| DUAL2 | 96 | 1 | Acceptable | Optimal | 8.07e-09 | 12 | 12 | 326.3ms | 6.4ms | 0.0x | PASS |
| DUAL4 | 75 | 1 | Acceptable | Optimal | 1.40e-07 | 11 | 12 | 193.5ms | 4.2ms | 0.0x | PASS |
| DUALC1 | 9 | 215 | Acceptable | Optimal | 6.76e-06 | 47 | 18 | 1.11s | 11.6ms | 0.0x | PASS |
| DUALC2 | 7 | 229 | Acceptable | Optimal | 1.62e-06 | 17 | 12 | 1.17s | 8.1ms | 0.0x | PASS |
| DUALC5 | 8 | 278 | Acceptable | Optimal | 1.83e-07 | 12 | 11 | 1.95s | 8.2ms | 0.0x | PASS |
| DUALC8 | 8 | 503 | Acceptable | Optimal | 1.77e-07 | 24 | 13 | 11.12s | 15.5ms | 0.0x | PASS |
| ECKERLE4 | 3 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 249us | 110us | 0.4x | BOTH_FAIL |
| ECKERLE4LS | 3 | 0 | Acceptable | Optimal | 6.98e-01 | 12 | 36 | 23us | 5.7ms | 245.2x | MISMATCH |
| EG1 | 3 | 0 | Optimal | Optimal | 2.07e-01 | 9 | 8 | 30.3ms | 1.7ms | 0.1x | MISMATCH |
| EGGCRATE | 2 | 0 | Optimal | Optimal | 5.62e-16 | 10 | 5 | 8us | 967us | 123.4x | PASS |
| EGGCRATEB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 12 | 6 | 18us | 1.3ms | 72.1x | PASS |
| EGGCRATENE | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 42us | 95us | 2.2x | BOTH_FAIL |
| ELATTAR | 7 | 102 | Acceptable | Optimal | 1.00e+00 | 1098 | 81 | 395.3ms | 35.2ms | 0.1x | MISMATCH |
| ELATVIDU | 2 | 0 | Optimal | Optimal | 9.69e-01 | 14 | 11 | 13us | 1.7ms | 130.7x | MISMATCH |
| ELATVIDUB | 2 | 0 | Acceptable | Optimal | 9.69e-01 | 18 | 11 | 38us | 1.8ms | 48.5x | MISMATCH |
| ELATVIDUNE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 48us | 97us | 2.0x | BOTH_FAIL |
| ENGVAL2 | 3 | 0 | Optimal | Optimal | 1.70e-20 | 31 | 21 | 18us | 2.9ms | 159.3x | PASS |
| ENGVAL2NE | 3 | 5 | Optimal | IpoptStatus( | N/A | 17 | 0 | 57us | 97us | 1.7x | ipopt_FAIL |
| ENSO | 9 | 168 | LocalInfeasi | IpoptStatus( | N/A | 31 | 0 | 8.6ms | 164us | 0.0x | BOTH_FAIL |
| ENSOLS | 9 | 0 | Acceptable | Optimal | 4.33e-16 | 27 | 7 | 1.8ms | 2.0ms | 1.1x | PASS |
| EQC | 9 | 3 | Acceptable | ErrorInStepC | N/A | 3998 | 15 | 21.6ms | 4.0ms | 0.2x | ipopt_FAIL |
| ERRINBAR | 18 | 9 | Acceptable | Optimal | 4.70e-07 | 68 | 37 | 1.0ms | 7.5ms | 7.2x | PASS |
| EXP2 | 2 | 0 | Optimal | Optimal | 1.58e-20 | 13 | 7 | 12us | 1.3ms | 112.8x | PASS |
| EXP2B | 2 | 0 | Optimal | Optimal | 3.41e-21 | 13 | 7 | 12us | 1.3ms | 102.1x | PASS |
| EXP2NE | 2 | 10 | Optimal | IpoptStatus( | N/A | 7 | 0 | 39us | 97us | 2.5x | ipopt_FAIL |
| EXPFIT | 2 | 0 | Optimal | Optimal | 1.39e-16 | 14 | 8 | 15us | 1.3ms | 88.5x | PASS |
| EXPFITA | 5 | 22 | Optimal | Optimal | 3.40e-08 | 40 | 13 | 269us | 2.6ms | 9.6x | PASS |
| EXPFITB | 5 | 102 | Optimal | Optimal | 2.00e-03 | 189 | 16 | 7.8ms | 5.2ms | 0.7x | MISMATCH |
| EXPFITC | 5 | 502 | Optimal | Optimal | 9.42e-03 | 17 | 18 | 7.6ms | 17.0ms | 2.2x | MISMATCH |
| EXPFITNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 5 | 0 | 61us | 94us | 1.5x | BOTH_FAIL |
| EXTRASIM | 2 | 1 | Optimal | Optimal | 1.82e-08 | 4 | 3 | 13us | 647us | 48.2x | PASS |
| FBRAIN | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 28.8ms | 272us | 0.0x | BOTH_FAIL |
| FBRAIN2 | 4 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 318.1ms | 456us | 0.0x | BOTH_FAIL |
| FBRAIN2LS | 4 | 0 | Optimal | Optimal | 5.11e-10 | 10 | 10 | 199.3ms | 9.2ms | 0.0x | PASS |
| FBRAIN2NE | 4 | 2211 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAIN3 | 6 | 2211 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAIN3LS | 6 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| FBRAINLS | 2 | 0 | Optimal | Optimal | 7.77e-16 | 14 | 7 | 3.8ms | 3.6ms | 0.9x | PASS |
| FBRAINNE | 2 | 2211 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 52.7ms | 265us | 0.0x | BOTH_FAIL |
| FCCU | 19 | 8 | Optimal | Optimal | 1.75e-15 | 10 | 9 | 94us | 1.8ms | 18.9x | PASS |
| FEEDLOC | 90 | 259 | Optimal | Optimal | 1.10e-08 | 75 | 23 | 53.1ms | 15.8ms | 0.3x | PASS |
| FLETCHER | 4 | 4 | Optimal | Optimal | 9.68e-09 | 38 | 28 | 435us | 5.2ms | 11.9x | PASS |
| FLT | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 21us | 922us | 43.6x | PASS |
| GAUSS1 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 3.76s | 174us | 0.0x | BOTH_FAIL |
| GAUSS1LS | 8 | 0 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 1.33s | 1.2ms | 0.0x | PASS |
| GAUSS2 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 3.24s | 136us | 0.0x | BOTH_FAIL |
| GAUSS2LS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 10 | 5 | 1.28s | 1.0ms | 0.0x | PASS |
| GAUSS3 | 8 | 250 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 3.24s | 139us | 0.0x | BOTH_FAIL |
| GAUSS3LS | 8 | 0 | Optimal | Optimal | 3.65e-16 | 7 | 11 | 1.00s | 2.5ms | 0.0x | PASS |
| GAUSSIAN | 3 | 0 | Optimal | Optimal | 6.24e-18 | 3 | 2 | 6us | 411us | 68.0x | PASS |
| GBRAIN | 2 | 2200 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 35.9ms | 260us | 0.0x | BOTH_FAIL |
| GBRAINLS | 2 | 0 | Optimal | Optimal | 6.23e-16 | 12 | 6 | 4.7ms | 3.2ms | 0.7x | PASS |
| GENHS28 | 10 | 8 | Optimal | Optimal | 1.22e-15 | 1 | 1 | 32us | 388us | 12.3x | PASS |
| GIGOMEZ1 | 3 | 3 | Optimal | Optimal | 2.87e-09 | 10 | 13 | 34us | 2.3ms | 68.0x | PASS |
| GIGOMEZ2 | 3 | 3 | Optimal | Optimal | 2.39e-02 | 29 | 7 | 137us | 1.3ms | 9.8x | MISMATCH |
| GIGOMEZ3 | 3 | 3 | Optimal | Optimal | 4.08e-09 | 10 | 8 | 30us | 1.4ms | 48.2x | PASS |
| GOFFIN | 51 | 50 | Optimal | Optimal | 9.97e-09 | 4 | 7 | 1.5ms | 3.5ms | 2.3x | PASS |
| GOTTFR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 5 | 26us | 899us | 34.4x | PASS |
| GOULDQP1 | 32 | 17 | Acceptable | Optimal | 4.28e-07 | 34 | 15 | 60.2ms | 3.0ms | 0.1x | PASS |
| GROUPING | 100 | 125 | Acceptable | IpoptStatus( | N/A | 7 | 0 | 4.6ms | 104us | 0.0x | ipopt_FAIL |
| GROWTH | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 42 | 0 | 228us | 90us | 0.4x | BOTH_FAIL |
| GROWTHLS | 3 | 0 | Optimal | Optimal | 1.00e+00 | 1 | 71 | 7us | 11.3ms | 1713.0x | MISMATCH |
| GULF | 3 | 0 | Optimal | Optimal | 1.82e-20 | 50 | 28 | 625us | 5.0ms | 8.1x | PASS |
| GULFNE | 3 | 99 | Optimal | IpoptStatus( | N/A | 22 | 0 | 1.2ms | 110us | 0.1x | ipopt_FAIL |
| HAHN1 | 7 | 236 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| HAHN1LS | 7 | 0 | Acceptable | Optimal | 7.55e-02 | 557 | 78 | 49.8ms | 16.8ms | 0.3x | MISMATCH |
| HAIFAM | 99 | 150 | Acceptable | Optimal | 8.33e-07 | 27 | 40 | 7.40s | 15.3ms | 0.0x | PASS |
| HAIFAS | 13 | 9 | Optimal | Optimal | 2.50e-01 | 188 | 16 | 2.0ms | 3.3ms | 1.6x | MISMATCH |
| HAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 53 | 62 | 37us | 10.2ms | 277.6x | PASS |
| HALDMADS | 6 | 42 | Acceptable | Optimal | 9.79e-01 | 134 | 8 | 26.6ms | 4.0ms | 0.1x | MISMATCH |
| HART6 | 6 | 0 | Optimal | Optimal | 0.00e+00 | 18 | 7 | 14us | 1.5ms | 105.5x | PASS |
| HATFLDA | 4 | 0 | Optimal | Optimal | 1.78e-19 | 49 | 13 | 27us | 2.1ms | 78.5x | PASS |
| HATFLDANE | 4 | 4 | Acceptable | Optimal | 0.00e+00 | 9 | 6 | 28us | 1.2ms | 42.2x | PASS |
| HATFLDB | 4 | 0 | Optimal | Optimal | 3.90e-09 | 9 | 8 | 22.3ms | 1.4ms | 0.1x | PASS |
| HATFLDBNE | 4 | 4 | MaxIteration | Infeasible | N/A | 2999 | 13 | 1.07s | 2.6ms | 0.0x | BOTH_FAIL |
| HATFLDC | 25 | 0 | Acceptable | Optimal | 3.35e-14 | 25 | 5 | 35us | 1.1ms | 30.8x | PASS |
| HATFLDCNE | 25 | 25 | Acceptable | Optimal | 0.00e+00 | 9 | 4 | 129us | 950us | 7.4x | PASS |
| HATFLDD | 3 | 0 | Optimal | Optimal | 1.89e-07 | 27 | 21 | 24us | 2.8ms | 117.6x | PASS |
| HATFLDDNE | 3 | 10 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 150us | 98us | 0.7x | BOTH_FAIL |
| HATFLDE | 3 | 0 | Optimal | Optimal | 2.22e-06 | 32 | 20 | 48us | 2.5ms | 52.2x | PASS |
| HATFLDENE | 3 | 21 | LocalInfeasi | IpoptStatus( | N/A | 19 | 0 | 223us | 104us | 0.5x | BOTH_FAIL |
| HATFLDF | 3 | 3 | Optimal | Optimal | 0.00e+00 | 18 | 135 | 52us | 23.8ms | 456.5x | PASS |
| HATFLDFL | 3 | 0 | Optimal | Optimal | 1.03e-08 | 494 | 1281 | 10.0ms | 196.1ms | 19.6x | PASS |
| HATFLDFLNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 15 | 15 | 799us | 2.8ms | 3.5x | PASS |
| HATFLDFLS | 3 | 0 | Optimal | Optimal | 3.78e-18 | 73 | 36 | 43us | 5.3ms | 122.9x | PASS |
| HATFLDG | 25 | 25 | Optimal | Optimal | 0.00e+00 | 13 | 7 | 323us | 1.5ms | 4.6x | PASS |
| HATFLDGLS | 25 | 0 | Acceptable | Optimal | 2.19e-13 | 44 | 14 | 67us | 2.5ms | 36.4x | PASS |
| HATFLDH | 4 | 7 | Acceptable | Optimal | 1.24e-04 | 96 | 17 | 14.3ms | 2.9ms | 0.2x | MISMATCH |
| HEART6 | 6 | 6 | Optimal | Optimal | 0.00e+00 | 775 | 22 | 2.9ms | 4.9ms | 1.7x | PASS |
| HEART6LS | 6 | 0 | Optimal | Optimal | 9.36e-23 | 907 | 875 | 781us | 140.8ms | 180.2x | PASS |
| HEART8 | 8 | 8 | Optimal | Optimal | 0.00e+00 | 66 | 12 | 350us | 2.2ms | 6.4x | PASS |
| HEART8LS | 8 | 0 | Acceptable | Optimal | 6.00e-18 | 551 | 106 | 552us | 17.9ms | 32.3x | PASS |
| HELIX | 3 | 0 | Optimal | Optimal | 9.65e-19 | 25 | 13 | 15us | 2.2ms | 143.0x | PASS |
| HELIXNE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 12 | 7 | 38us | 1.2ms | 30.9x | PASS |
| HET-Z | 2 | 1002 | Optimal | Optimal | 1.85e-03 | 12 | 11 | 13.5ms | 20.6ms | 1.5x | MISMATCH |
| HIELOW | 3 | 0 | Optimal | Optimal | 5.46e-15 | 7 | 8 | 51.2ms | 11.0ms | 0.2x | PASS |
| HIMMELBA | 2 | 2 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 12us | 300us | 24.8x | PASS |
| HIMMELBB | 2 | 0 | Optimal | Optimal | 1.39e-17 | 12 | 18 | 12us | 3.6ms | 296.0x | PASS |
| HIMMELBC | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 6 | 45us | 1.4ms | 30.3x | PASS |
| HIMMELBCLS | 2 | 0 | Optimal | Optimal | 1.82e-23 | 8 | 6 | 5us | 922us | 187.4x | PASS |
| HIMMELBD | 2 | 2 | RestorationF | Infeasible | N/A | 14 | 22 | 5.7ms | 4.5ms | 0.8x | BOTH_FAIL |
| HIMMELBE | 3 | 3 | Optimal | Optimal | 0.00e+00 | 5 | 2 | 17us | 388us | 22.8x | PASS |
| HIMMELBF | 4 | 0 | Acceptable | Optimal | 4.10e-15 | 56 | 75 | 135us | 10.7ms | 79.0x | PASS |
| HIMMELBFNE | 4 | 7 | LocalInfeasi | IpoptStatus( | N/A | 37 | 0 | 477us | 96us | 0.2x | BOTH_FAIL |
| HIMMELBG | 2 | 0 | Optimal | Optimal | 3.53e-18 | 10 | 6 | 6us | 1.0ms | 175.5x | PASS |
| HIMMELBH | 2 | 0 | Optimal | Optimal | 4.44e-16 | 7 | 4 | 5us | 869us | 178.3x | PASS |
| HIMMELBI | 100 | 12 | Optimal | Optimal | 5.80e-10 | 24 | 13 | 899us | 3.2ms | 3.6x | PASS |
| HIMMELBJ | 45 | 14 | Acceptable | ErrorInStepC | N/A | 28 | 580 | 34.3ms | 135.9ms | 4.0x | ipopt_FAIL |
| HIMMELBK | 24 | 14 | Optimal | Optimal | 4.89e-08 | 19 | 18 | 594us | 3.4ms | 5.8x | PASS |
| HIMMELP1 | 2 | 0 | Acceptable | Optimal | 1.83e-15 | 21 | 10 | 55us | 2.0ms | 36.3x | PASS |
| HIMMELP2 | 2 | 1 | Optimal | Optimal | 8.68e-01 | 12 | 17 | 30us | 3.2ms | 108.3x | MISMATCH |
| HIMMELP3 | 2 | 2 | Optimal | Optimal | 8.66e-01 | 0 | 11 | 10us | 1.9ms | 198.1x | MISMATCH |
| HIMMELP4 | 2 | 3 | Optimal | Optimal | 1.23e-08 | 9 | 23 | 31us | 4.0ms | 128.0x | PASS |
| HIMMELP5 | 2 | 3 | Optimal | Optimal | 7.50e-01 | 9 | 46 | 32us | 8.2ms | 254.6x | MISMATCH |
| HIMMELP6 | 2 | 5 | Optimal | Optimal | 7.50e-01 | 11 | 31 | 42us | 6.5ms | 152.6x | MISMATCH |
| HONG | 4 | 1 | Optimal | Optimal | 4.72e-16 | 9 | 7 | 29us | 1.2ms | 42.4x | PASS |
| HS1 | 2 | 0 | Optimal | Optimal | 9.07e-21 | 25 | 28 | 15us | 5.1ms | 340.8x | PASS |
| HS10 | 2 | 1 | Optimal | Optimal | 4.99e-09 | 19 | 12 | 36us | 2.0ms | 56.3x | PASS |
| HS100 | 7 | 4 | Optimal | Optimal | 1.36e-11 | 11 | 9 | 85us | 1.9ms | 22.1x | PASS |
| HS100LNP | 7 | 2 | Optimal | Optimal | 1.67e-16 | 6 | 20 | 32us | 2.4ms | 76.3x | PASS |
| HS100MOD | 7 | 4 | Optimal | Optimal | 1.73e-11 | 9 | 14 | 59us | 2.5ms | 42.8x | PASS |
| HS101 | 7 | 5 | Optimal | Optimal | 2.36e-08 | 30 | 39 | 306us | 9.1ms | 29.8x | PASS |
| HS102 | 7 | 5 | Optimal | Optimal | 4.27e-08 | 24 | 52 | 258us | 11.4ms | 44.3x | PASS |
| HS103 | 7 | 5 | Optimal | Optimal | 4.13e-08 | 26 | 21 | 363us | 4.2ms | 11.5x | PASS |
| HS104 | 8 | 5 | Optimal | Optimal | 1.50e-08 | 24 | 8 | 212us | 1.6ms | 7.3x | PASS |
| HS105 | 8 | 1 | Optimal | Optimal | 9.48e-12 | 30 | 23 | 3.2ms | 6.4ms | 2.0x | PASS |
| HS106 | 8 | 6 | Optimal | Optimal | 1.74e-08 | 18 | 18 | 76us | 3.1ms | 40.3x | PASS |
| HS107 | 9 | 6 | Acceptable | Optimal | 1.29e-07 | 52 | 7 | 7.1ms | 1.4ms | 0.2x | PASS |
| HS108 | 9 | 13 | Optimal | Optimal | 3.66e-01 | 184 | 11 | 40.3ms | 2.3ms | 0.1x | MISMATCH |
| HS109 | 9 | 10 | Acceptable | Optimal | 5.47e-03 | 3511 | 14 | 27.5ms | 2.6ms | 0.1x | MISMATCH |
| HS11 | 2 | 1 | Optimal | Optimal | 3.59e-09 | 6 | 6 | 20us | 1.1ms | 56.6x | PASS |
| HS111 | 10 | 3 | Optimal | Optimal | 1.05e-11 | 11 | 15 | 118us | 2.9ms | 24.2x | PASS |
| HS111LNP | 10 | 3 | Optimal | Optimal | 1.03e-10 | 11 | 15 | 94us | 2.4ms | 25.3x | PASS |
| HS112 | 10 | 3 | Optimal | Optimal | 6.49e-14 | 10 | 10 | 95us | 1.9ms | 20.2x | PASS |
| HS113 | 10 | 8 | Optimal | Optimal | 1.68e-09 | 14 | 9 | 126us | 1.9ms | 14.9x | PASS |
| HS114 | 10 | 11 | Optimal | Optimal | 1.04e-07 | 19 | 13 | 138us | 2.3ms | 16.9x | PASS |
| HS116 | 13 | 14 | Acceptable | Optimal | 4.43e-04 | 131 | 19 | 17.6ms | 3.6ms | 0.2x | MISMATCH |
| HS117 | 15 | 5 | Acceptable | Optimal | 3.56e-07 | 38 | 19 | 9.3ms | 3.5ms | 0.4x | PASS |
| HS118 | 15 | 17 | Optimal | Optimal | 1.39e-09 | 236 | 10 | 4.2ms | 2.0ms | 0.5x | PASS |
| HS119 | 16 | 8 | Acceptable | Optimal | 8.77e-05 | 21 | 17 | 16.3ms | 3.3ms | 0.2x | PASS |
| HS12 | 2 | 1 | Optimal | Optimal | 1.66e-10 | 7 | 6 | 23us | 1.1ms | 49.8x | PASS |
| HS13 | 2 | 1 | Acceptable | Optimal | 1.21e-03 | 53 | 47 | 1.5ms | 7.7ms | 5.1x | MISMATCH |
| HS14 | 2 | 2 | Optimal | Optimal | 1.32e-08 | 5 | 5 | 17us | 1.1ms | 61.2x | PASS |
| HS15 | 2 | 2 | Optimal | Optimal | 8.08e-08 | 11 | 13 | 36us | 2.5ms | 69.5x | PASS |
| HS16 | 2 | 2 | Optimal | Optimal | 9.89e-01 | 10 | 10 | 38us | 1.9ms | 50.0x | MISMATCH |
| HS17 | 2 | 2 | Optimal | Optimal | 2.01e-04 | 12 | 22 | 31us | 3.8ms | 123.0x | MISMATCH |
| HS18 | 2 | 2 | Optimal | Optimal | 9.03e-10 | 14 | 10 | 35us | 1.7ms | 48.9x | PASS |
| HS19 | 2 | 2 | Optimal | Optimal | 3.09e-09 | 20 | 12 | 47us | 2.3ms | 48.2x | PASS |
| HS1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 25 | 30 | 48us | 6.3ms | 129.4x | PASS |
| HS2 | 2 | 0 | Optimal | Optimal | 5.48e-09 | 10 | 10 | 18.4ms | 2.1ms | 0.1x | PASS |
| HS20 | 2 | 3 | Optimal | Optimal | 6.17e-08 | 18 | 5 | 57us | 1.0ms | 18.3x | PASS |
| HS21 | 2 | 1 | Optimal | Optimal | 1.38e-10 | 9 | 6 | 21us | 1.1ms | 53.4x | PASS |
| HS21MOD | 7 | 1 | Acceptable | Optimal | 4.30e-09 | 23 | 13 | 1.8ms | 2.3ms | 1.3x | PASS |
| HS22 | 2 | 2 | Optimal | Optimal | 1.16e-08 | 5 | 5 | 23us | 982us | 43.0x | PASS |
| HS23 | 2 | 5 | Optimal | Optimal | 7.89e-01 | 103 | 9 | 1.3ms | 1.6ms | 1.2x | MISMATCH |
| HS24 | 2 | 3 | Optimal | Optimal | 1.80e-08 | 7 | 14 | 29us | 2.7ms | 94.3x | PASS |
| HS25 | 3 | 0 | Acceptable | Optimal | 6.08e-05 | 16 | 27 | 650.2ms | 5.3ms | 0.0x | PASS |
| HS25NE | 3 | 99 | LocalInfeasi | IpoptStatus( | N/A | 14 | 0 | 577us | 99us | 0.2x | BOTH_FAIL |
| HS26 | 3 | 1 | Acceptable | Optimal | 2.94e-12 | 18 | 25 | 53us | 2.9ms | 55.4x | PASS |
| HS268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 21us | 2.3ms | 114.3x | MISMATCH |
| HS27 | 3 | 1 | Optimal | Optimal | 6.69e-13 | 17 | 57 | 47us | 8.1ms | 171.5x | PASS |
| HS28 | 3 | 1 | Optimal | Optimal | 9.24e-31 | 1 | 1 | 11us | 310us | 28.3x | PASS |
| HS29 | 3 | 1 | Optimal | Optimal | 3.10e-10 | 18 | 7 | 55us | 1.3ms | 23.6x | PASS |
| HS2NE | 2 | 2 | MaxIteration | Infeasible | N/A | 2999 | 12 | 168.1ms | 2.5ms | 0.0x | BOTH_FAIL |
| HS3 | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 49us | 795us | 16.2x | PASS |
| HS30 | 3 | 1 | Acceptable | Optimal | 2.36e-06 | 19 | 7 | 146us | 1.3ms | 9.0x | PASS |
| HS31 | 3 | 1 | Optimal | Optimal | 8.67e-09 | 11 | 6 | 30us | 1.1ms | 38.2x | PASS |
| HS32 | 3 | 2 | Acceptable | Optimal | 1.21e-06 | 19 | 15 | 1.8ms | 2.5ms | 1.4x | PASS |
| HS33 | 3 | 2 | Acceptable | Optimal | 8.21e-08 | 14 | 9 | 143us | 2.1ms | 14.7x | PASS |
| HS34 | 3 | 2 | Optimal | Optimal | 9.46e-09 | 13 | 7 | 32us | 1.3ms | 40.7x | PASS |
| HS35 | 3 | 1 | Optimal | Optimal | 3.81e-09 | 9 | 7 | 27us | 1.3ms | 48.4x | PASS |
| HS35I | 3 | 1 | Optimal | Optimal | 3.44e-09 | 9 | 7 | 27us | 1.4ms | 54.0x | PASS |
| HS35MOD | 3 | 1 | Optimal | Optimal | 3.62e-08 | 9 | 14 | 32us | 2.4ms | 75.9x | PASS |
| HS36 | 3 | 1 | Optimal | Optimal | 6.33e-09 | 16 | 11 | 40us | 2.1ms | 53.0x | PASS |
| HS37 | 3 | 2 | Optimal | Optimal | 4.17e-10 | 10 | 11 | 39us | 2.2ms | 57.3x | PASS |
| HS38 | 4 | 0 | Optimal | Optimal | 1.23e-21 | 24 | 39 | 17us | 6.7ms | 400.6x | PASS |
| HS39 | 4 | 2 | Optimal | Optimal | 1.15e-11 | 15 | 13 | 41us | 1.7ms | 41.6x | PASS |
| HS3MOD | 2 | 0 | Optimal | Optimal | 1.00e-08 | 5 | 4 | 59us | 734us | 12.4x | PASS |
| HS4 | 2 | 0 | Optimal | Optimal | 1.97e-08 | 6 | 4 | 16.3ms | 752us | 0.0x | PASS |
| HS40 | 4 | 3 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 24us | 526us | 21.8x | PASS |
| HS41 | 4 | 1 | Optimal | Optimal | 1.15e-09 | 11 | 7 | 33us | 1.3ms | 40.0x | PASS |
| HS42 | 4 | 2 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 20us | 806us | 40.5x | PASS |
| HS43 | 4 | 3 | Optimal | Optimal | 6.81e-10 | 8 | 8 | 38us | 1.5ms | 39.0x | PASS |
| HS44 | 4 | 6 | Acceptable | Optimal | 2.97e-06 | 28 | 24 | 1.9ms | 4.5ms | 2.4x | PASS |
| HS44NEW | 4 | 6 | Acceptable | Optimal | 1.15e-05 | 19 | 18 | 1.9ms | 3.6ms | 2.0x | PASS |
| HS45 | 5 | 0 | Acceptable | Optimal | 3.14e-04 | 31 | 11 | 15.6ms | 2.1ms | 0.1x | MISMATCH |
| HS46 | 5 | 2 | Acceptable | Optimal | 3.43e-11 | 17 | 19 | 72us | 2.3ms | 31.5x | PASS |
| HS47 | 5 | 3 | Acceptable | Optimal | 3.22e-12 | 17 | 19 | 85us | 2.5ms | 28.9x | PASS |
| HS48 | 5 | 2 | Optimal | Optimal | 4.44e-31 | 1 | 1 | 14us | 336us | 24.8x | PASS |
| HS49 | 5 | 2 | Acceptable | Optimal | 2.61e-10 | 17 | 19 | 63us | 2.3ms | 37.1x | PASS |
| HS5 | 2 | 0 | Optimal | Optimal | 2.32e-16 | 7 | 7 | 5us | 1.2ms | 254.0x | PASS |
| HS50 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 9 | 9 | 37us | 1.3ms | 33.6x | PASS |
| HS51 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 16us | 359us | 21.8x | PASS |
| HS52 | 5 | 3 | Optimal | Optimal | 1.17e-15 | 1 | 1 | 17us | 314us | 18.1x | PASS |
| HS53 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 8 | 6 | 30us | 1.2ms | 39.0x | PASS |
| HS54 | 6 | 1 | Optimal | Optimal | 7.51e-01 | 10 | 15 | 41us | 2.7ms | 64.4x | MISMATCH |
| HS55 | 6 | 6 | Optimal | Optimal | 2.04e-02 | 11 | 18 | 54us | 4.3ms | 79.4x | MISMATCH |
| HS56 | 7 | 4 | Optimal | Optimal | 1.99e-14 | 5 | 10 | 35us | 1.6ms | 44.5x | PASS |
| HS57 | 2 | 1 | Optimal | Optimal | 6.11e-15 | 10 | 10 | 40us | 1.5ms | 37.4x | PASS |
| HS59 | 2 | 3 | Optimal | Optimal | 1.58e-14 | 19 | 17 | 5.3ms | 3.3ms | 0.6x | PASS |
| HS6 | 2 | 1 | Optimal | Optimal | 4.93e-32 | 7 | 5 | 24us | 879us | 36.3x | PASS |
| HS60 | 3 | 1 | Optimal | Optimal | 1.19e-13 | 8 | 6 | 22us | 1.2ms | 52.7x | PASS |
| HS61 | 3 | 2 | Optimal | Optimal | 7.91e-16 | 9 | 10 | 29us | 1.4ms | 46.1x | PASS |
| HS62 | 3 | 1 | Optimal | Optimal | 1.38e-16 | 9 | 6 | 30us | 1.3ms | 43.4x | PASS |
| HS63 | 3 | 2 | Optimal | Optimal | 0.00e+00 | 9 | 5 | 27us | 940us | 34.6x | PASS |
| HS64 | 3 | 1 | Optimal | Optimal | 3.62e-09 | 18 | 16 | 43us | 3.0ms | 68.3x | PASS |
| HS65 | 3 | 1 | Optimal | Optimal | 7.39e-09 | 15 | 16 | 43us | 3.3ms | 77.4x | PASS |
| HS66 | 3 | 2 | Optimal | Optimal | 1.15e-08 | 10 | 10 | 28us | 1.6ms | 57.7x | PASS |
| HS67 | 3 | 14 | Optimal | Optimal | 2.99e-09 | 12 | 9 | 113us | 1.8ms | 16.3x | PASS |
| HS68 | 4 | 2 | Optimal | Optimal | 5.21e-11 | 18 | 16 | 69us | 2.8ms | 40.7x | PASS |
| HS69 | 4 | 2 | Optimal | Optimal | 2.85e-15 | 11 | 10 | 38us | 2.0ms | 53.2x | PASS |
| HS7 | 2 | 1 | Optimal | Optimal | 5.30e-12 | 10 | 27 | 24us | 4.1ms | 170.8x | PASS |
| HS70 | 4 | 1 | Optimal | Optimal | 1.80e-01 | 10 | 46 | 127us | 8.0ms | 63.1x | MISMATCH |
| HS71 | 4 | 2 | Optimal | Optimal | 1.01e-09 | 10 | 8 | 41us | 1.5ms | 37.3x | PASS |
| HS72 | 4 | 2 | Optimal | Optimal | 6.74e-07 | 13 | 16 | 37us | 2.6ms | 71.8x | PASS |
| HS73 | 4 | 3 | Optimal | Optimal | 5.36e-10 | 12 | 8 | 36us | 1.5ms | 42.0x | PASS |
| HS74 | 4 | 5 | Optimal | Optimal | 0.00e+00 | 15 | 8 | 66us | 1.6ms | 23.9x | PASS |
| HS75 | 4 | 5 | Optimal | Optimal | 5.37e-09 | 15 | 8 | 67us | 1.5ms | 22.0x | PASS |
| HS76 | 4 | 3 | Optimal | Optimal | 5.89e-09 | 9 | 7 | 31us | 1.4ms | 44.3x | PASS |
| HS76I | 4 | 3 | Optimal | Optimal | 3.59e-09 | 10 | 6 | 33us | 1.2ms | 35.0x | PASS |
| HS77 | 5 | 2 | Optimal | Optimal | 1.99e-11 | 9 | 11 | 35us | 1.5ms | 42.8x | PASS |
| HS78 | 5 | 3 | Optimal | Optimal | 1.52e-16 | 4 | 4 | 25us | 632us | 24.9x | PASS |
| HS79 | 5 | 3 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 28us | 698us | 24.5x | PASS |
| HS8 | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 23us | 739us | 32.7x | PASS |
| HS80 | 5 | 3 | Optimal | Optimal | 6.75e-13 | 9 | 5 | 39us | 1.1ms | 27.3x | PASS |
| HS81 | 5 | 3 | Optimal | Optimal | 3.46e-14 | 31 | 68 | 118us | 11.8ms | 100.5x | PASS |
| HS83 | 5 | 3 | Acceptable | Optimal | 2.47e-06 | 8991 | 9 | 236.8ms | 1.6ms | 0.0x | PASS |
| HS84 | 5 | 3 | Acceptable | Optimal | 2.55e-08 | 3019 | 9 | 103.8ms | 1.7ms | 0.0x | PASS |
| HS85 | 5 | 21 | Optimal | Optimal | 4.42e-01 | 555 | 13 | 281.9ms | 3.8ms | 0.0x | MISMATCH |
| HS86 | 5 | 10 | Optimal | Optimal | 2.21e-09 | 11 | 10 | 64us | 2.0ms | 31.3x | PASS |
| HS87 | 6 | 4 | Acceptable | MaxIteration | N/A | 5994 | 3000 | 206.8ms | 489.9ms | 2.4x | ipopt_FAIL |
| HS88 | 2 | 1 | Optimal | Optimal | 6.34e-06 | 24 | 18 | 1.1ms | 3.6ms | 3.4x | PASS |
| HS89 | 3 | 1 | Optimal | Optimal | 3.42e-06 | 45 | 15 | 4.6ms | 3.5ms | 0.8x | PASS |
| HS9 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 3 | 3 | 17us | 1.0ms | 58.3x | PASS |
| HS90 | 4 | 1 | Optimal | Optimal | 5.20e-06 | 18 | 16 | 1.3ms | 3.8ms | 2.9x | PASS |
| HS91 | 5 | 1 | Optimal | Optimal | 6.44e-06 | 18 | 16 | 1.8ms | 4.2ms | 2.3x | PASS |
| HS92 | 6 | 1 | Optimal | Optimal | 6.80e-06 | 16 | 35 | 1.9ms | 10.0ms | 5.3x | PASS |
| HS93 | 6 | 2 | Optimal | Optimal | 9.89e-09 | 14 | 7 | 77us | 1.4ms | 18.6x | PASS |
| HS95 | 6 | 4 | Optimal | Optimal | 1.83e-06 | 58 | 9 | 4.8ms | 1.7ms | 0.4x | PASS |
| HS96 | 6 | 4 | Acceptable | Optimal | 1.20e-04 | 58 | 8 | 383us | 1.5ms | 4.0x | MISMATCH |
| HS97 | 6 | 4 | Acceptable | Optimal | 1.05e-06 | 48 | 24 | 380us | 4.4ms | 11.5x | PASS |
| HS98 | 6 | 4 | Acceptable | Optimal | 1.38e-06 | 169 | 13 | 702us | 2.4ms | 3.4x | PASS |
| HS99 | 7 | 2 | Optimal | Optimal | 0.00e+00 | 8 | 5 | 44us | 973us | 22.2x | PASS |
| HS99EXP | 31 | 21 | Optimal | Optimal | 0.00e+00 | 9 | 17 | 434us | 3.1ms | 7.0x | PASS |
| HUBFIT | 2 | 1 | Optimal | Optimal | 2.24e-09 | 7 | 7 | 20us | 1.3ms | 64.1x | PASS |
| HUMPS | 2 | 0 | Acceptable | Optimal | 2.78e-14 | 132 | 1533 | 97us | 216.3ms | 2240.2x | PASS |
| HYDC20LS | 99 | 0 | Acceptable | Optimal | 2.98e-01 | 2999 | 639 | 1.04s | 207.7ms | 0.2x | MISMATCH |
| HYDCAR20 | 99 | 99 | Optimal | Optimal | 0.00e+00 | 11 | 9 | 211.3ms | 2.5ms | 0.0x | PASS |
| HYDCAR6 | 29 | 29 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 10.2ms | 1.0ms | 0.1x | PASS |
| HYDCAR6LS | 29 | 0 | Acceptable | Optimal | 8.94e-02 | 931 | 149 | 30.2ms | 33.1ms | 1.1x | MISMATCH |
| HYPCIR | 2 | 2 | Optimal | Optimal | 0.00e+00 | 6 | 5 | 22us | 909us | 42.0x | PASS |
| JENSMP | 2 | 0 | Optimal | Optimal | 9.38e-01 | 1 | 9 | 5us | 1.4ms | 293.1x | MISMATCH |
| JENSMPNE | 2 | 10 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 85us | 122us | 1.4x | BOTH_FAIL |
| JUDGE | 2 | 0 | Optimal | Optimal | 2.15e-01 | 17 | 9 | 19us | 1.4ms | 69.6x | MISMATCH |
| JUDGEB | 2 | 0 | Optimal | Optimal | 2.15e-01 | 17 | 9 | 42us | 1.6ms | 38.0x | MISMATCH |
| JUDGENE | 2 | 20 | LocalInfeasi | IpoptStatus( | N/A | 15 | 0 | 87us | 126us | 1.4x | BOTH_FAIL |
| KIRBY2 | 5 | 151 | LocalInfeasi | IpoptStatus( | N/A | 118 | 0 | 5.1ms | 129us | 0.0x | BOTH_FAIL |
| KIRBY2LS | 5 | 0 | Acceptable | Optimal | 1.21e-14 | 26 | 11 | 2.0ms | 2.2ms | 1.1x | PASS |
| KIWCRESC | 3 | 2 | Optimal | Optimal | 1.07e-08 | 13 | 8 | 41us | 1.5ms | 36.9x | PASS |
| KOEBHELB | 3 | 0 | Optimal | Optimal | 1.83e-16 | 164 | 71 | 1.1ms | 14.9ms | 13.0x | PASS |
| KOEBHELBNE | 3 | 156 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 6.1ms | 129us | 0.0x | BOTH_FAIL |
| KOWOSB | 4 | 0 | Optimal | Optimal | 4.11e-17 | 21 | 8 | 19us | 1.6ms | 85.2x | PASS |
| KOWOSBNE | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 95us | 90us | 1.0x | BOTH_FAIL |
| KSIP | 20 | 1001 | Optimal | Optimal | 9.79e-01 | 1 | 22 | 20.2ms | 71.8ms | 3.5x | MISMATCH |
| LAKES | 90 | 78 | Optimal | Optimal | 1.53e-13 | 13 | 11 | 1.0ms | 3.1ms | 3.0x | PASS |
| LANCZOS1 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 41 | 0 | 444us | 168us | 0.4x | ipopt_FAIL |
| LANCZOS1LS | 6 | 0 | Acceptable | Optimal | 4.35e-06 | 55 | 115 | 95us | 19.6ms | 206.7x | PASS |
| LANCZOS2 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 70 | 0 | 687us | 131us | 0.2x | ipopt_FAIL |
| LANCZOS2LS | 6 | 0 | Acceptable | Optimal | 4.37e-06 | 54 | 101 | 88us | 16.6ms | 189.0x | PASS |
| LANCZOS3 | 6 | 24 | Acceptable | IpoptStatus( | N/A | 30 | 0 | 313us | 99us | 0.3x | ipopt_FAIL |
| LANCZOS3LS | 6 | 0 | Optimal | Optimal | 4.33e-06 | 67 | 174 | 125us | 31.9ms | 254.6x | PASS |
| LAUNCH | 25 | 28 | Acceptable | Optimal | 3.93e-07 | 675 | 12 | 79.1ms | 2.6ms | 0.0x | PASS |
| LEVYMONE10 | 10 | 20 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 105us | 96us | 0.9x | BOTH_FAIL |
| LEVYMONE5 | 2 | 4 | LocalInfeasi | IpoptStatus( | N/A | 10 | 0 | 48us | 101us | 2.1x | BOTH_FAIL |
| LEVYMONE6 | 3 | 6 | Optimal | IpoptStatus( | N/A | 29 | 0 | 65us | 88us | 1.3x | ipopt_FAIL |
| LEVYMONE7 | 4 | 8 | LocalInfeasi | IpoptStatus( | N/A | 13 | 0 | 86us | 124us | 1.4x | BOTH_FAIL |
| LEVYMONE8 | 5 | 10 | LocalInfeasi | IpoptStatus( | N/A | 8 | 0 | 62us | 131us | 2.1x | BOTH_FAIL |
| LEVYMONE9 | 8 | 16 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 95us | 100us | 1.1x | BOTH_FAIL |
| LEVYMONT10 | 10 | 0 | Optimal | Optimal | 1.74e-16 | 14 | 4 | 19us | 866us | 46.1x | PASS |
| LEVYMONT5 | 2 | 0 | Acceptable | Optimal | 1.00e+00 | 10 | 10 | 19us | 1.9ms | 97.5x | MISMATCH |
| LEVYMONT6 | 3 | 0 | Optimal | Optimal | 1.00e+00 | 27 | 8 | 18us | 1.5ms | 84.8x | MISMATCH |
| LEVYMONT7 | 4 | 0 | Optimal | Optimal | 1.88e-01 | 14 | 7 | 13us | 1.6ms | 117.9x | MISMATCH |
| LEVYMONT8 | 5 | 0 | Optimal | Optimal | 9.71e-01 | 12 | 4 | 13us | 762us | 57.3x | MISMATCH |
| LEVYMONT9 | 8 | 0 | Optimal | Optimal | 6.77e-01 | 12 | 4 | 16us | 776us | 49.5x | MISMATCH |
| LEWISPOL | 6 | 9 | Acceptable | IpoptStatus( | N/A | 10 | 0 | 4.5ms | 92us | 0.0x | ipopt_FAIL |
| LHAIFAM | 99 | 150 | Optimal | InvalidNumbe | N/A | 0 | 0 | 979.6ms | 299us | 0.0x | ipopt_FAIL |
| LIN | 4 | 2 | Optimal | Optimal | 2.03e-03 | 9 | 7 | 39us | 1.2ms | 30.8x | MISMATCH |
| LINSPANH | 97 | 33 | Acceptable | Optimal | 5.18e-07 | 2999 | 24 | 271.3ms | 5.1ms | 0.0x | PASS |
| LOADBAL | 31 | 31 | Optimal | Optimal | 6.55e-09 | 15 | 13 | 839us | 2.8ms | 3.4x | PASS |
| LOGHAIRY | 2 | 0 | Optimal | Optimal | 0.00e+00 | 408 | 2747 | 291us | 415.6ms | 1427.3x | PASS |
| LOGROS | 2 | 0 | Optimal | Optimal | 0.00e+00 | 50 | 49 | 119us | 9.7ms | 81.3x | PASS |
| LOOTSMA | 3 | 2 | Optimal | Optimal | 8.76e-08 | 26 | 13 | 71us | 2.4ms | 34.2x | PASS |
| LOTSCHD | 12 | 7 | Optimal | Optimal | 4.44e-10 | 14 | 9 | 106us | 1.6ms | 14.9x | PASS |
| LRCOVTYPE | 54 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| LRIJCNN1 | 22 | 0 | Acceptable | Optimal | 5.92e-03 | 31 | 11 | 87.2ms | 183.3ms | 2.1x | MISMATCH |
| LSC1 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 55 | 0 | 117us | 101us | 0.9x | BOTH_FAIL |
| LSC1LS | 3 | 0 | Optimal | Optimal | 1.27e-15 | 49 | 16 | 32us | 2.7ms | 84.7x | PASS |
| LSC2 | 3 | 6 | LocalInfeasi | IpoptStatus( | N/A | 93 | 0 | 314us | 100us | 0.3x | BOTH_FAIL |
| LSC2LS | 3 | 0 | Acceptable | Optimal | 2.32e-04 | 32 | 38 | 28.8ms | 4.8ms | 0.2x | MISMATCH |
| LSNNODOC | 5 | 4 | Acceptable | Optimal | 3.42e-06 | 13 | 10 | 2.3ms | 1.9ms | 0.8x | PASS |
| LSQFIT | 2 | 1 | Optimal | Optimal | 4.24e-09 | 7 | 7 | 22us | 1.3ms | 59.5x | PASS |
| MADSEN | 3 | 6 | Optimal | Optimal | 8.02e-09 | 29 | 18 | 102us | 3.3ms | 32.3x | PASS |
| MAKELA1 | 3 | 2 | Optimal | Optimal | 2.00e+00 | 6 | 12 | 28us | 2.5ms | 89.9x | MISMATCH |
| MAKELA2 | 3 | 3 | Optimal | Optimal | 1.27e-01 | 3 | 6 | 18us | 1.2ms | 66.4x | MISMATCH |
| MAKELA3 | 21 | 20 | Optimal | Optimal | 1.00e-08 | 40 | 11 | 1.1ms | 2.4ms | 2.3x | PASS |
| MAKELA4 | 21 | 40 | Optimal | Optimal | 6.88e-01 | 1 | 5 | 79us | 1.2ms | 15.3x | MISMATCH |
| MARATOS | 2 | 1 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 14us | 628us | 46.2x | PASS |
| MARATOSB | 2 | 0 | Optimal | Optimal | 0.00e+00 | 1107 | 672 | 500us | 103.9ms | 207.6x | PASS |
| MATRIX2 | 6 | 2 | Acceptable | Optimal | 7.84e-10 | 12 | 42 | 3.4ms | 6.8ms | 2.0x | PASS |
| MAXLIKA | 8 | 0 | Acceptable | Optimal | 1.13e-02 | 20 | 23 | 1.99s | 6.8ms | 0.0x | MISMATCH |
| MCONCON | 15 | 11 | Acceptable | Optimal | 6.32e-08 | 31 | 7 | 7.9ms | 1.4ms | 0.2x | PASS |
| MDHOLE | 2 | 0 | Optimal | Optimal | 9.98e-09 | 35 | 42 | 185us | 8.2ms | 44.3x | PASS |
| MESH | 41 | 48 | Optimal | IpoptStatus( | N/A | 112 | 79 | 29.0ms | 20.2ms | 0.7x | ipopt_FAIL |
| METHANB8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 9 | 3 | 322us | 676us | 2.1x | PASS |
| METHANB8LS | 31 | 0 | Optimal | Optimal | 5.35e-26 | 9 | 8 | 8.4ms | 1.3ms | 0.2x | PASS |
| METHANL8 | 31 | 31 | Optimal | Optimal | 0.00e+00 | 4 | 4 | 10.0ms | 831us | 0.1x | PASS |
| METHANL8LS | 31 | 0 | Optimal | Optimal | 5.70e-17 | 977 | 40 | 30.5ms | 8.9ms | 0.3x | PASS |
| MEXHAT | 2 | 0 | Optimal | Optimal | 6.77e-10 | 45 | 26 | 26us | 3.5ms | 135.5x | PASS |
| MEYER3 | 3 | 0 | Acceptable | Optimal | 3.05e-12 | 720 | 194 | 56.1ms | 30.3ms | 0.5x | PASS |
| MEYER3NE | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 731 | 0 | 135.8ms | 113us | 0.0x | BOTH_FAIL |
| MGH09 | 4 | 11 | LocalInfeasi | IpoptStatus( | N/A | 59 | 0 | 442us | 98us | 0.2x | BOTH_FAIL |
| MGH09LS | 4 | 0 | Acceptable | Optimal | 7.12e-04 | 47 | 72 | 35us | 12.2ms | 351.0x | MISMATCH |
| MGH10 | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 107.9ms | 100us | 0.0x | BOTH_FAIL |
| MGH10LS | 3 | 0 | Acceptable | Optimal | 4.85e-12 | 2999 | 1828 | 49.2ms | 284.5ms | 5.8x | PASS |
| MGH10S | 3 | 16 | LocalInfeasi | IpoptStatus( | N/A | 7 | 0 | 2.05s | 107us | 0.0x | BOTH_FAIL |
| MGH10SLS | 3 | 0 | Acceptable | Optimal | 1.00e+00 | 1304 | 354 | 1.4ms | 55.6ms | 39.5x | MISMATCH |
| MGH17 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 933us | 101us | 0.1x | BOTH_FAIL |
| MGH17LS | 5 | 0 | Acceptable | Optimal | 2.43e-05 | 214 | 47 | 543us | 8.8ms | 16.3x | PASS |
| MGH17S | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 67 | 0 | 948us | 103us | 0.1x | BOTH_FAIL |
| MGH17SLS | 5 | 0 | Optimal | Optimal | 4.18e-07 | 71 | 41 | 201us | 7.7ms | 38.4x | PASS |
| MIFFLIN1 | 3 | 2 | Optimal | Optimal | 9.24e-09 | 6 | 5 | 26us | 1.0ms | 38.8x | PASS |
| MIFFLIN2 | 3 | 2 | Optimal | Optimal | 9.97e-09 | 15 | 11 | 38us | 2.1ms | 54.1x | PASS |
| MINMAXBD | 5 | 20 | Optimal | Optimal | 5.55e-11 | 469 | 25 | 4.3ms | 5.9ms | 1.4x | PASS |
| MINMAXRB | 3 | 4 | Optimal | Optimal | 9.85e-09 | 3 | 8 | 22us | 1.5ms | 67.1x | PASS |
| MINSURF | 64 | 0 | Optimal | Optimal | 0.00e+00 | 17 | 4 | 49us | 1.1ms | 22.5x | PASS |
| MISRA1A | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 30 | 0 | 495us | 114us | 0.2x | BOTH_FAIL |
| MISRA1ALS | 2 | 0 | Acceptable | Optimal | 1.91e-14 | 35 | 40 | 168us | 6.2ms | 36.7x | PASS |
| MISRA1B | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 273us | 103us | 0.4x | BOTH_FAIL |
| MISRA1BLS | 2 | 0 | Optimal | Optimal | 4.76e-14 | 25 | 34 | 117us | 4.9ms | 41.6x | PASS |
| MISRA1C | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 18 | 0 | 657us | 96us | 0.1x | BOTH_FAIL |
| MISRA1CLS | 2 | 0 | Acceptable | Optimal | 3.55e-14 | 59 | 14 | 161us | 2.3ms | 14.0x | PASS |
| MISRA1D | 2 | 14 | LocalInfeasi | IpoptStatus( | N/A | 24 | 0 | 108.8ms | 138us | 0.0x | BOTH_FAIL |
| MISRA1DLS | 2 | 0 | Acceptable | Optimal | 4.95e-15 | 56 | 30 | 157us | 4.6ms | 29.0x | PASS |
| MISTAKE | 9 | 13 | Optimal | Optimal | 3.43e-01 | 685 | 16 | 23.5ms | 3.4ms | 0.1x | MISMATCH |
| MRIBASIS | 36 | 55 | Acceptable | Optimal | 1.16e-08 | 684 | 15 | 106.3ms | 3.9ms | 0.0x | PASS |
| MSS1 | 90 | 73 | Acceptable | Optimal | 1.25e-01 | 74 | 95 | 4.61s | 49.8ms | 0.0x | MISMATCH |
| MUONSINE | 1 | 512 | LocalInfeasi | IpoptStatus( | N/A | 9 | 0 | 4.72s | 144us | 0.0x | BOTH_FAIL |
| MUONSINELS | 1 | 0 | Acceptable | Optimal | 1.42e-01 | 13 | 8 | 6.4ms | 1.5ms | 0.2x | MISMATCH |
| MWRIGHT | 5 | 3 | Optimal | Optimal | 9.48e-01 | 13 | 10 | 50us | 1.6ms | 31.1x | MISMATCH |
| NASH | 72 | 24 | RestorationF | Infeasible | N/A | 50 | 45 | 640.5ms | 12.2ms | 0.0x | BOTH_FAIL |
| NELSON | 3 | 128 | LocalInfeasi | IpoptStatus( | N/A | 573 | 0 | 102.4ms | 110us | 0.0x | BOTH_FAIL |
| NET1 | 48 | 57 | Acceptable | Optimal | 7.72e-08 | 53 | 26 | 1.28s | 5.8ms | 0.0x | PASS |
| NYSTROM5 | 18 | 20 | Optimal | IpoptStatus( | N/A | 21 | 0 | 431us | 110us | 0.3x | ipopt_FAIL |
| NYSTROM5C | 18 | 20 | Optimal | IpoptStatus( | N/A | 21 | 0 | 358us | 104us | 0.3x | ipopt_FAIL |
| ODFITS | 10 | 6 | Optimal | Optimal | 1.91e-16 | 11 | 8 | 59us | 1.5ms | 26.1x | PASS |
| OET1 | 3 | 1002 | Optimal | Optimal | 4.53e-02 | 53 | 33 | 22.0ms | 49.8ms | 2.3x | MISMATCH |
| OET2 | 3 | 1002 | Acceptable | Optimal | 4.46e-02 | 605 | 181 | 992.7ms | 281.5ms | 0.3x | MISMATCH |
| OET3 | 4 | 1002 | Optimal | Optimal | 5.47e-04 | 94 | 13 | 73.6ms | 20.9ms | 0.3x | MISMATCH |
| OET4 | 4 | 1002 | Optimal | Optimal | 7.66e-09 | 288 | 165 | 101.2ms | 258.7ms | 2.6x | PASS |
| OET5 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET6 | 5 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OET7 | 7 | 1002 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| OPTCNTRL | 32 | 20 | Optimal | Optimal | 2.03e-09 | 36 | 9 | 1.5ms | 1.8ms | 1.2x | PASS |
| OPTPRLOC | 30 | 30 | Optimal | Optimal | 1.22e-08 | 57 | 13 | 2.8ms | 3.1ms | 1.1x | PASS |
| ORTHREGB | 27 | 6 | Optimal | Optimal | 4.20e-19 | 2 | 2 | 68us | 714us | 10.6x | PASS |
| OSBORNE1 | 5 | 33 | LocalInfeasi | IpoptStatus( | N/A | 102 | 0 | 31.8ms | 112us | 0.0x | BOTH_FAIL |
| OSBORNE2 | 11 | 65 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 2.0ms | 121us | 0.1x | BOTH_FAIL |
| OSBORNEA | 5 | 0 | Acceptable | Optimal | 2.05e-18 | 99 | 64 | 307us | 10.5ms | 34.2x | PASS |
| OSBORNEB | 11 | 0 | Acceptable | Optimal | 2.34e-12 | 95 | 19 | 502us | 3.2ms | 6.4x | PASS |
| OSLBQP | 8 | 0 | Acceptable | Optimal | 7.24e-07 | 13 | 15 | 20.5ms | 2.4ms | 0.1x | PASS |
| PALMER1 | 4 | 0 | Acceptable | Optimal | 1.55e-16 | 34 | 13 | 88us | 2.3ms | 26.0x | PASS |
| PALMER1A | 6 | 0 | Acceptable | Optimal | 1.33e-14 | 173 | 48 | 629us | 9.7ms | 15.4x | PASS |
| PALMER1ANE | 6 | 35 | LocalInfeasi | IpoptStatus( | N/A | 187 | 0 | 3.9ms | 164us | 0.0x | BOTH_FAIL |
| PALMER1B | 4 | 0 | Acceptable | Optimal | 4.87e-14 | 51 | 17 | 332us | 3.0ms | 8.9x | PASS |
| PALMER1BNE | 4 | 35 | LocalInfeasi | IpoptStatus( | N/A | 17 | 0 | 2.6ms | 100us | 0.0x | BOTH_FAIL |
| PALMER1C | 8 | 0 | Optimal | Optimal | 1.79e-13 | 4 | 1 | 61.0ms | 281us | 0.0x | PASS |
| PALMER1D | 7 | 0 | Acceptable | Optimal | 1.65e-13 | 196 | 1 | 458us | 395us | 0.9x | PASS |
| PALMER1E | 8 | 0 | Acceptable | Optimal | 4.05e-13 | 675 | 55 | 1.8ms | 10.8ms | 6.2x | PASS |
| PALMER1ENE | 8 | 35 | LocalInfeasi | IpoptStatus( | N/A | 121 | 0 | 11.3ms | 104us | 0.0x | BOTH_FAIL |
| PALMER1NE | 4 | 31 | LocalInfeasi | IpoptStatus( | N/A | 25 | 0 | 840us | 114us | 0.1x | BOTH_FAIL |
| PALMER2 | 4 | 0 | Optimal | Optimal | 2.49e-16 | 24 | 28 | 41us | 6.1ms | 148.2x | PASS |
| PALMER2A | 6 | 0 | Optimal | Optimal | 1.78e-15 | 142 | 91 | 332us | 20.1ms | 60.3x | PASS |
| PALMER2ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 117 | 0 | 2.1ms | 101us | 0.0x | BOTH_FAIL |
| PALMER2B | 4 | 0 | Acceptable | Optimal | 3.44e-15 | 37 | 15 | 76us | 3.2ms | 42.0x | PASS |
| PALMER2BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 22 | 0 | 1.2ms | 100us | 0.1x | BOTH_FAIL |
| PALMER2C | 8 | 0 | Acceptable | Optimal | 1.30e-14 | 965 | 1 | 1.2ms | 241us | 0.2x | PASS |
| PALMER2E | 8 | 0 | Acceptable | Optimal | 2.09e-12 | 633 | 114 | 1.0ms | 25.8ms | 25.8x | PASS |
| PALMER2ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 109 | 0 | 6.3ms | 163us | 0.0x | BOTH_FAIL |
| PALMER2NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 29 | 0 | 837us | 121us | 0.1x | BOTH_FAIL |
| PALMER3 | 4 | 0 | Acceptable | Optimal | 6.25e-02 | 20 | 44 | 75.9ms | 7.5ms | 0.1x | MISMATCH |
| PALMER3A | 6 | 0 | Acceptable | Optimal | 1.72e-15 | 139 | 73 | 250us | 15.0ms | 60.2x | PASS |
| PALMER3ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 74 | 0 | 2.5ms | 95us | 0.0x | BOTH_FAIL |
| PALMER3B | 4 | 0 | Optimal | Optimal | 8.40e-16 | 31 | 15 | 47us | 3.1ms | 65.3x | PASS |
| PALMER3BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 32 | 0 | 568us | 98us | 0.2x | BOTH_FAIL |
| PALMER3C | 8 | 0 | Acceptable | Optimal | 3.99e-12 | 680 | 1 | 893us | 257us | 0.3x | PASS |
| PALMER3E | 8 | 0 | Acceptable | Optimal | 1.13e-13 | 304 | 32 | 517us | 5.8ms | 11.2x | PASS |
| PALMER3ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 340 | 0 | 3.5ms | 104us | 0.0x | BOTH_FAIL |
| PALMER3NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 218.0ms | 94us | 0.0x | BOTH_FAIL |
| PALMER4 | 4 | 0 | Optimal | Optimal | 5.72e-02 | 29 | 16 | 76.6ms | 3.6ms | 0.0x | MISMATCH |
| PALMER4A | 6 | 0 | Acceptable | Optimal | 3.23e-15 | 115 | 53 | 173us | 10.1ms | 58.6x | PASS |
| PALMER4ANE | 6 | 23 | LocalInfeasi | IpoptStatus( | N/A | 57 | 0 | 2.1ms | 100us | 0.0x | BOTH_FAIL |
| PALMER4B | 4 | 0 | Acceptable | Optimal | 1.43e-15 | 31 | 16 | 74us | 3.4ms | 45.7x | PASS |
| PALMER4BNE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 606us | 113us | 0.2x | BOTH_FAIL |
| PALMER4C | 8 | 0 | Acceptable | Optimal | 3.12e-11 | 607 | 1 | 857us | 254us | 0.3x | PASS |
| PALMER4E | 8 | 0 | Acceptable | Optimal | 3.97e-13 | 290 | 25 | 463us | 4.8ms | 10.4x | PASS |
| PALMER4ENE | 8 | 23 | LocalInfeasi | IpoptStatus( | N/A | 513 | 0 | 4.1ms | 113us | 0.0x | BOTH_FAIL |
| PALMER4NE | 4 | 23 | LocalInfeasi | IpoptStatus( | N/A | 23 | 0 | 452.2ms | 102us | 0.0x | BOTH_FAIL |
| PALMER5A | 8 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 16.0ms | 690.7ms | 43.2x | ipopt_FAIL |
| PALMER5ANE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 11.1ms | 133us | 0.0x | BOTH_FAIL |
| PALMER5B | 9 | 0 | Acceptable | Optimal | 5.49e-03 | 1132 | 113 | 1.2ms | 21.6ms | 18.3x | MISMATCH |
| PALMER5BNE | 9 | 12 | LocalInfeasi | IpoptStatus( | N/A | 65 | 0 | 2.9ms | 102us | 0.0x | BOTH_FAIL |
| PALMER5C | 6 | 0 | Optimal | Optimal | 8.35e-16 | 13 | 1 | 11us | 258us | 23.4x | PASS |
| PALMER5D | 4 | 0 | Acceptable | Optimal | 4.23e-15 | 21 | 1 | 57us | 265us | 4.7x | PASS |
| PALMER5E | 8 | 0 | Acceptable | MaxIteration | N/A | 14 | 3000 | 3.2ms | 502.2ms | 157.2x | ipopt_FAIL |
| PALMER5ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 7.1ms | 99us | 0.0x | BOTH_FAIL |
| PALMER6A | 6 | 0 | Acceptable | Optimal | 3.48e-15 | 229 | 105 | 317us | 19.2ms | 60.6x | PASS |
| PALMER6ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 104 | 0 | 1.3ms | 99us | 0.1x | BOTH_FAIL |
| PALMER6C | 8 | 0 | Acceptable | Optimal | 3.49e-12 | 499 | 1 | 402us | 242us | 0.6x | PASS |
| PALMER6E | 8 | 0 | Acceptable | Optimal | 2.53e-11 | 255 | 30 | 305us | 6.0ms | 19.8x | PASS |
| PALMER6ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 399 | 0 | 1.2ms | 100us | 0.1x | BOTH_FAIL |
| PALMER7A | 6 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 13.0ms | 533.4ms | 41.0x | ipopt_FAIL |
| PALMER7ANE | 6 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 20.5ms | 104us | 0.0x | BOTH_FAIL |
| PALMER7C | 8 | 0 | Optimal | Optimal | 2.59e-13 | 2 | 1 | 35.5ms | 345us | 0.0x | PASS |
| PALMER7E | 8 | 0 | Acceptable | MaxIteration | N/A | 1550 | 3000 | 10.4ms | 658.0ms | 63.2x | ipopt_FAIL |
| PALMER7ENE | 8 | 13 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 26.0ms | 227us | 0.0x | BOTH_FAIL |
| PALMER8A | 6 | 0 | Acceptable | Optimal | 2.22e-16 | 84 | 36 | 125us | 7.8ms | 62.4x | PASS |
| PALMER8ANE | 6 | 12 | LocalInfeasi | IpoptStatus( | N/A | 36 | 0 | 472us | 98us | 0.2x | BOTH_FAIL |
| PALMER8C | 8 | 0 | Acceptable | Optimal | 1.03e-14 | 525 | 1 | 458us | 266us | 0.6x | PASS |
| PALMER8E | 8 | 0 | Acceptable | Optimal | 9.31e-14 | 257 | 23 | 282us | 4.3ms | 15.2x | PASS |
| PALMER8ENE | 8 | 12 | LocalInfeasi | IpoptStatus( | N/A | 28 | 0 | 1.3ms | 102us | 0.1x | BOTH_FAIL |
| PARKCH | 15 | 0 | Acceptable | Optimal | 1.27e-14 | 19 | 17 | 11.41s | 4.03s | 0.4x | PASS |
| PENTAGON | 6 | 15 | Optimal | Optimal | 1.47e-05 | 13 | 19 | 77us | 4.0ms | 52.5x | PASS |
| PFIT1 | 3 | 3 | Optimal | Infeasible | N/A | 257 | 266 | 692us | 47.6ms | 68.8x | ipopt_FAIL |
| PFIT1LS | 3 | 0 | Optimal | Optimal | 1.37e-20 | 264 | 263 | 173us | 46.1ms | 266.6x | PASS |
| PFIT2 | 3 | 3 | Optimal | RestorationF | N/A | 111 | 247 | 293us | 48.3ms | 165.1x | ipopt_FAIL |
| PFIT2LS | 3 | 0 | Optimal | Optimal | 1.90e-20 | 529 | 82 | 351us | 16.0ms | 45.6x | PASS |
| PFIT3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 115 | 133 | 300us | 26.3ms | 87.9x | PASS |
| PFIT3LS | 3 | 0 | Optimal | Optimal | 2.58e-20 | 772 | 132 | 519us | 25.0ms | 48.1x | PASS |
| PFIT4 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 209 | 190 | 559us | 35.9ms | 64.2x | PASS |
| PFIT4LS | 3 | 0 | Optimal | Optimal | 6.66e-20 | 950 | 215 | 633us | 44.3ms | 70.0x | PASS |
| POLAK1 | 3 | 2 | Optimal | Optimal | 3.67e-09 | 9 | 5 | 30us | 1.2ms | 38.7x | PASS |
| POLAK2 | 11 | 2 | Optimal | Optimal | 3.59e-11 | 28 | 10 | 135us | 1.8ms | 13.4x | PASS |
| POLAK3 | 12 | 10 | Optimal | MaxIteration | N/A | 15 | 3000 | 336us | 695.7ms | 2073.3x | ipopt_FAIL |
| POLAK4 | 3 | 3 | Acceptable | Optimal | 4.53e-09 | 14 | 4 | 92us | 942us | 10.2x | PASS |
| POLAK5 | 3 | 2 | Acceptable | Optimal | 4.38e-11 | 84 | 31 | 315us | 5.1ms | 16.3x | PASS |
| POLAK6 | 5 | 4 | Optimal | MaxIteration | N/A | 13 | 3000 | 58us | 900.9ms | 15588.3x | ipopt_FAIL |
| PORTFL1 | 12 | 1 | Acceptable | Optimal | 1.28e-06 | 10 | 9 | 19.8ms | 2.0ms | 0.1x | PASS |
| PORTFL2 | 12 | 1 | Acceptable | Optimal | 4.37e-07 | 10 | 8 | 19.9ms | 1.7ms | 0.1x | PASS |
| PORTFL3 | 12 | 1 | Acceptable | Optimal | 1.50e-06 | 10 | 9 | 19.8ms | 1.9ms | 0.1x | PASS |
| PORTFL4 | 12 | 1 | Acceptable | Optimal | 3.83e-06 | 9 | 8 | 20.0ms | 1.6ms | 0.1x | PASS |
| PORTFL6 | 12 | 1 | Acceptable | Optimal | 4.96e-07 | 10 | 8 | 19.9ms | 1.7ms | 0.1x | PASS |
| POWELLBS | 2 | 2 | Optimal | Optimal | 0.00e+00 | 89 | 11 | 176us | 1.4ms | 7.9x | PASS |
| POWELLBSLS | 2 | 0 | Optimal | Optimal | 6.41e-26 | 159 | 91 | 82us | 12.9ms | 157.0x | PASS |
| POWELLSQ | 2 | 2 | Acceptable | Infeasible | N/A | 9 | 29 | 27us | 5.1ms | 190.4x | ipopt_FAIL |
| POWELLSQLS | 2 | 0 | Optimal | Optimal | 1.53e-14 | 13 | 10 | 9us | 1.6ms | 171.3x | PASS |
| PRICE3NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 11 | 7 | 27us | 1.1ms | 39.6x | PASS |
| PRICE4 | 2 | 0 | Acceptable | Optimal | 1.85e-11 | 56 | 8 | 30us | 1.1ms | 36.0x | PASS |
| PRICE4B | 2 | 0 | Acceptable | Optimal | 2.19e-11 | 59 | 8 | 31us | 1.5ms | 47.8x | PASS |
| PRICE4NE | 2 | 2 | Optimal | Acceptable | 0.00e+00 | 10 | 23 | 26us | 3.2ms | 121.4x | PASS |
| PRODPL0 | 60 | 29 | Acceptable | Optimal | 1.33e-06 | 20 | 15 | 723.6ms | 3.3ms | 0.0x | PASS |
| PRODPL1 | 60 | 29 | Acceptable | Optimal | 9.46e-05 | 231 | 28 | 32.0ms | 6.5ms | 0.2x | PASS |
| PSPDOC | 4 | 0 | Optimal | Optimal | 3.50e-09 | 7 | 5 | 2.8ms | 1.1ms | 0.4x | PASS |
| PT | 2 | 501 | Optimal | Optimal | 8.94e-03 | 4 | 106 | 2.0ms | 83.2ms | 41.7x | MISMATCH |
| QC | 9 | 4 | Acceptable | Optimal | 9.35e-02 | 3998 | 44 | 24.7ms | 8.3ms | 0.3x | MISMATCH |
| QCNEW | 9 | 3 | Acceptable | Optimal | 8.91e-03 | 3998 | 6 | 21.1ms | 1.2ms | 0.1x | MISMATCH |
| QPCBLEND | 83 | 74 | Optimal | Optimal | 5.56e-07 | 37 | 19 | 3.0ms | 5.0ms | 1.7x | PASS |
| QPNBLEND | 83 | 74 | Optimal | Optimal | 4.77e-07 | 50 | 18 | 4.7ms | 5.0ms | 1.1x | PASS |
| RAT42 | 3 | 9 | LocalInfeasi | IpoptStatus( | N/A | 21 | 0 | 113us | 95us | 0.8x | BOTH_FAIL |
| RAT42LS | 3 | 0 | Optimal | Optimal | 8.82e-16 | 22 | 28 | 79us | 4.1ms | 52.1x | PASS |
| RAT43 | 4 | 15 | LocalInfeasi | IpoptStatus( | N/A | 3 | 0 | 139us | 115us | 0.8x | BOTH_FAIL |
| RAT43LS | 4 | 0 | Optimal | Optimal | 9.92e-01 | 3 | 34 | 21us | 5.5ms | 256.5x | MISMATCH |
| RECIPE | 3 | 3 | Acceptable | Optimal | 0.00e+00 | 19 | 16 | 46us | 2.4ms | 52.0x | PASS |
| RECIPELS | 3 | 0 | Acceptable | Optimal | 1.17e-11 | 32 | 29 | 19us | 4.6ms | 243.9x | PASS |
| RES | 20 | 14 | Acceptable | Optimal | 0.00e+00 | 12 | 10 | 8.6ms | 1.8ms | 0.2x | PASS |
| RK23 | 17 | 11 | Acceptable | Optimal | 6.20e-05 | 17 | 10 | 4.6ms | 2.3ms | 0.5x | PASS |
| ROBOT | 14 | 2 | Optimal | IpoptStatus( | N/A | 13 | 18 | 130us | 3.8ms | 29.0x | ipopt_FAIL |
| ROSENBR | 2 | 0 | Optimal | Optimal | 1.65e-21 | 34 | 21 | 16us | 3.0ms | 185.7x | PASS |
| ROSENBRTU | 2 | 0 | Optimal | Optimal | 4.56e-26 | 90 | 87 | 43us | 14.2ms | 329.6x | PASS |
| ROSENMMX | 5 | 4 | Optimal | Optimal | 2.27e-10 | 10 | 13 | 58us | 2.7ms | 47.3x | PASS |
| ROSZMAN1 | 4 | 25 | LocalInfeasi | IpoptStatus( | N/A | 11 | 0 | 611us | 124us | 0.2x | BOTH_FAIL |
| ROSZMAN1LS | 4 | 0 | Acceptable | Optimal | 3.90e-02 | 118 | 28 | 183us | 4.5ms | 24.7x | MISMATCH |
| RSNBRNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 41us | 309us | 7.6x | PASS |
| S268 | 5 | 5 | Optimal | Optimal | 1.00e+00 | 1 | 14 | 21us | 2.4ms | 116.5x | MISMATCH |
| S308 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 13 | 9 | 9us | 1.4ms | 154.3x | PASS |
| S308NE | 2 | 3 | LocalInfeasi | IpoptStatus( | N/A | 12 | 0 | 39us | 129us | 3.3x | BOTH_FAIL |
| S316-322 | 2 | 1 | Optimal | Optimal | 0.00e+00 | 9 | 7 | 24us | 1.1ms | 46.9x | PASS |
| S365 | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 43.7ms | 642us | 0.0x | BOTH_FAIL |
| S365MOD | 7 | 5 | MaxIteration | RestorationF | N/A | 2999 | 1 | 44.5ms | 632us | 0.0x | BOTH_FAIL |
| SANTA | 21 | 23 | LocalInfeasi | IpoptStatus( | N/A | 34 | 0 | 1.2ms | 96us | 0.1x | BOTH_FAIL |
| SANTALS | 21 | 0 | Acceptable | Optimal | 2.69e-09 | 141 | 31 | 397us | 7.0ms | 17.7x | PASS |
| SIM2BQP | 2 | 0 | Optimal | Optimal | 7.86e-09 | 6 | 5 | 27us | 890us | 32.8x | PASS |
| SIMBQP | 2 | 0 | Optimal | Optimal | 8.29e-09 | 6 | 5 | 216us | 909us | 4.2x | PASS |
| SIMPLLPA | 2 | 2 | Optimal | Optimal | 1.18e-08 | 6 | 8 | 23us | 1.5ms | 65.8x | PASS |
| SIMPLLPB | 2 | 3 | Optimal | Optimal | 9.07e-09 | 3 | 7 | 22us | 1.4ms | 62.4x | PASS |
| SINEVAL | 2 | 0 | Optimal | Optimal | 1.03e-19 | 75 | 42 | 39us | 6.6ms | 170.3x | PASS |
| SINVALNE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 21 | 1 | 44us | 299us | 6.8x | PASS |
| SIPOW1 | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 81 | 38.8ms | 228.0ms | 5.9x | MISMATCH |
| SIPOW1M | 2 | 2000 | Optimal | Optimal | 1.43e+00 | 1 | 88 | 40.5ms | 244.1ms | 6.0x | MISMATCH |
| SIPOW2 | 2 | 2000 | Optimal | Optimal | 7.73e-03 | 11 | 69 | 10.8ms | 181.3ms | 16.7x | MISMATCH |
| SIPOW2M | 2 | 2000 | Optimal | Optimal | 7.62e-03 | 11 | 73 | 11.7ms | 188.6ms | 16.1x | MISMATCH |
| SIPOW3 | 4 | 2000 | Optimal | Optimal | 9.99e-01 | 273 | 12 | 274.1ms | 36.8ms | 0.1x | MISMATCH |
| SIPOW4 | 4 | 2000 | Optimal | Optimal | 8.69e-01 | 1 | 11 | 56.6ms | 34.2ms | 0.6x | MISMATCH |
| SISSER | 2 | 0 | Acceptable | Optimal | 1.22e-10 | 15 | 18 | 10us | 2.1ms | 205.4x | PASS |
| SISSER2 | 2 | 0 | Acceptable | Optimal | 1.98e-11 | 21 | 20 | 14us | 2.4ms | 167.9x | PASS |
| SNAIL | 2 | 0 | Optimal | Optimal | 3.13e-21 | 10 | 63 | 8us | 9.4ms | 1173.4x | PASS |
| SNAKE | 2 | 2 | Optimal | Optimal | 2.00e-04 | 5 | 8 | 33us | 1.7ms | 51.9x | MISMATCH |
| SPANHYD | 97 | 33 | Optimal | Optimal | 0.00e+00 | 65 | 20 | 5.2ms | 5.2ms | 1.0x | PASS |
| SPIRAL | 3 | 2 | Optimal | Infeasible | N/A | 121 | 370 | 1.7ms | 61.6ms | 36.3x | ipopt_FAIL |
| SSI | 3 | 0 | Acceptable | MaxIteration | N/A | 2999 | 3000 | 5.1ms | 459.6ms | 89.4x | ipopt_FAIL |
| SSINE | 3 | 2 | Acceptable | Optimal | 0.00e+00 | 2999 | 224 | 30.5ms | 35.0ms | 1.1x | PASS |
| STANCMIN | 3 | 2 | Optimal | Optimal | 9.32e-09 | 10 | 9 | 31us | 1.6ms | 53.1x | PASS |
| STRATEC | 10 | 0 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| STREG | 4 | 0 | Optimal | Optimal | 8.90e-02 | 71 | 13 | 41us | 2.2ms | 53.2x | MISMATCH |
| STREGNE | 4 | 2 | Optimal | Optimal | 0.00e+00 | 2 | 2 | 88us | 454us | 5.1x | PASS |
| SUPERSIM | 2 | 2 | Optimal | Optimal | 2.22e-16 | 7 | 1 | 70us | 842us | 11.9x | PASS |
| SWOPF | 83 | 92 | Optimal | Optimal | 9.19e-10 | 15 | 13 | 1.6ms | 3.7ms | 2.3x | PASS |
| SYNTHES1 | 6 | 6 | Optimal | Optimal | 2.39e-07 | 25 | 8 | 2.8ms | 1.6ms | 0.6x | PASS |
| SYNTHES2 | 11 | 14 | Optimal | Optimal | 8.24e-07 | 21 | 14 | 181us | 3.2ms | 17.7x | PASS |
| SYNTHES3 | 17 | 23 | Optimal | Optimal | 3.99e-08 | 36 | 13 | 707us | 2.6ms | 3.7x | PASS |
| TAME | 2 | 1 | Optimal | Optimal | 0.00e+00 | 5 | 5 | 14us | 856us | 59.5x | PASS |
| TAX13322 | 72 | 1261 | Optimal | N/A | N/A | 103 | 0 | 154.9ms | N/A | N/A | ipopt_FAIL |
| TAXR13322 | 72 | 1261 | Timeout | Timeout | N/A | 0 | 0 | 60.00s | 60.00s | 1.0x | BOTH_FAIL |
| TENBARS1 | 18 | 9 | Optimal | Optimal | 9.73e-09 | 192 | 39 | 3.1ms | 7.4ms | 2.4x | PASS |
| TENBARS2 | 18 | 8 | Optimal | Optimal | 1.00e-08 | 27 | 33 | 314us | 6.8ms | 21.7x | PASS |
| TENBARS3 | 18 | 8 | Optimal | Optimal | 1.01e-08 | 25 | 34 | 284us | 7.1ms | 25.1x | PASS |
| TENBARS4 | 18 | 9 | Optimal | Optimal | 1.01e-10 | 14 | 14 | 209us | 3.3ms | 15.7x | PASS |
| TFI1 | 3 | 101 | Optimal | Optimal | 1.26e-10 | 112 | 19 | 3.3ms | 7.5ms | 2.3x | PASS |
| TFI2 | 3 | 101 | Optimal | Optimal | 5.14e-04 | 83 | 8 | 2.1ms | 2.7ms | 1.3x | MISMATCH |
| TFI3 | 3 | 101 | Optimal | Optimal | 6.13e-09 | 89 | 13 | 2.0ms | 4.2ms | 2.1x | PASS |
| THURBER | 7 | 37 | LocalInfeasi | IpoptStatus( | N/A | 16 | 0 | 476.4ms | 153us | 0.0x | BOTH_FAIL |
| THURBERLS | 7 | 0 | Acceptable | Optimal | 1.61e-16 | 525 | 19 | 141.1ms | 3.2ms | 0.0x | PASS |
| TOINTGOR | 50 | 0 | Acceptable | Optimal | 5.19e-13 | 105 | 7 | 284us | 1.1ms | 4.0x | PASS |
| TOINTPSP | 50 | 0 | Acceptable | Optimal | 1.90e-14 | 91 | 20 | 200us | 5.1ms | 25.5x | PASS |
| TOINTQOR | 50 | 0 | Acceptable | Optimal | 9.67e-16 | 37 | 1 | 90us | 426us | 4.7x | PASS |
| TRIGGER | 7 | 6 | Acceptable | Optimal | 0.00e+00 | 299 | 15 | 1.2ms | 2.5ms | 2.0x | PASS |
| TRO3X3 | 30 | 13 | Acceptable | Optimal | 3.61e-03 | 383 | 47 | 71.3ms | 10.8ms | 0.2x | MISMATCH |
| TRO4X4 | 63 | 25 | Acceptable | IpoptStatus( | N/A | 567 | 157 | 364.1ms | 47.5ms | 0.1x | ipopt_FAIL |
| TRO6X2 | 45 | 21 | Acceptable | RestorationF | N/A | 3048 | 353 | 530.8ms | 102.6ms | 0.2x | ipopt_FAIL |
| TRUSPYR1 | 11 | 4 | Optimal | Optimal | 2.03e-08 | 30 | 10 | 243us | 2.2ms | 8.9x | PASS |
| TRUSPYR2 | 11 | 11 | Acceptable | Optimal | 5.18e-08 | 141 | 13 | 14.9ms | 2.6ms | 0.2x | PASS |
| TRY-B | 2 | 1 | Optimal | Optimal | 5.01e-19 | 13 | 23 | 32us | 4.6ms | 144.4x | PASS |
| TWOBARS | 2 | 2 | Optimal | Optimal | 7.00e-01 | 19 | 8 | 107us | 1.6ms | 14.6x | MISMATCH |
| VESUVIA | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 2999 | 0 | 18.27s | 289us | 0.0x | BOTH_FAIL |
| VESUVIALS | 8 | 0 | Acceptable | Optimal | 3.39e-01 | 330 | 48 | 116.9ms | 19.9ms | 0.2x | MISMATCH |
| VESUVIO | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 105 | 0 | 373.9ms | 394us | 0.0x | BOTH_FAIL |
| VESUVIOLS | 8 | 0 | Acceptable | Optimal | 0.00e+00 | 520 | 10 | 239.9ms | 5.0ms | 0.0x | PASS |
| VESUVIOU | 8 | 1025 | LocalInfeasi | IpoptStatus( | N/A | 141 | 0 | 69.6ms | 286us | 0.0x | BOTH_FAIL |
| VESUVIOULS | 8 | 0 | Acceptable | Optimal | 5.55e-16 | 22 | 8 | 25.7ms | 4.7ms | 0.2x | PASS |
| VIBRBEAM | 8 | 0 | Optimal | Optimal | 9.48e-01 | 115 | 58 | 8.9ms | 10.4ms | 1.2x | MISMATCH |
| VIBRBEAMNE | 8 | 30 | LocalInfeasi | IpoptStatus( | N/A | 35 | 0 | 13.0ms | 279us | 0.0x | BOTH_FAIL |
| WACHBIEG | 3 | 2 | Optimal | Infeasible | N/A | 12 | 15 | 114us | 4.1ms | 35.5x | ipopt_FAIL |
| WATER | 31 | 10 | Acceptable | Optimal | 9.14e-08 | 34 | 17 | 15.9ms | 5.2ms | 0.3x | PASS |
| WAYSEA1 | 2 | 0 | Optimal | Optimal | 2.69e-15 | 34 | 14 | 20us | 2.3ms | 113.1x | PASS |
| WAYSEA1B | 2 | 0 | Optimal | Optimal | 4.54e-16 | 31 | 14 | 21us | 2.7ms | 129.5x | PASS |
| WAYSEA1NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 15 | 7 | 35us | 977us | 28.1x | PASS |
| WAYSEA2 | 2 | 0 | Optimal | Optimal | 9.84e-18 | 19 | 22 | 12us | 2.8ms | 220.2x | PASS |
| WAYSEA2B | 2 | 0 | Optimal | Optimal | 9.15e-18 | 17 | 22 | 10us | 3.4ms | 333.3x | PASS |
| WAYSEA2NE | 2 | 2 | Optimal | Optimal | 0.00e+00 | 22 | 11 | 46us | 1.5ms | 32.6x | PASS |
| WEEDS | 3 | 0 | Acceptable | Optimal | 1.99e-14 | 225 | 28 | 463us | 6.9ms | 15.0x | PASS |
| WEEDSNE | 3 | 12 | LocalInfeasi | IpoptStatus( | N/A | 70 | 0 | 589us | 100us | 0.2x | BOTH_FAIL |
| WOMFLET | 3 | 3 | Optimal | Optimal | 1.00e+00 | 13 | 8 | 92us | 1.6ms | 17.3x | MISMATCH |
| YFIT | 3 | 0 | Optimal | Optimal | 1.33e-19 | 77 | 36 | 82us | 7.6ms | 92.9x | PASS |
| YFITNE | 3 | 17 | Acceptable | IpoptStatus( | N/A | 36 | 0 | 209us | 100us | 0.5x | ipopt_FAIL |
| YFITU | 3 | 0 | Optimal | Optimal | 5.31e-21 | 77 | 36 | 83us | 5.6ms | 68.2x | PASS |
| ZANGWIL2 | 2 | 0 | Optimal | Optimal | 0.00e+00 | 2 | 1 | 2us | 257us | 162.2x | PASS |
| ZANGWIL3 | 3 | 3 | Optimal | Optimal | 0.00e+00 | 1 | 1 | 10us | 349us | 36.5x | PASS |
| ZECEVIC2 | 2 | 2 | Optimal | Optimal | 5.57e-10 | 11 | 8 | 32us | 1.4ms | 45.7x | PASS |
| ZECEVIC3 | 2 | 2 | Optimal | Optimal | 8.24e-10 | 12 | 17 | 47us | 3.3ms | 70.5x | PASS |
| ZECEVIC4 | 2 | 2 | Optimal | Optimal | 2.62e-09 | 11 | 10 | 37us | 2.2ms | 61.4x | PASS |
| ZY2 | 3 | 2 | Acceptable | Optimal | 4.22e-05 | 23 | 14 | 1.4ms | 3.6ms | 2.6x | PASS |

## Performance Comparison (where both solve)

### Iteration Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 166.8 | 44.9 |
| Median | 18 | 13 |
| Total  | 92598 | 24914 |

- ripopt fewer iterations: 121/555
- Ipopt fewer iterations: 386/555
- Tied: 48/555

### Timing Comparison

| Metric | ripopt | Ipopt |
|--------|--------|-------|
| Mean   | 120.9ms | 19.5ms |
| Median | 117us | 2.4ms |
| Total  | 67.12s | 10.82s |

- Geometric mean speedup (Ipopt_time/ripopt_time): **8.22x**
  - \>1 means ripopt is faster, <1 means Ipopt is faster
- ripopt faster: 433/555 problems
- Ipopt faster: 122/555 problems
- Overall speedup (total time): 0.16x

## Failure Analysis

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
| EQC | 9 | 3 | ErrorInStepComputation | -8.367978e+02 |
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
| MESH | 41 | 48 | IpoptStatus(4) | -1.129753e+06 |
| NYSTROM5 | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| NYSTROM5C | 18 | 20 | IpoptStatus(-10) | 0.000000e+00 |
| PALMER5A | 8 | 0 | MaxIterations | 3.750621e-02 |
| PALMER5E | 8 | 0 | MaxIterations | 2.128087e+00 |
| PALMER7A | 6 | 0 | MaxIterations | 1.033779e+01 |
| PALMER7E | 8 | 0 | MaxIterations | 6.697960e+00 |
| PFIT1 | 3 | 3 | Infeasible | 0.000000e+00 |
| PFIT2 | 3 | 3 | RestorationFailed | 0.000000e+00 |
| POLAK3 | 12 | 10 | MaxIterations | 7.084571e+00 |
| POLAK6 | 5 | 4 | MaxIterations | -1.494339e+01 |
| POWELLSQ | 2 | 2 | Infeasible | 0.000000e+00 |
| ROBOT | 14 | 2 | IpoptStatus(3) | 6.593299e+00 |
| SPIRAL | 3 | 2 | Infeasible | 8.687279e-09 |
| SSI | 3 | 0 | MaxIterations | 1.510022e-09 |
| TAX13322 | 72 | 1261 | N/A | -3.726455e+03 |
| TRO4X4 | 63 | 25 | IpoptStatus(4) | 8.999846e+00 |
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
| RAT42 | 3 | 9 | LocalInfeasibility | IpoptStatus(-10) |
| RAT43 | 4 | 15 | LocalInfeasibility | IpoptStatus(-10) |
| ROSZMAN1 | 4 | 25 | LocalInfeasibility | IpoptStatus(-10) |
| S308NE | 2 | 3 | LocalInfeasibility | IpoptStatus(-10) |
| S365 | 7 | 5 | MaxIterations | RestorationFailed |
| S365MOD | 7 | 5 | MaxIterations | RestorationFailed |
| SANTA | 21 | 23 | LocalInfeasibility | IpoptStatus(-10) |
| STRATEC | 10 | 0 | Timeout | Timeout |
| TAXR13322 | 72 | 1261 | Timeout | Timeout |
| THURBER | 7 | 37 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIA | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIO | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VESUVIOU | 8 | 1025 | LocalInfeasibility | IpoptStatus(-10) |
| VIBRBEAMNE | 8 | 30 | LocalInfeasibility | IpoptStatus(-10) |
| WEEDSNE | 3 | 12 | LocalInfeasibility | IpoptStatus(-10) |

### Objective mismatches (101)

Both solvers converged but found different objective values (rel diff > 1e-4).

- **Different local minimum** (both Optimal): 63
- **Convergence gap** (one Acceptable): 38
- **Better objective found by**: ripopt 25, Ipopt 76

| Problem | ripopt obj | Ipopt obj | Rel Diff | r_status | i_status | Better |
|---------|-----------|-----------|----------|----------|----------|--------|
| MAKELA1 | 1.414214e+00 | -1.414214e+00 | 2.00e+00 | Optimal | Optimal | ipopt |
| SIPOW1 | 4.338956e-01 | -1.000000e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| SIPOW1M | 4.338960e-01 | -1.000001e+00 | 1.43e+00 | Optimal | Optimal | ipopt |
| LEVYMONT5 | 7.773811e+00 | 1.239502e-25 | 1.00e+00 | Acceptable | Optimal | ipopt |
| LEVYMONT6 | 1.073629e-22 | 1.250286e+01 | 1.00e+00 | Optimal | Optimal | ripopt |
| WOMFLET | 2.313759e-08 | 6.050000e+00 | 1.00e+00 | Optimal | Optimal | ripopt |
| MGH10SLS | 1.417866e+09 | 8.794586e+01 | 1.00e+00 | Acceptable | Optimal | ipopt |
| HS268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| S268 | 3.182746e+00 | 8.886855e-07 | 1.00e+00 | Optimal | Optimal | ipopt |
| DANIWOODLS | 1.039178e+02 | 4.317308e-03 | 1.00e+00 | Optimal | Optimal | ipopt |
| ELATTAR | 4.781305e+05 | 7.420618e+01 | 1.00e+00 | Acceptable | Optimal | ipopt |
| GROWTHLS | 3.542149e+03 | 1.004041e+00 | 1.00e+00 | Optimal | Optimal | ipopt |
| SIPOW3 | 5.582563e+02 | 5.346586e-01 | 9.99e-01 | Optimal | Optimal | ipopt |
| RAT43LS | 1.076462e+06 | 8.786405e+03 | 9.92e-01 | Optimal | Optimal | ipopt |
| HS16 | 2.314466e+01 | 2.500000e-01 | 9.89e-01 | Optimal | Optimal | ipopt |
| KSIP | 2.768715e+01 | 5.757979e-01 | 9.79e-01 | Optimal | Optimal | ipopt |
| HALDMADS | 4.764200e-02 | 2.218282e+00 | 9.79e-01 | Acceptable | Optimal | ripopt |
| LEVYMONT8 | 4.980481e+00 | 1.734956e+02 | 9.71e-01 | Optimal | Optimal | ripopt |
| ELATVIDU | 1.712780e+00 | 5.475112e+01 | 9.69e-01 | Optimal | Optimal | ripopt |
| ELATVIDUB | 1.712780e+00 | 5.475112e+01 | 9.69e-01 | Acceptable | Optimal | ripopt |
| MWRIGHT | 1.288383e+00 | 2.497881e+01 | 9.48e-01 | Optimal | Optimal | ripopt |
| VIBRBEAM | 6.346254e+00 | 3.322376e-01 | 9.48e-01 | Optimal | Optimal | ipopt |
| JENSMP | 2.020000e+03 | 1.243622e+02 | 9.38e-01 | Optimal | Optimal | ipopt |
| BT4 | -4.551055e+01 | -3.704768e+00 | 9.19e-01 | Optimal | Optimal | ripopt |
| BOXBODLS | 1.168009e+03 | 9.771500e+03 | 8.80e-01 | Acceptable | Optimal | ripopt |
| SIPOW4 | 2.080345e+00 | 2.723620e-01 | 8.69e-01 | Optimal | Optimal | ipopt |
| HIMMELP2 | -6.205394e+01 | -8.198044e+00 | 8.68e-01 | Optimal | Optimal | ripopt |
| HIMMELP3 | -7.913699e+00 | -5.901318e+01 | 8.66e-01 | Optimal | Optimal | ipopt |
| HS23 | 9.473230e+00 | 2.000000e+00 | 7.89e-01 | Optimal | Optimal | ipopt |
| HS54 | -1.566691e-01 | -9.080748e-01 | 7.51e-01 | Optimal | Optimal | ipopt |
| HIMMELP6 | -1.475339e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| HIMMELP5 | -1.475901e+01 | -5.901318e+01 | 7.50e-01 | Optimal | Optimal | ipopt |
| TWOBARS | 5.028823e+00 | 1.508652e+00 | 7.00e-01 | Optimal | Optimal | ipopt |
| ECKERLE4LS | 6.996961e-01 | 1.463589e-03 | 6.98e-01 | Acceptable | Optimal | ipopt |
| MAKELA4 | 6.877365e-01 | -9.600000e-09 | 6.88e-01 | Optimal | Optimal | ipopt |
| LEVYMONT9 | 5.367398e+01 | 1.661496e+02 | 6.77e-01 | Optimal | Optimal | ripopt |
| HS85 | -1.237175e+00 | -2.215605e+00 | 4.42e-01 | Optimal | Optimal | ipopt |
| HS108 | -8.660254e-01 | -5.000000e-01 | 3.66e-01 | Optimal | Optimal | ripopt |
| MISTAKE | -6.571639e-01 | -1.000000e+00 | 3.43e-01 | Optimal | Optimal | ipopt |
| VESUVIALS | 1.500440e+03 | 9.914100e+02 | 3.39e-01 | Acceptable | Optimal | ipopt |
| HYDC20LS | 2.977216e-01 | 2.967522e-15 | 2.98e-01 | Acceptable | Optimal | ipopt |
| HAIFAS | -2.000000e-01 | -4.500000e-01 | 2.50e-01 | Optimal | Optimal | ipopt |
| JUDGEB | 2.048234e+01 | 1.608173e+01 | 2.15e-01 | Optimal | Optimal | ipopt |
| JUDGE | 2.048234e+01 | 1.608173e+01 | 2.15e-01 | Optimal | Optimal | ipopt |
| EG1 | -1.132801e+00 | -1.429307e+00 | 2.07e-01 | Optimal | Optimal | ipopt |
| LEVYMONT7 | 1.015401e+01 | 1.251076e+01 | 1.88e-01 | Optimal | Optimal | ripopt |
| DENSCHNC | 1.833617e-01 | 2.177679e-20 | 1.83e-01 | Optimal | Optimal | ipopt |
| HS70 | 1.877383e-01 | 7.498464e-03 | 1.80e-01 | Optimal | Optimal | ipopt |
| ACOPR30 | 6.751849e+02 | 5.768924e+02 | 1.46e-01 | Acceptable | Optimal | ipopt |
| MUONSINELS | 5.114204e+04 | 4.387412e+04 | 1.42e-01 | Acceptable | Optimal | ipopt |
| MAKELA2 | 8.244898e+00 | 7.200000e+00 | 1.27e-01 | Optimal | Optimal | ipopt |
| MSS1 | -1.600000e+01 | -1.400000e+01 | 1.25e-01 | Acceptable | Optimal | ripopt |
| QC | -8.671253e+02 | -9.565379e+02 | 9.35e-02 | Acceptable | Optimal | ipopt |
| HYDCAR6LS | 8.942843e-02 | 2.863560e-18 | 8.94e-02 | Acceptable | Optimal | ipopt |
| STREG | 1.110999e-23 | 8.901950e-02 | 8.90e-02 | Optimal | Optimal | ripopt |
| HAHN1LS | 3.086398e+01 | 3.338424e+01 | 7.55e-02 | Acceptable | Optimal | ripopt |
| PALMER3 | 2.265958e+03 | 2.416980e+03 | 6.25e-02 | Acceptable | Optimal | ripopt |
| PALMER4 | 2.424016e+03 | 2.285383e+03 | 5.72e-02 | Optimal | Optimal | ipopt |
| OET1 | 5.834999e-01 | 5.382431e-01 | 4.53e-02 | Optimal | Optimal | ipopt |
| OET2 | 1.317369e-01 | 8.715962e-02 | 4.46e-02 | Acceptable | Optimal | ipopt |
| ROSZMAN1LS | 3.951903e-02 | 4.948485e-04 | 3.90e-02 | Acceptable | Optimal | ipopt |
| DEVGLA2 | 2.501146e-02 | 6.672171e-19 | 2.50e-02 | Acceptable | Optimal | ipopt |
| GIGOMEZ2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CHACONN1 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| CB2 | 2.000000e+00 | 1.952224e+00 | 2.39e-02 | Optimal | Optimal | ipopt |
| HS55 | 6.666667e+00 | 6.805833e+00 | 2.04e-02 | Optimal | Optimal | ripopt |
| MAXLIKA | 1.149346e+03 | 1.136307e+03 | 1.13e-02 | Acceptable | Optimal | ipopt |
| EXPFITC | 3.272073e-02 | 2.330262e-02 | 9.42e-03 | Optimal | Optimal | ipopt |
| PT | 1.873320e-01 | 1.783942e-01 | 8.94e-03 | Optimal | Optimal | ipopt |
| QCNEW | -8.137719e+02 | -8.065219e+02 | 8.91e-03 | Acceptable | Optimal | ripopt |
| SIPOW2 | -9.922730e-01 | -1.000000e+00 | 7.73e-03 | Optimal | Optimal | ipopt |
| SIPOW2M | -9.923838e-01 | -1.000005e+00 | 7.62e-03 | Optimal | Optimal | ipopt |
| CLIFF | 1.997866e-01 | 2.072380e-01 | 7.45e-03 | Optimal | Optimal | ripopt |
| LRIJCNN1 | 2.730781e-01 | 2.671576e-01 | 5.92e-03 | Acceptable | Optimal | ipopt |
| BIGGS5 | 5.655650e-03 | 1.088200e-19 | 5.66e-03 | Optimal | Optimal | ipopt |
| BIGGS6 | 5.655650e-03 | 2.316725e-22 | 5.66e-03 | Optimal | Optimal | ipopt |
| PALMER5B | 1.524120e-02 | 9.752496e-03 | 5.49e-03 | Acceptable | Optimal | ipopt |
| HS109 | 5.391541e+03 | 5.362069e+03 | 5.47e-03 | Acceptable | Optimal | ipopt |
| DGOSPEC | -9.887540e+02 | -9.933506e+02 | 4.63e-03 | Acceptable | Optimal | ipopt |
| TRO3X3 | 9.000000e+00 | 8.967478e+00 | 3.61e-03 | Acceptable | Optimal | ipopt |
| ACOPP30 | 5.788980e+02 | 5.768924e+02 | 3.46e-03 | Optimal | Optimal | ipopt |
| DEGENLPB | -3.069162e+01 | -3.076401e+01 | 2.35e-03 | Acceptable | Optimal | ipopt |
| LIN | -1.960628e-02 | -1.757754e-02 | 2.03e-03 | Optimal | Optimal | ripopt |
| EXPFITB | 7.016641e-03 | 5.019367e-03 | 2.00e-03 | Optimal | Optimal | ipopt |
| HET-Z | 1.001858e+00 | 1.000000e+00 | 1.85e-03 | Optimal | Optimal | ipopt |
| DEGENLPA | 3.060434e+00 | 3.054881e+00 | 1.81e-03 | Acceptable | Optimal | ipopt |
| HS13 | 9.933710e-01 | 9.945785e-01 | 1.21e-03 | Acceptable | Optimal | ripopt |
| DECONVC | 3.734409e-03 | 2.569475e-03 | 1.16e-03 | Optimal | Optimal | ipopt |
| MGH09LS | 1.019673e-03 | 3.075056e-04 | 7.12e-04 | Acceptable | Optimal | ipopt |
| OET3 | 5.052162e-03 | 4.505043e-03 | 5.47e-04 | Optimal | Optimal | ipopt |
| TFI2 | 6.495452e-01 | 6.490311e-01 | 5.14e-04 | Optimal | Optimal | ipopt |
| HS116 | 9.754427e+01 | 9.758747e+01 | 4.43e-04 | Acceptable | Optimal | ripopt |
| ANTWERP | 3.246598e+03 | 3.245241e+03 | 4.18e-04 | Acceptable | Optimal | ipopt |
| HS45 | 1.000314e+00 | 1.000000e+00 | 3.14e-04 | Acceptable | Optimal | ipopt |
| LSC2LS | 1.333749e+01 | 1.333439e+01 | 2.32e-04 | Acceptable | Optimal | ipopt |
| DENSCHND | 3.820110e-09 | 2.221899e-04 | 2.22e-04 | Acceptable | Optimal | ripopt |
| HS17 | 1.000201e+00 | 1.000000e+00 | 2.01e-04 | Optimal | Optimal | ipopt |
| SNAKE | -2.601155e-10 | -1.999999e-04 | 2.00e-04 | Optimal | Optimal | ipopt |
| ACOPR14 | 8.082839e+03 | 8.081526e+03 | 1.62e-04 | Acceptable | Optimal | ipopt |
| HATFLDH | -2.450304e+01 | -2.450000e+01 | 1.24e-04 | Acceptable | Optimal | ripopt |
| HS96 | 1.573793e-02 | 1.561775e-02 | 1.20e-04 | Acceptable | Optimal | ipopt |

---
*Generated by cutest_suite/compare.py*