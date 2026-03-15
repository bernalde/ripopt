# Changelog

## [0.4.0] - 2026-03-15

### Added
- **Mehrotra predictor-corrector** enabled by default (`mehrotra_pc=true`), with Gondzio centrality corrections (`gondzio_mcc_max=3`) for better centering and fewer iterations
- **Dense condensed KKT for tall-narrow problems**: when m >> n and n <= 100, uses an n x n dense solve instead of (n+m) x (n+m) sparse factorization (up to 845x faster on problems like EXPFITC)
- **SuiteSparse AMD ordering** for rmumps: replaced custom O(n^2) AMD with the `amd` crate, fixing 10s+ ordering times on 40K+ dimensional systems
- **Early stall detection** (`early_stall_timeout=10.0`): bails out fast when stuck in early iterations to trigger fallback strategies
- **CLI `--help` flag**: lists all 50+ solver options organized by category with types, defaults, and descriptions
- **Quality function mu oracle** (disabled by default, `mu_oracle_quality_function`): evaluates barrier KKT error for candidate mu values
- **Large-scale benchmark with Ipopt comparison**: both solvers receive the same NlpProblem via Rust trait, up to 100K variables
- **faer** sparse solver as optional feature (default enabled alongside rmumps)

### Changed
- Default sparse direct solver switched back to rmumps (from faer) after fixing AMD ordering — rmumps multifrontal factorization is 5x faster than faer SparseLdl
- Mehrotra centering parameter (sigma) now feeds into cross-iteration mu update via geometric blend
- Sufficient progress check tightened (`refs_red_fact`: 0.9999 -> 0.999) for earlier Free-to-Fixed mode switching
- Wall-time check runs every iteration in early phase (was every 10)
- Skip expensive NLP restoration when approaching early stall timeout

### Fixed
- Compiler warnings: unused assignment (`tried_compl_polish`), dead code (`HsProblemEntry`)
- Missing fields in `HsSolveResult` for ipopt_native benchmark macro
- Sparse condensed Schur complement with near-full bandwidth (> n/2) now falls back to augmented KKT instead of attempting dense-equivalent sparse factorization

### Performance
- SparseQP 100K: 25.4s -> 4.9s (SuiteSparse AMD + rmumps vs faer)
- EXPFITC (n=5, m=502): 10.1s -> 0.012s (dense condensed path)
- OET3 (n=4, m=1002): 1.4s -> 0.009s (dense condensed path)
- cho_parmest (n=21672, m=21660): first factorization 10.9s -> 0.066s (faer AMD)
- ACOPR14: 60s timeout -> 0.4s (early stall detection + fallback)
- CUTEst geometric mean speedup vs Ipopt: 7.8x -> 9.1x
- HS suite geometric mean speedup vs Ipopt: 15.0x -> 16.5x

## [0.3.0] - 2026-02-14

Initial public release with full IPM implementation, CUTEst/HS benchmarks, C API, AMPL interface, and Pyomo integration.
