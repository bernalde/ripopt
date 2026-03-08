# rmumps

A pure Rust multifrontal sparse symmetric indefinite (LDL^T) solver.

## Overview

rmumps factorizes sparse symmetric matrices A = LDL^T using a multifrontal method with Bunch-Kaufman pivoting. It provides inertia detection (counts of positive, negative, and zero eigenvalues in D), which is essential for KKT systems in interior point optimization.

The name reflects its inspiration from [MUMPS](http://mumps-solver.org/) (MUltifrontal Massively Parallel sparse direct Solver), a widely-used Fortran sparse solver. rmumps is a MUMPS-inspired, LLM-guided Rust implementation that follows the same algorithmic approach (multifrontal with supernodal elimination trees). The LLM (Claude Code + Opus 4.6) that guided the implementation was likely trained on MUMPS and related sparse solver codebases, so indirect influence should be assumed. rmumps inherits the MUMPS license as it could be considered a derivative work.

## Features

- Multifrontal LDL^T factorization with supernodal elimination tree
- Bunch-Kaufman pivoting (1x1 and 2x2 pivots) for symmetric indefinite matrices
- Inertia detection from the D factor
- Approximate Minimum Degree (AMD) fill-reducing ordering
- Cached symbolic analysis (analyze once, refactor with new values)
- Parallel level-set factorization via rayon
- Iterative refinement for improved solution accuracy
- Matrix scaling (diagonal equilibration, Ruiz iterative) for ill-conditioned systems
- `min_diagonal()` for direct inertia correction shortcuts
- `factor_csc()` to skip COO-to-CSC conversion during refactorization

## Usage

```rust
use rmumps::coo::CooMatrix;
use rmumps::solver::{Solver, SolverOptions};

// Create a 3x3 symmetric matrix (upper triangle, COO format)
let coo = CooMatrix::new(3,
    vec![0, 1, 2, 0, 1],  // rows
    vec![0, 1, 2, 1, 2],  // cols (col >= row)
    vec![4.0, 5.0, 6.0, 1.0, 2.0],  // values
).unwrap();

let mut solver = Solver::new(SolverOptions::default());
let inertia = solver.analyze_and_factor(&coo).unwrap();
assert_eq!(inertia.positive, 3);

let rhs = vec![1.0, 2.0, 3.0];
let mut solution = vec![0.0; 3];
solver.solve(&rhs, &mut solution).unwrap();
```

## Three-phase workflow

1. **Symbolic analysis** (`analyze`): computes AMD ordering, elimination tree, and supernodes from the sparsity pattern. Done once per pattern.
2. **Numeric factorization** (`factor` or `factor_csc`): assembles frontal matrices and performs partial LDL^T with Bunch-Kaufman pivoting bottom-up through the tree. Can be repeated with new values.
3. **Solve** (`solve`): forward/backward substitution with optional iterative refinement.

## Differences from Fortran MUMPS

We did not attempt to implement all the features of MUMPS, nor make it a direct translation line by line. We aimed for features that make ripopt comparable with ipopt for large-scale sparse optimization problems. We summarize the features and differences below.

| Aspect       | rmumps               | MUMPS                        |
|--------------|----------------------|------------------------------|
| Language     | Pure Rust            | Fortran 90 + C               |
| Pivoting     | Bunch-Kaufman        | Threshold partial pivoting   |
| Parallelism  | rayon level-set      | MPI + OpenMP                 |
| Ordering     | AMD (built-in)       | METIS, SCOTCH, AMD, PORD     |
| Scaling      | Diagonal / Ruiz      | MC64 row/column scaling      |
| Out-of-core  | Not supported        | Supported                    |
| Matrix types | Symmetric only       | General, symmetric, SPD      |
| Precision    | f64 only             | Single, double, complex      |
| Memory       | Rust ownership model | Manual Fortran allocation    |
| Dependencies | rayon only           | BLAS, LAPACK, ScaLAPACK, MPI |

## Limitations

- **Symmetric matrices only.** rmumps handles symmetric indefinite (LDL^T) but not general unsymmetric (LU) or SPD (Cholesky) factorizations.
- **No distributed parallelism.** Level-set parallelism via rayon works well on shared-memory systems but does not scale across nodes like MUMPS with MPI.
- **Basic matrix scaling only.** rmumps supports diagonal equilibration and Ruiz iterative scaling. MUMPS uses MC64 maximum-weight matching which can handle more pathological cases.
- **Single precision not supported.** Only f64 (double precision).
- **No out-of-core mode.** The entire factorization must fit in memory.
- **AMD ordering only.** MUMPS supports METIS and SCOTCH which can produce better orderings for some problems, reducing fill-in and factorization time.
- **Performance gap on large problems.** Fortran MUMPS with optimized BLAS (OpenBLAS, MKL) is faster on large frontal matrices due to highly-tuned dense linear algebra kernels. rmumps uses hand-written dense operations without SIMD intrinsics.

## Performance

rmumps is used as the default sparse solver in [ripopt](https://github.com/jkitchin/ripopt). On ripopt's benchmark suite (Apple Silicon, single-threaded):

| Problem         | n + m   | rmumps time |
|-----------------|---------|-------------|
| Bratu 1K        | 1,998   | 0.002s      |
| OptControl 2.5K | 3,749   | 0.006s      |
| Bratu 10K       | 19,998  | 0.136s      |
| OptControl 20K  | 29,999  | 0.200s      |
| Poisson 50K     | 74,892  | 3.9s        |
| SparseQP 100K   | 100,000 | 4.8s        |

For comparison, Ipopt with Fortran MUMPS solves OptControl 20K in 0.02s and SparseQP 100K in 0.35s, reflecting the performance gap from optimized BLAS kernels. The gap is smaller on banded problems where rmumps auto-delegates to a specialized banded solver.

## License

[CeCILL-C](LICENSE) (LGPL-compatible free software license).

This is a different license from the parent ripopt project (EPL-2.0). CeCILL-C was chosen because that is how MUMPS is licensed, and this is likely considered a derivative work of that.
