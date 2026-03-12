//! # rmumps
//!
//! A pure Rust multifrontal sparse symmetric indefinite (LDL^T) solver.
//!
//! rmumps factorizes sparse symmetric matrices using a multifrontal method with
//! Bunch-Kaufman pivoting, providing inertia detection (counts of positive, negative,
//! and zero eigenvalues in the D factor). It is designed for solving KKT systems
//! arising in interior point methods for nonlinear optimization.
//!
//! ## Features
//!
//! - **Multifrontal factorization**: organizes sparse LDL^T into dense operations on
//!   frontal matrices for better cache utilization than simplicial methods
//! - **Bunch-Kaufman pivoting**: 1x1 and 2x2 pivots for symmetric indefinite matrices
//! - **Inertia detection**: counts positive/negative/zero eigenvalues from D factor
//! - **AMD ordering**: approximate minimum degree fill-reducing permutation
//! - **Supernodal elimination tree**: groups columns into supernodes for efficiency
//! - **Iterative refinement**: configurable refinement steps for improved accuracy
//! - **Cached symbolic analysis**: analyze once, refactor with new values (same pattern)
//! - **Parallel level-set factorization**: independent supernodes factored via rayon
//!
//! ## Usage
//!
//! ```rust
//! use rmumps::coo::CooMatrix;
//! use rmumps::solver::{Solver, SolverOptions};
//!
//! // Create a 3x3 symmetric positive definite matrix (upper triangle, COO format)
//! let coo = CooMatrix::new(3,
//!     vec![0, 1, 2, 0, 1],  // rows
//!     vec![0, 1, 2, 1, 2],  // cols (col >= row)
//!     vec![4.0, 5.0, 6.0, 1.0, 2.0],  // values
//! ).unwrap();
//!
//! let mut solver = Solver::new(SolverOptions::default());
//! let inertia = solver.analyze_and_factor(&coo).unwrap();
//! assert_eq!(inertia.positive, 3);
//!
//! let rhs = vec![1.0, 2.0, 3.0];
//! let mut solution = vec![0.0; 3];
//! solver.solve(&rhs, &mut solution).unwrap();
//! ```
//!
//! ## Architecture
//!
//! The solver follows a three-phase approach:
//!
//! 1. **Symbolic analysis** (`symbolic.rs`): computes elimination tree, supernodes,
//!    and frontal matrix structure from the sparsity pattern
//! 2. **Numeric factorization** (`numeric.rs`): assembles and partially factors
//!    frontal matrices bottom-up through the elimination tree
//! 3. **Solve** (`solve.rs`): forward/backward substitution through the supernodal
//!    factorization, with optional iterative refinement
//!
//! Key modules:
//! - `coo`: COO (triplet) sparse matrix format
//! - `csc`: CSC (compressed sparse column) format with COO conversion
//! - `ordering`: fill-reducing orderings (AMD, natural)
//! - `etree`: elimination tree construction (Liu 1990)
//! - `symbolic`: supernodal symbolic factorization
//! - `frontal`: frontal matrix assembly and partial factorization
//! - `dense`: column-major dense matrix with BLAS-like operations
//! - `pivot`: Bunch-Kaufman pivoting for dense symmetric indefinite matrices
//! - `numeric`: multifrontal numeric factorization with inertia tracking
//! - `solve`: triangular solve with iterative refinement
//! - `solver`: high-level API combining all phases

/// Sparse symmetric matrix in COO (triplet) format.
pub mod coo;
/// Sparse symmetric matrix in CSC (compressed sparse column) format.
pub mod csc;
/// Column-major dense matrix with BLAS-like operations.
pub mod dense;
/// Elimination tree construction from sparse matrix structure.
pub mod etree;
/// Fill-reducing orderings (AMD, natural).
pub mod ordering;
/// Bunch-Kaufman pivoting for dense symmetric indefinite matrices.
pub mod pivot;
/// Frontal matrix assembly and partial factorization.
pub mod frontal;
/// Multifrontal numeric factorization with inertia tracking.
pub mod numeric;
/// Triangular solve with iterative refinement.
pub mod solve;
/// High-level solver API combining symbolic analysis, numeric factorization, and solve.
pub mod solver;
/// Matrix scaling (diagonal equilibration, Ruiz) for improved pivot quality.
pub mod scaling;
/// Supernodal symbolic factorization from sparsity pattern.
pub mod symbolic;
/// Preconditioner trait and basic implementations.
pub mod precond;
/// MINRES iterative solver for symmetric systems.
pub mod minres;
/// Incomplete LDL^T factorization for preconditioning.
pub mod incomplete;

use std::fmt;

/// Inertia of a symmetric matrix after LDL^T factorization.
/// Counts the signs of diagonal entries in D.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Inertia {
    /// Number of positive eigenvalues.
    pub positive: usize,
    /// Number of negative eigenvalues.
    pub negative: usize,
    /// Number of zero eigenvalues.
    pub zero: usize,
}

impl fmt::Display for Inertia {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(+{}, -{}, 0:{})", self.positive, self.negative, self.zero)
    }
}

/// Error from the solver.
#[derive(Debug, Clone)]
pub enum SolverError {
    /// Matrix is structurally singular.
    SingularMatrix,
    /// Numerical failure during factorization.
    NumericalFailure(String),
    /// Dimension mismatch.
    DimensionMismatch {
        /// Expected dimension.
        expected: usize,
        /// Actual dimension.
        got: usize,
    },
    /// Invalid input (e.g., bad COO indices).
    InvalidInput(String),
    /// Solver not in correct state (e.g., solve before factor).
    InvalidState(String),
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::SingularMatrix => write!(f, "singular matrix"),
            SolverError::NumericalFailure(msg) => write!(f, "numerical failure: {}", msg),
            SolverError::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: expected {}, got {}", expected, got)
            }
            SolverError::InvalidInput(msg) => write!(f, "invalid input: {}", msg),
            SolverError::InvalidState(msg) => write!(f, "invalid state: {}", msg),
        }
    }
}

impl std::error::Error for SolverError {}
