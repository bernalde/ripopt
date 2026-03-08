pub mod coo;
pub mod csc;
pub mod dense;
pub mod etree;
pub mod ordering;
pub mod pivot;
pub mod frontal;
pub mod numeric;
pub mod solve;
pub mod solver;
pub mod symbolic;

use std::fmt;

/// Inertia of a symmetric matrix after LDL^T factorization.
/// Counts the signs of diagonal entries in D.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Inertia {
    pub positive: usize,
    pub negative: usize,
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
    DimensionMismatch { expected: usize, got: usize },
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
