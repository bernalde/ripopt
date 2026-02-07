pub mod dense;

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
        write!(
            f,
            "(+{}, -{}, 0:{})",
            self.positive, self.negative, self.zero
        )
    }
}

/// Error from a linear solver.
#[derive(Debug, Clone)]
pub enum SolverError {
    /// Matrix is structurally singular (e.g., zero row).
    SingularMatrix,
    /// Numerical failure during factorization.
    NumericalFailure(String),
    /// Dimension mismatch.
    DimensionMismatch { expected: usize, got: usize },
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::SingularMatrix => write!(f, "singular matrix"),
            SolverError::NumericalFailure(msg) => write!(f, "numerical failure: {}", msg),
            SolverError::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: expected {}, got {}", expected, got)
            }
        }
    }
}

impl std::error::Error for SolverError {}

/// Symmetric matrix stored as dense lower triangle, column-major.
#[derive(Debug, Clone)]
pub struct SymmetricMatrix {
    /// Dimension of the matrix.
    pub n: usize,
    /// Lower triangle entries stored column-major.
    /// For column j, rows j..n are stored.
    /// Entry (i,j) where i >= j is at index: j*n - j*(j-1)/2 + (i-j)
    pub data: Vec<f64>,
}

impl SymmetricMatrix {
    /// Create a new zero symmetric matrix of dimension n.
    pub fn zeros(n: usize) -> Self {
        let nnz = n * (n + 1) / 2;
        Self {
            n,
            data: vec![0.0; nnz],
        }
    }

    /// Get the index into data for entry (i, j) where i >= j.
    #[inline]
    fn packed_index(n: usize, i: usize, j: usize) -> usize {
        debug_assert!(i >= j);
        debug_assert!(i < n);
        j * n - j * (j + 1) / 2 + i
    }

    /// Get element (i, j), automatically handling symmetry.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i >= j {
            self.data[Self::packed_index(self.n, i, j)]
        } else {
            self.data[Self::packed_index(self.n, j, i)]
        }
    }

    /// Set element (i, j), automatically handling symmetry.
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        if i >= j {
            self.data[Self::packed_index(self.n, i, j)] = val;
        } else {
            self.data[Self::packed_index(self.n, j, i)] = val;
        }
    }

    /// Add val to element (i, j), automatically handling symmetry.
    pub fn add(&mut self, i: usize, j: usize, val: f64) {
        if i >= j {
            self.data[Self::packed_index(self.n, i, j)] += val;
        } else {
            self.data[Self::packed_index(self.n, j, i)] += val;
        }
    }

    /// Add delta to all diagonal entries.
    pub fn add_diagonal(&mut self, delta: f64) {
        for i in 0..self.n {
            self.data[Self::packed_index(self.n, i, i)] += delta;
        }
    }

    /// Add delta to diagonal entries in range [start, end).
    pub fn add_diagonal_range(&mut self, start: usize, end: usize, delta: f64) {
        for i in start..end {
            self.data[Self::packed_index(self.n, i, i)] += delta;
        }
    }

    /// Compute y = A * x (symmetric matrix-vector product).
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        let n = self.n;
        for i in 0..n {
            y[i] = 0.0;
        }
        for j in 0..n {
            // Diagonal
            let ajj = self.data[Self::packed_index(n, j, j)];
            y[j] += ajj * x[j];
            // Off-diagonal (lower triangle)
            for i in (j + 1)..n {
                let aij = self.data[Self::packed_index(n, i, j)];
                y[i] += aij * x[j];
                y[j] += aij * x[i];
            }
        }
    }

    /// Convert to full dense matrix (row-major) for debugging.
    pub fn to_full(&self) -> Vec<Vec<f64>> {
        let mut m = vec![vec![0.0; self.n]; self.n];
        for i in 0..self.n {
            for j in 0..=i {
                let v = self.get(i, j);
                m[i][j] = v;
                m[j][i] = v;
            }
        }
        m
    }
}

impl SymmetricMatrix {
    /// Compute eigenvalues using the Jacobi eigenvalue algorithm.
    /// Only suitable for small matrices (n <= ~50).
    /// Returns eigenvalues sorted in ascending order.
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n;
        if n == 0 {
            return vec![];
        }
        // Work with full dense matrix
        let mut m = self.to_full();

        let max_sweeps = 100;
        for _sweep in 0..max_sweeps {
            // Find largest off-diagonal |m[p][q]|
            let mut max_val = 0.0f64;
            let mut p = 0;
            let mut q = 1;
            for i in 0..n {
                for j in (i + 1)..n {
                    if m[i][j].abs() > max_val {
                        max_val = m[i][j].abs();
                        p = i;
                        q = j;
                    }
                }
            }

            // Convergence check: off-diagonal is small relative to diagonal
            let diag_max = (0..n)
                .map(|i| m[i][i].abs())
                .fold(1e-300, f64::max);
            if max_val < 1e-12 * diag_max {
                break;
            }

            // Apply Jacobi rotation to zero m[p][q]
            let diff = m[q][q] - m[p][p];
            let t = if diff.abs() < 1e-20 * max_val {
                1.0
            } else {
                let phi = diff / (2.0 * m[p][q]);
                phi.signum() / (phi.abs() + (1.0 + phi * phi).sqrt())
            };
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = t * c;
            let tau = s / (1.0 + c);

            let apq = m[p][q];
            m[p][q] = 0.0;
            m[q][p] = 0.0;
            m[p][p] -= t * apq;
            m[q][q] += t * apq;

            for r in 0..n {
                if r != p && r != q {
                    let rp = m[r][p];
                    let rq = m[r][q];
                    m[r][p] = rp - s * (rq + tau * rp);
                    m[p][r] = m[r][p];
                    m[r][q] = rq + s * (rp - tau * rq);
                    m[q][r] = m[r][q];
                }
            }
        }

        let mut eigs: Vec<f64> = (0..n).map(|i| m[i][i]).collect();
        eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        eigs
    }
}

/// Trait for linear solvers used within the IPM.
pub trait LinearSolver {
    /// Factor the symmetric matrix. Returns inertia if the solver can compute it.
    fn factor(&mut self, matrix: &SymmetricMatrix) -> Result<Option<Inertia>, SolverError>;

    /// Solve the system using the most recent factorization.
    /// Reads rhs and writes the solution into `solution`.
    fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError>;

    /// Whether this solver can report inertia.
    fn provides_inertia(&self) -> bool;
}
