use crate::SolverError;

/// Sparse symmetric matrix in COO (triplet) format, upper triangle (row <= col).
/// Duplicate entries at the same position are summed during CSC conversion.
#[derive(Debug, Clone)]
pub struct CooMatrix {
    /// Matrix dimension (n x n).
    pub n: usize,
    /// Row indices (row <= col).
    pub rows: Vec<usize>,
    /// Column indices (col >= row).
    pub cols: Vec<usize>,
    /// Values.
    pub vals: Vec<f64>,
}

impl CooMatrix {
    /// Create a new COO matrix with validation.
    /// All entries must satisfy row <= col and indices < n.
    pub fn new(
        n: usize,
        rows: Vec<usize>,
        cols: Vec<usize>,
        vals: Vec<f64>,
    ) -> Result<Self, SolverError> {
        if rows.len() != cols.len() || rows.len() != vals.len() {
            return Err(SolverError::InvalidInput(format!(
                "rows ({}), cols ({}), and vals ({}) must have the same length",
                rows.len(),
                cols.len(),
                vals.len()
            )));
        }
        for k in 0..rows.len() {
            if rows[k] > cols[k] {
                return Err(SolverError::InvalidInput(format!(
                    "entry {} has row ({}) > col ({}); upper triangle required",
                    k, rows[k], cols[k]
                )));
            }
            if cols[k] >= n {
                return Err(SolverError::InvalidInput(format!(
                    "entry {} has index ({}, {}) out of bounds for n={}",
                    k, rows[k], cols[k], n
                )));
            }
        }
        Ok(Self { n, rows, cols, vals })
    }

    /// Create an empty COO matrix.
    pub fn empty(n: usize) -> Self {
        Self {
            n,
            rows: Vec::new(),
            cols: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// Number of stored entries (may include duplicates).
    pub fn nnz(&self) -> usize {
        self.rows.len()
    }

    /// Symmetric matrix-vector product y = A * x.
    /// Treats the upper-triangle COO as a full symmetric matrix.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        if x.len() != self.n || y.len() != self.n {
            return Err(SolverError::DimensionMismatch {
                expected: self.n,
                got: x.len().min(y.len()),
            });
        }
        for yi in y.iter_mut() {
            *yi = 0.0;
        }
        for k in 0..self.nnz() {
            let i = self.rows[k];
            let j = self.cols[k];
            let v = self.vals[k];
            y[i] += v * x[j];
            if i != j {
                y[j] += v * x[i];
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_matrix() {
        let m = CooMatrix::empty(5);
        assert_eq!(m.n, 5);
        assert_eq!(m.nnz(), 0);
    }

    #[test]
    fn test_valid_construction() {
        // 3x3 identity
        let m = CooMatrix::new(
            3,
            vec![0, 1, 2],
            vec![0, 1, 2],
            vec![1.0, 1.0, 1.0],
        );
        assert!(m.is_ok());
        assert_eq!(m.unwrap().nnz(), 3);
    }

    #[test]
    fn test_reject_lower_triangle() {
        // row=1, col=0 violates row <= col
        let m = CooMatrix::new(3, vec![1], vec![0], vec![1.0]);
        assert!(m.is_err());
    }

    #[test]
    fn test_reject_out_of_bounds() {
        let m = CooMatrix::new(3, vec![0], vec![3], vec![1.0]);
        assert!(m.is_err());
    }

    #[test]
    fn test_reject_mismatched_lengths() {
        let m = CooMatrix::new(3, vec![0, 1], vec![0], vec![1.0]);
        assert!(m.is_err());
    }

    #[test]
    fn test_matvec_identity() {
        let m = CooMatrix::new(3, vec![0, 1, 2], vec![0, 1, 2], vec![1.0, 1.0, 1.0]).unwrap();
        let x = [2.0, 3.0, 4.0];
        let mut y = [0.0; 3];
        m.matvec(&x, &mut y).unwrap();
        assert!((y[0] - 2.0).abs() < 1e-15);
        assert!((y[1] - 3.0).abs() < 1e-15);
        assert!((y[2] - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_matvec_symmetric() {
        // A = [[2, 1, 0], [1, 3, 1], [0, 1, 4]]
        // Upper triangle: (0,0)=2, (0,1)=1, (1,1)=3, (1,2)=1, (2,2)=4
        let m = CooMatrix::new(
            3,
            vec![0, 0, 1, 1, 2],
            vec![0, 1, 1, 2, 2],
            vec![2.0, 1.0, 3.0, 1.0, 4.0],
        )
        .unwrap();
        let x = [1.0, 2.0, 3.0];
        let mut y = [0.0; 3];
        m.matvec(&x, &mut y).unwrap();
        // y = [2*1+1*2, 1*1+3*2+1*3, 1*2+4*3] = [4, 10, 14]
        assert!((y[0] - 4.0).abs() < 1e-15);
        assert!((y[1] - 10.0).abs() < 1e-15);
        assert!((y[2] - 14.0).abs() < 1e-15);
    }

    #[test]
    fn test_matvec_with_duplicates() {
        // Two entries at (0,0): 1.0 + 2.0 = 3.0
        let m = CooMatrix::new(2, vec![0, 0, 1], vec![0, 0, 1], vec![1.0, 2.0, 4.0]).unwrap();
        let x = [1.0, 1.0];
        let mut y = [0.0; 2];
        m.matvec(&x, &mut y).unwrap();
        assert!((y[0] - 3.0).abs() < 1e-15);
        assert!((y[1] - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_matvec_dimension_mismatch() {
        let m = CooMatrix::new(3, vec![0], vec![0], vec![1.0]).unwrap();
        let x = [1.0, 2.0]; // wrong size
        let mut y = [0.0; 3];
        assert!(m.matvec(&x, &mut y).is_err());
    }
}
