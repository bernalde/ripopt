use crate::coo::CooMatrix;

/// Sparse symmetric matrix in CSC (compressed sparse column) format.
/// Stores upper triangle only (row <= col for each entry).
#[derive(Debug, Clone)]
pub struct CscMatrix {
    /// Matrix dimension (n x n).
    pub n: usize,
    /// Column pointers (length n+1). Column j has entries `col_ptr[j]..col_ptr[j+1]`.
    pub col_ptr: Vec<usize>,
    /// Row indices (sorted within each column).
    pub row_idx: Vec<usize>,
    /// Values corresponding to row_idx.
    pub vals: Vec<f64>,
}

impl CscMatrix {
    /// Convert a COO matrix to CSC, summing duplicate entries.
    pub fn from_coo(coo: &CooMatrix) -> Self {
        let n = coo.n;
        if coo.nnz() == 0 {
            return Self {
                n,
                col_ptr: vec![0; n + 1],
                row_idx: Vec::new(),
                vals: Vec::new(),
            };
        }

        // Count entries per column
        let mut col_count = vec![0usize; n];
        for &c in &coo.cols {
            col_count[c] += 1;
        }

        // Build column pointers
        let mut col_ptr = vec![0usize; n + 1];
        for j in 0..n {
            col_ptr[j + 1] = col_ptr[j] + col_count[j];
        }
        let nnz = col_ptr[n];

        // Place entries into CSC arrays
        let mut row_idx = vec![0usize; nnz];
        let mut vals = vec![0.0f64; nnz];
        let mut offset = col_ptr.clone(); // workspace: current write position per column

        for k in 0..coo.nnz() {
            let j = coo.cols[k];
            let pos = offset[j];
            row_idx[pos] = coo.rows[k];
            vals[pos] = coo.vals[k];
            offset[j] += 1;
        }

        // Sort each column by row index and sum duplicates
        for j in 0..n {
            let start = col_ptr[j];
            let end = col_ptr[j + 1];
            if start == end {
                continue;
            }

            // Sort by row index (insertion sort — columns are typically small)
            for i in (start + 1)..end {
                let key_row = row_idx[i];
                let key_val = vals[i];
                let mut pos = i;
                while pos > start && row_idx[pos - 1] > key_row {
                    row_idx[pos] = row_idx[pos - 1];
                    vals[pos] = vals[pos - 1];
                    pos -= 1;
                }
                row_idx[pos] = key_row;
                vals[pos] = key_val;
            }
        }

        // Sum duplicates in-place
        let mut new_row_idx = Vec::with_capacity(nnz);
        let mut new_vals = Vec::with_capacity(nnz);
        let mut new_col_ptr = vec![0usize; n + 1];

        for j in 0..n {
            let start = col_ptr[j];
            let end = col_ptr[j + 1];
            let col_start = new_row_idx.len();

            let mut i = start;
            while i < end {
                let r = row_idx[i];
                let mut v = vals[i];
                i += 1;
                while i < end && row_idx[i] == r {
                    v += vals[i];
                    i += 1;
                }
                new_row_idx.push(r);
                new_vals.push(v);
            }

            new_col_ptr[j + 1] = new_row_idx.len();
            debug_assert!(new_col_ptr[j] == col_start);
        }

        Self {
            n,
            col_ptr: new_col_ptr,
            row_idx: new_row_idx,
            vals: new_vals,
        }
    }

    /// Number of stored entries (after duplicate summing).
    pub fn nnz(&self) -> usize {
        self.col_ptr[self.n]
    }

    /// Symmetric matrix-vector product y = A * x.
    /// Treats upper-triangle CSC as full symmetric matrix.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        for yi in y.iter_mut() {
            *yi = 0.0;
        }
        for j in 0..self.n {
            for idx in self.col_ptr[j]..self.col_ptr[j + 1] {
                let i = self.row_idx[idx];
                let v = self.vals[idx];
                y[i] += v * x[j];
                if i != j {
                    y[j] += v * x[i];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_csc() {
        let coo = CooMatrix::empty(3);
        let csc = CscMatrix::from_coo(&coo);
        assert_eq!(csc.n, 3);
        assert_eq!(csc.nnz(), 0);
        assert_eq!(csc.col_ptr, vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_identity_csc() {
        let coo = CooMatrix::new(3, vec![0, 1, 2], vec![0, 1, 2], vec![1.0, 1.0, 1.0]).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        assert_eq!(csc.n, 3);
        assert_eq!(csc.nnz(), 3);
        assert_eq!(csc.col_ptr, vec![0, 1, 2, 3]);
        assert_eq!(csc.row_idx, vec![0, 1, 2]);
    }

    #[test]
    fn test_duplicate_summing() {
        // Two entries at (0,0): 1.0 + 2.0 = 3.0, one at (1,1): 4.0
        let coo = CooMatrix::new(2, vec![0, 0, 1], vec![0, 0, 1], vec![1.0, 2.0, 4.0]).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        assert_eq!(csc.nnz(), 2); // duplicates merged
        assert_eq!(csc.col_ptr, vec![0, 1, 2]);
        assert!((csc.vals[0] - 3.0).abs() < 1e-15); // summed
        assert!((csc.vals[1] - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_matvec_matches_coo() {
        // A = [[2, 1, 0], [1, 3, 1], [0, 1, 4]]
        let coo = CooMatrix::new(
            3,
            vec![0, 0, 1, 1, 2],
            vec![0, 1, 1, 2, 2],
            vec![2.0, 1.0, 3.0, 1.0, 4.0],
        )
        .unwrap();

        let x = [1.0, 2.0, 3.0];
        let mut y_coo = [0.0; 3];
        coo.matvec(&x, &mut y_coo).unwrap();

        let csc = CscMatrix::from_coo(&coo);
        let mut y_csc = [0.0; 3];
        csc.matvec(&x, &mut y_csc);

        for i in 0..3 {
            assert!(
                (y_coo[i] - y_csc[i]).abs() < 1e-15,
                "mismatch at {}: coo={} csc={}",
                i,
                y_coo[i],
                y_csc[i]
            );
        }
    }

    #[test]
    fn test_csc_sorted_rows() {
        // Insert entries out of row order within a column
        // Column 2: rows 2, 0 (out of order)
        let coo = CooMatrix::new(
            3,
            vec![2, 0, 1],
            vec![2, 2, 1],
            vec![4.0, 1.0, 3.0],
        )
        .unwrap();
        let csc = CscMatrix::from_coo(&coo);
        // Column 2 should have rows sorted: [0, 2]
        let start = csc.col_ptr[2];
        let end = csc.col_ptr[3];
        assert_eq!(end - start, 2);
        assert_eq!(csc.row_idx[start], 0);
        assert_eq!(csc.row_idx[start + 1], 2);
        assert!((csc.vals[start] - 1.0).abs() < 1e-15);
        assert!((csc.vals[start + 1] - 4.0).abs() < 1e-15);
    }
}
