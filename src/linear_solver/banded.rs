use super::{Inertia, KktMatrix, LinearSolver, SolverError};

/// Banded LDL^T solver for symmetric indefinite matrices with small bandwidth.
///
/// For a matrix of dimension n with half-bandwidth b, factorization is O(n*b²)
/// and solve is O(n*b). This is dramatically faster than general sparse solvers
/// for problems with banded structure (e.g., PDE discretizations).
///
/// Storage: lower triangle in column-oriented band format.
/// `band[k * n + j]` = A[j+k, j] for k = 0..=b, j = 0..n.
pub struct BandedLdl {
    n: usize,
    b: usize,
    band: Vec<f64>,
    d: Vec<f64>,
    factored: bool,
    zero_pivot_tol: f64,
}

impl BandedLdl {
    pub fn new() -> Self {
        Self {
            n: 0,
            b: 0,
            band: Vec::new(),
            d: Vec::new(),
            factored: false,
            zero_pivot_tol: 1e-20,
        }
    }

    /// Compute half-bandwidth from upper-triangle COO triplets.
    /// Returns max(col - row) over all entries.
    pub fn compute_bandwidth(rows: &[usize], cols: &[usize]) -> usize {
        let mut bw = 0usize;
        for (&r, &c) in rows.iter().zip(cols.iter()) {
            // Upper triangle: r <= c
            let diff = if c >= r { c - r } else { r - c };
            if diff > bw {
                bw = diff;
            }
        }
        bw
    }

    /// Scatter upper-triangle COO triplets into lower-triangle band storage.
    /// COO entry (r, c) with r <= c maps to lower entry A[c, r], stored at band[(c-r)*n + r].
    fn scatter_coo_to_band(&mut self, rows: &[usize], cols: &[usize], vals: &[f64]) {
        // Zero the band
        for v in self.band.iter_mut() {
            *v = 0.0;
        }
        let n = self.n;
        for i in 0..rows.len() {
            let r = rows[i];
            let c = cols[i];
            let v = vals[i];
            // Upper triangle: r <= c. Lower triangle equivalent: A[c, r]
            // band index: (c - r) * n + r
            let k = if c >= r { c - r } else { r - c };
            let j = r.min(c);
            debug_assert!(k <= self.b, "entry ({},{}) exceeds bandwidth {}", r, c, self.b);
            self.band[k * n + j] += v;
        }
    }

    /// Banded LDL^T factorization without pivoting.
    /// L overwrites the sub-diagonal entries of band (band[k*n+j] for k > 0).
    /// D is stored in self.d.
    fn factor_banded_ldlt(&mut self) {
        let n = self.n;
        let b = self.b;

        for j in 0..n {
            // d[j] = band[0][j] (current diagonal after updates)
            self.d[j] = self.band[j]; // band[0 * n + j]

            if self.d[j].abs() < self.zero_pivot_tol {
                // Near-zero pivot — leave L entries as zero, inertia will flag this
                for k in 1..=b.min(n - 1 - j) {
                    self.band[k * n + j] = 0.0;
                }
                continue;
            }

            let dj = self.d[j];

            // Compute L entries: L[j+k, j] = A[j+k, j] / d[j]
            for k in 1..=b.min(n - 1 - j) {
                self.band[k * n + j] /= dj;
            }

            // Rank-1 update of trailing submatrix
            for k in 1..=b.min(n - 1 - j) {
                let l_jk = self.band[k * n + j]; // L[j+k, j]
                let factor = dj * l_jk;
                for i in k..=b.min(n - 1 - j) {
                    // A[j+i, j+k] -= L[j+i, j] * d[j] * L[j+k, j]
                    // band[(i-k), j+k] -= band[i, j] * factor
                    self.band[(i - k) * n + (j + k)] -= self.band[i * n + j] * factor;
                }
            }
        }
    }

    fn compute_inertia(&self) -> Inertia {
        let mut pos = 0;
        let mut neg = 0;
        let mut zero = 0;
        for j in 0..self.n {
            if self.d[j].abs() < self.zero_pivot_tol {
                zero += 1;
            } else if self.d[j] > 0.0 {
                pos += 1;
            } else {
                neg += 1;
            }
        }
        Inertia {
            positive: pos,
            negative: neg,
            zero,
        }
    }

    /// Forward substitution: solve L * y = rhs (L is unit lower triangular banded).
    /// Then diagonal solve: y /= d.
    /// Then backward substitution: solve L^T * x = y.
    fn solve_banded_ldlt(&self, rhs: &[f64], solution: &mut [f64]) {
        let n = self.n;
        let b = self.b;

        // Forward: L * y = rhs
        solution[..n].copy_from_slice(&rhs[..n]);
        for j in 0..n {
            let yj = solution[j];
            for k in 1..=b.min(n - 1 - j) {
                solution[j + k] -= self.band[k * n + j] * yj;
            }
        }

        // Diagonal: y /= d
        for j in 0..n {
            if self.d[j].abs() < self.zero_pivot_tol {
                solution[j] = 0.0;
            } else {
                solution[j] /= self.d[j];
            }
        }

        // Backward: L^T * x = y
        for j in (0..n).rev() {
            for k in 1..=b.min(n - 1 - j) {
                solution[j] -= self.band[k * n + j] * solution[j + k];
            }
        }
    }
}

impl LinearSolver for BandedLdl {
    fn factor(&mut self, matrix: &KktMatrix) -> Result<Option<Inertia>, SolverError> {
        let (rows, cols, vals, n) = match matrix {
            KktMatrix::Sparse(s) => (
                &s.triplet_rows,
                &s.triplet_cols,
                &s.triplet_vals,
                s.n,
            ),
            KktMatrix::Dense(_) => {
                return Err(SolverError::NumericalFailure(
                    "BandedLdl requires KktMatrix::Sparse".into(),
                ));
            }
        };

        if self.n != n || self.band.is_empty() {
            self.n = n;
            self.b = Self::compute_bandwidth(rows, cols);
            self.band = vec![0.0; (self.b + 1) * n];
            self.d = vec![0.0; n];
        }

        self.scatter_coo_to_band(rows, cols, vals);
        self.factor_banded_ldlt();
        self.factored = true;

        Ok(Some(self.compute_inertia()))
    }

    fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        if !self.factored {
            return Err(SolverError::NumericalFailure(
                "BandedLdl: must call factor() before solve()".into(),
            ));
        }
        self.solve_banded_ldlt(rhs, solution);
        Ok(())
    }

    fn provides_inertia(&self) -> bool {
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        if !self.factored || self.n == 0 {
            return None;
        }
        Some(self.d.iter().copied().reduce(f64::min).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::SparseSymmetricMatrix;

    /// Build a SparseSymmetricMatrix from dense lower triangle entries.
    fn sparse_from_entries(n: usize, entries: &[(usize, usize, f64)]) -> SparseSymmetricMatrix {
        let mut m = SparseSymmetricMatrix::with_capacity(n, entries.len() + n);
        for &(r, c, v) in entries {
            m.add(r, c, v);
        }
        m
    }

    #[test]
    fn test_tridiagonal_positive_definite() {
        // 1D Laplacian: [2 -1 0; -1 2 -1; 0 -1 2]
        let n = 4;
        let mat = sparse_from_entries(n, &[
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 1, 2.0), (1, 2, -1.0),
            (2, 2, 2.0), (2, 3, -1.0),
            (3, 3, 2.0),
        ]);
        let kkt = KktMatrix::Sparse(mat);

        let mut solver = BandedLdl::new();
        let inertia = solver.factor(&kkt).unwrap().unwrap();
        assert_eq!(inertia.positive, 4);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);

        // Solve A*x = [1, 0, 0, 1]
        let rhs = [1.0, 0.0, 0.0, 1.0];
        let mut sol = vec![0.0; 4];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify A*x ≈ rhs
        let mut ax = vec![0.0; 4];
        ax[0] = 2.0 * sol[0] - sol[1];
        ax[1] = -sol[0] + 2.0 * sol[1] - sol[2];
        ax[2] = -sol[1] + 2.0 * sol[2] - sol[3];
        ax[3] = -sol[2] + 2.0 * sol[3];
        for i in 0..4 {
            assert!((ax[i] - rhs[i]).abs() < 1e-12, "ax[{}]={}, rhs[{}]={}", i, ax[i], i, rhs[i]);
        }
    }

    #[test]
    fn test_kkt_indefinite() {
        // KKT system: [H J^T; J 0] with H=diag(2,2), J=[1 1]
        // dim=3, bandwidth=2
        // Upper triangle: (0,0)=2, (1,1)=2, (0,2)=1, (1,2)=1, (2,2)=0
        let mat = sparse_from_entries(3, &[
            (0, 0, 2.0), (1, 1, 2.0),
            (0, 2, 1.0), (1, 2, 1.0),
        ]);
        let kkt = KktMatrix::Sparse(mat);

        let mut solver = BandedLdl::new();
        let inertia = solver.factor(&kkt).unwrap().unwrap();
        assert_eq!(inertia.positive, 2, "expected 2 positive, got {}", inertia.positive);
        assert_eq!(inertia.negative, 1, "expected 1 negative, got {}", inertia.negative);
        assert_eq!(inertia.zero, 0);

        // Solve
        let rhs = [1.0, 1.0, 0.0];
        let mut sol = vec![0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify A*x ≈ rhs
        let ax0 = 2.0 * sol[0] + sol[2];
        let ax1 = 2.0 * sol[1] + sol[2];
        let ax2 = sol[0] + sol[1];
        assert!((ax0 - rhs[0]).abs() < 1e-12);
        assert!((ax1 - rhs[1]).abs() < 1e-12);
        assert!((ax2 - rhs[2]).abs() < 1e-12);
    }

    #[test]
    fn test_diagonal_matrix() {
        let mat = sparse_from_entries(3, &[
            (0, 0, 4.0), (1, 1, -2.0), (2, 2, 3.0),
        ]);
        let kkt = KktMatrix::Sparse(mat);

        let mut solver = BandedLdl::new();
        let inertia = solver.factor(&kkt).unwrap().unwrap();
        assert_eq!(inertia.positive, 2);
        assert_eq!(inertia.negative, 1);

        let rhs = [8.0, -6.0, 9.0];
        let mut sol = vec![0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();
        assert!((sol[0] - 2.0).abs() < 1e-12);
        assert!((sol[1] - 3.0).abs() < 1e-12);
        assert!((sol[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_bandwidth_detection() {
        let rows = vec![0, 0, 1, 1, 2, 3];
        let cols = vec![0, 1, 1, 2, 2, 3];
        assert_eq!(BandedLdl::compute_bandwidth(&rows, &cols), 1);

        // Pentadiagonal
        let rows2 = vec![0, 0, 0, 1, 2];
        let cols2 = vec![0, 1, 2, 1, 2];
        assert_eq!(BandedLdl::compute_bandwidth(&rows2, &cols2), 2);
    }

    #[test]
    fn test_large_tridiagonal() {
        // Large tridiagonal: verify O(n) performance
        let n = 1000;
        let mut entries = Vec::with_capacity(2 * n);
        for i in 0..n {
            entries.push((i, i, 4.0));
            if i + 1 < n {
                entries.push((i, i + 1, -1.0));
            }
        }
        let mat = sparse_from_entries(n, &entries);
        let kkt = KktMatrix::Sparse(mat);

        let mut solver = BandedLdl::new();
        let inertia = solver.factor(&kkt).unwrap().unwrap();
        assert_eq!(inertia.positive, n);

        let rhs: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut sol = vec![0.0; n];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify residual
        for i in 0..n {
            let mut ax = 4.0 * sol[i];
            if i > 0 { ax -= sol[i - 1]; }
            if i + 1 < n { ax -= sol[i + 1]; }
            assert!((ax - rhs[i]).abs() < 1e-8, "residual at {}: {}", i, (ax - rhs[i]).abs());
        }
    }
}
