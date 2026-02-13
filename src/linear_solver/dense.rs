use super::{Inertia, LinearSolver, SolverError, SymmetricMatrix};

/// Dense LDL^T factorization with Bunch-Kaufman pivoting for symmetric indefinite matrices.
///
/// Factors A = P * L * D * L^T * P^T where:
/// - P is a permutation matrix
/// - L is unit lower triangular
/// - D is block diagonal (1x1 and 2x2 blocks)
///
/// This handles the indefinite KKT matrices that arise in interior point methods.
pub struct DenseLdl {
    /// Dimension of the factored matrix.
    n: usize,
    /// The L factor stored as full n×n (lower triangular with unit diagonal).
    l: Vec<f64>,
    /// The D factor: for 1x1 blocks, d[i] is the diagonal entry.
    /// For 2x2 blocks, we store in d and d_offdiag.
    d: Vec<f64>,
    /// Off-diagonal of 2x2 blocks in D. d_offdiag[i] != 0 means (i, i+1) is a 2x2 block.
    d_offdiag: Vec<f64>,
    /// Pivot permutation.
    perm: Vec<usize>,
    /// Inverse permutation.
    perm_inv: Vec<usize>,
    /// Whether the matrix has been factored.
    factored: bool,
    /// Tolerance for determining zero pivots.
    zero_pivot_tol: f64,
}

impl Default for DenseLdl {
    fn default() -> Self {
        Self::new()
    }
}

impl DenseLdl {
    pub fn new() -> Self {
        Self {
            n: 0,
            l: Vec::new(),
            d: Vec::new(),
            d_offdiag: Vec::new(),
            perm: Vec::new(),
            perm_inv: Vec::new(),
            factored: false,
            zero_pivot_tol: 1e-12,
        }
    }

    /// Access L[i][j] (full storage, row-major).
    #[inline]
    fn l_idx(n: usize, i: usize, j: usize) -> usize {
        i * n + j
    }

    /// Factor the matrix using Bunch-Kaufman pivoting.
    ///
    /// This implements the symmetric indefinite factorization A = P L D L^T P^T.
    /// We use the standard Bunch-Kaufman algorithm with alpha = (1 + sqrt(17)) / 8.
    fn bunch_kaufman_factor(&mut self, matrix: &SymmetricMatrix) -> Result<Inertia, SolverError> {
        let n = matrix.n;
        self.n = n;
        self.l = vec![0.0; n * n];
        self.d = vec![0.0; n];
        self.d_offdiag = vec![0.0; n];
        self.perm = (0..n).collect();
        self.perm_inv = (0..n).collect();

        if n == 0 {
            self.factored = true;
            return Ok(Inertia {
                positive: 0,
                negative: 0,
                zero: 0,
            });
        }

        // Work with a full dense copy (column-major) for easier pivoting
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = matrix.get(i, j);
                a[i * n + j] = v;
                a[j * n + i] = v;
            }
        }

        // Bunch-Kaufman alpha parameter
        let alpha = (1.0 + 17.0_f64.sqrt()) / 8.0;

        let mut k = 0;
        while k < n {
            // Find pivot
            let (pivot_type, p1, p2) =
                self.find_pivot(&a, n, k, alpha);

            if pivot_type == 1 {
                // 1x1 pivot
                if p1 != k {
                    // Swap rows/columns k and p1
                    self.swap_rows_cols(&mut a, n, k, p1);
                    self.perm.swap(k, p1);
                    // Also swap L entries for previously computed columns
                    for j in 0..k {
                        let idx_k = Self::l_idx(n, k, j);
                        let idx_p1 = Self::l_idx(n, p1, j);
                        self.l.swap(idx_k, idx_p1);
                    }
                }

                let akk = a[k * n + k];
                self.d[k] = akk;

                if akk.abs() > self.zero_pivot_tol {
                    // Compute L column k
                    for i in (k + 1)..n {
                        let idx = Self::l_idx(n, i, k);
                        self.l[idx] = a[i * n + k] / akk;
                    }

                    // Update remaining submatrix: A -= L[:,k] * D[k] * L[:,k]^T
                    for j in (k + 1)..n {
                        for i in j..n {
                            let lik = self.l[Self::l_idx(n, i, k)];
                            let ljk = self.l[Self::l_idx(n, j, k)];
                            a[i * n + j] -= lik * akk * ljk;
                            a[j * n + i] = a[i * n + j]; // symmetric
                        }
                    }
                }
                // Set diagonal of L to 1
                self.l[Self::l_idx(n, k, k)] = 1.0;
                k += 1;
            } else {
                // 2x2 pivot using rows/columns k and p2
                // First bring p2 to position k+1
                if p2 != k + 1 {
                    self.swap_rows_cols(&mut a, n, k + 1, p2);
                    self.perm.swap(k + 1, p2);
                    // Also swap L entries for previously computed columns
                    for j in 0..k {
                        let idx_k1 = Self::l_idx(n, k + 1, j);
                        let idx_p2 = Self::l_idx(n, p2, j);
                        self.l.swap(idx_k1, idx_p2);
                    }
                }
                // Now if p1 != k, swap
                if p1 != k {
                    self.swap_rows_cols(&mut a, n, k, p1);
                    self.perm.swap(k, p1);
                    // Also swap L entries for previously computed columns
                    for j in 0..k {
                        let idx_k = Self::l_idx(n, k, j);
                        let idx_p1 = Self::l_idx(n, p1, j);
                        self.l.swap(idx_k, idx_p1);
                    }
                }

                let akk = a[k * n + k];
                let ak1k = a[(k + 1) * n + k];
                let ak1k1 = a[(k + 1) * n + (k + 1)];

                self.d[k] = akk;
                self.d[k + 1] = ak1k1;
                self.d_offdiag[k] = ak1k;

                // D is the 2x2 block [[akk, ak1k], [ak1k, ak1k1]]
                let det = akk * ak1k1 - ak1k * ak1k;

                if det.abs() > self.zero_pivot_tol {
                    // D^{-1} = (1/det) * [[ak1k1, -ak1k], [-ak1k, akk]]
                    let d_inv_00 = ak1k1 / det;
                    let d_inv_01 = -ak1k / det;
                    let d_inv_11 = akk / det;

                    // Compute L columns k and k+1
                    for i in (k + 2)..n {
                        let aik = a[i * n + k];
                        let aik1 = a[i * n + (k + 1)];
                        self.l[Self::l_idx(n, i, k)] = aik * d_inv_00 + aik1 * d_inv_01;
                        self.l[Self::l_idx(n, i, k + 1)] = aik * d_inv_01 + aik1 * d_inv_11;
                    }

                    // Update remaining submatrix
                    for j in (k + 2)..n {
                        for i in j..n {
                            let lik = self.l[Self::l_idx(n, i, k)];
                            let lik1 = self.l[Self::l_idx(n, i, k + 1)];
                            let ljk = self.l[Self::l_idx(n, j, k)];
                            let ljk1 = self.l[Self::l_idx(n, j, k + 1)];

                            let update = lik * (akk * ljk + ak1k * ljk1)
                                + lik1 * (ak1k * ljk + ak1k1 * ljk1);
                            a[i * n + j] -= update;
                            a[j * n + i] = a[i * n + j];
                        }
                    }
                }

                self.l[Self::l_idx(n, k, k)] = 1.0;
                self.l[Self::l_idx(n, k + 1, k + 1)] = 1.0;
                k += 2;
            }
        }

        // Compute inverse permutation
        for i in 0..n {
            self.perm_inv[self.perm[i]] = i;
        }

        self.factored = true;

        Ok(self.compute_inertia())
    }

    /// Find the pivot for Bunch-Kaufman.
    /// Returns (pivot_type, p1, p2) where:
    /// - pivot_type = 1 for 1x1 pivot at position p1
    /// - pivot_type = 2 for 2x2 pivot at positions (p1, p2)
    fn find_pivot(&self, a: &[f64], n: usize, k: usize, alpha: f64) -> (usize, usize, usize) {
        if k == n - 1 {
            return (1, k, k);
        }

        let akk = a[k * n + k].abs();

        // Find largest off-diagonal in column k (below diagonal)
        let mut lambda = 0.0f64;
        let mut r = k;
        for i in (k + 1)..n {
            let v = a[i * n + k].abs();
            if v > lambda {
                lambda = v;
                r = i;
            }
        }

        if lambda == 0.0 && akk == 0.0 {
            // Zero column — take 1x1 pivot (will be zero pivot)
            return (1, k, k);
        }

        if akk >= alpha * lambda {
            // 1x1 pivot is good enough
            return (1, k, k);
        }

        // Find largest off-diagonal in row r (the row with largest off-diag in col k)
        let mut sigma = 0.0f64;
        for j in k..n {
            if j != r {
                let v = a[r * n + j].abs();
                if v > sigma {
                    sigma = v;
                }
            }
        }

        if akk * sigma >= alpha * lambda * lambda {
            // 1x1 pivot at k
            return (1, k, k);
        }

        let arr = a[r * n + r].abs();
        if arr >= alpha * sigma {
            // 1x1 pivot at r
            return (1, r, r);
        }

        // 2x2 pivot using k and r
        (2, k, r)
    }

    /// Swap rows and columns p and q in a full symmetric matrix.
    fn swap_rows_cols(&self, a: &mut [f64], n: usize, p: usize, q: usize) {
        if p == q {
            return;
        }
        // Swap rows p and q
        for j in 0..n {
            let idx_p = p * n + j;
            let idx_q = q * n + j;
            a.swap(idx_p, idx_q);
        }
        // Swap columns p and q
        for i in 0..n {
            let idx_p = i * n + p;
            let idx_q = i * n + q;
            a.swap(idx_p, idx_q);
        }
    }

    /// Compute inertia from the D factor.
    fn compute_inertia(&self) -> Inertia {
        let mut positive = 0;
        let mut negative = 0;
        let mut zero = 0;

        let mut k = 0;
        while k < self.n {
            if k + 1 < self.n && self.d_offdiag[k].abs() > self.zero_pivot_tol {
                // 2x2 block: eigenvalues from [[d[k], d_offdiag[k]], [d_offdiag[k], d[k+1]]]
                let a = self.d[k];
                let b = self.d_offdiag[k];
                let c = self.d[k + 1];

                let trace = a + c;
                let det = a * c - b * b;
                let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();

                let eig1 = (trace + disc) / 2.0;
                let eig2 = (trace - disc) / 2.0;

                for eig in [eig1, eig2] {
                    if eig > self.zero_pivot_tol {
                        positive += 1;
                    } else if eig < -self.zero_pivot_tol {
                        negative += 1;
                    } else {
                        zero += 1;
                    }
                }
                k += 2;
            } else {
                // 1x1 block
                let d = self.d[k];
                if d > self.zero_pivot_tol {
                    positive += 1;
                } else if d < -self.zero_pivot_tol {
                    negative += 1;
                } else {
                    zero += 1;
                }
                k += 1;
            }
        }

        Inertia {
            positive,
            negative,
            zero,
        }
    }

    /// Solve L * x = b (forward substitution with permutation and 2x2 blocks in D).
    fn solve_internal(&self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        if !self.factored {
            return Err(SolverError::NumericalFailure(
                "matrix not factored".to_string(),
            ));
        }

        let n = self.n;
        if rhs.len() != n || solution.len() != n {
            return Err(SolverError::DimensionMismatch {
                expected: n,
                got: rhs.len(),
            });
        }

        // Step 1: Apply permutation: y = P * rhs
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] = rhs[self.perm[i]];
        }

        // Step 2: Forward substitution: L * z = y
        for i in 0..n {
            let mut sum = y[i];
            #[allow(clippy::needless_range_loop)]
            for j in 0..i {
                sum -= self.l[Self::l_idx(n, i, j)] * y[j];
            }
            y[i] = sum;
        }

        // Step 3: Solve D * w = z
        let mut w = vec![0.0; n];
        let mut k = 0;
        while k < n {
            if k + 1 < n && self.d_offdiag[k].abs() > self.zero_pivot_tol {
                // 2x2 block
                let a = self.d[k];
                let b = self.d_offdiag[k];
                let c = self.d[k + 1];
                let det = a * c - b * b;

                if det.abs() < self.zero_pivot_tol {
                    return Err(SolverError::SingularMatrix);
                }

                w[k] = (c * y[k] - b * y[k + 1]) / det;
                w[k + 1] = (a * y[k + 1] - b * y[k]) / det;
                k += 2;
            } else {
                // 1x1 block
                if self.d[k].abs() < self.zero_pivot_tol {
                    return Err(SolverError::SingularMatrix);
                }
                w[k] = y[k] / self.d[k];
                k += 1;
            }
        }

        // Step 4: Backward substitution: L^T * v = w
        for i in (0..n).rev() {
            let mut sum = w[i];
            #[allow(clippy::needless_range_loop)]
            for j in (i + 1)..n {
                sum -= self.l[Self::l_idx(n, j, i)] * w[j];
            }
            w[i] = sum;
        }

        // Step 5: Apply inverse permutation: solution = P^T * v
        for i in 0..n {
            solution[self.perm[i]] = w[i];
        }

        Ok(())
    }
}

impl DenseLdl {
    /// Return the minimum diagonal entry of D after factorization.
    /// For 2x2 blocks, returns the minimum eigenvalue of the block.
    /// Returns None if not factored or n=0.
    pub fn min_diagonal(&self) -> Option<f64> {
        if !self.factored || self.n == 0 {
            return None;
        }
        let mut min_d = f64::INFINITY;
        let mut k = 0;
        while k < self.n {
            if k + 1 < self.n && self.d_offdiag[k].abs() > self.zero_pivot_tol {
                // 2x2 block: use eigenvalues
                let a = self.d[k];
                let b = self.d_offdiag[k];
                let c = self.d[k + 1];
                let trace = a + c;
                let det = a * c - b * b;
                let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
                let eig_min = (trace - disc) / 2.0;
                min_d = min_d.min(eig_min);
                k += 2;
            } else {
                min_d = min_d.min(self.d[k]);
                k += 1;
            }
        }
        Some(min_d)
    }
}

impl LinearSolver for DenseLdl {
    fn factor(&mut self, matrix: &SymmetricMatrix) -> Result<Option<Inertia>, SolverError> {
        let inertia = self.bunch_kaufman_factor(matrix)?;
        Ok(Some(inertia))
    }

    fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        self.solve_internal(rhs, solution)
    }

    fn provides_inertia(&self) -> bool {
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        DenseLdl::min_diagonal(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_positive_definite_3x3() {
        // A = [[4, 2, 1], [2, 5, 3], [1, 3, 6]]
        let mut mat = SymmetricMatrix::zeros(3);
        mat.set(0, 0, 4.0);
        mat.set(1, 0, 2.0);
        mat.set(1, 1, 5.0);
        mat.set(2, 0, 1.0);
        mat.set(2, 1, 3.0);
        mat.set(2, 2, 6.0);

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(inertia.positive, 3);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);

        // Solve Ax = b where b = [1, 2, 3]
        let rhs = [1.0, 2.0, 3.0];
        let mut sol = [0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify Ax = b
        let full = mat.to_full();
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i,
                ax_i,
                rhs[i]
            );
        }
    }

    #[test]
    fn test_indefinite_matrix() {
        // A = [[1, 2], [2, -1]] — eigenvalues: (1-1)/2 ± sqrt(4+1) = ±sqrt(5)
        // So one positive, one negative eigenvalue
        let mut mat = SymmetricMatrix::zeros(2);
        mat.set(0, 0, 1.0);
        mat.set(1, 0, 2.0);
        mat.set(1, 1, -1.0);

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(inertia.positive, 1);
        assert_eq!(inertia.negative, 1);
        assert_eq!(inertia.zero, 0);

        // Solve
        let rhs = [3.0, 1.0];
        let mut sol = [0.0; 2];
        solver.solve(&rhs, &mut sol).unwrap();

        let full = mat.to_full();
        for i in 0..2 {
            let ax_i: f64 = (0..2).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i,
                ax_i,
                rhs[i]
            );
        }
    }

    #[test]
    fn test_kkt_shaped_matrix() {
        // KKT matrix: [[H, J^T], [J, 0]]
        // H = [[2, 0], [0, 2]] (positive definite)
        // J = [[1, 1]] (1 constraint, 2 variables)
        // Full: [[2, 0, 1], [0, 2, 1], [1, 1, 0]]
        // Should have inertia (2, 1, 0)
        let mut mat = SymmetricMatrix::zeros(3);
        mat.set(0, 0, 2.0);
        mat.set(1, 1, 2.0);
        mat.set(2, 0, 1.0);
        mat.set(2, 1, 1.0);
        mat.set(2, 2, 0.0);

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(
            inertia.positive, 2,
            "Expected 2 positive, got {}",
            inertia.positive
        );
        assert_eq!(
            inertia.negative, 1,
            "Expected 1 negative, got {}",
            inertia.negative
        );
        assert_eq!(inertia.zero, 0, "Expected 0 zero, got {}", inertia.zero);

        // Solve the KKT system
        let rhs = [1.0, 2.0, 3.0];
        let mut sol = [0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        let full = mat.to_full();
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i,
                ax_i,
                rhs[i]
            );
        }
    }

    #[test]
    fn test_near_singular_matrix() {
        // A = [[1, 0], [0, 1e-15]] — should have one zero eigenvalue
        let mut mat = SymmetricMatrix::zeros(2);
        mat.set(0, 0, 1.0);
        mat.set(1, 1, 1e-15);

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(inertia.positive, 1);
        assert_eq!(inertia.zero, 1);
    }

    #[test]
    fn test_identity() {
        let mut mat = SymmetricMatrix::zeros(4);
        for i in 0..4 {
            mat.set(i, i, 1.0);
        }

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();
        assert_eq!(inertia.positive, 4);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);

        let rhs = [1.0, 2.0, 3.0, 4.0];
        let mut sol = [0.0; 4];
        solver.solve(&rhs, &mut sol).unwrap();

        for i in 0..4 {
            assert!((sol[i] - rhs[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_larger_kkt() {
        // 4 variables, 2 constraints
        // H = diag(1,2,3,4), J = [[1,0,1,0],[0,1,0,1]]
        // Full 6x6 matrix
        let n = 4;
        let m = 2;
        let dim = n + m;
        let mut mat = SymmetricMatrix::zeros(dim);

        // H diagonal
        for i in 0..n {
            mat.set(i, i, (i + 1) as f64);
        }

        // J entries
        mat.set(4, 0, 1.0); // J[0][0]
        mat.set(4, 2, 1.0); // J[0][2]
        mat.set(5, 1, 1.0); // J[1][1]
        mat.set(5, 3, 1.0); // J[1][3]

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(inertia.positive, 4, "Expected 4 positive eigenvalues");
        assert_eq!(inertia.negative, 2, "Expected 2 negative eigenvalues");
        assert_eq!(inertia.zero, 0);

        // Solve and verify
        let rhs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut sol = [0.0; 6];
        solver.solve(&rhs, &mut sol).unwrap();

        let full = mat.to_full();
        for i in 0..dim {
            let ax_i: f64 = (0..dim).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-9,
                "Row {}: Ax={}, b={}",
                i,
                ax_i,
                rhs[i]
            );
        }
    }

    #[test]
    fn test_negative_definite() {
        let mut mat = SymmetricMatrix::zeros(2);
        mat.set(0, 0, -3.0);
        mat.set(1, 0, -1.0);
        mat.set(1, 1, -4.0);

        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();

        assert_eq!(inertia.positive, 0);
        assert_eq!(inertia.negative, 2);

        let rhs = [1.0, 2.0];
        let mut sol = [0.0; 2];
        solver.solve(&rhs, &mut sol).unwrap();

        let full = mat.to_full();
        for i in 0..2 {
            let ax_i: f64 = (0..2).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i,
                ax_i,
                rhs[i]
            );
        }
    }

    #[test]
    fn test_1x1_matrix() {
        let mut mat = SymmetricMatrix::zeros(1);
        mat.set(0, 0, 5.0);
        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();
        assert_eq!(inertia.positive, 1);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);
        // Solve 5x = 3
        let rhs = [3.0];
        let mut sol = [0.0];
        solver.solve(&rhs, &mut sol).unwrap();
        assert!((sol[0] - 0.6).abs() < 1e-12);
    }

    #[test]
    fn test_zero_matrix() {
        let mat = SymmetricMatrix::zeros(2);
        let mut solver = DenseLdl::new();
        let inertia = solver.factor(&mat).unwrap().unwrap();
        assert_eq!(inertia.positive, 0);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 2);
    }

    #[test]
    fn test_solve_unfactored_error() {
        let mut solver = DenseLdl::new();
        let rhs = [1.0, 2.0];
        let mut sol = [0.0; 2];
        let result = solver.solve(&rhs, &mut sol);
        assert!(result.is_err());
    }

    #[test]
    fn test_dimension_mismatch() {
        let mut mat = SymmetricMatrix::zeros(3);
        for i in 0..3 {
            mat.set(i, i, 1.0);
        }
        let mut solver = DenseLdl::new();
        solver.factor(&mat).unwrap();
        let rhs = [1.0, 2.0]; // Wrong size (2 instead of 3)
        let mut sol = [0.0; 2];
        let result = solver.solve(&rhs, &mut sol);
        assert!(result.is_err());
    }

    #[test]
    fn test_round_trip_pivot_swap() {
        // Matrix with zero diagonal to force Bunch-Kaufman pivoting
        // A = [[0, 3, 1], [3, 1, 2], [1, 2, 4]]
        let mut mat = SymmetricMatrix::zeros(3);
        mat.set(0, 0, 0.0);
        mat.set(1, 0, 3.0);
        mat.set(1, 1, 1.0);
        mat.set(2, 0, 1.0);
        mat.set(2, 1, 2.0);
        mat.set(2, 2, 4.0);

        let mut solver = DenseLdl::new();
        solver.factor(&mat).unwrap();

        let rhs = [1.0, 2.0, 3.0];
        let mut sol = [0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify A * sol = rhs
        let full = mat.to_full();
        for i in 0..3 {
            let ax_i: f64 = (0..3).map(|j| full[i][j] * sol[j]).sum();
            assert!(
                (ax_i - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i, ax_i, rhs[i]
            );
        }
    }
}
