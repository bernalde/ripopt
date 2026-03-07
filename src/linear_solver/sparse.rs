use super::{Inertia, KktMatrix, LinearSolver, SolverError};
use faer::Index;
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky, LdltRef, SymbolicCholesky,
};
use faer::linalg::cholesky::ldlt_diagonal::compute::LdltRegularization;
use faer::Side;
use faer::dyn_stack::{PodStack, GlobalPodBuffer};
use faer::Parallelism;
use faer::Mat;
use faer::Conj;

/// Sparse LDL^T factorization using faer with AMD ordering.
///
/// Uses COO assembly → CSC conversion → symbolic factorization (cached) →
/// numeric LDLT factorization each iteration. Extracts D diagonal for inertia.
/// Automatically selects simplicial or supernodal mode based on fill-in analysis.
/// Supernodal mode uses dense BLAS operations on column groups for better cache
/// efficiency on large problems.
///
/// Performance: after the first factorization, the CSC structure, symbolic
/// factorization, COO→CSC mapping, and workspace buffers are all cached.
/// Subsequent calls only update CSC values and run numeric factorization.
pub struct SparseLdl {
    /// Dimension of the factored matrix.
    n: usize,
    /// Cached symbolic factorization (computed once, reused).
    symbolic: Option<SymbolicCholesky<usize>>,
    /// L values buffer (filled by numeric factorization).
    l_values: Vec<f64>,
    /// Whether the matrix has been factored.
    factored: bool,
    /// Tolerance for determining zero pivots.
    zero_pivot_tol: f64,
    /// Force supernodal mode (for testing).
    force_supernodal: bool,

    // --- Cached CSC structure and mapping (optimization #1 and #3) ---
    /// Cached CSC matrix (structure is fixed for NLP, values updated in-place).
    csc: Option<faer::sparse::SparseColMat<usize, f64>>,
    /// Mapping from COO triplet index → CSC value index.
    /// For each triplet k, coo_to_csc[k] is the index into csc_values where
    /// the triplet's value should be accumulated (duplicates sum).
    coo_to_csc: Vec<usize>,

    // --- Cached workspace buffers (optimization #2) ---
    /// Workspace for numeric factorization.
    factor_buf: Option<GlobalPodBuffer>,
    /// Workspace for solve.
    solve_buf: Option<GlobalPodBuffer>,
}

impl Default for SparseLdl {
    fn default() -> Self {
        Self::new()
    }
}

impl SparseLdl {
    pub fn new() -> Self {
        Self {
            n: 0,
            symbolic: None,
            l_values: Vec::new(),
            factored: false,
            zero_pivot_tol: 1e-12,
            force_supernodal: false,
            csc: None,
            coo_to_csc: Vec::new(),
            factor_buf: None,
            solve_buf: None,
        }
    }

    /// Extract inertia from D diagonal of LDLT factorization.
    ///
    /// Handles both simplicial mode (D stored as first entry per CSC column)
    /// and supernodal mode (D stored as diagonal of dense supernode blocks).
    fn compute_inertia(&self) -> Inertia {
        let symbolic = self.symbolic.as_ref().unwrap();
        let mut positive = 0;
        let mut negative = 0;
        let mut zero = 0;

        match symbolic.raw() {
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Simplicial(ref s) => {
                let col_ptrs = s.col_ptrs();
                for j in 0..self.n {
                    let d_j = self.l_values[col_ptrs[j]];
                    if d_j > self.zero_pivot_tol {
                        positive += 1;
                    } else if d_j < -self.zero_pivot_tol {
                        negative += 1;
                    } else {
                        zero += 1;
                    }
                }
            }
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Supernodal(ref s) => {
                let n_supernodes = s.n_supernodes();
                let sn_begin = s.supernode_begin();
                let sn_end_slice = s.supernode_end();
                let col_ptrs_ri = s.col_ptrs_for_row_indices();
                let col_ptrs_val = s.col_ptrs_for_values();
                for sn in 0..n_supernodes {
                    let sn_start = sn_begin[sn].zx();
                    let sn_end = sn_end_slice[sn].zx();
                    let sn_ncols = sn_end - sn_start;
                    let pattern_len = col_ptrs_ri[sn + 1].zx() - col_ptrs_ri[sn].zx();
                    let sn_nrows = pattern_len + sn_ncols;
                    let val_start = col_ptrs_val[sn].zx();

                    // D diagonal is at position (j, j) in the dense column-major block
                    for j in 0..sn_ncols {
                        let d_j = self.l_values[val_start + j * sn_nrows + j];
                        if d_j > self.zero_pivot_tol {
                            positive += 1;
                        } else if d_j < -self.zero_pivot_tol {
                            negative += 1;
                        } else {
                            zero += 1;
                        }
                    }
                }
            }
        }

        Inertia {
            positive,
            negative,
            zero,
        }
    }

    /// Build COO→CSC mapping for the given triplets and CSC structure.
    /// For each COO triplet (row, col, val), finds the CSC value index where
    /// the value should be accumulated. Duplicate COO entries mapping to the
    /// same CSC position will sum their values.
    fn build_coo_to_csc_mapping(
        triplet_rows: &[usize],
        triplet_cols: &[usize],
        csc: &faer::sparse::SparseColMat<usize, f64>,
    ) -> Vec<usize> {
        let nnz_coo = triplet_rows.len();
        let mut mapping = Vec::with_capacity(nnz_coo);
        let sym = csc.symbolic();
        let col_ptrs = sym.col_ptrs();
        let row_indices = sym.row_indices();

        for k in 0..nnz_coo {
            let row = triplet_rows[k];
            let col = triplet_cols[k];

            // Find the CSC position for (row, col) in upper triangle
            let col_start = col_ptrs[col];
            let col_end = col_ptrs[col + 1];

            // Binary search for row in this column's row indices
            let slice = &row_indices[col_start..col_end];
            let pos = slice.binary_search(&row)
                .unwrap_or_else(|_| panic!(
                    "COO entry ({}, {}) not found in CSC structure", row, col
                ));
            mapping.push(col_start + pos);
        }

        mapping
    }

    /// Update CSC values from COO triplets using the cached mapping.
    /// Zeros all values first, then accumulates from triplets.
    fn scatter_coo_to_csc(
        coo_to_csc: &[usize],
        triplet_vals: &[f64],
        csc_values: &mut [f64],
    ) {
        // Zero out CSC values
        for v in csc_values.iter_mut() {
            *v = 0.0;
        }
        // Scatter-add from COO triplets using the first triplet_vals.len() mapping entries
        for (k, &val) in triplet_vals.iter().enumerate() {
            csc_values[coo_to_csc[k]] += val;
        }
    }
}

impl SparseLdl {
    /// Force supernodal mode for testing. Must be called before first `factor()`.
    #[cfg(test)]
    fn force_supernodal(&mut self) {
        self.force_supernodal = true;
    }
}

impl LinearSolver for SparseLdl {
    fn factor(&mut self, matrix: &KktMatrix) -> Result<Option<Inertia>, SolverError> {
        let sparse = match matrix {
            KktMatrix::Sparse(s) => s,
            KktMatrix::Dense(_) => {
                return Err(SolverError::NumericalFailure(
                    "SparseLdl requires KktMatrix::Sparse".into(),
                ))
            }
        };

        self.n = sparse.n;

        if self.n == 0 {
            self.factored = true;
            return Ok(Some(Inertia {
                positive: 0,
                negative: 0,
                zero: 0,
            }));
        }

        let first_call = self.symbolic.is_none();

        if first_call {
            // First call: build CSC from COO, compute symbolic factorization,
            // build COO→CSC mapping, allocate workspace buffers.

            // Build CSC via faer's triplet converter
            let triplets: Vec<(usize, usize, f64)> = sparse.triplet_rows.iter()
                .zip(sparse.triplet_cols.iter())
                .zip(sparse.triplet_vals.iter())
                .map(|((&r, &c), &v)| (r, c, v))
                .collect();
            let csc = faer::sparse::SparseColMat::<usize, f64>::try_new_from_triplets(
                self.n, self.n, &triplets,
            ).map_err(|e| SolverError::NumericalFailure(format!("CSC conversion: {:?}", e)))?;

            // Build COO → CSC mapping
            self.coo_to_csc = Self::build_coo_to_csc_mapping(
                &sparse.triplet_rows,
                &sparse.triplet_cols,
                &csc,
            );

            // Symbolic factorization
            use faer::sparse::linalg::cholesky::CholeskySymbolicParams;
            use faer::sparse::linalg::SupernodalThreshold;
            // Force supernodal for large systems. faer's AUTO threshold (40.0) is too
            // high for KKT systems, causing simplicial mode to be selected even at 20K+
            // dimensions. Simplicial does O(n) random memory accesses per sparse update,
            // while supernodal uses dense BLAS on column groups with good cache locality.
            // This is the root cause of a 70,000x gap vs MUMPS at scale.
            let threshold = if self.force_supernodal || self.n >= 500 {
                SupernodalThreshold::FORCE_SUPERNODAL
            } else {
                SupernodalThreshold::AUTO
            };
            let symbolic = factorize_symbolic_cholesky(
                csc.symbolic(),
                Side::Upper,
                Default::default(), // AMD ordering
                CholeskySymbolicParams {
                    supernodal_flop_ratio_threshold: threshold,
                    ..Default::default()
                },
            )
            .map_err(|e| SolverError::NumericalFailure(format!("symbolic factorization: {:?}", e)))?;

            self.l_values = vec![0.0; symbolic.len_values()];

            // Pre-allocate workspace buffers
            let par = Parallelism::Rayon(0);
            let factor_req = symbolic
                .factorize_numeric_ldlt_req::<f64>(false, par)
                .map_err(|e| SolverError::NumericalFailure(format!("stack req: {:?}", e)))?;
            self.factor_buf = Some(GlobalPodBuffer::new(factor_req));

            let solve_req = symbolic
                .solve_in_place_req::<f64>(1)
                .map_err(|e| SolverError::NumericalFailure(format!("solve stack req: {:?}", e)))?;
            self.solve_buf = Some(GlobalPodBuffer::new(solve_req));

            self.symbolic = Some(symbolic);
            self.csc = Some(csc);
        } else {
            // Subsequent calls: update CSC values from COO triplets via cached mapping.
            let new_nnz = sparse.triplet_rows.len();
            if new_nnz > self.coo_to_csc.len() {
                // Extra entries added (e.g. add_diagonal_range for inertia correction).
                // Extend the mapping for the new triplets.
                let extra_mapping = Self::build_coo_to_csc_mapping(
                    &sparse.triplet_rows[self.coo_to_csc.len()..],
                    &sparse.triplet_cols[self.coo_to_csc.len()..],
                    self.csc.as_ref().unwrap(),
                );
                self.coo_to_csc.extend(extra_mapping);
            }
            let csc = self.csc.as_mut().unwrap();
            let csc_values = csc.values_mut();
            Self::scatter_coo_to_csc(&self.coo_to_csc, &sparse.triplet_vals, csc_values);
        }

        let csc = self.csc.as_ref().unwrap();
        let symbolic = self.symbolic.as_ref().unwrap();
        let par = Parallelism::Rayon(0);

        // No regularization — let the inertia correction loop handle perturbation
        let reg = LdltRegularization {
            dynamic_regularization_signs: None,
            dynamic_regularization_delta: 0.0,
            dynamic_regularization_epsilon: 0.0,
        };

        // Numeric LDLT factorization using cached workspace
        let stack = PodStack::new(self.factor_buf.as_mut().unwrap());
        let _ldlt = symbolic.factorize_numeric_ldlt::<f64>(
            &mut self.l_values,
            csc.as_ref(),
            Side::Upper,
            reg,
            par,
            stack,
        );

        self.factored = true;

        // Extract inertia from D diagonal
        let inertia = self.compute_inertia();
        Ok(Some(inertia))
    }

    fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        if !self.factored {
            return Err(SolverError::NumericalFailure(
                "matrix not factored".to_string(),
            ));
        }

        if rhs.len() != self.n || solution.len() != self.n {
            return Err(SolverError::DimensionMismatch {
                expected: self.n,
                got: rhs.len(),
            });
        }

        if self.n == 0 {
            return Ok(());
        }

        let symbolic = self.symbolic.as_ref().unwrap();
        let ldlt = LdltRef::new(symbolic, self.l_values.as_slice());

        // Copy rhs into faer::Mat (column vector)
        let mut rhs_mat = Mat::<f64>::from_fn(self.n, 1, |i, _| rhs[i]);

        // Solve using cached workspace buffer
        let stack = PodStack::new(self.solve_buf.as_mut().unwrap());
        ldlt.solve_in_place_with_conj(Conj::No, rhs_mat.as_mut(), Parallelism::Rayon(0), stack);

        for i in 0..self.n {
            solution[i] = rhs_mat[(i, 0)];
        }

        Ok(())
    }

    fn provides_inertia(&self) -> bool {
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        if !self.factored || self.n == 0 {
            return None;
        }

        let symbolic = self.symbolic.as_ref()?;
        let mut min_d = f64::INFINITY;
        match symbolic.raw() {
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Simplicial(ref s) => {
                let col_ptrs = s.col_ptrs();
                for j in 0..self.n {
                    min_d = min_d.min(self.l_values[col_ptrs[j]]);
                }
            }
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Supernodal(ref s) => {
                let n_supernodes = s.n_supernodes();
                let sn_begin = s.supernode_begin();
                let sn_end_slice = s.supernode_end();
                let col_ptrs_ri = s.col_ptrs_for_row_indices();
                let col_ptrs_val = s.col_ptrs_for_values();
                for sn in 0..n_supernodes {
                    let sn_start = sn_begin[sn].zx();
                    let sn_end = sn_end_slice[sn].zx();
                    let sn_ncols = sn_end - sn_start;
                    let pattern_len = col_ptrs_ri[sn + 1].zx() - col_ptrs_ri[sn].zx();
                    let sn_nrows = pattern_len + sn_ncols;
                    let val_start = col_ptrs_val[sn].zx();
                    for j in 0..sn_ncols {
                        min_d = min_d.min(self.l_values[val_start + j * sn_nrows + j]);
                    }
                }
            }
        }
        Some(min_d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_solver::SparseSymmetricMatrix;

    #[test]
    fn test_sparse_ldl_positive_definite() {
        // A = [[4, 2, 1], [2, 5, 3], [1, 3, 6]]
        let mut s = SparseSymmetricMatrix::zeros(3);
        s.add(0, 0, 4.0);
        s.add(1, 0, 2.0);
        s.add(1, 1, 5.0);
        s.add(2, 0, 1.0);
        s.add(2, 1, 3.0);
        s.add(2, 2, 6.0);

        let matrix = KktMatrix::Sparse(s);
        let mut solver = SparseLdl::new();
        let inertia = solver.factor(&matrix).unwrap().unwrap();

        assert_eq!(inertia.positive, 3, "Expected 3 positive, got {:?}", inertia);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);

        // Solve Ax = b where b = [1, 2, 3]
        let rhs = [1.0, 2.0, 3.0];
        let mut sol = [0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        // Verify Ax ≈ b using COO matvec
        let mut ax = [0.0; 3];
        matrix.matvec(&sol, &mut ax);
        for i in 0..3 {
            assert!(
                (ax[i] - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i, ax[i], rhs[i]
            );
        }
    }

    #[test]
    fn test_sparse_ldl_kkt_matrix() {
        // KKT matrix: [[2, 0, 1], [0, 2, 1], [1, 1, 0]]
        // Should have inertia (2, 1, 0)
        let mut s = SparseSymmetricMatrix::zeros(3);
        s.add(0, 0, 2.0);
        s.add(1, 1, 2.0);
        s.add(2, 0, 1.0);
        s.add(2, 1, 1.0);

        let matrix = KktMatrix::Sparse(s);
        let mut solver = SparseLdl::new();
        let inertia = solver.factor(&matrix).unwrap().unwrap();

        assert_eq!(inertia.positive, 2, "Expected 2 positive, got {:?}", inertia);
        assert_eq!(inertia.negative, 1, "Expected 1 negative, got {:?}", inertia);
        assert_eq!(inertia.zero, 0);

        // Solve
        let rhs = [1.0, 2.0, 3.0];
        let mut sol = [0.0; 3];
        solver.solve(&rhs, &mut sol).unwrap();

        let mut ax = [0.0; 3];
        matrix.matvec(&sol, &mut ax);
        for i in 0..3 {
            assert!(
                (ax[i] - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i, ax[i], rhs[i]
            );
        }
    }

    #[test]
    fn test_sparse_ldl_identity() {
        let mut s = SparseSymmetricMatrix::zeros(4);
        for i in 0..4 {
            s.add(i, i, 1.0);
        }

        let matrix = KktMatrix::Sparse(s);
        let mut solver = SparseLdl::new();
        let inertia = solver.factor(&matrix).unwrap().unwrap();
        assert_eq!(inertia.positive, 4);

        let rhs = [1.0, 2.0, 3.0, 4.0];
        let mut sol = [0.0; 4];
        solver.solve(&rhs, &mut sol).unwrap();

        for i in 0..4 {
            assert!((sol[i] - rhs[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_sparse_dense_inertia_match() {
        // Verify sparse and dense solvers report same inertia for a known KKT matrix
        use crate::linear_solver::dense::DenseLdl;
        use crate::linear_solver::SymmetricMatrix;

        // Build same matrix in both formats
        let mut dense = SymmetricMatrix::zeros(3);
        let mut sparse = SparseSymmetricMatrix::zeros(3);

        let entries = [(0, 0, 2.0), (1, 1, 2.0), (2, 0, 1.0), (2, 1, 1.0)];
        for &(i, j, v) in &entries {
            dense.set(i, j, v);
            sparse.add(i, j, v);
        }

        let mut dense_solver = DenseLdl::new();
        let dense_inertia = dense_solver
            .factor(&KktMatrix::Dense(dense))
            .unwrap()
            .unwrap();

        let mut sparse_solver = SparseLdl::new();
        let sparse_inertia = sparse_solver
            .factor(&KktMatrix::Sparse(sparse))
            .unwrap()
            .unwrap();

        assert_eq!(
            dense_inertia, sparse_inertia,
            "Dense inertia {:?} != Sparse inertia {:?}",
            dense_inertia, sparse_inertia
        );
    }

    #[test]
    fn test_supernodal_kkt_matrix() {
        // KKT-shaped matrix that exercises the supernodal code path.
        // Structure: [H, A^T; A, 0] with n=5 vars and m=3 constraints.
        let n = 5;
        let m = 3;
        let dim = n + m;
        let mut s = SparseSymmetricMatrix::zeros(dim);

        // H block: positive definite diagonal + some off-diagonal
        for i in 0..n {
            s.add(i, i, 10.0 + i as f64);
        }
        s.add(1, 0, 1.0);
        s.add(2, 0, 0.5);
        s.add(2, 1, 0.5);

        // A block (rows n..n+m, cols 0..n): constraint Jacobian
        s.add(n, 0, 1.0);
        s.add(n, 1, 2.0);
        s.add(n + 1, 1, 1.0);
        s.add(n + 1, 2, 3.0);
        s.add(n + 2, 3, 1.0);
        s.add(n + 2, 4, 1.0);

        // Small negative perturbation on (2,2) block (like IPM inertia correction)
        for i in n..dim {
            s.add(i, i, -1e-8);
        }

        let matrix = KktMatrix::Sparse(s.clone());

        // First verify with simplicial (known working)
        let mut simp_solver = SparseLdl::new();
        let simp_inertia = simp_solver.factor(&KktMatrix::Sparse(s)).unwrap().unwrap();

        // Solve with forced supernodal mode
        let mut solver = SparseLdl::new();
        solver.force_supernodal();
        let inertia = solver.factor(&matrix).unwrap().unwrap();

        // Supernodal inertia should match simplicial
        assert_eq!(inertia, simp_inertia,
            "Supernodal {:?} != Simplicial {:?}", inertia, simp_inertia);

        // Solve and verify
        let rhs: Vec<f64> = (1..=dim).map(|i| i as f64).collect();
        let mut sol = vec![0.0; dim];
        solver.solve(&rhs, &mut sol).unwrap();

        let mut ax = vec![0.0; dim];
        matrix.matvec(&sol, &mut ax);
        for i in 0..dim {
            assert!(
                (ax[i] - rhs[i]).abs() < 1e-4,
                "Row {}: Ax={}, b={}, err={}",
                i, ax[i], rhs[i], (ax[i] - rhs[i]).abs()
            );
        }
    }

    #[test]
    fn test_supernodal_positive_definite() {
        // Positive definite banded matrix, forced supernodal
        let n = 20;
        let mut s = SparseSymmetricMatrix::zeros(n);
        for i in 0..n {
            s.add(i, i, 4.0);
            if i > 0 { s.add(i, i - 1, -1.0); }
            if i > 1 { s.add(i, i - 2, -0.5); }
        }

        let matrix = KktMatrix::Sparse(s);
        let mut solver = SparseLdl::new();
        solver.force_supernodal();
        let inertia = solver.factor(&matrix).unwrap().unwrap();

        assert_eq!(inertia.positive, n);
        assert_eq!(inertia.negative, 0);
        assert_eq!(inertia.zero, 0);

        let rhs: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut sol = vec![0.0; n];
        solver.solve(&rhs, &mut sol).unwrap();

        let mut ax = vec![0.0; n];
        matrix.matvec(&sol, &mut ax);
        for i in 0..n {
            assert!(
                (ax[i] - rhs[i]).abs() < 1e-8,
                "Row {}: Ax={}, b={}",
                i, ax[i], rhs[i]
            );
        }
    }

    #[test]
    fn test_cached_refactorization() {
        // Test that factoring the same solver twice (simulating iteration 2)
        // gives correct results when using the cached CSC mapping.
        let mut s = SparseSymmetricMatrix::zeros(3);
        s.add(0, 0, 4.0);
        s.add(1, 0, 2.0);
        s.add(1, 1, 5.0);
        s.add(2, 0, 1.0);
        s.add(2, 1, 3.0);
        s.add(2, 2, 6.0);

        let matrix1 = KktMatrix::Sparse(s);
        let mut solver = SparseLdl::new();

        // First factorization
        solver.factor(&matrix1).unwrap();
        let rhs = [1.0, 2.0, 3.0];
        let mut sol1 = [0.0; 3];
        solver.solve(&rhs, &mut sol1).unwrap();

        // Second factorization with different values but same structure
        let mut s2 = SparseSymmetricMatrix::zeros(3);
        s2.add(0, 0, 10.0);
        s2.add(1, 0, 1.0);
        s2.add(1, 1, 10.0);
        s2.add(2, 0, 0.5);
        s2.add(2, 1, 1.0);
        s2.add(2, 2, 10.0);

        let matrix2 = KktMatrix::Sparse(s2.clone());
        solver.factor(&matrix2).unwrap();
        let mut sol2 = [0.0; 3];
        solver.solve(&rhs, &mut sol2).unwrap();

        // Verify Ax ≈ b for second solve
        let mut ax = [0.0; 3];
        matrix2.matvec(&sol2, &mut ax);
        for i in 0..3 {
            assert!(
                (ax[i] - rhs[i]).abs() < 1e-10,
                "Row {}: Ax={}, b={}",
                i, ax[i], rhs[i]
            );
        }

        // Solutions should differ (different matrices)
        assert!((sol1[0] - sol2[0]).abs() > 1e-6, "Solutions should differ");
    }
}
