use super::{Inertia, KktMatrix, LinearSolver, SolverError};
use faer::sparse::linalg::cholesky::{
    factorize_symbolic_cholesky, LdltRef, SymbolicCholesky,
};
use faer::linalg::cholesky::ldlt_diagonal::compute::LdltRegularization;
use faer::Side;
use faer::dyn_stack::{PodStack, GlobalPodBuffer};
use faer::Parallelism;
use faer::Mat;
use faer::Conj;

/// Sparse LDL^T factorization using faer's simplicial solver with AMD ordering.
///
/// Uses COO assembly → CSC conversion → symbolic factorization (cached) →
/// numeric LDLT factorization each iteration. Extracts D diagonal for inertia.
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
        }
    }

    /// Extract inertia from D diagonal of simplicial LDLT factorization.
    ///
    /// In faer's simplicial LDLT, the values array stores L entries in CSC format
    /// where the first entry of each column is the D diagonal value:
    ///   values[col_ptrs[j]] = D[j,j]
    fn compute_inertia(&self) -> Inertia {
        let symbolic = self.symbolic.as_ref().unwrap();
        let mut positive = 0;
        let mut negative = 0;
        let mut zero = 0;

        // Access the simplicial structure to get column pointers
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
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Supernodal(_) => {
                // Should not happen since we force simplicial mode.
                // Fall back to reporting unknown inertia.
                return Inertia {
                    positive: 0,
                    negative: 0,
                    zero: self.n,
                };
            }
        }

        Inertia {
            positive,
            negative,
            zero,
        }
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

        // Convert COO to CSC (upper triangle)
        let csc = sparse.to_upper_csc();

        // Compute symbolic factorization on first call (sparsity pattern is fixed for NLP)
        if self.symbolic.is_none() {
            let symbolic = factorize_symbolic_cholesky(
                csc.symbolic(),
                Side::Upper,
                Default::default(), // AMD ordering
                Default::default(), // CholeskySymbolicParams with simplicial threshold
            )
            .map_err(|e| SolverError::NumericalFailure(format!("symbolic factorization: {:?}", e)))?;

            self.l_values = vec![0.0; symbolic.len_values()];
            self.symbolic = Some(symbolic);
        }

        let symbolic = self.symbolic.as_ref().unwrap();

        // Compute required stack size for numeric factorization
        let req = symbolic
            .factorize_numeric_ldlt_req::<f64>(false, Parallelism::None)
            .map_err(|e| SolverError::NumericalFailure(format!("stack req: {:?}", e)))?;
        let mut mem = GlobalPodBuffer::new(req);
        let stack = PodStack::new(&mut mem);

        // No regularization — let the inertia correction loop handle perturbation
        let reg = LdltRegularization {
            dynamic_regularization_signs: None,
            dynamic_regularization_delta: 0.0,
            dynamic_regularization_epsilon: 0.0,
        };

        // Numeric LDLT factorization
        let _ldlt = symbolic.factorize_numeric_ldlt::<f64>(
            &mut self.l_values,
            csc.as_ref(),
            Side::Upper,
            reg,
            Parallelism::None,
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

        // Compute required stack size for solve
        let req = symbolic
            .solve_in_place_req::<f64>(1)
            .map_err(|e| SolverError::NumericalFailure(format!("solve stack req: {:?}", e)))?;
        let mut mem = GlobalPodBuffer::new(req);
        let stack = PodStack::new(&mut mem);

        ldlt.solve_in_place_with_conj(Conj::No, rhs_mat.as_mut(), Parallelism::None, stack);

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
        match symbolic.raw() {
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Simplicial(ref s) => {
                let col_ptrs = s.col_ptrs();
                let mut min_d = f64::INFINITY;
                for j in 0..self.n {
                    min_d = min_d.min(self.l_values[col_ptrs[j]]);
                }
                Some(min_d)
            }
            _ => None,
        }
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
}
