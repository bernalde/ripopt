use super::{Inertia, KktMatrix, LinearSolver, SolverError};

/// Iterative sparse symmetric solver using preconditioned MINRES.
///
/// Uses an incomplete LDL^T factorization as preconditioner and MINRES
/// as the Krylov solver. Suitable for large problems where the direct
/// multifrontal solver is too expensive.
///
/// Inertia from the incomplete factorization is returned as an approximation.
pub struct IterativeMinres {
    /// Cached CSC matrix (structure fixed, values updated in-place).
    csc: Option<rmumps::csc::CscMatrix>,
    /// Mapping from COO triplet index → CSC value index.
    coo_to_csc: Vec<usize>,
    /// Incomplete LDL^T preconditioner.
    precond: Option<rmumps::incomplete::IncompleteLdlt>,
    /// Inertia from the incomplete factorization (approximate).
    inertia: Option<Inertia>,
    /// Matrix dimension.
    n: usize,
    /// Whether factor() has been called successfully.
    factored: bool,
    /// MINRES tolerance.
    tol: f64,
    /// MINRES max iterations.
    max_iter: usize,
    /// Incomplete LDL^T drop tolerance.
    drop_tolerance: f64,
}

impl Default for IterativeMinres {
    fn default() -> Self {
        Self::new()
    }
}

impl IterativeMinres {
    pub fn new() -> Self {
        Self {
            csc: None,
            coo_to_csc: Vec::new(),
            precond: None,
            inertia: None,
            n: 0,
            factored: false,
            tol: 1e-10,
            max_iter: 2000,
            drop_tolerance: 0.01,
        }
    }

    /// Create with custom MINRES and preconditioner settings.
    pub fn with_options(tol: f64, max_iter: usize, drop_tolerance: f64) -> Self {
        Self {
            tol,
            max_iter,
            drop_tolerance,
            ..Self::new()
        }
    }

    fn convert_inertia(inertia: rmumps::Inertia) -> Inertia {
        Inertia {
            positive: inertia.positive,
            negative: inertia.negative,
            zero: inertia.zero,
        }
    }

    /// Build COO→CSC mapping for the given triplets and CSC structure.
    fn build_coo_to_csc_mapping(
        triplet_rows: &[usize],
        triplet_cols: &[usize],
        csc: &rmumps::csc::CscMatrix,
    ) -> Vec<usize> {
        let mut mapping = Vec::with_capacity(triplet_rows.len());
        for k in 0..triplet_rows.len() {
            let row = triplet_rows[k];
            let col = triplet_cols[k];
            let col_start = csc.col_ptr[col];
            let col_end = csc.col_ptr[col + 1];
            let slice = &csc.row_idx[col_start..col_end];
            let pos = slice
                .binary_search(&row)
                .unwrap_or_else(|_| {
                    panic!("COO entry ({}, {}) not found in CSC structure", row, col)
                });
            mapping.push(col_start + pos);
        }
        mapping
    }

    /// Update CSC values from COO triplets using the cached mapping.
    fn scatter_coo_to_csc(
        coo_to_csc: &[usize],
        triplet_vals: &[f64],
        csc_values: &mut [f64],
    ) {
        for v in csc_values.iter_mut() {
            *v = 0.0;
        }
        for (k, &val) in triplet_vals.iter().enumerate() {
            csc_values[coo_to_csc[k]] += val;
        }
    }
}

impl LinearSolver for IterativeMinres {
    fn factor(&mut self, matrix: &KktMatrix) -> Result<Option<Inertia>, SolverError> {
        let sparse = match matrix {
            KktMatrix::Sparse(s) => s,
            KktMatrix::Dense(_) => {
                return Err(SolverError::NumericalFailure(
                    "IterativeMinres requires KktMatrix::Sparse".into(),
                ))
            }
        };

        self.n = sparse.n;

        if self.n == 0 {
            self.factored = true;
            self.inertia = Some(Inertia { positive: 0, negative: 0, zero: 0 });
            return Ok(self.inertia);
        }

        let first_call = self.csc.is_none();

        if first_call {
            // Build COO → CSC, cache mapping
            let coo = rmumps::coo::CooMatrix::new(
                sparse.n,
                sparse.triplet_rows.clone(),
                sparse.triplet_cols.clone(),
                sparse.triplet_vals.clone(),
            )
            .map_err(|e| SolverError::NumericalFailure(format!("rmumps COO: {}", e)))?;

            let csc = rmumps::csc::CscMatrix::from_coo(&coo);
            self.coo_to_csc = Self::build_coo_to_csc_mapping(
                &sparse.triplet_rows,
                &sparse.triplet_cols,
                &csc,
            );
            self.csc = Some(csc);
        } else {
            // Update CSC values via cached mapping
            let new_nnz = sparse.triplet_rows.len();
            if new_nnz > self.coo_to_csc.len() {
                let extra = Self::build_coo_to_csc_mapping(
                    &sparse.triplet_rows[self.coo_to_csc.len()..],
                    &sparse.triplet_cols[self.coo_to_csc.len()..],
                    self.csc.as_ref().unwrap(),
                );
                self.coo_to_csc.extend(extra);
            }

            let csc = self.csc.as_mut().unwrap();
            Self::scatter_coo_to_csc(&self.coo_to_csc, &sparse.triplet_vals, &mut csc.vals);
        }

        // Build incomplete LDL^T preconditioner
        let csc = self.csc.as_ref().unwrap();
        let sym = rmumps::symbolic::SymbolicFactorization::from_csc(csc);
        let opts = rmumps::incomplete::IncompleteLdltOptions {
            drop_tolerance: self.drop_tolerance,
            pivot_threshold: 0.01,
        };
        let precond = rmumps::incomplete::IncompleteLdlt::new(csc, &sym, &opts);
        let inertia = Self::convert_inertia(precond.inertia());
        self.inertia = Some(inertia);
        self.precond = Some(precond);
        self.factored = true;

        Ok(self.inertia)
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

        let csc = self.csc.as_ref().unwrap();
        let precond = self.precond.as_ref().unwrap();

        // Build matvec closure from CSC (symmetric: upper triangle stored)
        let matvec = |x: &[f64], y: &mut [f64]| {
            y.iter_mut().for_each(|v| *v = 0.0);
            for j in 0..csc.n {
                for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                    let i = csc.row_idx[idx];
                    let v = csc.vals[idx];
                    y[i] += v * x[j];
                    if i != j {
                        y[j] += v * x[i];
                    }
                }
            }
        };

        solution.iter_mut().for_each(|v| *v = 0.0);

        let opts = rmumps::minres::MinresOptions {
            max_iter: self.max_iter,
            tol: self.tol,
        };

        let result = rmumps::minres::minres(
            self.n,
            matvec,
            Some(precond as &dyn rmumps::precond::Preconditioner),
            rhs,
            solution,
            &opts,
        );

        if result.converged {
            Ok(())
        } else {
            Err(SolverError::NumericalFailure(format!(
                "MINRES did not converge: {} iterations, residual = {:.2e}",
                result.iterations, result.residual_norm
            )))
        }
    }

    fn provides_inertia(&self) -> bool {
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        // Could extract from incomplete factorization, but it's approximate
        None
    }
}
