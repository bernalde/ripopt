use super::{Inertia, KktMatrix, LinearSolver, SolverError};

/// Multifrontal sparse symmetric indefinite solver backed by rmumps.
///
/// Wraps `rmumps::Solver` to implement ripopt's `LinearSolver` trait.
/// Uses AMD ordering and iterative refinement by default.
/// Caches symbolic factorization and COO→CSC mapping across calls
/// with the same sparsity pattern to avoid redundant allocation.
pub struct MultifrontalLdl {
    solver: rmumps::solver::Solver,
    n: usize,
    factored: bool,

    // --- Cached CSC structure and COO→CSC mapping ---
    /// Cached CSC matrix (structure fixed, values updated in-place).
    csc: Option<rmumps::csc::CscMatrix>,
    /// Mapping from COO triplet index → CSC value index.
    coo_to_csc: Vec<usize>,
}

impl Default for MultifrontalLdl {
    fn default() -> Self {
        Self::new()
    }
}

impl MultifrontalLdl {
    pub fn new() -> Self {
        Self {
            solver: rmumps::solver::Solver::new(rmumps::solver::SolverOptions::default()),
            n: 0,
            factored: false,
            csc: None,
            coo_to_csc: Vec::new(),
        }
    }

    /// Create a new solver configured for KKT systems with n_primal primal variables.
    /// Uses KKT matching ordering to place primal-dual pairs adjacent in the
    /// elimination order, enabling stable 2x2 pivots for zero-diagonal equality
    /// constraint rows within the FS-only pivot search.
    pub fn new_kkt(n_primal: usize) -> Self {
        let options = rmumps::solver::SolverOptions {
            n_primal: Some(n_primal),
            ordering: rmumps::ordering::Ordering::KktMatchingAmd,
            ..rmumps::solver::SolverOptions::default()
        };
        Self {
            solver: rmumps::solver::Solver::new(options),
            n: 0,
            factored: false,
            csc: None,
            coo_to_csc: Vec::new(),
        }
    }

    /// Convert rmumps Inertia to ripopt Inertia.
    fn convert_inertia(inertia: rmumps::Inertia) -> Inertia {
        Inertia {
            positive: inertia.positive,
            negative: inertia.negative,
            zero: inertia.zero,
        }
    }

    fn map_error(e: rmumps::SolverError) -> SolverError {
        match e {
            rmumps::SolverError::SingularMatrix => SolverError::SingularMatrix,
            rmumps::SolverError::NumericalFailure(msg) => SolverError::NumericalFailure(msg),
            rmumps::SolverError::DimensionMismatch { expected, got } => {
                SolverError::DimensionMismatch { expected, got }
            }
            other => SolverError::NumericalFailure(format!("{}", other)),
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

impl LinearSolver for MultifrontalLdl {
    fn factor(&mut self, matrix: &KktMatrix) -> Result<Option<Inertia>, SolverError> {
        let sparse = match matrix {
            KktMatrix::Sparse(s) => s,
            KktMatrix::Dense(_) => {
                return Err(SolverError::NumericalFailure(
                    "MultifrontalLdl requires KktMatrix::Sparse".into(),
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

        let first_call = self.csc.is_none();

        if first_call {
            // First call: build COO→CSC, cache mapping, do analyze+factor
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

            // analyze uses the CSC structure, factor uses the values
            self.solver.analyze(&coo).map_err(Self::map_error)?;
            let rmumps_inertia = self.solver.factor_csc(&csc).map_err(Self::map_error)?;

            self.csc = Some(csc);
            self.factored = true;
            Ok(Some(Self::convert_inertia(rmumps_inertia)))
        } else {
            // Subsequent calls: update CSC values via cached mapping, refactor
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

            let rmumps_inertia = self.solver.factor_csc(csc).map_err(Self::map_error)?;

            self.factored = true;
            Ok(Some(Self::convert_inertia(rmumps_inertia)))
        }
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

        self.solver.solve(rhs, solution).map_err(Self::map_error)
    }

    fn provides_inertia(&self) -> bool {
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        self.solver.min_diagonal()
    }
}
