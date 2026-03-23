use crate::coo::CooMatrix;
use crate::csc::CscMatrix;
use crate::numeric::{multifrontal_factor, multifrontal_factor_threshold, NumericFactorization};
use crate::ordering::{self, Ordering};
use crate::scaling::{self, Scaling, ScalingFactors};
use crate::solve::multifrontal_solve;
use crate::symbolic::SymbolicFactorization;
use crate::{Inertia, SolverError};

/// Configuration options for the solver.
#[derive(Debug, Clone)]
pub struct SolverOptions {
    /// Fill-reducing ordering method.
    pub ordering: Ordering,
    /// Number of iterative refinement steps (default: 10, adaptive stopping).
    pub refine_steps: usize,
    /// Matrix scaling method applied before factorization.
    pub scaling: Scaling,
    /// Pivot threshold for delayed pivoting (0.0 to 1.0).
    /// A value of 0.0 disables threshold pivoting (classic Bunch-Kaufman).
    /// Default: 0.01 (matches MA57/MUMPS CNTL(1)).
    pub pivot_threshold: f64,
    /// Number of primal variables for KKT-aware ordering.
    /// When set and ordering is KktMatchingAmd, the solver pairs primal-dual
    /// variables for numerically stable 2×2 pivots.
    pub n_primal: Option<usize>,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            ordering: Ordering::Amd,
            refine_steps: 10,
            scaling: Scaling::Ruiz { max_iter: 10 },
            pivot_threshold: crate::pivot::DEFAULT_PIVOT_THRESHOLD,
            n_primal: None,
        }
    }
}

/// Multifrontal sparse symmetric indefinite solver.
///
/// Usage:
/// 1. Create with `Solver::new(options)`
/// 2. Call `analyze(&coo)` to compute ordering and symbolic factorization (once per sparsity pattern)
/// 3. Call `factor(&coo)` to perform numeric factorization (can repeat with same pattern)
/// 4. Call `solve(rhs, solution)` to solve Ax = b
///
/// Or use `analyze_and_factor(&coo)` for convenience.
pub struct Solver {
    options: SolverOptions,
    /// Permutation from ordering.
    perm: Vec<usize>,
    perm_inv: Vec<usize>,
    /// Symbolic factorization (computed during analyze).
    symbolic: Option<SymbolicFactorization>,
    /// Numeric factorization (computed during factor).
    numeric: Option<NumericFactorization>,
    /// Permuted CSC matrix (stored for iterative refinement).
    permuted_csc: Option<CscMatrix>,
    /// Scaling factors (computed during factor, if scaling is enabled).
    scaling_factors: Option<ScalingFactors>,
}

impl Solver {
    /// Create a new solver with the given options.
    pub fn new(options: SolverOptions) -> Self {
        Self {
            options,
            perm: Vec::new(),
            perm_inv: Vec::new(),
            symbolic: None,
            numeric: None,
            permuted_csc: None,
            scaling_factors: None,
        }
    }

    /// Symbolic analysis: compute fill-reducing ordering and symbolic factorization.
    /// Call once per sparsity pattern.
    pub fn analyze(&mut self, matrix: &CooMatrix) -> Result<(), SolverError> {
        let csc = CscMatrix::from_coo(matrix);
        let (perm, perm_inv) = ordering::compute_ordering_with_kkt(
            &csc, self.options.ordering, self.options.n_primal,
        );
        let permuted_csc = ordering::permute_symmetric_csc(&csc, &perm, &perm_inv);
        let symbolic = SymbolicFactorization::from_csc(&permuted_csc);

        self.perm = perm;
        self.perm_inv = perm_inv;
        self.symbolic = Some(symbolic);
        self.numeric = None;
        Ok(())
    }

    /// Numeric factorization using the most recent symbolic analysis.
    /// The matrix must have the same sparsity pattern as the one passed to `analyze`.
    /// If scaling is enabled, computes scaling factors and applies them before factorization.
    /// Returns the inertia of the matrix.
    pub fn factor(&mut self, matrix: &CooMatrix) -> Result<Inertia, SolverError> {
        let sym = self.symbolic.as_ref().ok_or_else(|| {
            SolverError::InvalidState("must call analyze() before factor()".into())
        })?;

        let csc = CscMatrix::from_coo(matrix);

        // Compute and apply scaling if enabled
        // Use MC64 matching-based scaling for KKT systems (n_primal is set)
        let sf = if self.options.n_primal.is_some() {
            scaling::compute_scaling_kkt(&csc, self.options.scaling)
        } else {
            scaling::compute_scaling(&csc, self.options.scaling)
        };
        let mut permuted_csc = ordering::permute_symmetric_csc(&csc, &self.perm, &self.perm_inv);
        if let Some(ref sf) = sf {
            // Apply scaling in the permuted space: need to permute the scaling vector too
            let mut perm_d = vec![0.0; csc.n];
            for i in 0..csc.n {
                perm_d[self.perm_inv[i]] = sf.d[i];
            }
            let perm_sf = ScalingFactors { d: perm_d };
            perm_sf.scale_csc(&mut permuted_csc);
        }
        self.scaling_factors = sf;

        let numeric = multifrontal_factor_threshold(&permuted_csc, sym, self.options.pivot_threshold, self.options.n_primal);
        let inertia = numeric.inertia;
        self.numeric = Some(numeric);
        self.permuted_csc = Some(permuted_csc);
        Ok(inertia)
    }

    /// Solve Ax = b using the most recent factorization.
    /// If scaling was applied during factorization, the solve automatically
    /// handles scaling/unscaling: solves (DAD)y = Db, returns x = Dy.
    pub fn solve(&self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        let sym = self.symbolic.as_ref().ok_or_else(|| {
            SolverError::InvalidState("must call analyze() before solve()".into())
        })?;
        let num = self.numeric.as_ref().ok_or_else(|| {
            SolverError::InvalidState("must call factor() before solve()".into())
        })?;

        let n = sym.n;
        if rhs.len() != n || solution.len() != n {
            return Err(SolverError::DimensionMismatch {
                expected: n,
                got: rhs.len().min(solution.len()),
            });
        }

        // Scale RHS if scaling is active: b_scaled = D * b
        let scaled_rhs;
        let effective_rhs = if let Some(ref sf) = self.scaling_factors {
            scaled_rhs = {
                let mut sr = vec![0.0; n];
                sf.scale_rhs(rhs, &mut sr);
                sr
            };
            &scaled_rhs
        } else {
            rhs
        };

        // Apply permutation to RHS: permuted_rhs[new_i] = rhs[old_i]
        let mut permuted_rhs = vec![0.0; n];
        for i in 0..n {
            permuted_rhs[self.perm_inv[i]] = effective_rhs[i];
        }

        // Solve in permuted space
        let mut permuted_sol = vec![0.0; n];
        multifrontal_solve(num, sym, &permuted_rhs, &mut permuted_sol)?;

        // Adaptive iterative refinement in permuted space (using the scaled+permuted matrix).
        // Stops early when residual stagnates (ratio > 0.9) or is small enough.
        if self.options.refine_steps > 0 {
            if let Some(ref pcsc) = self.permuted_csc {
                let mut residual = vec![0.0; n];
                let mut correction = vec![0.0; n];
                let mut prev_res_norm = f64::INFINITY;
                for _ in 0..self.options.refine_steps {
                    // r = permuted_rhs - A_perm * permuted_sol
                    pcsc.matvec(&permuted_sol, &mut residual);
                    let mut res_norm: f64 = 0.0;
                    for i in 0..n {
                        residual[i] = permuted_rhs[i] - residual[i];
                        res_norm = res_norm.max(residual[i].abs());
                    }

                    // Stop if residual is small enough
                    if res_norm < 1e-14 {
                        break;
                    }

                    // Stop if refinement has stagnated
                    if res_norm > 0.9 * prev_res_norm {
                        break;
                    }
                    prev_res_norm = res_norm;

                    // Solve A_perm * e = r
                    multifrontal_solve(num, sym, &residual, &mut correction)?;
                    // x += e
                    for i in 0..n {
                        permuted_sol[i] += correction[i];
                    }
                }
            }
        }

        // Apply inverse permutation to solution
        let mut unpermuted = vec![0.0; n];
        for i in 0..n {
            unpermuted[i] = permuted_sol[self.perm_inv[i]];
        }

        // Unscale solution if scaling is active: x = D * y
        if let Some(ref sf) = self.scaling_factors {
            sf.unscale_solution(&unpermuted, solution);
        } else {
            solution.copy_from_slice(&unpermuted);
        }

        Ok(())
    }

    /// Convenience method: analyze + factor in one call.
    /// On first call, performs both symbolic and numeric factorization.
    /// On subsequent calls with the same sparsity pattern, only re-does numeric.
    pub fn analyze_and_factor(&mut self, matrix: &CooMatrix) -> Result<Inertia, SolverError> {
        if self.symbolic.is_none() {
            self.analyze(matrix)?;
        }
        self.factor(matrix)
    }

    /// Numeric factorization from a pre-built CSC matrix (skips COO→CSC conversion).
    /// The CSC must have the same sparsity pattern as the one used during `analyze`.
    /// This is the fast path for IPM inertia correction, where only values change.
    /// If scaling is enabled, computes and applies scaling before factorization.
    pub fn factor_csc(&mut self, csc: &CscMatrix) -> Result<Inertia, SolverError> {
        let sym = self.symbolic.as_ref().ok_or_else(|| {
            SolverError::InvalidState("must call analyze() before factor_csc()".into())
        })?;

        // Compute and apply scaling — use MC64 for KKT systems
        let sf = if self.options.n_primal.is_some() {
            scaling::compute_scaling_kkt(csc, self.options.scaling)
        } else {
            scaling::compute_scaling(csc, self.options.scaling)
        };
        let mut permuted_csc = ordering::permute_symmetric_csc(csc, &self.perm, &self.perm_inv);
        if let Some(ref sf) = sf {
            let mut perm_d = vec![0.0; csc.n];
            for i in 0..csc.n {
                perm_d[self.perm_inv[i]] = sf.d[i];
            }
            let perm_sf = ScalingFactors { d: perm_d };
            perm_sf.scale_csc(&mut permuted_csc);
        }
        self.scaling_factors = sf;

        let numeric = multifrontal_factor_threshold(&permuted_csc, sym, self.options.pivot_threshold, self.options.n_primal);
        let inertia = numeric.inertia;
        self.numeric = Some(numeric);
        self.permuted_csc = Some(permuted_csc);
        Ok(inertia)
    }

    /// Whether the solver has a valid factorization.
    pub fn is_factored(&self) -> bool {
        self.numeric.is_some()
    }

    /// Return the minimum eigenvalue of D across all supernodes.
    /// See [`NumericFactorization::min_diagonal`] for details.
    pub fn min_diagonal(&self) -> Option<f64> {
        self.numeric.as_ref()?.min_diagonal()
    }

    /// Reset the solver, clearing all cached state.
    pub fn reset(&mut self) {
        self.perm.clear();
        self.perm_inv.clear();
        self.symbolic = None;
        self.numeric = None;
        self.permuted_csc = None;
        self.scaling_factors = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_coo(n: usize, triplets: &[(usize, usize, f64)]) -> CooMatrix {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        CooMatrix::new(n, rows, cols, vals).unwrap()
    }

    fn check_residual(coo: &CooMatrix, x: &[f64], b: &[f64], tol: f64) {
        let mut ax = vec![0.0; coo.n];
        coo.matvec(x, &mut ax).unwrap();
        let max_resid: f64 = ax.iter().zip(b).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(max_resid < tol, "max residual {} exceeds {}", max_resid, tol);
    }

    #[test]
    fn test_solver_natural_ordering() {
        let coo = make_coo(3, &[
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ]);
        let mut solver = Solver::new(SolverOptions { ordering: Ordering::Natural, ..Default::default() });
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 3);

        let b = [8.0, 18.0, 25.0];
        let mut x = [0.0; 3];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_amd_ordering() {
        let coo = make_coo(3, &[
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ]);
        let mut solver = Solver::new(SolverOptions::default());
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 3);

        let b = [8.0, 18.0, 25.0];
        let mut x = [0.0; 3];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_kkt_5x5() {
        let coo = make_coo(5, &[
            (0, 0, 4.0), (0, 3, 1.0),
            (1, 1, 5.0), (1, 4, 1.0),
            (2, 2, 6.0), (2, 3, 1.0), (2, 4, 1.0),
            (3, 3, 0.0),
            (4, 4, 0.0),
        ]);
        let mut solver = Solver::new(SolverOptions::default());
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 3);
        assert_eq!(inertia.negative, 2);

        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = [0.0; 5];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_tridiagonal_amd() {
        let coo = make_coo(6, &[
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0), (2, 3, 1.0),
            (3, 3, 4.0), (3, 4, 1.0),
            (4, 4, 4.0), (4, 5, 1.0),
            (5, 5, 4.0),
        ]);
        let mut solver = Solver::new(SolverOptions::default());
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 6);

        let b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut x = [0.0; 6];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_arrow_amd() {
        // Arrow matrix: AMD should reorder for less fill
        let n = 6;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 10.0));
            if i < n - 1 {
                triplets.push((i, n - 1, 1.0));
            }
        }
        let coo = make_coo(n, &triplets);
        let mut solver = Solver::new(SolverOptions::default());
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 6);

        let b = vec![1.0; n];
        let mut x = vec![0.0; n];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_refactor_same_pattern() {
        // Factor, solve, then refactor with different values but same pattern
        let coo1 = make_coo(3, &[
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0),
        ]);
        let mut solver = Solver::new(SolverOptions::default());
        solver.analyze_and_factor(&coo1).unwrap();

        let b1 = [1.0, 2.0, 3.0];
        let mut x1 = [0.0; 3];
        solver.solve(&b1, &mut x1).unwrap();
        check_residual(&coo1, &x1, &b1, 1e-10);

        // Refactor with different values
        let coo2 = make_coo(3, &[
            (0, 0, 10.0), (0, 1, 2.0),
            (1, 1, 10.0), (1, 2, 2.0),
            (2, 2, 10.0),
        ]);
        solver.factor(&coo2).unwrap(); // reuses symbolic

        let b2 = [14.0, 24.0, 14.0];
        let mut x2 = [0.0; 3];
        solver.solve(&b2, &mut x2).unwrap();
        check_residual(&coo2, &x2, &b2, 1e-10);
    }

    #[test]
    fn test_solver_error_no_analyze() {
        let mut solver = Solver::new(SolverOptions::default());
        let coo = make_coo(2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        assert!(solver.factor(&coo).is_err());
    }

    #[test]
    fn test_solver_error_no_factor() {
        let mut solver = Solver::new(SolverOptions::default());
        let coo = make_coo(2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        solver.analyze(&coo).unwrap();
        let b = [1.0, 2.0];
        let mut x = [0.0; 2];
        assert!(solver.solve(&b, &mut x).is_err());
    }

    #[test]
    fn test_solver_8x8_kkt() {
        // Larger KKT: 5 variables, 3 constraints
        // H = diag(2,3,4,5,6), A = [[1,1,0,0,0],[0,0,1,1,0],[0,0,0,1,1]]
        let coo = make_coo(8, &[
            (0, 0, 2.0), (0, 5, 1.0), (0, 6, 0.0),
            (1, 1, 3.0), (1, 5, 1.0),
            (2, 2, 4.0), (2, 6, 1.0),
            (3, 3, 5.0), (3, 6, 1.0), (3, 7, 1.0),
            (4, 4, 6.0), (4, 7, 1.0),
            (5, 5, 0.0),
            (6, 6, 0.0),
            (7, 7, 0.0),
        ]);
        let mut solver = Solver::new(SolverOptions::default());
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 5);
        assert_eq!(inertia.negative, 3);

        let b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let mut x = [0.0; 8];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-10);
    }

    #[test]
    fn test_solver_with_ruiz_scaling() {
        // Ill-conditioned KKT: H = diag(1e6, 1e-6), A = [1, 1]
        let coo = make_coo(3, &[
            (0, 0, 1e6), (1, 1, 1e-6),
            (0, 2, 1.0), (1, 2, 1.0),
            (2, 2, 0.0),
        ]);
        let opts = SolverOptions {
            scaling: crate::scaling::Scaling::Ruiz { max_iter: 10 },
            ..Default::default()
        };
        let mut solver = Solver::new(opts);
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, 2);
        assert_eq!(inertia.negative, 1);

        let b = [1e6 + 1.0, 1e-6 + 1.0, 2.0]; // x = [1, 1, 0]
        let mut x = [0.0; 3];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-4);
    }

    #[test]
    fn test_solver_with_diagonal_scaling() {
        // Same problem with diagonal scaling
        let coo = make_coo(3, &[
            (0, 0, 1e4), (0, 1, 1.0),
            (1, 1, 1e-4), (1, 2, 1e-2),
            (2, 2, 1.0),
        ]);
        let opts = SolverOptions {
            scaling: crate::scaling::Scaling::Diagonal,
            ..Default::default()
        };
        let mut solver = Solver::new(opts);
        solver.analyze_and_factor(&coo).unwrap();

        let b = [1.0, 2.0, 3.0];
        let mut x = [0.0; 3];
        solver.solve(&b, &mut x).unwrap();
        check_residual(&coo, &x, &b, 1e-8);
    }

    #[test]
    fn test_solver_scaling_matches_unscaled() {
        // Well-conditioned problem: scaling should give same answer
        let coo = make_coo(3, &[
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ]);

        let b = [8.0, 18.0, 25.0];

        // Solve without scaling
        let mut solver1 = Solver::new(SolverOptions::default());
        solver1.analyze_and_factor(&coo).unwrap();
        let mut x1 = [0.0; 3];
        solver1.solve(&b, &mut x1).unwrap();

        // Solve with Ruiz scaling
        let opts = SolverOptions {
            scaling: crate::scaling::Scaling::Ruiz { max_iter: 10 },
            ..Default::default()
        };
        let mut solver2 = Solver::new(opts);
        solver2.analyze_and_factor(&coo).unwrap();
        let mut x2 = [0.0; 3];
        solver2.solve(&b, &mut x2).unwrap();

        for i in 0..3 {
            assert!((x1[i] - x2[i]).abs() < 1e-10,
                "x[{}]: unscaled={:.10e}, scaled={:.10e}", i, x1[i], x2[i]);
        }
    }
}
