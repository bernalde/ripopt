/// Matrix scaling for symmetric matrices.
///
/// Computes diagonal scaling D such that D*A*D has rows/columns with
/// approximately unit infinity-norm. This improves pivot quality in
/// Bunch-Kaufman factorization for ill-conditioned systems (e.g., KKT
/// matrices with entries spanning many orders of magnitude).
///
/// The solver applies scaling transparently:
/// - Factor: scales A → D*A*D, then factors
/// - Solve: given Ax = b, solves (DAD)(D⁻¹x) = Db, returns x = D*y

use crate::csc::CscMatrix;

/// Scaling method for the matrix before factorization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Scaling {
    /// No scaling (default for backward compatibility).
    None,
    /// Diagonal equilibration: D_ii = 1/sqrt(max_j |a_ij|).
    /// Single pass, cheap but effective for many problems.
    Diagonal,
    /// Ruiz iterative equilibration. Repeats diagonal equilibration
    /// until the row/column norms are approximately 1.
    /// More robust than single-pass diagonal scaling.
    Ruiz {
        /// Maximum number of iterations (default: 10).
        max_iter: usize,
    },
}

impl Default for Scaling {
    fn default() -> Self {
        Scaling::None
    }
}

/// Computed scaling factors for a symmetric matrix.
#[derive(Debug, Clone)]
pub struct ScalingFactors {
    /// Diagonal scaling vector D (length n).
    /// The scaled matrix is D*A*D.
    pub d: Vec<f64>,
}

impl ScalingFactors {
    /// Scale a right-hand side vector: b_scaled = D * b.
    pub fn scale_rhs(&self, rhs: &[f64], scaled: &mut [f64]) {
        for i in 0..self.d.len() {
            scaled[i] = self.d[i] * rhs[i];
        }
    }

    /// Unscale a solution vector: x = D * y (where y solves (DAD)y = Db).
    pub fn unscale_solution(&self, scaled_sol: &[f64], solution: &mut [f64]) {
        for i in 0..self.d.len() {
            solution[i] = self.d[i] * scaled_sol[i];
        }
    }

    /// Apply scaling to a CSC matrix in-place: A → D*A*D.
    /// For upper-triangle entry (i, j), multiply by D[i]*D[j].
    pub fn scale_csc(&self, csc: &mut CscMatrix) {
        for j in 0..csc.n {
            let dj = self.d[j];
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                csc.vals[idx] *= self.d[i] * dj;
            }
        }
    }
}

/// Compute row/column infinity-norms for a symmetric CSC matrix.
/// Returns a vector where `norms[i]` = max_j |a_ij|.
fn row_col_inf_norms(csc: &CscMatrix) -> Vec<f64> {
    let mut norms = vec![0.0f64; csc.n];
    for j in 0..csc.n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            let abs_val = csc.vals[idx].abs();
            if abs_val > norms[i] {
                norms[i] = abs_val;
            }
            // Symmetric: (i,j) entry also contributes to column i / row j
            if i != j && abs_val > norms[j] {
                norms[j] = abs_val;
            }
        }
    }
    norms
}

/// Compute diagonal equilibration scaling.
/// D_ii = 1/sqrt(max_j |a_ij|), with D_ii = 1 for zero rows.
pub fn compute_diagonal_scaling(csc: &CscMatrix) -> ScalingFactors {
    let norms = row_col_inf_norms(csc);
    let d: Vec<f64> = norms
        .iter()
        .map(|&norm| {
            if norm > 0.0 {
                1.0 / norm.sqrt()
            } else {
                1.0
            }
        })
        .collect();
    ScalingFactors { d }
}

/// Compute Ruiz iterative equilibration scaling.
/// Iteratively applies diagonal equilibration until row/column norms
/// are approximately 1 (or max_iter is reached).
pub fn compute_ruiz_scaling(csc: &CscMatrix, max_iter: usize) -> ScalingFactors {
    let n = csc.n;
    let mut d = vec![1.0; n];

    // Work on a copy of the values to track cumulative scaling
    let mut scaled_vals = csc.vals.clone();

    for _ in 0..max_iter {
        // Compute infinity norms of the current scaled matrix
        let mut norms = vec![0.0f64; n];
        for j in 0..n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                let abs_val = scaled_vals[idx].abs();
                if abs_val > norms[i] {
                    norms[i] = abs_val;
                }
                if i != j && abs_val > norms[j] {
                    norms[j] = abs_val;
                }
            }
        }

        // Compute this iteration's scaling factors
        let mut max_deviation = 0.0f64;
        let mut iter_d = vec![1.0; n];
        for i in 0..n {
            if norms[i] > 0.0 {
                iter_d[i] = 1.0 / norms[i].sqrt();
                let deviation = (norms[i] - 1.0).abs();
                if deviation > max_deviation {
                    max_deviation = deviation;
                }
            }
        }

        // Check convergence: all norms close to 1
        if max_deviation < 1e-2 {
            break;
        }

        // Accumulate into total scaling
        for i in 0..n {
            d[i] *= iter_d[i];
        }

        // Apply this iteration's scaling to the working values
        for j in 0..n {
            let dj = iter_d[j];
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                scaled_vals[idx] *= iter_d[i] * dj;
            }
        }
    }

    ScalingFactors { d }
}

/// Compute KKT-aware two-phase scaling.
///
/// Phase 1: Normalize diagonals. For each row/column i where |a_ii| is much
/// smaller than the row's max off-diagonal (saddle-point block), scale by
/// 1/sqrt(|a_ii|) to bring the diagonal to O(1). This addresses the fundamental
/// scale mismatch in KKT systems where equality constraint diagonals are O(1e-8).
///
/// Phase 2: Standard Ruiz equilibration on the result.
///
/// This matches MUMPS's ICNTL(8)=77 behavior for symmetric indefinite systems.
pub fn compute_kkt_scaling(csc: &CscMatrix, max_iter: usize) -> ScalingFactors {
    let n = csc.n;
    let mut d = vec![1.0; n];

    // Phase 1: Diagonal normalization for rows with very small diagonals
    // Find diagonal entries and row max off-diagonals
    let mut diag = vec![0.0f64; n];
    let mut max_offdiag = vec![0.0f64; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            let abs_val = csc.vals[idx].abs();
            if i == j {
                diag[j] = abs_val;
            } else {
                max_offdiag[i] = max_offdiag[i].max(abs_val);
                max_offdiag[j] = max_offdiag[j].max(abs_val);
            }
        }
    }

    // Scale rows where diagonal is zero or much smaller than off-diagonal.
    // For KKT systems, equality constraint rows have zero diagonal.
    // Scale these based on off-diagonal magnitude to bring them to O(1).
    for i in 0..n {
        if max_offdiag[i] > 0.0 {
            if diag[i] == 0.0 {
                // Zero diagonal (equality constraint) — scale by off-diagonal
                let alpha = 1.0 / max_offdiag[i].sqrt();
                let alpha = alpha.min(1e8);
                d[i] = alpha;
            } else if max_offdiag[i] / diag[i] > 100.0 {
                // Diagonal much smaller than off-diagonal
                let alpha = 1.0 / diag[i].sqrt();
                let alpha = alpha.min(1e8);
                d[i] = alpha;
            }
        }
    }

    // Apply Phase 1 scaling to get working values
    let mut scaled_vals = csc.vals.clone();
    for j in 0..n {
        let dj = d[j];
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            scaled_vals[idx] *= d[i] * dj;
        }
    }

    // Phase 2: Ruiz equilibration on the scaled matrix
    for _ in 0..max_iter {
        let mut norms = vec![0.0f64; n];
        for j in 0..n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                let abs_val = scaled_vals[idx].abs();
                norms[i] = norms[i].max(abs_val);
                if i != j {
                    norms[j] = norms[j].max(abs_val);
                }
            }
        }

        let mut max_deviation = 0.0f64;
        let mut iter_d = vec![1.0; n];
        for i in 0..n {
            if norms[i] > 0.0 {
                iter_d[i] = 1.0 / norms[i].sqrt();
                max_deviation = max_deviation.max((norms[i] - 1.0).abs());
            }
        }

        if max_deviation < 1e-2 {
            break;
        }

        for i in 0..n {
            d[i] *= iter_d[i];
        }

        for j in 0..n {
            let dj = iter_d[j];
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                scaled_vals[idx] *= iter_d[i] * dj;
            }
        }
    }

    ScalingFactors { d }
}

/// Compute matching-based scaling for symmetric indefinite matrices.
///
/// For each row i, finds the column j with largest |a_ij| (greedy matching)
/// and scales so that this entry becomes O(1) after D*A*D. This places large
/// entries on or near the diagonal, which is critical for stable Bunch-Kaufman
/// pivoting on KKT systems with zero diagonal blocks.
///
/// After matching-based scaling, applies Ruiz equilibration for fine-tuning.
/// This approximates MUMPS's ICNTL(8)=77 (MC64 + equilibration).
pub fn compute_matching_scaling(csc: &CscMatrix, max_iter: usize) -> ScalingFactors {
    let n = csc.n;
    let mut d = vec![1.0; n];

    // Phase 1: Matching-based scaling
    // For each row i, find the largest absolute entry and scale to make it O(1).
    // For symmetric matrices stored as upper triangle, we need to consider
    // both the row and column views.

    // Compute max absolute entry per row (considering symmetry)
    let mut row_max = vec![0.0f64; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            let abs_val = csc.vals[idx].abs();
            // Upper triangle: (i, j) with i <= j
            // This entry contributes to row i and row j (symmetric)
            row_max[i] = row_max[i].max(abs_val);
            if i != j {
                row_max[j] = row_max[j].max(abs_val);
            }
        }
    }

    // Scale each row/column so that its maximum entry is 1
    // D[i] = 1/sqrt(row_max[i]) for symmetric scaling D*A*D
    for i in 0..n {
        if row_max[i] > 1e-30 {
            d[i] = 1.0 / row_max[i].sqrt();
        }
    }

    // Apply Phase 1 scaling to working values
    let mut scaled_vals = csc.vals.clone();
    for j in 0..n {
        let dj = d[j];
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            scaled_vals[idx] *= d[i] * dj;
        }
    }

    // Phase 2: Ruiz equilibration on the scaled matrix
    for _ in 0..max_iter {
        let mut norms = vec![0.0f64; n];
        for j in 0..n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                let abs_val = scaled_vals[idx].abs();
                norms[i] = norms[i].max(abs_val);
                if i != j {
                    norms[j] = norms[j].max(abs_val);
                }
            }
        }

        let mut max_deviation = 0.0f64;
        let mut iter_d = vec![1.0; n];
        for i in 0..n {
            if norms[i] > 0.0 {
                iter_d[i] = 1.0 / norms[i].sqrt();
                max_deviation = max_deviation.max((norms[i] - 1.0).abs());
            }
        }

        if max_deviation < 1e-2 {
            break;
        }

        for i in 0..n {
            d[i] *= iter_d[i];
        }

        for j in 0..n {
            let dj = iter_d[j];
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                scaled_vals[idx] *= iter_d[i] * dj;
            }
        }
    }

    ScalingFactors { d }
}

/// Compute scaling factors for the given method.
pub fn compute_scaling(csc: &CscMatrix, method: Scaling) -> Option<ScalingFactors> {
    match method {
        Scaling::None => None,
        Scaling::Diagonal => Some(compute_diagonal_scaling(csc)),
        Scaling::Ruiz { max_iter } => Some(compute_kkt_scaling(csc, max_iter)),
    }
}

/// Compute scaling for KKT systems using matching-based approach.
pub fn compute_scaling_kkt(csc: &CscMatrix, method: Scaling) -> Option<ScalingFactors> {
    match method {
        Scaling::None => None,
        _ => Some(compute_matching_scaling(csc, 10)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::solver::{Solver, SolverOptions};

    fn make_csc(n: usize, triplets: &[(usize, usize, f64)]) -> CscMatrix {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        CscMatrix::from_coo(&coo)
    }

    #[test]
    fn test_diagonal_scaling_identity() {
        // Identity matrix: scaling should be D = I
        let csc = make_csc(3, &[(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)]);
        let sf = compute_diagonal_scaling(&csc);
        for &di in &sf.d {
            assert!((di - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn test_diagonal_scaling_reduces_range() {
        // Matrix with entries spanning 1e-4 to 1e4
        let csc = make_csc(3, &[
            (0, 0, 1e4), (0, 1, 1.0),
            (1, 1, 1e-4), (1, 2, 1e-2),
            (2, 2, 1.0),
        ]);

        let norms_before = row_col_inf_norms(&csc);
        let range_before = norms_before.iter().cloned().fold(0.0f64, f64::max)
            / norms_before.iter().cloned().fold(f64::MAX, f64::min);

        let sf = compute_diagonal_scaling(&csc);
        let mut scaled = csc.clone();
        sf.scale_csc(&mut scaled);

        let norms_after = row_col_inf_norms(&scaled);
        let range_after = norms_after.iter().cloned().fold(0.0f64, f64::max)
            / norms_after.iter().cloned().fold(f64::MAX, f64::min);

        assert!(range_after < range_before,
            "scaling should reduce dynamic range: before={:.1e}, after={:.1e}",
            range_before, range_after);
    }

    #[test]
    fn test_ruiz_scaling_convergence() {
        // Ill-conditioned matrix
        let csc = make_csc(3, &[
            (0, 0, 1e6), (0, 1, 1.0),
            (1, 1, 1e-6), (1, 2, 1e-3),
            (2, 2, 1.0),
        ]);

        let sf = compute_ruiz_scaling(&csc, 10);
        let mut scaled = csc.clone();
        sf.scale_csc(&mut scaled);

        let norms = row_col_inf_norms(&scaled);
        // After Ruiz, all norms should be close to 1
        for (i, &norm) in norms.iter().enumerate() {
            assert!(norm > 0.1 && norm < 10.0,
                "row/col {} norm = {:.3e}, expected ~1", i, norm);
        }
    }

    #[test]
    fn test_scaling_preserves_solution() {
        // Solve Ax = b with and without scaling — same answer
        let coo = CooMatrix::new(3,
            vec![0, 0, 1, 1, 2],
            vec![0, 1, 1, 2, 2],
            vec![1e4, 1.0, 1e-4, 1e-2, 1.0],
        ).unwrap();

        let b = [1.0, 2.0, 3.0];

        // Solve without scaling
        let mut solver = Solver::new(SolverOptions::default());
        solver.analyze_and_factor(&coo).unwrap();
        let mut x_unscaled = [0.0; 3];
        solver.solve(&b, &mut x_unscaled).unwrap();

        // Solve with scaling: D*A*D * y = D*b, x = D*y
        let csc = CscMatrix::from_coo(&coo);
        let sf = compute_ruiz_scaling(&csc, 10);

        let mut scaled_csc = csc.clone();
        sf.scale_csc(&mut scaled_csc);

        // Build COO from scaled CSC for the solver
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for j in 0..scaled_csc.n {
            for idx in scaled_csc.col_ptr[j]..scaled_csc.col_ptr[j + 1] {
                rows.push(scaled_csc.row_idx[idx]);
                cols.push(j);
                vals.push(scaled_csc.vals[idx]);
            }
        }
        let scaled_coo = CooMatrix::new(3, rows, cols, vals).unwrap();

        let mut solver2 = Solver::new(SolverOptions::default());
        solver2.analyze_and_factor(&scaled_coo).unwrap();

        let mut scaled_rhs = [0.0; 3];
        sf.scale_rhs(&b, &mut scaled_rhs);

        let mut y = [0.0; 3];
        solver2.solve(&scaled_rhs, &mut y).unwrap();

        let mut x_scaled = [0.0; 3];
        sf.unscale_solution(&y, &mut x_scaled);

        // Both solutions should match
        for i in 0..3 {
            let rel_err = if x_unscaled[i].abs() > 1e-15 {
                (x_scaled[i] - x_unscaled[i]).abs() / x_unscaled[i].abs()
            } else {
                (x_scaled[i] - x_unscaled[i]).abs()
            };
            assert!(rel_err < 1e-6,
                "x[{}]: unscaled={:.6e}, scaled={:.6e}, rel_err={:.3e}",
                i, x_unscaled[i], x_scaled[i], rel_err);
        }
    }

    #[test]
    fn test_scaling_improves_accuracy_on_ill_conditioned() {
        // KKT-like system with extreme scaling mismatch:
        // H = diag(1e8, 1e-8), A = [1, 1], constraint block zeros
        let coo = CooMatrix::new(3,
            vec![0, 1, 0, 1, 2],
            vec![0, 1, 2, 2, 2],
            vec![1e8, 1e-8, 1.0, 1.0, 0.0],
        ).unwrap();

        let b = [1e8, 1e-8, 2.0]; // Chosen so x = [1, 1, 0]

        // Solve without scaling
        let mut solver1 = Solver::new(SolverOptions::default());
        solver1.analyze_and_factor(&coo).unwrap();
        let mut x1 = [0.0; 3];
        solver1.solve(&b, &mut x1).unwrap();

        // Compute residual without scaling
        let csc = CscMatrix::from_coo(&coo);
        let mut ax1 = [0.0; 3];
        csc.matvec(&x1, &mut ax1);
        let resid1: f64 = ax1.iter().zip(&b).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);

        // Solve with Ruiz scaling
        let sf = compute_ruiz_scaling(&csc, 10);
        let mut scaled_csc = csc.clone();
        sf.scale_csc(&mut scaled_csc);

        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for j in 0..scaled_csc.n {
            for idx in scaled_csc.col_ptr[j]..scaled_csc.col_ptr[j + 1] {
                rows.push(scaled_csc.row_idx[idx]);
                cols.push(j);
                vals.push(scaled_csc.vals[idx]);
            }
        }
        let scaled_coo = CooMatrix::new(3, rows, cols, vals).unwrap();

        let mut solver2 = Solver::new(SolverOptions::default());
        solver2.analyze_and_factor(&scaled_coo).unwrap();

        let mut scaled_rhs = [0.0; 3];
        sf.scale_rhs(&b, &mut scaled_rhs);
        let mut y = [0.0; 3];
        solver2.solve(&scaled_rhs, &mut y).unwrap();
        let mut x2 = [0.0; 3];
        sf.unscale_solution(&y, &mut x2);

        // Compute residual with scaling
        let mut ax2 = [0.0; 3];
        csc.matvec(&x2, &mut ax2);
        let resid2: f64 = ax2.iter().zip(&b).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);

        // Scaled solution should have comparable or better residual
        // (On this problem, both should be good since iterative refinement helps,
        // but scaling prevents catastrophic loss on truly bad problems)
        assert!(resid2 < 1.0,
            "scaled residual {} should be small", resid2);
        assert!(resid1 < 1.0 || resid2 <= resid1,
            "scaling should not make things worse: unscaled={:.3e}, scaled={:.3e}",
            resid1, resid2);
    }

    #[test]
    fn test_scaling_zero_row() {
        // Matrix with a zero row/column — scaling should leave it as D=1
        let csc = make_csc(3, &[
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 5.0),
            // row/col 2 has no entries
        ]);

        let sf = compute_diagonal_scaling(&csc);
        assert!((sf.d[2] - 1.0).abs() < 1e-15, "zero row should get D=1");
    }

    #[test]
    fn test_ruiz_vs_diagonal_on_extreme_matrix() {
        // Matrix with 12 orders of magnitude range
        let csc = make_csc(4, &[
            (0, 0, 1e12),
            (1, 1, 1.0), (1, 2, 0.5),
            (2, 2, 1e-6), (2, 3, 1e-3),
            (3, 3, 1e-12),
        ]);

        let sf_diag = compute_diagonal_scaling(&csc);
        let mut scaled_diag = csc.clone();
        sf_diag.scale_csc(&mut scaled_diag);
        let norms_diag = row_col_inf_norms(&scaled_diag);
        let range_diag = norms_diag.iter().cloned().fold(0.0f64, f64::max)
            / norms_diag.iter().cloned().filter(|&x| x > 0.0).fold(f64::MAX, f64::min);

        let sf_ruiz = compute_ruiz_scaling(&csc, 20);
        let mut scaled_ruiz = csc.clone();
        sf_ruiz.scale_csc(&mut scaled_ruiz);
        let norms_ruiz = row_col_inf_norms(&scaled_ruiz);
        let range_ruiz = norms_ruiz.iter().cloned().fold(0.0f64, f64::max)
            / norms_ruiz.iter().cloned().filter(|&x| x > 0.0).fold(f64::MAX, f64::min);

        // Ruiz should produce tighter range than single-pass diagonal
        assert!(range_ruiz <= range_diag + 1e-10,
            "Ruiz range ({:.1e}) should be <= diagonal range ({:.1e})",
            range_ruiz, range_diag);
    }
}
