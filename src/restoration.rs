use crate::convergence;
use crate::linear_solver::dense::DenseLdl;
use crate::linear_solver::{LinearSolver, SymmetricMatrix};
use crate::options::SolverOptions;

/// State for the restoration phase.
///
/// When the filter line search fails completely, the restoration phase
/// attempts to find a point that is acceptable to the filter by
/// minimizing constraint violation using Gauss-Newton steps on ||violation||^2.
pub struct RestorationPhase {
    /// Maximum iterations in restoration.
    max_iter: usize,
    /// Whether restoration is currently active.
    active: bool,
}

impl RestorationPhase {
    pub fn new(max_iter: usize) -> Self {
        Self {
            max_iter,
            active: false,
        }
    }

    pub fn is_active(&self) -> bool {
        self.active
    }

    /// Attempt restoration: minimize constraint violation subject to bounds.
    ///
    /// Uses Gauss-Newton steps on 0.5*||violation||^2, which provides quadratic
    /// convergence for nonlinear equality constraints (unlike gradient descent
    /// which converges linearly). Falls back to gradient descent if the
    /// Gauss-Newton system is singular.
    ///
    /// Returns (new_x, success) where success indicates whether a point
    /// with sufficiently small constraint violation was found.
    #[allow(clippy::too_many_arguments)]
    pub fn restore(
        &mut self,
        x: &[f64],
        x_l: &[f64],
        x_u: &[f64],
        g_l: &[f64],
        g_u: &[f64],
        jac_rows: &[usize],
        jac_cols: &[usize],
        n: usize,
        m: usize,
        options: &SolverOptions,
        _is_acceptable_to_filter: &dyn Fn(f64, f64) -> bool,
        eval_constraints: &dyn Fn(&[f64], &mut [f64]),
        eval_jacobian: &dyn Fn(&[f64], &mut [f64]),
    ) -> (Vec<f64>, bool) {
        self.active = true;

        if m == 0 {
            // No constraints: nothing to restore.
            self.active = false;
            return (x.to_vec(), true);
        }

        let mut x_rest = x.to_vec();
        let mut g = vec![0.0; m];
        let jac_nnz = jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];

        for _iter in 0..self.max_iter {
            // Evaluate constraints at current point
            eval_constraints(&x_rest, &mut g);

            // Compute constraint violation
            let theta = convergence::primal_infeasibility(&g, g_l, g_u);

            if theta < options.tol {
                self.active = false;
                return (x_rest, true);
            }

            // Evaluate Jacobian at current point
            eval_jacobian(&x_rest, &mut jac_vals);

            // Compute violation vector for each constraint
            let mut violation = vec![0.0; m];
            let mut active = vec![false; m];
            for i in 0..m {
                let is_equality = g_l[i].is_finite()
                    && g_u[i].is_finite()
                    && (g_l[i] - g_u[i]).abs() < 1e-15;
                if is_equality {
                    violation[i] = g[i] - g_l[i];
                    active[i] = true;
                } else if g_l[i].is_finite() && g[i] < g_l[i] {
                    violation[i] = g[i] - g_l[i];
                    active[i] = true;
                } else if g_u[i].is_finite() && g[i] > g_u[i] {
                    violation[i] = g[i] - g_u[i];
                    active[i] = true;
                }
            }

            // Collect active constraint indices
            let active_indices: Vec<usize> = (0..m).filter(|&i| active[i]).collect();
            let m_active = active_indices.len();

            if m_active == 0 {
                break;
            }

            // Try Gauss-Newton step: dx = -J_a^T * (J_a * J_a^T + eps*I)^{-1} * v_a
            let step = self.gauss_newton_step(
                &jac_rows,
                &jac_cols,
                &jac_vals,
                &violation,
                &active_indices,
                n,
            );

            let step = match step {
                Some(s) => s,
                None => {
                    // Fall back to gradient descent: step = -J^T * violation
                    let mut grad_step = vec![0.0; n];
                    for (idx, (&row, &col)) in
                        jac_rows.iter().zip(jac_cols.iter()).enumerate()
                    {
                        if active[row] {
                            grad_step[col] -= jac_vals[idx] * violation[row];
                        }
                    }
                    grad_step
                }
            };

            // Normalize step if too large
            let step_norm: f64 = step.iter().map(|s| s * s).sum::<f64>().sqrt();
            if step_norm < 1e-20 {
                break;
            }
            let scale = if step_norm > 10.0 {
                10.0 / step_norm
            } else {
                1.0
            };

            // Backtracking line search on theta
            let mut alpha = scale;
            let mut x_trial = vec![0.0; n];
            let mut g_trial = vec![0.0; m];
            let mut found_decrease = false;

            for _ls in 0..30 {
                for i in 0..n {
                    x_trial[i] = x_rest[i] + alpha * step[i];
                    if x_l[i].is_finite() {
                        x_trial[i] = x_trial[i].max(x_l[i] + 1e-8);
                    }
                    if x_u[i].is_finite() {
                        x_trial[i] = x_trial[i].min(x_u[i] - 1e-8);
                    }
                }

                eval_constraints(&x_trial, &mut g_trial);
                let theta_trial = convergence::primal_infeasibility(&g_trial, g_l, g_u);

                if theta_trial < (1.0 - 1e-4 * alpha) * theta {
                    x_rest.copy_from_slice(&x_trial);
                    found_decrease = true;
                    break;
                }

                alpha *= 0.5;
            }

            if !found_decrease {
                // Gauss-Newton line search failed — try gradient descent as fallback
                // step_gd = -J_a^T * violation_a (steepest descent on 0.5*||violation||^2)
                let mut grad_step = vec![0.0; n];
                for (idx, (&row, &col)) in
                    jac_rows.iter().zip(jac_cols.iter()).enumerate()
                {
                    if active[row] {
                        grad_step[col] -= jac_vals[idx] * violation[row];
                    }
                }

                let gd_norm: f64 = grad_step.iter().map(|s| s * s).sum::<f64>().sqrt();
                if gd_norm < 1e-20 {
                    break;
                }
                let gd_scale = if gd_norm > 10.0 {
                    10.0 / gd_norm
                } else {
                    1.0
                };

                let mut gd_alpha = gd_scale;
                let mut gd_found = false;
                for _ls in 0..30 {
                    for i in 0..n {
                        x_trial[i] = x_rest[i] + gd_alpha * grad_step[i];
                        if x_l[i].is_finite() {
                            x_trial[i] = x_trial[i].max(x_l[i] + 1e-8);
                        }
                        if x_u[i].is_finite() {
                            x_trial[i] = x_trial[i].min(x_u[i] - 1e-8);
                        }
                    }

                    eval_constraints(&x_trial, &mut g_trial);
                    let theta_trial = convergence::primal_infeasibility(&g_trial, g_l, g_u);

                    if theta_trial < (1.0 - 1e-4 * gd_alpha) * theta {
                        x_rest.copy_from_slice(&x_trial);
                        gd_found = true;
                        break;
                    }

                    gd_alpha *= 0.5;
                }

                if !gd_found {
                    break;
                }
            }
        }

        // Check final constraint violation
        eval_constraints(&x_rest, &mut g);
        let theta_final = convergence::primal_infeasibility(&g, g_l, g_u);

        self.active = false;
        if theta_final < options.tol {
            (x_rest, true)
        } else {
            // Return the improved point even if not fully feasible,
            // as long as it's different from the input
            let moved: f64 = x_rest
                .iter()
                .zip(x.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            (x_rest, moved > 1e-14)
        }
    }

    /// Compute Gauss-Newton step: dx = -J_a^T * (J_a * J_a^T + eps*I)^{-1} * v_a
    ///
    /// where J_a is the Jacobian restricted to active (violated) constraints
    /// and v_a is the violation vector for active constraints.
    fn gauss_newton_step(
        &self,
        jac_rows: &[usize],
        jac_cols: &[usize],
        jac_vals: &[f64],
        violation: &[f64],
        active_indices: &[usize],
        n: usize,
    ) -> Option<Vec<f64>> {
        let m_active = active_indices.len();
        if m_active == 0 {
            return None;
        }

        // Map from original constraint index to active index
        let mut active_map = vec![usize::MAX; violation.len()];
        for (ai, &orig) in active_indices.iter().enumerate() {
            active_map[orig] = ai;
        }

        // Form J_a * J_a^T (m_active x m_active)
        // Group Jacobian entries by column for efficient J*J^T computation
        let mut col_entries: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
        for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
            if active_map[row] != usize::MAX {
                col_entries[col].push((active_map[row], jac_vals[idx]));
            }
        }

        let mut jjt = SymmetricMatrix::zeros(m_active);
        for col_ents in &col_entries {
            for &(ai, val_i) in col_ents {
                for &(aj, val_j) in col_ents {
                    if ai >= aj {
                        jjt.add(ai, aj, val_i * val_j);
                    }
                }
            }
        }

        // Add Levenberg-Marquardt regularization for numerical stability
        let jjt_diag_max = (0..m_active)
            .map(|i| jjt.get(i, i).abs())
            .fold(0.0f64, f64::max);
        let eps = 1e-8 * jjt_diag_max.max(1.0);
        jjt.add_diagonal(eps);

        // Solve (J_a * J_a^T + eps*I) * w = violation_a
        let mut solver = DenseLdl::new();
        if solver.factor(&jjt).is_err() {
            return None;
        }

        let v_active: Vec<f64> = active_indices.iter().map(|&i| violation[i]).collect();
        let mut w = vec![0.0; m_active];
        if solver.solve(&v_active, &mut w).is_err() {
            return None;
        }

        // step = -J_a^T * w
        let mut step = vec![0.0; n];
        for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
            let ai = active_map[row];
            if ai != usize::MAX {
                step[col] -= jac_vals[idx] * w[ai];
            }
        }

        Some(step)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::options::SolverOptions;

    fn default_opts() -> SolverOptions {
        SolverOptions {
            print_level: 0,
            ..SolverOptions::default()
        }
    }

    #[test]
    fn test_restoration_no_constraints() {
        let mut phase = RestorationPhase::new(50);
        let x = vec![1.0, 2.0];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let opts = default_opts();

        let (x_new, success) = phase.restore(
            &x, &x_l, &x_u, &[], &[],
            &[], &[], 2, 0, &opts,
            &|_theta, _phi| true,
            &|_x, _g| {},
            &|_x, _vals| {},
        );

        assert!(success, "No constraints → immediate success");
        assert!((x_new[0] - 1.0).abs() < 1e-15);
        assert!((x_new[1] - 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_restoration_already_feasible() {
        let mut phase = RestorationPhase::new(50);
        let x = vec![0.5, 0.5];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let g_l = vec![1.0]; // g = x0 + x1 = 1.0 = g_l → feasible
        let g_u = vec![1.0];
        let jac_rows = vec![0, 0];
        let jac_cols = vec![0, 1];
        let opts = default_opts();

        let (_, success) = phase.restore(
            &x, &x_l, &x_u, &g_l, &g_u,
            &jac_rows, &jac_cols, 2, 1, &opts,
            &|_theta, _phi| true,
            &|x, g| { g[0] = x[0] + x[1]; },
            &|_x, vals| { vals[0] = 1.0; vals[1] = 1.0; },
        );

        assert!(success, "Already feasible → success");
    }

    #[test]
    fn test_restoration_linear_equality() {
        // g(x) = x0 + x1 = 1, from x = (0, 0)
        let mut phase = RestorationPhase::new(100);
        let x = vec![0.0, 0.0];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let g_l = vec![1.0];
        let g_u = vec![1.0];
        let jac_rows = vec![0, 0];
        let jac_cols = vec![0, 1];
        let opts = default_opts();

        let (x_new, success) = phase.restore(
            &x, &x_l, &x_u, &g_l, &g_u,
            &jac_rows, &jac_cols, 2, 1, &opts,
            &|_theta, _phi| true,
            &|x, g| { g[0] = x[0] + x[1]; },
            &|_x, vals| { vals[0] = 1.0; vals[1] = 1.0; },
        );

        assert!(success, "Linear equality should be restored");
        let g_val = x_new[0] + x_new[1];
        assert!((g_val - 1.0).abs() < 1e-6,
            "Constraint should be satisfied: g = {}", g_val);
    }

    #[test]
    fn test_restoration_linear_inequality() {
        // g(x) = x0 >= 2, from x = (0)
        let mut phase = RestorationPhase::new(100);
        let x = vec![0.0];
        let x_l = vec![f64::NEG_INFINITY];
        let x_u = vec![f64::INFINITY];
        let g_l = vec![2.0];
        let g_u = vec![f64::INFINITY];
        let jac_rows = vec![0];
        let jac_cols = vec![0];
        let opts = default_opts();

        let (x_new, success) = phase.restore(
            &x, &x_l, &x_u, &g_l, &g_u,
            &jac_rows, &jac_cols, 1, 1, &opts,
            &|_theta, _phi| true,
            &|x, g| { g[0] = x[0]; },
            &|_x, vals| { vals[0] = 1.0; },
        );

        assert!(success, "Linear inequality should be restored");
        assert!(x_new[0] >= 2.0 - 1e-6,
            "Should satisfy x >= 2, got {}", x_new[0]);
    }

    #[test]
    fn test_restoration_with_bounds() {
        // g(x) = x0 + x1 = 5, from x = (0, 0), with bounds 0 <= xi <= 3
        let mut phase = RestorationPhase::new(100);
        let x = vec![0.5, 0.5];
        let x_l = vec![0.0, 0.0];
        let x_u = vec![3.0, 3.0];
        let g_l = vec![5.0];
        let g_u = vec![5.0];
        let jac_rows = vec![0, 0];
        let jac_cols = vec![0, 1];
        let opts = default_opts();

        let (x_new, _success) = phase.restore(
            &x, &x_l, &x_u, &g_l, &g_u,
            &jac_rows, &jac_cols, 2, 1, &opts,
            &|_theta, _phi| true,
            &|x, g| { g[0] = x[0] + x[1]; },
            &|_x, vals| { vals[0] = 1.0; vals[1] = 1.0; },
        );

        // Verify bounds are respected
        for i in 0..2 {
            assert!(x_new[i] >= x_l[i], "x[{}] = {} below lower bound", i, x_new[i]);
            assert!(x_new[i] <= x_u[i], "x[{}] = {} above upper bound", i, x_new[i]);
        }
    }

    #[test]
    fn test_restoration_active_flag() {
        let mut phase = RestorationPhase::new(50);
        assert!(!phase.is_active(), "Should not be active initially");

        let x = vec![1.0];
        let x_l = vec![f64::NEG_INFINITY];
        let x_u = vec![f64::INFINITY];
        let opts = default_opts();

        phase.restore(
            &x, &x_l, &x_u, &[], &[],
            &[], &[], 1, 0, &opts,
            &|_theta, _phi| true,
            &|_x, _g| {},
            &|_x, _vals| {},
        );

        assert!(!phase.is_active(), "Should not be active after restore completes");
    }
}
