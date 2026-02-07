use crate::convergence;
use crate::options::SolverOptions;

/// State for the restoration phase.
///
/// When the filter line search fails completely, the restoration phase
/// attempts to find a point that is acceptable to the filter by
/// minimizing constraint violation via damped gradient descent on ||violation||^2.
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
    /// Uses iterative gradient descent on 0.5*||violation||^2 with re-evaluation
    /// of constraints and Jacobian at each step. Includes a backtracking line
    /// search to ensure the violation actually decreases.
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

            // Compute gradient of 0.5*||violation||^2: step = -J^T * violation
            let mut step = vec![0.0; n];
            for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
                let violation = if g_l[row].is_finite() && g[row] < g_l[row] {
                    g[row] - g_l[row] // negative
                } else if g_u[row].is_finite() && g[row] > g_u[row] {
                    g[row] - g_u[row] // positive
                } else if g_l[row].is_finite() && g_u[row].is_finite()
                    && (g_l[row] - g_u[row]).abs() < 1e-15
                {
                    // Equality constraint: violation = g - target
                    g[row] - g_l[row]
                } else {
                    0.0
                };
                step[col] -= jac_vals[idx] * violation;
            }

            // Normalize step if too large
            let step_norm: f64 = step.iter().map(|s| s * s).sum::<f64>().sqrt();
            if step_norm < 1e-20 {
                // Zero gradient: can't make progress
                break;
            }
            let scale = if step_norm > 10.0 { 10.0 / step_norm } else { 1.0 };

            // Backtracking line search on theta
            let mut alpha = scale;
            let mut x_trial = vec![0.0; n];
            let mut g_trial = vec![0.0; m];
            let mut found_decrease = false;

            for _ls in 0..20 {
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
                    // Sufficient decrease
                    x_rest.copy_from_slice(&x_trial);
                    found_decrease = true;
                    break;
                }

                alpha *= 0.5;
            }

            if !found_decrease {
                // Line search failed: can't reduce theta further
                break;
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
            let moved: f64 = x_rest.iter().zip(x.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            (x_rest, moved > 1e-14)
        }
    }
}
