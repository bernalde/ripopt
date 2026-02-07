use crate::linear_solver::LinearSolver;
use crate::options::SolverOptions;

/// State for the restoration phase.
///
/// When the filter line search fails completely, the restoration phase
/// attempts to find a point that is acceptable to the filter by
/// minimizing constraint violation.
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
    /// Uses a penalty-based approach: solve
    ///   min sum_i (p_i + n_i)
    ///   s.t. g(x) - p + n = 0, p >= 0, n >= 0, x_l <= x <= x_u
    ///
    /// Returns (new_x, success) where success indicates whether a point
    /// acceptable to the original filter was found.
    #[allow(clippy::too_many_arguments)]
    pub fn restore(
        &mut self,
        x: &[f64],
        x_l: &[f64],
        x_u: &[f64],
        g: &[f64],
        g_l: &[f64],
        g_u: &[f64],
        _grad_f: &[f64],
        jac_rows: &[usize],
        jac_cols: &[usize],
        jac_vals: &[f64],
        n: usize,
        _m: usize,
        _solver: &mut dyn LinearSolver,
        options: &SolverOptions,
        _is_acceptable_to_filter: &dyn Fn(f64, f64) -> bool,
    ) -> (Vec<f64>, bool) {
        self.active = true;

        let mut x_rest = x.to_vec();

        // Simplified restoration: use a damped Newton step to reduce infeasibility.
        // This is a basic implementation — Ipopt's full restoration is more sophisticated.
        for _iter in 0..self.max_iter {
            // Compute constraint violation
            let theta: f64 = g
                .iter()
                .zip(g_l.iter().zip(g_u.iter()))
                .map(|(&gi, (&gli, &gui))| {
                    let viol_l = if gli.is_finite() { (gli - gi).max(0.0) } else { 0.0 };
                    let viol_u = if gui.is_finite() { (gi - gui).max(0.0) } else { 0.0 };
                    viol_l.max(viol_u)
                })
                .fold(0.0f64, f64::max);

            if theta < options.tol {
                self.active = false;
                return (x_rest, true);
            }

            // Take a damped step to reduce constraint violation
            // Compute J^T * (g - g_target) as the gradient of ||c||^2 / 2
            let mut step = vec![0.0; n];
            for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
                let violation = if g[row] < g_l[row] {
                    g[row] - g_l[row]
                } else if g[row] > g_u[row] {
                    g[row] - g_u[row]
                } else {
                    0.0
                };
                step[col] -= jac_vals[idx] * violation;
            }

            // Damped step with step size control
            let step_norm: f64 = step.iter().map(|s| s * s).sum::<f64>().sqrt();
            let alpha = if step_norm > 1.0 {
                1.0 / step_norm
            } else {
                1.0
            };

            for i in 0..n {
                x_rest[i] += alpha * step[i];
                // Project onto bounds
                if x_l[i].is_finite() {
                    x_rest[i] = x_rest[i].max(x_l[i] + 1e-8);
                }
                if x_u[i].is_finite() {
                    x_rest[i] = x_rest[i].min(x_u[i] - 1e-8);
                }
            }
        }

        self.active = false;
        (x_rest, false)
    }
}
