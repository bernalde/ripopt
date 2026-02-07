use crate::options::SolverOptions;

/// Initialize a warm-start iterate from a previous solution.
///
/// Adjusts the primal variables to be strictly interior to bounds,
/// and initializes multipliers based on the warm-start point.
pub struct WarmStartInitializer;

impl WarmStartInitializer {
    /// Adjust a warm-start point to be feasible for the barrier problem.
    ///
    /// - Push x away from bounds by warm_start_bound_push
    /// - Ensure z_l, z_u are positive (at least warm_start_mult_bound_push)
    /// - Compute initial mu from complementarity of warm-start point
    #[allow(clippy::too_many_arguments)]
    pub fn initialize(
        x: &mut [f64],
        z_l: &mut [f64],
        z_u: &mut [f64],
        x_l: &[f64],
        x_u: &[f64],
        options: &SolverOptions,
    ) -> f64 {
        let n = x.len();
        let push = options.warm_start_bound_push;
        let frac = options.warm_start_bound_frac;
        let mult_push = options.warm_start_mult_bound_push;

        for i in 0..n {
            // Push x away from bounds
            if x_l[i].is_finite() {
                let bound_dist = if x_u[i].is_finite() {
                    push.min(frac * (x_u[i] - x_l[i]))
                } else {
                    push
                };
                x[i] = x[i].max(x_l[i] + bound_dist);
            }
            if x_u[i].is_finite() {
                let bound_dist = if x_l[i].is_finite() {
                    push.min(frac * (x_u[i] - x_l[i]))
                } else {
                    push
                };
                x[i] = x[i].min(x_u[i] - bound_dist);
            }

            // Ensure multipliers are positive
            z_l[i] = z_l[i].max(mult_push);
            z_u[i] = z_u[i].max(mult_push);
        }

        // Compute initial mu from average complementarity
        let mut sum_compl = 0.0;
        let mut count = 0;
        for i in 0..n {
            if x_l[i].is_finite() {
                sum_compl += (x[i] - x_l[i]) * z_l[i];
                count += 1;
            }
            if x_u[i].is_finite() {
                sum_compl += (x_u[i] - x[i]) * z_u[i];
                count += 1;
            }
        }

        if count > 0 {
            (sum_compl / count as f64).max(options.mu_min)
        } else {
            options.mu_init
        }
    }
}
