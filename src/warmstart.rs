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

        // Ensure multipliers are positive before computing complementarity
        for i in 0..n {
            z_l[i] = z_l[i].max(mult_push);
            z_u[i] = z_u[i].max(mult_push);
        }

        // Compute initial mu from average complementarity at the original warm-start point,
        // BEFORE pushing x away from bounds. This gives a mu that reflects the actual
        // warm-start state rather than the artificially modified point.
        let mut sum_compl = 0.0;
        let mut count = 0;
        for i in 0..n {
            if x_l[i].is_finite() {
                let slack = (x[i] - x_l[i]).max(1e-20);
                sum_compl += slack * z_l[i];
                count += 1;
            }
            if x_u[i].is_finite() {
                let slack = (x_u[i] - x[i]).max(1e-20);
                sum_compl += slack * z_u[i];
                count += 1;
            }
        }

        // Now push x away from bounds
        for i in 0..n {
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
        }

        if count > 0 {
            (sum_compl / count as f64).max(options.mu_min)
        } else {
            options.mu_init
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_warmstart_no_bounds() {
        let mut x = vec![1.0, 2.0];
        let mut z_l = vec![0.5, 0.5];
        let mut z_u = vec![0.5, 0.5];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let opts = SolverOptions::default();

        let mu = WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        // x unchanged (no bounds to push from)
        assert!((x[0] - 1.0).abs() < 1e-15);
        assert!((x[1] - 2.0).abs() < 1e-15);
        // z floored to mult_push
        assert!(z_l[0] >= opts.warm_start_mult_bound_push);
        assert!(z_u[0] >= opts.warm_start_mult_bound_push);
        // No finite bounds → mu = mu_init
        assert!((mu - opts.mu_init).abs() < 1e-12);
    }

    #[test]
    fn test_warmstart_lower_bound_push() {
        let mut x = vec![0.0]; // At lower bound
        let mut z_l = vec![1.0];
        let mut z_u = vec![0.0];
        let x_l = vec![0.0];
        let x_u = vec![f64::INFINITY];
        let opts = SolverOptions::default();

        WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        // x should be pushed away from lower bound
        assert!(x[0] > x_l[0], "x should be pushed from lower bound");
        assert!(x[0] >= x_l[0] + opts.warm_start_bound_push);
    }

    #[test]
    fn test_warmstart_both_bounds_push() {
        let mut x = vec![0.0]; // At lower bound
        let mut z_l = vec![1.0];
        let mut z_u = vec![1.0];
        let x_l = vec![0.0];
        let x_u = vec![10.0];
        let opts = SolverOptions::default();

        WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        assert!(x[0] > x_l[0], "x pushed from lower bound");
        assert!(x[0] < x_u[0], "x stays below upper bound");
    }

    #[test]
    fn test_warmstart_narrow_bounds() {
        // Very narrow range: [0, 0.001]
        let mut x = vec![0.0];
        let mut z_l = vec![1.0];
        let mut z_u = vec![1.0];
        let x_l = vec![0.0];
        let x_u = vec![0.001];
        let opts = SolverOptions::default();

        WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        // Should be strictly between bounds
        assert!(x[0] > x_l[0]);
        assert!(x[0] < x_u[0]);
    }

    #[test]
    fn test_warmstart_multiplier_floor() {
        let mut x = vec![5.0];
        let mut z_l = vec![-1.0]; // Negative — should be floored
        let mut z_u = vec![-2.0]; // Negative — should be floored
        let x_l = vec![0.0];
        let x_u = vec![10.0];
        let opts = SolverOptions::default();

        WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        assert!(z_l[0] >= opts.warm_start_mult_bound_push,
            "z_l should be floored, got {}", z_l[0]);
        assert!(z_u[0] >= opts.warm_start_mult_bound_push,
            "z_u should be floored, got {}", z_u[0]);
    }

    #[test]
    fn test_warmstart_mu_from_complementarity() {
        let mut x = vec![1.0, 3.0];
        let mut z_l = vec![2.0, 1.0];
        let mut z_u = vec![0.5, 0.5];
        let x_l = vec![0.0, 0.0];
        let x_u = vec![5.0, 5.0];
        let opts = SolverOptions::default();

        let mu = WarmStartInitializer::initialize(&mut x, &mut z_l, &mut z_u, &x_l, &x_u, &opts);

        // mu should be computed from complementarity products
        // After push adjustments, mu = average of (x-x_l)*z_l + (x_u-x)*z_u
        assert!(mu > 0.0, "mu should be positive");
        assert!(mu >= opts.mu_min, "mu should be at least mu_min");
    }
}
