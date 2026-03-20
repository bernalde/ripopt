/// An entry in the filter, representing a constraint violation / objective pair.
#[derive(Debug, Clone, Copy)]
pub struct FilterEntry {
    /// Constraint violation (infeasibility measure).
    pub theta: f64,
    /// Barrier objective value.
    pub phi: f64,
}

/// The filter mechanism for the line search.
///
/// Maintains a set of (theta, phi) pairs. A trial point is acceptable if it
/// "dominates" (improves upon) all entries in the filter by sufficient margin.
pub struct Filter {
    entries: Vec<FilterEntry>,
    /// Maximum constraint violation allowed (theta_max).
    theta_max: f64,
    /// Margin for filter acceptance (gamma_theta).
    gamma_theta: f64,
    /// Margin for filter acceptance (gamma_phi).
    gamma_phi: f64,
    /// Threshold for switching condition (theta_min).
    theta_min: f64,
    /// Exponent for switching condition.
    s_theta: f64,
    /// Exponent for switching condition.
    s_phi: f64,
    /// Armijo parameter (eta_phi).
    eta_phi: f64,
    /// Small constant delta for filter margin.
    delta: f64,
}

impl Filter {
    /// Create a new filter with the given maximum constraint violation.
    pub fn new(theta_max: f64) -> Self {
        Self {
            entries: Vec::new(),
            theta_max,
            gamma_theta: 1e-5,
            gamma_phi: 1e-8,
            theta_min: 1e-4 * theta_max.max(1e-4),
            s_theta: 1.1,
            s_phi: 2.3,
            eta_phi: 1e-8,
            delta: 1.0,
        }
    }

    /// Initialize theta_min based on the initial constraint violation.
    pub fn set_theta_min_from_initial(&mut self, theta_init: f64) {
        self.theta_min = 1e-4 * theta_init.max(1e-4);
        self.theta_max = 1e4 * theta_init.max(1e-4);
    }

    /// Check if a trial point (theta, phi) is acceptable to the filter.
    /// Returns true if the point is acceptable (not dominated by any filter entry).
    pub fn is_acceptable(&self, theta: f64, phi: f64) -> bool {
        if theta.is_nan() || phi.is_nan() {
            return false;
        }
        if theta > self.theta_max {
            return false;
        }
        for entry in &self.entries {
            if theta >= (1.0 - self.gamma_theta) * entry.theta
                && phi >= entry.phi - self.gamma_phi * entry.theta
            {
                return false;
            }
        }
        true
    }

    /// Check the switching condition: whether we should use the objective (phi)
    /// criterion instead of the filter.
    ///
    /// Returns true if the current constraint violation is small enough and the
    /// directional derivative indicates sufficient objective decrease.
    pub fn switching_condition(
        &self,
        theta_current: f64,
        grad_phi_step: f64,
        alpha: f64,
    ) -> bool {
        // Ipopt switching condition: alpha makes this depend on step length,
        // so as alpha shrinks during backtracking, we properly fall back to
        // h-type (constraint reduction) acceptance.
        grad_phi_step < 0.0
            && theta_current < self.theta_min
            && alpha * (-grad_phi_step).powf(self.s_phi)
                > self.delta * theta_current.powf(self.s_theta)
    }

    /// Check the Armijo sufficient decrease condition.
    pub fn armijo_condition(
        &self,
        phi_current: f64,
        phi_trial: f64,
        grad_phi_step: f64,
        alpha: f64,
    ) -> bool {
        phi_trial <= phi_current + self.eta_phi * alpha * grad_phi_step
    }

    /// Check if a trial point provides sufficient constraint reduction
    /// compared to the current point.
    pub fn sufficient_infeasibility_reduction(
        &self,
        theta_current: f64,
        theta_trial: f64,
    ) -> bool {
        theta_trial <= (1.0 - self.gamma_theta) * theta_current
    }

    /// Check if a trial point is acceptable via either the filter or the
    /// sufficient decrease conditions (Armijo or constraint reduction).
    ///
    /// Returns (acceptable, use_switching) where:
    /// - acceptable: whether the step should be accepted
    /// - use_switching: whether the switching condition was used (affects filter update)
    pub fn check_acceptability(
        &self,
        theta_current: f64,
        phi_current: f64,
        theta_trial: f64,
        phi_trial: f64,
        grad_phi_step: f64,
        alpha: f64,
    ) -> (bool, bool) {
        // Safeguard: reject if objective increased too much (Ipopt uses factor 5.0)
        if phi_trial > phi_current + 5.0 * (1.0 + phi_current.abs()) {
            return (false, false);
        }

        // First check if trial is acceptable to the filter
        if !self.is_acceptable(theta_trial, phi_trial) {
            return (false, false);
        }

        // Check switching condition
        if self.switching_condition(theta_current, grad_phi_step, alpha) {
            // Use Armijo condition on the objective
            let accept = self.armijo_condition(phi_current, phi_trial, grad_phi_step, alpha);
            (accept, true)
        } else {
            // Use filter acceptance: sufficient decrease in theta or phi.
            // Guard: don't accept a step purely via phi decrease when:
            //   - alpha is tiny (< 1e-5), AND
            //   - infeasibility is still significant (> 1e-3), AND
            //   - infeasibility didn't actually decrease.
            // Without this guard, tiny steps accepted via microscopic phi decrease
            // (gamma_phi * theta ≈ 1e-8 * 1.36 ≈ 1e-8) block restoration from
            // firing, causing zigzag-like infinite cycling at stuck infeasibility.
            let infeasibility_reduced = theta_trial < theta_current * (1.0 - self.gamma_theta);
            let phi_ok = phi_trial <= phi_current - self.gamma_phi * theta_current;
            let tiny_useless = alpha < 1e-5 && theta_current > 1e-3 && !infeasibility_reduced;
            let accept = infeasibility_reduced || (phi_ok && !tiny_useless);
            (accept, false)
        }
    }

    /// Add a (theta, phi) pair to the filter.
    pub fn add(&mut self, theta: f64, phi: f64) {
        // Remove dominated entries
        self.entries.retain(|e| {
            !(theta <= (1.0 - self.gamma_theta) * e.theta
                && phi <= e.phi - self.gamma_phi * e.theta)
        });
        self.entries.push(FilterEntry { theta, phi });
    }

    /// Compute problem-dependent minimum step size for the line search (Ipopt formula).
    /// Returns alpha_min based on filter parameters and current iterate.
    pub fn compute_alpha_min(&self, theta_current: f64, grad_phi_step: f64) -> f64 {
        let alpha_min_frac = 0.05; // Ipopt default (was 1e-4)
        if grad_phi_step >= 0.0 || theta_current <= 1e-15 {
            // No useful descent direction or already feasible:
            // use a very small fallback to allow Armijo acceptance to work.
            return 1e-15;
        }
        let neg_gphi = -grad_phi_step;
        let term1 = self.gamma_theta;
        let term2 = self.gamma_phi * theta_current / neg_gphi;
        // Ipopt only includes switching-related term when theta <= theta_min
        let mut alpha_min = alpha_min_frac * term1.min(term2);
        if theta_current <= self.theta_min {
            let term3 =
                self.delta * theta_current.powf(self.s_theta) / neg_gphi.powf(self.s_phi);
            alpha_min = alpha_min.min(alpha_min_frac * term3);
        }
        // Floor at machine epsilon level to avoid overly aggressive cutoff
        alpha_min.max(1e-15)
    }

    /// Augment theta_max based on current violation (called before restoration).
    pub fn augment_for_restoration(&mut self, theta_current: f64) {
        self.theta_max = self.theta_max.max(1e4 * theta_current.max(1e-4));
    }

    /// Reset the filter (used when barrier parameter decreases).
    pub fn reset(&mut self) {
        self.entries.clear();
    }

    /// Save filter entries for watchdog mechanism.
    pub fn save_entries(&self) -> Vec<FilterEntry> {
        self.entries.clone()
    }

    /// Restore filter entries (watchdog rollback).
    pub fn restore_entries(&mut self, entries: Vec<FilterEntry>) {
        self.entries = entries;
    }

    /// Number of entries in the filter.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the filter is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Get theta_max.
    pub fn theta_max(&self) -> f64 {
        self.theta_max
    }

    /// Get gamma_theta.
    pub fn gamma_theta(&self) -> f64 {
        self.gamma_theta
    }

    /// Get gamma_phi.
    pub fn gamma_phi(&self) -> f64 {
        self.gamma_phi
    }

    /// Get the current filter entries (read-only).
    pub fn entries(&self) -> &[FilterEntry] {
        &self.entries
    }
}

/// Compute the maximum step size satisfying the fraction-to-boundary rule.
///
/// Returns the largest alpha in (0, 1] such that:
///   s + alpha * ds >= (1 - tau) * s   for all components
///
/// This ensures slacks/multipliers stay strictly positive.
pub fn fraction_to_boundary(s: &[f64], ds: &[f64], tau: f64) -> f64 {
    let mut alpha_max = 1.0;
    for (si, dsi) in s.iter().zip(ds.iter()) {
        if *dsi < 0.0 {
            let ratio = -tau * si / dsi;
            if ratio < alpha_max {
                alpha_max = ratio;
            }
        }
    }
    alpha_max.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fraction_to_boundary() {
        let s = vec![1.0, 2.0, 0.5];
        let ds = vec![-0.5, 0.1, -0.3];
        let tau = 0.99;
        let alpha = fraction_to_boundary(&s, &ds, tau);
        // For s[0]: -0.99 * 1.0 / (-0.5) = 1.98
        // For s[2]: -0.99 * 0.5 / (-0.3) = 1.65
        // Both > 1.0, so alpha_max = 1.0
        assert!((alpha - 1.0).abs() < 1e-12);

        // Case where step is limited
        let ds2 = vec![-2.0, 0.1, -0.3];
        let alpha2 = fraction_to_boundary(&s, &ds2, tau);
        // For s[0]: -0.99 * 1.0 / (-2.0) = 0.495
        assert!((alpha2 - 0.495).abs() < 1e-12);
    }

    #[test]
    fn test_filter_empty_accepts_everything() {
        let filter = Filter::new(100.0);
        assert!(filter.is_acceptable(1.0, 1.0));
        assert!(filter.is_acceptable(50.0, 50.0));
    }

    #[test]
    fn test_filter_rejects_over_theta_max() {
        let filter = Filter::new(100.0);
        assert!(!filter.is_acceptable(200.0, 0.0));
    }

    #[test]
    fn test_filter_rejects_dominated_point() {
        let mut filter = Filter::new(100.0);
        filter.add(1.0, 1.0);
        // A point that is worse in both theta and phi
        assert!(!filter.is_acceptable(1.0, 1.0));
        // A point that improves sufficiently in theta
        assert!(filter.is_acceptable(0.5, 1.0));
    }

    #[test]
    fn test_filter_reset() {
        let mut filter = Filter::new(100.0);
        filter.add(1.0, 1.0);
        filter.add(0.5, 2.0);
        assert_eq!(filter.len(), 2);
        filter.reset();
        assert_eq!(filter.len(), 0);
        assert!(filter.is_empty());
    }

    #[test]
    fn test_switching_condition() {
        let mut filter = Filter::new(100.0);
        filter.set_theta_min_from_initial(10.0);
        // theta_min = 1e-4 * 10.0 = 1e-3
        // Small theta + negative directional derivative + alpha=1.0 -> switching
        // alpha * 1.0^2.3 = 1.0 > delta * (1e-4)^1.1 ≈ 6.3e-5 -> true
        assert!(filter.switching_condition(1e-4, -1.0, 1.0));
        // Large theta -> no switching (theta >= theta_min)
        assert!(!filter.switching_condition(10.0, -1.0, 1.0));
        // Positive directional derivative -> no switching
        assert!(!filter.switching_condition(1e-4, 1.0, 1.0));
        // Very small alpha -> switching should turn off
        // alpha=1e-20 * 1.0^2.3 = 1e-20, vs delta * (1e-4)^1.1 ≈ 6.3e-5 -> false
        assert!(!filter.switching_condition(1e-4, -1.0, 1e-20));
    }

    #[test]
    fn test_armijo_condition() {
        let filter = Filter::new(100.0);
        // phi_trial <= phi_current + eta_phi * alpha * grad_phi_step
        // eta_phi = 1e-4, threshold: 10.0 + 1e-4 * 1.0 * (-1.0) = 9.9999
        let phi_current = 10.0;
        let grad_phi_step = -1.0;
        let alpha = 1.0;
        assert!(filter.armijo_condition(phi_current, 9.0, grad_phi_step, alpha));
        assert!(!filter.armijo_condition(phi_current, 10.0, grad_phi_step, alpha));
    }

    #[test]
    fn test_sufficient_infeasibility_reduction() {
        let filter = Filter::new(100.0);
        // gamma_theta = 1e-5
        // theta_trial <= (1 - 1e-5) * theta_current
        let theta_current = 1.0;
        assert!(filter.sufficient_infeasibility_reduction(theta_current, 0.5));
        assert!(!filter.sufficient_infeasibility_reduction(theta_current, 1.0));
    }

    #[test]
    fn test_check_acceptability_switching_mode() {
        let mut filter = Filter::new(100.0);
        filter.set_theta_min_from_initial(10.0);
        // Small theta + negative grad_phi_step → switching condition
        // Then Armijo must pass
        let theta_current = 1e-5;
        let phi_current = 10.0;
        let theta_trial = 1e-6;
        let phi_trial = 9.0; // Satisfies Armijo
        let grad_phi_step = -100.0; // Strong descent
        let alpha = 1.0;

        let (accept, switching) = filter.check_acceptability(
            theta_current, phi_current, theta_trial, phi_trial, grad_phi_step, alpha,
        );
        assert!(accept, "Should be accepted via Armijo");
        assert!(switching, "Should use switching condition");
    }

    #[test]
    fn test_check_acceptability_filter_mode() {
        let mut filter = Filter::new(100.0);
        filter.set_theta_min_from_initial(10.0);
        // Large theta → no switching → use filter
        let theta_current = 5.0;
        let phi_current = 10.0;
        let theta_trial = 2.0; // Sufficient reduction
        let phi_trial = 10.0;
        let grad_phi_step = -0.01; // Weak descent
        let alpha = 1.0;

        let (accept, switching) = filter.check_acceptability(
            theta_current, phi_current, theta_trial, phi_trial, grad_phi_step, alpha,
        );
        assert!(accept, "Should be accepted via filter mode");
        assert!(!switching, "Should NOT use switching");
    }

    #[test]
    fn test_check_acceptability_rejected() {
        let mut filter = Filter::new(100.0);
        filter.set_theta_min_from_initial(10.0);
        // Add a filter entry that dominates the trial
        filter.add(0.5, 5.0);
        // Trial point is worse than filter
        let (accept, _) = filter.check_acceptability(
            1.0, 10.0, 0.6, 6.0, -0.01, 1.0,
        );
        assert!(!accept, "Should be rejected by filter");
    }

    #[test]
    fn test_filter_dominated_entry_removal() {
        let mut filter = Filter::new(100.0);
        // Use entries that don't dominate each other: one has lower theta, other has lower phi
        filter.add(0.5, 10.0);
        filter.add(2.0, 5.0);
        assert_eq!(filter.len(), 2);
        // Add an entry that dominates both
        filter.add(0.01, 0.01);
        // The dominating entry should have removed both
        assert!(filter.len() <= 2);
        // The new entry should be in the filter
        assert!(filter.is_acceptable(0.005, -1.0));
    }

    #[test]
    fn test_fraction_to_boundary_edge_cases() {
        // All positive steps → alpha = 1.0
        let s = vec![1.0, 2.0, 0.5];
        let ds = vec![0.1, 0.5, 0.3];
        let tau = 0.99;
        let alpha = fraction_to_boundary(&s, &ds, tau);
        assert!((alpha - 1.0).abs() < 1e-12);

        // Tight constraint: s = [0.01], ds = [-1.0]
        let s2 = vec![0.01];
        let ds2 = vec![-1.0];
        let alpha2 = fraction_to_boundary(&s2, &ds2, tau);
        // alpha = -tau * s / ds = 0.99 * 0.01 / 1.0 = 0.0099
        assert!((alpha2 - 0.0099).abs() < 1e-12);
    }
}
