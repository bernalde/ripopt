use crate::options::SolverOptions;

/// Result of a convergence check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceStatus {
    /// Not converged.
    NotConverged,
    /// Converged within desired tolerance.
    Converged,
    /// Converged within acceptable tolerance.
    Acceptable,
    /// Diverging (objective growing unboundedly).
    Diverging,
}

/// Information needed to check convergence.
pub struct ConvergenceInfo {
    /// Primal infeasibility: max |c_i(x)| for violated constraints.
    pub primal_inf: f64,
    /// Dual infeasibility: ||grad_f - J^T lambda - z||_inf.
    pub dual_inf: f64,
    /// Complementarity: max |x_i * z_i - mu|.
    pub compl_inf: f64,
    /// Current barrier parameter.
    pub mu: f64,
    /// Current objective value.
    pub objective: f64,
    /// Sum of absolute values of all multipliers (y, z_l, z_u).
    /// Used for Ipopt-style dual scaling.
    pub multiplier_sum: f64,
    /// Total number of multiplier components (n + m for scaling denominator).
    pub multiplier_count: usize,
}

/// Check convergence of the IPM algorithm.
///
/// Returns the convergence status based on current optimality measures.
pub fn check_convergence(
    info: &ConvergenceInfo,
    options: &SolverOptions,
    consecutive_acceptable: usize,
) -> ConvergenceStatus {
    // Ipopt-style dual scaling: s_d = max(s_max, sum|mults| / count) / s_max
    // This scales the tolerance to account for large multiplier magnitudes.
    let s_max: f64 = 100.0;
    let s_d = if info.multiplier_count > 0 {
        (s_max.max(info.multiplier_sum / info.multiplier_count as f64)) / s_max
    } else {
        1.0
    };

    let primal_tol = options.tol;
    let dual_tol = options.tol * s_d;
    let compl_tol = options.tol * s_d;

    // Check strict convergence
    if info.primal_inf <= primal_tol
        && info.dual_inf <= dual_tol
        && info.compl_inf <= compl_tol
    {
        return ConvergenceStatus::Converged;
    }

    // Check acceptable convergence (also scaled)
    let acc_primal_tol = options.acceptable_tol;
    let acc_dual_tol = options.acceptable_tol * s_d;
    let acc_compl_tol = options.acceptable_tol * s_d;

    if info.primal_inf <= acc_primal_tol
        && info.dual_inf <= acc_dual_tol
        && info.compl_inf <= acc_compl_tol
        && consecutive_acceptable >= options.acceptable_iter
    {
        return ConvergenceStatus::Acceptable;
    }

    // Check divergence
    if info.objective.abs() > 1e20 {
        return ConvergenceStatus::Diverging;
    }

    ConvergenceStatus::NotConverged
}

/// Compute primal infeasibility (constraint violation).
/// Takes constraint values g(x), and constraint bounds g_l, g_u.
pub fn primal_infeasibility(g: &[f64], g_l: &[f64], g_u: &[f64]) -> f64 {
    let mut max_viol = 0.0f64;
    for i in 0..g.len() {
        if g[i] < g_l[i] {
            max_viol = max_viol.max(g_l[i] - g[i]);
        }
        if g[i] > g_u[i] {
            max_viol = max_viol.max(g[i] - g_u[i]);
        }
    }
    max_viol
}

/// Compute dual infeasibility: ||grad_f - J^T * lambda - z_l + z_u||_inf.
///
/// `grad_f`: gradient of objective
/// `jac_rows`, `jac_cols`, `jac_vals`: Jacobian in COO format
/// `lambda`: constraint multipliers
/// `z_l`, `z_u`: bound multipliers
/// `n`: number of variables
#[allow(clippy::too_many_arguments)]
pub fn dual_infeasibility(
    grad_f: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    lambda: &[f64],
    z_l: &[f64],
    z_u: &[f64],
    n: usize,
) -> f64 {
    let mut residual = vec![0.0; n];

    // Start with gradient of objective
    residual[..n].copy_from_slice(&grad_f[..n]);

    // Add J^T * lambda (Ipopt convention: L = f + y^T g)
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        residual[col] += jac_vals[idx] * lambda[row];
    }

    // Subtract z_l and add z_u (bound multipliers)
    for i in 0..n {
        residual[i] -= z_l[i];
        residual[i] += z_u[i];
    }

    residual.iter().map(|r| r.abs()).fold(0.0f64, f64::max)
}

/// Compute complementarity error for bound constraints.
/// compl = max_i |x_i * z_l_i| where x_i is near lower bound,
///         max_i |s_u_i * z_u_i| where x_i is near upper bound.
///
/// For the barrier method: complementarity = max(|(x-x_l)*z_l - mu|, |(x_u-x)*z_u - mu|).
pub fn complementarity_error(
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    z_l: &[f64],
    z_u: &[f64],
    mu: f64,
) -> f64 {
    let mut max_err = 0.0f64;
    let n = x.len();
    for i in 0..n {
        if x_l[i].is_finite() {
            let slack = x[i] - x_l[i];
            max_err = max_err.max((slack * z_l[i] - mu).abs());
        }
        if x_u[i].is_finite() {
            let slack = x_u[i] - x[i];
            max_err = max_err.max((slack * z_u[i] - mu).abs());
        }
    }
    max_err
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primal_infeasibility_feasible() {
        let g = vec![1.5, 3.0];
        let g_l = vec![1.0, 2.0];
        let g_u = vec![2.0, 4.0];
        assert_eq!(primal_infeasibility(&g, &g_l, &g_u), 0.0);
    }

    #[test]
    fn test_primal_infeasibility_violated() {
        let g = vec![0.5, 5.0];
        let g_l = vec![1.0, 2.0];
        let g_u = vec![2.0, 4.0];
        assert_eq!(primal_infeasibility(&g, &g_l, &g_u), 1.0);
    }

    #[test]
    fn test_convergence_optimal() {
        let info = ConvergenceInfo {
            primal_inf: 1e-10,
            dual_inf: 1e-10,
            compl_inf: 1e-10,
            mu: 1e-11,
            objective: 17.0,
            multiplier_sum: 0.0,
            multiplier_count: 0,
        };
        let opts = SolverOptions::default();
        assert_eq!(
            check_convergence(&info, &opts, 0),
            ConvergenceStatus::Converged
        );
    }

    #[test]
    fn test_convergence_not_converged() {
        let info = ConvergenceInfo {
            primal_inf: 1e-3,
            dual_inf: 1e-3,
            compl_inf: 1e-3,
            mu: 0.01,
            objective: 17.0,
            multiplier_sum: 0.0,
            multiplier_count: 0,
        };
        let opts = SolverOptions::default();
        assert_eq!(
            check_convergence(&info, &opts, 0),
            ConvergenceStatus::NotConverged
        );
    }
}
