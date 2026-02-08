/// Solver options matching Ipopt defaults.
#[derive(Debug, Clone)]
pub struct SolverOptions {
    /// Convergence tolerance for optimality.
    pub tol: f64,
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Acceptable convergence tolerance (less strict).
    pub acceptable_tol: f64,
    /// Number of consecutive acceptable iterations before declaring success.
    pub acceptable_iter: usize,
    /// Initial barrier parameter.
    pub mu_init: f64,
    /// Minimum barrier parameter.
    pub mu_min: f64,
    /// Fraction-to-boundary parameter minimum.
    pub tau_min: f64,
    /// Barrier parameter reduction factor (monotone mode).
    pub mu_linear_decrease_factor: f64,
    /// Barrier parameter superlinear decrease power.
    pub mu_superlinear_decrease_power: f64,
    /// Print level (0 = silent, 5 = verbose).
    pub print_level: u8,
    /// Bound push for initial point (kappa_1 in Ipopt).
    pub bound_push: f64,
    /// Bound fraction for initial point (kappa_2 in Ipopt).
    pub bound_frac: f64,
    /// Slack bound push.
    pub slack_bound_push: f64,
    /// Slack bound fraction.
    pub slack_bound_frac: f64,
    /// Constraint violation tolerance for convergence.
    pub constr_viol_tol: f64,
    /// Dual infeasibility tolerance for convergence.
    pub dual_inf_tol: f64,
    /// Complementarity tolerance for convergence.
    pub compl_inf_tol: f64,
    /// Use adaptive barrier parameter update (vs monotone).
    pub mu_strategy_adaptive: bool,
    /// Maximum number of second-order correction steps.
    pub max_soc: usize,
    /// Warm-start initialization enabled.
    pub warm_start: bool,
    /// Warm-start bound push.
    pub warm_start_bound_push: f64,
    /// Warm-start bound fraction.
    pub warm_start_bound_frac: f64,
    /// Warm-start multiplier initial value.
    pub warm_start_mult_bound_push: f64,
    /// Any bound less than this value is treated as -infinity (no bound).
    /// Set to a finite value to add artificial lower bounds on unbounded variables.
    pub nlp_lower_bound_inf: f64,
    /// Any bound greater than this value is treated as +infinity (no bound).
    /// Set to a finite value to add artificial upper bounds on unbounded variables.
    pub nlp_upper_bound_inf: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            max_iter: 3000,
            acceptable_tol: 1e-6,
            acceptable_iter: 15,
            mu_init: 0.1,
            mu_min: 1e-11,
            tau_min: 0.99,
            mu_linear_decrease_factor: 0.2,
            mu_superlinear_decrease_power: 1.5,
            print_level: 5,
            bound_push: 1e-2,
            bound_frac: 1e-2,
            slack_bound_push: 1e-2,
            slack_bound_frac: 1e-2,
            constr_viol_tol: 1e-4,
            dual_inf_tol: 1.0,
            compl_inf_tol: 1e-4,
            mu_strategy_adaptive: true,
            max_soc: 4,
            warm_start: false,
            warm_start_bound_push: 1e-3,
            warm_start_bound_frac: 1e-3,
            warm_start_mult_bound_push: 1e-3,
            nlp_lower_bound_inf: -1e19,
            nlp_upper_bound_inf: 1e19,
        }
    }
}
