/// Status of the solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStatus {
    /// Converged to optimal solution within tolerance.
    Optimal,
    /// Converged to acceptable solution (less strict tolerances).
    Acceptable,
    /// Problem is infeasible.
    Infeasible,
    /// Reached maximum number of iterations.
    MaxIterations,
    /// Numerical difficulties (e.g., singular KKT system).
    NumericalError,
    /// Problem appears unbounded below.
    Unbounded,
    /// Restoration phase failed.
    RestorationFailed,
    /// Internal error.
    InternalError,
}

/// Result of solving an NLP.
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// Optimal primal variables x*.
    pub x: Vec<f64>,
    /// Optimal objective value f(x*).
    pub objective: f64,
    /// Constraint multipliers (lambda).
    pub constraint_multipliers: Vec<f64>,
    /// Lower bound multipliers (z_L).
    pub bound_multipliers_lower: Vec<f64>,
    /// Upper bound multipliers (z_U).
    pub bound_multipliers_upper: Vec<f64>,
    /// Constraint values g(x*).
    pub constraint_values: Vec<f64>,
    /// Solve status.
    pub status: SolveStatus,
    /// Number of iterations performed.
    pub iterations: usize,
}
