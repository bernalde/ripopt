/// Status of the solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStatus {
    /// Converged to optimal solution within tolerance.
    Optimal,
    /// Converged to acceptable solution (less strict tolerances).
    Acceptable,
    /// Problem is infeasible.
    Infeasible,
    /// Local infeasibility detected: constraint violation is at a stationary
    /// point (gradient of violation ≈ 0) but violation is still large.
    /// For NE-to-LS reformulations, this means the system is inconsistent
    /// and x* is the best least-squares solution.
    LocalInfeasibility,
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

/// Structured diagnostic summary from a solve.
///
/// Captures counts of key solver events (restoration entries, barrier parameter
/// mode switches, filter rejects, etc.) and final convergence measures.
/// Useful for automated analysis and solver tuning.
#[derive(Debug, Clone, Default)]
pub struct SolverDiagnostics {
    /// Number of GN (Gauss-Newton) restoration entries.
    pub restoration_count: usize,
    /// Number of full NLP restoration entries.
    pub nlp_restoration_count: usize,
    /// Number of mu mode switches (Free↔Fixed).
    pub mu_mode_switches: usize,
    /// Number of filter rejects (line search exhausted backtracking).
    pub filter_rejects: usize,
    /// Number of watchdog activations.
    pub watchdog_activations: usize,
    /// Number of second-order corrections (SOC) applied.
    pub soc_corrections: usize,
    /// Final barrier parameter mu.
    pub final_mu: f64,
    /// Final primal infeasibility.
    pub final_primal_inf: f64,
    /// Final dual infeasibility (iterative z, used in unscaled gate).
    pub final_dual_inf: f64,
    /// Final dual infeasibility (z_opt, used in scaled gate).
    pub final_dual_inf_scaled: f64,
    /// Final complementarity error.
    pub final_compl: f64,
    /// Dual scaling factor s_d.
    pub final_s_d: f64,
    /// Total wall-clock time in seconds.
    pub wall_time_secs: f64,
    /// Fallback strategy used, if any.
    pub fallback_used: Option<String>,
}

impl SolverDiagnostics {
    /// Print a structured diagnostic summary to stderr.
    pub fn print_summary(&self, status: SolveStatus, iterations: usize) {
        eprintln!("\n--- ripopt diagnostics ---");
        eprintln!("status: {:?}", status);
        eprintln!("iterations: {}", iterations);
        eprintln!("wall_time: {:.3}s", self.wall_time_secs);
        eprintln!("final_mu: {:.2e}", self.final_mu);
        eprintln!("final_primal_inf: {:.2e}", self.final_primal_inf);
        eprintln!("final_dual_inf: {:.2e}", self.final_dual_inf);
        eprintln!("final_compl: {:.2e}", self.final_compl);
        eprintln!("restoration_count: {}", self.restoration_count);
        eprintln!("nlp_restoration_count: {}", self.nlp_restoration_count);
        eprintln!("mu_mode_switches: {}", self.mu_mode_switches);
        eprintln!("filter_rejects: {}", self.filter_rejects);
        eprintln!("watchdog_activations: {}", self.watchdog_activations);
        eprintln!("soc_corrections: {}", self.soc_corrections);
        if let Some(ref fb) = self.fallback_used {
            eprintln!("fallback_used: {}", fb);
        }
        eprintln!("--- end diagnostics ---");
    }
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
    /// Structured solver diagnostics.
    pub diagnostics: SolverDiagnostics,
}
