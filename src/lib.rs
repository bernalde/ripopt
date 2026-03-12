pub mod augmented_lagrangian;
pub mod c_api;
pub mod convergence;
pub mod filter;
pub mod ipm;
pub mod kkt;
pub mod lbfgs;
pub mod linearity;
pub mod linear_solver;
pub mod nl;
pub mod options;
pub mod preprocessing;
pub mod problem;
pub mod restoration;
pub mod restoration_nlp;
pub mod result;
pub mod sensitivity;
pub mod slack_formulation;
pub mod sqp;
pub mod warmstart;

pub use options::{SolverOptions, LinearSolverChoice};
pub use problem::NlpProblem;
pub use result::{SolveResult, SolverDiagnostics, SolveStatus};
pub use sensitivity::{ParametricNlpProblem, SensitivityContext, SensitivityResult};

/// Solve a nonlinear programming problem using the interior point method.
pub fn solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    ipm::solve(problem, options)
}

/// Solve and retain factored KKT for parametric sensitivity analysis.
pub fn solve_with_sensitivity<P: ParametricNlpProblem>(
    problem: &P,
    options: &SolverOptions,
) -> SensitivityContext {
    sensitivity::solve_with_sensitivity(problem, options)
}
