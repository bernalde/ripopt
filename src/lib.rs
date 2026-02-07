pub mod convergence;
pub mod filter;
pub mod ipm;
pub mod kkt;
pub mod linear_solver;
pub mod options;
pub mod problem;
pub mod restoration;
pub mod result;
pub mod warmstart;

pub use options::SolverOptions;
pub use problem::NlpProblem;
pub use result::{SolveResult, SolveStatus};

/// Solve a nonlinear programming problem using the interior point method.
pub fn solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    ipm::solve(problem, options)
}
