//! # ripopt — Rust Interior Point Optimizer
//!
//! ripopt is a primal-dual interior point solver for nonlinear programming (NLP) problems:
//!
//! ```text
//! min  f(x)
//!  x
//! s.t. g_l ≤ g(x) ≤ g_u
//!      x_l ≤  x   ≤ x_u
//! ```
//!
//! It closely follows the algorithm of [Ipopt](https://coin-or.github.io/Ipopt/) and achieves
//! comparable or better solve rates on standard benchmark suites (HS120, CUTEst).
//!
//! ## Quick Start
//!
//! Implement [`NlpProblem`] for your problem, then call [`solve`]:
//!
//! ```rust,no_run
//! use ripopt::{NlpProblem, SolveResult, SolverOptions, SolveStatus};
//!
//! struct MyProblem;
//!
//! impl NlpProblem for MyProblem {
//!     fn n(&self) -> usize { 2 }
//!     fn m(&self) -> usize { 0 }
//!     fn x_l(&self) -> Vec<f64> { vec![-1e20; 2] }
//!     fn x_u(&self) -> Vec<f64> { vec![ 1e20; 2] }
//!     fn g_l(&self) -> Vec<f64> { vec![] }
//!     fn g_u(&self) -> Vec<f64> { vec![] }
//!     fn x0(&self) -> Vec<f64> { vec![0.5, 0.5] }
//!     fn f(&self, x: &[f64]) -> f64 { (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0].powi(2)).powi(2) }
//!     fn grad_f(&self, x: &[f64]) -> Vec<f64> {
//!         vec![
//!             -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0].powi(2)),
//!              200.0 * (x[1] - x[0].powi(2)),
//!         ]
//!     }
//!     fn g(&self, _x: &[f64]) -> Vec<f64> { vec![] }
//!     fn jac_g(&self, _x: &[f64]) -> Vec<f64> { vec![] }
//!     fn jac_g_sparsity(&self) -> (Vec<usize>, Vec<usize>) { (vec![], vec![]) }
//!     fn hess_l(&self, x: &[f64], _sigma: f64, _lambda: &[f64]) -> Vec<f64> {
//!         let h00 = 2.0 - 400.0 * (x[1] - x[0].powi(2)) + 800.0 * x[0].powi(2);
//!         let h10 = -400.0 * x[0];
//!         let h11 = 200.0_f64;
//!         vec![h00, h10, h11]
//!     }
//!     fn hess_l_sparsity(&self) -> (Vec<usize>, Vec<usize>) {
//!         (vec![0, 1, 1], vec![0, 0, 1])
//!     }
//! }
//!
//! let result = solve(&MyProblem, &SolverOptions::default());
//! assert_eq!(result.status, SolveStatus::Optimal);
//! ```
//!
//! ## Key Types
//!
//! - [`NlpProblem`] — trait to implement for your problem
//! - [`SolverOptions`] — all solver tuning parameters
//! - [`SolveResult`] — solution, status, and diagnostics
//! - [`SolveStatus`] — outcome (`Optimal`, `LocalInfeasibility`, etc.)
//!
//! ## Algorithm
//!
//! The solver implements:
//! - Mehrotra predictor-corrector IPM with Gondzio corrections
//! - Filter line search with second-order corrections
//! - Gauss-Newton and NLP restoration phases
//! - Fallback cascade: IPM → L-BFGS → Augmented Lagrangian → SQP → slack reformulation
//! - Sparse (multifrontal LDL^T) and dense (Bunch-Kaufman LDL^T) linear solvers
//! - Parametric sensitivity analysis ([`solve_with_sensitivity`])

pub mod augmented_lagrangian;
pub(crate) mod logging;
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
