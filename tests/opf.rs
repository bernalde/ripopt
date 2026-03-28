//! AC Optimal Power Flow test suite for ripopt.
//!
//! Run with: cargo test opf -- --nocapture

#[path = "common/opf_problems.rs"]
mod problems;
use problems::*;

use ripopt::{NlpProblem, SolveStatus, SolverOptions};
use std::time::Instant;

fn default_options() -> SolverOptions {
    SolverOptions {
        tol: 1e-6,
        max_iter: 3000,
        print_level: 0,
        ..SolverOptions::default()
    }
}

/// Compute max constraint violation.
fn max_cv(problem: &dyn NlpProblem, g: &[f64]) -> f64 {
    let m = problem.num_constraints();
    if m == 0 {
        return 0.0;
    }
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    let mut cv = 0.0_f64;
    for i in 0..m {
        cv = cv.max((g_l[i] - g[i]).max(0.0)).max((g[i] - g_u[i]).max(0.0));
    }
    cv
}

#[test]
#[ignore = "known solver limitation: does not reach strict Optimal (previously hid as Acceptable)"]
fn opf_case3_lmbd() {
    let problem = case3_lmbd();
    let options = default_options();
    let start = Instant::now();
    let result = ripopt::solve(&problem, &options);
    let elapsed = start.elapsed();
    let cv = max_cv(&problem, &result.constraint_values);
    eprintln!(
        "case3: status={:?}, obj={:.2}, cv={:.2e}, iters={}, time={:.3}s",
        result.status, result.objective, cv, result.iterations, elapsed.as_secs_f64()
    );
    assert!(
        result.status == SolveStatus::Optimal,
        "Expected Optimal/Acceptable, got {:?}", result.status
    );
    assert!(cv < 1e-4, "cv={:.2e}", cv);
    // OPF is non-convex; any locally optimal feasible solution is acceptable.
}

#[test]
#[ignore = "known solver limitation: does not reach strict Optimal (previously hid as Acceptable)"]
fn opf_case5_pjm() {
    let problem = case5_pjm();
    let options = default_options();
    let start = Instant::now();
    let result = ripopt::solve(&problem, &options);
    let elapsed = start.elapsed();
    let cv = max_cv(&problem, &result.constraint_values);
    eprintln!(
        "case5: status={:?}, obj={:.2}, cv={:.2e}, iters={}, time={:.3}s",
        result.status, result.objective, cv, result.iterations, elapsed.as_secs_f64()
    );
    assert!(
        result.status == SolveStatus::Optimal,
        "Expected Optimal/Acceptable, got {:?}", result.status
    );
    assert!(cv < 1e-4, "cv={:.2e}", cv);
    // OPF is non-convex; any locally optimal feasible solution is acceptable.
}

#[test]
fn opf_case14_ieee() {
    let problem = case14_ieee();
    let options = default_options();
    let start = Instant::now();
    let result = ripopt::solve(&problem, &options);
    let elapsed = start.elapsed();
    let cv = max_cv(&problem, &result.constraint_values);
    eprintln!(
        "case14: status={:?}, obj={:.2}, cv={:.2e}, iters={}, time={:.3}s",
        result.status, result.objective, cv, result.iterations, elapsed.as_secs_f64()
    );
    assert!(
        result.status == SolveStatus::Optimal,
        "Expected Optimal/Acceptable, got {:?}", result.status
    );
    assert!(cv < 1e-4, "cv={:.2e}", cv);
    // OPF is non-convex; any locally optimal feasible solution is acceptable.
}

#[test]
fn opf_case30_ieee() {
    let problem = case30_ieee();
    let options = default_options();
    let start = Instant::now();
    let result = ripopt::solve(&problem, &options);
    let elapsed = start.elapsed();
    let cv = max_cv(&problem, &result.constraint_values);
    eprintln!(
        "case30: status={:?}, obj={:.2}, cv={:.2e}, iters={}, time={:.3}s",
        result.status, result.objective, cv, result.iterations, elapsed.as_secs_f64()
    );
    assert!(
        result.status == SolveStatus::Optimal,
        "Expected Optimal/Acceptable, got {:?}", result.status
    );
    assert!(cv < 1e-4, "cv={:.2e}", cv);
    // OPF is non-convex; any locally optimal feasible solution is acceptable.
}
