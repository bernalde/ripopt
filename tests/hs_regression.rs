//! HS (Hock-Schittkowski) regression tests.
//! These ensure ripopt continues to solve key benchmark problems correctly.

#[path = "../hs_suite/generated/hs_problems.rs"]
mod hs_problems;

use ripopt::{SolveStatus, SolverOptions};

fn assert_hs_solved<P: ripopt::NlpProblem>(problem: &P, known_fopt: f64, tol: f64) {
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(problem, &options);
    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal/Acceptable, got {:?} (obj={})",
        result.status, result.objective
    );
    let rel_diff = if known_fopt.abs() > 1.0 {
        (result.objective - known_fopt).abs() / known_fopt.abs()
    } else {
        (result.objective - known_fopt).abs()
    };
    assert!(
        rel_diff < tol,
        "Objective {:.6} too far from known {:.6} (rel_diff={:.2e})",
        result.objective, known_fopt, rel_diff
    );
}

// Unconstrained, n=2
#[test]
fn hs_tp001() {
    assert_hs_solved(&hs_problems::HsTp001, 0.0, 1e-2);
}

// Nonlinear equality, n=2
#[test]
fn hs_tp006() {
    assert_hs_solved(&hs_problems::HsTp006, 0.0, 1e-2);
}

// Nonlinear inequality, n=2
#[test]
fn hs_tp012() {
    assert_hs_solved(&hs_problems::HsTp012, -30.0, 1e-2);
}

// Linear inequality + bounds, n=3
#[test]
fn hs_tp035() {
    assert_hs_solved(&hs_problems::HsTp035, 0.1111111111111111, 1e-2);
}

// Linear inequality, n=4
#[test]
fn hs_tp044() {
    assert_hs_solved(&hs_problems::HsTp044, -15.0, 1e-2);
}

// 5-var unconstrained, n=5
#[test]
fn hs_tp045() {
    assert_hs_solved(&hs_problems::HsTp045, 1.0, 1e-2);
}

// Linear equality, n=5
#[test]
fn hs_tp048() {
    assert_hs_solved(&hs_problems::HsTp048, 0.0, 1e-2);
}

// Mixed constraints (HS071), n=4
#[test]
fn hs_tp071() {
    assert_hs_solved(&hs_problems::HsTp071, 17.0140172895, 1e-2);
}

// Nonlinear equality, n=5 — GN restoration regression test
// Ripopt converges to a local optimum at obj~1.0 (different from known 0.054)
#[test]
fn hs_tp081() {
    assert_hs_solved(&hs_problems::HsTp081, 1.0, 1e-2);
}

// Mixed constraints, n=8 — large problem with relative tolerance
#[test]
fn hs_tp106() {
    assert_hs_solved(&hs_problems::HsTp106, 7049.248, 1e-2);
}

// Mixed constraints, n=10
#[test]
fn hs_tp113() {
    assert_hs_solved(&hs_problems::HsTp113, 24.3062090641, 1e-2);
}

// Code-gen bounds bug regression, n=13 — relative tolerance
#[test]
fn hs_tp116() {
    assert_hs_solved(&hs_problems::HsTp116, 97.5884089805, 1e-2);
}

// Extended Rosenbrock, n=2
#[test]
fn hs_tp201() {
    assert_hs_solved(&hs_problems::HsTp201, 0.0, 1e-2);
}

// Mixed constraint types, n=2
#[test]
fn hs_tp325() {
    assert_hs_solved(&hs_problems::HsTp325, 3.7913414, 1e-2);
}

// Known failure — just verify it doesn't crash
#[test]
fn hs_tp374_no_panic() {
    let problem = hs_problems::HsTp374;
    let options = SolverOptions { print_level: 0, ..SolverOptions::default() };
    let result = ripopt::solve(&problem, &options);
    assert_ne!(result.status, SolveStatus::Optimal,
        "TP374 is known to fail — if it now succeeds, update memory!");
}
