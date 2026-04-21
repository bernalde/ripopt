//! HS (Hock-Schittkowski) regression tests.
//! These ensure ripopt continues to solve key benchmark problems correctly.

#[path = "hs_problems_subset.rs"]
mod hs_problems;

use ripopt::{SolveStatus, SolverOptions};

fn assert_hs_solved<P: ripopt::NlpProblem>(problem: &P) {
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(problem, &options);
    assert!(
        result.status == SolveStatus::Optimal,
        "Expected Optimal, got {:?} (obj={})",
        result.status, result.objective
    );

    // Check constraint feasibility.
    let m = problem.num_constraints();
    if m > 0 {
        let mut g_l = vec![f64::NEG_INFINITY; m];
        let mut g_u = vec![f64::INFINITY; m];
        problem.constraint_bounds(&mut g_l, &mut g_u);
        for i in 0..m {
            let gi = result.constraint_values[i];
            let viol = if gi < g_l[i] { g_l[i] - gi } else if gi > g_u[i] { gi - g_u[i] } else { 0.0 };
            assert!(
                viol < 1e-4,
                "Constraint {} violated by {:.2e}: g={:.6}, bounds=[{:.6},{:.6}]",
                i, viol, gi, g_l[i], g_u[i]
            );
        }
    }

    // Check variable bounds.
    let n = problem.num_variables();
    let mut x_l = vec![f64::NEG_INFINITY; n];
    let mut x_u = vec![f64::INFINITY; n];
    problem.bounds(&mut x_l, &mut x_u);
    for i in 0..n {
        if x_l[i].is_finite() {
            assert!(result.x[i] >= x_l[i] - 1e-4, "x[{}]={:.6} below bound {:.6}", i, result.x[i], x_l[i]);
        }
        if x_u[i].is_finite() {
            assert!(result.x[i] <= x_u[i] + 1e-4, "x[{}]={:.6} above bound {:.6}", i, result.x[i], x_u[i]);
        }
    }
}

// Unconstrained, n=2
#[test]
fn hs_tp001() {
    assert_hs_solved(&hs_problems::HsTp001);
}

// Nonlinear equality, n=2
#[test]
fn hs_tp006() {
    assert_hs_solved(&hs_problems::HsTp006);
}

// Nonlinear inequality, n=2
#[test]
fn hs_tp012() {
    assert_hs_solved(&hs_problems::HsTp012);
}

// Linear inequality + bounds, n=3
#[test]
fn hs_tp035() {
    assert_hs_solved(&hs_problems::HsTp035);
}

// Linear inequality, n=4
// TODO(z_opt-refactor): blocked by dz-step corruption at active bounds exposed by z_opt removal.
#[test]
#[ignore = "blocked by hs071 dz-step corruption exposed by z_opt removal"]
fn hs_tp044() {
    assert_hs_solved(&hs_problems::HsTp044);
}

// 5-var unconstrained, n=5
// Known failure: solver reaches near-tolerance but not strict Optimal.
#[test]
#[ignore = "known solver limitation: does not reach strict Optimal (previously hid as Acceptable)"]
fn hs_tp045() {
    assert_hs_solved(&hs_problems::HsTp045);
}

// Linear equality, n=5
#[test]
fn hs_tp048() {
    assert_hs_solved(&hs_problems::HsTp048);
}

// Mixed constraints (HS071), n=4
// TODO(z_opt-refactor): blocked by dz-step corruption at active bounds exposed by z_opt removal.
#[test]
#[ignore = "blocked by hs071 dz-step corruption exposed by z_opt removal"]
fn hs_tp071() {
    assert_hs_solved(&hs_problems::HsTp071);
}

// Nonlinear equality, n=5 — multiple local optima, non-convex.
#[test]
fn hs_tp081() {
    assert_hs_solved(&hs_problems::HsTp081);
}

// Mixed constraints, n=8
#[test]
fn hs_tp106() {
    assert_hs_solved(&hs_problems::HsTp106);
}

// Mixed constraints, n=10
#[test]
fn hs_tp113() {
    assert_hs_solved(&hs_problems::HsTp113);
}

// Code-gen bounds bug regression, n=13
// Known failure: solver reaches near-tolerance but not strict Optimal.
#[test]
#[ignore = "known solver limitation: does not reach strict Optimal (previously hid as Acceptable)"]
fn hs_tp116() {
    assert_hs_solved(&hs_problems::HsTp116);
}

// Extended Rosenbrock, n=2
#[test]
fn hs_tp201() {
    assert_hs_solved(&hs_problems::HsTp201);
}

// Mixed constraint types, n=2
#[test]
fn hs_tp325() {
    assert_hs_solved(&hs_problems::HsTp325);
}

// TP374: known solver limitation — does not reach strict Optimal.
#[test]
#[ignore = "known solver limitation: does not reach strict Optimal"]
fn hs_tp374() {
    assert_hs_solved(&hs_problems::HsTp374);
}
