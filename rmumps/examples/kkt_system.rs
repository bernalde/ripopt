//! Solving a KKT system from constrained optimization.
//!
//! KKT systems are symmetric indefinite:
//!
//!   [H   J^T] [dx]   [r_d]
//!   [J  -δI ] [dy] = [r_c]
//!
//! rmumps handles these via Bunch-Kaufman pivoting, which produces 1x1 and 2x2
//! pivots for symmetric indefinite matrices. The inertia (positive/negative/zero
//! eigenvalue counts) tells the optimizer whether the KKT matrix has the correct
//! structure for a local minimum.
//!
//! Run with: cargo run --example kkt_system

use rmumps::coo::CooMatrix;
use rmumps::solver::{Solver, SolverOptions};

fn main() {
    // A small KKT system: 3 variables, 2 constraints
    //
    //   H = [2 0 0]    J = [1 1 0]    delta = 1e-8
    //       [0 2 0]        [0 1 1]
    //       [0 0 2]
    //
    //   KKT = [2  0  0  1  0]
    //         [0  2  0  1  1]
    //         [0  0  2  0  1]
    //         [1  1  0 -δ  0]
    //         [0  1  1  0 -δ]
    //
    let n = 5;
    let delta = 1e-8;

    // Build upper triangle in COO format
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    // H block (diagonal)
    for i in 0..3 {
        rows.push(i);
        cols.push(i);
        vals.push(2.0);
    }
    // -delta*I block
    for i in 3..5 {
        rows.push(i);
        cols.push(i);
        vals.push(-delta);
    }
    // J^T block (upper triangle: row < col)
    // J = [[1,1,0],[0,1,1]] so J^T has entries at (0,3),(1,3),(1,4),(2,4)
    rows.push(0);
    cols.push(3);
    vals.push(1.0);
    rows.push(1);
    cols.push(3);
    vals.push(1.0);
    rows.push(1);
    cols.push(4);
    vals.push(1.0);
    rows.push(2);
    cols.push(4);
    vals.push(1.0);

    let coo = CooMatrix::new(n, rows, cols, vals).unwrap();

    let mut solver = Solver::new(SolverOptions::default());
    let inertia = solver.analyze_and_factor(&coo).unwrap();
    println!("Inertia: {}", inertia);
    // For a well-posed KKT: 3 positive (from H), 2 negative (from constraints)
    assert_eq!(inertia.positive, 3);
    assert_eq!(inertia.negative, 2);

    // Solve KKT system
    let rhs = vec![1.0, 2.0, 3.0, 0.0, 0.0]; // r_d = gradient, r_c = constraint violation
    let mut solution = vec![0.0; n];
    solver.solve(&rhs, &mut solution).unwrap();

    println!("dx = [{:.6}, {:.6}, {:.6}]", solution[0], solution[1], solution[2]);
    println!("dy = [{:.6}, {:.6}]", solution[3], solution[4]);

    // Verify residual
    let mut ax = vec![0.0; n];
    coo.matvec(&solution, &mut ax).unwrap();
    let residual: f64 = ax
        .iter()
        .zip(&rhs)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!("Max residual: {:.2e}", residual);
    assert!(residual < 1e-6);
}
