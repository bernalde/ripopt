//! Basic usage: create a sparse symmetric matrix and solve Ax = b.
//!
//! Run with: cargo run --example basic

use rmumps::coo::CooMatrix;
use rmumps::solver::{Solver, SolverOptions};

fn main() {
    // 3x3 symmetric positive definite matrix (upper triangle, COO format):
    //
    //     [4  1  0]
    // A = [1  5  2]
    //     [0  2  6]
    //
    let coo = CooMatrix::new(
        3,
        vec![0, 1, 2, 0, 1], // rows
        vec![0, 1, 2, 1, 2], // cols (col >= row for upper triangle)
        vec![4.0, 5.0, 6.0, 1.0, 2.0], // values
    )
    .unwrap();

    // Create solver with default options (AMD ordering, iterative refinement)
    let mut solver = Solver::new(SolverOptions::default());

    // Analyze sparsity pattern + factorize
    let inertia = solver.analyze_and_factor(&coo).unwrap();
    println!("Inertia: {}", inertia);
    assert_eq!(inertia.positive, 3); // all positive eigenvalues (SPD)

    // Solve Ax = b
    let b = vec![1.0, 2.0, 3.0];
    let mut x = vec![0.0; 3];
    solver.solve(&b, &mut x).unwrap();

    println!("Solution: {:?}", x);

    // Verify: compute residual ||Ax - b||
    let mut ax = vec![0.0; 3];
    coo.matvec(&x, &mut ax).unwrap();
    let residual: f64 = ax
        .iter()
        .zip(&b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!("Max residual: {:.2e}", residual);
    assert!(residual < 1e-14);
}
