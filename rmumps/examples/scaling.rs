//! Demonstrating matrix scaling for ill-conditioned systems.
//!
//! Scaling transforms A*x = b into (D*A*D)*y = D*b, then x = D*y, where D is
//! a diagonal matrix chosen to equilibrate row/column norms. This improves
//! numerical accuracy for ill-conditioned matrices.
//!
//! Run with: cargo run --example scaling

use rmumps::coo::CooMatrix;
use rmumps::scaling::Scaling;
use rmumps::solver::{Solver, SolverOptions};

fn main() {
    // An ill-conditioned 3x3 matrix:
    //     [1e8   1    0 ]
    // A = [1     1e-4 1 ]
    //     [0     1    1e6]
    //
    // Condition number ~1e12: small pivots can lose significant digits.
    let coo = CooMatrix::new(
        3,
        vec![0, 1, 2, 0, 1],
        vec![0, 1, 2, 1, 2],
        vec![1e8, 1e-4, 1e6, 1.0, 1.0],
    )
    .unwrap();

    let b = vec![1e8 + 1.0, 1.0 + 1e-4 + 1.0, 1.0 + 1e6]; // A * [1, 1, 1]
    let x_exact = vec![1.0, 1.0, 1.0];

    // Without scaling
    let mut solver = Solver::new(SolverOptions {
        scaling: Scaling::None,
        ..Default::default()
    });
    solver.analyze_and_factor(&coo).unwrap();
    let mut x_unscaled = vec![0.0; 3];
    solver.solve(&b, &mut x_unscaled).unwrap();
    let err_unscaled: f64 = x_unscaled
        .iter()
        .zip(&x_exact)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);

    // With Ruiz scaling
    let mut solver = Solver::new(SolverOptions {
        scaling: Scaling::Ruiz { max_iter: 10 },
        ..Default::default()
    });
    solver.analyze_and_factor(&coo).unwrap();
    let mut x_scaled = vec![0.0; 3];
    solver.solve(&b, &mut x_scaled).unwrap();
    let err_scaled: f64 = x_scaled
        .iter()
        .zip(&x_exact)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);

    // With diagonal scaling
    let mut solver = Solver::new(SolverOptions {
        scaling: Scaling::Diagonal,
        ..Default::default()
    });
    solver.analyze_and_factor(&coo).unwrap();
    let mut x_diag = vec![0.0; 3];
    solver.solve(&b, &mut x_diag).unwrap();
    let err_diag: f64 = x_diag
        .iter()
        .zip(&x_exact)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);

    println!("Ill-conditioned 3x3 system (exact solution = [1, 1, 1]):");
    println!("  No scaling:       error = {:.2e}", err_unscaled);
    println!("  Diagonal scaling: error = {:.2e}", err_diag);
    println!("  Ruiz scaling:     error = {:.2e}", err_scaled);
}
