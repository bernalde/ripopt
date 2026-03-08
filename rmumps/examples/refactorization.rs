//! Demonstrating the analyze-once, refactor-many workflow.
//!
//! In interior point methods, the KKT matrix changes values each iteration but
//! keeps the same sparsity pattern. rmumps caches the symbolic analysis (ordering,
//! elimination tree, supernodes) so subsequent factorizations only do numeric work.
//!
//! Run with: cargo run --example refactorization

use rmumps::coo::CooMatrix;
use rmumps::solver::{Solver, SolverOptions};
use std::time::Instant;

/// 2D Laplacian on an nx x nx grid with a scalar shift.
fn laplacian_2d(nx: usize, shift: f64) -> CooMatrix {
    let n = nx * nx;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for iy in 0..nx {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            rows.push(idx);
            cols.push(idx);
            vals.push(4.0 + shift);
            if ix + 1 < nx {
                rows.push(idx);
                cols.push(iy * nx + ix + 1);
                vals.push(-1.0);
            }
            if iy + 1 < nx {
                rows.push(idx);
                cols.push((iy + 1) * nx + ix);
                vals.push(-1.0);
            }
        }
    }
    CooMatrix::new(n, rows, cols, vals).unwrap()
}

fn main() {
    let nx = 50;
    let n = nx * nx; // 2500

    // Phase 1: Analyze (once)
    let coo0 = laplacian_2d(nx, 0.0);
    let mut solver = Solver::new(SolverOptions::default());

    let t = Instant::now();
    solver.analyze(&coo0).unwrap();
    println!("Analyze: {:.2}ms (n={})", t.elapsed().as_secs_f64() * 1000.0, n);

    // Phase 2: Factor with different shifts (simulating IPM iterations)
    for shift in [0.0, 0.1, 1.0, 10.0] {
        let coo = laplacian_2d(nx, shift);

        let t = Instant::now();
        let inertia = solver.factor(&coo).unwrap();
        let factor_ms = t.elapsed().as_secs_f64() * 1000.0;

        let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut x = vec![0.0; n];

        let t = Instant::now();
        solver.solve(&b, &mut x).unwrap();
        let solve_ms = t.elapsed().as_secs_f64() * 1000.0;

        // Residual
        let mut ax = vec![0.0; n];
        coo.matvec(&x, &mut ax).unwrap();
        let residual: f64 = ax
            .iter()
            .zip(&b)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        let norm_b: f64 = b.iter().map(|v| v.abs()).fold(0.0, f64::max);

        println!(
            "  shift={:>5.1}: factor {:.2}ms, solve {:.2}ms, inertia={}, rel_resid={:.2e}",
            shift,
            factor_ms,
            solve_ms,
            inertia,
            residual / norm_b
        );
    }
}
