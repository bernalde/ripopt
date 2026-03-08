//! Benchmark on larger problems: 2D Laplacian, KKT, and tridiagonal.
//!
//! Run with: cargo run --release --example large

use rmumps::coo::CooMatrix;
use rmumps::ordering::Ordering;
use rmumps::solver::{Solver, SolverOptions};
use std::time::Instant;

fn laplacian_2d(nx: usize) -> CooMatrix {
    let n = nx * nx;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for iy in 0..nx {
        for ix in 0..nx {
            let idx = iy * nx + ix;
            rows.push(idx);
            cols.push(idx);
            vals.push(4.0);
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

fn kkt_matrix(nvar: usize, ncon: usize) -> CooMatrix {
    let n = nvar + ncon;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..nvar {
        rows.push(i);
        cols.push(i);
        vals.push(4.0 + i as f64 * 0.01);
        if i + 1 < nvar {
            rows.push(i);
            cols.push(i + 1);
            vals.push(0.5);
        }
    }
    for c in 0..ncon {
        let row = nvar + c;
        for k in 0..3 {
            let var = (c * 3 + k) % nvar;
            let (r, col) = if var < row { (var, row) } else { (row, var) };
            rows.push(r);
            cols.push(col);
            vals.push(1.0);
        }
        rows.push(row);
        cols.push(row);
        vals.push(-1e-8);
    }
    CooMatrix::new(n, rows, cols, vals).unwrap()
}

fn bench(name: &str, coo: &CooMatrix) {
    let n = coo.n;
    println!("\n=== {} (n={}, nnz={}) ===", name, n, coo.nnz());

    let mut solver = Solver::new(SolverOptions {
        ordering: Ordering::Amd,
        ..Default::default()
    });

    let t = Instant::now();
    solver.analyze(coo).unwrap();
    let analyze_ms = t.elapsed().as_secs_f64() * 1000.0;
    println!("  Analyze:  {:>10.1}ms", analyze_ms);

    let t = Instant::now();
    solver.factor(coo).unwrap();
    let factor_ms = t.elapsed().as_secs_f64() * 1000.0;
    println!("  Factor:   {:>10.1}ms", factor_ms);

    let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let mut x = vec![0.0; n];
    let t = Instant::now();
    solver.solve(&b, &mut x).unwrap();
    let solve_ms = t.elapsed().as_secs_f64() * 1000.0;
    println!("  Solve:    {:>10.1}ms", solve_ms);

    // Residual
    let mut ax = vec![0.0; n];
    coo.matvec(&x, &mut ax).unwrap();
    let max_resid: f64 = ax
        .iter()
        .zip(&b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let norm_b: f64 = b.iter().map(|v| v.abs()).fold(0.0, f64::max);
    println!("  Residual: {:>14.2e}", max_resid / norm_b);
    println!("  Total:    {:>10.1}ms", analyze_ms + factor_ms + solve_ms);
}

fn main() {
    println!("rmumps large problem benchmark");
    println!("Rayon threads: {}", rayon::current_num_threads());

    bench("Laplacian 100x100", &laplacian_2d(100));
    bench("Laplacian 316x316", &laplacian_2d(316));
    bench("KKT 10000+5000", &kkt_matrix(10000, 5000));
    bench("KKT 70000+30000", &kkt_matrix(70000, 30000));
    bench("Tridiagonal 100k", &kkt_matrix(100000, 0));
}
