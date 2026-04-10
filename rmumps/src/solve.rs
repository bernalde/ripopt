use crate::numeric::NumericFactorization;
use crate::symbolic::SymbolicFactorization;
use crate::SolverError;

/// Dense GEMV on contiguous arrays: y[0..m] -= A[m x n, col-major] * x[0..n].
/// Column-oriented loop for cache efficiency and auto-vectorization.
#[inline]
fn gemv_sub(m: usize, n: usize, a: &[f64], lda: usize, x: &[f64], y: &mut [f64]) {
    for j in 0..n {
        let xj = x[j];
        if xj == 0.0 { continue; }
        let col = &a[j * lda..j * lda + m];
        for i in 0..m {
            y[i] -= col[i] * xj;
        }
    }
}

/// Dense transposed GEMV on contiguous arrays: y[0..n] -= A[m x n, col-major]^T * x[0..m].
/// Uses Kahan summation for each output element.
#[inline]
fn gemv_t_sub_kahan(m: usize, n: usize, a: &[f64], lda: usize, x: &[f64], y: &mut [f64]) {
    for j in 0..n {
        let col = &a[j * lda..j * lda + m];
        let mut sum = 0.0;
        let mut comp = 0.0;
        for i in 0..m {
            let term = col[i] * x[i] - comp;
            let new_sum = sum + term;
            comp = (new_sum - sum) - term;
            sum = new_sum;
        }
        y[j] -= sum;
    }
}

/// Solve A*x = b using the multifrontal factorization.
///
/// The factorization represents A = product of P*L*D*L^T*P^T blocks
/// organized by the supernodal elimination tree. The solve proceeds:
/// 1. Forward substitution (postorder): apply L^{-1} per supernode + L21 update
/// 2. Diagonal solve: apply D^{-1} per supernode
/// 3. Backward substitution (reverse postorder): apply L^{-T} per supernode
pub fn multifrontal_solve(
    num: &NumericFactorization,
    _sym: &SymbolicFactorization,
    rhs: &[f64],
    solution: &mut [f64],
) -> Result<(), SolverError> {
    let n = num.n;
    if rhs.len() != n || solution.len() != n {
        return Err(SolverError::DimensionMismatch {
            expected: n,
            got: rhs.len().min(solution.len()),
        });
    }

    let mut x = rhs.to_vec();
    let ns = num.num_snodes;

    // Pre-allocate workspace
    let max_nfs = num.node_factors.iter().map(|nf| nf.fs_indices.len()).max().unwrap_or(0);
    let max_ncb = num.node_factors.iter().map(|nf| nf.cb_indices.len()).max().unwrap_or(0);
    let mut fs_vals = vec![0.0; max_nfs];
    let mut cb_vals = vec![0.0; max_ncb];

    // Phase 1: Forward substitution (L^{-1}) + L21 scatter
    for s in 0..ns {
        let nf = &num.node_factors[s];
        let nfs = nf.fs_indices.len();
        let ncb = nf.cb_indices.len();

        if nfs == 1 {
            // Fast path for single-column supernodes (most common)
            let gi = nf.fs_indices[nf.bk.perm[0]];
            let val = x[gi]; // No L11 solve needed for nfs=1
            if ncb > 0 && val != 0.0 {
                let l21_data = &nf.l21.data;
                let cb_indices = &nf.cb_indices;
                for i in 0..ncb {
                    x[cb_indices[i]] -= l21_data[i] * val;
                }
            }
            continue;
        }

        if nfs == 0 { continue; }

        // Gather FS entries and apply BK permutation
        for (k, &perm_idx) in nf.bk.perm.iter().enumerate() {
            fs_vals[k] = x[nf.fs_indices[perm_idx]];
        }

        // Forward solve with L11
        let l_data = &nf.bk.l.data;
        for i in 0..nfs {
            let row = &l_data[i * nfs..i * nfs + i];
            let mut val = fs_vals[i];
            for j in 0..i {
                val -= row[j] * fs_vals[j];
            }
            fs_vals[i] = val;
        }

        // Write back
        for k in 0..nfs {
            x[nf.fs_indices[nf.bk.perm[k]]] = fs_vals[k];
        }

        // Update CB: x[cb] -= L21 * fs_vals (gather/GEMV/scatter)
        if ncb > 0 {
            let l21_data = &nf.l21.data;
            let cb_indices = &nf.cb_indices;
            for i in 0..ncb { cb_vals[i] = x[cb_indices[i]]; }
            gemv_sub(ncb, nfs, l21_data, ncb, &fs_vals[..nfs], &mut cb_vals[..ncb]);
            for i in 0..ncb { x[cb_indices[i]] = cb_vals[i]; }
        }
    }

    // Phase 2: Diagonal solve (D^{-1}) per supernode
    for s in 0..ns {
        let nf = &num.node_factors[s];
        let nfs = nf.fs_indices.len();
        if nfs == 0 { continue; }

        let mut k = 0;
        while k < nfs {
            let gi0 = nf.fs_indices[nf.bk.perm[k]];
            if k + 1 < nfs && nf.bk.d_offdiag[k].abs() > 1e-12 {
                let a = nf.bk.d_diag[k];
                let b = nf.bk.d_offdiag[k];
                let c = nf.bk.d_diag[k + 1];
                let det = a * c - b * b;
                if det.abs() < 1e-30 {
                    return Err(SolverError::SingularMatrix);
                }
                let gi1 = nf.fs_indices[nf.bk.perm[k + 1]];
                let r0 = x[gi0];
                let r1 = x[gi1];
                x[gi0] = (c * r0 - b * r1) / det;
                x[gi1] = (a * r1 - b * r0) / det;
                k += 2;
            } else {
                if nf.bk.d_diag[k].abs() < 1e-30 {
                    return Err(SolverError::SingularMatrix);
                }
                x[gi0] /= nf.bk.d_diag[k];
                k += 1;
            }
        }
    }

    // Phase 3: Backward substitution (reverse postorder = ns-1..0)
    for s in (0..ns).rev() {
        let nf = &num.node_factors[s];
        let nfs = nf.fs_indices.len();
        let ncb = nf.cb_indices.len();

        if nfs == 1 {
            // Fast path for single-column supernodes (gather + Kahan dot)
            let gi = nf.fs_indices[nf.bk.perm[0]];
            if ncb > 0 {
                let l21_data = &nf.l21.data;
                let cb_indices = &nf.cb_indices;
                for i in 0..ncb { cb_vals[i] = x[cb_indices[i]]; }
                let mut sum = 0.0;
                let mut comp = 0.0;
                for i in 0..ncb {
                    let term = l21_data[i] * cb_vals[i] - comp;
                    let new_sum = sum + term;
                    comp = (new_sum - sum) - term;
                    sum = new_sum;
                }
                x[gi] -= sum;
            }
            // No L11^T solve needed for nfs=1
            continue;
        }

        if nfs == 0 { continue; }

        // Gather FS vals
        for k in 0..nfs {
            fs_vals[k] = x[nf.fs_indices[nf.bk.perm[k]]];
        }

        // Update from CB: fs_vals -= L21^T * x[cb] (gather/GEMV^T/Kahan)
        if ncb > 0 {
            let l21_data = &nf.l21.data;
            let cb_indices = &nf.cb_indices;
            for i in 0..ncb { cb_vals[i] = x[cb_indices[i]]; }
            gemv_t_sub_kahan(ncb, nfs, l21_data, ncb, &cb_vals[..ncb], &mut fs_vals[..nfs]);
        }

        // Backward solve with L11^T (column-oriented for cache efficiency)
        let l_data = &nf.bk.l.data;
        for j in (1..nfs).rev() {
            let fj = fs_vals[j];
            if fj == 0.0 { continue; }
            let row = &l_data[j * nfs..j * nfs + j];
            for i in 0..j {
                fs_vals[i] -= row[i] * fj;
            }
        }

        // Write back
        for k in 0..nfs {
            x[nf.fs_indices[nf.bk.perm[k]]] = fs_vals[k];
        }
    }

    solution.copy_from_slice(&x);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::csc::CscMatrix;
    use crate::numeric::{multifrontal_factor, multifrontal_factor_threshold};
    use crate::symbolic::SymbolicFactorization;

    fn factor_and_solve(
        n: usize,
        triplets: &[(usize, usize, f64)],
        b: &[f64],
    ) -> Vec<f64> {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);
        let num = multifrontal_factor(&csc, &sym);

        let mut x = vec![0.0; n];
        multifrontal_solve(&num, &sym, b, &mut x).unwrap();
        x
    }

    fn factor_threshold_and_solve(
        n: usize,
        triplets: &[(usize, usize, f64)],
        b: &[f64],
        threshold: f64,
        n_primal: Option<usize>,
    ) -> Vec<f64> {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);
        let num = multifrontal_factor_threshold(&csc, &sym, threshold, n_primal);

        let mut x = vec![0.0; n];
        multifrontal_solve(&num, &sym, b, &mut x).unwrap();
        x
    }

    fn check_residual(
        n: usize,
        triplets: &[(usize, usize, f64)],
        x: &[f64],
        b: &[f64],
        tol: f64,
    ) {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let mut ax = vec![0.0; n];
        coo.matvec(x, &mut ax).unwrap();

        let mut max_resid = 0.0f64;
        for i in 0..n {
            max_resid = max_resid.max((ax[i] - b[i]).abs());
        }
        assert!(
            max_resid < tol,
            "max residual {} exceeds tolerance {}",
            max_resid, tol
        );
    }

    /// Compute componentwise backward error: max_i |b - Ax|_i / (|A||x|_i + |b|_i).
    fn backward_error(
        n: usize,
        triplets: &[(usize, usize, f64)],
        x: &[f64],
        b: &[f64],
    ) -> f64 {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows.clone(), cols.clone(), vals.clone()).unwrap();
        let mut ax = vec![0.0; n];
        coo.matvec(x, &mut ax).unwrap();

        // Compute |A|*|x| componentwise
        let mut abs_ax = vec![0.0; n];
        for (idx, &v) in vals.iter().enumerate() {
            let r = rows[idx];
            let c = cols[idx];
            abs_ax[r] += v.abs() * x[c].abs();
            if r != c {
                abs_ax[c] += v.abs() * x[r].abs();
            }
        }

        let mut max_berr = 0.0f64;
        for i in 0..n {
            let resid = (b[i] - ax[i]).abs();
            let denom = abs_ax[i] + b[i].abs();
            if denom > 0.0 {
                max_berr = max_berr.max(resid / denom);
            }
        }
        max_berr
    }

    #[test]
    fn test_solve_diagonal() {
        let triplets = [(0, 0, 2.0), (1, 1, 3.0), (2, 2, 5.0)];
        let b = [4.0, 9.0, 25.0];
        let x = factor_and_solve(3, &triplets, &b);
        check_residual(3, &triplets, &x, &b, 1e-10);
        assert!((x[0] - 2.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
        assert!((x[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_tridiagonal_spd() {
        let triplets = [
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0), (2, 3, 1.0),
            (3, 3, 4.0),
        ];
        let b = [1.0, 2.0, 3.0, 4.0];
        let x = factor_and_solve(4, &triplets, &b);
        check_residual(4, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_solve_dense_spd_3x3() {
        let triplets = [
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ];
        let b = [8.0, 18.0, 25.0];
        let x = factor_and_solve(3, &triplets, &b);
        check_residual(3, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_solve_indefinite_kkt_3x3() {
        // [[2, 0, 1], [0, 2, 1], [1, 1, 0]]
        let triplets = [
            (0, 0, 2.0), (0, 2, 1.0),
            (1, 1, 2.0), (1, 2, 1.0),
            (2, 2, 0.0),
        ];
        let b = [3.0, 5.0, 2.0];
        let x = factor_and_solve(3, &triplets, &b);
        check_residual(3, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_solve_5x5_kkt() {
        let triplets = [
            (0, 0, 4.0), (0, 3, 1.0),
            (1, 1, 5.0), (1, 4, 1.0),
            (2, 2, 6.0), (2, 3, 1.0), (2, 4, 1.0),
            (3, 3, 0.0),
            (4, 4, 0.0),
        ];
        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let x = factor_and_solve(5, &triplets, &b);
        check_residual(5, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_solve_arrow_4x4() {
        let triplets = [
            (0, 0, 10.0), (0, 3, 1.0),
            (1, 1, 10.0), (1, 3, 1.0),
            (2, 2, 10.0), (2, 3, 1.0),
            (3, 3, 10.0),
        ];
        let b = [11.0, 12.0, 13.0, 13.0];
        let x = factor_and_solve(4, &triplets, &b);
        check_residual(4, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_solve_multiple_rhs() {
        // Factor once, solve twice with different RHS
        let triplets_data = [
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0),
        ];
        let rows: Vec<usize> = triplets_data.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets_data.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets_data.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(3, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);
        let num = multifrontal_factor(&csc, &sym);

        let b1 = [1.0, 2.0, 3.0];
        let mut x1 = vec![0.0; 3];
        multifrontal_solve(&num, &sym, &b1, &mut x1).unwrap();
        check_residual(3, &triplets_data, &x1, &b1, 1e-10);

        let b2 = [5.0, 10.0, 15.0];
        let mut x2 = vec![0.0; 3];
        multifrontal_solve(&num, &sym, &b2, &mut x2).unwrap();
        check_residual(3, &triplets_data, &x2, &b2, 1e-10);
    }

    // ========== Threshold pivoting solve tests ==========
    // All tests below use multifrontal_factor_threshold with threshold > 0,
    // exercising the MUMPS-style LDLT code path (symmetric swaps, 2x2 pivots,
    // delayed pivots, between-block GEMM) that the above tests skip.

    #[test]
    fn test_threshold_solve_3x3_kkt() {
        // Same 3x3 KKT as above but with threshold pivoting
        // [[2, 0, 1], [0, 2, 1], [1, 1, 0]] — inertia (2, 1, 0)
        let triplets: Vec<(usize, usize, f64)> = vec![
            (0, 0, 2.0), (0, 2, 1.0),
            (1, 1, 2.0), (1, 2, 1.0),
            (2, 2, 0.0),
        ];
        let b = [3.0, 5.0, 2.0];
        let x = factor_threshold_and_solve(3, &triplets, &b, 0.01, None);
        check_residual(3, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_threshold_solve_5x5_kkt() {
        // 5x5 KKT: 3 primal vars + 2 constraints with zero diagonal
        let triplets: Vec<(usize, usize, f64)> = vec![
            (0, 0, 4.0), (0, 3, 1.0),
            (1, 1, 5.0), (1, 4, 1.0),
            (2, 2, 6.0), (2, 3, 1.0), (2, 4, 1.0),
            (3, 3, 0.0),
            (4, 4, 0.0),
        ];
        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let x = factor_threshold_and_solve(5, &triplets, &b, 0.01, Some(3));
        check_residual(5, &triplets, &x, &b, 1e-10);
    }

    #[test]
    fn test_threshold_solve_7x7_kkt_scale_disparity() {
        // 7x7 KKT with scale disparity: H = diag(1e4, 1e-1, 1e2), 2 constraints
        // Forces 2x2 pivots on zero-diagonal rows and tests pivot selection
        // with entries spanning 5 orders of magnitude.
        let n = 5; // 5 primal
        let m = 2; // 2 constraints
        let nm = n + m;
        let triplets: Vec<(usize, usize, f64)> = vec![
            // H block
            (0, 0, 1e4), (1, 1, 1e-1), (2, 2, 1e2),
            (3, 3, 1.0), (4, 4, 10.0),
            // H off-diagonal coupling
            (0, 1, 0.5), (2, 3, 0.3),
            // J^T block (rows 0..5, cols 5..7)
            (0, 5, 1.0), (1, 5, 1.0), (2, 5, 1.0),
            (3, 6, 1.0), (4, 6, 1.0),
            // Zero diagonal for dual block
            (5, 5, 0.0), (6, 6, 0.0),
        ];
        let b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let x = factor_threshold_and_solve(nm, &triplets, &b, 0.01, Some(n));
        check_residual(nm, &triplets, &x, &b, 1e-8);
        let berr = backward_error(nm, &triplets, &x, &b);
        assert!(berr < 1e-10, "backward error {berr:.2e} too large");
    }

    #[test]
    fn test_threshold_solve_10x10_kkt_with_coupling() {
        // 10x10 KKT: 6 primal + 4 constraints, off-diagonal coupling across
        // supernodes to create multi-level elimination tree and delayed pivots.
        let n = 6;
        let m = 4;
        let nm = n + m;
        let triplets: Vec<(usize, usize, f64)> = vec![
            // H block with varying scale
            (0, 0, 1e3), (1, 1, 1e-2), (2, 2, 50.0),
            (3, 3, 1.0), (4, 4, 1e4), (5, 5, 0.1),
            // H coupling (creates fill across supernodes)
            (0, 2, 1.0), (1, 3, 0.5), (4, 5, 2.0),
            // J^T block
            (0, 6, 1.0), (1, 6, 0.3),
            (2, 7, 1.0), (3, 7, 0.7),
            (4, 8, 1.0),
            (5, 9, 1.0), (3, 9, 0.2),
            // Zero dual diagonal
            (6, 6, 0.0), (7, 7, 0.0), (8, 8, 0.0), (9, 9, 0.0),
        ];
        let b: Vec<f64> = (1..=nm).map(|i| i as f64).collect();
        let x = factor_threshold_and_solve(nm, &triplets, &b, 0.01, Some(n));
        check_residual(nm, &triplets, &x, &b, 1e-8);
        let berr = backward_error(nm, &triplets, &x, &b);
        assert!(berr < 1e-10, "backward error {berr:.2e} too large");
    }

    #[test]
    fn test_threshold_solve_high_threshold() {
        // Same 10x10 KKT with threshold 0.1 to force delayed pivots.
        // Tests that delayed pivot promotion + solve is correct.
        let n = 6;
        let m = 4;
        let nm = n + m;
        let triplets: Vec<(usize, usize, f64)> = vec![
            (0, 0, 1e3), (1, 1, 1e-2), (2, 2, 50.0),
            (3, 3, 1.0), (4, 4, 1e4), (5, 5, 0.1),
            (0, 2, 1.0), (1, 3, 0.5), (4, 5, 2.0),
            (0, 6, 1.0), (1, 6, 0.3),
            (2, 7, 1.0), (3, 7, 0.7),
            (4, 8, 1.0),
            (5, 9, 1.0), (3, 9, 0.2),
            (6, 6, 0.0), (7, 7, 0.0), (8, 8, 0.0), (9, 9, 0.0),
        ];
        let b: Vec<f64> = (1..=nm).map(|i| i as f64).collect();
        // Test multiple thresholds: higher values force more delayed pivots
        for &thresh in &[0.01, 0.1, 0.2, 0.5] {
            let x = factor_threshold_and_solve(nm, &triplets, &b, thresh, Some(n));
            let berr = backward_error(nm, &triplets, &x, &b);
            assert!(
                berr < 1e-4,
                "threshold={thresh}: backward error {berr:.2e} too large"
            );
        }
    }

    #[test]
    fn test_threshold_solve_regularized_kkt() {
        // KKT with small negative diagonal on constraint block (as IPM does)
        // Tests that threshold pivoting handles near-quasidefinite structure.
        let n = 4;
        let m = 3;
        let nm = n + m;
        let delta_c = 1e-4;
        let triplets: Vec<(usize, usize, f64)> = vec![
            // H + Sigma block (positive definite)
            (0, 0, 5.0), (1, 1, 3.0), (2, 2, 7.0), (3, 3, 2.0),
            (0, 1, 0.5), (2, 3, 0.3),
            // J^T block
            (0, 4, 1.0), (1, 4, 1.0),
            (2, 5, 1.0), (3, 5, 0.5),
            (1, 6, 0.7), (3, 6, 1.0),
            // Regularized dual diagonal (-delta_c, not zero)
            (4, 4, -delta_c), (5, 5, -delta_c), (6, 6, -delta_c),
        ];
        let b: Vec<f64> = (1..=nm).map(|i| i as f64).collect();
        let x = factor_threshold_and_solve(nm, &triplets, &b, 0.01, Some(n));
        check_residual(nm, &triplets, &x, &b, 1e-8);
        let berr = backward_error(nm, &triplets, &x, &b);
        assert!(berr < 1e-10, "backward error {berr:.2e} too large");
    }

    #[test]
    fn test_threshold_solve_20x20_kkt_grid() {
        // Larger KKT from a discretized PDE-like problem: 12 primal (2D grid)
        // + 8 constraints. Tests multi-supernode tree with threshold pivoting.
        let n = 12;
        let m = 8;
        let nm = n + m;
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();

        // H block: 2D grid Laplacian (4x3) with diagonal dominance
        let nx = 4;
        let ny = 3;
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                triplets.push((i, i, 4.0 + 0.1 * i as f64)); // diagonal
                if ix + 1 < nx { triplets.push((i, i + 1, -1.0)); }
                if iy + 1 < ny { triplets.push((i, i + nx, -1.0)); }
            }
        }

        // J^T block: each constraint touches 2 adjacent primal vars
        for j in 0..m {
            let row1 = j % n;
            let row2 = (j + 1) % n;
            triplets.push((row1, n + j, 1.0));
            triplets.push((row2, n + j, 0.5));
            triplets.push((n + j, n + j, 0.0)); // zero dual diagonal
        }

        let b: Vec<f64> = (0..nm).map(|i| (i as f64 + 1.0).sin()).collect();
        let x = factor_threshold_and_solve(nm, &triplets, &b, 0.01, Some(n));
        check_residual(nm, &triplets, &x, &b, 1e-6);
        let berr = backward_error(nm, &triplets, &x, &b);
        assert!(berr < 1e-8, "backward error {berr:.2e} too large for 20x20 KKT");
    }

    #[test]
    fn test_threshold_vs_nothreshold_consistency() {
        // Verify that threshold=0.0 and threshold=0.01 produce solutions
        // with similar backward error on a well-conditioned system.
        let triplets: Vec<(usize, usize, f64)> = vec![
            (0, 0, 10.0), (0, 1, 1.0), (0, 2, 0.5),
            (1, 1, 10.0), (1, 2, 1.0),
            (2, 2, 10.0),
        ];
        let b = [1.0, 2.0, 3.0];
        let x_no = factor_and_solve(3, &triplets, &b);
        let x_th = factor_threshold_and_solve(3, &triplets, &b, 0.01, None);
        let berr_no = backward_error(3, &triplets, &x_no, &b);
        let berr_th = backward_error(3, &triplets, &x_th, &b);
        assert!(berr_no < 1e-14, "no-threshold berr {berr_no:.2e}");
        assert!(berr_th < 1e-14, "threshold berr {berr_th:.2e}");
    }
}
