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
    use crate::numeric::multifrontal_factor;
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
}
