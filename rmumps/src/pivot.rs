use crate::dense::DenseMat;
use crate::Inertia;

/// y[i] -= alpha * x[i], using SIMD when available.
#[cfg(target_arch = "aarch64")]
#[inline]
fn axpy_neg(alpha: f64, x: &[f64], y: &mut [f64]) {
    let n = x.len();
    debug_assert!(y.len() >= n);
    unsafe {
        use std::arch::aarch64::*;
        let av = vdupq_n_f64(alpha);
        let mut i = 0;
        while i + 8 <= n {
            let x0 = vld1q_f64(x.as_ptr().add(i));
            let x1 = vld1q_f64(x.as_ptr().add(i + 2));
            let x2 = vld1q_f64(x.as_ptr().add(i + 4));
            let x3 = vld1q_f64(x.as_ptr().add(i + 6));
            let y0 = vld1q_f64(y.as_ptr().add(i));
            let y1 = vld1q_f64(y.as_ptr().add(i + 2));
            let y2 = vld1q_f64(y.as_ptr().add(i + 4));
            let y3 = vld1q_f64(y.as_ptr().add(i + 6));
            vst1q_f64(y.as_mut_ptr().add(i),     vfmsq_f64(y0, x0, av));
            vst1q_f64(y.as_mut_ptr().add(i + 2), vfmsq_f64(y1, x1, av));
            vst1q_f64(y.as_mut_ptr().add(i + 4), vfmsq_f64(y2, x2, av));
            vst1q_f64(y.as_mut_ptr().add(i + 6), vfmsq_f64(y3, x3, av));
            i += 8;
        }
        while i + 2 <= n {
            let xv = vld1q_f64(x.as_ptr().add(i));
            let yv = vld1q_f64(y.as_ptr().add(i));
            vst1q_f64(y.as_mut_ptr().add(i), vfmsq_f64(yv, xv, av));
            i += 2;
        }
        if i < n { y[i] -= alpha * x[i]; }
    }
}

/// y[i] -= a0 * x0[i] + a1 * x1[i], rank-2 AXPY.
#[cfg(target_arch = "aarch64")]
#[inline]
fn axpy2_neg(a0: f64, a1: f64, x0: &[f64], x1: &[f64], y: &mut [f64]) {
    let n = x0.len();
    debug_assert!(x1.len() >= n && y.len() >= n);
    unsafe {
        use std::arch::aarch64::*;
        let av0 = vdupq_n_f64(a0);
        let av1 = vdupq_n_f64(a1);
        let mut i = 0;
        while i + 4 <= n {
            let x0a = vld1q_f64(x0.as_ptr().add(i));
            let x0b = vld1q_f64(x0.as_ptr().add(i + 2));
            let x1a = vld1q_f64(x1.as_ptr().add(i));
            let x1b = vld1q_f64(x1.as_ptr().add(i + 2));
            let ya = vld1q_f64(y.as_ptr().add(i));
            let yb = vld1q_f64(y.as_ptr().add(i + 2));
            let ya = vfmsq_f64(ya, x0a, av0);
            let yb = vfmsq_f64(yb, x0b, av0);
            let ya = vfmsq_f64(ya, x1a, av1);
            let yb = vfmsq_f64(yb, x1b, av1);
            vst1q_f64(y.as_mut_ptr().add(i), ya);
            vst1q_f64(y.as_mut_ptr().add(i + 2), yb);
            i += 4;
        }
        while i < n {
            y[i] -= a0 * x0[i] + a1 * x1[i];
            i += 1;
        }
    }
}

#[cfg(not(target_arch = "aarch64"))]
#[inline]
fn axpy_neg(alpha: f64, x: &[f64], y: &mut [f64]) {
    for i in 0..x.len() {
        y[i] -= alpha * x[i];
    }
}

#[cfg(not(target_arch = "aarch64"))]
#[inline]
fn axpy2_neg(a0: f64, a1: f64, x0: &[f64], x1: &[f64], y: &mut [f64]) {
    for i in 0..x0.len() {
        y[i] -= a0 * x0[i] + a1 * x1[i];
    }
}

/// Result of Bunch-Kaufman LDL^T factorization of a dense symmetric matrix.
#[derive(Debug, Clone)]
pub struct BunchKaufmanResult {
    /// Unit lower triangular factor L (stored in full n x n, row-major layout within DenseMat).
    pub l: DenseMat,
    /// Diagonal of D (1x1 blocks).
    pub d_diag: Vec<f64>,
    /// Off-diagonal of D (2x2 blocks). d_offdiag[k] != 0 means (k, k+1) is a 2x2 block.
    pub d_offdiag: Vec<f64>,
    /// Pivot permutation: original index -> factored position.
    pub perm: Vec<usize>,
    /// Inverse permutation.
    pub perm_inv: Vec<usize>,
    /// Inertia computed from D.
    pub inertia: Inertia,
}

/// Bunch-Kaufman alpha parameter: (1 + sqrt(17)) / 8.
const BK_ALPHA: f64 = 0.6404;

/// Tolerance for zero pivot detection.
const ZERO_PIVOT_TOL: f64 = 1e-12;

/// Find pivot for Bunch-Kaufman algorithm.
/// Returns (pivot_type, p1, p2):
/// - pivot_type=1: 1x1 pivot at row/col p1
/// - pivot_type=2: 2x2 pivot at rows/cols (p1, p2)
fn find_pivot(a: &[f64], n: usize, k: usize) -> (usize, usize, usize) {
    if k == n - 1 {
        return (1, k, k);
    }

    let akk = a[k * n + k].abs();

    // Find largest off-diagonal |a[i][k]| for i > k
    let mut lambda = 0.0f64;
    let mut r = k;
    for i in (k + 1)..n {
        let v = a[i * n + k].abs();
        if v > lambda {
            lambda = v;
            r = i;
        }
    }

    if lambda == 0.0 && akk == 0.0 {
        return (1, k, k); // zero column
    }

    if akk >= BK_ALPHA * lambda {
        return (1, k, k); // 1x1 pivot is good
    }

    // Find largest off-diagonal in row r
    let mut sigma = 0.0f64;
    for j in k..n {
        if j != r {
            let v = a[r * n + j].abs();
            if v > sigma {
                sigma = v;
            }
        }
    }

    if akk * sigma >= BK_ALPHA * lambda * lambda {
        return (1, k, k); // 1x1 pivot at k
    }

    let arr = a[r * n + r].abs();
    if arr >= BK_ALPHA * sigma {
        return (1, r, r); // 1x1 pivot at r
    }

    (2, k, r) // 2x2 pivot
}

/// Swap rows and columns p and q in the active submatrix [start..n, start..n].
/// Columns [0..start) have already been factored and stored in L separately,
/// so their entries in `a` are stale and don't need swapping.
#[inline]
fn swap_rows_cols_from(a: &mut [f64], n: usize, p: usize, q: usize, start: usize) {
    if p == q {
        return;
    }
    // Swap rows p and q (only columns start..n matter)
    for j in start..n {
        a.swap(p * n + j, q * n + j);
    }
    // Swap columns p and q (only rows start..n matter)
    for i in start..n {
        a.swap(i * n + p, i * n + q);
    }
}

/// Compute inertia from D diagonal and off-diagonal entries.
pub fn compute_inertia(d_diag: &[f64], d_offdiag: &[f64], n: usize) -> Inertia {
    let mut positive = 0;
    let mut negative = 0;
    let mut zero = 0;

    let mut k = 0;
    while k < n {
        if k + 1 < n && d_offdiag[k].abs() > ZERO_PIVOT_TOL {
            // 2x2 block eigenvalues
            let a = d_diag[k];
            let b = d_offdiag[k];
            let c = d_diag[k + 1];
            let trace = a + c;
            let det = a * c - b * b;
            let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
            let eig1 = (trace + disc) / 2.0;
            let eig2 = (trace - disc) / 2.0;

            for eig in [eig1, eig2] {
                if eig > ZERO_PIVOT_TOL {
                    positive += 1;
                } else if eig < -ZERO_PIVOT_TOL {
                    negative += 1;
                } else {
                    zero += 1;
                }
            }
            k += 2;
        } else {
            let d = d_diag[k];
            if d > ZERO_PIVOT_TOL {
                positive += 1;
            } else if d < -ZERO_PIVOT_TOL {
                negative += 1;
            } else {
                zero += 1;
            }
            k += 1;
        }
    }

    Inertia { positive, negative, zero }
}

/// Bunch-Kaufman LDL^T factorization of a dense symmetric indefinite matrix.
/// Input: symmetric matrix in `a` (full storage, row-major within DenseMat n x n).
/// The input matrix is consumed (overwritten).
pub fn dense_ldlt_bunch_kaufman(a: &mut DenseMat) -> BunchKaufmanResult {
    let n = a.nrows;
    debug_assert_eq!(a.ncols, n);

    let mut l = DenseMat::zeros(n, n);
    let mut d_diag = vec![0.0; n];
    let mut d_offdiag = vec![0.0; n];
    let mut perm: Vec<usize> = (0..n).collect();
    let mut work = vec![0.0f64; 2 * n];

    let aa = &mut a.data;

    let mut k = 0;
    while k < n {
        let (pivot_type, p1, p2) = find_pivot(aa, n, k);

        if pivot_type == 1 {
            if p1 != k {
                // Swap only the active part of aa (columns k..n)
                swap_rows_cols_from(aa, n, k, p1, k);
                perm.swap(k, p1);
                // Swap L entries for previously computed columns
                for j in 0..k {
                    l.data.swap(k * n + j, p1 * n + j);
                }
            }

            let akk = aa[k * n + k];
            d_diag[k] = akk;

            if akk.abs() > ZERO_PIVOT_TOL {
                let m = n - k - 1;
                for i in 0..m {
                    work[i] = aa[(k + 1 + i) * n + k] / akk;
                    l.data[(k + 1 + i) * n + k] = work[i];
                }
                for i in 0..m {
                    let si = work[i] * akk;
                    let base = (k + 1 + i) * n + (k + 1);
                    axpy_neg(si, &work[..m], &mut aa[base..base + m]);
                }
            }
            l.data[k * n + k] = 1.0;
            k += 1;
        } else {
            if p2 != k + 1 {
                swap_rows_cols_from(aa, n, k + 1, p2, k);
                perm.swap(k + 1, p2);
                for j in 0..k {
                    l.data.swap((k + 1) * n + j, p2 * n + j);
                }
            }
            if p1 != k {
                swap_rows_cols_from(aa, n, k, p1, k);
                perm.swap(k, p1);
                for j in 0..k {
                    l.data.swap(k * n + j, p1 * n + j);
                }
            }

            let akk = aa[k * n + k];
            let ak1k = aa[(k + 1) * n + k];
            let ak1k1 = aa[(k + 1) * n + (k + 1)];

            d_diag[k] = akk;
            d_diag[k + 1] = ak1k1;
            d_offdiag[k] = ak1k;

            let det = akk * ak1k1 - ak1k * ak1k;

            if det.abs() > ZERO_PIVOT_TOL {
                let d_inv_00 = ak1k1 / det;
                let d_inv_01 = -ak1k / det;
                let d_inv_11 = akk / det;

                let m = n - k - 2;
                for i in 0..m {
                    let aik = aa[(k + 2 + i) * n + k];
                    let aik1 = aa[(k + 2 + i) * n + (k + 1)];
                    work[i] = aik * d_inv_00 + aik1 * d_inv_01;
                    work[m + i] = aik * d_inv_01 + aik1 * d_inv_11;
                    l.data[(k + 2 + i) * n + k] = work[i];
                    l.data[(k + 2 + i) * n + (k + 1)] = work[m + i];
                }

                for i in 0..m {
                    let li0 = work[i];
                    let li1 = work[m + i];
                    let si0 = li0 * akk + li1 * ak1k;
                    let si1 = li0 * ak1k + li1 * ak1k1;
                    let base = (k + 2 + i) * n + (k + 2);
                    axpy2_neg(si0, si1, &work[..m], &work[m..m + m], &mut aa[base..base + m]);
                }
            }

            l.data[k * n + k] = 1.0;
            l.data[(k + 1) * n + (k + 1)] = 1.0;
            k += 2;
        }
    }

    let mut perm_inv = vec![0; n];
    for i in 0..n {
        perm_inv[perm[i]] = i;
    }

    let inertia = compute_inertia(&d_diag, &d_offdiag, n);

    BunchKaufmanResult { l, d_diag, d_offdiag, perm, perm_inv, inertia }
}

/// Solve A*x = b given the Bunch-Kaufman factorization A = P*L*D*L^T*P^T.
/// rhs is the right-hand side b, solution receives x.
pub fn bunch_kaufman_solve(bk: &BunchKaufmanResult, rhs: &[f64], solution: &mut [f64]) {
    let n = bk.l.nrows;

    // Step 1: y = P * b
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = rhs[bk.perm[i]];
    }

    // Step 2: Forward substitution L * z = y
    for i in 0..n {
        for j in 0..i {
            y[i] -= bk.l.data[i * n + j] * y[j];
        }
    }

    // Step 3: Solve D * w = z
    let mut w = vec![0.0; n];
    let mut k = 0;
    while k < n {
        if k + 1 < n && bk.d_offdiag[k].abs() > ZERO_PIVOT_TOL {
            let a = bk.d_diag[k];
            let b = bk.d_offdiag[k];
            let c = bk.d_diag[k + 1];
            let det = a * c - b * b;
            w[k] = (c * y[k] - b * y[k + 1]) / det;
            w[k + 1] = (a * y[k + 1] - b * y[k]) / det;
            k += 2;
        } else {
            w[k] = y[k] / bk.d_diag[k];
            k += 1;
        }
    }

    // Step 4: Backward substitution L^T * v = w
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            w[i] -= bk.l.data[j * n + i] * w[j];
        }
    }

    // Step 5: solution = P^T * v
    for i in 0..n {
        solution[bk.perm[i]] = w[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_full_symmetric(vals: &[&[f64]]) -> DenseMat {
        let n = vals.len();
        let mut m = DenseMat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                m.set(i, j, vals[i][j]);
            }
        }
        m
    }

    fn verify_factorization(orig: &DenseMat, bk: &BunchKaufmanResult) {
        let n = orig.nrows;
        // Reconstruct P*L*D*L^T*P^T and compare to original
        // First compute L*D*L^T
        let mut ldlt = DenseMat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut val = 0.0;
                let mut kk = 0;
                while kk < n {
                    if kk + 1 < n && bk.d_offdiag[kk].abs() > ZERO_PIVOT_TOL {
                        // 2x2 block
                        let lik0 = bk.l.data[i * n + kk];
                        let lik1 = bk.l.data[i * n + kk + 1];
                        let ljk0 = bk.l.data[j * n + kk];
                        let ljk1 = bk.l.data[j * n + kk + 1];
                        let d00 = bk.d_diag[kk];
                        let d01 = bk.d_offdiag[kk];
                        let d11 = bk.d_diag[kk + 1];
                        // (L*D)_{i,kk} = lik0*d00 + lik1*d01
                        // (L*D)_{i,kk+1} = lik0*d01 + lik1*d11
                        val += (lik0 * d00 + lik1 * d01) * ljk0
                            + (lik0 * d01 + lik1 * d11) * ljk1;
                        kk += 2;
                    } else {
                        let lik = bk.l.data[i * n + kk];
                        let ljk = bk.l.data[j * n + kk];
                        val += lik * bk.d_diag[kk] * ljk;
                        kk += 1;
                    }
                }
                ldlt.set(i, j, val);
            }
        }

        // Apply P: A_orig[perm[i], perm[j]] should == ldlt[i,j]
        for i in 0..n {
            for j in 0..n {
                let expected = orig.get(bk.perm[i], bk.perm[j]);
                let got = ldlt.get(i, j);
                assert!(
                    (expected - got).abs() < 1e-10,
                    "P*L*D*L^T*P^T mismatch at ({},{}): expected {} got {}",
                    i, j, expected, got
                );
            }
        }
    }

    #[test]
    fn test_bk_spd_3x3() {
        let mut a = make_full_symmetric(&[
            &[4.0, 2.0, 1.0],
            &[2.0, 5.0, 3.0],
            &[1.0, 3.0, 6.0],
        ]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);
        assert_eq!(bk.inertia, Inertia { positive: 3, negative: 0, zero: 0 });
        verify_factorization(&orig, &bk);
    }

    #[test]
    fn test_bk_indefinite_2x2() {
        // [[1, 2], [2, 1]] — eigenvalues 3 and -1
        let mut a = make_full_symmetric(&[&[1.0, 2.0], &[2.0, 1.0]]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);
        assert_eq!(bk.inertia, Inertia { positive: 1, negative: 1, zero: 0 });
        verify_factorization(&orig, &bk);
    }

    #[test]
    fn test_bk_kkt_like() {
        // KKT system: [[H, A^T], [A, 0]]
        // H = [[2, 0], [0, 2]], A = [[1, 1]]
        // Full: [[2, 0, 1], [0, 2, 1], [1, 1, 0]]
        // Expected inertia: 2 positive, 1 negative
        let mut a = make_full_symmetric(&[
            &[2.0, 0.0, 1.0],
            &[0.0, 2.0, 1.0],
            &[1.0, 1.0, 0.0],
        ]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);
        assert_eq!(bk.inertia.positive, 2);
        assert_eq!(bk.inertia.negative, 1);
        assert_eq!(bk.inertia.zero, 0);
        verify_factorization(&orig, &bk);
    }

    #[test]
    fn test_bk_solve_spd() {
        let mut a = make_full_symmetric(&[
            &[4.0, 2.0, 1.0],
            &[2.0, 5.0, 3.0],
            &[1.0, 3.0, 6.0],
        ]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);

        let b = [8.0, 18.0, 25.0];
        let mut x = [0.0; 3];
        bunch_kaufman_solve(&bk, &b, &mut x);

        // Verify A*x = b
        for i in 0..3 {
            let mut ax = 0.0;
            for j in 0..3 {
                ax += orig.get(i, j) * x[j];
            }
            assert!(
                (ax - b[i]).abs() < 1e-10,
                "residual at {}: {}",
                i,
                (ax - b[i]).abs()
            );
        }
    }

    #[test]
    fn test_bk_solve_indefinite() {
        // KKT: [[2, 0, 1], [0, 2, 1], [1, 1, 0]]
        let mut a = make_full_symmetric(&[
            &[2.0, 0.0, 1.0],
            &[0.0, 2.0, 1.0],
            &[1.0, 1.0, 0.0],
        ]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);

        let b = [3.0, 5.0, 2.0];
        let mut x = [0.0; 3];
        bunch_kaufman_solve(&bk, &b, &mut x);

        for i in 0..3 {
            let mut ax = 0.0;
            for j in 0..3 {
                ax += orig.get(i, j) * x[j];
            }
            assert!(
                (ax - b[i]).abs() < 1e-10,
                "residual at {}: {}",
                i,
                (ax - b[i]).abs()
            );
        }
    }

    #[test]
    fn test_bk_solve_larger_kkt() {
        // 5x5 KKT: H = 3x3 diagonal, A = 2x3
        // [[4,0,0,1,0], [0,5,0,0,1], [0,0,6,1,1], [1,0,1,0,0], [0,1,1,0,0]]
        let mut a = make_full_symmetric(&[
            &[4.0, 0.0, 0.0, 1.0, 0.0],
            &[0.0, 5.0, 0.0, 0.0, 1.0],
            &[0.0, 0.0, 6.0, 1.0, 1.0],
            &[1.0, 0.0, 1.0, 0.0, 0.0],
            &[0.0, 1.0, 1.0, 0.0, 0.0],
        ]);
        let orig = a.clone();
        let bk = dense_ldlt_bunch_kaufman(&mut a);

        assert_eq!(bk.inertia.positive, 3);
        assert_eq!(bk.inertia.negative, 2);

        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = [0.0; 5];
        bunch_kaufman_solve(&bk, &b, &mut x);

        for i in 0..5 {
            let mut ax = 0.0;
            for j in 0..5 {
                ax += orig.get(i, j) * x[j];
            }
            assert!(
                (ax - b[i]).abs() < 1e-10,
                "residual at {}: {}",
                i,
                (ax - b[i]).abs()
            );
        }
    }
}
