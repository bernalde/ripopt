/// Dense matrix in column-major storage for frontal matrix computations.
/// Entry (i, j) is at index j * nrows + i.
#[derive(Debug, Clone)]
pub struct DenseMat {
    pub nrows: usize,
    pub ncols: usize,
    /// Column-major data, length nrows * ncols.
    pub data: Vec<f64>,
}

impl DenseMat {
    /// Create a zero matrix.
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            data: vec![0.0; nrows * ncols],
        }
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.nrows + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.data[j * self.nrows + i] = val;
    }

    #[inline]
    pub fn add(&mut self, i: usize, j: usize, val: f64) {
        self.data[j * self.nrows + i] += val;
    }
}

/// In-place LDL^T factorization of an n x n symmetric matrix (diagonal pivoting only).
/// The matrix is stored in the top-left n x n block of `mat` (column-major).
/// On exit, L is unit lower triangular in the lower triangle, D is on the diagonal.
/// Returns the diagonal entries of D (length n).
///
/// Only suitable for symmetric positive definite matrices.
/// For indefinite matrices, use `dense_ldlt_bunch_kaufman` (step 4).
pub fn dense_ldlt_factor(mat: &mut DenseMat) -> Vec<f64> {
    let n = mat.nrows;
    debug_assert_eq!(mat.ncols, n);
    let mut d = vec![0.0; n];

    for j in 0..n {
        // d[j] = a[j][j] - sum_{k<j} L[j][k]^2 * d[k]
        let mut dj = mat.get(j, j);
        for k in 0..j {
            let ljk = mat.get(j, k);
            dj -= ljk * ljk * d[k];
        }
        d[j] = dj;

        if dj.abs() < 1e-30 {
            continue; // zero pivot, skip
        }

        // L[i][j] = (a[i][j] - sum_{k<j} L[i][k]*L[j][k]*d[k]) / d[j]
        for i in (j + 1)..n {
            let mut lij = mat.get(i, j);
            for k in 0..j {
                lij -= mat.get(i, k) * mat.get(j, k) * d[k];
            }
            mat.set(i, j, lij / dj);
        }
    }

    d
}

/// Forward solve: solve L * x = b where L is unit lower triangular stored in `mat`.
/// `x` is modified in-place (input as b, output as x).
/// Uses column-oriented access for cache efficiency with column-major storage.
pub fn dense_forward_solve(mat: &DenseMat, n: usize, x: &mut [f64]) {
    let data = &mat.data;
    let nrows = mat.nrows;
    for j in 0..n {
        let xj = x[j];
        let col = &data[j * nrows..(j * nrows + n)];
        // Contiguous read of column j, auto-vectorizable
        for i in (j + 1)..n {
            x[i] -= col[i] * xj;
        }
    }
}

/// Backward solve: solve L^T * x = b where L is unit lower triangular stored in `mat`.
/// `x` is modified in-place.
/// Uses column-oriented access for cache efficiency.
pub fn dense_backward_solve(mat: &DenseMat, n: usize, x: &mut [f64]) {
    let data = &mat.data;
    let nrows = mat.nrows;
    for j in (0..n).rev() {
        let col = &data[j * nrows..(j * nrows + n)];
        let mut sum = 0.0;
        for i in (j + 1)..n {
            sum += col[i] * x[i];
        }
        x[j] -= sum;
    }
}

/// Block size for cache-friendly tiling.
/// Chosen so that 3 blocks of BLOCK x BLOCK f64 fit in L1 cache (~32KB).
/// 3 * 64 * 64 * 8 = 98KB fits comfortably in L2; for L1 use 48.
/// 64 is a good balance for L1/L2 and SIMD alignment.
const BLOCK: usize = 64;

/// Dense matrix multiply: C += alpha * A * B
/// A is m x k, B is k x n, C is m x n. All column-major.
///
/// Uses cache-blocked algorithm with contiguous inner loops for auto-vectorization.
pub fn dense_gemm(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &DenseMat,
    b: &DenseMat,
    c: &mut DenseMat,
) {
    let a_data = &a.data;
    let b_data = &b.data;
    let c_data = &mut c.data;
    let a_rows = a.nrows;
    let b_rows = b.nrows;
    let c_rows = c.nrows;

    // Blocked GEMM: tile over j, kk, i in blocks
    let mut jj = 0;
    while jj < n {
        let j_end = (jj + BLOCK).min(n);
        let mut kk = 0;
        while kk < k {
            let k_end = (kk + BLOCK).min(k);
            let mut ii = 0;
            while ii < m {
                let i_end = (ii + BLOCK).min(m);

                // Micro-kernel: C[ii..i_end, jj..j_end] += alpha * A[ii..i_end, kk..k_end] * B[kk..k_end, jj..j_end]
                for j in jj..j_end {
                    for p in kk..k_end {
                        let bpj = alpha * b_data[j * b_rows + p];
                        let c_col = &mut c_data[j * c_rows + ii..j * c_rows + i_end];
                        let a_col = &a_data[p * a_rows + ii..p * a_rows + i_end];
                        // Contiguous inner loop — LLVM will auto-vectorize this
                        for idx in 0..c_col.len() {
                            c_col[idx] += a_col[idx] * bpj;
                        }
                    }
                }

                ii += BLOCK;
            }
            kk += BLOCK;
        }
        jj += BLOCK;
    }
}

/// Micro-kernel tile size for SIMD GEMM.
/// 8×4 on aarch64 (uses 24 of 32 NEON registers), 4×4 on x86_64.
#[cfg(target_arch = "aarch64")]
const MR: usize = 8;
#[cfg(not(target_arch = "aarch64"))]
const MR: usize = 4;
const NR: usize = 4;

/// Cache-blocked C -= A * B^T where A is m×k, B is n×k (both column-major).
/// C is m×n column-major. Used for Schur complement: S -= W * L21^T.
///
/// For large enough matrices, uses packed data + SIMD micro-kernels.
pub fn gemm_nt_sub(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    // For small matrices, use simple scalar code (packing overhead not worthwhile)
    if m < MR * 2 || n < NR * 2 || k < 4 {
        gemm_nt_sub_scalar(m, n, k, a, lda, b, ldb, c, ldc);
        return;
    }

    // Pack A and B for contiguous micro-kernel access
    let m_padded = (m + MR - 1) / MR * MR;
    let n_padded = (n + NR - 1) / NR * NR;
    let mut packed_a = vec![0.0f64; m_padded * k];
    let mut packed_b = vec![0.0f64; n_padded * k];

    // Pack A: for each MR-wide panel, layout packed_a[panel * MR * k + p * MR + i]
    for ii in (0..m).step_by(MR) {
        let panel = ii / MR;
        let ib = MR.min(m - ii);
        let base = panel * MR * k;
        for p in 0..k {
            for i in 0..ib {
                packed_a[base + p * MR + i] = a[p * lda + ii + i];
            }
            // Remaining entries already zero from allocation
        }
    }

    // Pack B: for each NR-wide panel, layout packed_b[panel * NR * k + p * NR + j]
    for jj in (0..n).step_by(NR) {
        let panel = jj / NR;
        let jb = NR.min(n - jj);
        let base = panel * NR * k;
        for p in 0..k {
            for j in 0..jb {
                packed_b[base + p * NR + j] = b[p * ldb + jj + j];
            }
        }
    }

    // GEMM using micro-kernels on packed data
    for jj in (0..n).step_by(NR) {
        let jb = NR.min(n - jj);
        let b_panel = (jj / NR) * NR * k;
        for ii in (0..m).step_by(MR) {
            let ib = MR.min(m - ii);
            let a_panel = (ii / MR) * MR * k;

            if ib == MR && jb == NR {
                // Full MR×NR micro-kernel
                unsafe {
                    microkernel_nt_sub(
                        k,
                        packed_a.as_ptr().add(a_panel),
                        packed_b.as_ptr().add(b_panel),
                        c.as_mut_ptr().add(jj * ldc + ii),
                        ldc,
                    );
                }
            } else {
                // Edge tile: scalar fallback
                for j in 0..jb {
                    for p in 0..k {
                        let bjp = packed_b[b_panel + p * NR + j];
                        for i in 0..ib {
                            c[(jj + j) * ldc + ii + i] -= packed_a[a_panel + p * MR + i] * bjp;
                        }
                    }
                }
            }
        }
    }
}

/// Scalar fallback for small GEMM-NT.
fn gemm_nt_sub_scalar(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) {
    for j in 0..n {
        for p in 0..k {
            let bjp = b[p * ldb + j];
            let c_col = &mut c[j * ldc..j * ldc + m];
            let a_col = &a[p * lda..p * lda + m];
            for i in 0..m {
                c_col[i] -= a_col[i] * bjp;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SIMD micro-kernels: C[0..MR, 0..NR] -= packed_a * packed_b^T
// packed_a layout: a[p * MR + i] for p in 0..k, i in 0..MR
// packed_b layout: b[p * NR + j] for p in 0..k, j in 0..NR
// ---------------------------------------------------------------------------

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn microkernel_nt_sub(
    k: usize,
    packed_a: *const f64,
    packed_b: *const f64,
    c: *mut f64,
    ldc: usize,
) {
    use std::arch::aarch64::*;

    // 16 accumulator registers for 8×4 tile (each 2×f64)
    // Row groups: [0..2], [2..4], [4..6], [6..8] × 4 columns
    let mut c00 = vdupq_n_f64(0.0); let mut c02 = vdupq_n_f64(0.0);
    let mut c04 = vdupq_n_f64(0.0); let mut c06 = vdupq_n_f64(0.0);
    let mut c10 = vdupq_n_f64(0.0); let mut c12 = vdupq_n_f64(0.0);
    let mut c14 = vdupq_n_f64(0.0); let mut c16 = vdupq_n_f64(0.0);
    let mut c20 = vdupq_n_f64(0.0); let mut c22 = vdupq_n_f64(0.0);
    let mut c24 = vdupq_n_f64(0.0); let mut c26 = vdupq_n_f64(0.0);
    let mut c30 = vdupq_n_f64(0.0); let mut c32 = vdupq_n_f64(0.0);
    let mut c34 = vdupq_n_f64(0.0); let mut c36 = vdupq_n_f64(0.0);

    for p in 0..k {
        let ap = packed_a.add(p * 8);
        let bp = packed_b.add(p * 4);

        // Load 8 elements of A (4 NEON registers)
        let a01 = vld1q_f64(ap);
        let a23 = vld1q_f64(ap.add(2));
        let a45 = vld1q_f64(ap.add(4));
        let a67 = vld1q_f64(ap.add(6));

        // Broadcast 4 elements of B
        let b0 = vdupq_n_f64(*bp);
        let b1 = vdupq_n_f64(*bp.add(1));
        let b2 = vdupq_n_f64(*bp.add(2));
        let b3 = vdupq_n_f64(*bp.add(3));

        // 16 FMAs (8 rows × 4 cols, 2 f64 per FMA = 32 FMAs total)
        c00 = vfmaq_f64(c00, a01, b0); c02 = vfmaq_f64(c02, a23, b0);
        c04 = vfmaq_f64(c04, a45, b0); c06 = vfmaq_f64(c06, a67, b0);
        c10 = vfmaq_f64(c10, a01, b1); c12 = vfmaq_f64(c12, a23, b1);
        c14 = vfmaq_f64(c14, a45, b1); c16 = vfmaq_f64(c16, a67, b1);
        c20 = vfmaq_f64(c20, a01, b2); c22 = vfmaq_f64(c22, a23, b2);
        c24 = vfmaq_f64(c24, a45, b2); c26 = vfmaq_f64(c26, a67, b2);
        c30 = vfmaq_f64(c30, a01, b3); c32 = vfmaq_f64(c32, a23, b3);
        c34 = vfmaq_f64(c34, a45, b3); c36 = vfmaq_f64(c36, a67, b3);
    }

    // Store: C -= accumulated products
    macro_rules! store_col {
        ($col:expr, $r01:expr, $r23:expr, $r45:expr, $r67:expr) => {{
            let ptr = c.add($col * ldc);
            vst1q_f64(ptr,        vsubq_f64(vld1q_f64(ptr),        $r01));
            vst1q_f64(ptr.add(2), vsubq_f64(vld1q_f64(ptr.add(2)), $r23));
            vst1q_f64(ptr.add(4), vsubq_f64(vld1q_f64(ptr.add(4)), $r45));
            vst1q_f64(ptr.add(6), vsubq_f64(vld1q_f64(ptr.add(6)), $r67));
        }};
    }
    store_col!(0, c00, c02, c04, c06);
    store_col!(1, c10, c12, c14, c16);
    store_col!(2, c20, c22, c24, c26);
    store_col!(3, c30, c32, c34, c36);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn microkernel_nt_sub(
    k: usize,
    packed_a: *const f64,
    packed_b: *const f64,
    c: *mut f64,
    ldc: usize,
) {
    use std::arch::x86_64::*;

    // AVX: 256-bit = 4×f64, so 4×4 tile needs 4 C registers
    let mut c0 = _mm256_setzero_pd();
    let mut c1 = _mm256_setzero_pd();
    let mut c2 = _mm256_setzero_pd();
    let mut c3 = _mm256_setzero_pd();

    for p in 0..k {
        let ap = packed_a.add(p * MR);
        let bp = packed_b.add(p * NR);

        let a_vec = _mm256_loadu_pd(ap);

        let b0 = _mm256_broadcast_sd(&*bp);
        let b1 = _mm256_broadcast_sd(&*bp.add(1));
        let b2 = _mm256_broadcast_sd(&*bp.add(2));
        let b3 = _mm256_broadcast_sd(&*bp.add(3));

        c0 = _mm256_fmadd_pd(a_vec, b0, c0);
        c1 = _mm256_fmadd_pd(a_vec, b1, c1);
        c2 = _mm256_fmadd_pd(a_vec, b2, c2);
        c3 = _mm256_fmadd_pd(a_vec, b3, c3);
    }

    // Store: C -= accumulated
    macro_rules! store_col {
        ($col:expr, $acc:expr) => {{
            let ptr = c.add($col * ldc);
            let cur = _mm256_loadu_pd(ptr);
            _mm256_storeu_pd(ptr, _mm256_sub_pd(cur, $acc));
        }};
    }
    store_col!(0, c0);
    store_col!(1, c1);
    store_col!(2, c2);
    store_col!(3, c3);
}

// Fallback for other architectures
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
unsafe fn microkernel_nt_sub(
    k: usize,
    packed_a: *const f64,
    packed_b: *const f64,
    c: *mut f64,
    ldc: usize,
) {
    for p in 0..k {
        for j in 0..NR {
            let bjp = *packed_b.add(p * NR + j);
            for i in 0..MR {
                *c.add(j * ldc + i) -= *packed_a.add(p * MR + i) * bjp;
            }
        }
    }
}

/// Symmetric rank-k update: C -= L21 * D11 * L21^T
/// where L21 is (m x k) and D11 is diagonal (length k).
/// C is m x m symmetric; both triangles are updated.
///
/// Uses a cache-friendly column-oriented approach: for each column p of L21,
/// scale by D[p] and do a rank-1 symmetric update.
pub fn dense_schur_complement(
    m: usize,
    k: usize,
    l21: &DenseMat,
    d11: &[f64],
    c: &mut DenseMat,
) {
    let l_data = &l21.data;
    let c_data = &mut c.data;
    let l_rows = l21.nrows;
    let c_rows = c.nrows;

    // Process in blocks of columns of L21 for cache locality
    let mut pp = 0;
    while pp < k {
        let p_end = (pp + BLOCK).min(k);

        for j in 0..m {
            for p in pp..p_end {
                let scaled_lj = d11[p] * l_data[p * l_rows + j];
                // Update column j of C: C[i,j] -= L[i,p] * d[p] * L[j,p] for i >= j
                let c_col = &mut c_data[j * c_rows..j * c_rows + m];
                let l_col = &l_data[p * l_rows..p * l_rows + m];
                // Lower triangle: i >= j
                for i in j..m {
                    c_col[i] -= l_col[i] * scaled_lj;
                }
            }
        }

        // Mirror lower to upper for full symmetric storage
        pp += BLOCK;
    }

    // Symmetrize: copy lower to upper
    for j in 0..m {
        for i in (j + 1)..m {
            let val = c_data[j * c_rows + i];
            c_data[i * c_rows + j] = val;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_symmetric(vals: &[&[f64]]) -> DenseMat {
        let n = vals.len();
        let mut m = DenseMat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                m.set(i, j, vals[i][j]);
            }
        }
        m
    }

    #[test]
    fn test_ldlt_3x3_spd() {
        // A = [[4, 2, 0], [2, 5, 1], [0, 1, 3]]
        let mut mat = make_symmetric(&[&[4.0, 2.0, 0.0], &[2.0, 5.0, 1.0], &[0.0, 1.0, 3.0]]);
        let orig = mat.clone();
        let d = dense_ldlt_factor(&mut mat);

        // Verify L*D*L^T = A
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                let mut val = 0.0;
                for k in 0..n {
                    let lik = if k == i { 1.0 } else if i > k { mat.get(i, k) } else { 0.0 };
                    let ljk = if k == j { 1.0 } else if j > k { mat.get(j, k) } else { 0.0 };
                    val += lik * d[k] * ljk;
                }
                assert!(
                    (val - orig.get(i, j)).abs() < 1e-12,
                    "L*D*L^T mismatch at ({},{}): {} vs {}",
                    i, j, val, orig.get(i, j)
                );
            }
        }
    }

    #[test]
    fn test_ldlt_solve_3x3() {
        // A = [[4, 2, 0], [2, 5, 1], [0, 1, 3]], b = [8, 13, 7]
        // Solution: x = [1, 2, 1+2/3] ... let's just verify A*x = b
        let mut mat = make_symmetric(&[&[4.0, 2.0, 0.0], &[2.0, 5.0, 1.0], &[0.0, 1.0, 3.0]]);
        let orig = mat.clone();
        let d = dense_ldlt_factor(&mut mat);

        let b = [8.0, 13.0, 7.0];
        let mut x = b;
        // Solve L*D*L^T*x = b: forward (L), diagonal (D), backward (L^T)
        dense_forward_solve(&mat, 3, &mut x);
        for i in 0..3 {
            x[i] /= d[i];
        }
        dense_backward_solve(&mat, 3, &mut x);

        // Verify A*x = b
        for i in 0..3 {
            let mut ax = 0.0;
            for j in 0..3 {
                ax += orig.get(i, j) * x[j];
            }
            assert!(
                (ax - b[i]).abs() < 1e-10,
                "residual at {}: |Ax-b| = {}",
                i,
                (ax - b[i]).abs()
            );
        }
    }

    #[test]
    fn test_ldlt_5x5_spd() {
        // Diagonally dominant 5x5
        let mut mat = DenseMat::zeros(5, 5);
        for i in 0..5 {
            mat.set(i, i, 10.0);
            if i + 1 < 5 {
                mat.set(i, i + 1, 1.0);
                mat.set(i + 1, i, 1.0);
            }
        }
        let orig = mat.clone();
        let d = dense_ldlt_factor(&mut mat);

        // All D entries should be positive (SPD)
        for &di in &d {
            assert!(di > 0.0, "D entry {} should be positive", di);
        }

        // Solve and verify
        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = b;
        dense_forward_solve(&mat, 5, &mut x);
        for i in 0..5 {
            x[i] /= d[i];
        }
        dense_backward_solve(&mat, 5, &mut x);

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

    #[test]
    fn test_gemm() {
        // A = [[1, 2], [3, 4], [5, 6]] (3x2)
        // B = [[7, 8, 9], [10, 11, 12]] (2x3)
        // C = A*B = [[27, 30, 33], [61, 68, 75], [95, 106, 117]]
        let mut a = DenseMat::zeros(3, 2);
        a.set(0, 0, 1.0); a.set(1, 0, 3.0); a.set(2, 0, 5.0);
        a.set(0, 1, 2.0); a.set(1, 1, 4.0); a.set(2, 1, 6.0);

        let mut b = DenseMat::zeros(2, 3);
        b.set(0, 0, 7.0); b.set(0, 1, 8.0); b.set(0, 2, 9.0);
        b.set(1, 0, 10.0); b.set(1, 1, 11.0); b.set(1, 2, 12.0);

        let mut c = DenseMat::zeros(3, 3);
        dense_gemm(3, 3, 2, 1.0, &a, &b, &mut c);

        let expected = [[27.0, 30.0, 33.0], [61.0, 68.0, 75.0], [95.0, 106.0, 117.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (c.get(i, j) - expected[i][j]).abs() < 1e-12,
                    "GEMM mismatch at ({},{})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_schur_complement() {
        // L21 = [[2], [3]], D11 = [4], C_orig = [[10, 5], [5, 20]]
        // Schur = C - L21 * D * L21^T = [[10-16, 5-24], [5-24, 20-36]] = [[-6, -19], [-19, -16]]
        let mut l21 = DenseMat::zeros(2, 1);
        l21.set(0, 0, 2.0);
        l21.set(1, 0, 3.0);
        let d = [4.0];

        let mut c = DenseMat::zeros(2, 2);
        c.set(0, 0, 10.0);
        c.set(0, 1, 5.0);
        c.set(1, 0, 5.0);
        c.set(1, 1, 20.0);

        dense_schur_complement(2, 1, &l21, &d, &mut c);
        assert!((c.get(0, 0) - (-6.0)).abs() < 1e-12);
        assert!((c.get(0, 1) - (-19.0)).abs() < 1e-12);
        assert!((c.get(1, 0) - (-19.0)).abs() < 1e-12);
        assert!((c.get(1, 1) - (-16.0)).abs() < 1e-12);
    }
}
