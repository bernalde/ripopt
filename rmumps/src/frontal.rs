use crate::dense::{DenseMat, gemm_nt_sub};
use crate::pivot::{dense_ldlt_bunch_kaufman, BunchKaufmanResult};

/// A frontal matrix in the multifrontal method.
///
/// Represents a dense submatrix indexed by a set of global indices.
/// The first `nfs` indices are "fully summed" (to be eliminated at this node).
/// The remaining indices form the "contribution block" passed to the parent.
#[derive(Debug, Clone)]
pub struct FrontalMatrix {
    /// The dense matrix (size = indices.len() x indices.len()).
    pub mat: DenseMat,
    /// Global row/column indices. The first `nfs` are fully summed.
    pub indices: Vec<usize>,
    /// Number of fully summed variables.
    pub nfs: usize,
}

/// Result of partially factoring a frontal matrix.
#[derive(Debug)]
pub struct PartialFactorResult {
    /// Bunch-Kaufman factorization of the fully-summed block.
    pub bk: BunchKaufmanResult,
    /// L21 block: (ncb x nfs) matrix, the sub-diagonal part of L for this front.
    pub l21: DenseMat,
    /// Contribution (Schur complement) block: (ncb x ncb), to be passed to parent.
    pub contrib: DenseMat,
    /// Global indices of the contribution block (subset of front indices, after nfs).
    pub contrib_indices: Vec<usize>,
    /// Global indices of the fully-summed variables (for storing L entries).
    pub fs_indices: Vec<usize>,
}

impl FrontalMatrix {
    /// Create a new zero frontal matrix.
    pub fn new(indices: Vec<usize>, nfs: usize) -> Self {
        let size = indices.len();
        Self {
            mat: DenseMat::zeros(size, size),
            indices,
            nfs,
        }
    }

    /// Total size of the front.
    pub fn size(&self) -> usize {
        self.indices.len()
    }

    /// Number of contribution block variables.
    pub fn ncb(&self) -> usize {
        self.indices.len() - self.nfs
    }

    /// Find the local index for a global index, or None.
    pub fn local_index(&self, global: usize) -> Option<usize> {
        // Indices are sorted, so use binary search
        self.indices.binary_search(&global).ok()
    }

    /// Assemble an original matrix entry (global_row, global_col, val) into this front.
    /// The entry should be in the lower triangle (global_row >= global_col) for the
    /// full symmetric storage used internally.
    pub fn assemble_entry(&mut self, global_row: usize, global_col: usize, val: f64) {
        if let (Some(li), Some(lj)) = (self.local_index(global_row), self.local_index(global_col)) {
            self.mat.add(li, lj, val);
            if li != lj {
                self.mat.add(lj, li, val); // symmetric
            }
        }
    }

    /// Extend-add: merge a child's contribution block into this front.
    pub fn extend_add(&mut self, contrib: &DenseMat, contrib_indices: &[usize]) {
        let ncb = contrib_indices.len();
        // Precompute local indices for the contribution
        let local_map: Vec<usize> = contrib_indices
            .iter()
            .map(|&gi| self.local_index(gi).expect("extend_add: index not found in parent"))
            .collect();

        let size = self.mat.nrows;
        for cj in 0..ncb {
            let lj = local_map[cj];
            let dst_base = lj * size;
            let src_base = cj * ncb;
            for ci in 0..ncb {
                self.mat.data[dst_base + local_map[ci]] += contrib.data[src_base + ci];
            }
        }
    }

    /// Partial factorization: factor the fully-summed block, compute L21 and Schur complement.
    ///
    /// After this, the front contains:
    /// - Top-left nfs x nfs: factored (L11, D11) via Bunch-Kaufman
    /// - Bottom-left ncb x nfs: L21 = A21 * D11^{-1} (through BK)
    /// - Bottom-right ncb x ncb: Schur complement S = A22 - L21 * D11 * L21^T
    pub fn partial_factor(self) -> PartialFactorResult {
        let nfs = self.nfs;
        let ncb = self.ncb();
        let size = self.size();

        let fs_indices = self.indices[..nfs].to_vec();
        let contrib_indices = self.indices[nfs..].to_vec();

        if nfs == size {
            // No contribution block — just factor the whole thing
            let mut full_mat = self.mat;
            let bk = dense_ldlt_bunch_kaufman(&mut full_mat);
            return PartialFactorResult {
                bk,
                l21: DenseMat::zeros(0, nfs),
                contrib: DenseMat::zeros(0, 0),
                contrib_indices,
                fs_indices,
            };
        }

        // Extract the fully-summed block A11 (nfs x nfs) and factor it
        let mut a11 = DenseMat::zeros(nfs, nfs);
        for i in 0..nfs {
            for j in 0..nfs {
                a11.set(i, j, self.mat.get(i, j));
            }
        }
        let bk = dense_ldlt_bunch_kaufman(&mut a11);

        // Compute L21: each row of A21 is solved via the BK factorization
        // L21 = A21 * (L11 * D11 * L11^T)^{-1} ... but we need to be careful.
        // Actually, L21[i, :] satisfies: A21[i, :] = L21[i, :] * D11 * L11^T
        // So L21[i, :] * D11 * L11^T = A21[i, :]
        // Let's compute L21 column by column. For the BK factorization P*L*D*L^T*P^T = A11,
        // the L21 block must satisfy: A21 * P = L21_perm * L11 * D11
        // where L21_perm is L21 with columns permuted by P.
        //
        // More precisely, for the full matrix:
        // [A11  A12] = [P 0] [L11  0 ] [D11  0 ] [L11^T  L21^T] [P^T 0]
        // [A21  A22]   [0 I] [L21  I ] [0    S ] [0      I    ] [0   I]
        //
        // From the (2,1) block: A21 * P = L21 * L11 * D11... no.
        // Let me think again. The block factorization of the permuted matrix is:
        // P^T * A * P (rows/cols of the FS block permuted by BK pivot perm).
        //
        // Simpler approach: solve for L21 row by row.
        // A21[i,:] = L21[i,:] * D * L^T, permuted appropriately.
        // For each row i of A21, solve L*D*L^T * x = P^T * A21[i,:]^T => L21[i,:] = x^T * P^T... complicated.
        //
        // Easiest correct approach: use the BK factorization directly.
        // L21_row solves: for each row i of the CB block,
        //   A21_permuted[i, :] = L21_permuted[i, :] * L11 * D11
        //   => L21_permuted[i, :] = A21_permuted[i, :] * D11^{-1} * L11^{-1}
        //
        // Step 1: permute A21 columns by BK perm
        // Step 2: forward solve with L11
        // Step 3: solve with D11

        // Compute L21 via batched column-oriented forward solve (TRSM-like).
        // Build A21 with permuted columns: a21 is ncb x nfs, column-major.
        let mut l21 = DenseMat::zeros(ncb, nfs);
        for j in 0..nfs {
            let src_col = bk.perm[j]; // permuted column
            for i in 0..ncb {
                l21.set(i, j, self.mat.get(nfs + i, src_col));
            }
        }

        // Column-oriented forward solve: L * Z = A21
        // For each column j of L, update all later columns of l21.
        // L is row-major in bk.l: L(col, j) = bk.l.data[col * nfs + j]
        {
            let l21d = &mut l21.data;
            for j in 0..nfs {
                for col in (j + 1)..nfs {
                    let l_val = bk.l.data[col * nfs + j];
                    if l_val != 0.0 {
                        let src_base = j * ncb;
                        let dst_base = col * ncb;
                        // Contiguous inner loop — auto-vectorizable
                        for i in 0..ncb {
                            l21d[dst_base + i] -= l_val * l21d[src_base + i];
                        }
                    }
                }
            }
        }

        // Column-oriented D solve: handle 1x1 and 2x2 blocks
        {
            let l21d = &mut l21.data;
            let mut k = 0;
            while k < nfs {
                if k + 1 < nfs && bk.d_offdiag[k].abs() > 1e-12 {
                    // 2x2 block
                    let a = bk.d_diag[k];
                    let b = bk.d_offdiag[k];
                    let c = bk.d_diag[k + 1];
                    let det = a * c - b * b;
                    let inv_det = 1.0 / det;
                    let base0 = k * ncb;
                    let base1 = (k + 1) * ncb;
                    for i in 0..ncb {
                        let r0 = l21d[base0 + i];
                        let r1 = l21d[base1 + i];
                        l21d[base0 + i] = (c * r0 - b * r1) * inv_det;
                        l21d[base1 + i] = (a * r1 - b * r0) * inv_det;
                    }
                    k += 2;
                } else {
                    if bk.d_diag[k].abs() > 1e-30 {
                        let inv_d = 1.0 / bk.d_diag[k];
                        let base = k * ncb;
                        for i in 0..ncb {
                            l21d[base + i] *= inv_d;
                        }
                    }
                    k += 1;
                }
            }
        }

        // Compute Schur complement: S = A22 - L21 * D * L21^T
        // where D is the block diagonal from BK.
        let mut contrib = DenseMat::zeros(ncb, ncb);
        for j in 0..ncb {
            let src_col = &self.mat.data[(nfs + j) * size + nfs..(nfs + j) * size + nfs + ncb];
            let dst_col = &mut contrib.data[j * ncb..j * ncb + ncb];
            dst_col.copy_from_slice(src_col);
        }

        // S -= L21 * D * L21^T using cache-blocked GEMM
        // Step 1: Compute W = L21 * D (scale columns by block-diagonal D)
        let l_data = &l21.data;
        let mut w_data = vec![0.0f64; ncb * nfs];
        {
            let mut k = 0;
            while k < nfs {
                if k + 1 < nfs && bk.d_offdiag[k].abs() > 1e-12 {
                    let d00 = bk.d_diag[k];
                    let d01 = bk.d_offdiag[k];
                    let d11 = bk.d_diag[k + 1];
                    let l0 = &l_data[k * ncb..(k + 1) * ncb];
                    let l1 = &l_data[(k + 1) * ncb..(k + 2) * ncb];
                    let w0 = &mut w_data[k * ncb..(k + 1) * ncb];
                    for i in 0..ncb {
                        w0[i] = d00 * l0[i] + d01 * l1[i];
                    }
                    let w1 = &mut w_data[(k + 1) * ncb..(k + 2) * ncb];
                    for i in 0..ncb {
                        w1[i] = d01 * l0[i] + d11 * l1[i];
                    }
                    k += 2;
                } else {
                    let dk = bk.d_diag[k];
                    let l_col = &l_data[k * ncb..(k + 1) * ncb];
                    let w_col = &mut w_data[k * ncb..(k + 1) * ncb];
                    for i in 0..ncb {
                        w_col[i] = dk * l_col[i];
                    }
                    k += 1;
                }
            }
        }

        // Step 2: S -= W * L21^T using cache-blocked GEMM-NT
        // W is ncb × nfs (col-major, stride ncb), L21 is ncb × nfs (col-major, stride ncb)
        gemm_nt_sub(
            ncb, ncb, nfs,
            &w_data, ncb,
            l_data, ncb,
            &mut contrib.data, ncb,
        );

        PartialFactorResult {
            bk,
            l21,
            contrib,
            contrib_indices,
            fs_indices,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frontal_assemble_and_factor_full() {
        // 2x2 front, both fully summed (like a leaf with no CB)
        // A = [[4, 2], [2, 5]]
        let mut front = FrontalMatrix::new(vec![0, 1], 2);
        front.assemble_entry(0, 0, 4.0);
        front.assemble_entry(0, 1, 2.0);
        front.assemble_entry(1, 1, 5.0);

        let result = front.partial_factor();
        assert_eq!(result.bk.inertia.positive, 2);
        assert_eq!(result.bk.inertia.negative, 0);
        assert_eq!(result.contrib_indices.len(), 0);
    }

    #[test]
    fn test_frontal_partial_factor() {
        // 3x3 matrix, front eliminates variable 0 (nfs=1), CB = {1, 2}
        // A = [[4, 2, 1], [2, 5, 3], [1, 3, 6]]
        let mut front = FrontalMatrix::new(vec![0, 1, 2], 1);
        front.assemble_entry(0, 0, 4.0);
        front.assemble_entry(0, 1, 2.0);
        front.assemble_entry(0, 2, 1.0);
        front.assemble_entry(1, 1, 5.0);
        front.assemble_entry(1, 2, 3.0);
        front.assemble_entry(2, 2, 6.0);

        let result = front.partial_factor();
        assert_eq!(result.bk.inertia.positive, 1); // D[0] = 4 > 0
        assert_eq!(result.contrib_indices, vec![1, 2]);

        // Schur complement: S = A22 - A21 * A11^{-1} * A12
        // A11 = 4, A21 = [2; 1], A12 = [2, 1], A22 = [[5,3],[3,6]]
        // S = [[5 - 4/4, 3 - 2/4], [3 - 2/4, 6 - 1/4]] = [[4, 2.5], [2.5, 5.75]]
        assert!((result.contrib.get(0, 0) - 4.0).abs() < 1e-10, "S[0,0] = {}", result.contrib.get(0, 0));
        assert!((result.contrib.get(0, 1) - 2.5).abs() < 1e-10, "S[0,1] = {}", result.contrib.get(0, 1));
        assert!((result.contrib.get(1, 0) - 2.5).abs() < 1e-10);
        assert!((result.contrib.get(1, 1) - 5.75).abs() < 1e-10, "S[1,1] = {}", result.contrib.get(1, 1));
    }

    #[test]
    fn test_extend_add() {
        // Parent front with indices {1, 2, 3}, child contrib has indices {2, 3}
        let mut parent = FrontalMatrix::new(vec![1, 2, 3], 1);
        // Some existing entries
        parent.assemble_entry(1, 1, 10.0);
        parent.assemble_entry(2, 2, 20.0);
        parent.assemble_entry(3, 3, 30.0);

        // Child contribution
        let mut contrib = DenseMat::zeros(2, 2);
        contrib.set(0, 0, 1.0); // index 2,2
        contrib.set(0, 1, 0.5); // index 2,3
        contrib.set(1, 0, 0.5); // index 3,2
        contrib.set(1, 1, 2.0); // index 3,3

        parent.extend_add(&contrib, &[2, 3]);

        assert!((parent.mat.get(1, 1) - 21.0).abs() < 1e-15); // 20 + 1
        assert!((parent.mat.get(1, 2) - 0.5).abs() < 1e-15);
        assert!((parent.mat.get(2, 1) - 0.5).abs() < 1e-15);
        assert!((parent.mat.get(2, 2) - 32.0).abs() < 1e-15); // 30 + 2
    }
}
