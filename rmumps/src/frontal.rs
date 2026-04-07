use crate::dense::{DenseMat, gemm_nt_sub};
use crate::pivot::{
    dense_ldlt_bunch_kaufman, BunchKaufmanResult, PivotResult,
    compute_inertia,
};

// Thread-local reusable buffer for W = L21 * D in Schur complement computation.
thread_local! {
    static W_BUF: std::cell::RefCell<Vec<f64>> = std::cell::RefCell::new(Vec::new());
}

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
    /// Number of FS columns that were actually eliminated (may be < original nfs due to delays).
    pub nfs_eliminated: usize,
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

    /// Promote delayed child FS columns from the CB portion to the FS portion.
    /// This moves the specified columns (by global index) from the CB set to the FS set,
    /// rearranging the front matrix accordingly. The promoted columns are appended
    /// after the current FS columns.
    pub fn promote_cb_to_fs(&mut self, delayed_cols: &[usize]) {
        if delayed_cols.is_empty() {
            return;
        }

        let n = self.size();
        let old_nfs = self.nfs;

        // Find local indices (in CB portion) of the delayed columns
        let mut cb_local_positions: Vec<usize> = Vec::new();
        for &gc in delayed_cols {
            if let Ok(pos) = self.indices[old_nfs..].binary_search(&gc) {
                cb_local_positions.push(old_nfs + pos);
            }
        }

        if cb_local_positions.is_empty() {
            return;
        }

        // Build new ordering: [old FS cols] [promoted cols] [remaining CB cols]
        let mut new_order: Vec<usize> = (0..old_nfs).collect();
        let cb_promote_set: std::collections::HashSet<usize> = cb_local_positions.iter().copied().collect();
        for &pos in &cb_local_positions {
            new_order.push(pos);
        }
        for i in old_nfs..n {
            if !cb_promote_set.contains(&i) {
                new_order.push(i);
            }
        }

        // Apply the permutation to indices and matrix
        let new_indices: Vec<usize> = new_order.iter().map(|&i| self.indices[i]).collect();

        // Permute the dense matrix: new[i,j] = old[new_order[i], new_order[j]]
        let mut new_data = vec![0.0; n * n];
        for new_j in 0..n {
            let old_j = new_order[new_j];
            for new_i in 0..n {
                let old_i = new_order[new_i];
                new_data[new_j * n + new_i] = self.mat.data[old_j * n + old_i];
            }
        }

        self.indices = new_indices;
        self.mat.data = new_data;
        self.nfs = old_nfs + cb_local_positions.len();
    }

    /// Find the local index for a global index, or None.
    pub fn local_index(&self, global: usize) -> Option<usize> {
        // Front indices have structure [FS cols (sorted) | CB cols (sorted)].
        // FS cols may have higher values than CB cols (due to delayed pivoting
        // expansion), so the full list is NOT globally sorted.
        // Search FS and CB portions separately.
        if let Ok(pos) = self.indices[..self.nfs].binary_search(&global) {
            return Some(pos);
        }
        if let Ok(pos) = self.indices[self.nfs..].binary_search(&global) {
            return Some(self.nfs + pos);
        }
        None
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
    /// If the contribution contains delayed child FS columns not in this front,
    /// the front is dynamically expanded to accommodate them.
    pub fn extend_add(&mut self, contrib: &DenseMat, contrib_indices: &[usize]) {
        let ncb = contrib_indices.len();

        // Check if any contrib indices are missing from the front
        let mut missing: Vec<usize> = Vec::new();
        for &gi in contrib_indices {
            if self.local_index(gi).is_none() {
                missing.push(gi);
            }
        }

        // Dynamically expand the front if needed (for delayed pivoting)
        if !missing.is_empty() {
            self.expand_for_delayed(&missing);
        }

        // Now all indices should be present
        let local_map: Vec<usize> = contrib_indices
            .iter()
            .map(|&gi| self.local_index(gi).expect("extend_add: index not found after expansion"))
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

    /// Expand the front to include additional indices in the CB portion.
    /// This is used for delayed pivoting when a child's delayed columns
    /// weren't predicted by the symbolic phase.
    fn expand_for_delayed(&mut self, new_indices: &[usize]) {
        let old_size = self.size();
        let new_size = old_size + new_indices.len();
        let nfs = self.nfs;

        // Build the expanded index array: [FS cols] [old CB + new indices, sorted]
        let mut new_front_indices = self.indices.clone();
        new_front_indices.extend_from_slice(new_indices);

        // Build sort permutation for CB portion so we can permute matrix data to match.
        // perm[new_pos - nfs] = old_pos, where old_pos < old_size means an existing
        // CB column and old_pos >= old_size means a newly added column.
        let cb_len = new_size - nfs;
        let mut cb_perm: Vec<usize> = (0..cb_len).collect();
        cb_perm.sort_unstable_by_key(|&i| new_front_indices[nfs + i]);
        // Apply the sort to the indices
        let sorted_cb: Vec<usize> = cb_perm.iter().map(|&i| new_front_indices[nfs + i]).collect();
        new_front_indices[nfs..].copy_from_slice(&sorted_cb);

        // Build full old-to-new position mapping:
        // FS columns (0..nfs) keep their positions.
        // CB columns are reordered according to cb_perm.
        let mut old_to_new = vec![usize::MAX; new_size];
        for i in 0..nfs {
            old_to_new[i] = i;
        }
        for (new_cb_pos, &old_cb_pos) in cb_perm.iter().enumerate() {
            old_to_new[nfs + old_cb_pos] = nfs + new_cb_pos;
        }

        // Expand and permute the dense matrix data.
        // Old entries at (old_i, old_j) go to (new_i, new_j).
        // New columns (old_pos >= old_size) are zero-initialized.
        let mut new_data = vec![0.0; new_size * new_size];
        for old_j in 0..old_size {
            let new_j = old_to_new[old_j];
            if new_j == usize::MAX { continue; }
            for old_i in 0..old_size {
                let new_i = old_to_new[old_i];
                if new_i == usize::MAX { continue; }
                new_data[new_j * new_size + new_i] = self.mat.data[old_j * old_size + old_i];
            }
        }

        self.indices = new_front_indices;
        self.mat = crate::dense::DenseMat {
            nrows: new_size,
            ncols: new_size,
            data: new_data,
        };
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
                nfs_eliminated: nfs,
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
        // Reuse thread-local buffer for W to avoid per-call allocation
        let l_data = &l21.data;
        W_BUF.with(|buf| {
            let mut w_buf = buf.borrow_mut();
            let w_len = ncb * nfs;
            if w_buf.len() < w_len {
                w_buf.resize(w_len, 0.0);
            } else {
                w_buf[..w_len].fill(0.0);
            }
            let w_data = &mut w_buf[..w_len];
            {
                let mut k = 0;
                while k < nfs {
                    if k + 1 < nfs && bk.d_offdiag[k].abs() > 1e-12 {
                        let d00 = bk.d_diag[k];
                        let d01 = bk.d_offdiag[k];
                        let d11 = bk.d_diag[k + 1];
                        let l0 = &l_data[k * ncb..(k + 1) * ncb];
                        let l1 = &l_data[(k + 1) * ncb..(k + 2) * ncb];
                        let (w0, w_rest) = w_data[k * ncb..].split_at_mut(ncb);
                        let w1 = &mut w_rest[..ncb];
                        for i in 0..ncb {
                            w0[i] = d00 * l0[i] + d01 * l1[i];
                        }
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
            gemm_nt_sub(
                ncb, ncb, nfs,
                w_data, ncb,
                l_data, ncb,
                &mut contrib.data, ncb,
            );
        });

        let nfs_eliminated = nfs;
        PartialFactorResult {
            bk,
            l21,
            contrib,
            contrib_indices,
            fs_indices,
            nfs_eliminated,
        }
    }

    /// Partial factorization with threshold pivoting and delayed pivots.
    ///
    /// Unlike `partial_factor`, this can reject pivots that fail the threshold test.
    /// Rejected FS columns are moved to the contribution block (delayed to parent).
    /// This is the key mechanism that makes MA57/MUMPS reliable on KKT systems.
    pub fn partial_factor_threshold(self, threshold: f64, _n_primal: Option<usize>) -> PartialFactorResult {
        let mut orig_nfs = self.nfs;
        let orig_size = self.size();

        if orig_nfs == 0 {
            // No FS columns — just pass through as contribution
            return PartialFactorResult {
                bk: BunchKaufmanResult {
                    l: DenseMat::zeros(0, 0),
                    d_diag: vec![],
                    d_offdiag: vec![],
                    perm: vec![],
                    perm_inv: vec![],
                    inertia: crate::Inertia { positive: 0, negative: 0, zero: 0 },
                },
                l21: DenseMat::zeros(orig_size, 0),
                contrib: self.mat,
                contrib_indices: self.indices,
                fs_indices: vec![],
                nfs_eliminated: 0,
            };
        }

        // Work with the full dense matrix. We'll perform threshold pivoting
        // on the FS block, rejecting columns that fail the threshold test.
        let mut a = self.mat.data.clone();
        let n = orig_size;

        // Track which FS columns are eliminated vs delayed
        // perm[k] = original column index in the front
        let mut perm: Vec<usize> = (0..n).collect();
        let mut nfs_elim = 0usize; // number of successfully eliminated columns
        let mut d_diag = vec![0.0; orig_nfs];
        let mut d_offdiag = vec![0.0; orig_nfs];

        // L is stored column-major: l[col * n + row]
        let mut l_data = vec![0.0; n * orig_nfs];
        let mut work = vec![0.0; 2 * n];

        let mut k = 0;
        while k < orig_nfs {
            // Only look for pivots among the remaining FS columns [k..orig_nfs]
            // But we also need to consider the active submatrix starting at index nfs_elim
            // in the permuted system.
            //
            // At this point, columns 0..nfs_elim have been eliminated.
            // Columns nfs_elim..orig_nfs are remaining FS candidates.
            // Columns orig_nfs..n are CB.

            let active_start = nfs_elim;

            // Find pivot in the active FS portion
            let fs_remaining = orig_nfs - k;

            // Find best pivot among remaining FS columns
            let pivot = find_pivot_in_fs(
                &a, n, active_start, fs_remaining, threshold,
            );

            match pivot {
                PivotResult::OneByOne(p) => {
                    // Swap p to position active_start
                    if p != active_start {
                        swap_full(&mut a, n, active_start, p);
                        perm.swap(active_start, p);
                        // Swap L entries for previously eliminated columns
                        for j in 0..nfs_elim {
                            l_data.swap(j * n + active_start, j * n + p);
                        }
                    }

                    let akk = a[active_start * n + active_start];
                    d_diag[nfs_elim] = akk;

                    if akk.abs() > 1e-30 {
                        let m = n - active_start - 1;
                        for i in 0..m {
                            work[i] = a[(active_start + 1 + i) * n + active_start] / akk;
                            l_data[nfs_elim * n + (active_start + 1 + i)] = work[i];
                        }
                        // Update trailing matrix
                        for i in 0..m {
                            let si = work[i] * akk;
                            let base = (active_start + 1 + i) * n + (active_start + 1);
                            for j in 0..m {
                                a[base + j] -= si * work[j];
                            }
                        }
                    }
                    l_data[nfs_elim * n + active_start] = 1.0;
                    nfs_elim += 1;
                    k += 1;
                }
                PivotResult::TwoByTwo(p1, p2) => {
                    // Need to bring p2 to active_start+1, p1 to active_start
                    if p2 != active_start + 1 {
                        swap_full(&mut a, n, active_start + 1, p2);
                        perm.swap(active_start + 1, p2);
                        for j in 0..nfs_elim {
                            l_data.swap(j * n + (active_start + 1), j * n + p2);
                        }
                    }
                    if p1 != active_start {
                        swap_full(&mut a, n, active_start, p1);
                        perm.swap(active_start, p1);
                        for j in 0..nfs_elim {
                            l_data.swap(j * n + active_start, j * n + p1);
                        }
                    }

                    let akk = a[active_start * n + active_start];
                    let ak1k = a[(active_start + 1) * n + active_start];
                    let ak1k1 = a[(active_start + 1) * n + (active_start + 1)];

                    d_diag[nfs_elim] = akk;
                    d_diag[nfs_elim + 1] = ak1k1;
                    d_offdiag[nfs_elim] = ak1k;

                    let det = akk * ak1k1 - ak1k * ak1k;

                    if det.abs() > 1e-30 {
                        let d_inv_00 = ak1k1 / det;
                        let d_inv_01 = -ak1k / det;
                        let d_inv_11 = akk / det;

                        let m = n - active_start - 2;
                        for i in 0..m {
                            let aik = a[(active_start + 2 + i) * n + active_start];
                            let aik1 = a[(active_start + 2 + i) * n + (active_start + 1)];
                            work[i] = aik * d_inv_00 + aik1 * d_inv_01;
                            work[m + i] = aik * d_inv_01 + aik1 * d_inv_11;
                            l_data[nfs_elim * n + (active_start + 2 + i)] = work[i];
                            l_data[(nfs_elim + 1) * n + (active_start + 2 + i)] = work[m + i];
                        }

                        // Update trailing matrix
                        for i in 0..m {
                            let li0 = work[i];
                            let li1 = work[m + i];
                            let si0 = li0 * akk + li1 * ak1k;
                            let si1 = li0 * ak1k + li1 * ak1k1;
                            let base = (active_start + 2 + i) * n + (active_start + 2);
                            for j in 0..m {
                                a[base + j] -= si0 * work[j] + si1 * work[m + j];
                            }
                        }
                    }

                    l_data[nfs_elim * n + active_start] = 1.0;
                    l_data[(nfs_elim + 1) * n + (active_start + 1)] = 1.0;
                    nfs_elim += 2;
                    k += 2;
                }
                PivotResult::Delayed => {
                    // No FS-only pivot passed threshold. Before genuinely delaying,
                    // search CB columns for a 2×2 partner (MA57-style).
                    // This pairs primal-dual variables across the FS-CB boundary.
                    let cb_partner = find_cb_pivot_partner(
                        &a, n, active_start, orig_nfs - k, orig_nfs, threshold,
                    );

                    if let Some((fs_pos, cb_pos)) = cb_partner {
                        // Found a good FS-CB 2×2 pivot. Promote the CB column
                        // into the FS range by swapping it to orig_nfs position
                        // (right after the last FS column).
                        //
                        // First, swap the CB column to position orig_nfs (end of FS range)
                        if cb_pos != orig_nfs {
                            swap_full(&mut a, n, cb_pos, orig_nfs);
                            perm.swap(cb_pos, orig_nfs);
                            for j in 0..nfs_elim {
                                l_data.swap(j * n + cb_pos, j * n + orig_nfs);
                            }
                        }
                        // Increase FS range to include the promoted column
                        orig_nfs += 1;
                        // Expand d_diag/d_offdiag/l_data to accommodate all FS columns
                        if orig_nfs > d_diag.len() {
                            d_diag.resize(orig_nfs, 0.0);
                            d_offdiag.resize(orig_nfs, 0.0);
                        }
                        if orig_nfs > l_data.len() / n {
                            l_data.resize(orig_nfs * n, 0.0);
                        }

                        // Now swap fs_pos to active_start, and orig_nfs-1 to active_start+1
                        // (the promoted column is now at orig_nfs - 1)
                        let promoted_pos = orig_nfs - 1;
                        if promoted_pos != active_start + 1 {
                            swap_full(&mut a, n, active_start + 1, promoted_pos);
                            perm.swap(active_start + 1, promoted_pos);
                            for j in 0..nfs_elim {
                                l_data.swap(j * n + (active_start + 1), j * n + promoted_pos);
                            }
                        }
                        if fs_pos != active_start {
                            swap_full(&mut a, n, active_start, fs_pos);
                            perm.swap(active_start, fs_pos);
                            for j in 0..nfs_elim {
                                l_data.swap(j * n + active_start, j * n + fs_pos);
                            }
                        }

                        // Now eliminate as 2×2 pivot (same code as TwoByTwo handler)
                        let akk = a[active_start * n + active_start];
                        let ak1k = a[(active_start + 1) * n + active_start];
                        let ak1k1 = a[(active_start + 1) * n + (active_start + 1)];

                        d_diag[nfs_elim] = akk;
                        d_diag[nfs_elim + 1] = ak1k1;
                        d_offdiag[nfs_elim] = ak1k;

                        let det = akk * ak1k1 - ak1k * ak1k;
                        if det.abs() > 1e-30 {
                            let d_inv_00 = ak1k1 / det;
                            let d_inv_01 = -ak1k / det;
                            let d_inv_11 = akk / det;
                            let m = n - active_start - 2;
                            for i in 0..m {
                                let aik = a[(active_start + 2 + i) * n + active_start];
                                let aik1 = a[(active_start + 2 + i) * n + (active_start + 1)];
                                work[i] = aik * d_inv_00 + aik1 * d_inv_01;
                                work[m + i] = aik * d_inv_01 + aik1 * d_inv_11;
                                l_data[nfs_elim * n + (active_start + 2 + i)] = work[i];
                                l_data[(nfs_elim + 1) * n + (active_start + 2 + i)] = work[m + i];
                            }
                            for i in 0..m {
                                let li0 = work[i];
                                let li1 = work[m + i];
                                let si0 = li0 * akk + li1 * ak1k;
                                let si1 = li0 * ak1k + li1 * ak1k1;
                                let base = (active_start + 2 + i) * n + (active_start + 2);
                                for j in 0..m {
                                    a[base + j] -= si0 * work[j] + si1 * work[m + j];
                                }
                            }
                        }
                        l_data[nfs_elim * n + active_start] = 1.0;
                        l_data[(nfs_elim + 1) * n + (active_start + 1)] = 1.0;
                        nfs_elim += 2;
                        k += 2; // consumed 1 original FS + 1 promoted CB
                    } else {
                        // No CB partner found — genuinely delay this column.
                        let last_fs = active_start + (orig_nfs - k) - 1;
                        if active_start != last_fs {
                            swap_full(&mut a, n, active_start, last_fs);
                            perm.swap(active_start, last_fs);
                            for j in 0..nfs_elim {
                                l_data.swap(j * n + active_start, j * n + last_fs);
                            }
                        }
                        k += 1;
                    }
                }
            }
        }

        // Now nfs_elim columns have been eliminated. The remaining columns
        // [nfs_elim..n] form the contribution block (including delayed FS columns).
        let ncb_new = n - nfs_elim;


        // Build the BK result from what we've computed
        let d_diag = d_diag[..nfs_elim].to_vec();
        let d_offdiag = d_offdiag[..nfs_elim].to_vec();

        let inertia = compute_inertia(&d_diag, &d_offdiag, nfs_elim);

        // Build perm/perm_inv for the eliminated block.
        // fs_indices is already built in factored order (fs_indices[k] = global
        // index of the k-th pivot), so the BK perm is the identity mapping.
        // The solve code uses bk.perm to index into fs_indices.
        let bk_perm: Vec<usize> = (0..nfs_elim).collect();
        let mut bk_perm_inv = vec![0usize; nfs_elim];
        for i in 0..nfs_elim {
            bk_perm_inv[i] = i;
        }

        // Build L factor for the eliminated block (nfs_elim x nfs_elim)
        // l_data is stored column-major: l_data[col * n + row] = L[row, col]
        // but the solve code expects row-major: l.data[row * nfs + col] = L[row, col]
        let mut l_factor = DenseMat::zeros(nfs_elim, nfs_elim);
        for col in 0..nfs_elim {
            for row in 0..nfs_elim {
                l_factor.data[row * nfs_elim + col] = l_data[col * n + row];
            }
        }

        // Build L21: rows [nfs_elim..n], cols [0..nfs_elim]
        let mut l21 = DenseMat::zeros(ncb_new, nfs_elim);
        for col in 0..nfs_elim {
            for row in 0..ncb_new {
                l21.data[col * ncb_new + row] = l_data[col * n + (nfs_elim + row)];
            }
        }

        // Extract contribution block from the trailing matrix
        let mut contrib = DenseMat::zeros(ncb_new, ncb_new);
        for j in 0..ncb_new {
            for i in 0..ncb_new {
                contrib.data[j * ncb_new + i] = a[(nfs_elim + j) * n + (nfs_elim + i)];
            }
        }

        // Map global indices
        let fs_indices: Vec<usize> = (0..nfs_elim).map(|i| self.indices[perm[i]]).collect();
        let contrib_indices: Vec<usize> = (nfs_elim..n).map(|i| self.indices[perm[i]]).collect();

        let bk = BunchKaufmanResult {
            l: l_factor,
            d_diag,
            d_offdiag,
            perm: bk_perm,
            perm_inv: bk_perm_inv,
            inertia,
        };

        PartialFactorResult {
            bk,
            l21,
            contrib,
            contrib_indices,
            fs_indices,
            nfs_eliminated: nfs_elim,
        }
    }
}

/// Find pivot among the first `fs_remaining` columns of the active submatrix
/// starting at `start` in an n x n matrix.
fn find_pivot_in_fs(
    a: &[f64],
    n: usize,
    start: usize,
    fs_remaining: usize,
    threshold: f64,
) -> PivotResult {
    if fs_remaining == 0 {
        return PivotResult::Delayed;
    }

    // Search for best 1x1 pivot among FS columns
    let fs_end = start + fs_remaining;

    // Try each FS column as a potential 1x1 pivot
    let mut best_1x1: Option<(usize, f64)> = None; // (col, ratio)

    for col in start..fs_end {
        let diag = a[col * n + col].abs();
        // Find max off-diagonal in column (full active submatrix)
        let mut max_offdiag = 0.0f64;
        for row in start..n {
            if row != col {
                max_offdiag = max_offdiag.max(a[row * n + col].abs());
            }
        }

        if max_offdiag == 0.0 && diag == 0.0 {
            continue; // skip zero columns
        }

        if max_offdiag == 0.0 {
            // Pure diagonal — always acceptable
            let ratio = f64::INFINITY;
            if best_1x1.map_or(true, |(_, r)| ratio > r) {
                best_1x1 = Some((col, ratio));
            }
            continue;
        }

        let ratio = diag / max_offdiag;
        if ratio >= threshold {
            if best_1x1.map_or(true, |(_, r)| ratio > r) {
                best_1x1 = Some((col, ratio));
            }
        }
    }

    if let Some((col, _)) = best_1x1 {
        return PivotResult::OneByOne(col);
    }

    // No 1x1 pivot passed threshold — try 2x2 pivots
    for i in start..fs_end {
        for j in (i + 1)..fs_end {
            let akk = a[i * n + i];
            let akj = a[j * n + i];
            let ajj = a[j * n + j];
            let det = (akk * ajj - akj * akj).abs();
            let max_elem = akk.abs().max(ajj.abs()).max(akj.abs());
            if max_elem > 1e-30 && det >= threshold * max_elem * max_elem {
                return PivotResult::TwoByTwo(i, j);
            }
        }
    }

    PivotResult::Delayed
}

/// Search CB columns for a 2×2 pivot partner for an FS column.
/// For KKT systems, pairs primal variables with dual variables across the FS-CB boundary.
/// Returns Some((fs_pos, cb_pos)) if a good pair is found.
fn find_cb_pivot_partner(
    a: &[f64],
    n: usize,
    active_start: usize,
    fs_remaining: usize,
    orig_nfs: usize,
    threshold: f64,
) -> Option<(usize, usize)> {
    let fs_end = active_start + fs_remaining;
    let mut best: Option<(usize, usize, f64)> = None; // (fs_pos, cb_pos, |det|)

    for fs_pos in active_start..fs_end {
        let diag_fs = a[fs_pos * n + fs_pos].abs();
        let is_zero_diag_fs = diag_fs < 1e-12;

        // Search CB columns for a partner with complementary diagonal
        // (pair zero-diagonal with non-zero-diagonal for good 2×2 pivots)
        for cb_pos in orig_nfs..n {
            let diag_cb = a[cb_pos * n + cb_pos].abs();
            let is_zero_diag_cb = diag_cb < 1e-12;

            // Only pair zero-diagonal with non-zero-diagonal
            if is_zero_diag_fs == is_zero_diag_cb {
                continue;
            }

            let akk = a[fs_pos * n + fs_pos];
            let akj = a[cb_pos * n + fs_pos]; // off-diagonal coupling
            let ajj = a[cb_pos * n + cb_pos];

            // Skip if no coupling
            if akj.abs() < 1e-30 {
                continue;
            }

            let det = (akk * ajj - akj * akj).abs();
            let max_elem = akk.abs().max(ajj.abs()).max(akj.abs());

            if max_elem > 1e-30 && det >= threshold * max_elem * max_elem {
                // Valid 2×2 pivot — track the best (largest determinant)
                if best.map_or(true, |(_, _, d)| det > d) {
                    best = Some((fs_pos, cb_pos, det));
                }
            }
        }
    }

    best.map(|(fs, cb, _)| (fs, cb))
}

/// Swap rows and columns p and q in the full n x n matrix stored in `a`.
fn swap_full(a: &mut [f64], n: usize, p: usize, q: usize) {
    if p == q {
        return;
    }
    for j in 0..n {
        a.swap(p * n + j, q * n + j);
    }
    for i in 0..n {
        a.swap(i * n + p, i * n + q);
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
