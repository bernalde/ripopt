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
        self.partial_factor_threshold_inner(threshold, _n_primal, 0)
    }

    /// Partial factorization with must-eliminate support for promoted columns.
    /// `n_must_eliminate` columns at the END of the FS range are force-eliminated
    /// via static pivoting if they fail the threshold test, preventing multi-hop
    /// delays that lose pivots in the elimination tree.
    pub fn partial_factor_threshold_with_must_eliminate(
        self, threshold: f64, _n_primal: Option<usize>, n_must_eliminate: usize,
    ) -> PartialFactorResult {
        self.partial_factor_threshold_inner(threshold, _n_primal, n_must_eliminate)
    }

    /// MUMPS-style LDLT frontal factorization with threshold pivoting and panel blocking.
    ///
    /// This implements the complete algorithm from MUMPS 5.8.2's dfac_front_LDLT_type1.F:
    /// - Two-level panel blocking (inner 16-32, outer 128)
    /// - First-acceptable pivot search with AMAX/RMAX separation
    /// - Triangular within-block + rectangular beyond-block MQ updates
    /// - U (originals) saved in pivot rows, L (D^{-1}-scaled) in columns
    /// - Between-block GEMM for deferred updates
    /// - MUMPS-style symmetric swap
    ///
    /// Reference: rmumps/MUMPS_LDLT_ALGORITHM.md
    fn partial_factor_threshold_inner(
        self, threshold: f64, _n_primal: Option<usize>, n_must_eliminate: usize,
    ) -> PartialFactorResult {
        let nass = self.nfs; // number of fully-summed columns (may include promoted)
        let nfront = self.size();

        if nass == 0 {
            return PartialFactorResult {
                bk: BunchKaufmanResult {
                    l: DenseMat::zeros(0, 0),
                    d_diag: vec![],
                    d_offdiag: vec![],
                    perm: vec![],
                    perm_inv: vec![],
                    inertia: crate::Inertia { positive: 0, negative: 0, zero: 0 },
                },
                l21: DenseMat::zeros(nfront, 0),
                contrib: self.mat,
                contrib_indices: self.indices,
                fs_indices: vec![],
                nfs_eliminated: 0,
            };
        }

        // Column-major working array: a[col * n + row] = entry(row, col)
        // DenseMat is already column-major: data[col * nrows + row]
        let mut a = self.mat.data.clone();
        let n = nfront;
        let uu = threshold;

        // Track permutation and D factor
        let mut perm: Vec<usize> = (0..n).collect();
        let mut npiv: usize = 0;
        let mut d_diag = vec![0.0; nass];
        let mut d_offdiag = vec![0.0; nass];
        let mut nneg: usize = 0;

        // Must-eliminate range for promoted columns (static pivoting at root)
        let must_eliminate_start = nass.saturating_sub(n_must_eliminate);

        // MUMPS panel blocking parameters
        let nbkjib = if nass > 96 { 32 } else if nass > 32 { 16 } else { nass };
        let nblr: usize = 128;

        let mut iend_blr: usize = 0;
        let mut iend_block: usize = 0;
        let mut last_panel = false;
        let mut ibeg_blr: usize = 0;

        // Outer panel loop
        while iend_blr < nass && !last_panel {
            ibeg_blr = npiv;
            iend_blr = (iend_blr + nblr).min(nass);

            // Inner block loop
            while iend_block < iend_blr && !last_panel {
                let ibeg_block = npiv;
                iend_block = (iend_block + nbkjib).min(iend_blr);

                // Pivot-by-pivot loop within this inner block
                #[cfg(debug_assertions)]
                if nass <= 10 {
                    eprintln!("=== Inner block: ibeg={} iend={} npiv={} nass={} nfront={}", ibeg_block, iend_block, npiv, nass, n);
                }
                loop {
                    if npiv >= iend_block { break; }

                    // STEP 1: Find pivot (MUMPS FAC_I_LDLT)
                    let (inopv, pivsiz) = find_pivot_ldlt(
                        &mut a, n, nass, npiv, iend_block, uu,
                        &mut perm, &mut nneg,
                        must_eliminate_start, nass,
                    );

                    match inopv {
                        1 => {
                            // No pivot in entire NASS — done
                            last_panel = true;
                            break;
                        }
                        2 => {
                            // No pivot in this inner block — advance to next block
                            break;
                        }
                        _ => {} // pivot found (inopv=0 normal, -1 static)
                    }

                    // STEP 2: Within-panel update (MUMPS FAC_MQ_LDLT)
                    let last_row = n; // PIVOT_OPTION=3: update all rows
                    let ifinb = update_within_block_ldlt(
                        &mut a, n, nass, npiv, pivsiz,
                        iend_block, last_row,
                        &mut d_diag, &mut d_offdiag,
                    );

                    npiv += pivsiz;

                    match ifinb {
                        0 => continue,  // more pivots in this block
                        1 => break,     // block done, panel continues
                        _ => { last_panel = true; break; } // -1: NASS done
                    }
                }

                // STEP 3: Between-block GEMM (within BLR panel)
                if iend_blr > iend_block && npiv > ibeg_block {
                    between_block_gemm_ldlt(
                        &mut a, n, ibeg_block, iend_block, npiv,
                        iend_blr, n,
                    );
                }
            }

            // STEP 4: Inter-panel GEMM (between BLR panels)
            // For PIVOT_OPTION=3, MQ already handled L scaling for all rows,
            // so only GEMM is needed (no TRSM). Update columns iend_blr..nass
            // and rows nass..nfront.
            if npiv > ibeg_blr {
                let nel1 = nass.saturating_sub(iend_blr);
                if nel1 > 0 {
                    between_block_gemm_ldlt(
                        &mut a, n, ibeg_blr, iend_blr, npiv,
                        nass, n,
                    );
                }
            }
        }

        // STEP 5: Tail update — apply all pivots to CB columns (MUMPS FAC_T_LDLT).
        // The MQ update within panels only covers columns up to NASS. The CB columns
        // (NASS..NFRONT) need the accumulated update from all pivots.
        if npiv > 0 && n > nass {
            between_block_gemm_ldlt(
                &mut a, n, 0, nass, npiv,
                n, n, // update columns nass..n and rows nass..n
            );
        }

        // Extract results in the format expected by the solve phase.
        let nfs_elim = npiv;
        let ncb_new = n - nfs_elim;

        // Build D factor
        let d_diag = d_diag[..nfs_elim].to_vec();
        let d_offdiag = d_offdiag[..nfs_elim].to_vec();
        let inertia = compute_inertia(&d_diag, &d_offdiag, nfs_elim);

        // Build L11 factor (nfs_elim x nfs_elim) in row-major for solve.
        // In the factored matrix, L is stored in columns below the diagonal
        // (D^{-1}-scaled for within-panel pivots). The L11 block is unit lower
        // triangular: L[i,j] for i > j is at a[j*n + i] (column j, row i).
        // For the solve phase, extract into row-major: l_factor.data[row*nfs+col].
        let bk_perm: Vec<usize> = (0..nfs_elim).collect();
        let bk_perm_inv: Vec<usize> = (0..nfs_elim).collect();
        let mut l_factor = DenseMat::zeros(nfs_elim, nfs_elim);
        let mut col = 0;
        while col < nfs_elim {
            let is_2x2 = col + 1 < nfs_elim && d_offdiag[col].abs() > 1e-30;
            l_factor.data[col * nfs_elim + col] = 1.0; // unit diagonal
            if is_2x2 {
                l_factor.data[(col + 1) * nfs_elim + (col + 1)] = 1.0;
                // L entries are at K1POS = a[row * n + col] (upper triangle), NOT
                // at a[col * n + row] (lower triangle, which has saved originals).
                for row in (col + 2)..nfs_elim {
                    l_factor.data[row * nfs_elim + col] = a[row * n + col];
                    l_factor.data[row * nfs_elim + (col + 1)] = a[row * n + (col + 1)];
                }
                col += 2;
            } else {
                for row in (col + 1)..nfs_elim {
                    l_factor.data[row * nfs_elim + col] = a[row * n + col];
                }
                col += 1;
            }
        }

        // Build L21 (ncb_new x nfs_elim) in column-major for solve.
        // L21[cb_row, pivot_col] is stored at K1POS = a[(nfs_elim+cb_row)*n + pivot_col]
        // (the CB column's entry at the pivot row), NOT at the symmetric position
        // a[pivot_col*n + (nfs_elim+cb_row)] which has the saved original (U storage).
        let mut l21 = DenseMat::zeros(ncb_new, nfs_elim);
        for col in 0..nfs_elim {
            for row in 0..ncb_new {
                l21.data[col * ncb_new + row] = a[(nfs_elim + row) * n + col];
            }
        }

        // Extract contribution block (already in-place in the trailing matrix).
        // Contribution at (i, j) is at a[(nfs_elim+j)*n + (nfs_elim+i)] (column-major).
        //
        // IMPORTANT: When threshold pivoting delays columns (nfs_elim < nass), the
        // MQ within-block update only maintains the upper triangle (col > row in the
        // factored matrix) for within-block columns. The lower triangle entries
        // involving delayed columns are stale. We extract from the upper triangle
        // only and mirror to the lower triangle to ensure symmetry.
        let mut contrib = DenseMat::zeros(ncb_new, ncb_new);
        for col in 0..ncb_new {
            for row in 0..=col {
                // Upper triangle (row <= col): a[(nfs_elim+col)*n + (nfs_elim+row)]
                // has factored_row = nfs_elim+row <= nfs_elim+col = factored_col
                let val = a[(nfs_elim + col) * n + (nfs_elim + row)];
                contrib.data[col * ncb_new + row] = val;
                contrib.data[row * ncb_new + col] = val; // mirror to lower triangle
            }
        }

        // Map global indices
        let fs_indices: Vec<usize> = (0..nfs_elim).map(|i| self.indices[perm[i]]).collect();
        let contrib_indices: Vec<usize> = (nfs_elim..n).map(|i| self.indices[perm[i]]).collect();

        PartialFactorResult {
            bk: BunchKaufmanResult {
                l: l_factor,
                d_diag,
                d_offdiag,
                perm: bk_perm,
                perm_inv: bk_perm_inv,
                inertia,
            },
            l21,
            contrib,
            contrib_indices,
            fs_indices,
            nfs_eliminated: nfs_elim,
        }
    }
}

// ============================================================================
// MUMPS-style helper functions for LDLT frontal factorization
// ============================================================================

/// Symmetric swap of rows/columns p and q in a column-major n×n matrix.
/// Matches MUMPS's DMUMPS_SWAP_LDLT: swaps the full symmetric structure
/// including the "bridge" region between p and q.
fn symmetric_swap_ldlt(a: &mut [f64], n: usize, p: usize, q: usize) {
    if p == q { return; }
    let (p, q) = if p < q { (p, q) } else { (q, p) };

    // 1. Swap UPPER triangle entries in already-factored columns: L entries at
    // rows 0..p-1 in columns p and q. In column-major a[col*n+row], the upper
    // triangle L entries are at a[p*n+row] and a[q*n+row] for row < p.
    // (MUMPS SWAP_LDLT line 2128: dswap(NPIVP1-1, A(col_p), 1, A(col_q), 1))
    for row in 0..p {
        a.swap(p * n + row, q * n + row);
    }

    // 2. Bridge: row p in cols p+1..q-1 ↔ col q in rows p+1..q-1
    for k in 1..(q - p) {
        let idx1 = (p + k) * n + p; // column p+k, row p
        let idx2 = q * n + (p + k); // column q, row p+k
        a.swap(idx1, idx2);
    }

    // 3. Swap diagonals
    a.swap(p * n + p, q * n + q);

    // 4. Swap columns q+1..n-1 at rows p and q
    for col in (q + 1)..n {
        a.swap(col * n + p, col * n + q);
    }
}

/// MUMPS-style pivot search (FAC_I_LDLT).
/// Searches for the FIRST ACCEPTABLE pivot within [npiv, iend_block).
/// Returns (inopv, pivsiz):
///   inopv: 0=found, 1=none in NASS, 2=none in this inner block
///   pivsiz: 1 for 1x1 pivot, 2 for 2x2 pivot
/// On success, the pivot is swapped into position npiv (and npiv+1 for 2x2).
fn find_pivot_ldlt(
    a: &mut [f64], n: usize, nass: usize, npiv: usize,
    iend_block: usize, uu: f64,
    perm: &mut [usize], nneg: &mut usize,
    must_eliminate_start: usize, must_eliminate_end: usize,
) -> (i32, usize) {
    let seuil: f64 = 0.0; // no static pivot threshold by default

    for ipiv in npiv..iend_block {
        let pivot = a[ipiv * n + ipiv]; // diagonal (column-major: a[col*n+row], col=ipiv, row=ipiv)

        // Compute AMAX: max off-diagonal in column ipiv within [npiv, iend_block)
        let mut amax = 0.0f64;
        let mut jmax: Option<usize> = None;
        // Scan rows npiv..ipiv-1 in column ipiv (below diagonal in column-major)
        for row in npiv..ipiv {
            let val = a[ipiv * n + row].abs();
            if val > amax { amax = val; jmax = Some(row); }
        }
        // Scan columns ipiv+1..iend_block-1 at row ipiv (by symmetry: column col, row ipiv)
        for col in (ipiv + 1)..iend_block {
            let val = a[col * n + ipiv].abs();
            if val > amax { amax = val; jmax = Some(col); }
        }

        // Compute RMAX: max off-diagonal outside inner block (PIVOT_OPTION=3)
        let mut rmax = 0.0f64;
        for col in iend_block..n {
            let val = a[col * n + ipiv].abs();
            rmax = rmax.max(val);
        }

        // Null check
        let col_max = amax.max(rmax).max(pivot.abs());
        if col_max <= 1e-30 { continue; }

        // 1x1 pivot test: |diag| >= uu * max(AMAX, RMAX) and |diag| > seuil
        #[cfg(debug_assertions)]
        if n <= 10 {
            eprintln!("  pivot search: ipiv={} diag={:.4} amax={:.4} rmax={:.4} jmax={:?} test={}",
                ipiv, pivot, amax, rmax, jmax, pivot.abs() >= uu * amax.max(rmax));
        }
        if pivot.abs() >= uu * amax.max(rmax)
            && pivot.abs() > seuil.max(f64::MIN_POSITIVE)
        {
            if pivot < 0.0 { *nneg += 1; }
            #[cfg(debug_assertions)]
            if n <= 10 { eprintln!("  -> 1x1 accepted at {}, nneg={}", ipiv, *nneg); }
            // Swap ipiv into position npiv
            if ipiv != npiv {
                symmetric_swap_ldlt(a, n, npiv, ipiv);
                perm.swap(npiv, ipiv);
            }
            return (0, 1);
        }

        // 2x2 pivot attempt
        if jmax.is_none() || npiv + 1 >= iend_block { continue; }
        let jmax_pos = jmax.unwrap();

        // Compute TMAX: max off-diagonal in column jmax, excluding ipiv
        let mut tmax = 0.0f64;
        for row in npiv..n {
            if row != ipiv && row != jmax_pos {
                let val = a[jmax_pos * n + row].abs();
                tmax = tmax.max(val);
            }
        }
        // Also check row jmax_pos in columns beyond jmax_pos
        for col in (jmax_pos + 1)..n {
            if col != ipiv {
                let val = a[col * n + jmax_pos].abs();
                tmax = tmax.max(val);
            }
        }
        tmax = tmax.max(seuil / uu.max(1e-30));

        let d_ii = a[ipiv * n + ipiv];
        let d_jj = a[jmax_pos * n + jmax_pos];
        let d_ij = if ipiv < jmax_pos {
            a[jmax_pos * n + ipiv]
        } else {
            a[ipiv * n + jmax_pos]
        };

        let detpiv = d_ii * d_jj - d_ij * d_ij;
        let abs_det = detpiv.abs();

        // 2x2 pivot test (modified Bunch-Kaufman)
        if abs_det <= 1e-30 { continue; }
        if (d_jj.abs() * rmax + amax * tmax) * uu > abs_det { continue; }
        if (d_ii.abs() * tmax + amax * rmax) * uu > abs_det { continue; }

        // 2x2 accepted! Count negative eigenvalues
        #[cfg(debug_assertions)]
        if n <= 10 {
            eprintln!("  -> 2x2 accepted: ({},{}) d_ii={:.4} d_jj={:.4} det={:.4} d_ij={:.4} nneg_before={}",
                ipiv, jmax_pos, d_ii, d_jj, detpiv, d_ij, *nneg);
        }
        if detpiv < 0.0 {
            *nneg += 1; // one positive + one negative
        } else if d_jj < 0.0 {
            *nneg += 2; // both negative
        }
        #[cfg(debug_assertions)]
        if n <= 10 { eprintln!("  -> nneg_after={}", *nneg); }

        // Swap: put min(ipiv, jmax) at npiv, max at npiv+1
        let first = ipiv.min(jmax_pos);
        let second = ipiv.max(jmax_pos);
        if first != npiv {
            symmetric_swap_ldlt(a, n, npiv, first);
            perm.swap(npiv, first);
        }
        if second != npiv + 1 {
            symmetric_swap_ldlt(a, n, npiv + 1, second);
            perm.swap(npiv + 1, second);
        }

        // Store DETPIV at sub-diagonal of pivot block: a[npiv, npiv+1] (col npiv, row npiv+1)
        a[npiv * n + (npiv + 1)] = detpiv;

        return (0, 2);
    }

    // No pivot found — check if ANY remaining column in [npiv, iend_block) is must-eliminate.
    // Must-eliminate columns are promoted from children or at the root supernode —
    // they cannot be delayed further. Swap the first must-eliminate column to npiv
    // and accept it via static pivoting.
    for pos in npiv..iend_block {
        if pos >= nass { break; }
        let orig_pos = perm[pos];
        if orig_pos >= must_eliminate_start && orig_pos < must_eliminate_end {
            // Swap this must-eliminate column to position npiv
            if pos != npiv {
                symmetric_swap_ldlt(a, n, npiv, pos);
                perm.swap(npiv, pos);
            }
            let pivot = a[npiv * n + npiv];
            if pivot < 0.0 { *nneg += 1; }
            return (-1, 1);
        }
    }

    #[cfg(debug_assertions)]
    if n <= 10 {
        eprintln!("  -> NO PIVOT: iend_block={} nass={}", iend_block, nass);
    }
    if iend_block >= nass {
        (1, 0) // exhausted all of NASS
    } else {
        (2, 0) // only this inner block exhausted
    }
}

/// MUMPS-style within-panel update (FAC_MQ_LDLT).
/// Performs rank-1 (1x1) or rank-2 (2x2) Schur complement update.
/// The update region is:
/// - Within inner block: triangular (JJ=1..I for column I)
/// - Beyond inner block: rectangular (JJ=1..NEL2 for column I > iend_block)
/// Original values are saved in pivot rows (U storage).
/// L entries (D^{-1}-scaled) are stored in pivot columns.
/// Returns IFINB: 0 (more in block), 1 (block done), -1 (NASS done).
fn update_within_block_ldlt(
    a: &mut [f64], n: usize, nass: usize, npiv: usize, pivsiz: usize,
    iend_block: usize, last_row: usize,
    d_diag: &mut [f64], d_offdiag: &mut [f64],
) -> i32 {
    let npiv_new = npiv + pivsiz;
    let nel2 = iend_block.saturating_sub(npiv_new); // remaining cols in inner block
    let ncb1 = last_row.saturating_sub(iend_block); // rows beyond inner block

    let ifinb = if nel2 == 0 {
        if iend_block >= nass { -1 } else { 1 }
    } else {
        0
    };

    if pivsiz == 1 {
        // 1x1 pivot at position npiv
        let apos_diag = npiv * n + npiv; // column npiv, row npiv
        let d = a[apos_diag];
        d_diag[npiv] = d;

        if d.abs() <= 1e-30 { return ifinb; }
        let valpiv = 1.0 / d;

        // Process each trailing column I = 1..nel2+ncb1
        for i in 1..=(nel2 + ncb1) {
            let col = npiv + i;
            let k1pos = col * n + npiv; // column col, row npiv (pivot row entry)

            // Save original to U storage (below diagonal in pivot column)
            a[npiv * n + (npiv + i)] = a[k1pos]; // row npiv+i in col npiv ← row npiv in col npiv+i

            // Wait: in column-major, a[col*n + row].
            // k1pos = col*n + npiv = entry(row=npiv, col=col) — the pivot ROW entry
            // a[npiv*n + (npiv+i)] = entry(row=npiv+i, col=npiv) — below diagonal in pivot COL
            // This saves the original row entry to the column (U storage).

            // Scale L entry: L[col, npiv] = original / d
            a[k1pos] *= valpiv;

            // Schur complement update
            let jmax = if i <= nel2 { i } else { nel2 };
            for jj in 1..=jmax {
                // a[col*n + (npiv+jj)] -= a[col*n + npiv] * a[npiv*n + (npiv+jj)]
                // = L_scaled[col] * U_original[jj]
                let target = col * n + (npiv + jj);
                let l_val = a[k1pos]; // L (scaled) for this column
                let u_val = a[npiv * n + (npiv + jj)]; // U (saved original) for row jj
                a[target] -= l_val * u_val;
            }
        }

    } else {
        // 2x2 pivot at positions npiv, npiv+1
        let d11 = a[npiv * n + npiv];
        let d22 = a[(npiv + 1) * n + (npiv + 1)];
        let detpiv = a[npiv * n + (npiv + 1)]; // stored by find_pivot_ldlt
        let d12_orig = a[(npiv + 1) * n + npiv]; // off-diagonal (lower triangle)

        d_diag[npiv] = d11;
        d_diag[npiv + 1] = d22;
        d_offdiag[npiv] = d12_orig;

        if detpiv.abs() <= 1e-30 { return ifinb; }

        // D^{-1} computation (MUMPS swaps indices: A11_inv = d22/det, A22_inv = d11/det)
        let a11_inv = d22 / detpiv;
        let a22_inv = d11 / detpiv;
        let a12_inv = -d12_orig / detpiv;

        // Fix storage: move off-diagonal to lower triangle, clear upper
        a[(npiv + 1) * n + npiv] = d12_orig; // lower: a[col=npiv, row=npiv+1]
        // Actually this is already there. Clear upper:
        // a[(npiv+1)*n + npiv] is col=npiv+1, row=npiv — wait, that's the upper position.
        // Let me re-check: a[col*n + row]. a[npiv*n + (npiv+1)] = col=npiv, row=npiv+1 = LOWER.
        // a[(npiv+1)*n + npiv] = col=npiv+1, row=npiv = UPPER.
        // The DETPIV was stored at a[npiv*n + (npiv+1)] (lower position).
        // The original d12 is at a[(npiv+1)*n + npiv] (upper position).
        // MUMPS wants: lower = d12, upper = 0
        a[npiv * n + (npiv + 1)] = d12_orig; // overwrite DETPIV with d12 in lower
        a[(npiv + 1) * n + npiv] = 0.0;       // clear upper

        // Process each trailing column
        for i in 1..=(nel2 + ncb1) {
            let col = npiv + 1 + i;
            let k1 = col * n + npiv;     // entry(row=npiv, col=col) — pivot row 1
            let k2 = col * n + (npiv + 1); // entry(row=npiv+1, col=col) — pivot row 2

            let a1_orig = a[k1];
            let a2_orig = a[k2];

            // Compute L*D^{-1} (negated for subtraction)
            let mult1 = -(a11_inv * a1_orig + a12_inv * a2_orig);
            let mult2 = -(a12_inv * a1_orig + a22_inv * a2_orig);

            // Save originals to U storage (pivot rows, below diagonal)
            // U row 1: a[npiv*n + col] = entry(row=col, col=npiv) — not right...
            // In column-major: to save in the pivot column below the diagonal:
            // Pivot 1 row storage: a[npiv*n + (npiv+1+i)] = entry(row=npiv+1+i, col=npiv)
            // Pivot 2 row storage: a[(npiv+1)*n + (npiv+1+i)] = entry(row=npiv+1+i, col=npiv+1)
            a[npiv * n + (npiv + 1 + i)] = a1_orig;
            a[(npiv + 1) * n + (npiv + 1 + i)] = a2_orig;

            // Schur complement update
            let jmax = if i <= nel2 { i } else { nel2 };
            for jj in 1..=jmax {
                let target = col * n + (npiv + 1 + jj); // entry(row=npiv+1+jj, col=col)
                let uk1 = a[npiv * n + (npiv + 1 + jj)];     // U row 1 (saved original)
                let uk2 = a[(npiv + 1) * n + (npiv + 1 + jj)]; // U row 2 (saved original)
                a[target] += mult1 * uk1 + mult2 * uk2;
            }

            // Store L*D^{-1} back in the column entries
            a[k1] = -mult1;
            a[k2] = -mult2;
        }
    }

    ifinb
}

/// MUMPS-style between-block GEMM update (FAC_SQ_LDLT).
/// Applies accumulated pivots from [ibeg_block, npiv) to trailing columns
/// [iend_block, last_col_gemm) and rows up to last_row_gemm.
/// U (original values) are in pivot rows, L (D^{-1}-scaled) in pivot columns.
///
/// With the `faer` feature, uses faer's optimized GEMM kernel (SIMD+FMA) for
/// precision matching MUMPS's BLAS-based between-block update.
fn between_block_gemm_ldlt(
    a: &mut [f64], n: usize,
    ibeg_block: usize, iend_block: usize, npiv: usize,
    last_col_gemm: usize, last_row_gemm: usize,
) {
    let npiv_block = npiv - ibeg_block;
    if npiv_block == 0 { return; }
    let nel1 = last_col_gemm.saturating_sub(iend_block);
    if nel1 == 0 { return; }

    #[cfg(feature = "faer")]
    {
        between_block_gemm_faer(a, n, ibeg_block, iend_block, npiv, last_col_gemm, last_row_gemm);
    }
    #[cfg(not(feature = "faer"))]
    {
        between_block_gemm_naive(a, n, ibeg_block, iend_block, npiv, last_col_gemm, last_row_gemm);
    }
}

/// Naive triple-loop GEMM fallback (no faer).
#[cfg(not(feature = "faer"))]
fn between_block_gemm_naive(
    a: &mut [f64], n: usize,
    ibeg_block: usize, iend_block: usize, npiv: usize,
    last_col_gemm: usize, last_row_gemm: usize,
) {
    let npiv_block = npiv - ibeg_block;

    // Symmetric part: update [iend_block, last_col_gemm) x [iend_block, last_col_gemm)
    for col in iend_block..last_col_gemm {
        for row in iend_block..last_col_gemm {
            let mut sum = 0.0;
            for k in ibeg_block..npiv {
                sum += a[k * n + row] * a[col * n + k];
            }
            a[col * n + row] -= sum;
        }
    }

    // Rectangular part: [iend_block, last_col_gemm) x [last_col_gemm, last_row_gemm)
    if last_row_gemm > last_col_gemm {
        for col in iend_block..last_col_gemm {
            for row in last_col_gemm..last_row_gemm {
                let mut sum = 0.0;
                for k in ibeg_block..npiv {
                    sum += a[k * n + row] * a[col * n + k];
                }
                a[col * n + row] -= sum;
            }
        }
    }
}

/// faer-accelerated GEMM for between-block Schur complement update.
/// Uses faer's optimized matmul kernel (SIMD + FMA) matching MUMPS's BLAS quality.
#[cfg(feature = "faer")]
fn between_block_gemm_faer(
    a: &mut [f64], n: usize,
    ibeg_block: usize, iend_block: usize, npiv: usize,
    last_col_gemm: usize, last_row_gemm: usize,
) {
    let npiv_block = npiv - ibeg_block;
    // "Symmetric" part: rows [iend_block, last_col_gemm) — must match naive version.
    // The rectangular part [last_col_gemm, last_row_gemm) is handled separately below.
    // Using last_row_gemm here would double-update the rectangular rows.
    let update_rows = last_col_gemm - iend_block;
    let update_cols = last_col_gemm - iend_block;
    if update_rows == 0 || update_cols == 0 { return; }

    // Copy U, L, target into faer Mats (column-major).
    // U: (update_rows x npiv_block), from pivot cols at target rows
    // L: (npiv_block x update_cols), from target cols at pivot rows
    let mut u_mat = faer::Mat::<f64>::zeros(update_rows, npiv_block);
    let mut l_mat = faer::Mat::<f64>::zeros(npiv_block, update_cols);
    let mut target = faer::Mat::<f64>::zeros(update_rows, update_cols);

    for k in 0..npiv_block {
        let src_col = ibeg_block + k;
        for i in 0..update_rows {
            u_mat[(i, k)] = a[src_col * n + (iend_block + i)];
        }
    }
    for j in 0..update_cols {
        let src_col = iend_block + j;
        for k in 0..npiv_block {
            l_mat[(k, j)] = a[src_col * n + (ibeg_block + k)];
        }
    }
    for j in 0..update_cols {
        for i in 0..update_rows {
            target[(i, j)] = a[(iend_block + j) * n + (iend_block + i)];
        }
    }

    // target = 1.0 * target + (-1.0) * U * L
    faer::linalg::matmul::matmul(
        target.as_mut(),
        u_mat.as_ref(),
        l_mat.as_ref(),
        Some(1.0),
        -1.0,
        faer::Parallelism::None,
    );

    // Write back
    for j in 0..update_cols {
        for i in 0..update_rows {
            a[(iend_block + j) * n + (iend_block + i)] = target[(i, j)];
        }
    }

    // Rectangular part: rows [last_col_gemm, last_row_gemm) x cols [iend_block, last_col_gemm)
    if last_row_gemm > last_col_gemm {
        let rect_rows = last_row_gemm - last_col_gemm;
        let nel1 = last_col_gemm - iend_block;

        let mut u_rect = faer::Mat::<f64>::zeros(rect_rows, npiv_block);
        let mut rect_target = faer::Mat::<f64>::zeros(rect_rows, nel1);

        for k in 0..npiv_block {
            for i in 0..rect_rows {
                u_rect[(i, k)] = a[(ibeg_block + k) * n + (last_col_gemm + i)];
            }
        }
        for j in 0..nel1 {
            for i in 0..rect_rows {
                rect_target[(i, j)] = a[(iend_block + j) * n + (last_col_gemm + i)];
            }
        }

        let l_rect = l_mat.get(.., ..nel1);
        faer::linalg::matmul::matmul(
            rect_target.as_mut(),
            u_rect.as_ref(),
            l_rect,
            Some(1.0),
            -1.0,
            faer::Parallelism::None,
        );

        for j in 0..nel1 {
            for i in 0..rect_rows {
                a[(iend_block + j) * n + (last_col_gemm + i)] = rect_target[(i, j)];
            }
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
