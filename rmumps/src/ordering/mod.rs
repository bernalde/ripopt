/// Approximate minimum degree ordering.
pub mod amd;
/// KKT-aware matching ordering for saddle-point systems.
pub mod kkt_matching;
/// Nested dissection ordering.
pub mod nested_dissection;
/// Identity (no reordering) permutation.
pub mod natural;

use crate::csc::CscMatrix;

/// Fill-reducing ordering methods.
#[derive(Debug, Clone, Copy)]
pub enum Ordering {
    /// Approximate Minimum Degree ordering (recommended).
    Amd,
    /// No reordering (identity permutation). Useful for testing.
    Natural,
    /// Nested dissection ordering (good for 2D/3D mesh problems).
    NestedDissection,
    /// KKT matching + AMD: pairs primal-dual variables for 2×2 pivots,
    /// then applies AMD. Requires `n_primal` to be set in SolverOptions.
    KktMatchingAmd,
}

/// Compute a fill-reducing permutation for the given symmetric CSC matrix.
/// Returns `(perm, perm_inv)` where `perm[new] = old` and `perm_inv[old] = new`.
pub fn compute_ordering(csc: &CscMatrix, ordering: Ordering) -> (Vec<usize>, Vec<usize>) {
    compute_ordering_with_kkt(csc, ordering, None)
}

/// Compute ordering with optional KKT structure info.
/// `n_primal` is required for `KktMatchingAmd` ordering.
pub fn compute_ordering_with_kkt(
    csc: &CscMatrix,
    ordering: Ordering,
    n_primal: Option<usize>,
) -> (Vec<usize>, Vec<usize>) {
    match ordering {
        Ordering::Amd => amd::amd_ordering(csc),
        Ordering::Natural => natural::natural_ordering(csc.n),
        Ordering::NestedDissection => nested_dissection::nd_ordering(csc),
        Ordering::KktMatchingAmd => {
            let np = n_primal.unwrap_or(0);
            if np > 0 && np < csc.n {
                kkt_matching::kkt_matching_ordering(csc, np)
            } else {
                amd::amd_ordering(csc)
            }
        }
    }
}

/// Apply a symmetric permutation to a CSC matrix: P * A * P^T.
/// `perm_inv[old]` = new position.
/// Returns a new CSC matrix with permuted indices.
pub fn permute_symmetric_csc(csc: &CscMatrix, _perm: &[usize], perm_inv: &[usize]) -> CscMatrix {
    let n = csc.n;

    // Collect all upper-triangle entries under the new ordering
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(csc.nnz());

    for old_col in 0..n {
        let new_col = perm_inv[old_col];
        for idx in csc.col_ptr[old_col]..csc.col_ptr[old_col + 1] {
            let old_row = csc.row_idx[idx];
            let new_row = perm_inv[old_row];
            let val = csc.vals[idx];

            // Ensure upper triangle (row <= col)
            let (r, c) = if new_row <= new_col {
                (new_row, new_col)
            } else {
                (new_col, new_row)
            };
            triplets.push((r, c, val));
        }
    }

    // Build CSC from the permuted triplets
    // Count per column
    let mut col_count = vec![0usize; n];
    for &(_, c, _) in &triplets {
        col_count[c] += 1;
    }
    let mut col_ptr = vec![0usize; n + 1];
    for j in 0..n {
        col_ptr[j + 1] = col_ptr[j] + col_count[j];
    }
    let nnz = col_ptr[n];
    let mut row_idx = vec![0usize; nnz];
    let mut vals = vec![0.0; nnz];
    let mut offset = col_ptr[..n].to_vec();

    for &(r, c, v) in &triplets {
        let pos = offset[c];
        row_idx[pos] = r;
        vals[pos] = v;
        offset[c] += 1;
    }

    // Sort rows within each column
    for j in 0..n {
        let start = col_ptr[j];
        let end = col_ptr[j + 1];
        // Simple insertion sort
        for i in (start + 1)..end {
            let key_row = row_idx[i];
            let key_val = vals[i];
            let mut pos = i;
            while pos > start && row_idx[pos - 1] > key_row {
                row_idx[pos] = row_idx[pos - 1];
                vals[pos] = vals[pos - 1];
                pos -= 1;
            }
            row_idx[pos] = key_row;
            vals[pos] = key_val;
        }
    }

    CscMatrix { n, col_ptr, row_idx, vals }
}
