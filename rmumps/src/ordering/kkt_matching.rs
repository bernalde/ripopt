//! KKT-aware matching ordering for saddle-point systems.
//!
//! For a KKT system with n primal variables and m dual variables:
//! ```text
//! [H + Σ    J^T]
//! [J       -D_c]
//! ```
//! This ordering pairs each dual variable y_i with the primal variable x_j
//! that has the largest |J_{ij}|. The paired columns are placed adjacent in
//! the elimination order, enabling Bunch-Kaufman to form 2×2 pivots that
//! capture the primal-dual coupling directly. This avoids the catastrophic
//! cancellation that occurs when 1×1 pivots are used on the tiny D_c entries.
//!
//! The ordering uses a compressed-graph AMD approach (like MUMPS ICNTL(12)):
//! matched pairs are merged into single super-nodes, AMD is run on the
//! compressed graph, then each super-node is expanded back into its two
//! columns (primal first, dual second). This preserves pair adjacency
//! while getting good fill-reduction from AMD.

use crate::csc::CscMatrix;

/// Compute a KKT matching + compressed AMD ordering.
///
/// `n_primal` is the number of primal variables (the (1,1) block size).
/// The matrix has dimension `n_primal + m_dual`.
///
/// Returns `(perm, perm_inv)` where `perm[new] = old`.
pub fn kkt_matching_ordering(csc: &CscMatrix, n_primal: usize) -> (Vec<usize>, Vec<usize>) {
    let dim = csc.n;
    let m_dual = dim - n_primal;

    if m_dual == 0 || n_primal == 0 {
        return super::amd::amd_ordering(csc);
    }

    // Step 1: Greedy maximum-weight matching
    let mut match_dual_to_primal: Vec<Option<usize>> = vec![None; m_dual];
    let mut match_primal_to_dual: Vec<Option<usize>> = vec![None; n_primal];

    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    for dual_idx in 0..m_dual {
        let col = n_primal + dual_idx;
        for idx in csc.col_ptr[col]..csc.col_ptr[col + 1] {
            let row = csc.row_idx[idx];
            if row < n_primal {
                let weight = csc.vals[idx].abs();
                if weight > 0.0 {
                    edges.push((dual_idx, row, weight));
                }
            }
        }
    }
    edges.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    let mut num_matched = 0usize;
    for &(dual_idx, primal_idx, _) in &edges {
        if match_dual_to_primal[dual_idx].is_none() && match_primal_to_dual[primal_idx].is_none() {
            match_dual_to_primal[dual_idx] = Some(primal_idx);
            match_primal_to_dual[primal_idx] = Some(dual_idx);
            num_matched += 1;
        }
    }

    // Step 2: Build compressed graph
    // Each matched pair becomes a single super-node.
    // Unmatched primal and dual variables are individual nodes.
    //
    // Compressed node numbering:
    //   0..num_matched: matched pairs (super-node k represents pair k)
    //   num_matched..num_matched+num_unmatched_primal: unmatched primals
    //   num_matched+num_unmatched_primal..: unmatched duals
    let num_unmatched_primal = n_primal - num_matched;
    let num_unmatched_dual = m_dual - num_matched;
    let compressed_n = num_matched + num_unmatched_primal + num_unmatched_dual;

    // Map original column → compressed node
    let mut orig_to_compressed = vec![0usize; dim];
    // Map compressed node → original columns (1 or 2)
    let mut compressed_to_orig: Vec<Vec<usize>> = Vec::with_capacity(compressed_n);

    // Matched pairs
    let mut pair_idx = 0;
    for dual_idx in 0..m_dual {
        if let Some(primal_idx) = match_dual_to_primal[dual_idx] {
            orig_to_compressed[primal_idx] = pair_idx;
            orig_to_compressed[n_primal + dual_idx] = pair_idx;
            compressed_to_orig.push(vec![primal_idx, n_primal + dual_idx]);
            pair_idx += 1;
        }
    }

    // Unmatched primals
    let mut unmatched_idx = num_matched;
    for j in 0..n_primal {
        if match_primal_to_dual[j].is_none() {
            orig_to_compressed[j] = unmatched_idx;
            compressed_to_orig.push(vec![j]);
            unmatched_idx += 1;
        }
    }

    // Unmatched duals
    for i in 0..m_dual {
        if match_dual_to_primal[i].is_none() {
            orig_to_compressed[n_primal + i] = unmatched_idx;
            compressed_to_orig.push(vec![n_primal + i]);
            unmatched_idx += 1;
        }
    }
    assert_eq!(compressed_to_orig.len(), compressed_n);

    // Build compressed adjacency: for each compressed node, its neighbor set
    let mut adj: Vec<std::collections::BTreeSet<usize>> = vec![std::collections::BTreeSet::new(); compressed_n];
    for j in 0..dim {
        let cj = orig_to_compressed[j];
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            let ci = orig_to_compressed[i];
            if ci != cj {
                adj[cj].insert(ci);
                adj[ci].insert(cj);
            }
        }
    }

    // Build full symmetric CSC for AMD (including self-loops on diagonal)
    let mut comp_col_ptr = vec![0i64; compressed_n + 1];
    for j in 0..compressed_n {
        comp_col_ptr[j + 1] = comp_col_ptr[j] + adj[j].len() as i64;
    }
    let comp_nnz = comp_col_ptr[compressed_n] as usize;
    let mut comp_row_idx = vec![0i64; comp_nnz];
    for j in 0..compressed_n {
        let start = comp_col_ptr[j] as usize;
        for (k, &neighbor) in adj[j].iter().enumerate() {
            comp_row_idx[start + k] = neighbor as i64;
        }
    }

    // Step 3: Run AMD on compressed graph
    let comp_perm = if compressed_n <= 2 || comp_nnz == 0 {
        // Too small for AMD — use natural ordering
        (0..compressed_n).collect::<Vec<_>>()
    } else {
        let control = amd::Control::default();
        match amd::order::<i64>(
            compressed_n as i64, &comp_col_ptr, &comp_row_idx, &control,
        ) {
            Ok((perm_i64, _, _)) => {
                perm_i64.iter().map(|&x| x as usize).collect::<Vec<_>>()
            }
            Err(_) => (0..compressed_n).collect(),
        }
    };

    // Step 4: Expand compressed ordering back to original columns
    // Each compressed super-node expands to its original columns (primal first for pairs)
    let mut perm = Vec::with_capacity(dim);
    for &comp_idx in &comp_perm {
        let orig_cols = &compressed_to_orig[comp_idx];
        perm.extend_from_slice(orig_cols);
    }
    assert_eq!(perm.len(), dim);

    let mut perm_inv = vec![0usize; dim];
    for (new, &old) in perm.iter().enumerate() {
        perm_inv[old] = new;
    }

    (perm, perm_inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;

    #[test]
    fn test_kkt_matching_simple() {
        // 3x3 KKT: n=2 primal, m=1 dual
        let coo = CooMatrix::new(3,
            vec![0, 1, 0, 1],
            vec![0, 1, 2, 2],
            vec![2.0, 2.0, 1.0, 1.0],
        ).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let (perm, perm_inv) = kkt_matching_ordering(&csc, 2);

        assert_eq!(perm.len(), 3);
        let mut sorted = perm.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2]);
        for i in 0..3 {
            assert_eq!(perm_inv[perm[i]], i);
        }
    }

    #[test]
    fn test_kkt_matching_preserves_adjacency() {
        // 5x5 KKT: n=3, m=2 (need enough structure for AMD)
        // H diagonal + Jacobian coupling each dual to one primal
        let coo = CooMatrix::new(5,
            vec![0, 1, 2, 0, 1],
            vec![0, 1, 2, 3, 4],
            vec![2.0, 2.0, 2.0, 3.0, 4.0],
        ).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let (perm, perm_inv) = kkt_matching_ordering(&csc, 3);

        assert_eq!(perm.len(), 5);
        let mut sorted = perm.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4]);

        // Matched pairs should be adjacent: (0,3) and (1,4)
        let pos0 = perm_inv[0];
        let pos3 = perm_inv[3];
        assert!((pos0 as isize - pos3 as isize).unsigned_abs() == 1,
            "Primal 0 and dual 0 should be adjacent: pos0={}, pos3={}", pos0, pos3);

        let pos1 = perm_inv[1];
        let pos4 = perm_inv[4];
        assert!((pos1 as isize - pos4 as isize).unsigned_abs() == 1,
            "Primal 1 and dual 1 should be adjacent: pos1={}, pos4={}", pos1, pos4);
    }
}
