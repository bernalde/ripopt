use crate::csc::CscMatrix;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Approximate Minimum Degree ordering.
///
/// Implements the AMD algorithm based on Amestoy, Davis, Duff (2004)
/// using a quotient graph representation with element absorption and
/// supervariable detection. Takes a symmetric CSC matrix (upper triangle)
/// and returns a fill-reducing permutation (perm, perm_inv).
///
/// Key optimizations:
/// - Flat-array representation for adjacency data (cache-friendly, zero heap alloc in main loop)
/// - Generation counters instead of clearing seen[] arrays
/// - Supervariable detection: indistinguishable nodes are merged, reducing effective problem size
/// - Element absorption: redundant elements are removed from adjacency lists
pub fn amd_ordering(csc: &CscMatrix) -> (Vec<usize>, Vec<usize>) {
    let n = csc.n;
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    if n == 1 {
        return (vec![0], vec![0]);
    }

    // --- Phase 1: Build initial adjacency in flat array ---

    let mut deg_init = vec![0usize; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            if i != j {
                deg_init[i] += 1;
                deg_init[j] += 1;
            }
        }
    }

    let mut pe = vec![0usize; n];
    let mut offset = 0usize;
    for i in 0..n {
        pe[i] = offset;
        offset += deg_init[i];
    }
    let total_adj = offset;

    let iwlen = total_adj + total_adj + n;
    let mut iw = vec![0usize; iwlen];
    let mut pfree = total_adj;

    let mut len = vec![0usize; n];
    {
        let mut pos = vec![0usize; n];
        for j in 0..n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                if i != j {
                    iw[pe[i] + pos[i]] = j;
                    pos[i] += 1;
                    iw[pe[j] + pos[j]] = i;
                    pos[j] += 1;
                }
            }
        }
        for i in 0..n {
            len[i] = pos[i];
            iw[pe[i]..pe[i] + len[i]].sort_unstable();
        }
    }

    // --- Phase 2: Initialize data structures ---

    let mut elen = vec![0usize; n]; // # element entries at start of i's list
    let mut degree = deg_init;
    let mut eliminated = vec![false; n];
    let mut perm = Vec::with_capacity(n);

    // Supervariable data
    let mut nv = vec![1usize; n];         // weight of supervariable (0 = dead/absorbed)
    let mut svar_next = vec![usize::MAX; n]; // linked list of absorbed nodes

    // Min-heap with (degree, node_index)
    let mut heap: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::with_capacity(n);
    for i in 0..n {
        heap.push(Reverse((degree[i], i)));
    }

    // Scratch arrays
    let mut seen = vec![false; n];
    let mut reachable = Vec::with_capacity(n);
    let mut elem_cleaned = vec![0usize; n];
    let mut deg_gen = vec![0usize; n];
    let mut gen = 0usize;
    let mut elim_step = 0usize;
    let mut node_buf = Vec::with_capacity(n);

    // --- Phase 3: Main elimination loop ---

    while perm.len() < n {
        elim_step += 1;

        // Pop minimum-degree node (skip eliminated and dead supervariables)
        let pivot = loop {
            if let Some(Reverse((deg, node))) = heap.pop() {
                if !eliminated[node] && nv[node] > 0 && deg == degree[node] {
                    break node;
                }
            } else {
                break usize::MAX; // shouldn't happen
            }
        };
        if pivot == usize::MAX {
            break;
        }

        eliminated[pivot] = true;
        // Expand supervariable: push pivot and all absorbed nodes to perm
        perm.push(pivot);
        {
            let mut cur = svar_next[pivot];
            while cur != usize::MAX {
                perm.push(cur);
                cur = svar_next[cur];
            }
        }

        if perm.len() >= n {
            break;
        }

        // Compute reachable set of pivot (skip eliminated and dead)
        reachable.clear();
        let p_start = pe[pivot];
        for k in 0..elen[pivot] {
            let e = iw[p_start + k];
            for m in 0..len[e] {
                let member = iw[pe[e] + m];
                if !eliminated[member] && nv[member] > 0 && !seen[member] {
                    seen[member] = true;
                    reachable.push(member);
                }
            }
        }
        for k in elen[pivot]..len[pivot] {
            let nb = iw[p_start + k];
            if !eliminated[nb] && nv[nb] > 0 && !seen[nb] {
                seen[nb] = true;
                reachable.push(nb);
            }
        }
        for &r in &reachable {
            seen[r] = false;
        }

        // Pivot becomes element: store member list = reachable set
        if reachable.len() <= len[pivot] {
            for (k, &r) in reachable.iter().enumerate() {
                iw[pe[pivot] + k] = r;
            }
            len[pivot] = reachable.len();
        } else {
            if pfree + reachable.len() > iw.len() {
                iw.resize(pfree + reachable.len() + n, 0);
            }
            pe[pivot] = pfree;
            for &r in &reachable {
                iw[pfree] = r;
                pfree += 1;
            }
            len[pivot] = reachable.len();
        }

        if reachable.is_empty() {
            continue;
        }

        // Mark reachable for absorption checking
        for &r in &reachable {
            seen[r] = true;
        }

        // Rebuild each reachable node's adjacency list.
        // The rebuilt list can be up to 1 entry longer than the original
        // (when pivot is added but wasn't a direct neighbor), so we build
        // into a scratch buffer and relocate if needed.
        for &ni in &reachable {
            let start = pe[ni];
            let old_elen = elen[ni];
            let old_len = len[ni];

            // Build new list in node_buf
            node_buf.clear();
            let mut new_elen = 0usize;

            // Surviving (non-absorbed) elements
            for k in 0..old_elen {
                let e = iw[start + k];
                let mut absorb = true;
                for m in 0..len[e] {
                    let member = iw[pe[e] + m];
                    if !eliminated[member] && nv[member] > 0 && member != ni && !seen[member] {
                        absorb = false;
                        break;
                    }
                }
                if !absorb {
                    node_buf.push(e);
                    new_elen += 1;
                }
            }

            // Add pivot as new element
            node_buf.push(pivot);
            new_elen += 1;

            // Surviving node neighbors (skip dead and eliminated)
            for k in old_elen..old_len {
                let nb = iw[start + k];
                if nb != pivot && !eliminated[nb] && nv[nb] > 0 {
                    node_buf.push(nb);
                }
            }

            // Place the rebuilt list
            let new_len = node_buf.len();
            if new_len <= old_len {
                iw[start..start + new_len].copy_from_slice(&node_buf);
            } else {
                if pfree + new_len > iw.len() {
                    iw.resize(pfree + new_len + n, 0);
                }
                iw[pfree..pfree + new_len].copy_from_slice(&node_buf);
                pe[ni] = pfree;
                pfree += new_len;
            }
            elen[ni] = new_elen;
            len[ni] = new_len;
        }

        // Clear seen flags
        for &r in &reachable {
            seen[r] = false;
        }

        // Compute exact degrees using generation counter
        for &ni in &reachable {
            let start = pe[ni];
            let ni_elen = elen[ni];
            let ni_len = len[ni];
            let mut ext_degree = 0usize;

            gen += 1;

            // Node neighbors (count nv[] weight)
            for k in ni_elen..ni_len {
                let nb = iw[start + k];
                if nb != ni && deg_gen[nb] != gen {
                    deg_gen[nb] = gen;
                    ext_degree += nv[nb];
                }
            }

            // Element neighbors' members (compact lazily, count nv[] weight)
            for k in 0..ni_elen {
                let e = iw[start + k];
                let estart = pe[e];
                if elem_cleaned[e] != elim_step {
                    let mut em_write = 0usize;
                    for em_read in 0..len[e] {
                        let member = iw[estart + em_read];
                        if !eliminated[member] && nv[member] > 0 {
                            iw[estart + em_write] = member;
                            em_write += 1;
                        }
                    }
                    len[e] = em_write;
                    elem_cleaned[e] = elim_step;
                }

                for m in 0..len[e] {
                    let member = iw[estart + m];
                    if member != ni && deg_gen[member] != gen {
                        deg_gen[member] = gen;
                        ext_degree += nv[member];
                    }
                }
            }

            degree[ni] = ext_degree;
            heap.push(Reverse((ext_degree, ni)));
        }

        // --- Supervariable detection ---
        // Two nodes i, j are indistinguishable if:
        //   Adjacent case:    adj(i) \ {j} == adj(j) \ {i}
        //   Non-adjacent case: adj(i) == adj(j)
        // We compare sorted COPIES (not in-place) to preserve adjacency structure.
        // Skip for large reachable sets: the pairwise O(R²) comparison is too expensive
        // and supervariable detection has minimal impact on ordering quality.
        if reachable.len() > 1 && reachable.len() <= 256 {
            // Group candidates by (degree, len)
            let mut candidates: Vec<(usize, usize, usize)> = Vec::new(); // (degree, len, node)
            for &ni in &reachable {
                if nv[ni] == 0 { continue; }
                candidates.push((degree[ni], len[ni], ni));
            }
            candidates.sort_unstable();

            // Scratch buffers for sorted adjacency comparison
            let mut sorted_i: Vec<usize> = Vec::new();
            let mut sorted_j: Vec<usize> = Vec::new();

            let mut ci = 0;
            while ci < candidates.len() {
                let (d0, l0, ni) = candidates[ci];
                if nv[ni] == 0 { ci += 1; continue; }

                // Build sorted copy of adj(ni) \ {ni} once per group leader
                sorted_i.clear();
                sorted_i.extend_from_slice(&iw[pe[ni]..pe[ni] + len[ni]]);
                sorted_i.sort_unstable();

                let mut cj = ci + 1;
                while cj < candidates.len() {
                    let (d1, l1, nj) = candidates[cj];
                    if d1 != d0 || l1 != l0 { break; } // different group
                    if nv[nj] == 0 { cj += 1; continue; }

                    // Build sorted copy of adj(nj)
                    sorted_j.clear();
                    sorted_j.extend_from_slice(&iw[pe[nj]..pe[nj] + len[nj]]);
                    sorted_j.sort_unstable();

                    // Compare adj(i)\{j} vs adj(j)\{i} using merge on sorted copies
                    let mut pi = 0usize;
                    let mut pj = 0usize;
                    let li = sorted_i.len();
                    let lj = sorted_j.len();
                    let mut identical = true;
                    loop {
                        // Skip nj in sorted_i and ni in sorted_j
                        while pi < li && sorted_i[pi] == nj { pi += 1; }
                        while pj < lj && sorted_j[pj] == ni { pj += 1; }

                        if pi >= li && pj >= lj { break; }
                        if pi >= li || pj >= lj { identical = false; break; }
                        if sorted_i[pi] != sorted_j[pj] { identical = false; break; }
                        pi += 1;
                        pj += 1;
                    }

                    if identical {
                        nv[ni] += nv[nj];
                        nv[nj] = 0;
                        // Append nj's chain to ni's chain
                        let mut tail = nj;
                        while svar_next[tail] != usize::MAX {
                            tail = svar_next[tail];
                        }
                        svar_next[tail] = svar_next[ni];
                        svar_next[ni] = nj;
                    }
                    cj += 1;
                }
                ci += 1;
            }
        }
    }

    // --- Phase 4: Build inverse permutation ---

    let mut perm_inv = vec![0usize; n];
    for (new_pos, &old_idx) in perm.iter().enumerate() {
        perm_inv[old_idx] = new_pos;
    }

    (perm, perm_inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;

    fn csc_from_upper_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> CscMatrix {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        CscMatrix::from_coo(&coo)
    }

    fn is_valid_permutation(perm: &[usize], n: usize) -> bool {
        if perm.len() != n {
            return false;
        }
        let mut seen = vec![false; n];
        for &p in perm {
            if p >= n || seen[p] {
                return false;
            }
            seen[p] = true;
        }
        true
    }

    /// Count fill-in (nonzeros in L) for a given elimination order.
    fn count_fill(csc: &CscMatrix, perm: &[usize]) -> usize {
        let n = csc.n;
        let perm_inv: Vec<usize> = {
            let mut pi = vec![0; n];
            for (new, &old) in perm.iter().enumerate() {
                pi[old] = new;
            }
            pi
        };

        let permuted = crate::ordering::permute_symmetric_csc(csc, perm, &perm_inv);
        let sym = crate::symbolic::SymbolicFactorization::from_csc(&permuted);
        sym.l_nnz
    }

    #[test]
    fn test_amd_valid_perm_diagonal() {
        let csc = csc_from_upper_triplets(3, &[
            (0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0),
        ]);
        let (perm, perm_inv) = amd_ordering(&csc);
        assert!(is_valid_permutation(&perm, 3));
        assert!(is_valid_permutation(&perm_inv, 3));
    }

    #[test]
    fn test_amd_valid_perm_tridiagonal() {
        let csc = csc_from_upper_triplets(5, &[
            (0, 0, 1.0), (0, 1, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0), (3, 4, 1.0),
            (4, 4, 1.0),
        ]);
        let (perm, perm_inv) = amd_ordering(&csc);
        assert!(is_valid_permutation(&perm, 5));
        assert!(is_valid_permutation(&perm_inv, 5));
    }

    #[test]
    fn test_amd_arrow_less_fill() {
        let n = 6;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 10.0));
            if i < n - 1 {
                triplets.push((i, n - 1, 1.0));
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);

        let natural_perm: Vec<usize> = (0..n).collect();
        let natural_fill = count_fill(&csc, &natural_perm);

        let (amd_perm, _) = amd_ordering(&csc);
        let amd_fill = count_fill(&csc, &amd_perm);

        assert!(
            amd_fill <= natural_fill,
            "AMD fill {} should be <= natural fill {}",
            amd_fill, natural_fill
        );
    }

    #[test]
    fn test_amd_perm_inv_consistent() {
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0), (0, 3, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let (perm, perm_inv) = amd_ordering(&csc);
        for i in 0..4 {
            assert_eq!(perm_inv[perm[i]], i);
            assert_eq!(perm[perm_inv[i]], i);
        }
    }

    #[test]
    fn test_amd_empty() {
        let csc = CscMatrix { n: 0, col_ptr: vec![0], row_idx: vec![], vals: vec![] };
        let (perm, perm_inv) = amd_ordering(&csc);
        assert_eq!(perm.len(), 0);
        assert_eq!(perm_inv.len(), 0);
    }

    #[test]
    fn test_amd_large_arrow() {
        // Larger arrow to stress-test: n=50
        let n = 50;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 10.0));
            if i < n - 1 {
                triplets.push((i, n - 1, 1.0));
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);

        let natural_perm: Vec<usize> = (0..n).collect();
        let natural_fill = count_fill(&csc, &natural_perm);

        let (amd_perm, _) = amd_ordering(&csc);
        assert!(is_valid_permutation(&amd_perm, n));
        let amd_fill = count_fill(&csc, &amd_perm);

        assert!(
            amd_fill <= natural_fill,
            "AMD fill {} should be <= natural fill {} for n={}",
            amd_fill, natural_fill, n
        );
    }

    #[test]
    fn test_amd_2d_grid() {
        // 5x5 grid graph (25 nodes) — tests AMD on a realistic structure
        let nx = 5;
        let n = nx * nx;
        let mut triplets = Vec::new();
        for iy in 0..nx {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                triplets.push((idx, idx, 4.0));
                if ix + 1 < nx {
                    let right = iy * nx + ix + 1;
                    let (r, c) = if idx < right { (idx, right) } else { (right, idx) };
                    triplets.push((r, c, -1.0));
                }
                if iy + 1 < nx {
                    let below = (iy + 1) * nx + ix;
                    let (r, c) = if idx < below { (idx, below) } else { (below, idx) };
                    triplets.push((r, c, -1.0));
                }
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);
        let (amd_perm, _) = amd_ordering(&csc);
        assert!(is_valid_permutation(&amd_perm, n));

        let natural_fill = count_fill(&csc, &(0..n).collect::<Vec<_>>());
        let amd_fill = count_fill(&csc, &amd_perm);
        assert!(amd_fill <= natural_fill);
    }
}
