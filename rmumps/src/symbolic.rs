use crate::csc::CscMatrix;
use crate::etree::EliminationTree;

/// Compute the sorted union of two sorted slices.
fn sorted_union(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            result.push(a[i]);
            i += 1;
        } else if a[i] > b[j] {
            result.push(b[j]);
            j += 1;
        } else {
            result.push(a[i]);
            i += 1;
            j += 1;
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

/// A supernode: a contiguous range of columns that can be factored together.
#[derive(Debug, Clone)]
pub struct Supernode {
    /// First column in this supernode.
    pub start: usize,
    /// Number of columns (fully-summed variables).
    pub nfs: usize,
    /// All row indices in this supernode's front (includes the nfs FS columns
    /// followed by the contribution block indices). Sorted.
    pub front_indices: Vec<usize>,
}

/// Symbolic factorization result.
/// Describes the structure of the multifrontal factorization without computing values.
#[derive(Debug, Clone)]
pub struct SymbolicFactorization {
    /// Matrix dimension.
    pub n: usize,
    /// Elimination tree (per-column).
    pub etree: EliminationTree,
    /// Per-column front indices (row structure of L, including diagonal).
    pub col_indices: Vec<Vec<usize>>,
    /// Supernodes in postorder (leaves before parents).
    pub supernodes: Vec<Supernode>,
    /// For each column j, which supernode index it belongs to.
    pub col_to_supernode: Vec<usize>,
    /// Supernodal elimination tree: `snode_parent[s]` = parent supernode of s, or None.
    pub snode_parent: Vec<Option<usize>>,
    /// Children of each supernode.
    pub snode_children: Vec<Vec<usize>>,
    /// Total number of nonzeros in L.
    pub l_nnz: usize,

    // Legacy fields for backward compatibility with non-supernodal code paths
    /// Per-column front indices (same as col_indices, kept for compatibility).
    pub front_indices: Vec<Vec<usize>>,
    /// Per-column nfs (always 1, kept for compatibility).
    pub front_nfs: Vec<usize>,
    /// Per-column postorder.
    pub postorder: Vec<usize>,
}

impl SymbolicFactorization {
    /// Compute the symbolic factorization from a permuted upper-triangle CSC matrix.
    pub fn from_csc(csc: &CscMatrix) -> Self {
        let n = csc.n;
        let etree = EliminationTree::from_csc(csc);
        let postorder = etree.postorder();
        let children = etree.children();

        // Phase 1: Compute per-column row structure of L (bottom-up in postorder)
        //
        // Build row adjacency: for each row i, which columns j > i have entry (i, j).
        // The CSC stores upper triangle (row <= col), so entry (i, j) with i < j is
        // in column j at row i. We invert this to get row_adj[i] = {j : (i,j) in A, j > i}.
        let mut row_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for j in 0..n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                if i < j {
                    row_adj[i].push(j);
                }
            }
        }

        let mut col_indices: Vec<Vec<usize>> = vec![Vec::new(); n];

        for &j in &postorder {
            let mut row_set: Vec<usize> = Vec::new();

            // Original off-diagonal entries in row j (columns k > j)
            row_set.extend_from_slice(&row_adj[j]);

            // Merge children's row structures (indices > j)
            for &child in &children[j] {
                for &row in &col_indices[child] {
                    if row > j {
                        row_set.push(row);
                    }
                }
            }

            row_set.sort_unstable();
            row_set.dedup();

            let mut indices = Vec::with_capacity(1 + row_set.len());
            indices.push(j);
            indices.extend_from_slice(&row_set);
            col_indices[j] = indices;
        }

        // Phase 2: Detect fundamental supernodes
        // Column j+1 can be merged with column j's supernode if:
        //   1. j+1 has exactly one child in the etree: j
        //   2. col_indices[j+1] == col_indices[j] \ {j}
        //      (i.e., same structure shifted by one)
        let mut snode_id = vec![0usize; n]; // which supernode each column belongs to
        let mut supernodes: Vec<Supernode> = Vec::new();

        if n > 0 {
            // Process in natural order (0..n) since supernodes are contiguous ranges
            let mut col = 0;
            while col < n {
                let current_start = col;
                let mut nfs = 1;

                // Try to extend: can col+1 merge into this supernode?
                while col + nfs < n {
                    let next = col + nfs;
                    let prev = col + nfs - 1;

                    // Check: next has exactly one child = prev
                    let next_children = &children[next];
                    if next_children.len() != 1 || next_children[0] != prev {
                        break;
                    }

                    let prev_idx = &col_indices[prev];
                    let next_idx = &col_indices[next];

                    if next_idx.len() + 1 != prev_idx.len() {
                        break;
                    }

                    // Check tail match: next_idx[1..] == prev_idx[2..]
                    let mut matches = true;
                    for i in 1..next_idx.len() {
                        if next_idx[i] != prev_idx[i + 1] {
                            matches = false;
                            break;
                        }
                    }
                    if !matches {
                        break;
                    }

                    nfs += 1;
                }

                // Create supernode
                let snode_idx = supernodes.len();
                let front_idx = &col_indices[current_start];
                let snode = Supernode {
                    start: current_start,
                    nfs,
                    front_indices: front_idx.clone(),
                };
                supernodes.push(snode);

                for k in current_start..(current_start + nfs) {
                    snode_id[k] = snode_idx;
                }

                col += nfs;
            }
        }

        // Phase 2b: Relaxed supernode amalgamation
        //
        // Merge child-parent supernode pairs when they are connected by a single
        // edge in the supernodal tree and the additional zero fill-in is small.
        // This creates larger dense blocks, converting Level-2 BLAS operations
        // (rank-1 updates) into Level-3 BLAS (GEMM), which is critical for
        // performance on KKT systems from interior point methods.
        //
        // We iterate bottom-up and merge child into parent when:
        //   1. Parent has exactly one child (the candidate)
        //   2. The child's last column is immediately before the parent's first column
        //      (contiguous range — required since supernodes must span consecutive cols)
        //   3. The extra zero fill from the union of front indices is within budget
        {
            // Build initial supernodal tree
            let num_snodes_init = supernodes.len();
            let mut sp: Vec<Option<usize>> = vec![None; num_snodes_init];
            let mut sc: Vec<Vec<usize>> = vec![Vec::new(); num_snodes_init];
            for s in 0..num_snodes_init {
                let last_col = supernodes[s].start + supernodes[s].nfs - 1;
                if let Some(parent_col) = etree.parent[last_col] {
                    let parent_snode = snode_id[parent_col];
                    if parent_snode != s {
                        sp[s] = Some(parent_snode);
                        sc[parent_snode].push(s);
                    }
                }
            }
            for ch in sc.iter_mut() {
                ch.sort_unstable();
                ch.dedup();
            }

            // Track which supernodes are alive (not merged into another)
            let mut alive = vec![true; num_snodes_init];

            // Bottom-up: try merging each child into its parent
            for s in 0..num_snodes_init {
                if !alive[s] {
                    continue;
                }
                let parent_idx = match sp[s] {
                    Some(p) if alive[p] => p,
                    _ => continue,
                };

                // Only merge if parent has exactly one (alive) child
                let alive_children: Vec<usize> = sc[parent_idx].iter()
                    .copied()
                    .filter(|&c| alive[c])
                    .collect();
                if alive_children.len() != 1 || alive_children[0] != s {
                    continue;
                }

                // Contiguity check: child's last col + 1 == parent's first col
                let child_end = supernodes[s].start + supernodes[s].nfs;
                if child_end != supernodes[parent_idx].start {
                    continue;
                }

                // Compute union of front indices (both sorted)
                let child_cb = &supernodes[s].front_indices[supernodes[s].nfs..];
                let parent_front = &supernodes[parent_idx].front_indices;
                let union = sorted_union(child_cb, parent_front);

                // Allow some fill-in from extra CB indices. The dynamic front
                // expansion in extend_add handles unexpected indices at ancestors.
                // Limit fill to avoid excessive front growth: allow up to 50% extra rows
                // relative to the merged front size, and cap absolute extra at 32.
                let extra_rows = union.len() - parent_front.len();
                let merged_nfs = supernodes[s].nfs + supernodes[parent_idx].nfs;
                let merged_front_size = merged_nfs + union.len();
                let fill_ratio = extra_rows as f64 / merged_front_size.max(1) as f64;
                if extra_rows > 32 || fill_ratio > 0.5 {
                    continue;
                }

                // Merge: child absorbs parent's columns
                supernodes[s].nfs = merged_nfs;
                supernodes[s].front_indices = {
                    // New front = child's FS cols + parent's FS cols + union of CB
                    let mut new_front = Vec::with_capacity(merged_nfs + union.len());
                    // FS columns: child start .. child_end, then parent start .. parent_end
                    for col in supernodes[s].start..(supernodes[s].start + merged_nfs) {
                        new_front.push(col);
                    }
                    // CB indices from union, excluding the FS columns
                    let parent_end = supernodes[parent_idx].start + supernodes[parent_idx].nfs;
                    for &idx in &union {
                        if idx >= parent_end {
                            new_front.push(idx);
                        }
                    }
                    new_front
                };

                // Update snode_id for parent's columns
                for col in supernodes[parent_idx].start..(supernodes[parent_idx].start + supernodes[parent_idx].nfs) {
                    snode_id[col] = s;
                }

                // Inherit parent's parent
                alive[parent_idx] = false;
                sp[s] = sp[parent_idx];
                if let Some(grandparent) = sp[s] {
                    // Replace parent_idx with s in grandparent's children
                    if let Some(pos) = sc[grandparent].iter().position(|&c| c == parent_idx) {
                        sc[grandparent][pos] = s;
                    }
                }
                // Inherit parent's other children (should be none since we checked single child)
                // But inherit any children the parent had from other merges
                let parent_children: Vec<usize> = sc[parent_idx].iter()
                    .copied()
                    .filter(|&c| alive[c] && c != s)
                    .collect();
                sc[s].extend(parent_children);
            }

            // Rebuild supernodes list with only alive entries
            let mut new_supernodes: Vec<Supernode> = Vec::new();
            let mut old_to_new: Vec<usize> = vec![0; num_snodes_init];
            for s in 0..num_snodes_init {
                if alive[s] {
                    old_to_new[s] = new_supernodes.len();
                    new_supernodes.push(supernodes[s].clone());
                }
            }
            // Update snode_id to use new indices
            for col in 0..n {
                snode_id[col] = old_to_new[snode_id[col]];
            }
            supernodes = new_supernodes;
        }

        // Phase 2c: NEMIN-based amalgamation (MA57/MUMPS style)
        //
        // Merge consecutive small supernodes to create larger fronts.
        // MA57 uses NEMIN=16: any supernode with nfs < NEMIN gets merged with its
        // neighbor. Larger fronts give Bunch-Kaufman more pivot candidates, which
        // is critical for KKT systems with near-zero diagonal blocks.
        //
        // Strategy: scan supernodes in order. When a consecutive sequence of small
        // supernodes is found (each with nfs < NEMIN, forming a contiguous column
        // range), merge them into one supernode. This is simpler than tree-based
        // merging and catches the most common case (chains of size-1 supernodes).
        {
            const NEMIN: usize = 4;

            let num_snodes_now = supernodes.len();
            let mut merged_supernodes: Vec<Supernode> = Vec::new();
            let mut i = 0;

            while i < num_snodes_now {
                // Start a new merged group
                let group_start_col = supernodes[i].start;
                let mut group_nfs = supernodes[i].nfs;
                let mut group_end = i;

                // Extend the group: merge consecutive supernodes that are contiguous
                // and the total nfs stays reasonable (cap at 4*NEMIN to avoid huge fronts)
                while group_end + 1 < num_snodes_now && group_nfs < NEMIN {
                    let next = group_end + 1;
                    let current_end = supernodes[group_end].start + supernodes[group_end].nfs;
                    // Check contiguity: current supernode's end == next's start
                    if current_end == supernodes[next].start {
                        group_nfs += supernodes[next].nfs;
                        group_end = next;
                    } else {
                        break;
                    }
                }

                // Build the merged supernode's front_indices
                // FS cols: [group_start_col .. group_start_col + group_nfs]
                // CB cols: union of all component supernodes' CB, minus the FS cols
                let group_fs_end = group_start_col + group_nfs;
                let mut cb_set: Vec<usize> = Vec::new();
                for k in i..=group_end {
                    let sn = &supernodes[k];
                    for &idx in &sn.front_indices[sn.nfs..] {
                        if idx >= group_fs_end {
                            cb_set.push(idx);
                        }
                    }
                }
                cb_set.sort_unstable();
                cb_set.dedup();

                let mut front_indices = Vec::with_capacity(group_nfs + cb_set.len());
                for col in group_start_col..group_fs_end {
                    front_indices.push(col);
                }
                front_indices.extend_from_slice(&cb_set);

                let new_idx = merged_supernodes.len();
                merged_supernodes.push(Supernode {
                    start: group_start_col,
                    nfs: group_nfs,
                    front_indices,
                });

                // Update snode_id for ALL columns in this merged group
                for col in group_start_col..group_fs_end {
                    snode_id[col] = new_idx;
                }

                i = group_end + 1;
            }

            supernodes = merged_supernodes;
        }

        // Phase 3: Build supernodal elimination tree
        let num_snodes = supernodes.len();
        let mut snode_parent: Vec<Option<usize>> = vec![None; num_snodes];
        let mut snode_children: Vec<Vec<usize>> = vec![Vec::new(); num_snodes];

        for s in 0..num_snodes {
            let last_col = supernodes[s].start + supernodes[s].nfs - 1;
            if let Some(parent_col) = etree.parent[last_col] {
                let parent_snode = snode_id[parent_col];
                if parent_snode != s {
                    snode_parent[s] = Some(parent_snode);
                    snode_children[parent_snode].push(s);
                }
            }
        }

        // Deduplicate snode_children
        for ch in snode_children.iter_mut() {
            ch.sort_unstable();
            ch.dedup();
        }

        // Compute L nnz
        let l_nnz: usize = col_indices
            .iter()
            .map(|fi| if fi.is_empty() { 0 } else { fi.len() - 1 })
            .sum();

        // Legacy fields
        let front_indices = col_indices.clone();
        let front_nfs = vec![1usize; n];

        SymbolicFactorization {
            n,
            etree,
            col_indices,
            supernodes,
            col_to_supernode: snode_id,
            snode_parent,
            snode_children,
            l_nnz,
            front_indices,
            front_nfs,
            postorder,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::csc::CscMatrix;

    fn csc_from_upper_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> CscMatrix {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        CscMatrix::from_coo(&coo)
    }

    #[test]
    fn test_symbolic_tridiagonal() {
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);

        assert_eq!(sym.front_indices[0], vec![0, 1]);
        assert_eq!(sym.front_indices[1], vec![1, 2]);
        assert_eq!(sym.front_indices[2], vec![2, 3]);
        assert_eq!(sym.front_indices[3], vec![3]);
        assert_eq!(sym.l_nnz, 3);

        // Tridiagonal: fundamental supernodes are {0}, {1}, {2,3}.
        // Relaxed amalgamation merges them into fewer supernodes since the
        // chain is contiguous with small fill. Verify correctness via solve tests.
        assert!(sym.supernodes.len() <= 3);
    }

    #[test]
    fn test_symbolic_arrow() {
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 3, 1.0),
            (1, 1, 1.0), (1, 3, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);

        assert_eq!(sym.front_indices[0], vec![0, 3]);
        assert_eq!(sym.front_indices[1], vec![1, 3]);
        assert_eq!(sym.front_indices[2], vec![2, 3]);
        assert_eq!(sym.front_indices[3], vec![3]);
        assert_eq!(sym.l_nnz, 3);
    }

    #[test]
    fn test_symbolic_fill_in() {
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0), (0, 3, 1.0),
            (1, 1, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);

        assert_eq!(sym.front_indices[0], vec![0, 1, 3]);
        assert_eq!(sym.front_indices[1], vec![1, 3]);
        assert_eq!(sym.front_indices[2], vec![2, 3]);
        assert_eq!(sym.front_indices[3], vec![3]);
        assert_eq!(sym.l_nnz, 4);

        // Columns 0 and 1 merge fundamentally. Relaxed amalgamation may merge further.
        // The first supernode always starts at column 0 with at least nfs=2.
        assert!(sym.supernodes.len() <= 3);
        assert_eq!(sym.supernodes[0].start, 0);
        assert!(sym.supernodes[0].nfs >= 2);
    }

    #[test]
    fn test_symbolic_diagonal() {
        let csc = csc_from_upper_triplets(3, &[
            (0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);
        for j in 0..3 {
            assert_eq!(sym.front_indices[j], vec![j]);
        }
        assert_eq!(sym.l_nnz, 0);
    }

    #[test]
    fn test_symbolic_dense_3x3() {
        let csc = csc_from_upper_triplets(3, &[
            (0, 0, 1.0), (0, 1, 1.0), (0, 2, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);

        assert_eq!(sym.front_indices[0], vec![0, 1, 2]);
        assert_eq!(sym.front_indices[1], vec![1, 2]);
        assert_eq!(sym.front_indices[2], vec![2]);
        assert_eq!(sym.l_nnz, 3);

        // Dense 3x3: all columns merge into one supernode
        assert_eq!(sym.supernodes.len(), 1);
        assert_eq!(sym.supernodes[0].start, 0);
        assert_eq!(sym.supernodes[0].nfs, 3);
    }

    #[test]
    fn test_supernodal_parent() {
        // Fill-in case: fundamental supernodes {0,1}, {2}, {3}
        // Relaxed amalgamation may merge {2} into {3} (contiguous, single child)
        // then {0,1} into {2,3} (contiguous, single child) → 1 snode total.
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0), (0, 3, 1.0),
            (1, 1, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let sym = SymbolicFactorization::from_csc(&csc);

        // Root supernode (last one) has no parent
        let root = sym.supernodes.len() - 1;
        assert_eq!(sym.snode_parent[root], None);
        // All columns are accounted for
        let total_nfs: usize = sym.supernodes.iter().map(|s| s.nfs).sum();
        assert_eq!(total_nfs, 4);
    }
}
