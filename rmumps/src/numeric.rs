use crate::csc::CscMatrix;
use crate::frontal::{FrontalMatrix, PartialFactorResult};
use crate::pivot::BunchKaufmanResult;
use crate::symbolic::SymbolicFactorization;
use crate::Inertia;
use rayon::prelude::*;
use std::cell::UnsafeCell;

/// A cell that is safe to share across threads when access is externally synchronized.
/// Level-set processing guarantees no concurrent access to the same slot.
struct SyncCell<T>(UnsafeCell<T>);
unsafe impl<T> Sync for SyncCell<T> {}
impl<T> SyncCell<T> {
    fn new(val: T) -> Self { Self(UnsafeCell::new(val)) }
    /// # Safety: caller must ensure no concurrent access to the same cell.
    unsafe fn get_mut(&self) -> &mut T { &mut *self.0.get() }
}

/// Result of the numeric multifrontal factorization.
#[derive(Debug)]
pub struct NumericFactorization {
    /// Per-supernode partial factorization results (indexed by supernode).
    pub node_factors: Vec<NodeFactor>,
    /// Overall inertia (sum across all supernodes).
    pub inertia: Inertia,
    /// Matrix dimension.
    pub n: usize,
    /// Number of supernodes.
    pub num_snodes: usize,
}

impl NumericFactorization {
    /// Return the minimum eigenvalue of D across all supernodes.
    ///
    /// For 1x1 blocks, this is the diagonal entry.
    /// For 2x2 blocks, this is the smaller eigenvalue.
    /// Used by IPM inertia correction to compute the minimum perturbation needed.
    pub fn min_diagonal(&self) -> Option<f64> {
        if self.n == 0 {
            return None;
        }
        let mut min_d = f64::INFINITY;
        for nf in &self.node_factors {
            let bk = &nf.bk;
            let n = bk.d_diag.len();
            let mut k = 0;
            while k < n {
                if k + 1 < n && bk.d_offdiag[k].abs() > 1e-12 {
                    // 2x2 block: compute smaller eigenvalue
                    let a = bk.d_diag[k];
                    let b = bk.d_offdiag[k];
                    let c = bk.d_diag[k + 1];
                    let trace = a + c;
                    let det = a * c - b * b;
                    let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
                    let eig_min = (trace - disc) / 2.0;
                    min_d = min_d.min(eig_min);
                    k += 2;
                } else {
                    min_d = min_d.min(bk.d_diag[k]);
                    k += 1;
                }
            }
        }
        Some(min_d)
    }
}

/// Factorization data for a single supernode.
#[derive(Debug)]
pub struct NodeFactor {
    /// Bunch-Kaufman factorization of the fully-summed block.
    pub bk: BunchKaufmanResult,
    /// L21 block (ncb x nfs), in BK-permuted column order.
    pub l21: crate::dense::DenseMat,
    /// Global indices of fully-summed variables.
    pub fs_indices: Vec<usize>,
    /// Global indices of contribution block variables.
    pub cb_indices: Vec<usize>,
}

/// Perform the numeric multifrontal factorization with parallel tree traversal.
///
/// Takes the CSC matrix (upper triangle, possibly permuted) and the symbolic factorization.
/// Uses supernodal factorization with level-set parallelism: supernodes at the same
/// tree depth are independent and factored concurrently via rayon.
pub fn multifrontal_factor(
    csc: &CscMatrix,
    sym: &SymbolicFactorization,
) -> NumericFactorization {
    let num_snodes = sym.supernodes.len();
    if num_snodes == 0 {
        return NumericFactorization {
            node_factors: vec![],
            inertia: Inertia { positive: 0, negative: 0, zero: 0 },
            n: sym.n,
            num_snodes: 0,
        };
    }

    // Compute topological levels: leaves = 0, level[s] = 1 + max(children levels)
    let mut level = vec![0usize; num_snodes];
    let mut max_level = 0usize;
    for s in 0..num_snodes {
        for &child in &sym.snode_children[s] {
            level[s] = level[s].max(level[child] + 1);
        }
        max_level = max_level.max(level[s]);
    }

    // Group supernodes by level (leaves first)
    let mut level_sets: Vec<Vec<usize>> = vec![Vec::new(); max_level + 1];
    for s in 0..num_snodes {
        level_sets[level[s]].push(s);
    }

    // Lock-free storage: level-set processing guarantees no concurrent access to same slot
    let node_factors: Vec<SyncCell<Option<NodeFactor>>> =
        (0..num_snodes).map(|_| SyncCell::new(None)).collect();
    let contributions: Vec<SyncCell<Option<(crate::dense::DenseMat, Vec<usize>)>>> =
        (0..num_snodes).map(|_| SyncCell::new(None)).collect();

    // Process levels bottom-up; within each level, supernodes are independent.
    // Use parallel execution only when there's enough work to justify the overhead.
    for level_nodes in &level_sets {
        // Estimate total work at this level
        let total_front_size: usize = level_nodes.iter()
            .map(|&s| {
                let f = sym.supernodes[s].front_indices.len();
                f * f  // approximate O(f²) work
            })
            .sum();

        // Use parallel only if enough work (> ~64KB of dense ops)
        if total_front_size > 4096 {
            level_nodes.par_iter().for_each(|&s| {
                // SAFETY: within a level, each thread writes to a unique s.
                // Children are at lower levels (already completed), so reads are safe.
                unsafe {
                    factor_supernode(s, csc, sym, &node_factors, &contributions);
                }
            });
        } else {
            for &s in level_nodes {
                unsafe {
                    factor_supernode(s, csc, sym, &node_factors, &contributions);
                }
            }
        }
    }

    // Compute inertia from completed factorizations (no locks needed)
    let mut total_inertia = Inertia { positive: 0, negative: 0, zero: 0 };
    let node_factors: Vec<NodeFactor> = node_factors
        .into_iter()
        .map(|cell| {
            let nf = cell.0.into_inner().unwrap();
            total_inertia.positive += nf.bk.inertia.positive;
            total_inertia.negative += nf.bk.inertia.negative;
            total_inertia.zero += nf.bk.inertia.zero;
            nf
        })
        .collect();

    NumericFactorization {
        node_factors,
        inertia: total_inertia,
        n: sym.n,
        num_snodes,
    }
}

/// Factor a single supernode: assemble entries, extend-add children, partial factor.
///
/// # Safety
/// Caller must ensure no concurrent access to `node_factors[s]` or `contributions[s]`.
/// Children's contributions must be fully written before this is called.
unsafe fn factor_supernode(
    s: usize,
    csc: &CscMatrix,
    sym: &SymbolicFactorization,
    node_factors: &[SyncCell<Option<NodeFactor>>],
    contributions: &[SyncCell<Option<(crate::dense::DenseMat, Vec<usize>)>>],
) {
    let snode = &sym.supernodes[s];
    let nfs = snode.nfs;
    let mut front = FrontalMatrix::new(snode.front_indices.clone(), nfs);

    // Assemble original matrix entries for all FS columns in this supernode.
    // Build a local index map for O(1) lookup instead of binary search.
    let fs_end = snode.start + nfs;

    // Build global-to-local index map for this front
    // For small fronts, binary search in assemble_entry is fine.
    // For larger fronts, a hash or direct-mapped array is faster.
    // Since indices are sorted and we know the range, we can use partition_point.

    for offset in 0..nfs {
        let col = snode.start + offset;
        let local_col = offset; // FS columns are first in front_indices

        for idx in csc.col_ptr[col]..csc.col_ptr[col + 1] {
            let row = csc.row_idx[idx];
            let val = csc.vals[idx];
            // Find local row index via binary search on front_indices
            if let Ok(local_row) = snode.front_indices.binary_search(&row) {
                let size = front.mat.nrows;
                front.mat.data[local_col * size + local_row] += val;
                if local_row != local_col {
                    front.mat.data[local_row * size + local_col] += val;
                }
            }
        }
    }

    // Off-diagonal entries: for columns gi > fs_end in the front
    for (fi, &gi) in snode.front_indices[nfs..].iter().enumerate() {
        let local_col = nfs + fi;
        let col_start = csc.col_ptr[gi];
        let col_end = csc.col_ptr[gi + 1];
        let rows = &csc.row_idx[col_start..col_end];
        let lo = rows.partition_point(|&r| r < snode.start);
        for k in lo..rows.len() {
            let row = rows[k];
            if row >= fs_end {
                break;
            }
            let local_row = row - snode.start; // FS rows have direct offset
            let size = front.mat.nrows;
            let val = csc.vals[col_start + k];
            front.mat.data[local_col * size + local_row] += val;
            front.mat.data[local_row * size + local_col] += val;
        }
    }

    // Extend-add contributions from children (already computed at a lower level)
    for &child_s in &sym.snode_children[s] {
        if let Some((contrib, contrib_indices)) = contributions[child_s].get_mut().take() {
            front.extend_add(&contrib, &contrib_indices);
        }
    }

    // Partial factorization
    let result = front.partial_factor();

    let PartialFactorResult { bk, l21, contrib, contrib_indices, fs_indices } = result;

    *node_factors[s].get_mut() = Some(NodeFactor {
        bk,
        l21,
        fs_indices,
        cb_indices: contrib_indices.clone(),
    });

    if !contrib_indices.is_empty() {
        *contributions[s].get_mut() = Some((contrib, contrib_indices));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::csc::CscMatrix;
    use crate::symbolic::SymbolicFactorization;

    fn factor_from_upper_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> (CscMatrix, NumericFactorization) {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);
        let num = multifrontal_factor(&csc, &sym);
        (csc, num)
    }

    #[test]
    fn test_diagonal_3x3() {
        let (_, num) = factor_from_upper_triplets(3, &[
            (0, 0, 2.0), (1, 1, 3.0), (2, 2, 5.0),
        ]);
        assert_eq!(num.inertia, Inertia { positive: 3, negative: 0, zero: 0 });
    }

    #[test]
    fn test_tridiagonal_spd() {
        // Tridiagonal SPD: diag=4, off-diag=1
        let (_, num) = factor_from_upper_triplets(4, &[
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0), (2, 3, 1.0),
            (3, 3, 4.0),
        ]);
        assert_eq!(num.inertia.positive, 4);
        assert_eq!(num.inertia.negative, 0);
    }

    #[test]
    fn test_indefinite_3x3() {
        // [[2, 0, 1], [0, 2, 1], [1, 1, 0]] — KKT-like, inertia (2, 1, 0)
        let (_, num) = factor_from_upper_triplets(3, &[
            (0, 0, 2.0), (0, 2, 1.0),
            (1, 1, 2.0), (1, 2, 1.0),
            (2, 2, 0.0),
        ]);
        assert_eq!(num.inertia.positive, 2);
        assert_eq!(num.inertia.negative, 1);
        assert_eq!(num.inertia.zero, 0);
    }

    #[test]
    fn test_5x5_kkt() {
        // 5x5 KKT: H=diag(4,5,6), A=[[1,0,1],[0,1,1]]
        // [[4,0,0,1,0],[0,5,0,0,1],[0,0,6,1,1],[1,0,1,0,0],[0,1,1,0,0]]
        let (_, num) = factor_from_upper_triplets(5, &[
            (0, 0, 4.0), (0, 3, 1.0),
            (1, 1, 5.0), (1, 4, 1.0),
            (2, 2, 6.0), (2, 3, 1.0), (2, 4, 1.0),
            (3, 3, 0.0),
            (4, 4, 0.0),
        ]);
        assert_eq!(num.inertia.positive, 3);
        assert_eq!(num.inertia.negative, 2);
    }

    #[test]
    fn test_dense_spd_3x3() {
        // Full 3x3 SPD: [[4, 2, 1], [2, 5, 3], [1, 3, 6]]
        let (_, num) = factor_from_upper_triplets(3, &[
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ]);
        assert_eq!(num.inertia.positive, 3);
        assert_eq!(num.inertia.negative, 0);
    }

    #[test]
    fn test_arrow_4x4() {
        // Arrow: diag=10, last col/row=1
        let (_, num) = factor_from_upper_triplets(4, &[
            (0, 0, 10.0), (0, 3, 1.0),
            (1, 1, 10.0), (1, 3, 1.0),
            (2, 2, 10.0), (2, 3, 1.0),
            (3, 3, 10.0),
        ]);
        assert_eq!(num.inertia.positive, 4);
        assert_eq!(num.inertia.negative, 0);
    }
}
