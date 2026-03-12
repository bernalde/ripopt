//! Incomplete LDL^T factorization for preconditioning.
//!
//! Builds an approximate multifrontal LDL^T factorization by applying dropping
//! rules to L21 and the Schur complement (contribution block) after each
//! supernode's partial factorization. This produces a sparser factorization
//! suitable for use as a preconditioner with iterative solvers like MINRES.

use crate::coo::CooMatrix;
use crate::csc::CscMatrix;
use crate::dense::DenseMat;
use crate::frontal::{FrontalMatrix, PartialFactorResult};
use crate::numeric::{NodeFactor, NumericFactorization};
use crate::ordering::{self, Ordering};
use crate::precond::Preconditioner;
use crate::solve::multifrontal_solve;
use crate::symbolic::SymbolicFactorization;
use crate::Inertia;

/// Options for the incomplete LDL^T factorization.
pub struct IncompleteLdltOptions {
    /// Drop tolerance: entries smaller than `drop_tolerance * column_norm` are zeroed.
    /// Default: 0.01.
    pub drop_tolerance: f64,
    /// Pivot threshold for Bunch-Kaufman pivoting. Default: 0.01.
    pub pivot_threshold: f64,
}

impl Default for IncompleteLdltOptions {
    fn default() -> Self {
        Self {
            drop_tolerance: 0.01,
            pivot_threshold: 0.01,
        }
    }
}

/// Incomplete LDL^T factorization that can be used as a preconditioner.
pub struct IncompleteLdlt {
    numeric: NumericFactorization,
    symbolic: SymbolicFactorization,
}

impl IncompleteLdlt {
    /// Build incomplete factorization from a permuted CSC matrix and symbolic factorization.
    pub fn new(
        csc: &CscMatrix,
        sym: &SymbolicFactorization,
        opts: &IncompleteLdltOptions,
    ) -> Self {
        let numeric = multifrontal_factor_incomplete(csc, sym, opts);
        Self {
            numeric,
            symbolic: sym.clone(),
        }
    }

    /// Return the inertia from the incomplete factorization.
    /// This is an approximation of the true inertia.
    pub fn inertia(&self) -> Inertia {
        self.numeric.inertia
    }

    /// Convenience: build from COO matrix with ordering.
    pub fn from_coo(coo: &CooMatrix, ordering: Ordering, opts: &IncompleteLdltOptions) -> Self {
        let csc = CscMatrix::from_coo(coo);
        let (perm, perm_inv) = ordering::compute_ordering(&csc, ordering);
        let permuted_csc = ordering::permute_symmetric_csc(&csc, &perm, &perm_inv);
        let sym = SymbolicFactorization::from_csc(&permuted_csc);
        Self::new(&permuted_csc, &sym, opts)
    }
}

impl Preconditioner for IncompleteLdlt {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        if multifrontal_solve(&self.numeric, &self.symbolic, r, z).is_err() {
            // Fall back to identity if solve fails
            z.copy_from_slice(r);
        }
    }
}

/// Sequential incomplete multifrontal factorization with dropping.
fn multifrontal_factor_incomplete(
    csc: &CscMatrix,
    sym: &SymbolicFactorization,
    opts: &IncompleteLdltOptions,
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

    let mut node_factors: Vec<Option<NodeFactor>> = (0..num_snodes).map(|_| None).collect();
    let mut contributions: Vec<Option<(DenseMat, Vec<usize>)>> =
        (0..num_snodes).map(|_| None).collect();

    // Process supernodes in postorder (leaves first)
    for s in 0..num_snodes {
        factor_supernode_incomplete(
            s, csc, sym, &mut node_factors, &mut contributions, opts,
        );
    }

    let mut total_inertia = Inertia { positive: 0, negative: 0, zero: 0 };
    let node_factors: Vec<NodeFactor> = node_factors
        .into_iter()
        .map(|nf| {
            let nf = nf.unwrap();
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

/// Factor a single supernode with dropping applied after partial factorization.
fn factor_supernode_incomplete(
    s: usize,
    csc: &CscMatrix,
    sym: &SymbolicFactorization,
    node_factors: &mut [Option<NodeFactor>],
    contributions: &mut [Option<(DenseMat, Vec<usize>)>],
    opts: &IncompleteLdltOptions,
) {
    let snode = &sym.supernodes[s];
    let nfs = snode.nfs;
    let mut front = FrontalMatrix::new(snode.front_indices.clone(), nfs);

    // Assemble original matrix entries (same as numeric.rs)
    let fs_end = snode.start + nfs;
    let size = front.mat.nrows;
    let front_indices = &snode.front_indices;

    // Build global-to-local index map
    let mut index_map = vec![usize::MAX; csc.n];
    for (local, &global) in front_indices.iter().enumerate() {
        index_map[global] = local;
    }

    for offset in 0..nfs {
        let col = snode.start + offset;
        let local_col = offset;

        for idx in csc.col_ptr[col]..csc.col_ptr[col + 1] {
            let row = csc.row_idx[idx];
            let val = csc.vals[idx];
            let local_row = index_map[row];
            if local_row != usize::MAX {
                front.mat.data[local_col * size + local_row] += val;
                if local_row != local_col {
                    front.mat.data[local_row * size + local_col] += val;
                }
            }
        }
    }

    // Off-diagonal entries from CB columns
    for (fi, &gi) in front_indices[nfs..].iter().enumerate() {
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
            let local_row = row - snode.start;
            let val = csc.vals[col_start + k];
            front.mat.data[local_col * size + local_row] += val;
            front.mat.data[local_row * size + local_col] += val;
        }
    }

    // Clean up index map
    for &global in front_indices.iter() {
        index_map[global] = usize::MAX;
    }

    // Extend-add contributions from children
    for &child_s in &sym.snode_children[s] {
        if let Some((contrib, contrib_indices)) = contributions[child_s].take() {
            front.extend_add(&contrib, &contrib_indices);
        }
    }

    // Partial factorization with threshold pivoting
    let mut result = if opts.pivot_threshold > 0.0 {
        front.partial_factor_threshold(opts.pivot_threshold)
    } else {
        front.partial_factor()
    };

    // Apply dropping to L21
    let drop_tol = opts.drop_tolerance;
    if drop_tol > 0.0 {
        let ncb = result.l21.nrows;
        let nfs_elim = result.l21.ncols;
        let l21_data = &mut result.l21.data;

        for j in 0..nfs_elim {
            let col_start = j * ncb;
            // Compute column norm (max absolute value)
            let col_norm: f64 = l21_data[col_start..col_start + ncb]
                .iter()
                .map(|v| v.abs())
                .fold(0.0_f64, f64::max);
            let threshold = drop_tol * col_norm;
            for i in 0..ncb {
                if l21_data[col_start + i].abs() < threshold {
                    l21_data[col_start + i] = 0.0;
                }
            }
        }

        // Apply dropping to contribution block (Schur complement)
        let contrib_data = &mut result.contrib.data;
        let ncb_contrib = result.contrib.nrows;
        if ncb_contrib > 0 {
            let ncb2 = ncb_contrib * ncb_contrib;
            let max_val: f64 = contrib_data[..ncb2]
                .iter()
                .map(|v| v.abs())
                .fold(0.0_f64, f64::max);
            let threshold = drop_tol * max_val;
            for v in &mut contrib_data[..ncb2] {
                if v.abs() < threshold {
                    *v = 0.0;
                }
            }
        }
    }

    let PartialFactorResult { bk, l21, contrib, contrib_indices, fs_indices, .. } = result;

    node_factors[s] = Some(NodeFactor {
        bk,
        l21,
        fs_indices,
        cb_indices: contrib_indices.clone(),
    });

    if !contrib_indices.is_empty() {
        contributions[s] = Some((contrib, contrib_indices));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::csc::CscMatrix;
    use crate::numeric::multifrontal_factor_threshold;
    use crate::solve::multifrontal_solve;
    use crate::symbolic::SymbolicFactorization;

    fn make_csc_and_sym(n: usize, triplets: &[(usize, usize, f64)]) -> (CscMatrix, SymbolicFactorization) {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);
        (csc, sym)
    }

    #[test]
    fn test_drop_tolerance_zero_matches_exact() {
        // With drop_tolerance=0, incomplete should match exact factorization
        let triplets = [
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ];
        let (csc, sym) = make_csc_and_sym(3, &triplets);

        // Exact factorization
        let exact = multifrontal_factor_threshold(&csc, &sym, 0.01);
        let mut x_exact = vec![0.0; 3];
        let b = [1.0, 2.0, 3.0];
        multifrontal_solve(&exact, &sym, &b, &mut x_exact).unwrap();

        // Incomplete with zero drop tolerance
        let opts = IncompleteLdltOptions { drop_tolerance: 0.0, pivot_threshold: 0.01 };
        let incomplete = IncompleteLdlt::new(&csc, &sym, &opts);
        let mut x_inc = vec![0.0; 3];
        multifrontal_solve(&incomplete.numeric, &incomplete.symbolic, &b, &mut x_inc).unwrap();

        for i in 0..3 {
            assert!((x_exact[i] - x_inc[i]).abs() < 1e-10,
                "x_exact[{}]={} != x_inc[{}]={}", i, x_exact[i], i, x_inc[i]);
        }
    }

    #[test]
    fn test_tridiagonal_exact_regardless() {
        // Tridiagonal has no fill, so incomplete should be exact for any drop tolerance
        let triplets = [
            (0, 0, 4.0), (0, 1, 1.0),
            (1, 1, 4.0), (1, 2, 1.0),
            (2, 2, 4.0), (2, 3, 1.0),
            (3, 3, 4.0),
        ];
        let (csc, sym) = make_csc_and_sym(4, &triplets);

        let exact = multifrontal_factor_threshold(&csc, &sym, 0.01);
        let mut x_exact = vec![0.0; 4];
        let b = [1.0, 2.0, 3.0, 4.0];
        multifrontal_solve(&exact, &sym, &b, &mut x_exact).unwrap();

        let opts = IncompleteLdltOptions { drop_tolerance: 0.5, pivot_threshold: 0.01 };
        let incomplete = IncompleteLdlt::new(&csc, &sym, &opts);
        let mut x_inc = vec![0.0; 4];
        multifrontal_solve(&incomplete.numeric, &incomplete.symbolic, &b, &mut x_inc).unwrap();

        for i in 0..4 {
            assert!((x_exact[i] - x_inc[i]).abs() < 1e-10,
                "mismatch at {}: {} vs {}", i, x_exact[i], x_inc[i]);
        }
    }

    #[test]
    fn test_preconditioner_reduces_residual() {
        let triplets = [
            (0, 0, 4.0), (0, 1, 2.0), (0, 2, 1.0),
            (1, 1, 5.0), (1, 2, 3.0),
            (2, 2, 6.0),
        ];
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(3, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);

        let opts = IncompleteLdltOptions { drop_tolerance: 0.0, pivot_threshold: 0.01 };
        let precond = IncompleteLdlt::new(&csc, &sym, &opts);

        let r = [1.0, 2.0, 3.0];
        let mut z = vec![0.0; 3];
        precond.apply(&r, &mut z);

        // Compute A*z
        let mut az = vec![0.0; 3];
        coo.matvec(&z, &mut az).unwrap();

        // Residual ‖r - A*M^{-1}*r‖ should be small (exact precond => 0)
        let resid: f64 = (0..3).map(|i| (r[i] - az[i]).powi(2)).sum::<f64>().sqrt();
        let r_norm: f64 = (0..3).map(|i| r[i].powi(2)).sum::<f64>().sqrt();
        assert!(resid < 1e-8 * r_norm, "residual {} too large (r_norm={})", resid, r_norm);
    }

    #[test]
    fn test_spd_grid_with_minres() {
        // Build a 1D Laplacian (tridiagonal: diag=2, offdiag=-1), size 10
        let n = 10;
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n {
            rows.push(i); cols.push(i); vals.push(2.0);
            if i + 1 < n {
                rows.push(i); cols.push(i + 1); vals.push(-1.0);
            }
        }
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        let csc = CscMatrix::from_coo(&coo);
        let sym = SymbolicFactorization::from_csc(&csc);

        let opts = IncompleteLdltOptions { drop_tolerance: 0.01, pivot_threshold: 0.01 };
        let precond = IncompleteLdlt::new(&csc, &sym, &opts);

        // RHS: b = [1, 0, 0, ..., 0]
        let mut b = vec![0.0; n];
        b[0] = 1.0;

        let mut x = vec![0.0; n];
        let minres_opts = crate::minres::MinresOptions { max_iter: 100, tol: 1e-10 };
        let mv = |xv: &[f64], yv: &mut [f64]| { coo.matvec(xv, yv).unwrap(); };
        let res = crate::minres::minres(n, mv, Some(&precond), &b, &mut x, &minres_opts);

        assert!(res.converged, "MINRES did not converge: iters={}, resid={}", res.iterations, res.residual_norm);

        // Verify solution
        let mut ax = vec![0.0; n];
        coo.matvec(&x, &mut ax).unwrap();
        let resid: f64 = (0..n).map(|i| (ax[i] - b[i]).powi(2)).sum::<f64>().sqrt();
        assert!(resid < 1e-8, "solution residual {} too large", resid);
    }

    #[test]
    fn test_full_pipeline_coo_to_minres() {
        // Full pipeline: COO → incomplete LDL^T → preconditioned MINRES
        // Use a denser SPD matrix: 5x5 with some fill
        let triplets = [
            (0, 0, 10.0), (0, 1, 1.0), (0, 4, 2.0),
            (1, 1, 10.0), (1, 2, 1.0),
            (2, 2, 10.0), (2, 3, 1.0),
            (3, 3, 10.0), (3, 4, 1.0),
            (4, 4, 10.0),
        ];
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(5, rows, cols, vals).unwrap();

        let opts = IncompleteLdltOptions { drop_tolerance: 0.01, pivot_threshold: 0.01 };
        let precond = IncompleteLdlt::from_coo(&coo, Ordering::Natural, &opts);

        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = vec![0.0; 5];
        let minres_opts = crate::minres::MinresOptions { max_iter: 100, tol: 1e-8 };
        let mv = |xv: &[f64], yv: &mut [f64]| { coo.matvec(xv, yv).unwrap(); };
        let res = crate::minres::minres(5, mv, Some(&precond), &b, &mut x, &minres_opts);

        assert!(res.converged, "MINRES did not converge: iters={}, resid={}", res.iterations, res.residual_norm);

        let mut ax = vec![0.0; 5];
        coo.matvec(&x, &mut ax).unwrap();
        let resid: f64 = (0..5).map(|i| (ax[i] - b[i]).powi(2)).sum::<f64>().sqrt();
        assert!(resid < 1e-8, "solution residual {} too large", resid);
    }
}
