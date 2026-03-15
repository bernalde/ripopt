use crate::csc::CscMatrix;

/// Approximate Minimum Degree ordering using the SuiteSparse AMD algorithm.
///
/// Takes a symmetric CSC matrix (upper triangle) and returns a fill-reducing
/// permutation (perm, perm_inv) where perm[new] = old.
pub fn amd_ordering(csc: &CscMatrix) -> (Vec<usize>, Vec<usize>) {
    let n = csc.n;
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    if n == 1 {
        return (vec![0], vec![0]);
    }

    // Build full symmetric adjacency from upper triangle CSC.
    // The amd crate expects the complete adjacency structure (both (i,j) and (j,i)).
    // Count entries per column in the full matrix
    let mut full_col_count = vec![0usize; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            full_col_count[j] += 1;
            if i != j {
                full_col_count[i] += 1; // transpose entry
            }
        }
    }

    let mut full_col_ptr = vec![0i64; n + 1];
    for j in 0..n {
        full_col_ptr[j + 1] = full_col_ptr[j] + full_col_count[j] as i64;
    }
    let nnz = full_col_ptr[n] as usize;
    let mut full_row_idx = vec![0i64; nnz];
    let mut pos = vec![0usize; n];

    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            let p = full_col_ptr[j] as usize + pos[j];
            full_row_idx[p] = i as i64;
            pos[j] += 1;
            if i != j {
                let p = full_col_ptr[i] as usize + pos[i];
                full_row_idx[p] = j as i64;
                pos[i] += 1;
            }
        }
    }

    // Sort row indices within each column (required by amd)
    for j in 0..n {
        let start = full_col_ptr[j] as usize;
        let end = full_col_ptr[j + 1] as usize;
        full_row_idx[start..end].sort_unstable();
    }

    // Call SuiteSparse AMD
    let control = amd::Control::default();
    match amd::order::<i64>(n as i64, &full_col_ptr, &full_row_idx, &control) {
        Ok((perm_i64, perm_inv_i64, _info)) => {
            let perm: Vec<usize> = perm_i64.iter().map(|&x| x as usize).collect();
            let perm_inv: Vec<usize> = perm_inv_i64.iter().map(|&x| x as usize).collect();
            (perm, perm_inv)
        }
        Err(_) => {
            // Fallback: natural ordering
            let perm: Vec<usize> = (0..n).collect();
            let perm_inv: Vec<usize> = (0..n).collect();
            (perm, perm_inv)
        }
    }
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
        if perm.len() != n { return false; }
        let mut seen = vec![false; n];
        for &p in perm {
            if p >= n || seen[p] { return false; }
            seen[p] = true;
        }
        true
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
        let n = 50;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 10.0));
            if i < n - 1 {
                triplets.push((i, n - 1, 1.0));
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);
        let (amd_perm, _) = amd_ordering(&csc);
        assert!(is_valid_permutation(&amd_perm, n));
    }
}
