/// Natural (identity) ordering — no permutation.
/// Returns (perm, perm_inv) both as identity.
pub fn natural_ordering(n: usize) -> (Vec<usize>, Vec<usize>) {
    let perm: Vec<usize> = (0..n).collect();
    let perm_inv = perm.clone();
    (perm, perm_inv)
}
