//! Preconditioner trait and basic implementations for iterative solvers.

use crate::csc::CscMatrix;

/// Trait for preconditioners: compute z = M^{-1} * r.
pub trait Preconditioner {
    fn apply(&self, r: &[f64], z: &mut [f64]);
}

/// Identity preconditioner (no preconditioning).
pub struct IdentityPrecond;

impl Preconditioner for IdentityPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        z.copy_from_slice(r);
    }
}

/// Diagonal (Jacobi) preconditioner: z_i = r_i / a_ii.
pub struct DiagonalPrecond {
    pub inv_diag: Vec<f64>,
}

impl DiagonalPrecond {
    /// Build from a CSC matrix (upper triangle). Extracts diagonal and inverts.
    /// Zero diagonals are replaced with 1.0 (no scaling for that row).
    pub fn from_csc(csc: &CscMatrix) -> Self {
        let mut inv_diag = vec![1.0; csc.n];
        for j in 0..csc.n {
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                if csc.row_idx[idx] == j {
                    let d = csc.vals[idx];
                    if d.abs() > 1e-30 {
                        inv_diag[j] = 1.0 / d;
                    }
                    break;
                }
            }
        }
        Self { inv_diag }
    }
}

impl Preconditioner for DiagonalPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for i in 0..r.len() {
            z[i] = r[i] * self.inv_diag[i];
        }
    }
}
