//! Preconditioned MINRES iterative solver for symmetric (possibly indefinite) systems.
//!
//! Implements the Paige & Saunders (1975) MINRES algorithm following the
//! Choi-Saunders reference implementation, with optional SPD preconditioning.

use crate::precond::Preconditioner;

/// Options for the MINRES solver.
pub struct MinresOptions {
    /// Maximum number of iterations (default: 1000).
    pub max_iter: usize,
    /// Relative residual tolerance (default: 1e-8).
    pub tol: f64,
}

impl Default for MinresOptions {
    fn default() -> Self {
        Self {
            max_iter: 1000,
            tol: 1e-8,
        }
    }
}

/// Result of a MINRES solve.
pub struct MinresResult {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual norm.
    pub residual_norm: f64,
    /// Whether the solver converged within tolerance.
    pub converged: bool,
}

/// Solve Ax = b using preconditioned MINRES.
///
/// `matvec` computes y = A * x. `precond` is an optional SPD preconditioner M
/// such that z = M^{-1} * r is computed. If `precond` is None, identity is used.
///
/// The initial guess x is used as x₀; on return x holds the solution.
pub fn minres<F>(
    n: usize,
    matvec: F,
    precond: Option<&dyn Preconditioner>,
    b: &[f64],
    x: &mut [f64],
    opts: &MinresOptions,
) -> MinresResult
where
    F: Fn(&[f64], &mut [f64]),
{
    let bnorm = norm(b);
    if bnorm == 0.0 {
        x.iter_mut().for_each(|v| *v = 0.0);
        return MinresResult { iterations: 0, residual_norm: 0.0, converged: true };
    }

    // r0 = b - A*x0
    let mut r2 = vec![0.0; n]; // tracks unnormalized Lanczos vectors
    if x.iter().all(|&v| v == 0.0) {
        r2.copy_from_slice(b);
    } else {
        let mut tmp = vec![0.0; n];
        matvec(x, &mut tmp);
        for i in 0..n { r2[i] = b[i] - tmp[i]; }
    }

    // Apply preconditioner
    let mut y = vec![0.0; n];
    apply_precond(precond, &r2, &mut y);

    let rtz = dot(&r2, &y);
    if rtz < 0.0 {
        // Preconditioner not SPD; fall back to unpreconditioned
        return minres(n, matvec, None, b, x, opts);
    }

    let beta1 = rtz.sqrt();
    if beta1 / bnorm < opts.tol {
        return MinresResult { iterations: 0, residual_norm: beta1, converged: true };
    }

    // Following Choi-Saunders reference implementation
    let mut r1 = vec![0.0; n];
    let mut v = vec![0.0; n];

    let mut oldb: f64 = 0.0;
    let mut beta: f64 = beta1;
    let mut dbar: f64 = 0.0;
    let mut epsln: f64 = 0.0;
    let mut phibar: f64 = beta1;

    // Givens rotation: initial cs=-1, sn=0 (Saunders convention)
    let mut cs: f64 = -1.0;
    let mut sn: f64 = 0.0;

    // Solution update vectors
    let mut w = vec![0.0; n];
    let mut w1 = vec![0.0; n];
    let mut w2 = vec![0.0; n];

    for itn in 1..=opts.max_iter {
        // Normalized Lanczos vector
        let s = 1.0 / beta;
        for i in 0..n { v[i] = s * y[i]; }

        // y = A * v
        matvec(&v, &mut y);

        // Lanczos 3-term recurrence (Saunders variant)
        if itn >= 2 {
            let coeff = beta / oldb;
            for i in 0..n { y[i] -= coeff * r1[i]; }
        }

        let alfa = dot(&v, &y);
        let coeff = alfa / beta;
        for i in 0..n { y[i] -= coeff * r2[i]; }

        // Shift r1, r2
        std::mem::swap(&mut r1, &mut r2);
        r2.copy_from_slice(&y);

        // Apply preconditioner
        apply_precond(precond, &r2, &mut y);

        let rtz = dot(&r2, &y);
        if rtz < 0.0 {
            return minres(n, matvec, None, b, x, opts);
        }

        oldb = beta;
        beta = rtz.sqrt();

        // QR factorization update
        // Apply previous rotation [cs sn; sn -cs] to [dbar; alfa]
        let oldeps = epsln;
        let delta = cs * dbar + sn * alfa;
        let gbar = sn * dbar - cs * alfa;
        epsln = sn * beta;
        dbar = -cs * beta;

        // New rotation to eliminate beta against gbar
        let gamma = (gbar * gbar + beta * beta).sqrt();
        let gamma = if gamma < 1e-30 { 1e-30 } else { gamma };
        cs = gbar / gamma;
        sn = beta / gamma;
        let phi = cs * phibar;
        phibar *= sn;

        // Update w vectors: w_new = (v - oldeps*w1 - delta*w2) / gamma
        let denom = 1.0 / gamma;
        std::mem::swap(&mut w1, &mut w2);
        std::mem::swap(&mut w2, &mut w);
        for i in 0..n {
            w[i] = (v[i] - oldeps * w1[i] - delta * w2[i]) * denom;
        }

        // Update solution
        for i in 0..n {
            x[i] += phi * w[i];
        }

        // Check convergence using residual norm estimate
        let rnorm = phibar.abs();
        if rnorm / bnorm < opts.tol {
            return MinresResult {
                iterations: itn,
                residual_norm: rnorm,
                converged: true,
            };
        }
    }

    MinresResult {
        iterations: opts.max_iter,
        residual_norm: phibar.abs(),
        converged: false,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

fn apply_precond(precond: Option<&dyn Preconditioner>, r: &[f64], z: &mut [f64]) {
    match precond {
        Some(p) => p.apply(r, z),
        None => z.copy_from_slice(r),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::precond::{DiagonalPrecond, IdentityPrecond};

    /// Symmetric matvec from dense upper triangle stored row-major.
    fn dense_sym_matvec(n: usize, upper: &[f64]) -> impl Fn(&[f64], &mut [f64]) + '_ {
        move |x: &[f64], y: &mut [f64]| {
            y.iter_mut().for_each(|v| *v = 0.0);
            for i in 0..n {
                for j in i..n {
                    let a = upper[i * n + j];
                    y[i] += a * x[j];
                    if i != j {
                        y[j] += a * x[i];
                    }
                }
            }
        }
    }

    #[test]
    fn test_spd_3x3() {
        // A = [[4, 1, 0], [1, 3, 1], [0, 1, 4]]
        let a = [
            4.0, 1.0, 0.0,
            0.0, 3.0, 1.0,
            0.0, 0.0, 4.0,
        ];
        let mv = dense_sym_matvec(3, &a);
        let b = [1.0, 2.0, 3.0];
        let mut x = vec![0.0; 3];
        let res = minres(3, mv, None, &b, &mut x, &MinresOptions::default());
        assert!(res.converged, "did not converge in {} iters, resid={}", res.iterations, res.residual_norm);
        assert!(res.iterations <= 3, "took {} iterations for 3x3", res.iterations);

        // Check Ax ≈ b
        let mv2 = dense_sym_matvec(3, &a);
        let mut ax = vec![0.0; 3];
        mv2(&x, &mut ax);
        for i in 0..3 {
            assert!((ax[i] - b[i]).abs() < 1e-6, "residual[{}] = {}", i, (ax[i] - b[i]).abs());
        }
    }

    #[test]
    fn test_indefinite_kkt_5x5() {
        // KKT: H=diag(4,5,6), A=[[1,0,1],[0,1,1]]
        let a = [
            4.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 5.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 6.0, 1.0, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mv = dense_sym_matvec(5, &a);
        let b = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = vec![0.0; 5];
        let opts = MinresOptions { max_iter: 100, tol: 1e-8 };
        let res = minres(5, mv, None, &b, &mut x, &opts);
        assert!(res.converged, "did not converge in {} iters, resid={}", res.iterations, res.residual_norm);

        // Check residual
        let mv2 = dense_sym_matvec(5, &a);
        let mut ax = vec![0.0; 5];
        mv2(&x, &mut ax);
        let resid: f64 = (0..5).map(|i| (ax[i] - b[i]).powi(2)).sum::<f64>().sqrt();
        assert!(resid < 1e-6, "residual = {}", resid);
    }

    #[test]
    fn test_identity_precond_same_as_unpreconditioned() {
        let a = [
            4.0, 1.0, 0.0,
            0.0, 3.0, 1.0,
            0.0, 0.0, 4.0,
        ];
        let b = [1.0, 2.0, 3.0];

        let mut x1 = vec![0.0; 3];
        let mv1 = dense_sym_matvec(3, &a);
        let r1 = minres(3, mv1, None, &b, &mut x1, &MinresOptions::default());

        let mut x2 = vec![0.0; 3];
        let mv2 = dense_sym_matvec(3, &a);
        let id = IdentityPrecond;
        let r2 = minres(3, mv2, Some(&id), &b, &mut x2, &MinresOptions::default());

        assert_eq!(r1.iterations, r2.iterations);
    }

    #[test]
    fn test_diagonal_precond() {
        // Ill-conditioned SPD: diag(100, 1, 0.01) + off-diag
        let a = [
            100.0, 1.0, 0.0,
            0.0,   1.0, 0.1,
            0.0,   0.0, 0.01,
        ];
        let b = [101.0, 1.1, 0.11];
        let opts = MinresOptions { max_iter: 100, tol: 1e-10 };

        let mut x_no = vec![0.0; 3];
        let mv1 = dense_sym_matvec(3, &a);
        let r_no = minres(3, mv1, None, &b, &mut x_no, &opts);

        // Build diagonal preconditioner manually
        let diag_precond = DiagonalPrecond {
            inv_diag: vec![1.0 / 100.0, 1.0 / 1.0, 1.0 / 0.01],
        };
        let mut x_diag = vec![0.0; 3];
        let mv2 = dense_sym_matvec(3, &a);
        let r_diag = minres(3, mv2, Some(&diag_precond), &b, &mut x_diag, &opts);

        assert!(r_diag.converged);
        assert!(r_no.converged);
        // With diagonal preconditioning, should converge in fewer or equal iterations
        assert!(r_diag.iterations <= r_no.iterations,
            "preconditioned {} > unpreconditioned {}", r_diag.iterations, r_no.iterations);
    }

    #[test]
    fn test_zero_rhs() {
        let a = [4.0, 1.0, 0.0, 3.0];
        let mv = dense_sym_matvec(2, &a);
        let b = [0.0, 0.0];
        let mut x = vec![0.0; 2];
        let res = minres(2, mv, None, &b, &mut x, &MinresOptions::default());
        assert!(res.converged);
        assert_eq!(res.iterations, 0);
        assert_eq!(x, vec![0.0, 0.0]);
    }

    #[test]
    fn test_already_solved() {
        let a = [
            4.0, 1.0,
            0.0, 3.0,
        ];
        let mv = dense_sym_matvec(2, &a);
        // x = [1, 1] => Ax = [5, 4]
        let b = [5.0, 4.0];
        let mut x = vec![1.0, 1.0];
        let res = minres(2, mv, None, &b, &mut x, &MinresOptions::default());
        assert!(res.converged);
        assert_eq!(res.iterations, 0);
    }
}
