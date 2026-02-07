use crate::linear_solver::{LinearSolver, SymmetricMatrix};

/// Information about the KKT system structure.
pub struct KktSystem {
    /// Dimension of the full KKT matrix (n + m).
    pub dim: usize,
    /// Number of primal variables.
    pub n: usize,
    /// Number of constraints.
    pub m: usize,
    /// The assembled KKT matrix.
    pub matrix: SymmetricMatrix,
    /// Right-hand side vector.
    pub rhs: Vec<f64>,
}

/// Assemble the augmented KKT matrix:
/// ```text
/// [H + Sigma + delta_w*I,  J^T ]   [dx]   [r_d]
/// [J,                -delta_c*I ] * [dy] = [r_p]
/// ```
///
/// Where:
/// - H is the Hessian of the Lagrangian (n x n, lower triangle in COO)
/// - J is the constraint Jacobian (m x n, in COO)
/// - Sigma contains the barrier diagonal terms: Sigma_ii = z_l_i/(x_i - x_l_i) + z_u_i/(x_u_i - x_i)
/// - r_d is the dual residual
/// - r_p is the primal residual
///
/// # Arguments
/// - `n`: number of variables
/// - `m`: number of constraints
/// - `hess_rows`, `hess_cols`, `hess_vals`: Hessian lower triangle in COO
/// - `jac_rows`, `jac_cols`, `jac_vals`: Jacobian in COO
/// - `sigma`: barrier diagonal (length n)
/// - `grad_f`: gradient of objective (length n)
/// - `g`: constraint values (length m)
/// - `g_l`, `g_u`: constraint bounds
/// - `y`: current constraint multipliers (length m)
/// - `z_l`, `z_u`: bound multipliers
/// - `x`, `x_l`, `x_u`: current point and bounds
/// - `mu`: barrier parameter
#[allow(clippy::too_many_arguments)]
pub fn assemble_kkt(
    n: usize,
    m: usize,
    hess_rows: &[usize],
    hess_cols: &[usize],
    hess_vals: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    sigma: &[f64],
    grad_f: &[f64],
    g: &[f64],
    g_l: &[f64],
    g_u: &[f64],
    y: &[f64],
    z_l: &[f64],
    z_u: &[f64],
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    mu: f64,
) -> KktSystem {
    let dim = n + m;
    let mut matrix = SymmetricMatrix::zeros(dim);
    let mut rhs = vec![0.0; dim];

    // (1,1) block: H + Sigma
    for (idx, (&row, &col)) in hess_rows.iter().zip(hess_cols.iter()).enumerate() {
        matrix.add(row, col, hess_vals[idx]);
    }

    // Add barrier diagonal Sigma for variable bounds
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        matrix.add(i, i, sigma[i]);
    }

    // (2,1) block: J
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        matrix.add(n + row, col, jac_vals[idx]);
    }

    // RHS: dual residual r_d (first n entries)
    // Ipopt convention (L = f + y^T g): stationarity is ∇f + J^T y - z_l + z_u = 0
    // Newton RHS: r_d = -(∇f + J^T y - z_l + z_u) + barrier correction
    for i in 0..n {
        let mut rd = -grad_f[i];
        rd += z_l[i];
        rd -= z_u[i];

        if x_l[i].is_finite() {
            rd += mu / (x[i] - x_l[i]);
        }
        if x_u[i].is_finite() {
            rd -= mu / (x_u[i] - x[i]);
        }

        rhs[i] = rd;
    }

    // Subtract J^T * y contribution from r_d
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        rhs[col] -= jac_vals[idx] * y[row];
    }

    // RHS: primal residual r_p (last m entries)
    // For equality constraints: rhs = -(g - g_l)
    // For inequality constraints with slack elimination:
    //   When feasible (slack > 0): rhs = -(g - g_l) - mu/y (barrier correction)
    //   When infeasible (slack <= 0): rhs = -(g - g_l) (no barrier, drive to feasibility)
    for i in 0..m {
        let is_equality = g_l[i].is_finite() && g_u[i].is_finite() && (g_l[i] - g_u[i]).abs() < 1e-15;
        if is_equality {
            rhs[n + i] = -(g[i] - g_l[i]);
        } else if g_l[i].is_finite() {
            rhs[n + i] = -(g[i] - g_l[i]);
        } else if g_u[i].is_finite() {
            rhs[n + i] = -(g[i] - g_u[i]);
        }
    }

    // (2,2) block: slack contribution for inequality constraints.
    // After implicit slack elimination, D_i = -y_i / s_i where s_i is the
    // constraint slack. The matrix entry is -D_i.
    // When the constraint is violated (s <= 0), skip the barrier contribution
    // and let the Newton step drive toward feasibility.
    for i in 0..m {
        let is_equality = g_l[i].is_finite() && g_u[i].is_finite() && (g_l[i] - g_u[i]).abs() < 1e-15;
        if is_equality {
            continue;
        }

        let mut d_ii = 0.0;
        if g_l[i].is_finite() {
            let slack = g[i] - g_l[i];
            if slack > 1e-10 && y[i] < -1e-20 {
                // Feasible: D = -y/s > 0 (y < 0, s > 0)
                d_ii += -y[i] / slack;
            }
        }
        if g_u[i].is_finite() {
            let slack = g_u[i] - g[i];
            if slack > 1e-10 && y[i] > 1e-20 {
                // Feasible: D = y/s > 0 (y > 0, s > 0)
                d_ii += y[i] / slack;
            }
        }
        if d_ii > 0.0 {
            // Matrix entry is -D (negative for correct KKT inertia)
            matrix.add(n + i, n + i, -d_ii);
        }
    }

    KktSystem {
        dim,
        n,
        m,
        matrix,
        rhs,
    }
}

/// Compute the barrier diagonal Sigma.
///
/// Sigma_ii = z_l_i / (x_i - x_l_i) + z_u_i / (x_u_i - x_i)
pub fn compute_sigma(
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    z_l: &[f64],
    z_u: &[f64],
) -> Vec<f64> {
    let n = x.len();
    let mut sigma = vec![0.0; n];
    for i in 0..n {
        if x_l[i].is_finite() {
            let slack = (x[i] - x_l[i]).max(1e-20);
            sigma[i] += z_l[i] / slack;
        }
        if x_u[i].is_finite() {
            let slack = (x_u[i] - x[i]).max(1e-20);
            sigma[i] += z_u[i] / slack;
        }
    }
    sigma
}

/// Parameters for inertia correction.
pub struct InertiaCorrectionParams {
    /// Initial primal regularization.
    pub delta_w_init: f64,
    /// Initial constraint regularization.
    pub delta_c: f64,
    /// Growth factor for delta_w.
    pub delta_w_growth: f64,
    /// Maximum number of correction attempts.
    pub max_attempts: usize,
    /// Last successful delta_w (for warm-starting perturbation).
    pub delta_w_last: f64,
}

impl Default for InertiaCorrectionParams {
    fn default() -> Self {
        Self {
            delta_w_init: 1e-4,
            delta_c: 1e-8,
            delta_w_growth: 8.0,
            max_attempts: 10,
            delta_w_last: 0.0,
        }
    }
}

/// Perform KKT factorization with inertia correction.
///
/// Factor the KKT matrix and check inertia. If inertia is wrong
/// (should be (n, m, 0) for an n-variable, m-constraint problem),
/// add regularization and re-factor.
///
/// Returns the factored solver and the regularization used.
pub fn factor_with_inertia_correction(
    kkt: &mut KktSystem,
    solver: &mut dyn LinearSolver,
    params: &mut InertiaCorrectionParams,
) -> Result<(f64, f64), crate::linear_solver::SolverError> {
    let n = kkt.n;
    let m = kkt.m;

    // First attempt: factor without perturbation
    let inertia = solver.factor(&kkt.matrix)?;

    if let Some(inertia) = inertia {
        if inertia.positive == n && inertia.negative == m && inertia.zero == 0 {
            params.delta_w_last = 0.0;
            return Ok((0.0, 0.0));
        }
    }

    // Inertia is wrong — apply perturbation and re-factor
    let mut delta_w = if params.delta_w_last == 0.0 {
        params.delta_w_init
    } else {
        (params.delta_w_last / params.delta_w_growth).max(params.delta_w_init)
    };
    let delta_c = params.delta_c;

    for attempt in 0..params.max_attempts {
        // Create perturbed matrix
        let mut perturbed = kkt.matrix.clone();
        perturbed.add_diagonal_range(0, n, delta_w);
        if m > 0 {
            perturbed.add_diagonal_range(n, n + m, -delta_c);
        }

        let inertia = solver.factor(&perturbed)?;

        if let Some(inertia) = inertia {
            if inertia.positive == n && inertia.negative == m && inertia.zero == 0 {
                // Update the KKT matrix to the perturbed version
                kkt.matrix = perturbed;
                params.delta_w_last = delta_w;
                return Ok((delta_w, delta_c));
            }
        }

        // Increase perturbation
        delta_w *= params.delta_w_growth;

        log::debug!(
            "Inertia correction attempt {}: delta_w = {:.2e}, inertia = {:?}",
            attempt + 1,
            delta_w,
            inertia
        );
    }

    Err(crate::linear_solver::SolverError::NumericalFailure(
        "inertia correction failed after maximum attempts".to_string(),
    ))
}

/// Solve the KKT system for the search direction, given a factored solver.
///
/// Returns (dx, dy) where dx is the primal step and dy is the dual step.
/// Bound multiplier steps dz_l, dz_u are recovered from complementarity.
pub fn solve_for_direction(
    kkt: &KktSystem,
    solver: &mut dyn LinearSolver,
) -> Result<(Vec<f64>, Vec<f64>), crate::linear_solver::SolverError> {
    let mut solution = vec![0.0; kkt.dim];
    solver.solve(&kkt.rhs, &mut solution)?;

    let dx = solution[..kkt.n].to_vec();
    let dy = solution[kkt.n..].to_vec();

    Ok((dx, dy))
}

/// Recover bound multiplier steps from complementarity.
///
/// dz_l_i = (mu - z_l_i * (x_i - x_l_i) - z_l_i * dx_i) / (x_i - x_l_i)
///        = (mu / (x_i - x_l_i)) - z_l_i - (z_l_i / (x_i - x_l_i)) * dx_i
///        = sigma_l_i * dx_i ... (simplified from complementarity)
///
/// More precisely:
/// dz_l_i = (mu - z_l_i * s_l_i) / s_l_i - z_l_i * dx_i / s_l_i
///        where s_l_i = x_i - x_l_i
pub fn recover_dz(
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    z_l: &[f64],
    z_u: &[f64],
    dx: &[f64],
    mu: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = x.len();
    let mut dz_l = vec![0.0; n];
    let mut dz_u = vec![0.0; n];

    for i in 0..n {
        if x_l[i].is_finite() {
            let s_l = (x[i] - x_l[i]).max(1e-20);
            dz_l[i] = (mu - z_l[i] * s_l) / s_l - (z_l[i] / s_l) * dx[i];
        }
        if x_u[i].is_finite() {
            let s_u = (x_u[i] - x[i]).max(1e-20);
            dz_u[i] = (mu - z_u[i] * s_u) / s_u + (z_u[i] / s_u) * dx[i];
        }
    }

    (dz_l, dz_u)
}
