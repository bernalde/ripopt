use crate::convergence::is_equality_constraint;
use crate::linear_solver::{KktMatrix, LinearSolver, SolverError, SymmetricMatrix};

/// Information about the KKT system structure.
pub struct KktSystem {
    /// Dimension of the full KKT matrix (n + m).
    pub dim: usize,
    /// Number of primal variables.
    pub n: usize,
    /// Number of constraints.
    pub m: usize,
    /// The assembled KKT matrix (dense or sparse).
    pub matrix: KktMatrix,
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
    use_sparse: bool,
    _v_l: &[f64],
    _v_u: &[f64],
) -> KktSystem {
    let dim = n + m;
    let capacity = hess_rows.len() + jac_rows.len() + n + m;
    let mut matrix = if use_sparse {
        KktMatrix::zeros_sparse(dim, capacity)
    } else {
        KktMatrix::zeros_dense(dim)
    };
    let mut rhs = vec![0.0; dim];

    // (1,1) block: H + Sigma
    for (idx, (&row, &col)) in hess_rows.iter().zip(hess_cols.iter()).enumerate() {
        let v = hess_vals[idx];
        if v.is_nan() || v.is_infinite() {
            log::warn!("NaN/Inf in Hessian at ({}, {}): {}", row, col, v);
        }
        matrix.add(row, col, v);
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
    //
    // After eliminating dz from the full Newton system, the z_l and z_u terms
    // cancel algebraically: correct condensed RHS = -∇f - J^T*y + μ/s_l - μ/s_u.
    // However, keeping the z terms (r_d = -∇f + z_l - z_u + μ/s_l - μ/s_u - J^T*y)
    // provides better convergence in practice by tracking the dual residual,
    // especially when z deviates from μ/s due to safeguarding (kappa_sigma).
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

    // RHS: primal residual r_p (last m entries) and (2,2) block for inequality constraints.
    //
    // After condensing the slack variables from the KKT system, the condensed system is:
    //   [H + Σ_x    J^T         ] [Δx]   [r_d                            ]
    //   [J          -Σ_s^{-1}   ] [Δy] = [Σ_s^{-1} * (y + μ/s_l - μ/s_u)]
    //
    // where Σ_s = z_sl/s_l + z_su/s_u is the barrier contribution from constraint slacks,
    // z_sl, z_su are the slack bound multipliers, and s_l = g - g_l, s_u = g_u - g.
    //
    // For equality constraints: no slack, (2,2) = 0, r_c = -(g - g_l).
    // For infeasible inequality constraints: no barrier, r_c = -(g - bound).
    for i in 0..m {
        if is_equality_constraint(g_l[i], g_u[i]) {
            rhs[n + i] = -(g[i] - g_l[i]);
            continue;
        }

        // Compute Σ_s and the RHS correction term (y + μ/s_l - μ/s_u)
        let mut sigma_s = 0.0;
        let mut rhs_correction = y[i]; // starts with y
        let mut any_feasible = false;
        let mut rhs_infeasible = 0.0;

        if g_l[i].is_finite() {
            let slack = g[i] - g_l[i];
            if slack >= -1e-8 {
                // Feasible or at bound: use barrier with safeguarded slack
                let safe_slack = slack.max(mu.max(1e-10));
                // Heuristic: use |y| when y has correct sign, else barrier estimate mu/s
                let z_sl = if y[i] < -1e-20 {
                    -y[i]
                } else {
                    mu / safe_slack
                };
                sigma_s += z_sl / safe_slack;
                rhs_correction += mu / safe_slack;
                any_feasible = true;
            } else {
                // Truly infeasible: drive toward feasibility
                rhs_infeasible += -(g[i] - g_l[i]);
            }
        }
        if g_u[i].is_finite() {
            let slack = g_u[i] - g[i];
            if slack >= -1e-8 {
                // Feasible or at bound: use barrier with safeguarded slack
                let safe_slack = slack.max(mu.max(1e-10));
                // Heuristic: use |y| when y has correct sign, else barrier estimate mu/s
                let z_su = if y[i] > 1e-20 {
                    y[i]
                } else {
                    mu / safe_slack
                };
                sigma_s += z_su / safe_slack;
                rhs_correction -= mu / safe_slack;
                any_feasible = true;
            } else {
                // Truly infeasible: drive toward feasibility
                rhs_infeasible += -(g[i] - g_u[i]);
            }
        }

        if any_feasible && sigma_s > 1e-20 {
            let sigma_s_inv = (1.0 / sigma_s).min(1e20);
            // (2,2) block: -Σ_s^{-1} (always negative, correct for KKT inertia)
            matrix.add(n + i, n + i, -sigma_s_inv);
            // RHS: Σ_s^{-1} * (y + μ/s_l - μ/s_u) + infeasible contributions
            rhs[n + i] = sigma_s_inv * rhs_correction + rhs_infeasible;
        } else {
            // All infeasible: just drive toward feasibility
            rhs[n + i] = rhs_infeasible;
        }
    }

    // Debug: check for NaN in matrix and RHS
    if rhs.iter().any(|v| v.is_nan() || v.is_infinite()) {
        log::warn!("NaN/Inf in KKT RHS!");
        for (i, v) in rhs.iter().enumerate() {
            if v.is_nan() || v.is_infinite() {
                log::warn!("  rhs[{}] = {}", i, v);
            }
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
    /// Base constraint regularization.
    pub delta_c_base: f64,
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
            delta_c_base: 1e-8,
            delta_w_growth: 4.0,
            max_attempts: 15,
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

    // For unconstrained problems (m=0), try direct delta_w from min diagonal.
    // This avoids the exponential growth that overshoots on indefinite Hessians,
    // producing gradient-like steps instead of Newton steps.
    // Only use this when the required perturbation is moderate (< 1e4) —
    // extreme indefiniteness needs the standard exponential growth strategy.
    if m == 0 {
        if let Some(min_d) = solver.min_diagonal() {
            if min_d < 0.0 {
                let delta_w_direct = -min_d + 1e-8;
                let mut perturbed = kkt.matrix.clone();
                perturbed.add_diagonal_range(0, n, delta_w_direct);
                let inertia = solver.factor(&perturbed)?;
                if let Some(inertia) = inertia {
                    if inertia.positive == n && inertia.negative == 0 && inertia.zero == 0 {
                        kkt.matrix = perturbed;
                        params.delta_w_last = delta_w_direct;
                        return Ok((delta_w_direct, 0.0));
                    }
                }
            }
        }
    }

    // Inertia is wrong — apply perturbation and re-factor
    let mut delta_w = if params.delta_w_last == 0.0 {
        params.delta_w_init
    } else {
        (params.delta_w_last / params.delta_w_growth).max(params.delta_w_init)
    };
    let mut best_delta_w = delta_w;

    for attempt in 0..params.max_attempts {
        let delta_c = params.delta_c_base;

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

        best_delta_w = delta_w;

        // Increase perturbation
        delta_w *= params.delta_w_growth;

        log::debug!(
            "Inertia correction attempt {}: delta_w = {:.2e}, delta_c = {:.2e}, inertia = {:?}",
            attempt + 1,
            delta_w,
            delta_c,
            inertia
        );
    }

    // Inertia correction failed — use last perturbed matrix and proceed.
    // The line search will reject bad steps, and restoration can recover.
    let delta_c = params.delta_c_base;

    log::warn!(
        "Inertia correction failed after {} attempts (delta_w={:.2e}, delta_c={:.2e}), proceeding with approximate factorization",
        params.max_attempts, best_delta_w, delta_c
    );
    let mut perturbed = kkt.matrix.clone();
    perturbed.add_diagonal_range(0, n, best_delta_w);
    if m > 0 {
        perturbed.add_diagonal_range(n, n + m, -delta_c);
    }
    solver.factor(&perturbed)?;
    kkt.matrix = perturbed;
    params.delta_w_last = best_delta_w;
    Ok((best_delta_w, delta_c))
}

/// Solve the KKT system for the search direction, given a factored solver.
///
/// Returns (dx, dy) where dx is the primal step and dy is the dual step.
/// Bound multiplier steps dz_l, dz_u are recovered from complementarity.
///
/// Uses iterative refinement to improve solution accuracy for ill-conditioned
/// systems (e.g., near-singular Hessians in equality-constrained problems).
pub fn solve_for_direction(
    kkt: &KktSystem,
    solver: &mut dyn LinearSolver,
) -> Result<(Vec<f64>, Vec<f64>), crate::linear_solver::SolverError> {
    let dim = kkt.dim;

    // NaN guard on RHS — if the RHS has NaN, the problem evaluation is broken
    if kkt.rhs.iter().any(|v| v.is_nan() || v.is_infinite()) {
        return Err(crate::linear_solver::SolverError::NumericalFailure(
            "KKT RHS contains NaN/Inf".to_string(),
        ));
    }

    let mut solution = vec![0.0; dim];
    solver.solve(&kkt.rhs, &mut solution)?;

    // NaN guard on solution — factorization may produce NaN for ill-conditioned systems
    if solution.iter().any(|v| v.is_nan() || v.is_infinite()) {
        return Err(crate::linear_solver::SolverError::NumericalFailure(
            "KKT solution contains NaN/Inf".to_string(),
        ));
    }

    // Iterative refinement: correct the solution using the residual
    let max_refinements = 3;
    let mut residual = vec![0.0; dim];
    for _ref_iter in 0..max_refinements {
        // Compute residual: r = b - A*x
        kkt.matrix.matvec(&solution, &mut residual);
        let mut res_norm: f64 = 0.0;
        for i in 0..dim {
            residual[i] = kkt.rhs[i] - residual[i];
            res_norm = res_norm.max(residual[i].abs());
        }

        if res_norm < 1e-12 {
            break;
        }

        // Solve A * correction = residual
        let mut correction = vec![0.0; dim];
        if solver.solve(&residual, &mut correction).is_err() {
            break;
        }

        // Update solution
        for i in 0..dim {
            solution[i] += correction[i];
        }
    }

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

/// Condensed KKT system for m >> n problems (Schur complement).
///
/// Instead of factoring the full (n+m)×(n+m) KKT system, we condense to n×n:
///   S = H + Σ + δ_w·I + J^T · D_c^{-1} · J
///   S · dx = r_d + J^T · D_c^{-1} · r_p
///   dy = D_c^{-1} · (J · dx - r_p)
///
/// where D_c is the (2,2) block diagonal (negative for inequalities).
/// Cost: O(n²·m + n³) instead of O((n+m)³).
pub struct CondensedKktSystem {
    /// Condensed matrix S (n × n).
    pub matrix: SymmetricMatrix,
    /// Condensed RHS (n-vector).
    pub rhs: Vec<f64>,
    /// Number of primal variables.
    pub n: usize,
    /// Number of constraints.
    pub m: usize,
    /// D_c diagonal (m-vector, from the (2,2) block).
    pub d_c: Vec<f64>,
    /// Original primal RHS (n-vector).
    pub rhs_primal: Vec<f64>,
    /// Original constraint RHS (m-vector).
    pub rhs_constraint: Vec<f64>,
    /// Jacobian in COO format.
    pub jac_rows: Vec<usize>,
    pub jac_cols: Vec<usize>,
    pub jac_vals: Vec<f64>,
}

/// Assemble the condensed (Schur complement) KKT system.
///
/// Takes the same inputs as `assemble_kkt` but produces an n×n system
/// instead of (n+m)×(n+m).
#[allow(clippy::too_many_arguments)]
pub fn assemble_condensed_kkt(
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
    _v_l: &[f64],
    _v_u: &[f64],
) -> CondensedKktSystem {
    // Build the condensed system directly from problem data without assembling
    // the full (n+m)×(n+m) KKT matrix. This saves O((n+m)^2) memory and work.

    // --- (1,1) block: H + Sigma (n×n dense symmetric) ---
    let mut matrix = SymmetricMatrix::zeros(n);

    // Hessian entries
    for (idx, (&row, &col)) in hess_rows.iter().zip(hess_cols.iter()).enumerate() {
        matrix.add(row, col, hess_vals[idx]);
    }

    // Barrier diagonal Sigma for variable bounds
    for i in 0..n {
        matrix.add(i, i, sigma[i]);
    }

    // --- RHS: dual residual r_d (n-vector) ---
    let mut rhs_primal = vec![0.0; n];
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
        rhs_primal[i] = rd;
    }
    // Subtract J^T * y
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        rhs_primal[col] -= jac_vals[idx] * y[row];
    }

    // --- (2,2) block diagonal D_c and constraint RHS r_p (m-vectors) ---
    let mut d_c = vec![0.0; m];
    let mut rhs_constraint = vec![0.0; m];

    for i in 0..m {
        if is_equality_constraint(g_l[i], g_u[i]) {
            rhs_constraint[i] = -(g[i] - g_l[i]);
            // d_c[i] = 0.0 for equalities (no (2,2) block entry)
            continue;
        }

        let mut sigma_s = 0.0;
        let mut rhs_correction = y[i];
        let mut any_feasible = false;
        let mut rhs_infeasible = 0.0;

        if g_l[i].is_finite() {
            let slack = g[i] - g_l[i];
            if slack >= -1e-8 {
                let safe_slack = slack.max(mu.max(1e-10));
                let z_sl = if y[i] < -1e-20 { -y[i] } else { mu / safe_slack };
                sigma_s += z_sl / safe_slack;
                rhs_correction += mu / safe_slack;
                any_feasible = true;
            } else {
                rhs_infeasible += -(g[i] - g_l[i]);
            }
        }
        if g_u[i].is_finite() {
            let slack = g_u[i] - g[i];
            if slack >= -1e-8 {
                let safe_slack = slack.max(mu.max(1e-10));
                let z_su = if y[i] > 1e-20 { y[i] } else { mu / safe_slack };
                sigma_s += z_su / safe_slack;
                rhs_correction -= mu / safe_slack;
                any_feasible = true;
            } else {
                rhs_infeasible += -(g[i] - g_u[i]);
            }
        }

        if any_feasible && sigma_s > 1e-20 {
            let sigma_s_inv = (1.0 / sigma_s).min(1e20);
            d_c[i] = -sigma_s_inv;
            rhs_constraint[i] = sigma_s_inv * rhs_correction + rhs_infeasible;
        } else {
            rhs_constraint[i] = rhs_infeasible;
        }
    }

    // --- Build condensed matrix: S = (1,1) block + J^T · (-D_c)^{-1} · J ---
    let mut j_dense = vec![0.0; m * n];
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        j_dense[row * n + col] += jac_vals[idx];
    }

    for i in 0..m {
        let d_c_eff = if d_c[i].abs() < 1e-20 {
            -1e-16  // equality: very stiff spring (inv = -1e16)
        } else {
            d_c[i]
        };
        let inv_neg_dc = 1.0 / (-d_c_eff);
        for p in 0..n {
            let jp = j_dense[i * n + p];
            if jp == 0.0 {
                continue;
            }
            for q in 0..=p {
                let jq = j_dense[i * n + q];
                if jq != 0.0 {
                    matrix.add(p, q, inv_neg_dc * jp * jq);
                }
            }
        }
    }

    // --- Build condensed RHS: r_d + J^T · (-D_c)^{-1} · r_p ---
    let mut rhs = rhs_primal.clone();
    for i in 0..m {
        let d_c_eff = if d_c[i].abs() < 1e-20 {
            -1e-16
        } else {
            d_c[i]
        };
        let inv_neg_dc = 1.0 / (-d_c_eff);
        let scaled_rp = inv_neg_dc * rhs_constraint[i];
        for p in 0..n {
            let jp = j_dense[i * n + p];
            if jp != 0.0 {
                rhs[p] += jp * scaled_rp;
            }
        }
    }

    CondensedKktSystem {
        matrix,
        rhs,
        n,
        m,
        d_c,
        rhs_primal,
        rhs_constraint,
        jac_rows: jac_rows.to_vec(),
        jac_cols: jac_cols.to_vec(),
        jac_vals: jac_vals.to_vec(),
    }
}

/// Solve the condensed KKT system: compute dx from condensed, recover dy.
pub fn solve_condensed(
    condensed: &CondensedKktSystem,
    solver: &mut dyn LinearSolver,
) -> Result<(Vec<f64>, Vec<f64>), SolverError> {
    let n = condensed.n;
    let m = condensed.m;

    // Solve S · dx = rhs_condensed
    let mut dx = vec![0.0; n];
    solver.solve(&condensed.rhs, &mut dx)?;

    // Recover dy = (-D_c)^{-1} · (J · dx - r_p)
    // First compute J · dx
    let mut jdx = vec![0.0; m];
    for (idx, (&row, &col)) in condensed
        .jac_rows
        .iter()
        .zip(condensed.jac_cols.iter())
        .enumerate()
    {
        jdx[row] += condensed.jac_vals[idx] * dx[col];
    }

    let mut dy = vec![0.0; m];
    for i in 0..m {
        let d_c_eff = if condensed.d_c[i].abs() < 1e-20 {
            -1e-16  // equality: consistent with assembly
        } else {
            condensed.d_c[i]
        };
        dy[i] = (jdx[i] - condensed.rhs_constraint[i]) / (-d_c_eff);
    }

    Ok((dx, dy))
}

/// Solve the condensed system with a modified constraint residual (for SOC).
///
/// Instead of using the original rhs_constraint, uses -c_soc as the constraint
/// residual. This avoids rebuilding the full (n+m)×(n+m) KKT system for SOC.
pub fn solve_condensed_soc(
    condensed: &CondensedKktSystem,
    solver: &mut dyn LinearSolver,
    c_soc: &[f64],
) -> Result<Vec<f64>, SolverError> {
    let n = condensed.n;
    let m = condensed.m;

    // Build per-constraint scaling: (-c_soc[i]) / (-d_c[i])
    let mut scaled = vec![0.0; m];
    for i in 0..m {
        let d_c_eff = if condensed.d_c[i].abs() < 1e-20 {
            -1e-16 // equality: consistent with assembly
        } else {
            condensed.d_c[i]
        };
        scaled[i] = (-c_soc[i]) / (-d_c_eff);
    }

    // Build modified condensed RHS: rhs_primal + J^T · scaled
    let mut rhs = condensed.rhs_primal.clone();
    for (idx, (&row, &col)) in condensed
        .jac_rows
        .iter()
        .zip(condensed.jac_cols.iter())
        .enumerate()
    {
        rhs[col] += condensed.jac_vals[idx] * scaled[row];
    }

    // Solve S · dx = modified_rhs
    let mut dx = vec![0.0; n];
    solver.solve(&rhs, &mut dx)?;

    Ok(dx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_solver::dense::DenseLdl;

    #[test]
    fn test_compute_sigma_no_bounds() {
        let x = vec![1.0, 2.0];
        let x_l = vec![f64::NEG_INFINITY, f64::NEG_INFINITY];
        let x_u = vec![f64::INFINITY, f64::INFINITY];
        let z_l = vec![0.0, 0.0];
        let z_u = vec![0.0, 0.0];
        let sigma = compute_sigma(&x, &x_l, &x_u, &z_l, &z_u);
        assert!((sigma[0]).abs() < 1e-15);
        assert!((sigma[1]).abs() < 1e-15);
    }

    #[test]
    fn test_compute_sigma_lower_bound_only() {
        let x = vec![1.5];
        let x_l = vec![1.0];
        let x_u = vec![f64::INFINITY];
        let z_l = vec![2.0];
        let z_u = vec![0.0];
        let sigma = compute_sigma(&x, &x_l, &x_u, &z_l, &z_u);
        // sigma = z_l / (x - x_l) = 2.0 / 0.5 = 4.0
        assert!((sigma[0] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_sigma_both_bounds() {
        let x = vec![1.5];
        let x_l = vec![1.0];
        let x_u = vec![2.0];
        let z_l = vec![2.0];
        let z_u = vec![3.0];
        let sigma = compute_sigma(&x, &x_l, &x_u, &z_l, &z_u);
        // sigma = 2.0/0.5 + 3.0/0.5 = 4.0 + 6.0 = 10.0
        assert!((sigma[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_compute_sigma_at_bound_clamped() {
        let x = vec![1.0]; // At lower bound
        let x_l = vec![1.0];
        let x_u = vec![f64::INFINITY];
        let z_l = vec![1.0];
        let z_u = vec![0.0];
        let sigma = compute_sigma(&x, &x_l, &x_u, &z_l, &z_u);
        // slack = max(0, 1e-20) = 1e-20, sigma = 1.0/1e-20 = 1e20
        assert!(sigma[0] > 1e19);
    }

    #[test]
    fn test_assemble_kkt_unconstrained() {
        // 2 vars, no constraints
        // H = [[2, 0], [0, 3]]
        let n = 2;
        let m = 0;
        let hess_rows = vec![0, 1];
        let hess_cols = vec![0, 1];
        let hess_vals = vec![2.0, 3.0];
        let sigma = vec![1.0, 2.0];
        let grad_f = vec![0.5, 0.5];
        let x = vec![1.0, 2.0];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let z_l = vec![0.0; 2];
        let z_u = vec![0.0; 2];

        let kkt = assemble_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &[], &[], &[], &sigma, &grad_f,
            &[], &[], &[], &[], &z_l, &z_u,
            &x, &x_l, &x_u, 0.1, false, &[], &[],
        );

        assert_eq!(kkt.dim, 2);
        // (1,1) block: H + Sigma = [[3, 0], [0, 5]]
        assert!((kkt.matrix.get(0, 0) - 3.0).abs() < 1e-12);
        assert!((kkt.matrix.get(1, 1) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_assemble_kkt_equality_constraint() {
        // 2 vars, 1 equality constraint: x0 + x1 = 1
        let n = 2;
        let m = 1;
        let hess_rows = vec![0, 1];
        let hess_cols = vec![0, 1];
        let hess_vals = vec![2.0, 2.0];
        let jac_rows = vec![0, 0];
        let jac_cols = vec![0, 1];
        let jac_vals = vec![1.0, 1.0];
        let sigma = vec![0.0; 2];
        let grad_f = vec![1.0, 1.0];
        let g = vec![0.7]; // current constraint value
        let g_l = vec![1.0];
        let g_u = vec![1.0];
        let y = vec![0.5];
        let x = vec![0.3, 0.4];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let z_l = vec![0.0; 2];
        let z_u = vec![0.0; 2];

        let v_l = vec![0.0; m];
        let v_u = vec![0.0; m];
        let kkt = assemble_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &jac_rows, &jac_cols, &jac_vals, &sigma, &grad_f,
            &g, &g_l, &g_u, &y, &z_l, &z_u,
            &x, &x_l, &x_u, 0.1, false, &v_l, &v_u,
        );

        assert_eq!(kkt.dim, 3);
        // Verify J block: matrix[2,0] and matrix[2,1] should be 1.0
        assert!((kkt.matrix.get(2, 0) - 1.0).abs() < 1e-12);
        assert!((kkt.matrix.get(2, 1) - 1.0).abs() < 1e-12);
        // Equality constraint: (2,2) block should be 0
        assert!((kkt.matrix.get(2, 2)).abs() < 1e-12);
        // Primal residual: -(g - g_l) = -(0.7 - 1.0) = 0.3
        assert!((kkt.rhs[2] - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_assemble_kkt_rhs_sign_convention() {
        // Regression: J^T*y must be subtracted from r_d
        let n = 1;
        let m = 1;
        let hess_rows = vec![0];
        let hess_cols = vec![0];
        let hess_vals = vec![1.0];
        let jac_rows = vec![0];
        let jac_cols = vec![0];
        let jac_vals = vec![2.0]; // J = [2]
        let sigma = vec![0.0];
        let grad_f = vec![3.0];
        let g = vec![1.0];
        let g_l = vec![1.0];
        let g_u = vec![1.0];
        let y = vec![1.0];
        let x = vec![1.0];
        let x_l = vec![f64::NEG_INFINITY];
        let x_u = vec![f64::INFINITY];
        let z_l = vec![0.0];
        let z_u = vec![0.0];

        let v_l = vec![0.0; m];
        let v_u = vec![0.0; m];
        let kkt = assemble_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &jac_rows, &jac_cols, &jac_vals, &sigma, &grad_f,
            &g, &g_l, &g_u, &y, &z_l, &z_u,
            &x, &x_l, &x_u, 0.1, false, &v_l, &v_u,
        );

        // r_d = -grad_f + z_l - z_u = -3.0 + 0 - 0 = -3.0
        // Then subtract J^T * y: -3.0 - 2.0*1.0 = -5.0
        assert!((kkt.rhs[0] - (-5.0)).abs() < 1e-12,
            "RHS sign convention: expected -5.0, got {}", kkt.rhs[0]);
    }

    #[test]
    fn test_assemble_kkt_inequality_constraint() {
        // Feasible inequality: g(x) = 2.0, g_l = 1.0, g_u = INF
        let n = 1;
        let m = 1;
        let hess_rows = vec![0];
        let hess_cols = vec![0];
        let hess_vals = vec![1.0];
        let jac_rows = vec![0];
        let jac_cols = vec![0];
        let jac_vals = vec![1.0];
        let sigma = vec![0.0];
        let grad_f = vec![0.0];
        let g = vec![2.0]; // feasible: 2.0 > 1.0
        let g_l = vec![1.0];
        let g_u = vec![f64::INFINITY];
        let y = vec![0.0];
        let x = vec![2.0];
        let x_l = vec![f64::NEG_INFINITY];
        let x_u = vec![f64::INFINITY];
        let z_l = vec![0.0];
        let z_u = vec![0.0];
        let mu = 0.1;

        let v_l = vec![0.0; m];
        let v_u = vec![0.0; m];
        let kkt = assemble_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &jac_rows, &jac_cols, &jac_vals, &sigma, &grad_f,
            &g, &g_l, &g_u, &y, &z_l, &z_u,
            &x, &x_l, &x_u, mu, false, &v_l, &v_u,
        );

        // (2,2) block should be negative (from -Σ_s^{-1})
        assert!(kkt.matrix.get(1, 1) < 0.0,
            "Inequality (2,2) block should be negative, got {}", kkt.matrix.get(1, 1));
    }

    #[test]
    fn test_factor_with_inertia_correction_good() {
        // KKT matrix with correct inertia (2, 1, 0) — no perturbation needed
        let n = 2;
        let m = 1;
        let mut matrix = SymmetricMatrix::zeros(3);
        matrix.set(0, 0, 2.0);
        matrix.set(1, 1, 2.0);
        matrix.set(2, 0, 1.0);
        matrix.set(2, 1, 1.0);

        let mut kkt = KktSystem {
            dim: 3, n, m,
            matrix: KktMatrix::Dense(matrix),
            rhs: vec![1.0, 2.0, 3.0],
        };
        let mut solver = DenseLdl::new();
        let mut params = InertiaCorrectionParams::default();

        let (delta_w, delta_c) = factor_with_inertia_correction(&mut kkt, &mut solver, &mut params).unwrap();
        assert!((delta_w).abs() < 1e-15, "Good inertia should need no delta_w");
        assert!((delta_c).abs() < 1e-15, "Good inertia should need no delta_c");
    }

    #[test]
    fn test_factor_with_inertia_correction_needs_perturbation() {
        // Matrix with wrong inertia: all positive (3, 0, 0) instead of (2, 1, 0)
        let n = 2;
        let m = 1;
        let mut matrix = SymmetricMatrix::zeros(3);
        matrix.set(0, 0, 2.0);
        matrix.set(1, 1, 2.0);
        matrix.set(2, 2, 1.0); // Positive instead of 0
        matrix.set(2, 0, 0.1);
        matrix.set(2, 1, 0.1);

        let mut kkt = KktSystem {
            dim: 3, n, m,
            matrix: KktMatrix::Dense(matrix),
            rhs: vec![1.0, 2.0, 3.0],
        };
        let mut solver = DenseLdl::new();
        let mut params = InertiaCorrectionParams::default();

        let (delta_w, _delta_c) = factor_with_inertia_correction(&mut kkt, &mut solver, &mut params).unwrap();
        assert!(delta_w > 0.0, "Wrong inertia should require delta_w > 0");
    }

    #[test]
    fn test_factor_with_inertia_correction_warm_start() {
        // Use delta_w_last to warm-start perturbation
        let n = 2;
        let m = 1;
        let mut matrix = SymmetricMatrix::zeros(3);
        matrix.set(0, 0, 2.0);
        matrix.set(1, 1, 2.0);
        matrix.set(2, 2, 1.0);
        matrix.set(2, 0, 0.1);
        matrix.set(2, 1, 0.1);

        let mut kkt = KktSystem {
            dim: 3, n, m,
            matrix: KktMatrix::Dense(matrix),
            rhs: vec![1.0, 2.0, 3.0],
        };
        let mut solver = DenseLdl::new();
        let mut params = InertiaCorrectionParams::default();
        params.delta_w_last = 1.0; // Warm-start from previous

        let (delta_w, _) = factor_with_inertia_correction(&mut kkt, &mut solver, &mut params).unwrap();
        // Should start from delta_w_last / growth = 1.0 / 8.0 = 0.125
        assert!(delta_w >= 0.125 - 1e-10, "Warm-start should begin from delta_w_last/growth");
    }

    #[test]
    fn test_solve_for_direction_simple() {
        // Create a simple 2-var, 1-constraint KKT system and solve
        let n = 2;
        let m = 1;
        let mut matrix = SymmetricMatrix::zeros(3);
        matrix.set(0, 0, 2.0);
        matrix.set(1, 1, 2.0);
        matrix.set(2, 0, 1.0);
        matrix.set(2, 1, 1.0);

        let rhs = vec![1.0, 2.0, 0.5];
        let kkt = KktSystem {
            dim: 3, n, m,
            matrix: KktMatrix::Dense(matrix.clone()),
            rhs: rhs.clone(),
        };

        let mut solver = DenseLdl::new();
        solver.factor(&KktMatrix::Dense(matrix.clone())).unwrap();

        let (dx, dy) = solve_for_direction(&kkt, &mut solver).unwrap();
        assert_eq!(dx.len(), 2);
        assert_eq!(dy.len(), 1);

        // Verify KKT * [dx; dy] ≈ rhs
        let mut sol = vec![0.0; 3];
        sol[..2].copy_from_slice(&dx);
        sol[2] = dy[0];
        let mut ax = vec![0.0; 3];
        matrix.matvec(&sol, &mut ax);
        for i in 0..3 {
            assert!((ax[i] - rhs[i]).abs() < 1e-8,
                "KKT*solution mismatch at {}: {} vs {}", i, ax[i], rhs[i]);
        }
    }

    #[test]
    fn test_recover_dz_no_bounds() {
        let x = vec![1.0, 2.0];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let z_l = vec![0.0; 2];
        let z_u = vec![0.0; 2];
        let dx = vec![0.1, 0.2];
        let (dz_l, dz_u) = recover_dz(&x, &x_l, &x_u, &z_l, &z_u, &dx, 0.1);
        for i in 0..2 {
            assert!((dz_l[i]).abs() < 1e-15);
            assert!((dz_u[i]).abs() < 1e-15);
        }
    }

    #[test]
    fn test_recover_dz_lower_bound() {
        let x = vec![1.5];
        let x_l = vec![1.0];
        let x_u = vec![f64::INFINITY];
        let z_l = vec![2.0];
        let z_u = vec![0.0];
        let dx = vec![0.1];
        let mu = 0.1;
        let (dz_l, _) = recover_dz(&x, &x_l, &x_u, &z_l, &z_u, &dx, mu);
        // s_l = 0.5
        // dz_l = (mu - z_l*s_l)/s_l - (z_l/s_l)*dx
        //      = (0.1 - 2.0*0.5)/0.5 - (2.0/0.5)*0.1
        //      = (0.1 - 1.0)/0.5 - 4.0*0.1
        //      = -0.9/0.5 - 0.4
        //      = -1.8 - 0.4 = -2.2
        assert!((dz_l[0] - (-2.2)).abs() < 1e-12);
    }

    #[test]
    fn test_condensed_kkt_matches_full() {
        // 2 variables, 3 inequality constraints (m > n)
        // min x0^2 + x1^2 s.t. x0 + x1 >= 1, x0 >= 0.2, x1 >= 0.3
        let n = 2;
        let m = 3;
        let hess_rows = vec![0, 1];
        let hess_cols = vec![0, 1];
        let hess_vals = vec![2.0, 2.0];
        let jac_rows = vec![0, 0, 1, 2];
        let jac_cols = vec![0, 1, 0, 1];
        let jac_vals = vec![1.0, 1.0, 1.0, 1.0];
        let x = vec![0.6, 0.7];
        let x_l = vec![f64::NEG_INFINITY; 2];
        let x_u = vec![f64::INFINITY; 2];
        let z_l = vec![0.0; 2];
        let z_u = vec![0.0; 2];
        let sigma = compute_sigma(&x, &x_l, &x_u, &z_l, &z_u);
        let grad_f = vec![1.2, 1.4];
        let g = vec![1.3, 0.6, 0.7];
        let g_l = vec![1.0, 0.2, 0.3];
        let g_u = vec![f64::INFINITY; 3];
        let y = vec![0.1, 0.05, 0.05];
        let mu = 0.01;

        // Solve with full KKT
        let v_l = vec![0.0; m];
        let v_u = vec![0.0; m];
        let mut full_kkt = assemble_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &jac_rows, &jac_cols, &jac_vals, &sigma, &grad_f,
            &g, &g_l, &g_u, &y, &z_l, &z_u, &x, &x_l, &x_u, mu, false,
            &v_l, &v_u,
        );
        let mut full_solver = DenseLdl::new();
        let mut params = InertiaCorrectionParams::default();
        let _ = factor_with_inertia_correction(&mut full_kkt, &mut full_solver, &mut params);
        let (dx_full, dy_full) = solve_for_direction(&full_kkt, &mut full_solver).unwrap();

        // Solve with condensed KKT
        let condensed = assemble_condensed_kkt(
            n, m, &hess_rows, &hess_cols, &hess_vals,
            &jac_rows, &jac_cols, &jac_vals, &sigma, &grad_f,
            &g, &g_l, &g_u, &y, &z_l, &z_u, &x, &x_l, &x_u, mu,
            &v_l, &v_u,
        );
        let mut cond_solver = DenseLdl::new();
        cond_solver.factor(&KktMatrix::Dense(condensed.matrix.clone())).unwrap();
        let (dx_cond, dy_cond) = solve_condensed(&condensed, &mut cond_solver).unwrap();

        // Compare solutions
        for i in 0..n {
            assert!(
                (dx_full[i] - dx_cond[i]).abs() < 1e-6,
                "dx mismatch at {}: full={}, condensed={}", i, dx_full[i], dx_cond[i]
            );
        }
        for i in 0..m {
            assert!(
                (dy_full[i] - dy_cond[i]).abs() < 1e-6,
                "dy mismatch at {}: full={}, condensed={}", i, dy_full[i], dy_cond[i]
            );
        }
    }
}
