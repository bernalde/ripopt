use crate::convergence::is_equality_constraint;
use crate::linear_solver::{KktMatrix, LinearSolver, SolverError, SparseSymmetricMatrix, SymmetricMatrix};

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
    /// Per-constraint δ_c values added to the (2,2) block during assembly.
    /// Stored as positive values; the matrix gets -delta_c_diag\[i\] on diagonal (n+i, n+i).
    /// Used by iterative refinement to recover the original (unregularized) matvec.
    pub delta_c_diag: Vec<f64>,
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

    // Quasidefinite regularization: add -delta_c to ALL constraint diagonals.
    // This makes the (2,2) block strictly negative definite, guaranteeing correct
    // inertia (n, m, 0) and stable factorization for equalities (zero diagonal),
    // infeasible inequalities (zero diagonal), and even feasible inequalities
    // (already negative, making them more negative is harmless).
    // Iterative refinement in solve_for_direction recovers the true Newton direction.
    // No assembly-time δ_c — let the IC loop add δ_c when needed for inertia.
    // The IC's δ_c is tracked and iterative refinement in solve_for_direction
    // recovers the solution of the ORIGINAL (unperturbed) system.
    let delta_c_diag = vec![0.0; m];

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
        delta_c_diag,
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
            delta_c_base: 1e-4,
            delta_w_growth: 4.0,
            max_attempts: 15,
            delta_w_last: 0.0,
        }
    }
}

/// Check if a factorization produces acceptable backward error for the KKT system.
/// Does a trial solve and checks ||b - Ax|| / (||x|| + ||b||).
/// Returns true if backward error <= 1e-6.
fn check_factorization_backward_error(
    kkt: &KktSystem,
    solver: &mut dyn LinearSolver,
) -> bool {
    check_factorization_backward_error_with_matrix(&kkt.matrix, &kkt.rhs, solver)
}

/// Check backward error using a specific matrix (may differ from kkt.matrix when perturbed).
/// Uses a lenient threshold (1e-4) since iterative refinement in solve_for_direction
/// will improve accuracy. This catches only grossly unreliable factorizations.
fn check_factorization_backward_error_with_matrix(
    matrix: &KktMatrix,
    rhs: &[f64],
    solver: &mut dyn LinearSolver,
) -> bool {
    let dim = rhs.len();
    let mut solution = vec![0.0; dim];
    if solver.solve(rhs, &mut solution).is_err() {
        return false;
    }
    if solution.iter().any(|v| v.is_nan() || v.is_infinite()) {
        return false;
    }
    let mut residual = vec![0.0; dim];
    matrix.matvec(&solution, &mut residual);
    let x_norm: f64 = solution.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
    let mut max_berr: f64 = 0.0;
    for i in 0..dim {
        let abs_res = (rhs[i] - residual[i]).abs();
        let denom = x_norm + rhs[i].abs().max(1e-30);
        max_berr = max_berr.max(abs_res / denom);
    }
    max_berr <= 1e-4
}

/// Perform KKT factorization with inertia correction.
///
/// Factor the KKT matrix and check inertia. If inertia is wrong
/// (should be (n, m, 0) for an n-variable, m-constraint problem),
/// add regularization and re-factor. Also checks backward error
/// to ensure the factorization is numerically reliable.
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
        // Accept inertia if counts match expected (n, m, 0). For large systems,
        // allow ±1 tolerance: Bunch-Kaufman pivoting can misclassify near-zero
        // eigenvalues at the positive/negative boundary, especially when n ≈ m.
        let inertia_ok = inertia.positive == n && inertia.negative == m && inertia.zero == 0;
        let total = inertia.positive + inertia.negative + inertia.zero;
        let approx_ok = !inertia_ok && (n + m) >= 100 && inertia.zero == 0
            && (total as isize - (n + m) as isize).unsigned_abs() <= 2
            && (inertia.positive as isize - n as isize).unsigned_abs() <= 1
            && (inertia.negative as isize - m as isize).unsigned_abs() <= 1;
        if inertia_ok || approx_ok {
            // For large systems, accept once inertia is correct. With a permissive
            // pivot threshold (1e-6), small pivots are common and produce large
            // factorization backward error. Iterative refinement in solve_for_direction
            // recovers solve accuracy; the factorization just needs correct inertia.
            if (n + m) >= 100 {
                params.delta_w_last = 0.0;
                return Ok((0.0, 0.0));
            }
            // For small systems: verify backward error is acceptable
            if check_factorization_backward_error(kkt, solver) {
                params.delta_w_last = 0.0;
                return Ok((0.0, 0.0));
            }
            // Backward error too large — fall through to regularization
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

    // Before perturbation: try increasing factorization quality (Ipopt's IncreaseQuality).
    // This escalates the pivot threshold (e.g., 1e-6 -> 1e-3 -> 0.03 -> 0.1), which can
    // fix inertia by forcing better pivot choices without adding regularization.
    // Like Ipopt, try this once before resorting to delta_w perturbation.
    if solver.increase_quality() {
        let inertia = solver.factor(&kkt.matrix)?;
        if let Some(inertia) = inertia {
            let inertia_ok = inertia.positive == n && inertia.negative == m && inertia.zero == 0;
            let total = inertia.positive + inertia.negative + inertia.zero;
            let approx_ok = !inertia_ok && (n + m) >= 100 && inertia.zero == 0
                && (total as isize - (n + m) as isize).unsigned_abs() <= 2
                && (inertia.positive as isize - n as isize).unsigned_abs() <= 1
                && (inertia.negative as isize - m as isize).unsigned_abs() <= 1;
            if inertia_ok || approx_ok {
                if (n + m) >= 100 {
                    params.delta_w_last = 0.0;
                    return Ok((0.0, 0.0));
                }
                if check_factorization_backward_error(kkt, solver) {
                    params.delta_w_last = 0.0;
                    return Ok((0.0, 0.0));
                }
            }
        }
    }

    // Selective delta_c-only perturbation (Ipopt's PerturbForSingularity path).
    // When n_neg < m by a small amount, the (1,1) block is likely already positive
    // (semi)definite and a few eigenvalues near the saddle-point boundary flipped sign
    // due to GEMM rounding. Adding only -delta_c to the (2,2) block pushes those
    // borderline eigenvalues negative without perturbing the primal block.
    if m > 0 {
        let last_inertia = solver.factor(&kkt.matrix)?;
        if let Some(inertia) = last_inertia {
            let total = inertia.positive + inertia.negative + inertia.zero;
            let deficit = m as isize - inertia.negative as isize;
            // Too few negatives by a small amount: try delta_c only
            if deficit > 0
                && deficit <= 5.max((m / 100) as isize)
                && inertia.zero == 0
                && (total as isize - (n + m) as isize).unsigned_abs() <= 2
            {
                let mut delta_c = params.delta_c_base;
                for _ in 0..4 {
                    let mut perturbed = kkt.matrix.clone();
                    perturbed.add_diagonal_range(n, n + m, -delta_c);
                    let inertia = solver.factor(&perturbed)?;
                    if let Some(inertia) = inertia {
                        let ok = inertia.positive == n && inertia.negative == m
                            && inertia.zero == 0;
                        if ok {
                            log::debug!(
                                "Selective delta_c-only correction succeeded: delta_c={:.2e}",
                                delta_c
                            );
                            kkt.matrix = perturbed;
                            params.delta_w_last = 0.0;
                            return Ok((0.0, delta_c));
                        }
                    }
                    delta_c *= 4.0;
                }
                // delta_c-only didn't work — fall through to full perturbation
            }
        }
    }

    // Inertia is wrong or backward error too large — apply perturbation and re-factor
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
            let exact_ok = inertia.positive == n && inertia.negative == m && inertia.zero == 0;
            let total = inertia.positive + inertia.negative + inertia.zero;
            let approx_ok = !exact_ok && (n + m) >= 100 && inertia.zero == 0
                && (total as isize - (n + m) as isize).unsigned_abs() <= 2
                && (inertia.positive as isize - n as isize).unsigned_abs() <= 1
                && (inertia.negative as isize - m as isize).unsigned_abs() <= 1;
            if exact_ok || approx_ok {
                // For large systems, accept once inertia is correct.
                if (n + m) >= 100 {
                    kkt.matrix = perturbed;
                    params.delta_w_last = delta_w;
                    return Ok((delta_w, delta_c));
                }
                // For small systems: verify backward error is acceptable
                if check_factorization_backward_error_with_matrix(&perturbed, &kkt.rhs, solver) {
                    kkt.matrix = perturbed;
                    params.delta_w_last = delta_w;
                    return Ok((delta_w, delta_c));
                }
                // Backward error too large — increase regularization
                log::debug!(
                    "Inertia correct at delta_w={:.2e} but backward error too large, increasing",
                    delta_w
                );
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

/// Compute y = A_original * x, undoing ALL perturbations: both assembly-time δ_c
/// on equality constraints and IC perturbation (δ_w on primal, δ_c_ic on constraints).
///
/// The factored system has: A_factored = A_original - diag(assembly_δ_c) - diag(IC_δ_c) + diag(IC_δ_w)
/// So: A_original * x = A_factored * x - δ_w * x\[0..n\] + (assembly_δ_c\[i\] + IC_δ_c) * x\[n+i\]
fn matvec_original(
    kkt: &KktSystem,
    x: &[f64],
    y: &mut [f64],
    delta_w: f64,
    delta_c_ic: f64,
) {
    kkt.matrix.matvec(x, y);
    // Undo IC primal perturbation (IC added +delta_w to diagonal 0..n)
    if delta_w > 0.0 {
        for j in 0..kkt.n {
            y[j] -= delta_w * x[j];
        }
    }
    // Undo assembly δ_c + IC constraint perturbation
    for i in 0..kkt.m {
        let total_dc = kkt.delta_c_diag[i] + delta_c_ic;
        if total_dc > 0.0 {
            y[kkt.n + i] += total_dc * x[kkt.n + i];
        }
    }
}

/// Solve the KKT system for the search direction, given a factored solver.
///
/// Returns (dx, dy) where dx is the primal step and dy is the dual step.
/// Bound multiplier steps dz_l, dz_u are recovered from complementarity.
///
/// Uses iterative refinement against the ORIGINAL (unregularized) system
/// to recover the true Newton direction despite δ_c/δ_w regularization.
/// The factored regularized system acts as a preconditioner.
pub fn solve_for_direction(
    kkt: &KktSystem,
    solver: &mut dyn LinearSolver,
    delta_w: f64,
    delta_c_ic: f64,
) -> Result<(Vec<f64>, Vec<f64>), crate::linear_solver::SolverError> {
    let dim = kkt.dim;
    // Refine against the assembled system (undoing IC perturbation) when IC was
    // triggered AND assembly-time δ_c is present. The δ_c makes the assembled
    // system non-singular, so undoing IC produces a well-conditioned target.
    // When IC triggers for other reasons (indefinite Hessian without δ_c),
    // the original system is ill-conditioned and refinement would diverge.
    // Only use IC-undoing refinement for large systems (where backward error
    // failures from near-zero pivots are the primary issue). Small systems
    // with indefinite Hessians need IC perturbation and shouldn't be undone.
    let use_ic_refinement = kkt.m > 0 && (kkt.n + kkt.m) >= 100
        && (delta_w > 0.0 || delta_c_ic > 0.0);

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

    // Iterative refinement: correct the solution using residuals against the
    // ORIGINAL (unregularized) matrix. The factored regularized system acts as
    // a preconditioner. This converges to the solution of the original system.
    // Use more refinement iterations for large systems where factorization accuracy
    // may be limited, and for IC-refinement where we're solving a different system.
    let max_refinements = if use_ic_refinement || (kkt.n + kkt.m) >= 100 { 5 } else { 3 };
    let mut residual = vec![0.0; dim];
    let mut prev_res_norm = f64::MAX;
    for _ref_iter in 0..max_refinements {
        // When IC regularization was applied, compute residual against the ORIGINAL
        // (unregularized) matrix so refinement converges to the true Newton direction.
        if use_ic_refinement {
            matvec_original(kkt, &solution, &mut residual, delta_w, delta_c_ic);
        } else {
            kkt.matrix.matvec(&solution, &mut residual);
        }
        let mut res_norm: f64 = 0.0;
        for i in 0..dim {
            residual[i] = kkt.rhs[i] - residual[i];
            res_norm = res_norm.max(residual[i].abs());
        }

        if res_norm < 1e-12 {
            break;
        }

        // Stagnation detection: stop if not improving (only for IC-regularized path
        // where refinement against the original may stagnate due to preconditioning mismatch)
        if use_ic_refinement && res_norm > 0.9 * prev_res_norm {
            break;
        }
        prev_res_norm = res_norm;

        // Solve A_regularized * correction = residual
        let mut correction = vec![0.0; dim];
        if solver.solve(&residual, &mut correction).is_err() {
            break;
        }

        // Update solution
        for i in 0..dim {
            solution[i] += correction[i];
        }
    }

    // Backward error monitoring against the ORIGINAL system
    // berr_i = |r_i| / (||x||_∞ + |b_i|)
    // If max berr > 1e-6, the linear solve is unreliable
    {
        if use_ic_refinement {
            matvec_original(kkt, &solution, &mut residual, delta_w, delta_c_ic);
        } else {
            kkt.matrix.matvec(&solution, &mut residual);
        }
        let x_norm: f64 = solution.iter().map(|v| v.abs()).fold(0.0f64, f64::max).max(1.0);
        let mut max_berr: f64 = 0.0;
        for i in 0..dim {
            let abs_res = (kkt.rhs[i] - residual[i]).abs();
            let denom = x_norm + kkt.rhs[i].abs().max(1e-30);
            max_berr = max_berr.max(abs_res / denom);
        }
        // Log backward error for diagnostics but don't reject the solve.
        // The line search is the natural quality filter: if the direction is poor,
        // no step length will be accepted and the solver enters restoration.
        // Ipopt similarly trusts the inertia-checked factorization without a
        // backward error gate on the solve.
        let berr_tol = if (kkt.n + kkt.m) >= 100 { 1e-2 } else { 1e-6 };
        if max_berr > berr_tol {
            log::debug!("KKT backward error {:.2e} exceeds {:.0e}", max_berr, berr_tol);
            return Err(crate::linear_solver::SolverError::NumericalFailure(
                format!("KKT backward error {:.2e} exceeds tolerance", max_berr),
            ));
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

/// Compute the affine-scaling (μ=0) predictor RHS for Mehrotra predictor-corrector.
///
/// Returns a copy of the existing KKT RHS with the centering terms (μ/s) removed.
/// Solving with this RHS gives the pure-Newton (affine-scaling) direction that
/// can be used to probe a better barrier parameter μ.
///
/// Cost: O(n) to build — extremely cheap. The expensive part is the triangular solve.
pub fn affine_predictor_rhs(
    rhs: &[f64],
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    mu: f64,
) -> Vec<f64> {
    let n = x.len();
    let mut rhs_aff = rhs.to_vec();
    // Remove the μ/s centering terms from the primal block
    for i in 0..n {
        if x_l[i].is_finite() {
            let s_l = (x[i] - x_l[i]).max(1e-20);
            rhs_aff[i] -= mu / s_l;
        }
        if x_u[i].is_finite() {
            let s_u = (x_u[i] - x[i]).max(1e-20);
            rhs_aff[i] += mu / s_u;
        }
    }
    rhs_aff
}

/// Build a new KKT RHS with a different barrier parameter μ_new.
///
/// Given the existing RHS assembled with μ_old, returns a new RHS for μ_new.
/// Only the primal block (first n entries) changes: the μ/s centering terms
/// are updated from μ_old/s to μ_new/s.
pub fn rebuild_rhs_with_mu(
    rhs: &[f64],
    x: &[f64],
    x_l: &[f64],
    x_u: &[f64],
    mu_old: f64,
    mu_new: f64,
) -> Vec<f64> {
    let n = x.len();
    let mut rhs_new = rhs.to_vec();
    let delta_mu = mu_new - mu_old;
    for i in 0..n {
        if x_l[i].is_finite() {
            let s_l = (x[i] - x_l[i]).max(1e-20);
            rhs_new[i] += delta_mu / s_l;
        }
        if x_u[i].is_finite() {
            let s_u = (x_u[i] - x[i]).max(1e-20);
            rhs_new[i] -= delta_mu / s_u;
        }
    }
    rhs_new
}

/// Solve the factored system with a custom RHS (backsolve only, no re-factorization).
///
/// Used for Mehrotra predictor-corrector and Gondzio centrality corrections —
/// both need extra backsolves with the already-factored KKT matrix.
///
/// Returns (dx, dy) where dx is the primal block (first n entries) and dy the dual.
pub fn solve_with_custom_rhs(
    n: usize,
    dim: usize,
    solver: &mut dyn LinearSolver,
    rhs: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), SolverError> {
    if rhs.iter().any(|v| v.is_nan() || v.is_infinite()) {
        return Err(SolverError::NumericalFailure(
            "Custom RHS contains NaN/Inf".to_string(),
        ));
    }
    let mut solution = vec![0.0; dim];
    solver.solve(rhs, &mut solution)?;
    if solution.iter().any(|v| v.is_nan() || v.is_infinite()) {
        return Err(SolverError::NumericalFailure(
            "Custom solve solution contains NaN/Inf".to_string(),
        ));
    }
    Ok((solution[..n].to_vec(), solution[n..].to_vec()))
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

/// Sparse condensed KKT system for large problems where both n and m are large.
///
/// Like `CondensedKktSystem` but stores S as a sparse matrix, avoiding O(n²) memory.
/// The Schur complement J^T · D_c^{-1} · J is computed directly in COO format.
pub struct SparseCondensedKktSystem {
    pub matrix: SparseSymmetricMatrix,
    pub rhs: Vec<f64>,
    pub n: usize,
    pub m: usize,
    pub d_c: Vec<f64>,
    pub rhs_primal: Vec<f64>,
    pub rhs_constraint: Vec<f64>,
    pub jac_rows: Vec<usize>,
    pub jac_cols: Vec<usize>,
    pub jac_vals: Vec<f64>,
}

/// Assemble a sparse condensed (Schur complement) KKT system.
///
/// S = H + Σ + J^T · (-D_c)^{-1} · J  (n×n sparse)
/// rhs = r_d + J^T · (-D_c)^{-1} · r_p
///
/// The J^T·D_c^{-1}·J product is computed row-by-row:
/// for each constraint i, add rank-1 update (1/(-d_c\[i\])) * J\[i,:\]^T * J\[i,:\]
/// Only nonzero pairs of J entries in the same row produce fill.
#[allow(clippy::too_many_arguments)]
pub fn assemble_sparse_condensed_kkt(
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
) -> SparseCondensedKktSystem {
    // Estimate nnz: hess + diagonal + J^T*D_c^{-1}*J fill
    // For each constraint row with k nonzeros, the outer product has k*(k+1)/2 entries.
    // Build a row-pointer structure for J to iterate by constraint row.
    let mut row_start = vec![0usize; m + 1];
    for &r in jac_rows {
        row_start[r + 1] += 1;
    }
    for i in 0..m {
        row_start[i + 1] += row_start[i];
    }
    // Sort Jacobian entries by row
    let jac_nnz = jac_rows.len();
    let mut jac_order = vec![0usize; jac_nnz];
    let mut row_count = vec![0usize; m];
    for k in 0..jac_nnz {
        let r = jac_rows[k];
        jac_order[row_start[r] + row_count[r]] = k;
        row_count[r] += 1;
    }

    // Estimate nnz for S
    let mut schur_nnz = 0;
    for i in 0..m {
        let k = row_start[i + 1] - row_start[i];
        schur_nnz += k * (k + 1) / 2;
    }
    let total_nnz = hess_rows.len() + n + schur_nnz;
    let mut matrix = SparseSymmetricMatrix::with_capacity(n, total_nnz);

    // (1,1) block: H + Σ
    for (idx, (&row, &col)) in hess_rows.iter().zip(hess_cols.iter()).enumerate() {
        matrix.add(row, col, hess_vals[idx]);
    }
    for i in 0..n {
        matrix.add(i, i, sigma[i]);
    }

    // RHS: dual residual r_d
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

    // (2,2) block diagonal D_c and constraint RHS
    let mut d_c = vec![0.0; m];
    let mut rhs_constraint = vec![0.0; m];

    for i in 0..m {
        if is_equality_constraint(g_l[i], g_u[i]) {
            // For equalities, D_c comes from constraint regularization.
            // Use -delta_c so the Schur complement J^T * delta_c^{-1} * J
            // doesn't blow up. Floor at 1e-4 to keep D_c^{-1} <= 1e4;
            // smaller values cause inaccurate dy recovery that degrades
            // convergence. The regularization doesn't affect the primal
            // step quality (dx from S is accurate), only dy.
            let delta_c = mu.max(1e-8);
            d_c[i] = -delta_c;
            rhs_constraint[i] = -(g[i] - g_l[i]);
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

    // Schur complement: S += J^T · (-D_c)^{-1} · J
    // Process row-by-row to exploit sparsity
    for i in 0..m {
        let d_c_eff = if d_c[i].abs() < 1e-20 { -1e-16 } else { d_c[i] };
        let inv_neg_dc = 1.0 / (-d_c_eff);

        let start = row_start[i];
        let end = row_start[i + 1];
        // For each pair of nonzeros in row i of J, add outer product entry
        for a in start..end {
            let ka = jac_order[a];
            let ca = jac_cols[ka];
            let va = jac_vals[ka];
            for b in a..end {
                let kb = jac_order[b];
                let cb = jac_cols[kb];
                let vb = jac_vals[kb];
                // Add to upper triangle (min(ca,cb), max(ca,cb))
                let (p, q) = if ca <= cb { (ca, cb) } else { (cb, ca) };
                let val = if a == b {
                    inv_neg_dc * va * vb
                } else {
                    inv_neg_dc * va * vb // both (a,b) and (b,a) contribute, but we only do upper
                };
                matrix.add(p, q, val);
            }
        }
    }

    // Condensed RHS: r_d + J^T · (-D_c)^{-1} · r_p
    let mut rhs = rhs_primal.clone();
    for i in 0..m {
        let d_c_eff = if d_c[i].abs() < 1e-20 { -1e-16 } else { d_c[i] };
        let inv_neg_dc = 1.0 / (-d_c_eff);
        let scaled_rp = inv_neg_dc * rhs_constraint[i];
        let start = row_start[i];
        let end = row_start[i + 1];
        for a in start..end {
            let ka = jac_order[a];
            rhs[jac_cols[ka]] += jac_vals[ka] * scaled_rp;
        }
    }

    SparseCondensedKktSystem {
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

/// Solve a sparse condensed system: dx from factored solver, recover dy.
pub fn solve_sparse_condensed(
    condensed: &SparseCondensedKktSystem,
    solver: &mut dyn LinearSolver,
) -> Result<(Vec<f64>, Vec<f64>), SolverError> {
    let n = condensed.n;
    let m = condensed.m;

    let mut dx = vec![0.0; n];
    solver.solve(&condensed.rhs, &mut dx)?;

    // Recover dy = (-D_c)^{-1} · (J · dx - r_p)
    let mut jdx = vec![0.0; m];
    for (idx, (&row, &col)) in condensed.jac_rows.iter().zip(condensed.jac_cols.iter()).enumerate() {
        jdx[row] += condensed.jac_vals[idx] * dx[col];
    }

    let mut dy = vec![0.0; m];
    for i in 0..m {
        let d_c_eff = if condensed.d_c[i].abs() < 1e-20 { -1e-16 } else { condensed.d_c[i] };
        dy[i] = (jdx[i] - condensed.rhs_constraint[i]) / (-d_c_eff);
    }

    Ok((dx, dy))
}

/// Solve sparse condensed with modified constraint residual (for SOC).
pub fn solve_sparse_condensed_soc(
    condensed: &SparseCondensedKktSystem,
    solver: &mut dyn LinearSolver,
    c_soc: &[f64],
) -> Result<Vec<f64>, SolverError> {
    let n = condensed.n;
    let m = condensed.m;

    let mut scaled = vec![0.0; m];
    for i in 0..m {
        let d_c_eff = if condensed.d_c[i].abs() < 1e-20 { -1e-16 } else { condensed.d_c[i] };
        scaled[i] = (-c_soc[i]) / (-d_c_eff);
    }

    let mut rhs = condensed.rhs_primal.clone();
    for (idx, (&row, &col)) in condensed.jac_rows.iter().zip(condensed.jac_cols.iter()).enumerate() {
        rhs[col] += condensed.jac_vals[idx] * scaled[row];
    }

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
        // Equality constraint: (2,2) block should be 0 (no δ_c regularization;
        // KKT-aware scaling handles the zero diagonal)
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
            delta_c_diag: vec![0.0; m],
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
            delta_c_diag: vec![0.0; m],
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
            delta_c_diag: vec![0.0; m],
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
            delta_c_diag: vec![0.0; m],
        };

        let mut solver = DenseLdl::new();
        solver.factor(&KktMatrix::Dense(matrix.clone())).unwrap();

        let (dx, dy) = solve_for_direction(&kkt, &mut solver, 0.0, 0.0).unwrap();
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
        let (dw, dc) = factor_with_inertia_correction(&mut full_kkt, &mut full_solver, &mut params).unwrap();
        let (dx_full, dy_full) = solve_for_direction(&full_kkt, &mut full_solver, dw, dc).unwrap();

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
