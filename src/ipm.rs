use crate::convergence::{
    self, check_convergence, complementarity_error, ConvergenceInfo, ConvergenceStatus,
};
use crate::filter::{self, Filter};
use crate::kkt::{self, InertiaCorrectionParams};
use crate::linear_solver::dense::DenseLdl;
use crate::linear_solver::{LinearSolver, SymmetricMatrix};
use crate::options::SolverOptions;
use crate::problem::NlpProblem;
use crate::restoration::RestorationPhase;
use crate::result::{SolveResult, SolveStatus};
use crate::warmstart::WarmStartInitializer;

/// Central state struct for the IPM solver.
pub(crate) struct SolverState {
    /// Current primal variables.
    pub x: Vec<f64>,
    /// Current constraint multipliers (lambda/y).
    pub y: Vec<f64>,
    /// Lower bound multipliers.
    pub z_l: Vec<f64>,
    /// Upper bound multipliers.
    pub z_u: Vec<f64>,
    /// Search direction: primal.
    pub dx: Vec<f64>,
    /// Search direction: constraint multipliers.
    pub dy: Vec<f64>,
    /// Search direction: lower bound multipliers.
    pub dz_l: Vec<f64>,
    /// Search direction: upper bound multipliers.
    pub dz_u: Vec<f64>,
    /// Barrier parameter.
    pub mu: f64,
    /// Primal step size.
    pub alpha_primal: f64,
    /// Dual step size.
    pub alpha_dual: f64,
    /// Iteration counter.
    pub iter: usize,
    /// Variable lower bounds.
    pub x_l: Vec<f64>,
    /// Variable upper bounds.
    pub x_u: Vec<f64>,
    /// Constraint lower bounds.
    pub g_l: Vec<f64>,
    /// Constraint upper bounds.
    pub g_u: Vec<f64>,
    /// Number of variables.
    pub n: usize,
    /// Number of constraints.
    pub m: usize,
    /// Current objective value.
    pub obj: f64,
    /// Current gradient.
    pub grad_f: Vec<f64>,
    /// Current constraint values.
    pub g: Vec<f64>,
    /// Jacobian structure and values.
    pub jac_rows: Vec<usize>,
    pub jac_cols: Vec<usize>,
    pub jac_vals: Vec<f64>,
    /// Hessian structure and values.
    pub hess_rows: Vec<usize>,
    pub hess_cols: Vec<usize>,
    pub hess_vals: Vec<f64>,
    /// Consecutive acceptable iterations.
    pub consecutive_acceptable: usize,
}

impl SolverState {
    /// Initialize from an NLP problem.
    fn new<P: NlpProblem>(problem: &P, options: &SolverOptions) -> Self {
        let n = problem.num_variables();
        let m = problem.num_constraints();

        let mut x_l = vec![0.0; n];
        let mut x_u = vec![0.0; n];
        problem.bounds(&mut x_l, &mut x_u);

        let mut g_l = vec![0.0; m];
        let mut g_u = vec![0.0; m];
        problem.constraint_bounds(&mut g_l, &mut g_u);

        let mut x = vec![0.0; n];
        problem.initial_point(&mut x);

        // Push initial point away from bounds
        for i in 0..n {
            if x_l[i].is_finite() && x_u[i].is_finite() {
                let range = x_u[i] - x_l[i];
                let push = options.bound_push.min(options.bound_frac * range);
                x[i] = x[i].max(x_l[i] + push).min(x_u[i] - push);
            } else if x_l[i].is_finite() {
                x[i] = x[i].max(x_l[i] + options.bound_push);
            } else if x_u[i].is_finite() {
                x[i] = x[i].min(x_u[i] - options.bound_push);
            }
        }

        // Initialize bound multipliers
        let mut z_l = vec![0.0; n];
        let mut z_u = vec![0.0; n];
        for i in 0..n {
            if x_l[i].is_finite() {
                let slack = (x[i] - x_l[i]).max(1e-20);
                z_l[i] = options.mu_init / slack;
            }
            if x_u[i].is_finite() {
                let slack = (x_u[i] - x[i]).max(1e-20);
                z_u[i] = options.mu_init / slack;
            }
        }

        let (jac_rows, jac_cols) = problem.jacobian_structure();
        let jac_nnz = jac_rows.len();
        let (hess_rows, hess_cols) = problem.hessian_structure();
        let hess_nnz = hess_rows.len();

        // Initialize constraint multipliers via least-squares estimate if enabled.
        // Solves min ||∇f + J^T y||^2  ⟹  (J J^T) y = -J ∇f
        let y = if options.least_squares_mult_init && m > 0 {
            let mut grad_f_init = vec![0.0; n];
            problem.gradient(&x, &mut grad_f_init);

            let mut jac_vals_init = vec![0.0; jac_nnz];
            problem.jacobian_values(&x, &mut jac_vals_init);

            // Compute b = -J * grad_f  (m-vector)
            let mut b = vec![0.0; m];
            for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
                b[row] -= jac_vals_init[idx] * grad_f_init[col];
            }

            // Compute A = J * J^T  (m x m dense symmetric matrix)
            // Build dense J (m x n) first for small m
            let mut j_dense = vec![0.0; m * n];
            for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
                j_dense[row * n + col] = jac_vals_init[idx];
            }
            let mut a_mat = SymmetricMatrix::zeros(m);
            for i in 0..m {
                for j in 0..=i {
                    let mut dot = 0.0;
                    for k in 0..n {
                        dot += j_dense[i * n + k] * j_dense[j * n + k];
                    }
                    a_mat.set(i, j, dot);
                }
            }

            // Solve (J J^T) y = b using DenseLdl
            let mut ls_solver = DenseLdl::new();
            let mut y_ls = vec![0.0; m];
            let factored = ls_solver.factor(&a_mat);
            let solved = factored.is_ok() && ls_solver.solve(&b, &mut y_ls).is_ok();

            if solved {
                let max_abs = y_ls.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
                if max_abs <= options.constr_mult_init_max {
                    y_ls
                } else {
                    vec![0.0; m]
                }
            } else {
                vec![0.0; m]
            }
        } else {
            vec![0.0; m]
        };

        Self {
            x,
            y,
            z_l,
            z_u,
            dx: vec![0.0; n],
            dy: vec![0.0; m],
            dz_l: vec![0.0; n],
            dz_u: vec![0.0; n],
            mu: options.mu_init,
            alpha_primal: 0.0,
            alpha_dual: 0.0,
            iter: 0,
            x_l,
            x_u,
            g_l,
            g_u,
            n,
            m,
            obj: 0.0,
            grad_f: vec![0.0; n],
            g: vec![0.0; m],
            jac_rows,
            jac_cols,
            jac_vals: vec![0.0; jac_nnz],
            hess_rows,
            hess_cols,
            hess_vals: vec![0.0; hess_nnz],
            consecutive_acceptable: 0,
        }
    }

    /// Evaluate all functions at the current point.
    fn evaluate<P: NlpProblem>(&mut self, problem: &P, obj_factor: f64) {
        self.obj = problem.objective(&self.x);
        problem.gradient(&self.x, &mut self.grad_f);
        if self.m > 0 {
            problem.constraints(&self.x, &mut self.g);
            problem.jacobian_values(&self.x, &mut self.jac_vals);
        }
        problem.hessian_values(&self.x, obj_factor, &self.y, &mut self.hess_vals);
    }

    /// Compute the barrier objective:
    /// f(x) - mu * sum(ln(x_i - x_l_i) + ln(x_u_i - x_i))
    /// Optionally includes constraint slack log-barriers when enabled.
    fn barrier_objective(&self, options: &SolverOptions) -> f64 {
        let mut phi = self.obj;
        for i in 0..self.n {
            if self.x_l[i].is_finite() {
                let slack = (self.x[i] - self.x_l[i]).max(1e-20);
                phi -= self.mu * slack.ln();
            }
            if self.x_u[i].is_finite() {
                let slack = (self.x_u[i] - self.x[i]).max(1e-20);
                phi -= self.mu * slack.ln();
            }
        }
        if options.constraint_slack_barrier {
            for i in 0..self.m {
                // Skip equality constraints (g_l == g_u): slack is zero by definition
                let is_eq = self.g_l[i].is_finite() && self.g_u[i].is_finite()
                    && (self.g_l[i] - self.g_u[i]).abs() < 1e-15;
                if is_eq {
                    continue;
                }
                if self.g_l[i].is_finite() {
                    let slack = self.g[i] - self.g_l[i];
                    if slack > self.mu * 1e-2 {
                        phi -= self.mu * slack.ln();
                    }
                }
                if self.g_u[i].is_finite() {
                    let slack = self.g_u[i] - self.g[i];
                    if slack > self.mu * 1e-2 {
                        phi -= self.mu * slack.ln();
                    }
                }
            }
        }
        phi
    }

    /// Compute constraint violation (theta).
    fn constraint_violation(&self) -> f64 {
        convergence::primal_infeasibility(&self.g, &self.g_l, &self.g_u)
    }

    /// Compute the directional derivative of the barrier objective along the search direction.
    ///
    /// ∇φ·dx = (∇f - μ/(x-x_l) + μ/(x_u-x))·dx
    /// Optionally includes constraint slack derivative terms when enabled.
    fn barrier_directional_derivative(&self, options: &SolverOptions) -> f64 {
        let mut grad_phi_dx = 0.0;
        for i in 0..self.n {
            let mut grad_phi_i = self.grad_f[i];
            if self.x_l[i].is_finite() {
                let slack = (self.x[i] - self.x_l[i]).max(1e-20);
                grad_phi_i -= self.mu / slack;
            }
            if self.x_u[i].is_finite() {
                let slack = (self.x_u[i] - self.x[i]).max(1e-20);
                grad_phi_i += self.mu / slack;
            }
            grad_phi_dx += grad_phi_i * self.dx[i];
        }
        if options.constraint_slack_barrier && self.m > 0 {
            // Compute J * dx (directional change in constraints)
            let mut jdx = vec![0.0; self.m];
            for (idx, (&row, &col)) in
                self.jac_rows.iter().zip(self.jac_cols.iter()).enumerate()
            {
                jdx[row] += self.jac_vals[idx] * self.dx[col];
            }
            for i in 0..self.m {
                let is_eq = self.g_l[i].is_finite() && self.g_u[i].is_finite()
                    && (self.g_l[i] - self.g_u[i]).abs() < 1e-15;
                if is_eq {
                    continue;
                }
                if self.g_l[i].is_finite() {
                    let slack = self.g[i] - self.g_l[i];
                    if slack > self.mu * 1e-2 {
                        grad_phi_dx -= self.mu * jdx[i] / slack;
                    }
                }
                if self.g_u[i].is_finite() {
                    let slack = self.g_u[i] - self.g[i];
                    if slack > self.mu * 1e-2 {
                        grad_phi_dx += self.mu * jdx[i] / slack;
                    }
                }
            }
        }
        grad_phi_dx
    }
}

/// Solve the NLP using the interior point method.
pub fn solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    let mut state = SolverState::new(problem, options);
    let n = state.n;
    let m = state.m;

    // Handle warm-start
    if options.warm_start {
        state.mu = WarmStartInitializer::initialize(
            &mut state.x,
            &mut state.z_l,
            &mut state.z_u,
            &state.x_l,
            &state.x_u,
            options,
        );
    }

    // Initialize linear solver
    let mut lin_solver = DenseLdl::new();
    let mut inertia_params = InertiaCorrectionParams::default();
    let mut restoration = RestorationPhase::new(50);

    // Initialize filter
    let mut filter = Filter::new(1e4);

    // Stall detection: count consecutive iterations with alpha_primal ≈ 0.
    // Progressive recovery: first reset filter only, then bump mu on repeated stalls.
    let mut consecutive_zero_alpha: usize = 0;
    let stall_threshold: usize = 5;
    let mut stall_recovery_count: usize = 0;

    // Initial evaluation
    state.evaluate(problem, 1.0);

    // Set filter parameters based on initial constraint violation
    let theta_init = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_init);

    // Print header
    if options.print_level >= 5 {
        log::info!(
            "{:>4} {:>14} {:>10} {:>10} {:>10} {:>10} {:>8} {:>8}",
            "iter",
            "objective",
            "inf_pr",
            "inf_du",
            "compl",
            "mu",
            "alpha_p",
            "alpha_d"
        );
    }

    // Main IPM loop
    for iteration in 0..options.max_iter {
        state.iter = iteration;

        // Compute optimality measures.
        let primal_inf = state.constraint_violation();

        // Compute z_opt from stationarity for the scaled convergence check.
        // At optimality, grad_f + J^T y - z_l + z_u = 0. For active bounds,
        // z_opt captures the true bound multiplier. For inactive bounds, z=0.
        //
        // Complementarity gate: only use z_opt when z_opt * slack is consistent
        // with the barrier problem (z*s ~ mu). If z_opt * slack >> mu, the point
        // is not a barrier-optimal point and z_opt would hide a true infeasibility.
        let (z_l_opt, z_u_opt) = {
            let mut grad_jty = state.grad_f.clone();
            for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
                grad_jty[col] += state.jac_vals[idx] * state.y[row];
            }
            let mut zl = vec![0.0; n];
            let mut zu = vec![0.0; n];
            let kappa_compl = 1e10; // allow z*s up to 1e10 * mu (matches kappa_sigma)
            for i in 0..n {
                if grad_jty[i] > 0.0 && state.x_l[i].is_finite() {
                    let s_l = (state.x[i] - state.x_l[i]).max(1e-20);
                    if grad_jty[i] * s_l <= kappa_compl * state.mu.max(1e-20) {
                        zl[i] = grad_jty[i];
                    }
                } else if grad_jty[i] < 0.0 && state.x_u[i].is_finite() {
                    let s_u = (state.x_u[i] - state.x[i]).max(1e-20);
                    if (-grad_jty[i]) * s_u <= kappa_compl * state.mu.max(1e-20) {
                        zu[i] = -grad_jty[i];
                    }
                }
            }
            (zl, zu)
        };

        let dual_inf = convergence::dual_infeasibility(
            &state.grad_f,
            &state.jac_rows,
            &state.jac_cols,
            &state.jac_vals,
            &state.y,
            &z_l_opt,
            &z_u_opt,
            n,
        );
        let compl_inf =
            complementarity_error(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0);

        if options.print_level >= 5 {
            log::info!(
                "{:>4} {:>14.7e} {:>10.2e} {:>10.2e} {:>10.2e} {:>10.2e} {:>8.2e} {:>8.2e}",
                iteration,
                state.obj,
                primal_inf,
                dual_inf,
                compl_inf,
                state.mu,
                state.alpha_primal,
                state.alpha_dual,
            );
        }

        // Compute multiplier scaling for convergence check
        let multiplier_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let multiplier_count = m + 2 * n;

        // Check convergence
        let conv_info = ConvergenceInfo {
            primal_inf,
            dual_inf,
            compl_inf,
            mu: state.mu,
            objective: state.obj,
            multiplier_sum,
            multiplier_count,
        };

        match check_convergence(&conv_info, options, state.consecutive_acceptable) {
            ConvergenceStatus::Converged => {
                return make_result(&state, SolveStatus::Optimal);
            }
            ConvergenceStatus::Acceptable => {
                return make_result(&state, SolveStatus::Acceptable);
            }
            ConvergenceStatus::Diverging => {
                return make_result(&state, SolveStatus::Unbounded);
            }
            ConvergenceStatus::NotConverged => {}
        }

        // Track acceptable iterations
        let acc_primal = primal_inf <= options.acceptable_tol;
        let acc_dual = dual_inf <= options.acceptable_tol;
        let acc_compl = compl_inf <= options.acceptable_tol;
        if acc_primal && acc_dual && acc_compl {
            state.consecutive_acceptable += 1;
        } else {
            state.consecutive_acceptable = 0;
        }

        // Compute sigma (barrier diagonal)
        let sigma = kkt::compute_sigma(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u);

        // Assemble and factor KKT system
        let mut kkt_system = kkt::assemble_kkt(
            n,
            m,
            &state.hess_rows,
            &state.hess_cols,
            &state.hess_vals,
            &state.jac_rows,
            &state.jac_cols,
            &state.jac_vals,
            &sigma,
            &state.grad_f,
            &state.g,
            &state.g_l,
            &state.g_u,
            &state.y,
            &state.z_l,
            &state.z_u,
            &state.x,
            &state.x_l,
            &state.x_u,
            state.mu,
        );

        // Factor with inertia correction
        let inertia_result =
            kkt::factor_with_inertia_correction(&mut kkt_system, &mut lin_solver, &mut inertia_params);

        if let Err(e) = inertia_result {
            log::warn!("KKT factorization failed: {}", e);
            // Try restoration instead of giving up
            let (x_rest, success) = restoration.restore(
                &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
                &state.jac_rows, &state.jac_cols, n, m, options,
                &|theta, phi| filter.is_acceptable(theta, phi),
                &|x_eval, g_out| problem.constraints(x_eval, g_out),
                &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
            );
            if success {
                state.x = x_rest;
                state.alpha_primal = 0.0;
                state.evaluate(problem, 1.0);
                continue;
            }
            return make_result(&state, SolveStatus::NumericalError);
        }

        // Solve for search direction
        let dir_result = kkt::solve_for_direction(&kkt_system, &mut lin_solver);
        let (dx, dy) = match dir_result {
            Ok(d) => d,
            Err(e) => {
                log::warn!("KKT solve failed: {}", e);
                // Try restoration instead of giving up
                let (x_rest, success) = restoration.restore(
                    &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
                    &state.jac_rows, &state.jac_cols, n, m, options,
                    &|theta, phi| filter.is_acceptable(theta, phi),
                    &|x_eval, g_out| problem.constraints(x_eval, g_out),
                    &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
                );
                if success {
                    state.x = x_rest;
                    state.alpha_primal = 0.0;
                    state.evaluate(problem, 1.0);
                    continue;
                }
                return make_result(&state, SolveStatus::NumericalError);
            }
        };

        // Recover bound multiplier steps
        let (dz_l, dz_u) =
            kkt::recover_dz(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, &dx, state.mu);

        state.dx = dx;
        state.dy = dy;
        state.dz_l = dz_l;
        state.dz_u = dz_u;

        // Compute maximum step sizes using fraction-to-boundary rule
        let tau = (1.0 - state.mu).max(options.tau_min);

        // Primal step: ensure x + alpha*dx stays within variable bounds
        let mut alpha_primal_max: f64 = 1.0;
        for i in 0..n {
            if state.x_l[i].is_finite() && state.dx[i] < 0.0 {
                let slack = state.x[i] - state.x_l[i];
                let ratio = -tau * slack / state.dx[i];
                alpha_primal_max = alpha_primal_max.min(ratio);
            }
            if state.x_u[i].is_finite() && state.dx[i] > 0.0 {
                let slack = state.x_u[i] - state.x[i];
                let ratio = tau * slack / state.dx[i];
                alpha_primal_max = alpha_primal_max.min(ratio);
            }
        }

        alpha_primal_max = alpha_primal_max.clamp(0.0, 1.0);

        // Dual step: ensure z + alpha*dz > 0
        let alpha_dual_max_l = filter::fraction_to_boundary(&state.z_l, &state.dz_l, tau);
        let alpha_dual_max_u = filter::fraction_to_boundary(&state.z_u, &state.dz_u, tau);
        let alpha_dual_max = alpha_dual_max_l.min(alpha_dual_max_u);

        // Line search
        let theta_current = primal_inf;
        let phi_current = state.barrier_objective(options);
        let grad_phi_step = state.barrier_directional_derivative(options);

        let mut alpha = alpha_primal_max;
        let mut step_accepted = false;
        let min_alpha = 1e-12;

        for _ls_iter in 0..40 {
            if alpha < min_alpha {
                break;
            }

            // Compute trial point
            let mut x_trial = vec![0.0; n];
            #[allow(clippy::needless_range_loop)]
            for i in 0..n {
                x_trial[i] = state.x[i] + alpha * state.dx[i];
                // Safeguard: ensure strictly within bounds
                if state.x_l[i].is_finite() {
                    x_trial[i] = x_trial[i].max(state.x_l[i] + 1e-14);
                }
                if state.x_u[i].is_finite() {
                    x_trial[i] = x_trial[i].min(state.x_u[i] - 1e-14);
                }
            }

            // Evaluate at trial point
            let obj_trial = problem.objective(&x_trial);
            let mut g_trial = vec![0.0; m];
            if m > 0 {
                problem.constraints(&x_trial, &mut g_trial);
            }

            let theta_trial =
                convergence::primal_infeasibility(&g_trial, &state.g_l, &state.g_u);

            // Compute barrier objective at trial
            let mut phi_trial = obj_trial;
            #[allow(clippy::needless_range_loop)]
            for i in 0..n {
                if state.x_l[i].is_finite() {
                    let slack = (x_trial[i] - state.x_l[i]).max(1e-20);
                    phi_trial -= state.mu * slack.ln();
                }
                if state.x_u[i].is_finite() {
                    let slack = (state.x_u[i] - x_trial[i]).max(1e-20);
                    phi_trial -= state.mu * slack.ln();
                }
            }
            if options.constraint_slack_barrier {
                for i in 0..m {
                    let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                        && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
                    if is_eq {
                        continue;
                    }
                    if state.g_l[i].is_finite() {
                        let slack = g_trial[i] - state.g_l[i];
                        if slack > state.mu * 1e-2 {
                            phi_trial -= state.mu * slack.ln();
                        }
                    }
                    if state.g_u[i].is_finite() {
                        let slack = state.g_u[i] - g_trial[i];
                        if slack > state.mu * 1e-2 {
                            phi_trial -= state.mu * slack.ln();
                        }
                    }
                }
            }

            // Check acceptability
            let (acceptable, _used_switching) = filter.check_acceptability(
                theta_current,
                phi_current,
                theta_trial,
                phi_trial,
                grad_phi_step,
                alpha,
            );

            if acceptable {
                // Accept step
                state.x = x_trial;
                state.obj = obj_trial;
                state.g = g_trial;
                state.alpha_primal = alpha;
                step_accepted = true;

                // Add to filter if not using switching condition
                if !_used_switching {
                    filter.add(theta_current, phi_current);
                }
                break;
            }

            // Second-order correction (SOC)
            if _ls_iter == 0 && theta_trial > theta_current && options.max_soc > 0 {
                let soc_accepted = attempt_soc(
                    &state,
                    problem,
                    &x_trial,
                    &g_trial,
                    &mut lin_solver,
                    &kkt_system,
                    &filter,
                    theta_current,
                    phi_current,
                    grad_phi_step,
                    alpha,
                    options,
                );

                if let Some((x_soc, obj_soc, g_soc, alpha_soc)) = soc_accepted {
                    state.x = x_soc;
                    state.obj = obj_soc;
                    state.g = g_soc;
                    state.alpha_primal = alpha_soc;
                    step_accepted = true;
                    filter.add(theta_current, phi_current);
                    break;
                }
            }

            // Backtrack
            alpha *= 0.5;
        }

        if !step_accepted {
            // Try restoration phase
            log::debug!("Line search failed at iteration {}, entering restoration", iteration);

            let (x_rest, success) = restoration.restore(
                &state.x,
                &state.x_l,
                &state.x_u,
                &state.g_l,
                &state.g_u,
                &state.jac_rows,
                &state.jac_cols,
                n,
                m,
                options,
                &|theta, phi| filter.is_acceptable(theta, phi),
                &|x_eval, g_out| problem.constraints(x_eval, g_out),
                &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
            );

            if success {
                // Check if restoration actually moved the point
                let move_dist: f64 = x_rest.iter().zip(state.x.iter())
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0f64, f64::max);

                state.x = x_rest;
                state.alpha_primal = 0.0;
                // Re-evaluate at restoration point
                state.evaluate(problem, 1.0);

                if move_dist < 1e-14 && primal_inf < options.tol {
                    // Restoration returned same point AND constraints already satisfied.
                    // This happens for unconstrained problems or when the filter blocks
                    // progress on a fully feasible iterate (e.g., concave objectives).
                    consecutive_zero_alpha += 1;
                    if consecutive_zero_alpha >= stall_threshold {
                        stall_recovery_count += 1;
                        filter.reset();
                        inertia_params.delta_w_last = 0.0;
                        consecutive_zero_alpha = 0;

                        if stall_recovery_count >= 2 {
                            // Repeated stalls: filter reset alone isn't enough.
                            // Bump mu to change the barrier landscape.
                            let mu_new = (state.mu * 100.0).max(options.mu_init);
                            log::debug!(
                                "Repeated stall (recovery #{}), mu: {:.2e} -> {:.2e}",
                                stall_recovery_count, state.mu, mu_new
                            );
                            state.mu = mu_new;
                            stall_recovery_count = 0;
                        } else {
                            log::debug!("Stall detected, resetting filter (recovery #{})", stall_recovery_count);
                        }
                    }
                } else {
                    consecutive_zero_alpha = 0;
                }

                continue;
            } else {
                log::warn!("Restoration failed at iteration {}", iteration);
                return make_result(&state, SolveStatus::RestorationFailed);
            }
        }

        // Step was accepted — reset stall counter (any accepted step is progress)
        consecutive_zero_alpha = 0;

        // Update dual variables
        let alpha_d = alpha_dual_max;
        for i in 0..m {
            state.y[i] += alpha_d * state.dy[i];
        }
        // Ipopt kappa_sigma safeguard: keep z*s in [mu/kappa_sigma, kappa_sigma*mu]
        let kappa_sigma = 1e10;
        for i in 0..n {
            if state.x_l[i].is_finite() {
                let z_new = (state.z_l[i] + alpha_d * state.dz_l[i]).max(1e-20);
                let s_l = (state.x[i] - state.x_l[i]).max(1e-20);
                let z_lo = state.mu / (kappa_sigma * s_l);
                let z_hi = kappa_sigma * state.mu / s_l;
                state.z_l[i] = z_new.clamp(z_lo, z_hi);
            }
            if state.x_u[i].is_finite() {
                let z_new = (state.z_u[i] + alpha_d * state.dz_u[i]).max(1e-20);
                let s_u = (state.x_u[i] - state.x[i]).max(1e-20);
                let z_lo = state.mu / (kappa_sigma * s_u);
                let z_hi = kappa_sigma * state.mu / s_u;
                state.z_u[i] = z_new.clamp(z_lo, z_hi);
            }
        }

        state.alpha_dual = alpha_d;

        // Re-evaluate at new point
        state.evaluate(problem, 1.0);

        // Update barrier parameter
        state.mu = update_barrier_parameter(&state, options);
    }

    make_result(&state, SolveStatus::MaxIterations)
}

/// Attempt a second-order correction step.
///
/// If the trial point has worse constraint violation than the current point,
/// try to correct by solving a modified system targeting the trial constraint values.
#[allow(clippy::too_many_arguments)]
fn attempt_soc<P: NlpProblem>(
    state: &SolverState,
    problem: &P,
    _x_trial: &[f64],
    g_trial: &[f64],
    solver: &mut dyn LinearSolver,
    kkt: &kkt::KktSystem,
    filter: &Filter,
    theta_current: f64,
    phi_current: f64,
    grad_phi_step: f64,
    alpha: f64,
    options: &SolverOptions,
) -> Option<(Vec<f64>, f64, Vec<f64>, f64)> {
    let n = state.n;
    let m = state.m;

    if m == 0 {
        return None;
    }

    // Compute constraint residual at trial point, respecting constraint type
    let mut c_soc = vec![0.0; m];
    for i in 0..m {
        let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_equality || state.g_l[i].is_finite() {
            c_soc[i] = g_trial[i] - state.g_l[i];
        } else if state.g_u[i].is_finite() {
            c_soc[i] = g_trial[i] - state.g_u[i];
        }
    }

    for _soc_iter in 0..options.max_soc {
        // Modify RHS for SOC: replace primal residual with trial constraint residual
        let mut rhs_soc = kkt.rhs.clone();
        for i in 0..m {
            rhs_soc[n + i] = -c_soc[i];
        }

        // Solve with same factored matrix
        let mut sol_soc = vec![0.0; n + m];
        if solver.solve(&rhs_soc, &mut sol_soc).is_err() {
            return None;
        }

        let dx_soc = &sol_soc[..n];

        // Compute SOC trial point
        #[allow(clippy::needless_range_loop)]
        let mut x_soc = vec![0.0; n];
        for i in 0..n {
            x_soc[i] = state.x[i] + alpha * dx_soc[i];
            if state.x_l[i].is_finite() {
                x_soc[i] = x_soc[i].max(state.x_l[i] + 1e-14);
            }
            if state.x_u[i].is_finite() {
                x_soc[i] = x_soc[i].min(state.x_u[i] - 1e-14);
            }
        }

        let obj_soc = problem.objective(&x_soc);
        let mut g_soc = vec![0.0; m];
        problem.constraints(&x_soc, &mut g_soc);

        let theta_soc = convergence::primal_infeasibility(&g_soc, &state.g_l, &state.g_u);

        // Compute barrier objective
        let mut phi_soc = obj_soc;
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            if state.x_l[i].is_finite() {
                let slack = (x_soc[i] - state.x_l[i]).max(1e-20);
                phi_soc -= state.mu * slack.ln();
            }
            if state.x_u[i].is_finite() {
                let slack = (state.x_u[i] - x_soc[i]).max(1e-20);
                phi_soc -= state.mu * slack.ln();
            }
        }
        if options.constraint_slack_barrier {
            for i in 0..m {
                let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                    && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
                if is_eq {
                    continue;
                }
                if state.g_l[i].is_finite() {
                    let slack = g_soc[i] - state.g_l[i];
                    if slack > state.mu * 1e-2 {
                        phi_soc -= state.mu * slack.ln();
                    }
                }
                if state.g_u[i].is_finite() {
                    let slack = state.g_u[i] - g_soc[i];
                    if slack > state.mu * 1e-2 {
                        phi_soc -= state.mu * slack.ln();
                    }
                }
            }
        }

        let (acceptable, _) = filter.check_acceptability(
            theta_current,
            phi_current,
            theta_soc,
            phi_soc,
            grad_phi_step,
            alpha,
        );

        if acceptable {
            return Some((x_soc, obj_soc, g_soc, alpha));
        }

        // Update c_soc for next SOC iteration (respecting constraint type)
        for i in 0..m {
            let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_equality || state.g_l[i].is_finite() {
                c_soc[i] = g_soc[i] - state.g_l[i];
            } else if state.g_u[i].is_finite() {
                c_soc[i] = g_soc[i] - state.g_u[i];
            }
        }
    }

    None
}

/// Update the barrier parameter mu.
fn update_barrier_parameter(state: &SolverState, options: &SolverOptions) -> f64 {
    let mu = state.mu;

    if options.mu_strategy_adaptive {
        // Adaptive (Loqo) rule: mu = (x^T z) / (n * kappa)
        let mut sum_compl = 0.0;
        let mut count = 0;
        for i in 0..state.n {
            if state.x_l[i].is_finite() {
                let slack = (state.x[i] - state.x_l[i]).max(1e-20);
                sum_compl += slack * state.z_l[i];
                count += 1;
            }
            if state.x_u[i].is_finite() {
                let slack = (state.x_u[i] - state.x[i]).max(1e-20);
                sum_compl += slack * state.z_u[i];
                count += 1;
            }
        }

        if count == 0 {
            return (mu * options.mu_linear_decrease_factor).max(options.mu_min);
        }

        let avg_compl = sum_compl / count as f64;

        // Kappa adaptive: reduce faster when progress is good
        let mu_new = (avg_compl / options.kappa).max(options.mu_min);

        // Allow mu to increase if the option is set (helps after restoration/stall)
        if options.mu_allow_increase {
            mu_new
        } else {
            mu_new.min(mu)
        }
    } else {
        // Monotone decrease
        let mu_new = mu * options.mu_linear_decrease_factor;
        mu_new.max(options.mu_min)
    }
}

/// Build the final solve result.
/// Computes z from stationarity for more accurate output multipliers.
fn make_result(state: &SolverState, status: SolveStatus) -> SolveResult {
    // Compute optimal z from stationarity: ∇f + J^T y - z_l + z_u = 0
    let n = state.n;
    let mut grad_jty = state.grad_f.clone();
    for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        grad_jty[col] += state.jac_vals[idx] * state.y[row];
    }
    let mut z_l_out = vec![0.0; n];
    let mut z_u_out = vec![0.0; n];
    for i in 0..n {
        if grad_jty[i] > 0.0 && state.x_l[i].is_finite() {
            z_l_out[i] = grad_jty[i];
        } else if grad_jty[i] < 0.0 && state.x_u[i].is_finite() {
            z_u_out[i] = -grad_jty[i];
        }
    }

    SolveResult {
        x: state.x.clone(),
        objective: state.obj,
        constraint_multipliers: state.y.clone(),
        bound_multipliers_lower: z_l_out,
        bound_multipliers_upper: z_u_out,
        constraint_values: state.g.clone(),
        status,
        iterations: state.iter,
    }
}
