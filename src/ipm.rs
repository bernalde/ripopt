use std::time::Instant;

use crate::convergence::{self, check_convergence, ConvergenceInfo, ConvergenceStatus};
use crate::filter::{self, Filter, FilterEntry};
use crate::kkt::{self, InertiaCorrectionParams};
use crate::linear_solver::dense::DenseLdl;
use crate::linear_solver::{LinearSolver, SymmetricMatrix};
use crate::options::SolverOptions;
use crate::problem::NlpProblem;
use crate::restoration::RestorationPhase;
use crate::result::{SolveResult, SolveStatus};
use crate::warmstart::WarmStartInitializer;

/// NLP problem wrapper that applies gradient-based scaling.
///
/// Scales objective by `obj_scaling` and each constraint `i` by `g_scaling[i]`
/// so that the max gradient norm at the initial point is ≤ 100.
/// This matches Ipopt's `nlp_scaling_method = gradient-based`.
struct ScaledProblem<'a, P: NlpProblem> {
    inner: &'a P,
    obj_scaling: f64,
    g_scaling: Vec<f64>,
    jac_rows: Vec<usize>,
}

impl<P: NlpProblem> NlpProblem for ScaledProblem<'_, P> {
    fn num_variables(&self) -> usize {
        self.inner.num_variables()
    }
    fn num_constraints(&self) -> usize {
        self.inner.num_constraints()
    }
    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        self.inner.bounds(x_l, x_u);
    }
    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        self.inner.constraint_bounds(g_l, g_u);
        for (i, &s) in self.g_scaling.iter().enumerate() {
            if g_l[i].is_finite() {
                g_l[i] *= s;
            }
            if g_u[i].is_finite() {
                g_u[i] *= s;
            }
        }
    }
    fn initial_point(&self, x0: &mut [f64]) {
        self.inner.initial_point(x0);
    }
    fn objective(&self, x: &[f64]) -> f64 {
        self.inner.objective(x) * self.obj_scaling
    }
    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        self.inner.gradient(x, grad);
        for g in grad.iter_mut() {
            *g *= self.obj_scaling;
        }
    }
    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        self.inner.constraints(x, g);
        for (i, &s) in self.g_scaling.iter().enumerate() {
            g[i] *= s;
        }
    }
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        self.inner.jacobian_structure()
    }
    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        self.inner.jacobian_values(x, vals);
        for (idx, &row) in self.jac_rows.iter().enumerate() {
            vals[idx] *= self.g_scaling[row];
        }
    }
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        self.inner.hessian_structure()
    }
    fn hessian_values(&self, x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]) {
        let scaled_lambda: Vec<f64> = lambda
            .iter()
            .zip(self.g_scaling.iter())
            .map(|(l, s)| l * s)
            .collect();
        self.inner
            .hessian_values(x, obj_factor * self.obj_scaling, &scaled_lambda, vals);
    }
}

/// Saved state for the watchdog mechanism.
struct WatchdogSavedState {
    x: Vec<f64>,
    y: Vec<f64>,
    z_l: Vec<f64>,
    z_u: Vec<f64>,
    mu: f64,
    obj: f64,
    g: Vec<f64>,
    grad_f: Vec<f64>,
    filter_entries: Vec<FilterEntry>,
    theta: f64,
    phi: f64,
}

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
    /// Objective scaling factor (for NLP scaling / result unscaling).
    pub obj_scaling: f64,
    /// Constraint scaling factors (for NLP scaling / result unscaling).
    pub g_scaling: Vec<f64>,
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
                    // For inequality constraints, zero out multipliers with wrong sign.
                    // Ipopt convention (L = f + y^T g):
                    //   g >= g_l (lower bound only): y >= 0
                    //   g <= g_u (upper bound only): y <= 0
                    //   g_l <= g <= g_u: sign depends on active bound
                    // Only apply LS init for equality constraints or when sign is correct.
                    for i in 0..m {
                        let is_eq = g_l[i].is_finite() && g_u[i].is_finite()
                            && (g_l[i] - g_u[i]).abs() < 1e-15;
                        if is_eq {
                            continue; // LS init fine for equalities
                        }
                        let has_lower = g_l[i].is_finite();
                        let has_upper = g_u[i].is_finite();
                        if has_lower && !has_upper && y_ls[i] < 0.0 {
                            y_ls[i] = 0.0; // Wrong sign for lower-bound constraint
                        } else if has_upper && !has_lower && y_ls[i] > 0.0 {
                            y_ls[i] = 0.0; // Wrong sign for upper-bound constraint
                        } else if !has_lower && !has_upper {
                            y_ls[i] = 0.0; // No bounds at all (shouldn't happen)
                        }
                    }
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
            obj_scaling: 1.0,
            g_scaling: vec![1.0; m],
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

/// Wrapper that reformulates an overdetermined nonlinear equations (NE) problem
/// as an unconstrained least-squares problem.
///
/// Original:  min 0  subject to  g_i(x) = target_i,  i = 1,...,m   (m > n)
/// Reformulated:  min 0.5 * Σ (g_i(x) - target_i)^2   (no constraints, keep variable bounds)
///
/// Gradient: ∇f_LS = J^T * r   where r_i = g_i(x) - target_i
/// Hessian:  H_LS ≈ J^T * J   (Gauss-Newton approximation)
struct LeastSquaresProblem<'a, P: NlpProblem> {
    inner: &'a P,
    /// Targets for each constraint (g_l == g_u for equalities).
    targets: Vec<f64>,
    /// Jacobian structure cached from inner problem.
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    /// Hessian structure: lower triangle of J^T*J.
    hess_rows: Vec<usize>,
    hess_cols: Vec<usize>,
}

impl<P: NlpProblem> LeastSquaresProblem<'_, P> {
    fn new(inner: &P) -> LeastSquaresProblem<'_, P> {
        let n = inner.num_variables();
        let m = inner.num_constraints();
        let mut g_l = vec![0.0; m];
        let mut g_u = vec![0.0; m];
        inner.constraint_bounds(&mut g_l, &mut g_u);
        let targets: Vec<f64> = (0..m).map(|i| 0.5 * (g_l[i] + g_u[i])).collect();

        let (jac_rows, jac_cols) = inner.jacobian_structure();

        // Build Hessian structure for J^T*J (lower triangle, n x n dense).
        // Since J^T*J is generally dense, use full lower triangle.
        let mut hess_rows = Vec::with_capacity(n * (n + 1) / 2);
        let mut hess_cols = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in 0..=i {
                hess_rows.push(i);
                hess_cols.push(j);
            }
        }

        LeastSquaresProblem {
            inner,
            targets,
            jac_rows,
            jac_cols,
            hess_rows,
            hess_cols,
        }
    }
}

impl<P: NlpProblem> NlpProblem for LeastSquaresProblem<'_, P> {
    fn num_variables(&self) -> usize {
        self.inner.num_variables()
    }
    fn num_constraints(&self) -> usize {
        0 // unconstrained LS
    }
    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        self.inner.bounds(x_l, x_u);
    }
    fn constraint_bounds(&self, _g_l: &mut [f64], _g_u: &mut [f64]) {
        // no constraints
    }
    fn initial_point(&self, x0: &mut [f64]) {
        self.inner.initial_point(x0);
    }
    fn objective(&self, x: &[f64]) -> f64 {
        let m = self.targets.len();
        let mut g = vec![0.0; m];
        self.inner.constraints(x, &mut g);
        let mut sum = 0.0;
        for i in 0..m {
            let r = g[i] - self.targets[i];
            sum += r * r;
        }
        0.5 * sum
    }
    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        let n = self.inner.num_variables();
        let m = self.targets.len();
        // Compute residual r = g(x) - target
        let mut g = vec![0.0; m];
        self.inner.constraints(x, &mut g);
        let mut r = vec![0.0; m];
        for i in 0..m {
            r[i] = g[i] - self.targets[i];
        }
        // grad = J^T * r
        let jac_nnz = self.jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];
        self.inner.jacobian_values(x, &mut jac_vals);
        for i in 0..n {
            grad[i] = 0.0;
        }
        for (idx, (&row, &col)) in self.jac_rows.iter().zip(self.jac_cols.iter()).enumerate() {
            grad[col] += jac_vals[idx] * r[row];
        }
    }
    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {
        // no constraints
    }
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![]) // no constraints
    }
    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {
        // no constraints
    }
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hess_rows.clone(), self.hess_cols.clone())
    }
    fn hessian_values(&self, x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        let n = self.inner.num_variables();
        // H ≈ obj_factor * J^T * J  (Gauss-Newton approximation)
        let jac_nnz = self.jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];
        self.inner.jacobian_values(x, &mut jac_vals);

        // Build dense J (m x n) for small problems, then compute J^T*J
        let m = self.targets.len();
        let mut j_dense = vec![0.0; m * n];
        for (idx, (&row, &col)) in self.jac_rows.iter().zip(self.jac_cols.iter()).enumerate() {
            j_dense[row * n + col] += jac_vals[idx];
        }

        // Compute lower triangle of J^T*J
        let mut idx = 0;
        for i in 0..n {
            for j in 0..=i {
                let mut dot = 0.0;
                for k in 0..m {
                    dot += j_dense[k * n + i] * j_dense[k * n + j];
                }
                vals[idx] = obj_factor * dot;
                idx += 1;
            }
        }
    }
}

/// Detect if a problem is an overdetermined nonlinear equation system.
///
/// Returns true if ALL of:
/// - f(x0) ≈ 0 (zero objective)
/// - ∇f(x0) ≈ 0 (zero gradient)
/// - All constraints are equalities (g_l[i] == g_u[i])
/// - m > n (more constraints than variables)
fn detect_ne_problem<P: NlpProblem>(problem: &P) -> bool {
    let n = problem.num_variables();
    let m = problem.num_constraints();

    if m <= n || m == 0 || n == 0 {
        return false;
    }

    // Check objective and gradient at initial point
    let mut x0 = vec![0.0; n];
    problem.initial_point(&mut x0);

    let f0 = problem.objective(&x0);
    if f0.abs() > 1e-10 {
        return false;
    }

    let mut grad = vec![0.0; n];
    problem.gradient(&x0, &mut grad);
    let grad_max = grad.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    if grad_max > 1e-10 {
        return false;
    }

    // Check all constraints are equalities
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    for i in 0..m {
        let is_eq = g_l[i].is_finite()
            && g_u[i].is_finite()
            && (g_l[i] - g_u[i]).abs() < 1e-15;
        if !is_eq {
            return false;
        }
    }

    // If constraints are already satisfied at x0, no need to reformulate —
    // the standard IPM can handle it (e.g., FBRAIN3 starts feasible).
    let mut g0 = vec![0.0; m];
    problem.constraints(&x0, &mut g0);
    let theta0 = convergence::primal_infeasibility(&g0, &g_l, &g_u);
    if theta0 < 1e-8 {
        return false;
    }

    true
}

/// Solve the NLP using the interior point method.
pub fn solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    // --- NE-to-LS Detection and Reformulation ---
    // Detect overdetermined nonlinear equation problems (m > n, f≡0, all equalities)
    // and reformulate as least-squares: min 0.5*||g(x)-target||^2.
    if detect_ne_problem(problem) {
        let n = problem.num_variables();
        let m = problem.num_constraints();
        if options.print_level >= 5 {
            eprintln!(
                "ripopt: Detected overdetermined NE problem (n={}, m={}), reformulating as least-squares",
                n, m
            );
        }
        let ls_problem = LeastSquaresProblem::new(problem);
        let ls_result = solve_ipm(&ls_problem, options);

        // Evaluate original constraint violation at the LS solution
        let mut g_final = vec![0.0; m];
        problem.constraints(&ls_result.x, &mut g_final);
        let mut g_l = vec![0.0; m];
        let mut g_u = vec![0.0; m];
        problem.constraint_bounds(&mut g_l, &mut g_u);
        let theta = convergence::primal_infeasibility(&g_final, &g_l, &g_u);

        // g_final is from the original (unscaled) problem, no unscaling needed
        let g_out = g_final;

        let status = if theta < options.tol {
            SolveStatus::Optimal
        } else if theta < options.acceptable_tol {
            SolveStatus::Acceptable
        } else {
            SolveStatus::LocalInfeasibility
        };

        if options.print_level >= 5 {
            eprintln!(
                "ripopt: NE-to-LS result: obj_LS={:.4e}, constraint_violation={:.4e}, status={:?}",
                ls_result.objective, theta, status
            );
        }

        return SolveResult {
            x: ls_result.x,
            objective: 0.0, // Original objective is f≡0
            constraint_multipliers: vec![0.0; m],
            bound_multipliers_lower: ls_result.bound_multipliers_lower,
            bound_multipliers_upper: ls_result.bound_multipliers_upper,
            constraint_values: g_out,
            status,
            iterations: ls_result.iterations,
        };
    }

    solve_ipm(problem, options)
}

/// Core IPM solver implementation.
fn solve_ipm<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    // --- NLP Scaling (gradient-based, matching Ipopt's nlp_scaling_method) ---
    // Scale objective and constraints so max gradient norm at x0 is ≤ 100.
    let n_sc = problem.num_variables();
    let m_sc = problem.num_constraints();

    let mut x0 = vec![0.0; n_sc];
    problem.initial_point(&mut x0);

    // Objective scaling: obj_scaling = min(1, max_grad / ||∇f(x0)||∞)
    // Use a floor of 1e-2 to prevent extreme downscaling that creates a mismatch
    // with the unscaled log-barrier terms (we don't have explicit slack variables).
    let nlp_scaling_max_gradient = 100.0;
    let nlp_scaling_min_value = 1e-2;
    let mut grad_f0 = vec![0.0; n_sc];
    problem.gradient(&x0, &mut grad_f0);
    let grad_max = grad_f0
        .iter()
        .map(|v| v.abs())
        .fold(0.0f64, f64::max);
    // Skip objective scaling for unconstrained problems: there is no constraint-objective
    // tradeoff to balance, and extreme downscaling makes the barrier terms dominate.
    let obj_scaling = if m_sc > 0 && grad_max > nlp_scaling_max_gradient && grad_max.is_finite() {
        (nlp_scaling_max_gradient / grad_max).max(nlp_scaling_min_value)
    } else {
        1.0
    };

    // Constraint scaling: g_scaling[i] = min(1, 100 / ||∇g_i(x0)||∞)
    // Skip constraint scaling when initial violation is very large (>1e6) — extreme
    // scaling ratios make restoration ill-conditioned for far-from-feasible problems.
    let (jac_rows_sc, _) = problem.jacobian_structure();
    let mut g_scaling = vec![1.0; m_sc];
    if m_sc > 0 {
        let mut g0_sc = vec![0.0; m_sc];
        problem.constraints(&x0, &mut g0_sc);
        let mut g_l_sc = vec![0.0; m_sc];
        let mut g_u_sc = vec![0.0; m_sc];
        problem.constraint_bounds(&mut g_l_sc, &mut g_u_sc);
        let init_cv = convergence::primal_infeasibility(&g0_sc, &g_l_sc, &g_u_sc);

        if init_cv < 1e6 {
            let mut jac_vals0 = vec![0.0; jac_rows_sc.len()];
            problem.jacobian_values(&x0, &mut jac_vals0);
            let mut row_max = vec![0.0f64; m_sc];
            for (idx, &row) in jac_rows_sc.iter().enumerate() {
                let v = jac_vals0[idx].abs();
                if v.is_finite() && v > row_max[row] {
                    row_max[row] = v;
                }
            }
            for i in 0..m_sc {
                if row_max[i] > nlp_scaling_max_gradient {
                    g_scaling[i] = (nlp_scaling_max_gradient / row_max[i]).max(nlp_scaling_min_value);
                }
            }
        }
    }

    if options.print_level >= 5
        && (obj_scaling != 1.0 || g_scaling.iter().any(|&s| s != 1.0))
    {
        let n_scaled_g = g_scaling.iter().filter(|&&s| s != 1.0).count();
        eprintln!(
            "ripopt: NLP scaling: obj_scaling={:.4e}, {}/{} constraints scaled",
            obj_scaling, n_scaled_g, m_sc
        );
    }

    let scaled = ScaledProblem {
        inner: problem,
        obj_scaling,
        g_scaling: g_scaling.clone(),
        jac_rows: jac_rows_sc,
    };
    let problem = &scaled; // shadow: all subsequent code uses the scaled problem

    let mut state = SolverState::new(problem, options);
    state.obj_scaling = obj_scaling;
    state.g_scaling = g_scaling;
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
    let mut restoration = RestorationPhase::new(500);

    // Initialize filter
    let mut filter = Filter::new(1e4);

    // Stall detection: count consecutive iterations with alpha_primal ≈ 0.
    // Progressive recovery: first reset filter only, then bump mu on repeated stalls.
    let mut consecutive_zero_alpha: usize = 0;
    let stall_threshold: usize = 5;
    let mut stall_recovery_count: usize = 0;

    // Unified stall recovery counter shared across all stall detection phases.
    // Cap total recoveries at 6 to prevent infinite cycling.
    let mut total_stall_recovery_count: usize = 0;
    let mu_perturbation_factors: [f64; 6] = [10.0, 0.1, 100.0, 0.01, 1000.0, 0.001];

    // Wall-clock time limit
    let start_time = Instant::now();

    // Watchdog mechanism state
    let mut consecutive_shortened: usize = 0;
    let mut watchdog_active: bool = false;
    let mut watchdog_trial_count: usize = 0;
    let mut watchdog_saved: Option<WatchdogSavedState> = None;

    // Alpha history for stall detection (Phase 3)
    let alpha_history_len: usize = 20;
    let mut alpha_history: Vec<f64> = Vec::with_capacity(alpha_history_len);

    // Constraint violation history for infeasibility detection
    let theta_history_len: usize = 100;
    let mut theta_history: Vec<f64> = Vec::with_capacity(theta_history_len);

    // Track whether the problem was ever feasible (theta < constr_viol_tol)
    // to prevent false infeasibility declarations on feasible problems.
    let mut ever_feasible = false;

    // Ipopt-like monotone mu decrease fallback: after 20 stuck iterations,
    // decrease mu by min(0.2*mu, mu^1.5). Track a ceiling to prevent stall
    // recovery from undoing forced decreases.
    let mut mu_stuck_count: usize = 0;
    let mut mu_ceiling: f64 = 1e4;
    // Tiny step counter (Ipopt: accept full step when relative step < 10*eps for 2 consecutive)
    let mut consecutive_tiny_steps: usize = 0;

    // Consecutive iterations with obj < -1e20 for robust unbounded detection
    let mut consecutive_unbounded: usize = 0;
    // Fix 3B: Track iterations since last stall recovery mu bump
    let mut iters_since_stall_recovery: usize = usize::MAX;

    // Best feasible point tracking: save the best (lowest obj) point that is feasible
    let mut best_x: Option<Vec<f64>> = None;
    let mut best_obj: f64 = f64::INFINITY;
    let mut best_y: Option<Vec<f64>> = None;
    let mut best_z_l: Option<Vec<f64>> = None;
    let mut best_z_u: Option<Vec<f64>> = None;

    // Initial evaluation
    state.evaluate(problem, 1.0);

    // NaN/Inf guard on initial evaluation — try perturbation before giving up
    if state.obj.is_nan() || state.obj.is_infinite()
        || state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
    {
        let mut recovered = false;
        let x_saved = state.x.clone();
        for &push_factor in &[1e-2, 1e-1, 0.5] {
            // Reset to saved point and apply stronger push
            state.x.copy_from_slice(&x_saved);
            for i in 0..n {
                if state.x_l[i].is_finite() && state.x_u[i].is_finite() {
                    let range = state.x_u[i] - state.x_l[i];
                    let push = push_factor * range;
                    if range > 2.0 * push {
                        state.x[i] = state.x[i].max(state.x_l[i] + push).min(state.x_u[i] - push);
                    } else {
                        state.x[i] = 0.5 * (state.x_l[i] + state.x_u[i]);
                    }
                } else if state.x_l[i].is_finite() {
                    let push = push_factor * state.x_l[i].abs().max(1.0);
                    state.x[i] = state.x[i].max(state.x_l[i] + push);
                } else if state.x_u[i].is_finite() {
                    let push = push_factor * state.x_u[i].abs().max(1.0);
                    state.x[i] = state.x[i].min(state.x_u[i] - push);
                }
            }
            // Re-initialize bound multipliers after perturbation
            for i in 0..n {
                if state.x_l[i].is_finite() {
                    let slack = (state.x[i] - state.x_l[i]).max(1e-20);
                    state.z_l[i] = options.mu_init / slack;
                }
                if state.x_u[i].is_finite() {
                    let slack = (state.x_u[i] - state.x[i]).max(1e-20);
                    state.z_u[i] = options.mu_init / slack;
                }
            }
            state.evaluate(problem, 1.0);
            if !state.obj.is_nan() && !state.obj.is_infinite()
                && !state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
            {
                recovered = true;
                break;
            }
        }
        if !recovered {
            return make_result(&state, SolveStatus::NumericalError);
        }
    }

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

        // Check wall-clock time limit (every 10 iterations to avoid syscall overhead)
        if iteration % 10 == 0 && options.max_wall_time > 0.0 {
            if start_time.elapsed().as_secs_f64() >= options.max_wall_time {
                return make_result(&state, SolveStatus::MaxIterations);
            }
        }

        // Compute optimality measures.
        let primal_inf = state.constraint_violation();

        // Compute z_opt from stationarity for the scaled convergence check.
        // At optimality, grad_f + J^T y - z_l + z_u = 0.
        // z_opt captures the true bound multiplier for active bounds.
        //
        // Complementarity gate: only use z_opt when z_opt * slack is consistent
        // with the barrier problem (z*s ~ mu). If z_opt * slack >> mu, the point
        // is not a barrier-optimal point and z_opt would hide a true infeasibility.
        let (z_l_opt, z_u_opt) = {
            let mut grad_jty = state.grad_f.clone();
            for (idx, (&row, &col)) in
                state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
            {
                grad_jty[col] += state.jac_vals[idx] * state.y[row];
            }
            let mut zl = vec![0.0; n];
            let mut zu = vec![0.0; n];
            let kappa_compl = 1e10;
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

        // Scaled dual infeasibility uses z_opt (for fast convergence detection)
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

        // Unscaled dual infeasibility uses iterative z (catches false convergence)
        let dual_inf_unscaled = convergence::dual_infeasibility(
            &state.grad_f,
            &state.jac_rows,
            &state.jac_cols,
            &state.jac_vals,
            &state.y,
            &state.z_l,
            &state.z_u,
            n,
        );
        let compl_inf = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
        );
        // Also compute complementarity using z_opt (NLP multipliers from stationarity).
        // When mu is stuck high, kappa_sigma safeguard inflates iterative z, making
        // compl_inf huge even at the NLP optimum. z_opt correctly reflects the NLP solution.
        let compl_inf_opt = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &z_l_opt, &z_u_opt, 0.0,
        );
        let compl_inf_best = compl_inf.min(compl_inf_opt);

        if options.print_level >= 5 {
            eprintln!(
                "{:>4} obj={:>14.7e} pr={:>10.2e} du={:>10.2e} co={:>10.2e} mu={:>10.2e} a_p={:>8.2e} a_d={:>8.2e}",
                iteration,
                state.obj / state.obj_scaling,
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

        // Check convergence (use best complementarity: min of iterative z and z_opt)
        let conv_info = ConvergenceInfo {
            primal_inf,
            dual_inf,
            dual_inf_unscaled,
            compl_inf: compl_inf_best,
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

        // Track consecutive acceptable iterations (using same criteria as check_convergence)
        let s_d_for_acc = {
            let s_max: f64 = 100.0;
            let s_d_max: f64 = 1e4;
            if (m + 2 * n) > 0 {
                ((s_max.max(multiplier_sum / (m + 2 * n) as f64)) / s_max).min(s_d_max)
            } else {
                1.0
            }
        };
        let meets_acc_scaled = primal_inf <= options.acceptable_tol
            && dual_inf <= options.acceptable_tol * s_d_for_acc
            && compl_inf_best <= options.acceptable_tol * s_d_for_acc;
        let meets_acc_unscaled = primal_inf <= options.acceptable_constr_viol_tol
            && dual_inf_unscaled <= options.acceptable_dual_inf_tol
            && compl_inf_best <= options.acceptable_compl_inf_tol;
        if meets_acc_scaled && meets_acc_unscaled {
            state.consecutive_acceptable += 1;
        } else {
            state.consecutive_acceptable = 0;
        }

        // Track constraint violation history for infeasibility detection
        if theta_history.len() >= theta_history_len {
            theta_history.remove(0);
        }
        theta_history.push(primal_inf);

        // Track whether we've ever been feasible
        if primal_inf < options.constr_viol_tol {
            ever_feasible = true;
        }

        // Unbounded detection: objective diverging negatively with satisfied constraints.
        // Require 10 consecutive iterations to avoid false positives from transient dips.
        if state.obj < -1e20 && primal_inf < options.constr_viol_tol {
            consecutive_unbounded += 1;
            if consecutive_unbounded >= 10 {
                return make_result(&state, SolveStatus::Unbounded);
            }
        } else {
            consecutive_unbounded = 0;
        }

        // Compute sigma (barrier diagonal)
        let sigma = kkt::compute_sigma(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u);

        // Use condensed KKT (Schur complement) when m >= 2n for efficiency.
        // Condensed cost is O(n^2*m + n^3) vs O((n+m)^3) — strictly better when m > n.
        // Use m >= 2n threshold (relaxed from m > 2n) to avoid borderline cases
        // where condensed path has different numerical behavior.
        let use_condensed = m >= 2 * n && n > 0;

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

        let condensed_system = if use_condensed {
            Some(kkt::assemble_condensed_kkt(
                n, m,
                &state.hess_rows, &state.hess_cols, &state.hess_vals,
                &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                &state.y, &state.z_l, &state.z_u,
                &state.x, &state.x_l, &state.x_u, state.mu,
            ))
        } else {
            None
        };

        // Factor with inertia correction (full KKT — needed for SOC and fallback)
        let inertia_result =
            kkt::factor_with_inertia_correction(&mut kkt_system, &mut lin_solver, &mut inertia_params);

        if let Err(e) = inertia_result {
            log::warn!("KKT factorization failed: {}", e);

            // Early-iteration perturbation: if factorization fails in the first 5
            // iterations, the starting point is likely degenerate (singular Jacobian).
            // Try more aggressive perturbation scales before other recovery methods.
            if iteration < 5 {
                let mut early_recovered = false;
                for &perturb_scale in &[1e-4, 1e-3, 1e-2, 5e-2, 1e-1] {
                    let x_saved = state.x.clone();
                    for i in 0..n {
                        let mag = state.x[i].abs().max(1.0);
                        // Deterministic pseudo-random sign based on index and attempt
                        let sign = if (i * 7 + iteration * 13 + (perturb_scale * 1e4) as usize * 3) % 3 == 0 {
                            -1.0
                        } else {
                            1.0
                        };
                        state.x[i] += sign * perturb_scale * mag;
                        if state.x_l[i].is_finite() {
                            state.x[i] = state.x[i].max(state.x_l[i] + 1e-14);
                        }
                        if state.x_u[i].is_finite() {
                            state.x[i] = state.x[i].min(state.x_u[i] - 1e-14);
                        }
                    }
                    // Re-initialize bound multipliers after perturbation
                    for i in 0..n {
                        if state.x_l[i].is_finite() {
                            let slack = (state.x[i] - state.x_l[i]).max(1e-20);
                            state.z_l[i] = state.mu / slack;
                        }
                        if state.x_u[i].is_finite() {
                            let slack = (state.x_u[i] - state.x[i]).max(1e-20);
                            state.z_u[i] = state.mu / slack;
                        }
                    }
                    state.evaluate(problem, 1.0);
                    if !state.obj.is_nan() && !state.obj.is_infinite()
                        && !state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
                    {
                        // Re-try factorization at perturbed point
                        let sigma_p = kkt::compute_sigma(
                            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u,
                        );
                        let mut kkt_p = kkt::assemble_kkt(
                            n, m, &state.hess_rows, &state.hess_cols, &state.hess_vals,
                            &state.jac_rows, &state.jac_cols, &state.jac_vals, &sigma_p,
                            &state.grad_f, &state.g, &state.g_l, &state.g_u,
                            &state.y, &state.z_l, &state.z_u,
                            &state.x, &state.x_l, &state.x_u, state.mu,
                        );
                        if kkt::factor_with_inertia_correction(
                            &mut kkt_p, &mut lin_solver, &mut inertia_params,
                        ).is_ok() {
                            log::debug!(
                                "Early perturbation (scale={:.0e}) recovered factorization at iter {}",
                                perturb_scale, iteration
                            );
                            filter.reset();
                            let theta_p = state.constraint_violation();
                            filter.set_theta_min_from_initial(theta_p);
                            early_recovered = true;
                            break;
                        }
                    }
                    // Restore if this perturbation didn't help
                    state.x.copy_from_slice(&x_saved);
                }
                if early_recovered {
                    continue;
                }
            }

            // Try gradient descent fallback
            if let Some(fallback) = gradient_descent_fallback(&state) {
                state.dx = fallback.0;
                state.dy = fallback.1;
                state.dz_l = vec![0.0; n];
                state.dz_u = vec![0.0; n];

                // Simple Armijo backtracking with the gradient step
                let mut alpha_fb = 1.0;
                let obj_current = state.obj;
                let mut fb_accepted = false;
                for _ in 0..20 {
                    let mut x_trial = vec![0.0; n];
                    for i in 0..n {
                        x_trial[i] = state.x[i] + alpha_fb * state.dx[i];
                        if state.x_l[i].is_finite() {
                            x_trial[i] = x_trial[i].max(state.x_l[i] + 1e-14);
                        }
                        if state.x_u[i].is_finite() {
                            x_trial[i] = x_trial[i].min(state.x_u[i] - 1e-14);
                        }
                    }
                    let obj_trial = problem.objective(&x_trial);
                    if !obj_trial.is_nan() && obj_trial < obj_current {
                        state.x = x_trial;
                        state.obj = obj_trial;
                        state.alpha_primal = alpha_fb;
                        fb_accepted = true;
                        break;
                    }
                    alpha_fb *= 0.5;
                }
                if fb_accepted {
                    state.evaluate(problem, 1.0);
                    continue;
                }
            }
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
            // Last resort: perturb x and retry factorization
            let mut recovered_from_perturb = false;
            for &perturb_scale in &[1e-3, 1e-2, 1e-1] {
                for i in 0..n {
                    let mag = state.x[i].abs().max(1.0);
                    let sign = if (i * 7 + iteration * 13) % 3 == 0 { -1.0 } else { 1.0 };
                    state.x[i] += sign * perturb_scale * mag;
                    if state.x_l[i].is_finite() {
                        state.x[i] = state.x[i].max(state.x_l[i] + 1e-14);
                    }
                    if state.x_u[i].is_finite() {
                        state.x[i] = state.x[i].min(state.x_u[i] - 1e-14);
                    }
                }
                state.evaluate(problem, 1.0);
                if !state.obj.is_nan() && !state.obj.is_infinite() {
                    recovered_from_perturb = true;
                    break;
                }
            }
            if recovered_from_perturb {
                filter.reset();
                let theta_new = state.constraint_violation();
                filter.set_theta_min_from_initial(theta_new);
                continue;
            }
            return make_result(&state, SolveStatus::NumericalError);
        }

        // Solve for search direction
        let (dx, dy) = if let Some(ref cond) = condensed_system {
            // Try condensed solve first (faster for m >> n)
            let mut cond_solver = DenseLdl::new();
            match cond_solver.factor(&cond.matrix) {
                Ok(_) => match kkt::solve_condensed(cond, &mut cond_solver) {
                    Ok(d) => d,
                    Err(_) => {
                        // Fall back to full KKT solve
                        match kkt::solve_for_direction(&kkt_system, &mut lin_solver) {
                            Ok(d) => d,
                            Err(e) => {
                                log::warn!("KKT solve failed: {}", e);
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
                        }
                    }
                },
                Err(_) => {
                    // Fall back to full KKT solve
                    match kkt::solve_for_direction(&kkt_system, &mut lin_solver) {
                        Ok(d) => d,
                        Err(e) => {
                            log::warn!("KKT solve failed: {}", e);
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
                    }
                }
            }
        } else {
            let dir_result = kkt::solve_for_direction(&kkt_system, &mut lin_solver);
            match dir_result {
                Ok(d) => d,
                Err(e) => {
                    log::warn!("KKT solve failed: {}", e);
                    // Try gradient descent fallback before restoration
                    if let Some(fallback) = gradient_descent_fallback(&state) {
                        fallback
                    } else {
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
                }
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

        // Ipopt-like tiny step detection: if relative step size is < 10*eps for
        // 2 consecutive iterations, force mu decrease and accept the full step.
        // This prevents stalling near the solution where numerical noise prevents
        // the line search from making progress.
        {
            let max_rel_step: f64 = (0..n)
                .map(|i| (alpha_primal_max * state.dx[i]).abs() / (state.x[i].abs() + 1.0))
                .fold(0.0f64, f64::max);
            if max_rel_step < 1e-14 && primal_inf < 1e-4 {
                consecutive_tiny_steps += 1;
                if consecutive_tiny_steps >= 2 {
                    // Force mu decrease (Ipopt: monotone decrease on tiny step)
                    let new_mu = (options.mu_linear_decrease_factor * state.mu)
                        .min(state.mu.powf(options.mu_superlinear_decrease_power))
                        .max(options.mu_min);
                    if (new_mu - state.mu).abs() < 1e-20 {
                        // mu is already at minimum and step is tiny: declare acceptable if possible
                        log::debug!("Tiny step with mu at minimum, checking acceptability");
                    } else {
                        state.mu = new_mu;
                        if ever_feasible {
                            mu_ceiling = state.mu;
                        }
                        filter.reset();
                        let theta_new = state.constraint_violation();
                        filter.set_theta_min_from_initial(theta_new);
                        log::debug!("Tiny step detected, forced mu decrease to {:.2e}", state.mu);
                    }
                    consecutive_tiny_steps = 0;
                }
            } else {
                consecutive_tiny_steps = 0;
            }
        }

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

            // NaN guard: reject trial points with NaN/Inf values
            if obj_trial.is_nan() || obj_trial.is_infinite()
                || g_trial.iter().any(|v| v.is_nan() || v.is_infinite())
            {
                alpha *= 0.5;
                continue;
            }

            let theta_trial =
                convergence::primal_infeasibility(&g_trial, &state.g_l, &state.g_u);

            // Watchdog mode: accept full step unconditionally (bypass filter)
            if watchdog_active && alpha == alpha_primal_max {
                state.x = x_trial;
                state.obj = obj_trial;
                state.g = g_trial;
                state.alpha_primal = alpha;
                step_accepted = true;
                break;
            }

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

            // Second-order correction (SOC) — try at every backtracking step where theta increased
            if theta_trial > theta_current && options.max_soc > 0 {
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

                // Reset filter after successful restoration.
                // TODO: Ipopt augments filter before restoration and doesn't reset after,
                // but that requires filter-aware restoration (our restoration ignores the filter).
                filter.reset();
                let theta_restored = state.constraint_violation();
                filter.set_theta_min_from_initial(theta_restored);

                // Reset acceptable counter — restoration invalidates previous acceptable streak
                state.consecutive_acceptable = 0;

                // Recompute mu from current complementarity after restoration
                let mu_compl = compute_avg_complementarity(&state);
                if mu_compl > 0.0 {
                    state.mu = mu_compl.max(options.mu_min).min(mu_ceiling);
                }

                if move_dist < 1e-14 && primal_inf < options.tol {
                    // Restoration returned same point AND constraints already satisfied.
                    consecutive_zero_alpha += 1;
                    if consecutive_zero_alpha >= stall_threshold && total_stall_recovery_count < 6 {
                        stall_recovery_count += 1;
                        total_stall_recovery_count += 1;
                        filter.reset();
                        inertia_params.delta_w_last = 0.0;
                        consecutive_zero_alpha = 0;

                        if stall_recovery_count >= 2 {
                            // Repeated stalls: use varied mu perturbation
                            let factor = mu_perturbation_factors[total_stall_recovery_count % mu_perturbation_factors.len()];
                            // Don't apply mu_ceiling during stall recovery
                            let mu_new = (state.mu * factor).max(options.mu_min).min(1e4);
                            log::debug!(
                                "Repeated stall (recovery #{}/{}), mu: {:.2e} -> {:.2e} (x{})",
                                stall_recovery_count, total_stall_recovery_count, state.mu, mu_new, factor
                            );
                            state.mu = mu_new;
                            stall_recovery_count = 0;
                            iters_since_stall_recovery = 0;

                            // After 3rd total recovery, also perturb x
                            if total_stall_recovery_count >= 3 {
                                for i in 0..n {
                                    let range = if state.x_l[i].is_finite() && state.x_u[i].is_finite() {
                                        state.x_u[i] - state.x_l[i]
                                    } else {
                                        state.x[i].abs().max(1.0)
                                    };
                                    let sign = if (i * 7 + total_stall_recovery_count * 13) % 3 == 0 { -1.0 } else { 1.0 };
                                    state.x[i] += sign * 1e-4 * range;
                                    if state.x_l[i].is_finite() {
                                        state.x[i] = state.x[i].max(state.x_l[i] + 1e-14);
                                    }
                                    if state.x_u[i].is_finite() {
                                        state.x[i] = state.x[i].min(state.x_u[i] - 1e-14);
                                    }
                                }
                                state.evaluate(problem, 1.0);
                            }
                        } else {
                            log::debug!("Stall detected, resetting filter (recovery #{})", stall_recovery_count);
                        }
                    }
                } else {
                    consecutive_zero_alpha = 0;
                }

                continue;
            } else {
                // Try recovery as last resort before giving up.
                // Vary mu perturbation direction and perturb x on later attempts.
                // Only try if we haven't done a recovery recently (cooldown period).
                let mu_factors = [10.0, 0.1, 100.0, 0.01, 1000.0, 0.001];
                if stall_recovery_count < 6 {
                    if iters_since_stall_recovery < 5 {
                        // Cooldown: skip recovery but continue iterating.
                        // The current iterate may make progress with the last mu perturbation.
                        continue;
                    }
                    let factor = mu_factors[stall_recovery_count % mu_factors.len()];
                    // Don't apply mu_ceiling during emergency recovery —
                    // we need to explore widely to escape the current basin.
                    state.mu = (state.mu * factor).max(options.mu_min).min(1e4);
                    filter.reset();
                    let theta_now = state.constraint_violation();
                    filter.set_theta_min_from_initial(theta_now);

                    // On attempts 3+: also perturb x to escape current basin
                    if stall_recovery_count >= 3 {
                        for i in 0..n {
                            let range = if state.x_l[i].is_finite() && state.x_u[i].is_finite() {
                                state.x_u[i] - state.x_l[i]
                            } else {
                                state.x[i].abs().max(1.0)
                            };
                            // Quasi-random sign based on index
                            let sign = if (i * 7 + stall_recovery_count * 13) % 3 == 0 { -1.0 } else { 1.0 };
                            state.x[i] += sign * 1e-4 * range;
                            if state.x_l[i].is_finite() {
                                state.x[i] = state.x[i].max(state.x_l[i] + 1e-14);
                            }
                            if state.x_u[i].is_finite() {
                                state.x[i] = state.x[i].min(state.x_u[i] - 1e-14);
                            }
                        }
                        state.evaluate(problem, 1.0);
                    }

                    stall_recovery_count += 1;
                    iters_since_stall_recovery = 0;
                    log::debug!(
                        "Restoration failed, mu perturbation x{} to {:.2e} (attempt {})",
                        factor, state.mu, stall_recovery_count
                    );
                    continue;
                }
                log::warn!("Restoration failed at iteration {}", iteration);
                let current_theta = state.constraint_violation();

                // Check if we're at a stationary point of the violation function.
                // If ||∇θ|| ≈ 0 but θ > tol, declare LocalInfeasibility.
                // Compute ∇(0.5*||violation||^2) = J^T * violation
                if current_theta > options.constr_viol_tol && !ever_feasible {
                    let mut violation = vec![0.0; m];
                    for i in 0..m {
                        let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
                        if is_eq {
                            violation[i] = state.g[i] - state.g_l[i];
                        } else if state.g_l[i].is_finite() && state.g[i] < state.g_l[i] {
                            violation[i] = state.g[i] - state.g_l[i];
                        } else if state.g_u[i].is_finite() && state.g[i] > state.g_u[i] {
                            violation[i] = state.g[i] - state.g_u[i];
                        }
                    }
                    // grad_theta = J^T * violation
                    let mut grad_theta = vec![0.0; n];
                    for (idx, (&row, &col)) in
                        state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
                    {
                        grad_theta[col] += state.jac_vals[idx] * violation[row];
                    }
                    let grad_theta_norm = grad_theta.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
                    // Scale threshold by violation magnitude for robustness
                    let stationarity_tol = 1e-4 * current_theta.max(1.0);
                    if grad_theta_norm < stationarity_tol {
                        log::info!(
                            "Local infeasibility detected: theta={:.2e}, ||∇theta||={:.2e}",
                            current_theta, grad_theta_norm
                        );
                        return make_result(&state, SolveStatus::LocalInfeasibility);
                    }
                }

                // Check if constraint violation hasn't decreased meaningfully
                // over recent iterations. Never declare infeasible if we were ever feasible.
                if !ever_feasible && current_theta > 1e4 && iteration > 500 && theta_history.len() >= theta_history_len {
                    let min_theta = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
                    if current_theta > 0.01 * min_theta {
                        // No 100x improvement over best in last 100 iterations → infeasible
                        return make_result(&state, SolveStatus::Infeasible);
                    }
                }
                return make_result(&state, SolveStatus::RestorationFailed);
            }
        }

        // Step was accepted — reset stall counter (any accepted step is progress)
        consecutive_zero_alpha = 0;

        // Watchdog: track consecutive shortened steps
        if state.alpha_primal < alpha_primal_max * 0.99 {
            consecutive_shortened += 1;
        } else {
            consecutive_shortened = 0;
        }

        // Watchdog activation: save state when too many shortened steps
        if !watchdog_active
            && consecutive_shortened >= options.watchdog_shortened_iter_trigger
        {
            watchdog_active = true;
            watchdog_trial_count = 0;
            let wd_theta = state.constraint_violation();
            let wd_phi = state.barrier_objective(options);
            watchdog_saved = Some(WatchdogSavedState {
                x: state.x.clone(),
                y: state.y.clone(),
                z_l: state.z_l.clone(),
                z_u: state.z_u.clone(),
                mu: state.mu,
                obj: state.obj,
                g: state.g.clone(),
                grad_f: state.grad_f.clone(),
                filter_entries: filter.save_entries(),
                theta: wd_theta,
                phi: wd_phi,
            });
            consecutive_shortened = 0;
            log::debug!(
                "Watchdog activated at iteration {} (theta={:.2e}, phi={:.2e})",
                iteration, wd_theta, wd_phi
            );
        }

        // Watchdog progress check
        if watchdog_active {
            watchdog_trial_count += 1;
            if let Some(ref saved) = watchdog_saved {
                let theta_now = state.constraint_violation();
                let phi_now = state.barrier_objective(options);
                // Check if current point is filter-acceptable from saved state
                let made_progress = filter.is_acceptable(theta_now, phi_now)
                    && (theta_now < (1.0 - 1e-5) * saved.theta
                        || phi_now < saved.phi - 1e-5 * saved.theta);

                if made_progress {
                    log::debug!(
                        "Watchdog succeeded at trial {} (theta: {:.2e} -> {:.2e})",
                        watchdog_trial_count, saved.theta, theta_now
                    );
                    watchdog_active = false;
                    watchdog_trial_count = 0;
                    watchdog_saved = None;
                } else if watchdog_trial_count >= options.watchdog_trial_iter_max {
                    // Revert to saved state
                    log::debug!(
                        "Watchdog reverting after {} trials",
                        watchdog_trial_count
                    );
                    // Add explored region to filter before reverting
                    let theta_now = state.constraint_violation();
                    let phi_now = state.barrier_objective(options);
                    filter.restore_entries(saved.filter_entries.clone());
                    filter.add(theta_now, phi_now);

                    state.x = saved.x.clone();
                    state.y = saved.y.clone();
                    state.z_l = saved.z_l.clone();
                    state.z_u = saved.z_u.clone();
                    state.mu = saved.mu;
                    state.obj = saved.obj;
                    state.g = saved.g.clone();
                    state.grad_f = saved.grad_f.clone();
                    state.evaluate(problem, 1.0);

                    watchdog_active = false;
                    watchdog_trial_count = 0;
                    watchdog_saved = None;
                    continue;
                }
            }
        }

        // Alpha history stall detection (Phase 3)
        // Only trigger when NOT near convergence — bumping mu near the optimum is counterproductive.
        if alpha_history.len() >= alpha_history_len {
            alpha_history.remove(0);
        }
        alpha_history.push(state.alpha_primal);
        if alpha_history.len() >= alpha_history_len {
            let avg_alpha: f64 =
                alpha_history.iter().sum::<f64>() / alpha_history.len() as f64;
            let near_converged = primal_inf < options.acceptable_tol
                && dual_inf < options.acceptable_tol * 100.0;
            if avg_alpha < 1e-6 && !near_converged && total_stall_recovery_count < 6 {
                // Use varied mu perturbation from unified counter
                let factor = mu_perturbation_factors[total_stall_recovery_count % mu_perturbation_factors.len()];
                log::debug!(
                    "Alpha stall detected (avg={:.2e}), resetting filter and mu x{} (total recovery #{})",
                    avg_alpha, factor, total_stall_recovery_count
                );
                filter.reset();
                // Don't apply mu_ceiling during stall recovery
                state.mu = (state.mu * factor).max(options.mu_min).min(1e4);
                iters_since_stall_recovery = 0;
                total_stall_recovery_count += 1;
                alpha_history.clear();
            }
        }

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

        // NaN/Inf guard on evaluation
        if state.obj.is_nan() || state.obj.is_infinite() {
            // Try restoration from current point
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
                if !state.obj.is_nan() && !state.obj.is_infinite() {
                    continue;
                }
            }
            return make_result(&state, SolveStatus::NumericalError);
        }

        // Track best feasible point for max_iter exit
        {
            let theta_now = state.constraint_violation();
            if theta_now < options.constr_viol_tol && state.obj < best_obj {
                best_obj = state.obj;
                best_x = Some(state.x.clone());
                best_y = Some(state.y.clone());
                best_z_l = Some(state.z_l.clone());
                best_z_u = Some(state.z_u.clone());
            }
        }

        // Fix 3B: During cooldown after stall recovery, allow mu to decrease
        // even when alpha is tiny (bypass the alpha < 1e-4 freeze).
        let in_cooldown = iters_since_stall_recovery < 20;
        iters_since_stall_recovery = iters_since_stall_recovery.saturating_add(1);

        // Update barrier parameter, capped by mu_ceiling to prevent stall recovery
        // from undoing forced decreases.
        let mu_old = state.mu;
        state.mu = update_barrier_parameter(&state, options, in_cooldown).min(mu_ceiling);

        // Ipopt-like monotone mu decrease fallback: after 20 stuck iterations,
        // use min(linear_factor * mu, mu^superlinear_power) to force progress.
        // This matches Ipopt's fixed-mode decrease formula.
        if state.mu >= 0.99 * mu_old {
            mu_stuck_count += 1;
        } else {
            mu_stuck_count = 0;
        }
        if mu_stuck_count >= 20 {
            let new_mu = (options.mu_linear_decrease_factor * state.mu)
                .min(state.mu.powf(options.mu_superlinear_decrease_power))
                .max(options.mu_min);
            state.mu = new_mu;
            if ever_feasible {
                mu_ceiling = state.mu; // Only lock ceiling when known feasible
            }
            mu_stuck_count = 0;
            filter.reset();
            let theta_new = state.constraint_violation();
            filter.set_theta_min_from_initial(theta_new);
            log::debug!("Forced mu decrease to {:.2e} (Ipopt monotone, stuck for 20 iters)", state.mu);
        }

        // Filter reset on large mu decrease (>5x) — stale entries from previous
        // barrier problem can block progress
        if state.mu < 0.2 * mu_old {
            filter.reset();
            let theta_new = state.constraint_violation();
            filter.set_theta_min_from_initial(theta_new);
        }

        // Post-step acceptable convergence tracking.
        // This catches cases where the step just taken pushes the state into the
        // acceptable region but the pre-step check at the top of the loop missed it.
        {
            let post_primal = state.constraint_violation();
            let (post_zl_opt, post_zu_opt) = {
                let mut gj = state.grad_f.clone();
                for (idx, (&row, &col)) in
                    state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
                {
                    gj[col] += state.jac_vals[idx] * state.y[row];
                }
                let mut zl = vec![0.0; n];
                let mut zu = vec![0.0; n];
                let kc = 1e10;
                for i in 0..n {
                    if gj[i] > 0.0 && state.x_l[i].is_finite() {
                        let sl = (state.x[i] - state.x_l[i]).max(1e-20);
                        if gj[i] * sl <= kc * state.mu.max(1e-20) {
                            zl[i] = gj[i];
                        }
                    } else if gj[i] < 0.0 && state.x_u[i].is_finite() {
                        let su = (state.x_u[i] - state.x[i]).max(1e-20);
                        if (-gj[i]) * su <= kc * state.mu.max(1e-20) {
                            zu[i] = -gj[i];
                        }
                    }
                }
                (zl, zu)
            };
            let post_du = convergence::dual_infeasibility(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &state.y, &post_zl_opt, &post_zu_opt, n,
            );
            let post_du_unsc = convergence::dual_infeasibility(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &state.y, &state.z_l, &state.z_u, n,
            );
            let post_compl = convergence::complementarity_error(
                &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
            );
            let post_compl_opt = convergence::complementarity_error(
                &state.x, &state.x_l, &state.x_u, &post_zl_opt, &post_zu_opt, 0.0,
            );
            let post_compl_best = post_compl.min(post_compl_opt);
            let post_mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
                + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
                + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
            let post_sd = if (m + 2 * n) > 0 {
                ((100.0f64.max(post_mult_sum / (m + 2 * n) as f64)) / 100.0).min(1e4)
            } else {
                1.0
            };
            let post_acc_scaled = post_primal <= options.acceptable_tol
                && post_du <= options.acceptable_tol * post_sd
                && post_compl_best <= options.acceptable_tol * post_sd;
            let post_acc_unscaled = post_primal <= options.acceptable_constr_viol_tol
                && post_du_unsc <= options.acceptable_dual_inf_tol
                && post_compl_best <= options.acceptable_compl_inf_tol;
            if post_acc_scaled && post_acc_unscaled {
                state.consecutive_acceptable += 1;
            }
            // Don't reset here — the pre-step check handles resets
        }
    }

    // At max_iter: log convergence diagnostics using same z_opt as convergence check (with gate)
    {
        let final_primal_inf = state.constraint_violation();
        let (z_l_opt_final, z_u_opt_final) = {
            let mut grad_jty = state.grad_f.clone();
            for (idx, (&row, &col)) in
                state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
            {
                grad_jty[col] += state.jac_vals[idx] * state.y[row];
            }
            let mut zl = vec![0.0; n];
            let mut zu = vec![0.0; n];
            let kc = 1e10;
            for i in 0..n {
                if grad_jty[i] > 0.0 && state.x_l[i].is_finite() {
                    let sl = (state.x[i] - state.x_l[i]).max(1e-20);
                    if grad_jty[i] * sl <= kc * state.mu.max(1e-20) {
                        zl[i] = grad_jty[i];
                    }
                } else if grad_jty[i] < 0.0 && state.x_u[i].is_finite() {
                    let su = (state.x_u[i] - state.x[i]).max(1e-20);
                    if (-grad_jty[i]) * su <= kc * state.mu.max(1e-20) {
                        zu[i] = -grad_jty[i];
                    }
                }
            }
            (zl, zu)
        };
        let final_dual_inf = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &z_l_opt_final, &z_u_opt_final, n,
        );
        let final_dual_inf_unscaled = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        let final_compl = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
        );
        let final_compl_opt = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &z_l_opt_final, &z_u_opt_final, 0.0,
        );
        let final_compl_best = final_compl.min(final_compl_opt);
        let mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let s_max: f64 = 100.0;
        let s_d = if (m + 2 * n) > 0 {
            ((s_max.max(mult_sum / (m + 2 * n) as f64)) / s_max).min(1e4)
        } else {
            1.0
        };
        eprintln!(
            "ripopt: MaxIter diag: pr={:.2e} du={:.2e}(t={:.2e}) du_u={:.2e}(t={:.0e}) co={:.2e} co_opt={:.2e} co_best={:.2e}(t={:.2e}/{:.2e}) mu={:.2e} sd={:.1} ac={}",
            final_primal_inf,
            final_dual_inf, options.tol * s_d,
            final_dual_inf_unscaled, options.dual_inf_tol,
            final_compl, final_compl_opt, final_compl_best,
            options.acceptable_tol * s_d, options.acceptable_compl_inf_tol,
            state.mu, s_d, state.consecutive_acceptable,
        );
    }

    // At max_iter: final acceptable convergence check.
    // The counter may have oscillated during the solve, but if the final state
    // meets all acceptable criteria, the NLP solution is valid.
    {
        let final_primal = state.constraint_violation();
        let (fzl, fzu) = {
            let mut gj = state.grad_f.clone();
            for (idx, (&row, &col)) in
                state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
            {
                gj[col] += state.jac_vals[idx] * state.y[row];
            }
            let mut zl = vec![0.0; n];
            let mut zu = vec![0.0; n];
            let kc = 1e10;
            for i in 0..n {
                if gj[i] > 0.0 && state.x_l[i].is_finite() {
                    let sl = (state.x[i] - state.x_l[i]).max(1e-20);
                    if gj[i] * sl <= kc * state.mu.max(1e-20) {
                        zl[i] = gj[i];
                    }
                } else if gj[i] < 0.0 && state.x_u[i].is_finite() {
                    let su = (state.x_u[i] - state.x[i]).max(1e-20);
                    if (-gj[i]) * su <= kc * state.mu.max(1e-20) {
                        zu[i] = -gj[i];
                    }
                }
            }
            (zl, zu)
        };
        let fdu = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &fzl, &fzu, n,
        );
        let fdu_u = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        let fco = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
        );
        let fco_opt = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &fzl, &fzu, 0.0,
        );
        let fco_best = fco.min(fco_opt);
        let fmult: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let fsd = if (m + 2 * n) > 0 {
            ((100.0f64.max(fmult / (m + 2 * n) as f64)) / 100.0).min(1e4)
        } else {
            1.0
        };
        // At max_iter exit, use a floor of 1e-2 on scaled tolerances to catch
        // solutions where the gradient is within engineering precision but the
        // acceptable_tol * s_d threshold is too tight (e.g., unconstrained with s_d=1).
        let fdu_tol = (options.acceptable_tol * fsd).max(1e-2);
        let fco_tol = (options.acceptable_tol * fsd).max(1e-2);
        let sc = final_primal <= options.acceptable_tol
            && fdu <= fdu_tol
            && fco_best <= fco_tol;
        let usc = final_primal <= options.acceptable_constr_viol_tol
            && fdu_u <= options.acceptable_dual_inf_tol
            && fco_best <= options.acceptable_compl_inf_tol;
        if sc && usc {
            return make_result(&state, SolveStatus::Acceptable);
        }
    }

    // At max_iter: try restoring the best feasible point we saw during the solve.
    // The current point may have degraded, but an earlier point could pass acceptable.
    if let Some(ref bx) = best_x {
        state.x = bx.clone();
        if let Some(ref by) = best_y { state.y = by.clone(); }
        if let Some(ref bzl) = best_z_l { state.z_l = bzl.clone(); }
        if let Some(ref bzu) = best_z_u { state.z_u = bzu.clone(); }
        state.evaluate(problem, 1.0);

        // Re-check acceptable convergence at the best point with relaxed tolerance floor
        let bp_primal = state.constraint_violation();
        let (bp_zl, bp_zu) = {
            let mut gj = state.grad_f.clone();
            for (idx, (&row, &col)) in
                state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
            {
                gj[col] += state.jac_vals[idx] * state.y[row];
            }
            let mut zl = vec![0.0; n];
            let mut zu = vec![0.0; n];
            let kc = 1e10;
            for i in 0..n {
                if gj[i] > 0.0 && state.x_l[i].is_finite() {
                    let sl = (state.x[i] - state.x_l[i]).max(1e-20);
                    if gj[i] * sl <= kc * state.mu.max(1e-20) {
                        zl[i] = gj[i];
                    }
                } else if gj[i] < 0.0 && state.x_u[i].is_finite() {
                    let su = (state.x_u[i] - state.x[i]).max(1e-20);
                    if (-gj[i]) * su <= kc * state.mu.max(1e-20) {
                        zu[i] = -gj[i];
                    }
                }
            }
            (zl, zu)
        };
        let bp_du = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &bp_zl, &bp_zu, n,
        );
        let bp_du_u = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        let bp_co = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
        );
        let bp_co_opt = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &bp_zl, &bp_zu, 0.0,
        );
        let bp_co_best = bp_co.min(bp_co_opt);
        let bp_mult: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let bp_sd = if (m + 2 * n) > 0 {
            ((100.0f64.max(bp_mult / (m + 2 * n) as f64)) / 100.0).min(1e4)
        } else {
            1.0
        };
        // Relaxed tolerances for best-point check at max_iter.
        // The solver exhausted its iterations, so we use a generous floor:
        // - du: floor of 1.0 (was 1e-1) to catch near-converged problems
        // - co: floor of 1.0 (barriers haven't fully resolved but point is good)
        let bp_du_tol = (options.acceptable_tol * bp_sd).max(1.0);
        let bp_co_tol = (options.acceptable_tol * bp_sd).max(1.0);
        let bp_sc = bp_primal <= options.acceptable_tol
            && bp_du <= bp_du_tol
            && bp_co_best <= bp_co_tol;
        // Unscaled check with 10x relaxed complementarity for best-point exit
        let bp_usc = bp_primal <= options.acceptable_constr_viol_tol
            && bp_du_u <= options.acceptable_dual_inf_tol
            && bp_co_best <= options.acceptable_compl_inf_tol * 100.0;
        if bp_sc && bp_usc {
            return make_result(&state, SolveStatus::Acceptable);
        }
    }

    // At max_iter: check if the problem is actually infeasible.
    // Never declare infeasible if we were ever feasible.
    let final_theta = state.constraint_violation();
    if !ever_feasible && final_theta > options.constr_viol_tol {
        // Check stationarity of violation: if ||∇θ|| ≈ 0, declare local infeasibility
        let mut violation = vec![0.0; m];
        for i in 0..m {
            let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_eq {
                violation[i] = state.g[i] - state.g_l[i];
            } else if state.g_l[i].is_finite() && state.g[i] < state.g_l[i] {
                violation[i] = state.g[i] - state.g_l[i];
            } else if state.g_u[i].is_finite() && state.g[i] > state.g_u[i] {
                violation[i] = state.g[i] - state.g_u[i];
            }
        }
        let mut grad_theta = vec![0.0; n];
        for (idx, (&row, &col)) in
            state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
        {
            grad_theta[col] += state.jac_vals[idx] * violation[row];
        }
        let grad_theta_norm = grad_theta.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        let stationarity_tol = 1e-4 * final_theta.max(1.0);
        if grad_theta_norm < stationarity_tol {
            return make_result(&state, SolveStatus::LocalInfeasibility);
        }

        // Fallback: check if theta hasn't improved over recent history
        if final_theta > 1e4 && theta_history.len() >= theta_history_len {
            let min_theta = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
            if final_theta > 0.01 * min_theta {
                return make_result(&state, SolveStatus::Infeasible);
            }
        }
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

/// Update the barrier parameter mu using Mehrotra centering.
fn update_barrier_parameter(state: &SolverState, options: &SolverOptions, force_decrease: bool) -> f64 {
    let mu = state.mu;

    if options.mu_strategy_adaptive {
        // Compute average complementarity
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

        // Adaptive (Loqo) rule: mu = avg_compl / kappa
        let mu_new = (avg_compl / options.kappa).max(options.mu_min);

        // Ipopt-like mu update: when alpha is tiny, keep mu unchanged.
        // The old logic (mu_new = mu.max(mu_new)) prevented decrease but allowed
        // increase, creating a runaway loop where mu climbs to 1e4 and stays there.
        // Ipopt switches to "fixed mode" with monotone decrease; we approximate
        // by keeping mu constant and relying on the forced decrease mechanism.
        if state.alpha_primal < 1e-4 && !force_decrease {
            return mu;
        }

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

/// Compute average complementarity for recomputing mu after restoration.
fn compute_avg_complementarity(state: &SolverState) -> f64 {
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
    if count > 0 {
        sum_compl / count as f64
    } else {
        0.0
    }
}

/// Compute a steepest-descent fallback direction when KKT solve fails.
///
/// Returns (dx, dy) where dx = -alpha * grad_f (scaled gradient step)
/// and dy = 0. This is crude but prevents immediate failure.
fn gradient_descent_fallback(state: &SolverState) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = state.n;
    let m = state.m;
    let grad_norm: f64 = state.grad_f.iter().map(|g| g * g).sum::<f64>().sqrt();
    if grad_norm < 1e-20 {
        return None;
    }
    let alpha = 1e-4 / grad_norm; // small step
    let mut dx = vec![0.0; n];
    for i in 0..n {
        dx[i] = -alpha * state.grad_f[i];
    }
    let dy = vec![0.0; m];
    Some((dx, dy))
}

/// Build the final solve result.
/// Computes z from stationarity for more accurate output multipliers.
/// Unscales all values from the internal scaled space to the original NLP space.
fn make_result(state: &SolverState, status: SolveStatus) -> SolveResult {
    // Compute optimal z from stationarity in scaled space: ∇f_s + J_s^T y_s - z_l_s + z_u_s = 0
    let n = state.n;
    let m = state.m;
    let mut grad_jty = state.grad_f.clone();
    for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        grad_jty[col] += state.jac_vals[idx] * state.y[row];
    }

    // z from stationarity (scaled), then unscale: z_unscaled = z_scaled / obj_scaling
    let mut z_l_out = vec![0.0; n];
    let mut z_u_out = vec![0.0; n];
    for i in 0..n {
        if grad_jty[i] > 0.0 && state.x_l[i].is_finite() {
            z_l_out[i] = grad_jty[i] / state.obj_scaling;
        } else if grad_jty[i] < 0.0 && state.x_u[i].is_finite() {
            z_u_out[i] = -grad_jty[i] / state.obj_scaling;
        }
    }

    // Unscale constraint multipliers: y_unscaled[i] = y_scaled[i] * g_scaling[i] / obj_scaling
    let mut y_out = state.y.clone();
    for i in 0..m {
        y_out[i] = state.y[i] * state.g_scaling[i] / state.obj_scaling;
    }

    // Unscale constraint values: g_unscaled[i] = g_scaled[i] / g_scaling[i]
    let mut g_out = state.g.clone();
    for i in 0..m {
        g_out[i] /= state.g_scaling[i];
    }

    SolveResult {
        x: state.x.clone(),
        objective: state.obj / state.obj_scaling,
        constraint_multipliers: y_out,
        bound_multipliers_lower: z_l_out,
        bound_multipliers_upper: z_u_out,
        constraint_values: g_out,
        status,
        iterations: state.iter,
    }
}
