use std::time::{Duration, Instant};

use crate::convergence::{self, check_convergence, ConvergenceInfo, ConvergenceStatus};
use crate::filter::{self, Filter, FilterEntry};
use crate::kkt::{self, InertiaCorrectionParams};
use crate::linear_solver::banded::BandedLdl;
use crate::linear_solver::dense::DenseLdl;
#[cfg(feature = "faer")]
use crate::linear_solver::sparse::SparseLdl;
#[cfg(feature = "rmumps")]
use crate::linear_solver::multifrontal::MultifrontalLdl;
use crate::linear_solver::{KktMatrix, LinearSolver, SymmetricMatrix};

/// Create a new sparse linear solver using the best available backend.
/// Prefers rmumps (multifrontal) when available, falls back to faer (SparseLdl).
fn new_sparse_solver() -> Box<dyn LinearSolver> {
    #[cfg(feature = "rmumps")]
    { return Box::new(MultifrontalLdl::new()); }
    #[cfg(all(not(feature = "rmumps"), feature = "faer"))]
    { return Box::new(SparseLdl::new()); }
    #[cfg(not(any(feature = "rmumps", feature = "faer")))]
    { return Box::new(DenseLdl::new()); }
}

/// Create the appropriate linear solver for a fallback KKT system.
/// Uses sparse solver when `use_sparse` is true, dense otherwise.
fn new_fallback_solver(use_sparse: bool) -> Box<dyn LinearSolver> {
    if use_sparse {
        new_sparse_solver()
    } else {
        Box::new(DenseLdl::new())
    }
}
use crate::options::SolverOptions;
use crate::problem::NlpProblem;
use crate::restoration::RestorationPhase;
use crate::restoration_nlp::RestorationNlp;
use crate::result::{SolveResult, SolverDiagnostics, SolveStatus};
use crate::slack_formulation::SlackFormulation;
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
    v_l: Vec<f64>,
    v_u: Vec<f64>,
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
    /// Constraint slack lower-bound multipliers (Ipopt's v_L).
    /// v_l[i] > 0 for inequality constraints with finite g_l[i], 0 otherwise.
    pub v_l: Vec<f64>,
    /// Constraint slack upper-bound multipliers (Ipopt's v_U).
    /// v_u[i] > 0 for inequality constraints with finite g_u[i], 0 otherwise.
    pub v_u: Vec<f64>,
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
    /// Accumulated solver diagnostics.
    pub diagnostics: SolverDiagnostics,
}

/// Barrier parameter mode (Ipopt's adaptive mu strategy).
#[derive(Debug, Clone, Copy, PartialEq)]
enum MuMode {
    /// Free mode: mu chosen by oracle each iteration. Can increase or decrease.
    Free,
    /// Fixed mode: monotone mu decrease. Subproblem solved to barrier_tol_factor * mu.
    Fixed,
}

/// State for the free/fixed mu update strategy.
struct MuState {
    mode: MuMode,
    /// Sliding window of KKT error values for progress tracking.
    ref_vals: Vec<f64>,
    /// Maximum reference values to keep.
    num_refs_max: usize,
    /// Required reduction factor (sufficient progress if error < refs_red_fact * any ref).
    refs_red_fact: f64,
    /// Flag for tiny step detection.
    tiny_step: bool,
    /// Flag to indicate first iteration after mode switch.
    first_iter_in_mode: bool,
    /// Count of consecutive restoration failures for giving up.
    consecutive_restoration_failures: usize,
}

impl MuState {
    fn new() -> Self {
        Self {
            mode: MuMode::Free,
            ref_vals: Vec::with_capacity(8),
            num_refs_max: 4,
            refs_red_fact: 0.9999,
            tiny_step: false,
            first_iter_in_mode: true,
            consecutive_restoration_failures: 0,
        }
    }

    /// Check if sufficient progress is being made (KKT error reference check).
    fn check_sufficient_progress(&self, kkt_error: f64) -> bool {
        if self.ref_vals.len() < self.num_refs_max {
            return true; // Not enough history yet
        }
        // Sufficient if current error < refs_red_fact * any reference
        self.ref_vals.iter().any(|&r| kkt_error <= self.refs_red_fact * r)
    }

    /// Remember an accepted KKT error value.
    fn remember_accepted(&mut self, kkt_error: f64) {
        if self.ref_vals.len() >= self.num_refs_max {
            self.ref_vals.remove(0);
        }
        self.ref_vals.push(kkt_error);
    }
}

/// L-BFGS Hessian approximation state for use inside the IPM loop.
///
/// Maintains curvature pairs (s_k, y_k) from Lagrangian gradient differences
/// and forms an explicit dense B_k matrix for the KKT system.
pub struct LbfgsIpmState {
    /// Number of primal variables.
    n: usize,
    /// Maximum number of stored pairs.
    m_max: usize,
    /// Stored s_k vectors (x_{k+1} - x_k).
    s_store: Vec<Vec<f64>>,
    /// Stored y_k vectors (∇L_{k+1} - ∇L_k, evaluated at new multipliers).
    y_store: Vec<Vec<f64>>,
    /// Previous iterate x_k.
    prev_x: Vec<f64>,
    /// Previous Lagrangian gradient ∇_x L(x_k, λ_{k+1}).
    prev_lag_grad: Vec<f64>,
    /// Whether we have a previous iterate (skip update on first call).
    has_prev: bool,
    /// Initial Hessian scaling factor gamma (H0 = gamma * I).
    gamma: f64,
}

impl LbfgsIpmState {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            m_max: 10,
            s_store: Vec::with_capacity(10),
            y_store: Vec::with_capacity(10),
            prev_x: vec![0.0; n],
            prev_lag_grad: vec![0.0; n],
            has_prev: false,
            gamma: 1.0,
        }
    }

    /// Compute ∇_x L = ∇f + J^T λ
    fn compute_lagrangian_gradient(
        grad_f: &[f64],
        jac_rows: &[usize],
        jac_cols: &[usize],
        jac_vals: &[f64],
        lambda: &[f64],
        n: usize,
    ) -> Vec<f64> {
        let mut lag_grad = grad_f.to_vec();
        for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
            // J^T λ: column `col` of J^T gets contribution from row `row`
            if col < n {
                lag_grad[col] += jac_vals[idx] * lambda[row];
            }
        }
        lag_grad
    }

    /// Update L-BFGS pairs after a step has been accepted.
    /// Uses Powell damping to ensure positive curvature (s^T y > 0).
    pub fn update(
        &mut self,
        new_x: &[f64],
        new_lag_grad: &[f64],
    ) {
        if !self.has_prev {
            self.prev_x.copy_from_slice(new_x);
            self.prev_lag_grad.copy_from_slice(new_lag_grad);
            self.has_prev = true;
            return;
        }

        let n = self.n;
        let mut s_k = vec![0.0; n];
        let mut y_k = vec![0.0; n];
        for i in 0..n {
            s_k[i] = new_x[i] - self.prev_x[i];
            y_k[i] = new_lag_grad[i] - self.prev_lag_grad[i];
        }

        let ss: f64 = s_k.iter().map(|v| v * v).sum();
        if ss < 1e-30 {
            // Step too small, skip update
            self.prev_x.copy_from_slice(new_x);
            self.prev_lag_grad.copy_from_slice(new_lag_grad);
            return;
        }

        let sy: f64 = s_k.iter().zip(y_k.iter()).map(|(s, y)| s * y).sum();

        // Compute B_k * s_k for Powell damping
        let bs = self.multiply_bk(&s_k);
        let sbs: f64 = s_k.iter().zip(bs.iter()).map(|(s, b)| s * b).sum();

        // Powell damping: ensure s^T y >= 0.2 * s^T B s
        if sy >= 0.2 * sbs {
            // Use y_k as-is
        } else {
            let theta = if (sbs - sy).abs() < 1e-30 {
                1.0
            } else {
                0.8 * sbs / (sbs - sy)
            };
            for i in 0..n {
                y_k[i] = theta * y_k[i] + (1.0 - theta) * bs[i];
            }
        }

        // Verify positive curvature after damping
        let sy_damped: f64 = s_k.iter().zip(y_k.iter()).map(|(s, y)| s * y).sum();
        if sy_damped <= 1e-20 {
            self.prev_x.copy_from_slice(new_x);
            self.prev_lag_grad.copy_from_slice(new_lag_grad);
            return;
        }

        // Update gamma = s^T y / y^T y
        let yy: f64 = y_k.iter().map(|v| v * v).sum();
        if yy > 1e-30 {
            self.gamma = sy_damped / yy;
        }

        // Store pair
        if self.s_store.len() == self.m_max {
            self.s_store.remove(0);
            self.y_store.remove(0);
        }
        self.s_store.push(s_k);
        self.y_store.push(y_k);

        self.prev_x.copy_from_slice(new_x);
        self.prev_lag_grad.copy_from_slice(new_lag_grad);
    }

    /// Compute B_k * v using the L-BFGS compact representation.
    /// B_k = gamma^{-1} I - ... (inverse Hessian formulation, then invert).
    /// Instead, we directly build B_k * v using the recursive formula.
    pub fn multiply_bk(&self, v: &[f64]) -> Vec<f64> {
        let n = self.n;
        let k = self.s_store.len();

        if k == 0 {
            // B_0 = (1/gamma) * I
            let scale = 1.0 / self.gamma.max(1e-12);
            return v.iter().map(|&vi| scale * vi).collect();
        }

        // Use the explicit B_k formation and multiply
        // B_k = (1/gamma) I + sum of rank-2 updates
        // It's easier to just form the full matrix and multiply for correctness
        let mut result = vec![0.0; n];
        let bk = self.form_dense_bk();
        for i in 0..n {
            for j in 0..n {
                let (r, c) = if i >= j { (i, j) } else { (j, i) };
                let idx = r * (r + 1) / 2 + c;
                result[i] += bk[idx] * v[j];
            }
        }
        result
    }

    /// Form explicit dense B_k matrix in lower-triangle format.
    /// Uses the L-BFGS compact representation:
    ///   B_k = B_0 - [B_0 S_k  Y_k] * M^{-1} * [S_k^T B_0; Y_k^T]
    /// where B_0 = (1/gamma) I.
    pub fn form_dense_bk(&self) -> Vec<f64> {
        let n = self.n;
        let k = self.s_store.len();
        let nnz = n * (n + 1) / 2;

        // Start with B_0 = (1/gamma) * I
        let b0_diag = 1.0 / self.gamma.max(1e-12);
        let mut bk = vec![0.0; nnz];
        for i in 0..n {
            let idx = i * (i + 1) / 2 + i;
            bk[idx] = b0_diag;
        }

        if k == 0 {
            return bk;
        }

        // Build B_k iteratively using rank-2 updates (BFGS formula):
        // B_{k+1} = B_k - (B_k s_k s_k^T B_k) / (s_k^T B_k s_k) + (y_k y_k^T) / (s_k^T y_k)
        for p in 0..k {
            let s = &self.s_store[p];
            let y = &self.y_store[p];

            // Compute B_k * s
            let mut bs = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    let (r, c) = if i >= j { (i, j) } else { (j, i) };
                    let idx = r * (r + 1) / 2 + c;
                    bs[i] += bk[idx] * s[j];
                }
            }

            let sbs: f64 = s.iter().zip(bs.iter()).map(|(si, bsi)| si * bsi).sum();
            let sy: f64 = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum();

            if sbs.abs() < 1e-30 || sy.abs() < 1e-30 {
                continue;
            }

            // Rank-2 update: B -= (bs bs^T) / sbs + (y y^T) / sy
            for i in 0..n {
                for j in 0..=i {
                    let idx = i * (i + 1) / 2 + j;
                    bk[idx] += -bs[i] * bs[j] / sbs + y[i] * y[j] / sy;
                }
            }
        }

        bk
    }

    /// Fill the hess_vals buffer with the dense B_k matrix.
    fn fill_hessian(&self, hess_vals: &mut [f64]) {
        let bk = self.form_dense_bk();
        hess_vals[..bk.len()].copy_from_slice(&bk);
    }
}

/// Generate dense lower-triangle sparsity pattern for n variables.
fn dense_lower_triangle_pattern(n: usize) -> (Vec<usize>, Vec<usize>) {
    let nnz = n * (n + 1) / 2;
    let mut rows = Vec::with_capacity(nnz);
    let mut cols = Vec::with_capacity(nnz);
    for i in 0..n {
        for j in 0..=i {
            rows.push(i);
            cols.push(j);
        }
    }
    (rows, cols)
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

        // Relax fixed variables: when x_l == x_u, the variable is fixed.
        // Interior-point methods require strictly interior starting points,
        // so we relax the bounds slightly (Ipopt's relax_bounds approach).
        for i in 0..n {
            if x_l[i].is_finite() && x_u[i].is_finite() && (x_u[i] - x_l[i]).abs() < 1e-10 {
                let center = (x_l[i] + x_u[i]) / 2.0;
                let relax = 1e-8 * center.abs().max(1.0);
                x_l[i] = center - relax;
                x_u[i] = center + relax;
            }
        }

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
        let (hess_rows, hess_cols) = if options.hessian_approximation_lbfgs {
            dense_lower_triangle_pattern(n)
        } else {
            problem.hessian_structure()
        };
        let hess_nnz = hess_rows.len();

        // Initialize constraint multipliers via least-squares estimate if enabled.
        // Solves min ||∇f + J^T y||^2  ⟹  (J J^T) y = -J ∇f
        // Gate LS mult init on problem size: dense J*J^T is O(m^2*n), too slow for large problems
        let ls_init_dim_limit = 500;
        let y = if options.least_squares_mult_init && m > 0 && (m + n) <= ls_init_dim_limit {
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
            let factored = ls_solver.bunch_kaufman_factor(&a_mat);
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
                        if convergence::is_equality_constraint(g_l[i], g_u[i]) {
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
            v_l: vec![0.0; m],
            v_u: vec![0.0; m],
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
            diagnostics: SolverDiagnostics::default(),
        }
    }

    /// Evaluate all functions, zeroing Hessian lambda for linear constraints.
    /// When `skip_hessian` is true (L-BFGS mode), the Hessian evaluation is skipped.
    fn evaluate_with_linear<P: NlpProblem>(
        &mut self,
        problem: &P,
        obj_factor: f64,
        linear_constraints: Option<&[bool]>,
        skip_hessian: bool,
    ) {
        self.obj = problem.objective(&self.x);
        problem.gradient(&self.x, &mut self.grad_f);
        if self.m > 0 {
            problem.constraints(&self.x, &mut self.g);
            problem.jacobian_values(&self.x, &mut self.jac_vals);
        }
        if skip_hessian {
            return;
        }
        if let Some(flags) = linear_constraints {
            let mut lambda_for_hess = self.y.clone();
            for (i, &is_lin) in flags.iter().enumerate() {
                if is_lin {
                    lambda_for_hess[i] = 0.0;
                }
            }
            problem.hessian_values(&self.x, obj_factor, &lambda_for_hess, &mut self.hess_vals);
        } else {
            problem.hessian_values(&self.x, obj_factor, &self.y, &mut self.hess_vals);
        }
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
    /// Hessian structure: lower triangle of J^T*J + ∑ r_i ∇²g_i.
    hess_rows: Vec<usize>,
    hess_cols: Vec<usize>,
    /// Mapping from inner hessian entries to our dense lower triangle index.
    /// inner_hess_map[k] = index into our vals[] for inner hessian entry k.
    inner_hess_map: Vec<usize>,
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

        // Build Hessian structure for J^T*J + ∑r_i∇²g_i (lower triangle, n x n dense).
        // Since J^T*J is generally dense, use full lower triangle.
        let mut hess_rows = Vec::with_capacity(n * (n + 1) / 2);
        let mut hess_cols = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in 0..=i {
                hess_rows.push(i);
                hess_cols.push(j);
            }
        }

        // Build mapping from inner hessian entries to our dense lower triangle.
        // Our layout: for (i,j) with i >= j, index = i*(i+1)/2 + j
        let (inner_hess_rows, inner_hess_cols) = inner.hessian_structure();
        let mut inner_hess_map = Vec::with_capacity(inner_hess_rows.len());
        for k in 0..inner_hess_rows.len() {
            let (r, c) = (inner_hess_rows[k], inner_hess_cols[k]);
            // Ensure lower triangle (r >= c)
            let (i, j) = if r >= c { (r, c) } else { (c, r) };
            let idx = i * (i + 1) / 2 + j;
            inner_hess_map.push(idx);
        }

        LeastSquaresProblem {
            inner,
            targets,
            jac_rows,
            jac_cols,
            hess_rows,
            hess_cols,
            inner_hess_map,
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
        let m = self.targets.len();

        // Compute residual r = g(x) - target
        let mut g = vec![0.0; m];
        self.inner.constraints(x, &mut g);
        let mut r = vec![0.0; m];
        for i in 0..m {
            r[i] = g[i] - self.targets[i];
        }

        // Part 1: J^T * J (Gauss-Newton term)
        let jac_nnz = self.jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];
        self.inner.jacobian_values(x, &mut jac_vals);

        let mut j_dense = vec![0.0; m * n];
        for (idx, (&row, &col)) in self.jac_rows.iter().zip(self.jac_cols.iter()).enumerate() {
            j_dense[row * n + col] += jac_vals[idx];
        }

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

        // Part 2: ∑ r_i * ∇²g_i (second-order correction)
        // inner.hessian_values(x, 0.0, r, hess) gives ∑ r_i * ∇²g_i
        let inner_hess_nnz = self.inner_hess_map.len();
        if inner_hess_nnz > 0 {
            let mut inner_hess_vals = vec![0.0; inner_hess_nnz];
            self.inner.hessian_values(x, 0.0, &r, &mut inner_hess_vals);
            for (k, &v) in inner_hess_vals.iter().enumerate() {
                vals[self.inner_hess_map[k]] += obj_factor * v;
            }
        }
    }
}

/// Detect if a problem is a nonlinear equation system (square or overdetermined).
///
/// Returns true if ALL of:
/// - f(x0) ≈ 0 (zero objective)
/// - ∇f(x0) ≈ 0 (zero gradient)
/// - All constraints are equalities (g_l[i] == g_u[i])
/// - m >= n (square or more constraints than variables)
fn detect_ne_problem<P: NlpProblem>(problem: &P) -> bool {
    let n = problem.num_variables();
    let m = problem.num_constraints();

    if m < n || m == 0 || n == 0 {
        return false;
    }

    // Square systems (m == n) are better solved by direct constrained IPM
    // (Newton on g(x)=0), which typically converges in a few iterations.
    // LS reformulation min 0.5*||g||^2 has a harder landscape for square systems.
    if m == n {
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
        if !convergence::is_equality_constraint(g_l[i], g_u[i]) {
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
    let solve_start = Instant::now();

    // --- Preprocessing: eliminate fixed variables and redundant constraints ---
    if options.enable_preprocessing {
        let prep = crate::preprocessing::PreprocessedProblem::new(problem as &dyn NlpProblem);
        if prep.did_reduce() {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: Preprocessing reduced problem: {} fixed vars, {} redundant constraints ({}x{} -> {}x{})",
                    prep.num_fixed(), prep.num_redundant(),
                    problem.num_variables(), problem.num_constraints(),
                    prep.num_variables(), prep.num_constraints(),
                );
            }
            let mut prep_opts = options.clone();
            prep_opts.enable_preprocessing = false; // prevent re-preprocessing
            let reduced_result = solve(&prep, &prep_opts);
            let result = prep.unmap_solution(&reduced_result);
            // If preprocessing made things worse, fall back to solving without it
            if matches!(result.status, SolveStatus::Optimal | SolveStatus::Acceptable) {
                return result;
            }
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: Preprocessed solve failed ({:?}), retrying without preprocessing",
                    result.status
                );
            }
            // Fall through to solve without preprocessing
        }
    }

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
        // For square systems (m == n), cap LS iterations — they often fail the LS
        // approach because there are no "extra" equations to drive residuals down.
        let mut ls_opts = options.clone();
        if m == n {
            ls_opts.max_iter = (options.max_iter / 10).min(100);
        }
        let ls_result = solve_ipm(&ls_problem, &ls_opts);

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

        // Fall back to constrained IPM when LS reports infeasibility:
        // - For square systems (m==n): LS may have found a non-root critical point
        //   of ||g||^2; constrained IPM can find the actual root.
        // - For non-square systems: only fall back if LS didn't converge
        //   (MaxIterations/RestorationFailed); if LS converged with high theta,
        //   the system is genuinely inconsistent.
        let ls_converged = matches!(
            ls_result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );
        if status == SolveStatus::LocalInfeasibility && (m == n || !ls_converged) {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: LS reformulation reports infeasibility (theta={:.4e}, ls_status={:?}), falling back to constrained IPM",
                    theta, ls_result.status
                );
            }
            let mut fallback_opts = options.clone();
            if options.max_wall_time > 0.0 {
                let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
                if remaining <= 0.1 {
                    return SolveResult {
                        x: ls_result.x,
                        objective: 0.0,
                        constraint_multipliers: vec![0.0; m],
                        bound_multipliers_lower: ls_result.bound_multipliers_lower,
                        bound_multipliers_upper: ls_result.bound_multipliers_upper,
                        constraint_values: g_out,
                        status: SolveStatus::MaxIterations,
                        iterations: ls_result.iterations,
                        diagnostics: SolverDiagnostics::default(),
                    };
                }
                fallback_opts.max_wall_time = remaining;
            }
            let ipm_result = solve_ipm(problem, &fallback_opts);
            if matches!(ipm_result.status, SolveStatus::Optimal | SolveStatus::Acceptable) {
                return ipm_result;
            }
            // IPM fallback failed too — try AL fallback for square NE systems
            // (e.g. HEART6 where both LS and constrained IPM fail)
            if options.enable_al_fallback {
                if options.print_level >= 5 {
                    eprintln!(
                        "ripopt: NE constrained IPM fallback failed ({:?}), trying AL",
                        ipm_result.status
                    );
                }
                let mut al_opts = options.clone();
                al_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
                if options.max_wall_time > 0.0 {
                    let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
                    if remaining <= 0.1 {
                        return ipm_result;
                    }
                    al_opts.max_wall_time = remaining;
                }
                let al_result = crate::augmented_lagrangian::solve(problem, &al_opts);
                if matches!(al_result.status, SolveStatus::Optimal | SolveStatus::Acceptable) {
                    if options.print_level >= 5 {
                        eprintln!(
                            "ripopt: NE AL fallback succeeded ({:?}, obj={:.6e})",
                            al_result.status, al_result.objective
                        );
                    }
                    return al_result;
                }
            }
            return ipm_result;
        }

        // L-BFGS retry on LS problem when IPM found a local min with nonzero residual
        let (final_x, final_status, final_g, final_iters, final_zl, final_zu) =
            if status == SolveStatus::LocalInfeasibility && options.enable_lbfgs_fallback {
                if options.print_level >= 5 {
                    eprintln!(
                        "ripopt: NE-to-LS LocalInfeasibility (theta={:.4e}), trying L-BFGS on LS",
                        theta
                    );
                }
                let lbfgs_ls = crate::lbfgs::solve(&ls_problem, options);
                let mut g_lb = vec![0.0; m];
                problem.constraints(&lbfgs_ls.x, &mut g_lb);
                let theta_lb = convergence::primal_infeasibility(&g_lb, &g_l, &g_u);

                if theta_lb < theta {
                    let new_status = if theta_lb < options.tol {
                        SolveStatus::Optimal
                    } else if theta_lb < options.acceptable_tol {
                        SolveStatus::Acceptable
                    } else {
                        SolveStatus::LocalInfeasibility
                    };
                    if options.print_level >= 5 {
                        eprintln!(
                            "ripopt: L-BFGS improved NE-to-LS (theta: {:.4e} -> {:.4e}, status={:?})",
                            theta, theta_lb, new_status
                        );
                    }
                    (lbfgs_ls.x, new_status, g_lb, lbfgs_ls.iterations,
                     lbfgs_ls.bound_multipliers_lower, lbfgs_ls.bound_multipliers_upper)
                } else {
                    if options.print_level >= 5 {
                        eprintln!(
                            "ripopt: L-BFGS did not improve NE-to-LS (theta_lb={:.4e} >= theta={:.4e})",
                            theta_lb, theta
                        );
                    }
                    (ls_result.x, status, g_out, ls_result.iterations,
                     ls_result.bound_multipliers_lower, ls_result.bound_multipliers_upper)
                }
            } else {
                (ls_result.x, status, g_out, ls_result.iterations,
                 ls_result.bound_multipliers_lower, ls_result.bound_multipliers_upper)
            };

        return SolveResult {
            x: final_x,
            objective: 0.0, // Original objective is f≡0
            constraint_multipliers: vec![0.0; m],
            bound_multipliers_lower: final_zl,
            bound_multipliers_upper: final_zu,
            constraint_values: final_g,
            status: final_status,
            iterations: final_iters,
            diagnostics: SolverDiagnostics::default(),
        };
    }

    // For unconstrained problems, try L-BFGS first (O(n·m) vs O(n³) per iteration)
    // then fall back to IPM if needed.
    let mut result = if options.enable_lbfgs_fallback && problem.num_constraints() == 0 {
        let lbfgs_result = crate::lbfgs::solve(problem, options);
        if matches!(lbfgs_result.status, SolveStatus::Optimal | SolveStatus::Acceptable) {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: L-BFGS solved unconstrained problem ({:?}, obj={:.6e})",
                    lbfgs_result.status, lbfgs_result.objective
                );
            }
            lbfgs_result
        } else {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: L-BFGS failed ({:?}, obj={:.6e}), trying IPM",
                    lbfgs_result.status, lbfgs_result.objective
                );
            }
            let ipm_result = solve_ipm(problem, options);
            if matches!(ipm_result.status, SolveStatus::Optimal | SolveStatus::Acceptable) {
                ipm_result
            } else if lbfgs_result.objective < ipm_result.objective {
                lbfgs_result
            } else {
                ipm_result
            }
        }
    } else {
        solve_ipm(problem, options)
    };

    // L-BFGS Hessian fallback: retry with L-BFGS Hessian approximation when
    // the exact-Hessian IPM failed. Skip if already in limited-memory mode.
    // Skip for large problems: L-BFGS forms an explicit dense n×n matrix (O(n²) memory).
    let n_lbfgs = problem.num_variables();
    if options.enable_lbfgs_hessian_fallback
        && !options.hessian_approximation_lbfgs
        && n_lbfgs <= 5000
        && matches!(
            result.status,
            SolveStatus::MaxIterations | SolveStatus::NumericalError | SolveStatus::RestorationFailed
        )
    {
        if options.print_level >= 5 {
            eprintln!(
                "ripopt: Exact-Hessian IPM failed ({:?}), retrying with L-BFGS Hessian approximation",
                result.status
            );
        }
        let mut lbfgs_opts = options.clone();
        lbfgs_opts.hessian_approximation_lbfgs = true;
        lbfgs_opts.enable_lbfgs_hessian_fallback = false; // prevent recursion
        // Cap fallback iterations
        lbfgs_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
        if options.max_wall_time > 0.0 {
            let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
            if remaining <= 0.1 {
                return result;
            }
            lbfgs_opts.max_wall_time = remaining;
        }
        let lbfgs_result = solve_ipm(problem, &lbfgs_opts);
        let lbfgs_better = matches!(
            lbfgs_result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );
        let current_solved = matches!(
            result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );
        if lbfgs_better && !current_solved {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: L-BFGS Hessian fallback succeeded ({:?}, obj={:.6e})",
                    lbfgs_result.status, lbfgs_result.objective
                );
            }
            result = lbfgs_result;
            result.diagnostics.fallback_used = Some("lbfgs_hessian".into());
        } else if options.print_level >= 5 {
            eprintln!(
                "ripopt: L-BFGS Hessian fallback did not improve ({:?})",
                lbfgs_result.status
            );
        }
    }

    // Augmented Lagrangian fallback for constrained problems
    // Skip for large problems: AL creates penalized subproblems that each run a full IPM solve,
    // which is prohibitively expensive at scale.
    let kkt_dim_al = problem.num_variables() + problem.num_constraints();
    if options.enable_al_fallback
        && problem.num_constraints() > 0
        && kkt_dim_al <= 5000
        && matches!(
            result.status,
            SolveStatus::MaxIterations | SolveStatus::NumericalError
        )
    {
        if options.print_level >= 5 {
            eprintln!(
                "ripopt: IPM failed ({:?}) on equality-only problem, trying AL fallback",
                result.status
            );
        }
        let mut al_opts = options.clone();
        // Cap fallback iterations — if it hasn't converged in 500 iters, it likely won't in 3000
        al_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
        if options.max_wall_time > 0.0 {
            let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
            if remaining <= 0.1 {
                return result;
            }
            al_opts.max_wall_time = remaining;
        }
        let al_result = crate::augmented_lagrangian::solve(problem, &al_opts);
        let al_better = matches!(
            al_result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );

        if al_better {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: AL fallback succeeded ({:?}, obj={:.6e})",
                    al_result.status, al_result.objective
                );
            }
            result = al_result;
            result.diagnostics.fallback_used = Some("augmented_lagrangian".into());
        } else if options.print_level >= 5 {
            eprintln!(
                "ripopt: AL fallback did not improve ({:?})",
                al_result.status
            );
        }
    }

    // SQP fallback: try SQP when IPM didn't reach Optimal and problem has constraints
    // Skip for large problems — SQP uses dense QP subproblems, O(n³) per iteration
    let kkt_dim = problem.num_variables() + problem.num_constraints();
    if options.enable_sqp_fallback
        && problem.num_constraints() > 0
        && kkt_dim <= 1000
        && !matches!(result.status, SolveStatus::Optimal)
    {
        if options.print_level >= 5 {
            eprintln!(
                "ripopt: Trying SQP fallback (result was {:?})",
                result.status
            );
        }
        let mut sqp_opts = options.clone();
        sqp_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
        if options.max_wall_time > 0.0 {
            let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
            if remaining <= 0.1 {
                return result;
            }
            sqp_opts.max_wall_time = remaining;
        }
        let sqp_result = crate::sqp::solve(problem, &sqp_opts);
        let sqp_better = matches!(
            sqp_result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );
        let current_solved = matches!(
            result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );
        let sqp_strictly_better = sqp_better
            && (!current_solved
                || matches!(sqp_result.status, SolveStatus::Optimal)
                    && !matches!(result.status, SolveStatus::Optimal)
                || sqp_result.status == result.status && sqp_result.objective < result.objective);
        if sqp_strictly_better {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: SQP fallback succeeded ({:?}, obj={:.6e})",
                    sqp_result.status, sqp_result.objective
                );
            }
            result = sqp_result;
            result.diagnostics.fallback_used = Some("sqp".into());
        } else if options.print_level >= 5 {
            eprintln!(
                "ripopt: SQP fallback did not improve ({:?})",
                sqp_result.status
            );
        }
    }

    // Slack variable fallback: if the solve didn't reach Optimal and the problem
    // has inequality constraints, retry with explicit slacks (g(x)-s=0, bounds on s).
    // Also try when only Acceptable was achieved — slack formulation may find a better solution
    // (e.g. ACOPR14 where AL gives Acceptable at a suboptimal local minimum).
    // For Acceptable results, only try slack if the problem is small enough
    // that the slack formulation (n + n_ineq variables) won't be too expensive.
    // For failure statuses, always try (the original behavior).
    let is_failure = matches!(
        result.status,
        SolveStatus::RestorationFailed | SolveStatus::MaxIterations | SolveStatus::NumericalError
    );
    let n_ineq = if !is_failure {
        let m = problem.num_constraints();
        let mut g_l = vec![0.0; m];
        let mut g_u = vec![0.0; m];
        problem.constraint_bounds(&mut g_l, &mut g_u);
        (0..m).filter(|&i| (g_l[i] - g_u[i]).abs() > 0.0).count()
    } else { 0 };
    let slack_too_large = !is_failure && problem.num_variables() + n_ineq > 200;

    if options.enable_slack_fallback
        && has_inequality_constraints(problem)
        && !matches!(result.status, SolveStatus::Optimal)
        && !slack_too_large
    {
        if options.print_level >= 5 {
            eprintln!(
                "ripopt: Trying slack fallback (result was {:?})",
                result.status
            );
        }

        let slack_prob = SlackFormulation::new(problem, &result.x);
        let mut slack_opts = options.clone();
        slack_opts.enable_slack_fallback = false; // prevent recursion
        // Cap fallback iterations
        slack_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
        if options.max_wall_time > 0.0 {
            let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
            if remaining <= 0.1 {
                return result;
            }
            slack_opts.max_wall_time = remaining;
        }

        let slack_result = solve_ipm(&slack_prob, &slack_opts);

        let slack_improved = matches!(
            slack_result.status,
            SolveStatus::Optimal | SolveStatus::Acceptable
        );

        // Only use slack result if it's strictly better than current result:
        // - Slack solved (Optimal/Acceptable) when current failed, OR
        // - Slack is Optimal when current is only Acceptable, OR
        // - Both are same status but slack has lower objective
        let current_solved = matches!(result.status, SolveStatus::Optimal | SolveStatus::Acceptable);
        let slack_strictly_better = slack_improved && (
            !current_solved  // any solve beats a failure
            || matches!(slack_result.status, SolveStatus::Optimal) && !matches!(result.status, SolveStatus::Optimal)
            || slack_result.status == result.status && slack_result.objective < result.objective
        );
        if slack_strictly_better {
            let n = problem.num_variables();
            let m = problem.num_constraints();

            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: Slack fallback succeeded ({:?}, obj={:.6e})",
                    slack_result.status, slack_result.objective
                );
            }

            // Extract original x from [x, s]
            let x_out = slack_result.x[..n].to_vec();

            // Evaluate original constraints at solution
            let mut g_out = vec![0.0; m];
            problem.constraints(&x_out, &mut g_out);

            // y from slack solve = original constraint multipliers
            let y_out = slack_result.constraint_multipliers;

            // z_l[..n], z_u[..n] = original bound multipliers
            let z_l_out = slack_result.bound_multipliers_lower[..n].to_vec();
            let z_u_out = slack_result.bound_multipliers_upper[..n].to_vec();

            let mut diag = slack_result.diagnostics;
            diag.fallback_used = Some("slack".into());
            diag.wall_time_secs = solve_start.elapsed().as_secs_f64();
            let slack_status = slack_result.status;
            let slack_iters = result.iterations + slack_result.iterations;
            if options.print_level >= 5 {
                diag.print_summary(slack_status, slack_iters);
            }
            return SolveResult {
                x: x_out,
                objective: slack_result.objective,
                constraint_multipliers: y_out,
                bound_multipliers_lower: z_l_out,
                bound_multipliers_upper: z_u_out,
                constraint_values: g_out,
                status: slack_status,
                iterations: slack_iters,
                diagnostics: diag,
            };
        } else if options.print_level >= 5 {
            eprintln!(
                "ripopt: Slack fallback did not improve ({:?}), returning original result",
                slack_result.status
            );
        }
    }

    result.diagnostics.wall_time_secs = solve_start.elapsed().as_secs_f64();
    if options.print_level >= 5 {
        result.diagnostics.print_summary(result.status, result.iterations);
    }
    result
}

/// Check if a problem has any inequality constraints (g_l[i] != g_u[i]).
fn has_inequality_constraints<P: NlpProblem>(problem: &P) -> bool {
    let m = problem.num_constraints();
    if m == 0 {
        return false;
    }
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    (0..m).any(|i| (g_l[i] - g_u[i]).abs() > 0.0)
}

/// Accumulates wall-clock time spent in each phase of the IPM loop.
/// Printed as a summary table at the end of `solve_ipm` when `print_level >= 5`.
struct PhaseTimings {
    problem_eval: Duration,
    kkt_assembly: Duration,
    factorization: Duration,
    direction_solve: Duration,
    line_search: Duration,
}

impl PhaseTimings {
    fn new() -> Self {
        PhaseTimings {
            problem_eval: Duration::ZERO,
            kkt_assembly: Duration::ZERO,
            factorization: Duration::ZERO,
            direction_solve: Duration::ZERO,
            line_search: Duration::ZERO,
        }
    }

    fn print_summary(&self, iterations: usize, total: Duration) {
        let total_secs = total.as_secs_f64();
        let phases = [
            ("Problem eval", self.problem_eval),
            ("KKT assembly", self.kkt_assembly),
            ("Factorization", self.factorization),
            ("Direction solve", self.direction_solve),
            ("Line search", self.line_search),
        ];
        let accounted: Duration = phases.iter().map(|(_, d)| *d).sum();
        let other = total.saturating_sub(accounted);

        eprintln!("\nPhase breakdown ({} iterations):", iterations);
        for (name, dur) in &phases {
            let secs = dur.as_secs_f64();
            let pct = if total_secs > 0.0 { 100.0 * secs / total_secs } else { 0.0 };
            eprintln!("  {:<20} {:>8.3}s ({:>5.1}%)", name, secs, pct);
        }
        let other_secs = other.as_secs_f64();
        let other_pct = if total_secs > 0.0 { 100.0 * other_secs / total_secs } else { 0.0 };
        eprintln!("  {:<20} {:>8.3}s ({:>5.1}%)", "Other", other_secs, other_pct);
        eprintln!("  {:<20} {:>8.3}s", "Total", total_secs);
    }
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

    // --- Linear constraint detection (on original unscaled problem for accuracy) ---
    let linear_constraints: Option<Vec<bool>> = if options.detect_linear_constraints && m_sc > 0 {
        let flags = crate::linearity::detect_linear_constraints(problem, &x0);
        let n_linear = flags.iter().filter(|&&f| f).count();
        if n_linear > 0 {
            if options.print_level >= 5 {
                eprintln!(
                    "ripopt: Detected {}/{} linear constraints (Hessian contribution skipped)",
                    n_linear, m_sc
                );
            }
            Some(flags)
        } else {
            None
        }
    } else {
        None
    };

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

    // L-BFGS-in-IPM mode
    let lbfgs_mode = options.hessian_approximation_lbfgs;
    let mut lbfgs_state = if lbfgs_mode {
        if options.print_level >= 5 {
            eprintln!("ripopt: Using L-BFGS Hessian approximation (limited-memory mode)");
        }
        Some(LbfgsIpmState::new(n))
    } else {
        None
    };

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

    // Initialize linear solver — use sparse for large KKT systems
    let use_sparse = (n + m) >= options.sparse_threshold;
    let mut lin_solver: Box<dyn LinearSolver> = if use_sparse {
        new_sparse_solver()
    } else {
        Box::new(DenseLdl::new())
    };
    let mut inertia_params = InertiaCorrectionParams::default();
    let mut restoration = RestorationPhase::new(500);

    // Estimate Schur complement density from Jacobian structure.
    // If J^T·D·J would be denser than the full augmented KKT system,
    // disable sparse condensed and use the full (n+m)×(n+m) system instead.
    let disable_sparse_condensed = if use_sparse && m > 0 {
        let (jac_rows_est, _) = problem.jacobian_structure();
        // Build row counts
        let mut row_nnz = vec![0usize; m];
        for &r in &jac_rows_est {
            row_nnz[r] += 1;
        }
        // Estimate Schur complement nnz: Σ k_i*(k_i+1)/2 (before dedup)
        let schur_nnz_upper: usize = row_nnz.iter().map(|&k| k * (k + 1) / 2).sum();
        // Augmented KKT nnz: hess_nnz + jac_nnz + n (diagonal)
        let (hess_rows_est, _) = problem.hessian_structure();
        let augmented_nnz = hess_rows_est.len() + jac_rows_est.len() + n;
        let disable = schur_nnz_upper > 2 * augmented_nnz;
        if disable && options.print_level >= 3 {
            eprintln!(
                "ripopt: Disabling sparse condensed KKT: Schur complement nnz estimate ({}) > 2× augmented KKT nnz ({})",
                schur_nnz_upper, augmented_nnz
            );
        }
        disable
    } else {
        false
    };

    // Initialize filter
    let mut filter = Filter::new(1e4);

    // Free/fixed mu mode state (replaces ad-hoc stall recovery)
    let mut mu_state = MuState::new();

    // Wall-clock time limit
    let start_time = Instant::now();

    // Phase timing instrumentation
    let mut timings = PhaseTimings::new();
    let ipm_start = Instant::now();

    // Watchdog mechanism state
    let mut consecutive_shortened: usize = 0;
    let mut watchdog_active: bool = false;
    let mut watchdog_trial_count: usize = 0;
    let mut watchdog_saved: Option<WatchdogSavedState> = None;

    // Constraint violation history for infeasibility detection
    let theta_history_len: usize = 100;
    let mut theta_history: Vec<f64> = Vec::with_capacity(theta_history_len);

    // Track whether the problem was ever feasible (theta < constr_viol_tol)
    // to prevent false infeasibility declarations on feasible problems.
    let mut ever_feasible = false;

    // Tiny step counter (Ipopt: accept full step when relative step < 10*eps for 2 consecutive)
    let mut consecutive_tiny_steps: usize = 0;

    // Overall progress stall detection: if neither primal nor dual infeasibility
    // improves by at least 1% over many consecutive iterations, terminate early.
    let mut stall_best_pr: f64 = f64::INFINITY;
    let mut stall_best_du: f64 = f64::INFINITY;
    let mut stall_no_progress_count: usize = 0;

    // Consecutive iterations with obj < -1e20 for robust unbounded detection
    let mut consecutive_unbounded: usize = 0;

    // Consecutive iterations where theta (primal infeasibility) stagnated.
    // Used by proactive infeasibility detection to exit early.
    let mut theta_stall_count: usize = 0;

    // Best feasible point tracking: save the best (lowest obj) point that is feasible
    let mut best_x: Option<Vec<f64>> = None;
    let mut best_obj: f64 = f64::INFINITY;
    let mut best_y: Option<Vec<f64>> = None;
    let mut best_z_l: Option<Vec<f64>> = None;
    let mut best_z_u: Option<Vec<f64>> = None;

    // Best-du point tracking
    let mut best_du_x: Option<Vec<f64>> = None;
    let mut best_du_val: f64 = f64::INFINITY;
    let mut best_du_y: Option<Vec<f64>> = None;
    let mut best_du_zl: Option<Vec<f64>> = None;
    let mut best_du_zu: Option<Vec<f64>> = None;

    // Dual stagnation detection: track best du improvement.
    // If du hasn't improved significantly over many iterations and we have a
    // best feasible point, restore it and restart with fresh parameters.
    let mut dual_stall_last_good_du: f64 = f64::INFINITY;
    let mut dual_stall_last_good_iter: usize = 0;
    let mut dual_stall_triggered: bool = false;

    // Initial evaluation
    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }

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
            state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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

    // Initialize constraint slack barrier multipliers v_l, v_u (Ipopt's v_L, v_U).
    // For each inequality constraint side: v = mu_init / slack.
    // This matches Ipopt's IpDefaultIterateInitializer.cpp.
    for i in 0..m {
        let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_eq {
            continue;
        }
        if state.g_l[i].is_finite() {
            let slack = (state.g[i] - state.g_l[i]).max(1e-20);
            state.v_l[i] = options.mu_init / slack;
        }
        if state.g_u[i].is_finite() {
            let slack = (state.g_u[i] - state.g[i]).max(1e-20);
            state.v_u[i] = options.mu_init / slack;
        }
    }

    // Set filter parameters based on initial constraint violation
    let theta_init = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_init);

    // Print header
    if options.print_level >= 5 {
        eprintln!(
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

    if options.print_level >= 5 {
        eprintln!("ripopt: Starting main loop (n={}, m={})", n, m);
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

        // --- Dual stagnation detection (runs every iteration, including restoration) ---
        // Track best du seen. If du hasn't improved for 500+ iterations and we have a
        // best feasible point, restore it with fresh filter/mu.
        // This catches problems like ACOPR14 where restoration cycling pushes the solver
        // off a good region and it gets stuck for thousands of iterations.
        if iteration > 0 {
            let current_du = convergence::dual_infeasibility(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &state.y, &state.z_l, &state.z_u, n,
            );
            if current_du < 0.5 * dual_stall_last_good_du {
                dual_stall_last_good_du = current_du;
                dual_stall_last_good_iter = iteration;
            }

            let stall_iters = iteration.saturating_sub(dual_stall_last_good_iter);
            if stall_iters >= 500
                && !dual_stall_triggered
                && current_du > options.acceptable_tol
                && best_x.is_some()
            {
                // Dual stagnation detected. Restore the best-du point (which had
                // du=best_du_val with stored x, y, z). This point was near-converged
                // but got disrupted by restoration cycling.
                if let Some(ref bdx) = best_du_x {
                    log::debug!(
                        "Dual stagnation at iter {}: du={:.2e}, restoring best-du point (du={:.2e} at iter {})",
                        iteration, current_du, dual_stall_last_good_du, dual_stall_last_good_iter
                    );
                    state.x.copy_from_slice(bdx);
                    if let Some(ref bdy) = best_du_y { state.y.copy_from_slice(bdy); }
                    if let Some(ref bdzl) = best_du_zl { state.z_l.copy_from_slice(bdzl); }
                    if let Some(ref bdzu) = best_du_zu { state.z_u.copy_from_slice(bdzu); }
                    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }

                    // Reset filter and bump mu for a fresh start from the good point.
                    filter.reset();
                    let theta_restart = state.constraint_violation();
                    filter.set_theta_min_from_initial(theta_restart);
                    state.mu = (state.mu * 100.0).max(1e-4).min(1e-1);
                    mu_state.mode = MuMode::Free;
                    mu_state.first_iter_in_mode = true;
                    mu_state.consecutive_restoration_failures = 0;
                    inertia_params.delta_w_last = 0.0;

                    // Check if the restored point meets acceptable convergence.
                    // The point had du=best_du_val which may already be excellent.
                    let rest_pr = state.constraint_violation();
                    let rest_du = convergence::dual_infeasibility(
                        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
                        &state.y, &state.z_l, &state.z_u, n,
                    );
                    let rest_co = convergence::complementarity_error(
                        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
                    );
                    // Use relaxed tolerances (acceptable level)
                    let s_max = 100.0_f64;
                    let mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
                        + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
                        + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
                    let s_d = if (m + 2 * n) > 0 {
                        (s_max.max(mult_sum / (m + 2 * n) as f64) / s_max).min(1e4)
                    } else { 1.0 };
                    let du_tol = (options.acceptable_tol * s_d).max(1e-2);
                    let co_tol = (options.acceptable_tol * s_d).max(1e-2);
                    let pr_tol = options.acceptable_tol.max(options.acceptable_constr_viol_tol);
                    if rest_pr <= pr_tol && rest_du <= du_tol && rest_co <= co_tol {
                        log::debug!(
                            "Restored best-du point passes acceptable (pr={:.2e}, du={:.2e}, co={:.2e})",
                            rest_pr, rest_du, rest_co
                        );
                        return make_result(&state, SolveStatus::Acceptable);
                    }

                    dual_stall_triggered = true;
                    // Fall through to normal iteration from good point
                }
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
                "{:>4} {:>14.7e} {:>10.2e} {:>10.2e} {:>10.2e} {:>10.2e} {:>8.2e} {:>8.2e}",
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
                if options.print_level >= 5 {
                    timings.print_summary(iteration + 1, ipm_start.elapsed());
                }
                return make_result(&state, SolveStatus::Optimal);
            }
            ConvergenceStatus::Acceptable => {
                if options.print_level >= 5 {
                    timings.print_summary(iteration + 1, ipm_start.elapsed());
                }
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

        // Track best-du point for cycling/stall detection at max_iter exit.
        if dual_inf < best_du_val {
            best_du_val = dual_inf;
            best_du_x = Some(state.x.clone());
            best_du_y = Some(state.y.clone());
            best_du_zl = Some(state.z_l.clone());
            best_du_zu = Some(state.z_u.clone());
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

        // Proactive infeasibility detection: if θ has stagnated for many consecutive
        // iterations AND the gradient of the violation is near-zero, declare infeasibility
        // earlier rather than burning iterations until restoration eventually fires.
        if options.proactive_infeasibility_detection
            && !ever_feasible
            && m > 0
            && iteration >= 50
            && primal_inf > options.constr_viol_tol
            && theta_history.len() >= theta_history_len
        {
            let theta_min_h = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
            let theta_max_h = theta_history.iter().cloned().fold(0.0f64, f64::max);
            // "Stagnated" = less than 1% relative variation over the history window
            if theta_max_h > 0.0 && (theta_max_h - theta_min_h) < 0.01 * primal_inf {
                theta_stall_count += 1;
            } else {
                theta_stall_count = 0;
            }
            // After 10 consecutive stagnation windows, check stationarity of ∇θ
            if theta_stall_count >= 10 {
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
                let stationarity_tol = 1e-3 * primal_inf.max(1.0);
                if grad_theta_norm < stationarity_tol {
                    log::info!(
                        "Proactive infeasibility at iter {}: θ stagnated at {:.2e}, ‖∇θ‖={:.2e}",
                        iteration, primal_inf, grad_theta_norm
                    );
                    return make_result(&state, SolveStatus::LocalInfeasibility);
                }
                // Stationarity not met — reset counter to check again in another window
                theta_stall_count = 0;
            }
        } else if ever_feasible {
            theta_stall_count = 0;
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

        // Overall progress stall detection: terminate when the solver is stuck
        // with no improvement in either primal or dual infeasibility.
        // Two triggers: (1) tiny steps for 15 iterations, (2) no metric improvement
        // for 30 iterations regardless of step size.
        // Only activate after 50 iterations to avoid tripping during early phases.
        if iteration > 50 {
            let pr_improved = primal_inf < 0.99 * stall_best_pr;
            let du_improved = dual_inf < 0.99 * stall_best_du;
            if pr_improved {
                stall_best_pr = primal_inf;
            }
            if du_improved {
                stall_best_du = dual_inf;
            }
            if pr_improved || du_improved {
                stall_no_progress_count = 0;
            } else {
                stall_no_progress_count += 1;
                let tiny_alpha = state.alpha_primal < 1e-8 && state.alpha_dual < 1e-4;
                // Terminate after 15 iters with truly negligible steps,
                // or 30 iters with no metric improvement regardless of step size
                let stall_limit = if tiny_alpha { 15 } else { 30 };
                if stall_no_progress_count >= stall_limit {
                    if options.print_level >= 3 {
                        eprintln!(
                            "ripopt: Stalled for {} iterations without progress (alpha_p={:.2e}, pr={:.2e}, du={:.2e}), terminating",
                            stall_no_progress_count, state.alpha_primal, primal_inf, dual_inf
                        );
                    }
                    return make_result(&state, SolveStatus::NumericalError);
                }
            }
        }

        // Compute sigma (barrier diagonal)
        let sigma = kkt::compute_sigma(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u);

        // Use dense condensed KKT (Schur complement) when m >= 2n and n is small.
        // Condensed cost is O(n^2*m + n^3) vs O((n+m)^3) — strictly better when m > n.
        let use_condensed = m >= 2 * n && n > 0 && !use_sparse;

        // Use sparse condensed KKT when problem is large and has constraints,
        // but only if the Schur complement is actually sparser than the augmented system.
        let use_sparse_condensed = use_sparse && m > 0 && !use_condensed && !disable_sparse_condensed;

        let t_kkt = Instant::now();
        let condensed_system = if use_condensed {
            Some(kkt::assemble_condensed_kkt(
                n, m,
                &state.hess_rows, &state.hess_cols, &state.hess_vals,
                &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                &state.y, &state.z_l, &state.z_u,
                &state.x, &state.x_l, &state.x_u, state.mu,
                &state.v_l, &state.v_u,
            ))
        } else {
            None
        };

        let sparse_condensed_system = if use_sparse_condensed {
            Some(kkt::assemble_sparse_condensed_kkt(
                n, m,
                &state.hess_rows, &state.hess_cols, &state.hess_vals,
                &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                &state.y, &state.z_l, &state.z_u,
                &state.x, &state.x_l, &state.x_u, state.mu,
                &state.v_l, &state.v_u,
            ))
        } else {
            None
        };

        // Build full KKT only when not using any condensed path
        let mut kkt_system_opt: Option<kkt::KktSystem> = if !use_condensed && !use_sparse_condensed {
            Some(kkt::assemble_kkt(
                n, m,
                &state.hess_rows, &state.hess_cols, &state.hess_vals,
                &state.jac_rows, &state.jac_cols, &state.jac_vals,
                &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                &state.y, &state.z_l, &state.z_u,
                &state.x, &state.x_l, &state.x_u, state.mu,
                use_sparse, &state.v_l, &state.v_u,
            ))
        } else {
            None
        };
        timings.kkt_assembly += t_kkt.elapsed();

        // On first iteration with sparse condensed, detect bandwidth for the condensed system
        if iteration == 0 && use_sparse_condensed {
            if let Some(ref sc) = sparse_condensed_system {
                let bw = BandedLdl::compute_bandwidth(&sc.matrix.triplet_rows, &sc.matrix.triplet_cols);
                if bw * bw <= n {
                    if options.print_level >= 5 {
                        eprintln!("ripopt: Sparse condensed S has bandwidth {} for n={}, using banded solver", bw, n);
                    }
                    lin_solver = Box::new(BandedLdl::new());
                } else if options.print_level >= 5 {
                    eprintln!("ripopt: Sparse condensed S has bandwidth {} for n={}, using sparse solver", bw, n);
                }
            }
        }

        // Factor with inertia correction (only for non-condensed path)
        if let Some(ref mut kkt_system) = kkt_system_opt {
        let t_fact = Instant::now();
        let inertia_result =
            kkt::factor_with_inertia_correction(kkt_system, lin_solver.as_mut(), &mut inertia_params);
        timings.factorization += t_fact.elapsed();

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
                    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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
                            use_sparse, &state.v_l, &state.v_u,
                        );
                        if kkt::factor_with_inertia_correction(
                            &mut kkt_p, lin_solver.as_mut(), &mut inertia_params,
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
                    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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
                Some(&|x_eval: &[f64]| problem.objective(x_eval)),
            );
            if success {
                state.x = x_rest;
                state.alpha_primal = 0.0;
                state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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
                state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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
        } // end: if let Some(ref mut kkt_system) = kkt_system_opt

        // Solve for search direction
        let t_dir = Instant::now();
        let mut cond_solver_for_soc: Option<DenseLdl> = None;
        let (dx, dy) = if let Some(ref cond) = condensed_system {
            // Try condensed solve first (faster for m >> n)
            let mut cond_solver = DenseLdl::new();
            let cond_ok = cond_solver.bunch_kaufman_factor(&cond.matrix).is_ok();
            let cond_result = if cond_ok {
                kkt::solve_condensed(cond, &mut cond_solver).ok()
            } else {
                None
            };

            if let Some(d) = cond_result {
                cond_solver_for_soc = Some(cond_solver);
                d
            } else {
                // Condensed failed — build full KKT on demand
                let mut kkt = kkt::assemble_kkt(
                    n, m,
                    &state.hess_rows, &state.hess_cols, &state.hess_vals,
                    &state.jac_rows, &state.jac_cols, &state.jac_vals,
                    &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                    &state.y, &state.z_l, &state.z_u,
                    &state.x, &state.x_l, &state.x_u, state.mu,
                    use_sparse, &state.v_l, &state.v_u,
                );
                if kkt::factor_with_inertia_correction(
                    &mut kkt, lin_solver.as_mut(), &mut inertia_params,
                ).is_err() {
                    let (x_rest, success) = restoration.restore(
                        &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
                        &state.jac_rows, &state.jac_cols, n, m, options,
                        &|theta, phi| filter.is_acceptable(theta, phi),
                        &|x_eval, g_out| problem.constraints(x_eval, g_out),
                        &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
                        Some(&|x_eval: &[f64]| problem.objective(x_eval)),
                    );
                    if success {
                        state.x = x_rest;
                        state.alpha_primal = 0.0;
                        state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
                        continue;
                    }
                    return make_result(&state, SolveStatus::NumericalError);
                }
                match kkt::solve_for_direction(&kkt, lin_solver.as_mut()) {
                    Ok(d) => {
                        kkt_system_opt = Some(kkt);
                        d
                    },
                    Err(e) => {
                        log::warn!("KKT solve failed: {}", e);
                        let (x_rest, success) = restoration.restore(
                            &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
                            &state.jac_rows, &state.jac_cols, n, m, options,
                            &|theta, phi| filter.is_acceptable(theta, phi),
                            &|x_eval, g_out| problem.constraints(x_eval, g_out),
                            &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
                            Some(&|x_eval: &[f64]| problem.objective(x_eval)),
                        );
                        if success {
                            state.x = x_rest;
                            state.alpha_primal = 0.0;
                            state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
                            continue;
                        }
                        return make_result(&state, SolveStatus::NumericalError);
                    }
                }
            }
        } else if let Some(ref sc) = sparse_condensed_system {
            // Sparse condensed path: factor S = H + Σ + J^T·D_c^{-1}·J with banded/sparse solver
            let kkt_sc = KktMatrix::Sparse(sc.matrix.clone());
            let factor_ok = lin_solver.factor(&kkt_sc).is_ok();
            if factor_ok {
                match kkt::solve_sparse_condensed(sc, lin_solver.as_mut()) {
                    Ok(d) => d,
                    Err(_) => {
                        // Fall back to full KKT
                        let mut kkt = kkt::assemble_kkt(
                            n, m,
                            &state.hess_rows, &state.hess_cols, &state.hess_vals,
                            &state.jac_rows, &state.jac_cols, &state.jac_vals,
                            &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                            &state.y, &state.z_l, &state.z_u,
                            &state.x, &state.x_l, &state.x_u, state.mu,
                            use_sparse, &state.v_l, &state.v_u,
                        );
                        let mut fallback_solver = new_fallback_solver(use_sparse);
                        if kkt::factor_with_inertia_correction(
                            &mut kkt, fallback_solver.as_mut(), &mut inertia_params,
                        ).is_ok() {
                            kkt::solve_for_direction(&kkt, fallback_solver.as_mut())
                                .unwrap_or_else(|_| gradient_descent_fallback(&state)
                                    .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m])))
                        } else {
                            gradient_descent_fallback(&state)
                                .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m]))
                        }
                    }
                }
            } else {
                // Factor failed — try full KKT
                let mut kkt = kkt::assemble_kkt(
                    n, m,
                    &state.hess_rows, &state.hess_cols, &state.hess_vals,
                    &state.jac_rows, &state.jac_cols, &state.jac_vals,
                    &sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
                    &state.y, &state.z_l, &state.z_u,
                    &state.x, &state.x_l, &state.x_u, state.mu,
                    use_sparse, &state.v_l, &state.v_u,
                );
                let mut fallback_solver = new_fallback_solver(use_sparse);
                if kkt::factor_with_inertia_correction(
                    &mut kkt, fallback_solver.as_mut(), &mut inertia_params,
                ).is_ok() {
                    kkt::solve_for_direction(&kkt, fallback_solver.as_mut())
                        .unwrap_or_else(|_| gradient_descent_fallback(&state)
                            .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m])))
                } else {
                    gradient_descent_fallback(&state)
                        .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m]))
                }
            }
        } else {
            // Mehrotra predictor-corrector (PC): probe the affine-scaling direction to
            // estimate a better barrier parameter μ before solving the main step.
            //
            // Algorithm:
            //   1. Solve affine predictor (μ=0 in RHS) — same factored matrix, new RHS.
            //   2. Compute max step α_aff for the predictor using fraction-to-boundary.
            //   3. Compute μ_aff = average complementarity after the affine step.
            //   4. Set σ = (μ_aff / μ)³  (centering parameter).
            //   5. Update KKT RHS to use μ_new = σ·μ (≤ μ → more aggressive decrease).
            //   6. Solve the main corrector with the improved RHS.
            //
            // Cost: one extra triangular solve (no re-factorization).
            // Reference: Mehrotra (1992, SIAM J. Optim.); Nocedal/Wächter (2006).
            if options.mehrotra_pc {
                let has_bounds = (0..n).any(|i| state.x_l[i].is_finite() || state.x_u[i].is_finite());
                if has_bounds {
                    // Scope the immutable borrow so it ends before the mutable update below.
                    let pc_rhs: Option<Vec<f64>> = {
                        let kkt = kkt_system_opt.as_ref().unwrap();
                        let rhs_aff = kkt::affine_predictor_rhs(
                            &kkt.rhs, &state.x, &state.x_l, &state.x_u, state.mu,
                        );
                        if let Ok((dx_aff, _)) = kkt::solve_with_custom_rhs(
                            kkt.n, kkt.dim, lin_solver.as_mut(), &rhs_aff,
                        ) {
                            // Complementarity steps for the affine predictor (μ=0)
                            let (dz_l_aff, dz_u_aff) = kkt::recover_dz(
                                &state.x, &state.x_l, &state.x_u,
                                &state.z_l, &state.z_u, &dx_aff, 0.0,
                            );
                            // Compute α_aff = max step along affine direction
                            let tau_aff = 1.0 - 1e-3;
                            let aff_zl = filter::fraction_to_boundary(&state.z_l, &dz_l_aff, tau_aff);
                            let aff_zu = filter::fraction_to_boundary(&state.z_u, &dz_u_aff, tau_aff);
                            let mut alpha_aff = aff_zl.min(aff_zu).min(1.0);
                            for i in 0..n {
                                if state.x_l[i].is_finite() && dx_aff[i] < 0.0 {
                                    let s = (state.x[i] - state.x_l[i]).max(1e-20);
                                    alpha_aff = alpha_aff.min(tau_aff * s / (-dx_aff[i]));
                                }
                                if state.x_u[i].is_finite() && dx_aff[i] > 0.0 {
                                    let s = (state.x_u[i] - state.x[i]).max(1e-20);
                                    alpha_aff = alpha_aff.min(tau_aff * s / dx_aff[i]);
                                }
                            }
                            alpha_aff = alpha_aff.clamp(0.0, 1.0);
                            // Compute μ_aff = average complementarity after affine step
                            let mut mu_aff_sum = 0.0_f64;
                            let mut nb: usize = 0;
                            for i in 0..n {
                                if state.x_l[i].is_finite() {
                                    let s = (state.x[i] + alpha_aff * dx_aff[i]
                                        - state.x_l[i]).max(1e-20);
                                    let z = (state.z_l[i] + alpha_aff * dz_l_aff[i]).max(1e-20);
                                    mu_aff_sum += s * z;
                                    nb += 1;
                                }
                                if state.x_u[i].is_finite() {
                                    let s = (state.x_u[i] - state.x[i]
                                        - alpha_aff * dx_aff[i]).max(1e-20);
                                    let z = (state.z_u[i] + alpha_aff * dz_u_aff[i]).max(1e-20);
                                    mu_aff_sum += s * z;
                                    nb += 1;
                                }
                            }
                            if nb > 0 {
                                let mu_aff = mu_aff_sum / nb as f64;
                                let sigma = (mu_aff / state.mu).powi(3).clamp(0.0, 1.0);
                                let mu_pc = (sigma * state.mu).max(options.mu_min);
                                // Apply PC only when the centering parameter suggests
                                // a meaningful μ decrease (σ < 0.95) and the probe is
                                // not degenerate. Skip early iterations to avoid
                                // amplifying noise at poorly-scaled starting points.
                                // Also skip when sigma is very small (near convergence,
                                // α_aff → 1) to avoid over-aggressive barrier decrease.
                                let sigma_skip_min = 0.05_f64;
                                if mu_pc < state.mu * 0.95 && sigma >= sigma_skip_min && iteration >= 2 {
                                    log::debug!(
                                        "Mehrotra PC iter {}: σ={:.4} α_aff={:.4} μ: {:.2e}→{:.2e}",
                                        iteration, sigma, alpha_aff, state.mu, mu_pc
                                    );
                                    Some(kkt::rebuild_rhs_with_mu(
                                        &kkt.rhs, &state.x, &state.x_l, &state.x_u,
                                        state.mu, mu_pc,
                                    ))
                                } else {
                                    None
                                }
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    }; // immutable borrow of kkt_system_opt ends here

                    // Apply the improved RHS to the KKT system
                    if let Some(new_rhs) = pc_rhs {
                        kkt_system_opt.as_mut().unwrap().rhs = new_rhs;
                    }
                }
            }

            let dir_result = kkt::solve_for_direction(kkt_system_opt.as_ref().unwrap(), lin_solver.as_mut());
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
                            Some(&|x_eval: &[f64]| problem.objective(x_eval)),
                        );
                        if success {
                            state.x = x_rest;
                            state.alpha_primal = 0.0;
                            state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
                            continue;
                        }
                        return make_result(&state, SolveStatus::NumericalError);
                    }
                }
            }
        };

        timings.direction_solve += t_dir.elapsed();

        // Recover bound multiplier steps
        let (dz_l, dz_u) =
            kkt::recover_dz(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, &dx, state.mu);

        state.dx = dx;
        state.dy = dy;
        state.dz_l = dz_l;
        state.dz_u = dz_u;

        // Gondzio multiple centrality corrections (MCC).
        //
        // After computing the main search direction (possibly Mehrotra-corrected),
        // perform up to `gondzio_mcc_max` additional centrality corrections.
        // Each correction uses the SAME factored KKT matrix (one extra backsolve each)
        // to drive complementarity pairs that are far from μ back toward the central path.
        //
        // Acceptance criterion: the correction is accepted only if it does not reduce
        // the maximum step length by more than 10%.
        //
        // Reference: Gondzio (1994, Comput. Optim. Appl.); Gondzio (2007).
        if options.gondzio_mcc_max > 0 {
            if let Some(ref kkt) = kkt_system_opt {
                // Compute a preliminary max step for the current direction
                let tau_mcc = if mu_state.mode == MuMode::Free {
                    let nlp_error = primal_inf + dual_inf + compl_inf_best;
                    (1.0 - nlp_error).max(options.tau_min)
                } else {
                    (1.0 - state.mu).max(options.tau_min)
                };
                let mcc_zl = filter::fraction_to_boundary(&state.z_l, &state.dz_l, tau_mcc);
                let mcc_zu = filter::fraction_to_boundary(&state.z_u, &state.dz_u, tau_mcc);
                let mut alpha_mcc = mcc_zl.min(mcc_zu).min(1.0);
                for i in 0..n {
                    if state.x_l[i].is_finite() && state.dx[i] < 0.0 {
                        let s = state.x[i] - state.x_l[i];
                        alpha_mcc = alpha_mcc.min(tau_mcc * s / (-state.dx[i]));
                    }
                    if state.x_u[i].is_finite() && state.dx[i] > 0.0 {
                        let s = state.x_u[i] - state.x[i];
                        alpha_mcc = alpha_mcc.min(tau_mcc * s / state.dx[i]);
                    }
                }
                alpha_mcc = alpha_mcc.clamp(0.0, 1.0);

                let mu_target = state.mu;
                let beta_min = 0.01_f64;  // centrality lower bound: z·s ≥ β_min·μ
                let beta_max = 100.0_f64; // centrality upper bound: z·s ≤ β_max·μ

                for _mcc_iter in 0..options.gondzio_mcc_max {
                    // Build centrality correction RHS: target z·s → μ for outliers
                    let mut rhs_mcc = vec![0.0_f64; kkt.dim];
                    let mut needs_correction = false;

                    for i in 0..n {
                        if state.x_l[i].is_finite() {
                            let s_t = (state.x[i] + alpha_mcc * state.dx[i]
                                - state.x_l[i]).max(1e-20);
                            let z_t = (state.z_l[i] + alpha_mcc * state.dz_l[i]).max(1e-20);
                            let c = z_t * s_t;
                            if c < beta_min * mu_target || c > beta_max * mu_target {
                                rhs_mcc[i] += (mu_target - c) / s_t;
                                needs_correction = true;
                            }
                        }
                        if state.x_u[i].is_finite() {
                            let s_t = (state.x_u[i] - state.x[i]
                                - alpha_mcc * state.dx[i]).max(1e-20);
                            let z_t = (state.z_u[i] + alpha_mcc * state.dz_u[i]).max(1e-20);
                            let c = z_t * s_t;
                            if c < beta_min * mu_target || c > beta_max * mu_target {
                                rhs_mcc[i] -= (mu_target - c) / s_t;
                                needs_correction = true;
                            }
                        }
                    }

                    if !needs_correction {
                        break;
                    }

                    // Solve for the centrality correction direction
                    match kkt::solve_with_custom_rhs(kkt.n, kkt.dim, lin_solver.as_mut(), &rhs_mcc) {
                        Ok((ddx, ddy)) => {
                            // Compute bound-multiplier corrections from the Newton step:
                            //   S_l · ddz_l + Z_l · ddx = 0  (no centering in correction)
                            //   ddz_l[i] = -(z_l[i] / s_l[i]) * ddx[i]
                            //   ddz_u[i] =  (z_u[i] / s_u[i]) * ddx[i]
                            // NOTE: do NOT use recover_dz(mu=0) here — that adds the
                            // affine centering term (-z_l[i]) which would drive z_l to zero.
                            let mut ddz_l = vec![0.0_f64; n];
                            let mut ddz_u = vec![0.0_f64; n];
                            for i in 0..n {
                                if state.x_l[i].is_finite() {
                                    let s_l = (state.x[i] - state.x_l[i]).max(1e-20);
                                    ddz_l[i] = -(state.z_l[i] / s_l) * ddx[i];
                                }
                                if state.x_u[i].is_finite() {
                                    let s_u = (state.x_u[i] - state.x[i]).max(1e-20);
                                    ddz_u[i] = (state.z_u[i] / s_u) * ddx[i];
                                }
                            }

                            // Tentatively update direction
                            let mut dx_c: Vec<f64> = state.dx.iter().zip(ddx.iter()).map(|(a, b)| a + b).collect();
                            let dy_c: Vec<f64> = state.dy.iter().zip(ddy.iter()).map(|(a, b)| a + b).collect();
                            let dz_l_c: Vec<f64> = state.dz_l.iter().zip(ddz_l.iter()).map(|(a, b)| a + b).collect();
                            let dz_u_c: Vec<f64> = state.dz_u.iter().zip(ddz_u.iter()).map(|(a, b)| a + b).collect();

                            // Compute new alpha for the corrected direction
                            let new_zl = filter::fraction_to_boundary(&state.z_l, &dz_l_c, tau_mcc);
                            let new_zu = filter::fraction_to_boundary(&state.z_u, &dz_u_c, tau_mcc);
                            let mut alpha_new = new_zl.min(new_zu).min(1.0);
                            for i in 0..n {
                                if state.x_l[i].is_finite() && dx_c[i] < 0.0 {
                                    let s = state.x[i] - state.x_l[i];
                                    alpha_new = alpha_new.min(tau_mcc * s / (-dx_c[i]));
                                }
                                if state.x_u[i].is_finite() && dx_c[i] > 0.0 {
                                    let s = state.x_u[i] - state.x[i];
                                    alpha_new = alpha_new.min(tau_mcc * s / dx_c[i]);
                                }
                            }
                            alpha_new = alpha_new.clamp(0.0, 1.0);

                            // Accept only if the correction doesn't shrink alpha by >10%
                            if alpha_new >= 0.9 * alpha_mcc {
                                // dx_c needs explicit type annotation for vec addition
                                let _ = &mut dx_c; // silence unused_mut warning
                                state.dx = dx_c;
                                state.dy = dy_c;
                                state.dz_l = dz_l_c;
                                state.dz_u = dz_u_c;
                                alpha_mcc = alpha_new;
                                log::debug!(
                                    "Gondzio MCC iter {}: correction accepted, α_mcc={:.4}",
                                    iteration, alpha_mcc
                                );
                            } else {
                                break;
                            }
                        }
                        Err(_) => break,
                    }
                }
            }
        }

        // Compute maximum step sizes using fraction-to-boundary rule.
        // Free mode: tau based on NLP error. Fixed mode: tau based on mu.
        let tau = if mu_state.mode == MuMode::Free {
            let nlp_error = primal_inf + dual_inf + compl_inf_best;
            (1.0 - nlp_error).max(options.tau_min)
        } else {
            (1.0 - state.mu).max(options.tau_min)
        };

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
        {
            let max_rel_step: f64 = (0..n)
                .map(|i| (alpha_primal_max * state.dx[i]).abs() / (state.x[i].abs() + 1.0))
                .fold(0.0f64, f64::max);
            if max_rel_step < 1e-14 && primal_inf < 1e-4 {
                consecutive_tiny_steps += 1;
                mu_state.tiny_step = true;
                if consecutive_tiny_steps >= 2 {
                    // Force mu decrease (Ipopt: monotone decrease on tiny step)
                    let new_mu = (options.mu_linear_decrease_factor * state.mu)
                        .min(state.mu.powf(options.mu_superlinear_decrease_power))
                        .max(options.mu_min);
                    if (new_mu - state.mu).abs() < 1e-20 {
                        log::debug!("Tiny step with mu at minimum, checking acceptability");
                    } else {
                        state.mu = new_mu;
                        filter.reset();
                        let theta_new = state.constraint_violation();
                        filter.set_theta_min_from_initial(theta_new);
                        log::debug!("Tiny step detected, forced mu decrease to {:.2e}", state.mu);
                    }
                    consecutive_tiny_steps = 0;
                }
            } else {
                consecutive_tiny_steps = 0;
                mu_state.tiny_step = false;
            }
        }

        // Line search
        let t_ls = Instant::now();
        let theta_current = primal_inf;
        let phi_current = state.barrier_objective(options);
        let grad_phi_step = state.barrier_directional_derivative(options);

        let mut alpha = alpha_primal_max;
        let mut step_accepted = false;
        let min_alpha = filter.compute_alpha_min(theta_current, grad_phi_step);

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

            // Second-order correction (SOC) — try on every backtracking step where theta increases
            if theta_trial > theta_current && options.max_soc > 0 {
                let soc_accepted = if let (Some(ref cond), Some(ref mut cs)) = (&condensed_system, &mut cond_solver_for_soc) {
                    // Use condensed SOC (avoids building full KKT)
                    attempt_soc_condensed(
                        &state, problem, &g_trial, cs, cond, &filter,
                        theta_current, phi_current, grad_phi_step, alpha, options,
                    )
                } else if let Some(ref sc) = sparse_condensed_system {
                    // Use sparse condensed SOC
                    attempt_soc_sparse_condensed(
                        &state, problem, &g_trial, lin_solver.as_mut(), sc, &filter,
                        theta_current, phi_current, grad_phi_step, alpha, options,
                    )
                } else if let Some(ref kkt) = kkt_system_opt {
                    attempt_soc(
                        &state, problem, &x_trial, &g_trial,
                        lin_solver.as_mut(), kkt, &filter,
                        theta_current, phi_current, grad_phi_step, alpha, options,
                    )
                } else {
                    None
                };

                if let Some((x_soc, obj_soc, g_soc, alpha_soc)) = soc_accepted {
                    state.diagnostics.soc_corrections += 1;
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
            state.diagnostics.filter_rejects += 1;

            // Add current point to filter before entering restoration (Ipopt convention).
            filter.add(theta_current, phi_current);
            filter.augment_for_restoration(theta_current);

            // Phase 1: Fast GN restoration
            log::debug!("Line search failed at iteration {}, entering restoration", iteration);

            let (x_rest, gn_success) = restoration.restore(
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
                Some(&|x_eval: &[f64]| problem.objective(x_eval)),
            );

            if gn_success {
                state.diagnostics.restoration_count += 1;
                // GN restoration succeeded — apply standard restoration success handling
                apply_restoration_success(
                    &mut state, &mut filter, &mut mu_state, options, n, m, problem, &x_rest,
                    linear_constraints.as_deref(), lbfgs_mode, &mut lbfgs_state,
                );
                continue;
            }

            // GN restoration failed — recovery logic with NLP restoration as last resort
            {
                mu_state.consecutive_restoration_failures += 1;
                let fail_count = mu_state.consecutive_restoration_failures;

                // At fail_count == 2: try full NLP restoration early.
                // The NLP restoration is the most robust approach (Ipopt's primary method).
                // Try it before exhausting simpler recovery strategies.
                // Skip for large problems: NLP restoration doubles the problem size,
                // making it prohibitively expensive for n+m > 10000.
                let kkt_dim = n + m;
                if fail_count == 2 && !options.disable_nlp_restoration && kkt_dim <= 10000 {
                    state.diagnostics.nlp_restoration_count += 1;
                    let (x_nlp, outcome) = attempt_nlp_restoration(
                        problem, &state, &filter, options, theta_current,
                    );
                    match outcome {
                        RestorationOutcome::Success => {
                            apply_restoration_success(
                                &mut state, &mut filter, &mut mu_state, options, n, m,
                                problem, &x_nlp,
                                linear_constraints.as_deref(), lbfgs_mode, &mut lbfgs_state,
                            );
                            continue;
                        }
                        RestorationOutcome::LocalInfeasibility
                        | RestorationOutcome::Failed => {
                            // Fall through to continue recovery.
                            // Don't immediately return LocalInfeasibility — the existing
                            // infeasibility detection at fail_count > 6 uses stationarity
                            // checks which are more reliable.
                        }
                    }
                }

                // For large problems (no NLP restoration), give up sooner
                let max_restore_attempts = if kkt_dim > 10000 { 3 } else { 6 };
                if fail_count > max_restore_attempts {
                    // Exhausted recovery attempts: check infeasibility and give up
                    log::warn!("Restoration failed at iteration {} (attempt #{})", iteration, fail_count);
                    let current_theta = state.constraint_violation();

                    // Check stationarity of violation
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
                        let mut grad_theta = vec![0.0; n];
                        for (idx, (&row, &col)) in
                            state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
                        {
                            grad_theta[col] += state.jac_vals[idx] * violation[row];
                        }
                        let grad_theta_norm = grad_theta.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
                        let stationarity_tol = 1e-4 * current_theta.max(1.0);
                        if grad_theta_norm < stationarity_tol {
                            log::info!(
                                "Local infeasibility detected: theta={:.2e}, ||∇theta||={:.2e}",
                                current_theta, grad_theta_norm
                            );
                            return make_result(&state, SolveStatus::LocalInfeasibility);
                        }
                    }

                    if !ever_feasible && current_theta > 1e4 && iteration > 500 && theta_history.len() >= theta_history_len {
                        let min_theta = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
                        if current_theta > 0.01 * min_theta {
                            return make_result(&state, SolveStatus::Infeasible);
                        }
                    }
                    return make_result(&state, SolveStatus::RestorationFailed);
                }

                // Recovery strategies: cycle through mode switches and mu perturbations
                log::debug!("Restoration failed (attempt #{}), trying recovery", fail_count);
                let mu_factors: [f64; 6] = [10.0, 0.1, 100.0, 0.01, 1000.0, 0.001];

                match fail_count {
                    1 => {
                        // First failure: switch mode
                        if mu_state.mode == MuMode::Free {
                            state.diagnostics.mu_mode_switches += 1;
                            mu_state.mode = MuMode::Fixed;
                            mu_state.first_iter_in_mode = true;
                            let avg_compl = compute_avg_complementarity(&state);
                            if avg_compl > 0.0 {
                                state.mu = (options.adaptive_mu_monotone_init_factor * avg_compl)
                                    .clamp(options.mu_min, 1e5);
                            }
                        } else {
                            // Force mu decrease
                            let new_mu = (options.mu_linear_decrease_factor * state.mu)
                                .min(state.mu.powf(options.mu_superlinear_decrease_power))
                                .max(options.mu_min);
                            state.mu = new_mu;
                        }
                    }
                    _ => {
                        // Subsequent failures: try varied mu perturbation
                        let factor = mu_factors[(fail_count - 2) % mu_factors.len()];
                        state.mu = (state.mu * factor).max(options.mu_min).min(1e5);
                    }
                }
                filter.reset();
                let theta_now = state.constraint_violation();
                filter.set_theta_min_from_initial(theta_now);
                inertia_params.delta_w_last = 0.0;

                // On attempts 3+: also perturb x to escape current basin
                if fail_count >= 3 {
                    for i in 0..n {
                        let range = if state.x_l[i].is_finite() && state.x_u[i].is_finite() {
                            state.x_u[i] - state.x_l[i]
                        } else {
                            state.x[i].abs().max(1.0)
                        };
                        let sign = if (i * 7 + fail_count * 13) % 3 == 0 { -1.0 } else { 1.0 };
                        state.x[i] += sign * 1e-4 * range;
                        if state.x_l[i].is_finite() {
                            state.x[i] = state.x[i].max(state.x_l[i] + 1e-14);
                        }
                        if state.x_u[i].is_finite() {
                            state.x[i] = state.x[i].min(state.x_u[i] - 1e-14);
                        }
                    }
                    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
                }
                continue;
            }
        }

        // Step was accepted — reset consecutive restoration failure counter
        mu_state.consecutive_restoration_failures = 0;

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
            state.diagnostics.watchdog_activations += 1;
            watchdog_active = true;
            watchdog_trial_count = 0;
            let wd_theta = state.constraint_violation();
            let wd_phi = state.barrier_objective(options);
            watchdog_saved = Some(WatchdogSavedState {
                x: state.x.clone(),
                y: state.y.clone(),
                z_l: state.z_l.clone(),
                z_u: state.z_u.clone(),
                v_l: state.v_l.clone(),
                v_u: state.v_u.clone(),
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
                    state.v_l = saved.v_l.clone();
                    state.v_u = saved.v_u.clone();
                    state.mu = saved.mu;
                    state.obj = saved.obj;
                    state.g = saved.g.clone();
                    state.grad_f = saved.grad_f.clone();
                    state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }

                    watchdog_active = false;
                    watchdog_trial_count = 0;
                    watchdog_saved = None;
                    continue;
                }
            }
        }

        timings.line_search += t_ls.elapsed();

        // Update dual variables
        let alpha_d = alpha_dual_max;
        for i in 0..m {
            state.y[i] += alpha_d * state.dy[i];
        }
        // Ipopt kappa_sigma safeguard: keep z*s in [mu_ks/kappa_sigma, kappa_sigma*mu_ks]
        // In free mode, use avg_compl (clamped to [mu, 1e3]) for more stable z bounds
        let kappa_sigma = 1e10;
        let mu_ks = if mu_state.mode == MuMode::Free {
            compute_avg_complementarity(&state).max(state.mu).min(1e3)
        } else {
            state.mu
        };
        for i in 0..n {
            if state.x_l[i].is_finite() {
                let z_new = (state.z_l[i] + alpha_d * state.dz_l[i]).max(1e-20);
                let s_l = (state.x[i] - state.x_l[i]).max(1e-20);
                let z_lo = mu_ks / (kappa_sigma * s_l);
                let z_hi = kappa_sigma * mu_ks / s_l;
                state.z_l[i] = z_new.clamp(z_lo, z_hi);
            }
            if state.x_u[i].is_finite() {
                let z_new = (state.z_u[i] + alpha_d * state.dz_u[i]).max(1e-20);
                let s_u = (state.x_u[i] - state.x[i]).max(1e-20);
                let z_lo = mu_ks / (kappa_sigma * s_u);
                let z_hi = kappa_sigma * mu_ks / s_u;
                state.z_u[i] = z_new.clamp(z_lo, z_hi);
            }
        }

        state.alpha_dual = alpha_d;

        // Re-evaluate at new point
        let t_eval = Instant::now();
        state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
        timings.problem_eval += t_eval.elapsed();

        // Reset v_l, v_u from barrier equilibrium v = mu_ks / slack.
        // Simple reset rather than Newton update (our dv is approximate since we
        // lack explicit slacks, and FTB on v can restrict alpha_d too much).
        for i in 0..m {
            if state.v_l[i] > 0.0 && state.g_l[i].is_finite() {
                let slack = (state.g[i] - state.g_l[i]).max(1e-20);
                state.v_l[i] = mu_ks / slack;
            }
            if state.v_u[i] > 0.0 && state.g_u[i].is_finite() {
                let slack = (state.g_u[i] - state.g[i]).max(1e-20);
                state.v_u[i] = mu_ks / slack;
            }
        }

        // NaN/Inf guard on evaluation
        if state.obj.is_nan() || state.obj.is_infinite() {
            // Try restoration from current point
            let (x_rest, success) = restoration.restore(
                &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
                &state.jac_rows, &state.jac_cols, n, m, options,
                &|theta, phi| filter.is_acceptable(theta, phi),
                &|x_eval, g_out| problem.constraints(x_eval, g_out),
                &|x_eval, jac_out| problem.jacobian_values(x_eval, jac_out),
                Some(&|x_eval: &[f64]| problem.objective(x_eval)),
            );
            if success {
                state.x = x_rest;
                state.alpha_primal = 0.0;
                state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
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

        // --- Barrier parameter update (free/fixed mode) ---
        // Special case: no bound constraints → no barrier, force mu to mu_min (Ipopt convention)
        let has_bounds = (0..n).any(|i| state.x_l[i].is_finite() || state.x_u[i].is_finite());
        if !has_bounds {
            state.mu = options.mu_min;
        } else {
            let kkt_error = {
                let pi = state.constraint_violation();
                let di = convergence::dual_infeasibility(
                    &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
                    &state.y, &state.z_l, &state.z_u, n,
                );
                let ci = convergence::complementarity_error(
                    &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
                );
                pi * pi + di * di + ci * ci
            };

            let sufficient = mu_state.check_sufficient_progress(kkt_error);

            match mu_state.mode {
                MuMode::Free => {
                    if sufficient && !mu_state.tiny_step {
                        mu_state.remember_accepted(kkt_error);
                        // Adaptive mu from complementarity (Loqo-style)
                        let avg_compl = compute_avg_complementarity(&state);
                        if avg_compl > 0.0 {
                            state.mu = (avg_compl / options.kappa).clamp(options.mu_min, 1e5);
                        } else {
                            // No complementarity products (all bounds inactive) → decrease mu
                            state.mu = (options.mu_linear_decrease_factor * state.mu)
                                .max(options.mu_min);
                        }
                        // In free mode: reset filter each iteration
                        filter.reset();
                        let theta_new = state.constraint_violation();
                        filter.set_theta_min_from_initial(theta_new);
                    } else {
                        // Switch to fixed mode
                        log::debug!("Switching to fixed mu mode (insufficient progress or tiny step)");
                        state.diagnostics.mu_mode_switches += 1;
                        mu_state.mode = MuMode::Fixed;
                        mu_state.first_iter_in_mode = true;
                        let avg_compl = compute_avg_complementarity(&state);
                        if avg_compl > 0.0 {
                            state.mu = (options.adaptive_mu_monotone_init_factor * avg_compl)
                                .clamp(options.mu_min, 1e5);
                        } else {
                            state.mu = (options.mu_linear_decrease_factor * state.mu)
                                .max(options.mu_min);
                        }
                        filter.reset();
                        let theta_new = state.constraint_violation();
                        filter.set_theta_min_from_initial(theta_new);
                    }
                }
                MuMode::Fixed => {
                    if sufficient && !mu_state.tiny_step && !mu_state.first_iter_in_mode {
                        // Switch back to free mode
                        log::debug!("Switching back to free mu mode (sufficient progress)");
                        state.diagnostics.mu_mode_switches += 1;
                        mu_state.mode = MuMode::Free;
                        mu_state.remember_accepted(kkt_error);
                        mu_state.first_iter_in_mode = true;
                    } else {
                        mu_state.first_iter_in_mode = false;
                        // Check if subproblem is solved (barrier error small enough)
                        let barrier_err = compute_barrier_error(&state);
                        if barrier_err <= options.barrier_tol_factor * state.mu || mu_state.tiny_step {
                            let new_mu = (options.mu_linear_decrease_factor * state.mu)
                                .min(state.mu.powf(options.mu_superlinear_decrease_power))
                                .max(options.mu_min);
                            if !(mu_state.tiny_step && (new_mu - state.mu).abs() < 1e-20) {
                                state.mu = new_mu;
                                filter.reset();
                                let theta_new = state.constraint_violation();
                                filter.set_theta_min_from_initial(theta_new);
                                log::debug!("Fixed mode: mu decreased to {:.2e}", state.mu);
                            }
                        }
                    }
                }
            }
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
        let fpr_tol = options.acceptable_tol.max(options.acceptable_constr_viol_tol);
        let sc = final_primal <= fpr_tol
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
        state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }

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

    // At max_iter: if the best du seen during the solve is below the acceptable floor,
    // restore the best-du point and declare Acceptable. This catches cycling problems
    // (e.g. HAHN1LS) where du oscillates but periodically reaches good values.
    let fdu_floor = (options.acceptable_tol * {
        let s_max = 100.0_f64;
        let fm: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        if (m + 2 * n) > 0 {
            (s_max.max(fm / (m + 2 * n) as f64) / s_max).min(1e4)
        } else { 1.0 }
    }).max(1e-2);
    if best_du_val < fdu_floor {
        if let Some(ref bdx) = best_du_x {
            state.x.copy_from_slice(bdx);
            if let Some(ref bdy) = best_du_y { state.y.copy_from_slice(bdy); }
            if let Some(ref bdzl) = best_du_zl { state.z_l.copy_from_slice(bdzl); }
            if let Some(ref bdzu) = best_du_zu { state.z_u.copy_from_slice(bdzu); }
            state.evaluate_with_linear(problem, 1.0, linear_constraints.as_deref(), lbfgs_mode);
        if let Some(ref mut lbfgs) = lbfgs_state {
            let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
                &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
            );
            lbfgs.update(&state.x, &lag_grad);
            lbfgs.fill_hessian(&mut state.hess_vals);
        }
            let bd_pr = state.constraint_violation();
            let bd_pr_tol = options.acceptable_tol.max(options.acceptable_constr_viol_tol);
            if bd_pr <= bd_pr_tol {
                if options.print_level >= 5 {
                    eprintln!(
                        "ripopt: best-du point passes acceptable (du={:.2e}, pr={:.2e})",
                        best_du_val, bd_pr
                    );
                }
                return make_result(&state, SolveStatus::Acceptable);
            }
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
    if options.print_level >= 5 {
        timings.print_summary(options.max_iter, ipm_start.elapsed());
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

    let kappa_soc = 0.99;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

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

        // Stop SOC iterations if theta isn't decreasing sufficiently
        if theta_soc >= kappa_soc * theta_prev_soc {
            return None;
        }
        theta_prev_soc = theta_soc;

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

/// Attempt a second-order correction step using the condensed KKT system.
///
/// Same logic as `attempt_soc` but uses the condensed system to avoid building
/// the full (n+m)×(n+m) KKT matrix. Uses `solve_condensed_soc` which rebuilds
/// only the n-dimensional condensed RHS with the modified constraint residual.
#[allow(clippy::too_many_arguments)]
fn attempt_soc_condensed<P: NlpProblem>(
    state: &SolverState,
    problem: &P,
    g_trial: &[f64],
    solver: &mut DenseLdl,
    condensed: &kkt::CondensedKktSystem,
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

    let kappa_soc = 0.99;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

    for _soc_iter in 0..options.max_soc {
        // Solve condensed system with modified constraint residual
        let dx_soc = match kkt::solve_condensed_soc(condensed, solver, &c_soc) {
            Ok(d) => d,
            Err(_) => return None,
        };

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

        // Stop SOC iterations if theta isn't decreasing sufficiently
        if theta_soc >= kappa_soc * theta_prev_soc {
            return None;
        }
        theta_prev_soc = theta_soc;

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

/// SOC using sparse condensed KKT system.
fn attempt_soc_sparse_condensed<P: NlpProblem>(
    state: &SolverState,
    problem: &P,
    g_trial: &[f64],
    solver: &mut dyn LinearSolver,
    condensed: &kkt::SparseCondensedKktSystem,
    filter: &Filter,
    theta_current: f64,
    phi_current: f64,
    grad_phi_step: f64,
    alpha: f64,
    options: &SolverOptions,
) -> Option<(Vec<f64>, f64, Vec<f64>, f64)> {
    let n = state.n;
    let m = state.m;
    if m == 0 { return None; }

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

    let kappa_soc = 0.99;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

    for _soc_iter in 0..options.max_soc {
        let dx_soc = match kkt::solve_sparse_condensed_soc(condensed, solver, &c_soc) {
            Ok(d) => d,
            Err(_) => return None,
        };

        let mut x_soc = vec![0.0; n];
        for i in 0..n {
            x_soc[i] = state.x[i] + alpha * dx_soc[i];
            if state.x_l[i].is_finite() { x_soc[i] = x_soc[i].max(state.x_l[i] + 1e-14); }
            if state.x_u[i].is_finite() { x_soc[i] = x_soc[i].min(state.x_u[i] - 1e-14); }
        }

        let obj_soc = problem.objective(&x_soc);
        let mut g_soc = vec![0.0; m];
        problem.constraints(&x_soc, &mut g_soc);

        let theta_soc = convergence::primal_infeasibility(&g_soc, &state.g_l, &state.g_u);
        if theta_soc >= kappa_soc * theta_prev_soc { return None; }
        theta_prev_soc = theta_soc;

        let mut phi_soc = obj_soc;
        for i in 0..n {
            if state.x_l[i].is_finite() {
                phi_soc -= state.mu * (x_soc[i] - state.x_l[i]).max(1e-20).ln();
            }
            if state.x_u[i].is_finite() {
                phi_soc -= state.mu * (state.x_u[i] - x_soc[i]).max(1e-20).ln();
            }
        }
        if options.constraint_slack_barrier {
            for i in 0..m {
                let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                    && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
                if is_eq { continue; }
                if state.g_l[i].is_finite() {
                    let slack = g_soc[i] - state.g_l[i];
                    if slack > state.mu * 1e-2 { phi_soc -= state.mu * slack.ln(); }
                }
                if state.g_u[i].is_finite() {
                    let slack = state.g_u[i] - g_soc[i];
                    if slack > state.mu * 1e-2 { phi_soc -= state.mu * slack.ln(); }
                }
            }
        }

        let (acceptable, _) = filter.check_acceptability(
            theta_current, phi_current, theta_soc, phi_soc, grad_phi_step, alpha,
        );
        if acceptable { return Some((x_soc, obj_soc, g_soc, alpha)); }

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

/// Apply post-restoration success handling: update state, reset multipliers, filter, and mu.
fn apply_restoration_success<P: NlpProblem>(
    state: &mut SolverState,
    filter: &mut Filter,
    mu_state: &mut MuState,
    options: &SolverOptions,
    n: usize,
    m: usize,
    problem: &P,
    x_new: &[f64],
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
    lbfgs_state: &mut Option<LbfgsIpmState>,
) {
    state.x.copy_from_slice(x_new);
    state.alpha_primal = 0.0;
    state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    if let Some(ref mut lbfgs) = lbfgs_state {
        let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
        );
        lbfgs.update(&state.x, &lag_grad);
        lbfgs.fill_hessian(&mut state.hess_vals);
    }

    // Reset multipliers after restoration (Ipopt-style).
    let bound_mult_reset_threshold = 1000.0;
    let mu_for_reset = state.mu;
    let mut any_large = false;
    for i in 0..n {
        if state.x_l[i].is_finite() {
            let slack = (state.x[i] - state.x_l[i]).max(1e-12);
            let z_new = mu_for_reset / slack;
            if z_new > bound_mult_reset_threshold {
                any_large = true;
            }
            state.z_l[i] = z_new.min(bound_mult_reset_threshold);
        }
        if state.x_u[i].is_finite() {
            let slack = (state.x_u[i] - state.x[i]).max(1e-12);
            let z_new = mu_for_reset / slack;
            if z_new > bound_mult_reset_threshold {
                any_large = true;
            }
            state.z_u[i] = z_new.min(bound_mult_reset_threshold);
        }
    }
    if any_large {
        for i in 0..n {
            if state.x_l[i].is_finite() {
                let slack = (state.x[i] - state.x_l[i]).max(1e-12);
                state.z_l[i] = (mu_for_reset / slack).min(bound_mult_reset_threshold);
            } else {
                state.z_l[i] = 0.0;
            }
            if state.x_u[i].is_finite() {
                let slack = (state.x_u[i] - state.x[i]).max(1e-12);
                state.z_u[i] = (mu_for_reset / slack).min(bound_mult_reset_threshold);
            } else {
                state.z_u[i] = 0.0;
            }
        }
    }
    // Reset constraint multipliers to zero.
    for i in 0..m {
        state.y[i] = 0.0;
    }

    // Reset constraint slack barrier multipliers v_l, v_u from mu/slack.
    let mu_r = state.mu;
    for i in 0..m {
        let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_eq {
            state.v_l[i] = 0.0;
            state.v_u[i] = 0.0;
            continue;
        }
        if state.g_l[i].is_finite() {
            let slack = (state.g[i] - state.g_l[i]).max(1e-12);
            state.v_l[i] = (mu_r / slack).min(bound_mult_reset_threshold);
        } else {
            state.v_l[i] = 0.0;
        }
        if state.g_u[i].is_finite() {
            let slack = (state.g_u[i] - state.g[i]).max(1e-12);
            state.v_u[i] = (mu_r / slack).min(bound_mult_reset_threshold);
        } else {
            state.v_u[i] = 0.0;
        }
    }

    // Reset filter and re-initialize from restored point
    filter.reset();
    let theta_restored = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_restored);
    state.consecutive_acceptable = 0;

    // Recompute mu from current complementarity after restoration
    let mu_compl = compute_avg_complementarity(state);
    if mu_compl > 0.0 {
        state.mu = mu_compl.max(options.mu_min).min(1e5);
    }

    // Reset mu_state to free mode (restoration is a fresh start)
    mu_state.mode = MuMode::Free;
    mu_state.first_iter_in_mode = true;
    mu_state.ref_vals.clear();
    mu_state.consecutive_restoration_failures = 0;
}

/// Outcome of the NLP restoration attempt.
enum RestorationOutcome {
    /// Restoration found a point with improved feasibility.
    Success,
    /// Local infeasibility: restoration converged but constraints not improved.
    LocalInfeasibility,
    /// Restoration failed (inner solve did not converge).
    Failed,
}

/// Attempt full NLP restoration by solving a restoration subproblem with the IPM.
///
/// Formulates a restoration NLP with p/n slack variables and solves it using
/// the same IPM engine (with `disable_nlp_restoration=true` to prevent recursion).
fn attempt_nlp_restoration<P: NlpProblem>(
    problem: &P,
    state: &SolverState,
    filter: &Filter,
    options: &SolverOptions,
    theta_current: f64,
) -> (Vec<f64>, RestorationOutcome) {
    let n = state.n;
    let m = state.m;

    if options.print_level >= 5 {
        eprintln!(
            "ripopt: Entering NLP restoration (theta={:.2e}, mu={:.2e})",
            theta_current, state.mu
        );
    }

    // Build restoration NLP
    let resto_nlp = RestorationNlp::new(problem, &state.x, state.mu, 1000.0, 1.0);

    // Configure inner solver options
    let mut inner_opts = options.clone();
    inner_opts.max_iter = options.restoration_max_iter;
    inner_opts.disable_nlp_restoration = true; // prevent recursion
    inner_opts.print_level = if options.print_level >= 5 { 3 } else { 0 };
    inner_opts.mu_init = state.mu.max(1e-2);
    // Relax convergence tolerances — we just need feasibility, not optimality
    inner_opts.tol = 1e-7;
    inner_opts.acceptable_tol = 1e-3;
    inner_opts.acceptable_iter = 5;

    // Solve the restoration NLP
    let result = solve_ipm(&resto_nlp, &inner_opts);

    // Extract x_orig from the restoration solution
    let x_nlp: Vec<f64> = result.x[..n].to_vec();

    // Evaluate original constraints at the restored point
    let mut g_new = vec![0.0; m];
    problem.constraints(&x_nlp, &mut g_new);
    let theta_new = convergence::primal_infeasibility(&g_new, &state.g_l, &state.g_u);

    // Evaluate original objective at the restored point
    let phi_new = problem.objective(&x_nlp);

    if options.print_level >= 5 {
        eprintln!(
            "ripopt: NLP restoration result: status={:?}, theta_new={:.2e} (was {:.2e}), phi_new={:.2e}",
            result.status, theta_new, theta_current, phi_new
        );
    }

    let inner_converged =
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable;

    // Check success criteria — require meaningful improvement
    if theta_new < options.constr_viol_tol {
        // Achieved feasibility
        return (x_nlp, RestorationOutcome::Success);
    }

    // Require >=50% reduction for non-feasible improvement (stricter than GN's 10%)
    // to avoid marginal "success" that prevents recovery mechanisms from engaging.
    if theta_new <= 0.5 * theta_current {
        return (x_nlp, RestorationOutcome::Success);
    }

    // Check if acceptable to outer filter AND has meaningful reduction
    if theta_new < 0.9 * theta_current {
        let filter_acceptable = {
            let entries = filter.entries();
            let theta_max = filter.theta_max();
            let gamma_theta = filter.gamma_theta();
            let gamma_phi = filter.gamma_phi();

            if theta_new.is_nan() || phi_new.is_nan() || theta_new > theta_max {
                false
            } else {
                let mut ok = true;
                for entry in entries {
                    if theta_new >= (1.0 - gamma_theta) * entry.theta
                        && phi_new >= entry.phi - gamma_phi * entry.theta
                    {
                        ok = false;
                        break;
                    }
                }
                ok
            }
        };

        if filter_acceptable {
            return (x_nlp, RestorationOutcome::Success);
        }
    }

    // Inner solve converged but didn't improve feasibility → locally infeasible
    if inner_converged {
        return (x_nlp, RestorationOutcome::LocalInfeasibility);
    }

    (x_nlp, RestorationOutcome::Failed)
}

/// Compute average complementarity for recomputing mu after restoration.
fn compute_avg_complementarity(state: &SolverState) -> f64 {
    let mut sum_compl = 0.0;
    let mut count = 0;
    // Variable bound complementarity: z_l * (x - x_l), z_u * (x_u - x)
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
    // If no variable bounds exist but inequality constraints do, include
    // constraint slack complementarity v_l*(g-g_l), v_u*(g_u-g) as fallback.
    // This prevents avg_compl=0 for problems with only inequality constraints
    // (e.g., OET2/6/7 with m=1002 inequalities and no variable bounds),
    // which otherwise causes mu to collapse to mu_min prematurely.
    //
    // When variable bounds exist, their z*slack products already drive mu;
    // adding v*slack (which ≈ mu since v = mu/slack) would bias avg_compl
    // and slow convergence (causes TP044/TP116 regressions).
    if count == 0 {
        for i in 0..state.m {
            if state.v_l[i] > 0.0 {
                let slack = (state.g[i] - state.g_l[i]).max(1e-20);
                sum_compl += state.v_l[i] * slack;
                count += 1;
            }
            if state.v_u[i] > 0.0 {
                let slack = (state.g_u[i] - state.g[i]).max(1e-20);
                sum_compl += state.v_u[i] * slack;
                count += 1;
            }
        }
    }
    if count > 0 {
        sum_compl / count as f64
    } else {
        0.0
    }
}

/// Compute barrier error for fixed-mode subproblem convergence check.
/// This is the optimality error of the current barrier subproblem (for fixed mu).
fn compute_barrier_error(state: &SolverState) -> f64 {
    let n = state.n;

    // Dual infeasibility of barrier problem:
    // grad_f + J^T y - z_l + z_u
    let mut grad_lag = state.grad_f.clone();
    for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        grad_lag[col] += state.jac_vals[idx] * state.y[row];
    }
    for i in 0..n {
        if state.x_l[i].is_finite() {
            grad_lag[i] -= state.z_l[i];
        }
        if state.x_u[i].is_finite() {
            grad_lag[i] += state.z_u[i];
        }
    }

    let sd = n.max(1) as f64;
    let dual_err = grad_lag.iter().map(|v| v.abs()).sum::<f64>() / sd;

    // Complementarity error (relative to mu)
    let mut compl_err = 0.0;
    let mut count = 0;
    for i in 0..n {
        if state.x_l[i].is_finite() {
            let slack = (state.x[i] - state.x_l[i]).max(1e-20);
            compl_err += (slack * state.z_l[i] - state.mu).abs();
            count += 1;
        }
        if state.x_u[i].is_finite() {
            let slack = (state.x_u[i] - state.x[i]).max(1e-20);
            compl_err += (slack * state.z_u[i] - state.mu).abs();
            count += 1;
        }
    }
    if count > 0 {
        compl_err /= count as f64;
    }

    // Primal infeasibility
    let primal_err = state.constraint_violation();

    dual_err.max(compl_err).max(primal_err)
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
    let n = state.n;
    let m = state.m;

    // Fill in final convergence measures for diagnostics
    let mut diag = state.diagnostics.clone();
    diag.final_mu = state.mu;
    diag.final_primal_inf = state.constraint_violation();
    diag.final_dual_inf = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    diag.final_compl = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );

    // Compute optimal z from stationarity in scaled space: ∇f_s + J_s^T y_s - z_l_s + z_u_s = 0
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
        diagnostics: diag,
    }
}
