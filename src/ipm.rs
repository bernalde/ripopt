use std::time::{Duration, Instant};

use crate::convergence::{self, check_convergence, ConvergenceInfo, ConvergenceStatus};
use crate::filter::{self, Filter, FilterEntry};
use crate::kkt::{self, InertiaCorrectionParams};
use crate::linear_solver::banded::BandedLdl;
use crate::linear_solver::dense::DenseLdl;
#[cfg(all(feature = "faer", not(feature = "rmumps")))]
use crate::linear_solver::sparse::SparseLdl;
#[cfg(feature = "rmumps")]
use crate::linear_solver::multifrontal::MultifrontalLdl;
#[cfg(feature = "rmumps")]
use crate::linear_solver::iterative::IterativeMinres;
#[cfg(feature = "rmumps")]
use crate::linear_solver::hybrid::HybridSolver;
use crate::linear_solver::{KktMatrix, LinearSolver, SymmetricMatrix};
use crate::options::LinearSolverChoice;

/// Create a new sparse linear solver using the best available backend.
/// Prefers rmumps (multifrontal) when available, falls back to faer (SparseLdl).
fn new_sparse_solver() -> Box<dyn LinearSolver> {
    new_sparse_solver_with_choice(LinearSolverChoice::Direct)
}

/// Create a sparse linear solver with the specified choice.
fn new_sparse_solver_with_choice(choice: LinearSolverChoice) -> Box<dyn LinearSolver> {
    match choice {
        LinearSolverChoice::Direct => {
            #[cfg(feature = "rmumps")]
            { return Box::new(MultifrontalLdl::new()); }
            #[cfg(all(not(feature = "rmumps"), feature = "faer"))]
            { return Box::new(SparseLdl::new()); }
            #[cfg(not(any(feature = "rmumps", feature = "faer")))]
            { return Box::new(DenseLdl::new()); }
        }
        LinearSolverChoice::Iterative => {
            #[cfg(feature = "rmumps")]
            { return Box::new(IterativeMinres::new()); }
            #[cfg(not(feature = "rmumps"))]
            {
                log::warn!("Iterative solver requires rmumps feature; falling back to direct");
                return new_sparse_solver_with_choice(LinearSolverChoice::Direct);
            }
        }
        LinearSolverChoice::Hybrid => {
            #[cfg(feature = "rmumps")]
            { return Box::new(HybridSolver::new()); }
            #[cfg(not(feature = "rmumps"))]
            {
                log::warn!("Hybrid solver requires rmumps feature; falling back to direct");
                return new_sparse_solver_with_choice(LinearSolverChoice::Direct);
            }
        }
    }
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
use crate::trace;
use crate::restoration_nlp::RestorationNlp;
use crate::result::{SolveResult, SolverDiagnostics, SolveStatus};
use crate::slack_formulation::SlackFormulation;
use crate::warmstart::WarmStartInitializer;
use crate::logging::rip_log;

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
    fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
        if !self.inner.objective(x, _new_x, obj) { return false; }
        *obj *= self.obj_scaling;
        true
    }
    fn gradient(&self, x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
        if !self.inner.gradient(x, _new_x, grad) { return false; }
        for g in grad.iter_mut() {
            *g *= self.obj_scaling;
        }
        true
    }
    fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
        if !self.inner.constraints(x, _new_x, g) { return false; }
        for (i, &s) in self.g_scaling.iter().enumerate() {
            g[i] *= s;
        }
        true
    }
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        self.inner.jacobian_structure()
    }
    fn jacobian_values(&self, x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
        if !self.inner.jacobian_values(x, _new_x, vals) { return false; }
        for (idx, &row) in self.jac_rows.iter().enumerate() {
            vals[idx] *= self.g_scaling[row];
        }
        true
    }
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        self.inner.hessian_structure()
    }
    fn hessian_values(&self, x: &[f64], _new_x: bool, obj_factor: f64, lambda: &[f64], vals: &mut [f64]) -> bool {
        let scaled_lambda: Vec<f64> = lambda
            .iter()
            .zip(self.g_scaling.iter())
            .map(|(l, s)| l * s)
            .collect();
        self.inner
            .hessian_values(x, _new_x, obj_factor * self.obj_scaling, &scaled_lambda, vals)
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
    /// Last point at which evaluations were performed (for new_x tracking).
    /// Initialized to NaN so the first evaluation always gets new_x = true.
    x_last_eval: Vec<f64>,
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
    /// Count of consecutive insufficient-progress iterations in Free mode.
    consecutive_insufficient: usize,
    /// Sliding window of dual infeasibility values for stagnation detection.
    dual_inf_window: Vec<f64>,
}

impl MuState {
    fn new() -> Self {
        Self {
            mode: MuMode::Free,
            ref_vals: Vec::with_capacity(8),
            num_refs_max: 4,
            refs_red_fact: 0.999,
            tiny_step: false,
            first_iter_in_mode: true,
            consecutive_restoration_failures: 0,
            consecutive_insufficient: 0,
            dual_inf_window: Vec::with_capacity(4),
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

        // Apply nlp_lower/upper_bound_inf: treat very large bounds as ±infinity.
        // This is the standard NLP convention (Ipopt uses the same option).
        // The GAMS/AMPL links map their "infinity" value to ±1e30 (which is finite
        // in f64). Without this conversion, z_l ≈ mu/1e30 and slack ≈ 1e30, giving
        // slack * z_l ≈ mu ≠ 0 as a spurious complementarity contribution that blocks
        // convergence detection even when the NLP is solved.
        for i in 0..n {
            if x_l[i] <= options.nlp_lower_bound_inf {
                x_l[i] = f64::NEG_INFINITY;
            }
            if x_u[i] >= options.nlp_upper_bound_inf {
                x_u[i] = f64::INFINITY;
            }
        }
        for i in 0..m {
            if g_l[i] <= options.nlp_lower_bound_inf {
                g_l[i] = f64::NEG_INFINITY;
            }
            if g_u[i] >= options.nlp_upper_bound_inf {
                g_u[i] = f64::INFINITY;
            }
        }

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
        // LS mult init: use dense path for small problems only.
        // Large problems skip init LS (sparse J*J^T is expensive and initial point
        // is often far from feasible where LS multipliers are meaningless).
        // The sparse LS path is available for post-restoration recovery.
        let mut y = if options.least_squares_mult_init && m > 0 && m <= 500 {
            let mut grad_f_init = vec![0.0; n];
            let grad_ok = problem.gradient(&x, true, &mut grad_f_init);

            let mut jac_vals_init = vec![0.0; jac_nnz];
            let jac_ok = problem.jacobian_values(&x, false, &mut jac_vals_init);
            if !grad_ok || !jac_ok {
                // Evaluation failed during LS mult init; skip and use default multipliers
                vec![0.0; m]
            } else {
                compute_ls_multiplier_estimate(
                    &grad_f_init, &jac_rows, &jac_cols, &jac_vals_init,
                    &g_l, &g_u, n, m, options.constr_mult_init_max,
                ).unwrap_or_else(|| vec![0.0; m])
            }
        } else {
            vec![0.0; m]
        };

        // Warm-start override: when warm_start is enabled and the problem
        // supplies initial multipliers, overwrite the defaults computed
        // above. WarmStartInitializer::initialize (called later in
        // `solve()`) will still floor z_l/z_u and recompute mu from
        // complementarity of the provided point.
        if options.warm_start {
            let mut ws_lam_g = vec![0.0; m];
            let mut ws_z_l = vec![0.0; n];
            let mut ws_z_u = vec![0.0; n];
            if problem.initial_multipliers(&mut ws_lam_g, &mut ws_z_l, &mut ws_z_u) {
                if m > 0 {
                    y.copy_from_slice(&ws_lam_g);
                }
                z_l.copy_from_slice(&ws_z_l);
                z_u.copy_from_slice(&ws_z_u);
            }
        }

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
            x_last_eval: vec![f64::NAN; n],
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
    ) -> bool {
        let new_x = self.x != self.x_last_eval;
        if !problem.objective(&self.x, new_x, &mut self.obj) { return false; }
        if !self.obj.is_finite() { return false; }
        if !problem.gradient(&self.x, false, &mut self.grad_f) { return false; }
        // NB: 42f4015 added element-wise is_finite checks on grad_f and g
        // here. Benchmarking showed they caused regressions on CUTEst
        // problems (BIGGS6NE, CERI651*, HS84/89/92, MAKELA3, OPTCNTRL, ...)
        // where a transient NaN/Inf in an individual grad or constraint
        // element at a boundary previously propagated but did not prevent
        // convergence. Matches Ipopt: IpOrigIpoptNLP only checks Nrm2 of
        // grad_f and its default `check_derivatives_for_naninf = false`
        // means it does NOT check Jacobian element-wise either. Keep only
        // the obj finite check (essential for convergence checks) and the
        // user callback bool return.
        if self.m > 0 {
            if !problem.constraints(&self.x, false, &mut self.g) { return false; }
            if !problem.jacobian_values(&self.x, false, &mut self.jac_vals) { return false; }
        }
        self.x_last_eval.copy_from_slice(&self.x);
        if skip_hessian {
            return true;
        }
        if let Some(flags) = linear_constraints {
            let mut lambda_for_hess = self.y.clone();
            for (i, &is_lin) in flags.iter().enumerate() {
                if is_lin {
                    lambda_for_hess[i] = 0.0;
                }
            }
            if !problem.hessian_values(&self.x, false, obj_factor, &lambda_for_hess, &mut self.hess_vals) { return false; }
        } else {
            if !problem.hessian_values(&self.x, false, obj_factor, &self.y, &mut self.hess_vals) { return false; }
        }
        true
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
    fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
        let m = self.targets.len();
        let mut g = vec![0.0; m];
        if !self.inner.constraints(x, _new_x, &mut g) { return false; }
        let mut sum = 0.0;
        for i in 0..m {
            let r = g[i] - self.targets[i];
            sum += r * r;
        }
        *obj = 0.5 * sum;
        true
    }
    fn gradient(&self, x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
        let n = self.inner.num_variables();
        let m = self.targets.len();
        let mut g = vec![0.0; m];
        if !self.inner.constraints(x, _new_x, &mut g) { return false; }
        let mut r = vec![0.0; m];
        for i in 0..m {
            r[i] = g[i] - self.targets[i];
        }
        let jac_nnz = self.jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];
        if !self.inner.jacobian_values(x, _new_x, &mut jac_vals) { return false; }
        for i in 0..n {
            grad[i] = 0.0;
        }
        for (idx, (&row, &col)) in self.jac_rows.iter().zip(self.jac_cols.iter()).enumerate() {
            grad[col] += jac_vals[idx] * r[row];
        }
        true
    }
    fn constraints(&self, _x: &[f64], _new_x: bool, _g: &mut [f64]) -> bool { true }
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }
    fn jacobian_values(&self, _x: &[f64], _new_x: bool, _vals: &mut [f64]) -> bool { true }
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hess_rows.clone(), self.hess_cols.clone())
    }
    fn hessian_values(&self, x: &[f64], _new_x: bool, obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) -> bool {
        let n = self.inner.num_variables();
        let m = self.targets.len();

        let mut g = vec![0.0; m];
        if !self.inner.constraints(x, _new_x, &mut g) { return false; }
        let mut r = vec![0.0; m];
        for i in 0..m {
            r[i] = g[i] - self.targets[i];
        }

        let jac_nnz = self.jac_rows.len();
        let mut jac_vals = vec![0.0; jac_nnz];
        if !self.inner.jacobian_values(x, _new_x, &mut jac_vals) { return false; }

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

        let inner_hess_nnz = self.inner_hess_map.len();
        if inner_hess_nnz > 0 {
            let mut inner_hess_vals = vec![0.0; inner_hess_nnz];
            if !self.inner.hessian_values(x, _new_x, 0.0, &r, &mut inner_hess_vals) { return false; }
            for (k, &v) in inner_hess_vals.iter().enumerate() {
                vals[self.inner_hess_map[k]] += obj_factor * v;
            }
        }
        true
    }
}

/// Detect if a problem is a nonlinear equation system (square or overdetermined).
///
/// Returns true if ALL of:
/// - f(x0) ≈ 0 (zero objective)
/// - ∇f(x0) ≈ 0 (zero gradient)
/// - All constraints are equalities (g_l[i] == g_u[i])
/// - m >= n (square or more constraints than variables)
/// Dense Cholesky solve: solve A*x = b where A is n×n symmetric positive definite.
/// Returns None if A is not positive definite (factorization fails).
fn dense_cholesky_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    // Cholesky: A = L * L^T
    let mut l = vec![0.0; n * n];
    for j in 0..n {
        let mut sum = 0.0;
        for k in 0..j {
            sum += l[j * n + k] * l[j * n + k];
        }
        let diag = a[j * n + j] - sum;
        if diag <= 0.0 {
            return None;
        }
        l[j * n + j] = diag.sqrt();

        for i in (j + 1)..n {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = (a[i * n + j] - sum) / l[j * n + j];
        }
    }

    // Forward solve: L * y = b
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut sum = 0.0;
        for k in 0..i {
            sum += l[i * n + k] * y[k];
        }
        y[i] = (b[i] - sum) / l[i * n + i];
    }

    // Back solve: L^T * x = y
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for k in (i + 1)..n {
            sum += l[k * n + i] * x[k];
        }
        x[i] = (y[i] - sum) / l[i * n + i];
    }

    Some(x)
}

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

    let mut f0 = 0.0;
    if !problem.objective(&x0, true, &mut f0) {
        return false;
    }
    if f0.abs() > 1e-10 {
        return false;
    }

    let mut grad = vec![0.0; n];
    if !problem.gradient(&x0, false, &mut grad) {
        return false;
    }
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
    if !problem.constraints(&x0, false, &mut g0) {
        return false;
    }
    let theta0 = convergence::primal_infeasibility(&g0, &g_l, &g_u);
    if theta0 < 1e-8 {
        return false;
    }

    true
}

/// Diagnosis of why the IPM failed, used to select targeted recovery strategies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FailureDiagnosis {
    /// Constraint violation is large and not decreasing.
    StallAtInfeasibility,
    /// All metrics are small but strict tolerances not met.
    StallNearOptimal,
    /// Factorization failed, NaN/Inf, or other numerical breakdown.
    NumericalBreakdown,
    /// Making progress but hit max iterations.
    SlowConvergence,
    /// Dual variables growing unboundedly (ill-conditioned Hessian).
    DualDivergence,
}

/// Classify the failure mode from diagnostics to select the best recovery strategy.
fn diagnose_failure(result: &SolveResult) -> FailureDiagnosis {
    let d = &result.diagnostics;
    match result.status {
        SolveStatus::NumericalError => {
            if d.final_dual_inf > 1e4 {
                FailureDiagnosis::DualDivergence
            } else {
                FailureDiagnosis::NumericalBreakdown
            }
        }
        SolveStatus::RestorationFailed => FailureDiagnosis::StallAtInfeasibility,
        SolveStatus::MaxIterations => {
            // Check dual divergence first — large dual infeasibility indicates
            // the Hessian or problem scaling is the root issue, not convergence speed.
            if d.final_dual_inf > 1e4 {
                FailureDiagnosis::DualDivergence
            } else if d.final_primal_inf > 1e-2 {
                FailureDiagnosis::StallAtInfeasibility
            } else if d.final_primal_inf < 1e-6 && d.final_compl < 1e-4 {
                FailureDiagnosis::StallNearOptimal
            } else {
                FailureDiagnosis::SlowConvergence
            }
        }
        _ => FailureDiagnosis::SlowConvergence,
    }
}

/// Check if `candidate` is strictly better than `current`.
///
/// A candidate is "better" when it is Optimal **and** either:
/// 1. `current` is clearly bad (infeasible or objective unusable), or
/// 2. `candidate.objective < current.objective` (strict improvement).
///
/// If `current` is NumericalError/MaxIterations but reached a feasible point
/// (`pr ≤ constr_viol_tol` and a finite objective), its objective is meaningful
/// and the candidate should only replace it if it finds a strictly lower objective.
/// Without this guard, a fallback that converges to a worse local minimum can
/// silently replace a near-optimal main-IPM iterate.
fn is_strictly_better(current: &SolveResult, candidate: &SolveResult) -> bool {
    let candidate_solved = matches!(candidate.status, SolveStatus::Optimal);
    if !candidate_solved {
        return false;
    }
    let current_solved = matches!(current.status, SolveStatus::Optimal);
    // Treat `current` as having a meaningful objective if it is feasible,
    // has a finite objective, and primal infeasibility is within the tolerance.
    let current_has_good_point = current.objective.is_finite()
        && current.diagnostics.final_primal_inf <= 1e-4;
    if current_solved || current_has_good_point {
        // Require strict objective improvement (with tiny tolerance for FP noise).
        let tol = 1e-8 * current.objective.abs().max(1.0);
        candidate.objective < current.objective - tol
    } else {
        // current is clearly bad — any Optimal candidate is an improvement.
        true
    }
}

/// Prepare options for a fallback solve: cap iterations and set remaining time budget.
/// Returns `None` if there is no time budget remaining.
fn prepare_fallback_opts(options: &SolverOptions, solve_start: &Instant) -> Option<SolverOptions> {
    let mut opts = options.clone();
    opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
    if options.max_wall_time > 0.0 {
        let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
        if remaining <= 0.1 {
            return None;
        }
        opts.max_wall_time = remaining;
    }
    Some(opts)
}

/// L-BFGS Hessian fallback: re-run IPM with `hessian_approximation_lbfgs`
/// enabled. Triggered for `NumericalBreakdown` and as a follow-up after a
/// failed plain-IPM retry. No-op when the user-provided Hessian is already
/// disabled or the fallback is opt-out.
fn try_lbfgs_hessian_fallback<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
) {
    if !options.enable_lbfgs_hessian_fallback || options.hessian_approximation_lbfgs {
        return;
    }
    let Some(mut opts) = prepare_fallback_opts(options, &solve_start) else { return };
    opts.hessian_approximation_lbfgs = true;
    opts.enable_lbfgs_hessian_fallback = false;
    opts.stall_iter_limit = 0;
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying L-BFGS Hessian fallback ({:?})", diagnosis);
    }
    let candidate = solve_ipm(problem, &opts);
    if is_strictly_better(result, &candidate) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: L-BFGS Hessian fallback succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        *result = candidate;
        result.diagnostics.fallback_used = Some("lbfgs_hessian".into());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: L-BFGS Hessian fallback did not improve ({:?})", candidate.status);
    }
}

/// Augmented Lagrangian fallback: solve via `crate::augmented_lagrangian`.
/// Only fires for constrained problems; primarily used for
/// `SlowConvergence` failures.
fn try_al_fallback<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
    has_constraints: bool,
) {
    if !options.enable_al_fallback || !has_constraints {
        return;
    }
    let Some(opts) = prepare_fallback_opts(options, &solve_start) else { return };
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying AL fallback ({:?})", diagnosis);
    }
    let candidate = crate::augmented_lagrangian::solve(problem, &opts);
    if is_strictly_better(result, &candidate) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: AL fallback succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        *result = candidate;
        result.diagnostics.fallback_used = Some("augmented_lagrangian".into());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: AL fallback did not improve ({:?})", candidate.status);
    }
}

/// SQP fallback: solve via `crate::sqp`. Only fires for constrained
/// problems. Used for `StallAtInfeasibility`, `SlowConvergence`, and
/// `StallNearOptimal` (where SQP refines a near-optimal IPM iterate).
fn try_sqp_fallback<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
    has_constraints: bool,
) {
    if !options.enable_sqp_fallback || !has_constraints {
        return;
    }
    let Some(opts) = prepare_fallback_opts(options, &solve_start) else { return };
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying SQP fallback ({:?})", diagnosis);
    }
    let candidate = crate::sqp::solve(problem, &opts);
    if is_strictly_better(result, &candidate) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: SQP fallback succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        *result = candidate;
        result.diagnostics.fallback_used = Some("sqp".into());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: SQP fallback did not improve ({:?})", candidate.status);
    }
}

/// Slack reformulation fallback: re-run IPM on `SlackFormulation::new`,
/// adding explicit slack variables for inequality constraints. Returns
/// `Some(SolveResult)` (with x truncated back to the original variable
/// space) when the slack solve strictly improves on the current result;
/// `None` otherwise.
fn try_slack_fallback<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
    has_inequalities: bool,
) -> Option<SolveResult> {
    if !options.enable_slack_fallback || !has_inequalities {
        return None;
    }
    let mut opts = prepare_fallback_opts(options, &solve_start)?;
    opts.enable_slack_fallback = false;
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying slack fallback ({:?})", diagnosis);
    }
    let slack_prob = SlackFormulation::new(problem, &result.x);
    let candidate = solve_ipm(&slack_prob, &opts);
    if is_strictly_better(result, &candidate) {
        let n = problem.num_variables();
        let m = problem.num_constraints();
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: Slack fallback succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        let x_out = candidate.x[..n].to_vec();
        let mut g_out = vec![0.0; m];
        let _ = problem.constraints(&x_out, true, &mut g_out);
        let mut diag = candidate.diagnostics;
        diag.fallback_used = Some("slack".into());
        diag.wall_time_secs = solve_start.elapsed().as_secs_f64();
        return Some(SolveResult {
            x: x_out,
            objective: candidate.objective,
            constraint_multipliers: candidate.constraint_multipliers,
            bound_multipliers_lower: candidate.bound_multipliers_lower[..n].to_vec(),
            bound_multipliers_upper: candidate.bound_multipliers_upper[..n].to_vec(),
            constraint_values: g_out,
            status: candidate.status,
            iterations: result.iterations + candidate.iterations,
            diagnostics: diag,
        });
    } else if options.print_level >= 5 {
        rip_log!("ripopt: Slack fallback did not improve ({:?})", candidate.status);
    }
    None
}

/// Plain-IPM retry for `DualDivergence`: re-run IPM with Gondzio MCC,
/// Mehrotra PC, and stall detection disabled. Large dual infeasibility
/// often means the advanced Newton corrections steered the solver into
/// a bad basin.
fn try_plain_ipm_retry<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
) {
    let Some(mut opts) = prepare_fallback_opts(options, &solve_start) else { return };
    opts.gondzio_mcc_max = 0;
    opts.mehrotra_pc = false;
    opts.stall_iter_limit = 0;
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying plain IPM retry (no corrections) for DualDivergence");
    }
    let candidate = solve_ipm(problem, &opts);
    if is_strictly_better(result, &candidate) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: Plain IPM retry succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        *result = candidate;
        result.diagnostics.fallback_used = Some("plain_ipm".into());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: Plain IPM retry did not improve ({:?})", candidate.status);
    }
}

/// Dispatch failure-recovery fallbacks based on the diagnosis. Returns
/// `Some(SolveResult)` if the slack fallback fires for
/// `StallAtInfeasibility` and produces a result we want to return early.
/// Otherwise updates `result` in place and returns `None`.
fn dispatch_failure_recovery<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
    has_constraints: bool,
    has_inequalities: bool,
) -> Option<SolveResult> {
    match diagnosis {
        FailureDiagnosis::StallAtInfeasibility => {
            if let Some(slack_result) = try_slack_fallback(
                result, problem, options, solve_start, diagnosis, has_inequalities,
            ) {
                return Some(slack_result);
            }
            try_sqp_fallback(result, problem, options, solve_start, diagnosis, has_constraints);
        }
        FailureDiagnosis::NumericalBreakdown => {
            try_lbfgs_hessian_fallback(result, problem, options, solve_start, diagnosis);
        }
        FailureDiagnosis::DualDivergence => {
            try_plain_ipm_retry(result, problem, options, solve_start);
            if !matches!(result.status, SolveStatus::Optimal) {
                try_lbfgs_hessian_fallback(result, problem, options, solve_start, diagnosis);
            }
        }
        FailureDiagnosis::SlowConvergence => {
            try_al_fallback(result, problem, options, solve_start, diagnosis, has_constraints);
            if !matches!(result.status, SolveStatus::Optimal) {
                try_sqp_fallback(result, problem, options, solve_start, diagnosis, has_constraints);
            }
        }
        FailureDiagnosis::StallNearOptimal => {
            try_sqp_fallback(result, problem, options, solve_start, diagnosis, has_constraints);
        }
    }
    None
}

/// Slow-optimal slack fallback: if the initial IPM was Optimal but
/// started from a feasible point and the objective worsened (or didn't
/// improve) while consuming > 5% of the wall-time budget, it likely
/// converged to a bad local minimum — try slack reformulation.
fn try_slow_optimal_slack_fallback<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
    diagnosis: FailureDiagnosis,
    has_inequalities: bool,
    initial_feasible: bool,
    initial_obj: f64,
) -> Option<SolveResult> {
    if !(matches!(result.status, SolveStatus::Optimal)
        && has_inequalities
        && options.enable_slack_fallback
        && options.max_wall_time > 0.0)
    {
        return None;
    }
    let time_used = solve_start.elapsed().as_secs_f64();
    let worsened_from_feasible = initial_feasible
        && initial_obj.is_finite()
        && result.objective > initial_obj - 1e-3 * initial_obj.abs().max(1.0);
    if !(time_used > 0.05 * options.max_wall_time && worsened_from_feasible) {
        return None;
    }
    if options.print_level >= 5 {
        rip_log!(
            "ripopt: Slow-optimal detected (obj={:.4e}, init_obj={:.4e}, time={:.1}s/{:.1}s), trying slack fallback",
            result.objective, initial_obj, time_used, options.max_wall_time
        );
    }
    try_slack_fallback(result, problem, options, solve_start, diagnosis, has_inequalities)
}

/// Run preprocessing (fixed-variable and redundant-constraint elimination)
/// and, if it reduces the problem, recursively `solve` the smaller problem
/// then unmap the solution. Returns `Some(result)` only when the
/// preprocessed solve reaches `Optimal`; otherwise returns `None` so the
/// caller falls through to the unpreprocessed path.
///
/// The preprocessed solve is capped at 50% of the remaining wall-time budget
/// so the unpreprocessed retry has time left when preprocessing leads the
/// solver into a bad basin (e.g. ganges.gms).
fn try_preprocessed_solve<P: NlpProblem>(
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
) -> Option<SolveResult> {
    if !options.enable_preprocessing {
        return None;
    }
    let prep = crate::preprocessing::PreprocessedProblem::new(problem as &dyn NlpProblem, options.bound_push);
    if !prep.did_reduce() {
        return None;
    }
    if options.print_level >= 5 {
        rip_log!(
            "ripopt: Preprocessing reduced problem: {} fixed vars, {} redundant constraints ({}x{} -> {}x{})",
            prep.num_fixed(), prep.num_redundant(),
            problem.num_variables(), problem.num_constraints(),
            prep.num_variables(), prep.num_constraints(),
        );
    }
    let mut prep_opts = options.clone();
    prep_opts.enable_preprocessing = false;
    if options.max_wall_time > 0.0 {
        let elapsed = solve_start.elapsed().as_secs_f64();
        let remaining = (options.max_wall_time - elapsed).max(1.0);
        prep_opts.max_wall_time = remaining * 0.5;
    }
    let reduced_result = solve(&prep, &prep_opts);
    let result = prep.unmap_solution(&reduced_result);
    if matches!(result.status, SolveStatus::Optimal) {
        return Some(result);
    }
    if options.print_level >= 5 {
        rip_log!(
            "ripopt: Preprocessed solve failed ({:?}), retrying without preprocessing",
            result.status
        );
    }
    None
}

/// Detect an overdetermined nonlinear equation problem (f ≡ 0, all equalities,
/// m ≥ n) and, if detected, solve it by reformulating as the unconstrained
/// least-squares problem `min 0.5·||g(x) − target||²` via
/// `LeastSquaresProblem`. Post-solve, polish with Gauss-Newton on the
/// original system, falling back to constrained IPM (and AL for square
/// systems) when the LS residual is above `options.tol`.
///
/// Returns `Some(SolveResult)` when the NE detection fires, `None`
/// otherwise (caller should fall through to the normal IPM path).
fn try_ne_to_ls_reformulation<P: NlpProblem>(
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
) -> Option<SolveResult> {
    if !detect_ne_problem(problem) {
        return None;
    }
    let n = problem.num_variables();
    let m = problem.num_constraints();
    if options.print_level >= 5 {
        rip_log!(
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
    if !problem.constraints(&ls_result.x, true, &mut g_final) {
        let mut diag = ls_result.diagnostics.clone();
        diag.fallback_used = Some("ne-to-ls".into());
        return Some(SolveResult {
            x: ls_result.x,
            objective: ls_result.objective,
            constraint_multipliers: vec![0.0; m],
            bound_multipliers_lower: ls_result.bound_multipliers_lower,
            bound_multipliers_upper: ls_result.bound_multipliers_upper,
            constraint_values: g_final,
            status: SolveStatus::EvaluationError,
            iterations: ls_result.iterations,
            diagnostics: diag,
        });
    }
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    let mut theta = convergence::primal_infeasibility(&g_final, &g_l, &g_u);

    // Newton polish: if theta is close but not quite at tol, try a few
    // Gauss-Newton steps on the original system g(x) = target to drive
    // constraint violation below tol.
    let mut polished_x = ls_result.x.clone();
    if theta > options.tol && theta < 1e-2 {
        let mut x_l_var = vec![0.0; n];
        let mut x_u_var = vec![0.0; n];
        problem.bounds(&mut x_l_var, &mut x_u_var);
        let (jac_rows, jac_cols) = problem.jacobian_structure();
        let nnz = jac_rows.len();

        let target: Vec<f64> = (0..m).map(|i| {
            if (g_u[i] - g_l[i]).abs() < 1e-15 { g_l[i] } else { 0.5 * (g_l[i] + g_u[i]) }
        }).collect();

        let max_newton_iters = 20;
        for newton_iter in 0..max_newton_iters {
            let r: Vec<f64> = (0..m).map(|i| g_final[i] - target[i]).collect();

            let mut jac_vals = vec![0.0; nnz];
            if !problem.jacobian_values(&polished_x, true, &mut jac_vals) {
                break;
            }

            let mut j_dense = vec![0.0; m * n];
            for k in 0..nnz {
                j_dense[jac_rows[k] * n + jac_cols[k]] += jac_vals[k];
            }

            // Solve for dx using normal equations: (J^T J) dx = -J^T r
            let mut jtj = vec![0.0; n * n];
            let mut jtr = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    let mut s = 0.0;
                    for k in 0..m {
                        s += j_dense[k * n + i] * j_dense[k * n + j];
                    }
                    jtj[i * n + j] = s;
                }
                let mut s = 0.0;
                for k in 0..m {
                    s += j_dense[k * n + i] * r[k];
                }
                jtr[i] = s;
            }

            for i in 0..n {
                jtj[i * n + i] += 1e-14;
            }

            let dx = match dense_cholesky_solve(&jtj, &jtr, n) {
                Some(dx) => dx,
                None => break,
            };

            // Fraction-to-boundary line search on variable bounds
            let mut alpha = 1.0;
            let tau = 0.995;
            for i in 0..n {
                if dx[i] < 0.0 && x_l_var[i].is_finite() {
                    let max_step = -tau * (polished_x[i] - x_l_var[i]) / dx[i];
                    if max_step < alpha { alpha = max_step; }
                }
                if dx[i] > 0.0 && x_u_var[i].is_finite() {
                    let max_step = tau * (x_u_var[i] - polished_x[i]) / dx[i];
                    if max_step < alpha { alpha = max_step; }
                }
            }
            alpha = alpha.max(0.0).min(1.0);

            let mut trial_x = vec![0.0; n];
            let mut trial_g = vec![0.0; m];
            let mut best_alpha = alpha;
            let mut best_theta = theta;
            for _ in 0..10 {
                for i in 0..n {
                    trial_x[i] = polished_x[i] - best_alpha * dx[i];
                }
                if !problem.constraints(&trial_x, true, &mut trial_g) {
                    best_alpha *= 0.5;
                    continue;
                }
                let trial_theta = convergence::primal_infeasibility(&trial_g, &g_l, &g_u);
                if trial_theta < theta {
                    best_theta = trial_theta;
                    break;
                }
                best_alpha *= 0.5;
            }

            if best_theta >= theta * 0.999 {
                break;
            }

            for i in 0..n {
                polished_x[i] -= best_alpha * dx[i];
            }
            if !problem.constraints(&polished_x, true, &mut g_final) {
                break;
            }
            theta = best_theta;

            if options.print_level >= 5 {
                rip_log!(
                    "ripopt: Newton polish iter {}: theta={:.2e}, alpha={:.4}",
                    newton_iter + 1, theta, best_alpha,
                );
            }

            if theta < options.tol {
                break;
            }
        }
    }

    let g_out = g_final;

    let status = if theta < options.tol {
        SolveStatus::Optimal
    } else {
        SolveStatus::LocalInfeasibility
    };

    if options.print_level >= 5 {
        rip_log!(
            "ripopt: NE-to-LS result: obj_LS={:.4e}, constraint_violation={:.4e}, status={:?}",
            ls_result.objective, theta, status
        );
    }

    // Fall back to constrained IPM when LS reports infeasibility.
    let ls_converged = matches!(ls_result.status, SolveStatus::Optimal);
    if status == SolveStatus::LocalInfeasibility && (m == n || !ls_converged) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: LS reformulation reports infeasibility (theta={:.4e}, ls_status={:?}), falling back to constrained IPM",
                theta, ls_result.status
            );
        }
        let mut fallback_opts = options.clone();
        if options.max_wall_time > 0.0 {
            let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
            if remaining <= 0.1 {
                return Some(SolveResult {
                    x: ls_result.x,
                    objective: 0.0,
                    constraint_multipliers: vec![0.0; m],
                    bound_multipliers_lower: ls_result.bound_multipliers_lower,
                    bound_multipliers_upper: ls_result.bound_multipliers_upper,
                    constraint_values: g_out,
                    status: SolveStatus::MaxIterations,
                    iterations: ls_result.iterations,
                    diagnostics: SolverDiagnostics::default(),
                });
            }
            fallback_opts.max_wall_time = remaining;
        }
        let ipm_result = solve_ipm(problem, &fallback_opts);
        if matches!(ipm_result.status, SolveStatus::Optimal) {
            return Some(ipm_result);
        }
        // IPM fallback failed — try AL for square NE systems
        if options.enable_al_fallback {
            if options.print_level >= 5 {
                rip_log!(
                    "ripopt: NE constrained IPM fallback failed ({:?}), trying AL",
                    ipm_result.status
                );
            }
            let mut al_opts = options.clone();
            al_opts.max_iter = options.max_iter.min(500).max(options.max_iter / 3);
            if options.max_wall_time > 0.0 {
                let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
                if remaining <= 0.1 {
                    return Some(ipm_result);
                }
                al_opts.max_wall_time = remaining;
            }
            let al_result = crate::augmented_lagrangian::solve(problem, &al_opts);
            if matches!(al_result.status, SolveStatus::Optimal) {
                if options.print_level >= 5 {
                    rip_log!(
                        "ripopt: NE AL fallback succeeded ({:?}, obj={:.6e})",
                        al_result.status, al_result.objective
                    );
                }
                return Some(al_result);
            }
        }
        return Some(ipm_result);
    }

    // L-BFGS retry on LS problem when IPM found a local min with nonzero residual
    let (final_x, final_status, final_g, final_iters, final_zl, final_zu) =
        if status == SolveStatus::LocalInfeasibility && options.enable_lbfgs_fallback {
            if options.print_level >= 5 {
                rip_log!(
                    "ripopt: NE-to-LS LocalInfeasibility (theta={:.4e}), trying L-BFGS on LS",
                    theta
                );
            }
            let lbfgs_ls = crate::lbfgs::solve(&ls_problem, options);
            let mut g_lb = vec![0.0; m];
            let theta_lb = if problem.constraints(&lbfgs_ls.x, true, &mut g_lb) {
                convergence::primal_infeasibility(&g_lb, &g_l, &g_u)
            } else {
                f64::INFINITY
            };

            if theta_lb < theta {
                let new_status = if theta_lb < options.tol {
                    SolveStatus::Optimal
                } else {
                    SolveStatus::LocalInfeasibility
                };
                if options.print_level >= 5 {
                    rip_log!(
                        "ripopt: L-BFGS improved NE-to-LS (theta: {:.4e} -> {:.4e}, status={:?})",
                        theta, theta_lb, new_status
                    );
                }
                (lbfgs_ls.x, new_status, g_lb, lbfgs_ls.iterations,
                 lbfgs_ls.bound_multipliers_lower, lbfgs_ls.bound_multipliers_upper)
            } else {
                if options.print_level >= 5 {
                    rip_log!(
                        "ripopt: L-BFGS did not improve NE-to-LS (theta_lb={:.4e} >= theta={:.4e})",
                        theta_lb, theta
                    );
                }
                (polished_x, status, g_out, ls_result.iterations,
                 ls_result.bound_multipliers_lower, ls_result.bound_multipliers_upper)
            }
        } else {
            (polished_x, status, g_out, ls_result.iterations,
             ls_result.bound_multipliers_lower, ls_result.bound_multipliers_upper)
        };

    Some(SolveResult {
        x: final_x,
        objective: 0.0,
        constraint_multipliers: vec![0.0; m],
        bound_multipliers_lower: final_zl,
        bound_multipliers_upper: final_zu,
        constraint_values: final_g,
        status: final_status,
        iterations: final_iters,
        diagnostics: SolverDiagnostics::default(),
    })
}

/// Solve the NLP using the interior point method.
pub fn solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    let solve_start = Instant::now();

    // Capture initial objective and feasibility for slow-optimal detection.
    // NOTE: disabled -- extra problem evaluations here change CUTEst FP state and cause regressions.
    let (initial_obj, initial_feasible) = (f64::INFINITY, false);

    if let Some(result) = try_preprocessed_solve(problem, options, solve_start) {
        return result;
    }

    if let Some(result) = try_ne_to_ls_reformulation(problem, options, solve_start) {
        return result;
    }

    let mut result = run_initial_solve(problem, options);

    let diagnosis = diagnose_failure(&result);
    let has_constraints = problem.num_constraints() > 0;
    let has_inequalities = has_inequality_constraints(problem);

    if options.print_level >= 5 && !matches!(result.status, SolveStatus::Optimal) {
        rip_log!("ripopt: Failure diagnosis: {:?}", diagnosis);
    }

    try_conservative_ipm_retry(&mut result, problem, options, solve_start);

    if !matches!(result.status, SolveStatus::Optimal) {
        if let Some(slack_result) = dispatch_failure_recovery(
            &mut result, problem, options, solve_start, diagnosis,
            has_constraints, has_inequalities,
        ) {
            return slack_result;
        }
    }

    if let Some(slack_result) = try_slow_optimal_slack_fallback(
        &mut result, problem, options, solve_start, diagnosis,
        has_inequalities, initial_feasible, initial_obj,
    ) {
        return slack_result;
    }

    apply_late_optimality_promotion(&mut result, options);

    result.diagnostics.wall_time_secs = solve_start.elapsed().as_secs_f64();
    result
}

/// Run the initial solve, picking a method based on problem structure:
///
///   * **Unconstrained (`m == 0`) and L-BFGS fallback enabled**: try
///     L-BFGS first; if L-BFGS converges to Optimal, return its result.
///     Otherwise run IPM and return the better of the two.
///   * **Constrained with a wall-time budget**: cap the initial IPM at
///     50% of `max_wall_time` so SQP/slack/AL fallbacks have time left.
///   * **Constrained without a wall-time budget**: full IPM with the
///     unmodified options.
fn run_initial_solve<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    if options.enable_lbfgs_fallback && problem.num_constraints() == 0 {
        let lbfgs_result = crate::lbfgs::solve(problem, options);
        if matches!(lbfgs_result.status, SolveStatus::Optimal) {
            if options.print_level >= 5 {
                rip_log!(
                    "ripopt: L-BFGS solved unconstrained problem ({:?}, obj={:.6e})",
                    lbfgs_result.status, lbfgs_result.objective
                );
            }
            return lbfgs_result;
        }
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: L-BFGS failed ({:?}, obj={:.6e}), trying IPM",
                lbfgs_result.status, lbfgs_result.objective
            );
        }
        let ipm_result = solve_ipm(problem, options);
        if matches!(ipm_result.status, SolveStatus::Optimal) {
            ipm_result
        } else if lbfgs_result.objective < ipm_result.objective {
            lbfgs_result
        } else {
            ipm_result
        }
    } else if options.max_wall_time > 0.0 && problem.num_constraints() > 0 {
        let mut main_opts = options.clone();
        main_opts.max_wall_time = options.max_wall_time * 0.5;
        solve_ipm(problem, &main_opts)
    } else {
        solve_ipm(problem, options)
    }
}

/// Conservative IPM retry: revert v0.4.0 algorithmic changes (Gondzio MCC,
/// Mehrotra PC, stall detection) to recover the pre-regression trajectory.
/// This is the most reliable recovery for problems sensitive to Newton
/// direction changes (TRO3X3, STRATEC, MGH10LS, ACOPR30). Only fires for
/// `n ≤ 200` and a non-Optimal current result; consumes 70% of the
/// remaining wall-time budget when one is set.
fn try_conservative_ipm_retry<P: NlpProblem>(
    result: &mut SolveResult,
    problem: &P,
    options: &SolverOptions,
    solve_start: Instant,
) {
    let n_problem = problem.num_variables();
    if n_problem > 200 || matches!(result.status, SolveStatus::Optimal) {
        return;
    }
    let Some(mut opts) = prepare_fallback_opts(options, &solve_start) else {
        return;
    };
    opts.gondzio_mcc_max = 0;
    opts.mehrotra_pc = false;
    opts.stall_iter_limit = 0;
    opts.proactive_infeasibility_detection = true;
    opts.max_iter = options.max_iter;
    if options.max_wall_time > 0.0 {
        let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
        opts.max_wall_time = remaining * 0.7;
    }
    if options.print_level >= 5 {
        rip_log!("ripopt: Trying conservative IPM retry (no Gondzio/Mehrotra, no stall detection)");
    }
    let candidate = solve_ipm(problem, &opts);
    if is_strictly_better(result, &candidate) {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: Conservative retry succeeded ({:?}, obj={:.6e})",
                candidate.status, candidate.objective
            );
        }
        *result = candidate;
        result.diagnostics.fallback_used = Some("conservative_ipm".into());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: Conservative retry did not improve ({:?})", candidate.status);
    }
}

/// Promote a stalled result (NumericalError/MaxIterations/Acceptable) to
/// Optimal or Acceptable when the final iterative-z residuals already meet
/// the KKT gates. Applied AFTER all fallbacks so conservative retries can
/// fix wrong local minima before we accept them.
///
/// Uses Ipopt-style residual checks:
///   * Optimal if pr ≤ `constr_viol_tol`, co ≤ `compl_inf_tol`, du strict.
///   * Acceptable (if not already) when du ≤ 1e-6, pr ≤ 1e-2, co ≤ 1e-2
///     (Ipopt's `acceptable_tol`).
fn apply_late_optimality_promotion(result: &mut SolveResult, options: &SolverOptions) {
    if !matches!(result.status, SolveStatus::NumericalError | SolveStatus::MaxIterations | SolveStatus::Acceptable) {
        return;
    }
    let d = &result.diagnostics;
    let pr_ok = d.final_primal_inf <= options.constr_viol_tol;
    let co_ok = d.final_compl <= options.compl_inf_tol;
    let du_strict_ok = d.final_dual_inf <= options.dual_inf_tol
        && d.final_dual_inf <= options.tol * 1000.0;
    if pr_ok && co_ok && du_strict_ok {
        if options.print_level >= 5 {
            rip_log!(
                "ripopt: Late-optimal (pr={:.2e}, du={:.2e}, co={:.2e}), returning Optimal",
                d.final_primal_inf, d.final_dual_inf, d.final_compl
            );
        }
        result.status = SolveStatus::Optimal;
    } else if !matches!(result.status, SolveStatus::Acceptable) {
        let du_acc_ok = d.final_dual_inf <= 1e-6
            && d.final_primal_inf <= 1e-2
            && d.final_compl <= 1e-2;
        if du_acc_ok {
            if options.print_level >= 5 {
                rip_log!(
                    "ripopt: Late-acceptable (pr={:.2e}, du={:.2e}, co={:.2e}), returning Acceptable",
                    d.final_primal_inf, d.final_dual_inf, d.final_compl
                );
            }
            result.status = SolveStatus::Acceptable;
        }
    }
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

        rip_log!("\nPhase breakdown ({} iterations):", iterations);
        for (name, dur) in &phases {
            let secs = dur.as_secs_f64();
            let pct = if total_secs > 0.0 { 100.0 * secs / total_secs } else { 0.0 };
            rip_log!("  {:<20} {:>8.3}s ({:>5.1}%)", name, secs, pct);
        }
        let other_secs = other.as_secs_f64();
        let other_pct = if total_secs > 0.0 { 100.0 * other_secs / total_secs } else { 0.0 };
        rip_log!("  {:<20} {:>8.3}s ({:>5.1}%)", "Other", other_secs, other_pct);
        rip_log!("  {:<20} {:>8.3}s", "Total", total_secs);
    }
}

/// Dump a KKT system to disk for external solver benchmarking.
///
/// Writes two files to `dir`:
/// - `<name>_<iter:04>.mtx` — Matrix Market symmetric format, lower triangle, 1-indexed
/// - `<name>_<iter:04>.json` — Metadata: problem_name, iteration, n, m, rhs, inertia, status
///
/// All IO errors are logged as warnings and never propagate to the caller.
fn dump_kkt_matrix(
    dir: &std::path::Path,
    name: &str,
    iteration: usize,
    kkt: &kkt::KktSystem,
    inertia: Option<(usize, usize, usize)>,
    delta_w: f64,
    delta_c: f64,
) {
    use std::io::Write;

    if let Err(e) = std::fs::create_dir_all(dir) {
        log::warn!("kkt_dump: cannot create directory {}: {}", dir.display(), e);
        return;
    }

    let stem = format!("{}_{:04}", name, iteration);
    let dim = kkt.dim;

    // --- Matrix Market (.mtx) ---
    let mtx_path = dir.join(format!("{}.mtx", stem));
    let write_mtx = || -> std::io::Result<()> {
        // Collect lower-triangle entries (1-indexed for Matrix Market).
        let entries: Vec<(usize, usize, f64)> = match &kkt.matrix {
            KktMatrix::Dense(d) => {
                let mut v = Vec::with_capacity(dim * (dim + 1) / 2);
                for j in 0..dim {
                    for i in j..dim {
                        let val = d.get(i, j);
                        if val != 0.0 {
                            v.push((i + 1, j + 1, val));
                        }
                    }
                }
                v
            }
            KktMatrix::Sparse(s) => {
                // Triplets are upper triangle (row <= col). Flip each entry to lower
                // triangle by swapping indices, then aggregate duplicates.
                let mut map: std::collections::HashMap<(usize, usize), f64> =
                    std::collections::HashMap::with_capacity(s.triplet_rows.len());
                for k in 0..s.triplet_rows.len() {
                    let r = s.triplet_rows[k]; // r <= c (upper tri)
                    let c = s.triplet_cols[k];
                    // Lower-triangle key: larger index first → (c, r) with c >= r.
                    *map.entry((c, r)).or_insert(0.0) += s.triplet_vals[k];
                }
                let mut v: Vec<(usize, usize, f64)> = map
                    .into_iter()
                    .filter(|(_, val)| *val != 0.0)
                    .map(|((i, j), val)| (i + 1, j + 1, val))
                    .collect();
                // Sort column-major for reader convenience.
                v.sort_unstable_by_key(|&(i, j, _)| (j, i));
                v
            }
        };

        let mut file = std::fs::File::create(&mtx_path)?;
        writeln!(file, "%%MatrixMarket matrix coordinate real symmetric")?;
        writeln!(file, "{} {} {}", dim, dim, entries.len())?;
        for (i, j, v) in &entries {
            writeln!(file, "{} {} {:.17e}", i, j, v)?;
        }
        Ok(())
    };
    if let Err(e) = write_mtx() {
        log::warn!("kkt_dump: failed to write {}.mtx: {}", stem, e);
        return;
    }

    // --- JSON sidecar (.json) ---
    let json_path = dir.join(format!("{}.json", stem));
    let (pos, neg, zer) = inertia.unwrap_or((0, 0, 0));
    let write_json = || -> std::io::Result<()> {
        let meta = serde_json::json!({
            "problem_name": name,
            "iteration": iteration,
            "n": kkt.n,
            "m": kkt.m,
            "rhs": kkt.rhs,
            "inertia": { "positive": pos, "negative": neg, "zero": zer },
            "delta_w": delta_w,
            "delta_c": delta_c,
            "status": "ongoing"
        });
        let mut file = std::fs::File::create(&json_path)?;
        write!(file, "{}", meta)?;
        Ok(())
    };
    if let Err(e) = write_json() {
        log::warn!("kkt_dump: failed to write {}.json: {}", stem, e);
    }
}

/// Update the L-BFGS Hessian approximation (no-op when `lbfgs_state` is `None`).
///
/// Recomputes the Lagrangian gradient at the current iterate and feeds
/// `(x, grad_L)` into the L-BFGS memory, then fills `state.hess_vals` with
/// the dense approximate Hessian `B_k`. When exact Hessian is used
/// (`lbfgs_state == None`) this is a no-op — the Hessian was populated
/// directly by `state.evaluate_with_linear`.
///
/// Ipopt parallel: `HessianUpdater` strategy object
/// (`IpHessianUpdater.{hpp,cpp}` with `LimMemQuasiNewtonUpdater` /
/// `ExactHessianUpdater` subclasses).
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2). Replaces ~16 copies of the same 6-line pattern.
fn update_lbfgs_hessian(
    lbfgs_state: &mut Option<LbfgsIpmState>,
    state: &mut SolverState,
) {
    if let Some(ref mut lbfgs) = lbfgs_state {
        let lag_grad = LbfgsIpmState::compute_lagrangian_gradient(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals, &state.y, state.n,
        );
        lbfgs.update(&state.x, &lag_grad);
        lbfgs.fill_hessian(&mut state.hess_vals);
    }
}

/// Bundled output of the per-iteration KKT assembly phase.
///
/// Holds the barrier diagonal `sigma` and whichever of the three KKT
/// representations the dispatch selected. Exactly one of
/// `condensed_system`, `sparse_condensed_system`, `kkt_system_opt`
/// is `Some` on exit.
struct AssembledKkt {
    sigma: Vec<f64>,
    use_sparse_condensed: bool,
    condensed_system: Option<kkt::CondensedKktSystem>,
    sparse_condensed_system: Option<kkt::SparseCondensedKktSystem>,
    kkt_system_opt: Option<kkt::KktSystem>,
}

/// Build the KKT system(s) for the current iterate.
///
/// Three-way dispatch:
/// - Dense condensed Schur complement when `m >= 2n` and `n` is small.
///   Cost `O(n^2 m + n^3)` beats `O((n+m)^3)` when `m >> n`.
/// - Sparse condensed Schur when the problem is large with constraints
///   and dense condensed is not already selected.
/// - Full augmented (n+m)x(n+m) KKT otherwise.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition. Does not touch `timings` — the caller records the
/// elapsed duration.
/// Outcome of the backtracking line-search trial loop.
enum LineSearchOutcome {
    StepAccepted,
    Rejected,
    Return(SolveResult),
}

/// Run the Ipopt-style backtracking line search on the current direction.
///
/// Up to 40 α-halving trials starting from `alpha_primal_max`. Each trial:
///  1. form `x_trial = x + α·dx` (clamped strictly inside bounds),
///  2. evaluate f, g at the trial; reject NaN/Inf,
///  3. watchdog shortcut: accept the full step unconditionally,
///  4. compute φ_trial with the (optional) constraint-slack barrier,
///  5. run the filter acceptability test (switching + Armijo + sufficient
///     infeasibility reduction),
///  6. on failure, attempt a second-order correction (SOC) on the first
///     trial only (Maratos-effect guard),
///  7. otherwise backtrack.
///
/// Returns `Return(NumericalError)` if the early-iteration stall timeout
/// trips, `StepAccepted` on any accepted trial / SOC / watchdog full
/// step, `Rejected` if α falls below `min_alpha` / loop budget exhausts.
#[allow(clippy::too_many_arguments)]
/// Compute the barrier-augmented objective phi(x, g) = obj - mu*Σ ln(slack)
/// summed over (a) all finite variable bounds via x and (b) when
/// `constraint_slack_barrier` is on, all finite inequality-constraint
/// bounds via g whose slack exceeds mu*1e-2 (the small-slack guard
/// matches the line-search loop and the SOC routines).
fn compute_barrier_phi(
    obj: f64,
    x: &[f64],
    g: &[f64],
    state: &SolverState,
    n: usize,
    m: usize,
    constraint_slack_barrier: bool,
) -> f64 {
    let mut phi = obj;
    for i in 0..n {
        if state.x_l[i].is_finite() {
            let slack = (x[i] - state.x_l[i]).max(1e-20);
            phi -= state.mu * slack.ln();
        }
        if state.x_u[i].is_finite() {
            let slack = (state.x_u[i] - x[i]).max(1e-20);
            phi -= state.mu * slack.ln();
        }
    }
    if constraint_slack_barrier {
        for i in 0..m {
            let is_eq = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_eq {
                continue;
            }
            if state.g_l[i].is_finite() {
                let slack = g[i] - state.g_l[i];
                if slack > state.mu * 1e-2 {
                    phi -= state.mu * slack.ln();
                }
            }
            if state.g_u[i].is_finite() {
                let slack = state.g_u[i] - g[i];
                if slack > state.mu * 1e-2 {
                    phi -= state.mu * slack.ln();
                }
            }
        }
    }
    phi
}

/// Dispatch a Second-Order Correction attempt to the appropriate solver
/// path: dense condensed (cond_solver_for_soc + condensed_system), sparse
/// condensed, or full augmented KKT. Returns the SOC trial tuple
/// `(x, obj, g, alpha)` on acceptance, `None` otherwise. Caller checks
/// `theta_trial > theta_current` and `*ls_steps == 0` before calling.
#[allow(clippy::too_many_arguments)]
fn dispatch_soc_attempt<P: NlpProblem>(
    state: &SolverState,
    problem: &P,
    x_trial: &[f64],
    g_trial: &[f64],
    condensed_system: &Option<kkt::CondensedKktSystem>,
    cond_solver_for_soc: &mut Option<DenseLdl>,
    sparse_condensed_system: &Option<kkt::SparseCondensedKktSystem>,
    kkt_system_opt: &Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    filter: &Filter,
    theta_current: f64,
    phi_current: f64,
    grad_phi_step: f64,
    alpha: f64,
    options: &SolverOptions,
) -> Option<(Vec<f64>, f64, Vec<f64>, f64)> {
    if let (Some(cond), Some(cs)) = (condensed_system.as_ref(), cond_solver_for_soc.as_mut()) {
        attempt_soc_condensed(
            state, problem, g_trial, cs, cond, filter,
            theta_current, phi_current, grad_phi_step, alpha, options,
        )
    } else if let Some(sc) = sparse_condensed_system.as_ref() {
        attempt_soc_sparse_condensed(
            state, problem, g_trial, lin_solver, sc, filter,
            theta_current, phi_current, grad_phi_step, alpha, options,
        )
    } else if let Some(kkt) = kkt_system_opt.as_ref() {
        attempt_soc(
            state, problem, x_trial, g_trial,
            lin_solver, kkt, filter,
            theta_current, phi_current, grad_phi_step, alpha, options,
        )
    } else {
        None
    }
}

fn run_line_search_loop<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    filter: &mut Filter,
    condensed_system: &Option<kkt::CondensedKktSystem>,
    cond_solver_for_soc: &mut Option<DenseLdl>,
    sparse_condensed_system: &Option<kkt::SparseCondensedKktSystem>,
    kkt_system_opt: &Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    alpha_primal_max: f64,
    theta_current: f64,
    phi_current: f64,
    grad_phi_step: f64,
    min_alpha: f64,
    force_restoration: bool,
    watchdog_active: bool,
    iteration: usize,
    n: usize,
    m: usize,
    start_time: Instant,
    early_timeout: f64,
    last_soc_accepted: &mut bool,
    ls_steps: &mut usize,
) -> LineSearchOutcome {
    let mut alpha = alpha_primal_max;
    let mut step_accepted = false;
    *ls_steps = 0;

    for _ls_iter in 0..40 {
        if force_restoration {
            break;
        }
        // Intra-iteration early stall check (scaled by problem size)
        if iteration < 3 && options.early_stall_timeout > 0.0
            && start_time.elapsed().as_secs_f64() > early_timeout
        {
            return LineSearchOutcome::Return(make_result(state, SolveStatus::NumericalError));
        }
        if alpha < min_alpha {
            break;
        }

        // Compute trial point
        let mut x_trial = vec![0.0; n];
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            x_trial[i] = state.x[i] + alpha * state.dx[i];
            if state.x_l[i].is_finite() {
                x_trial[i] = x_trial[i].max(state.x_l[i] + 1e-14);
            }
            if state.x_u[i].is_finite() {
                x_trial[i] = x_trial[i].min(state.x_u[i] - 1e-14);
            }
        }

        // Evaluate at trial point
        let mut obj_trial = f64::INFINITY;
        let obj_ok = problem.objective(&x_trial, true, &mut obj_trial);
        let mut g_trial = vec![0.0; m];
        let constr_ok = if m > 0 {
            problem.constraints(&x_trial, true, &mut g_trial)
        } else {
            true
        };

        if !obj_ok || !constr_ok || obj_trial.is_nan() || obj_trial.is_infinite()
            || g_trial.iter().any(|v| v.is_nan() || v.is_infinite())
        {
            alpha *= 0.5;
            *ls_steps += 1;
            continue;
        }

        let theta_trial =
            convergence::primal_infeasibility(&g_trial, &state.g_l, &state.g_u);

        // Watchdog: accept full step unconditionally (bypass filter)
        if watchdog_active && alpha == alpha_primal_max {
            state.x = x_trial;
            state.obj = obj_trial;
            state.g = g_trial;
            state.alpha_primal = alpha;
            step_accepted = true;
            break;
        }

        // Barrier objective at trial
        let phi_trial = compute_barrier_phi(
            obj_trial, &x_trial, &g_trial, state, n, m, options.constraint_slack_barrier,
        );

        let (acceptable, used_switching) = filter.check_acceptability(
            theta_current,
            phi_current,
            theta_trial,
            phi_trial,
            grad_phi_step,
            alpha,
        );

        if !acceptable && options.print_level >= 7 && *ls_steps < 5 {
            let sw = filter.switching_condition(theta_current, grad_phi_step, alpha);
            let ar = filter.armijo_condition(phi_current, phi_trial, grad_phi_step, alpha);
            let sr = filter.sufficient_infeasibility_reduction(theta_current, theta_trial);
            let fa = filter.is_acceptable(theta_trial, phi_trial);
            rip_log!("  LS reject: alpha={:.2e} theta_t={:.2e} phi_t={:.2e} switch={} armijo={} suff_red={} filter_ok={}",
                alpha, theta_trial, phi_trial, sw, ar, sr, fa);
        }

        if acceptable {
            state.x = x_trial;
            state.obj = obj_trial;
            state.g = g_trial;
            state.alpha_primal = alpha;
            step_accepted = true;
            if !used_switching {
                filter.add(theta_current, phi_current);
            }
            break;
        }

        // SOC on the first trial only, if full step increased theta
        if theta_trial > theta_current && options.max_soc > 0 && *ls_steps == 0 {
            let soc_accepted = dispatch_soc_attempt(
                state, problem, &x_trial, &g_trial, condensed_system,
                cond_solver_for_soc, sparse_condensed_system, kkt_system_opt,
                lin_solver, filter, theta_current, phi_current, grad_phi_step,
                alpha, options,
            );
            if let Some((x_soc, obj_soc, g_soc, alpha_soc)) = soc_accepted {
                state.diagnostics.soc_corrections += 1;
                *last_soc_accepted = true;
                state.x = x_soc;
                state.obj = obj_soc;
                state.g = g_soc;
                state.alpha_primal = alpha_soc;
                step_accepted = true;
                filter.add(theta_current, phi_current);
                break;
            }
        }

        alpha *= 0.5;
        *ls_steps += 1;
    }

    if step_accepted {
        LineSearchOutcome::StepAccepted
    } else {
        LineSearchOutcome::Rejected
    }
}

fn assemble_kkt_systems(
    state: &SolverState,
    n: usize,
    m: usize,
    use_sparse: bool,
    disable_sparse_condensed: bool,
) -> AssembledKkt {
    let sigma = kkt::compute_sigma(&state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u);

    let use_condensed = m >= 2 * n && n > 0 && (!use_sparse || n <= 100);
    let use_sparse_condensed = use_sparse && m > 0 && !use_condensed && !disable_sparse_condensed;

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

    let kkt_system_opt = if !use_condensed && !use_sparse_condensed {
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

    AssembledKkt {
        sigma,
        use_sparse_condensed,
        condensed_system,
        sparse_condensed_system,
        kkt_system_opt,
    }
}

/// Apply the accepted step to the dual variables (y, z_l, z_u) and
/// enforce Ipopt's kappa_sigma safeguard on bound multipliers.
///
/// - `y` (constraint multipliers) use `alpha_for_y = primal`: the
///   accepted primal step length (Ipopt default). Near convergence
///   (`consecutive_acceptable >= 1`), persistent sign changes in
///   `dy_i` across ≥3 iterations trigger 50% damping to suppress
///   oscillation.
/// - `z_l`, `z_u` use `alpha_d` (fraction-to-boundary on bound
///   multipliers). After the update, each is clamped to keep `z*s`
///   inside `[mu_ks/kappa_sigma, kappa_sigma*mu_ks]` (Ipopt
///   `IpIpoptCalculatedQuantities::ComputePDSystem`-style safeguard
///   with `kappa_sigma = 1e10`).
///
/// Writes `state.alpha_dual = alpha_d`. Returns `mu_ks`, which the
/// caller reuses when resetting slack-constraint multipliers `v_l`,
/// `v_u` after the post-step re-evaluation.
fn update_dual_variables(
    state: &mut SolverState,
    mu_state: &MuState,
    alpha_dual_max: f64,
    prev_dy: &mut Option<Vec<f64>>,
    dy_sign_change_count: &mut [u8],
) -> f64 {
    let n = state.n;
    let m = state.m;
    let alpha_y = state.alpha_primal;
    let alpha_d = alpha_dual_max;
    let near_convergence = state.consecutive_acceptable >= 1;
    for i in 0..m {
        let sign_change = if let Some(ref pdy) = prev_dy {
            pdy[i] * state.dy[i] < 0.0
        } else {
            false
        };
        if near_convergence && sign_change {
            dy_sign_change_count[i] = dy_sign_change_count[i].saturating_add(1);
        } else if !sign_change {
            dy_sign_change_count[i] = 0;
        }
        let dy_i = if near_convergence && dy_sign_change_count[i] >= 3 {
            0.5 * state.dy[i]
        } else {
            state.dy[i]
        };
        state.y[i] += alpha_y * dy_i;
    }
    *prev_dy = Some(state.dy.clone());

    let kappa_sigma = 1e10;
    let mu_ks = if mu_state.mode == MuMode::Free {
        compute_avg_complementarity(state)
            .max(state.mu)
            .min(1e3)
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
    mu_ks
}

/// Apply Gondzio multiple centrality corrections (MCC).
///
/// After the main Newton direction has been computed, perform up to
/// `options.gondzio_mcc_max` additional centrality corrections using
/// the SAME factored KKT matrix (one extra backsolve per correction)
/// to drive complementarity pairs far from μ back toward the central
/// path.
///
/// Guards:
/// - Reject corrections whose solution magnitude exceeds `1e10`× the
///   RHS magnitude (null-space blow-ups on rank-deficient systems).
/// - Dampen corrections that deflect the direction by more than 30°
///   (basin-switching on nonconvex problems).
/// - Accept each correction only if it does not shrink α by more
///   than 10%.
///
/// No-op when `options.gondzio_mcc_max == 0` or when the KKT is
/// condensed (this helper only runs on the full-augmented path).
///
/// Reference: Gondzio (1994, Comput. Optim. Appl.); Gondzio (2007).
/// Combine the Gondzio MCC corrector step `(ddx, ddy)` with the
/// current Newton direction stored in `state.{dx, dy, dz_l, dz_u}`,
/// recovering the bound-multiplier corrections via
/// `dz_L = -(z_L/s_L)·ddx`, `dz_U = (z_U/s_U)·ddx` (no centering
/// term). If the resulting `dx_c` deflects more than ~45° from the
/// original `state.dx` (cos(angle) < 0.7), damp the correction with
/// blending factor `alpha_damp = 0.3` to keep the corrected step
/// close to the Newton direction.
///
/// Returns `(dx_c, dy_c, dz_l_c, dz_u_c)`.
fn build_mcc_corrected_direction(
    state: &SolverState,
    ddx: &[f64],
    ddy: &[f64],
    iteration: usize,
    n: usize,
    m: usize,
    dx_norm_orig: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
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

    let mut dx_c: Vec<f64> = state.dx.iter().zip(ddx.iter()).map(|(a, b)| a + b).collect();
    let mut dy_c: Vec<f64> = state.dy.iter().zip(ddy.iter()).map(|(a, b)| a + b).collect();
    let mut dz_l_c: Vec<f64> = state.dz_l.iter().zip(ddz_l.iter()).map(|(a, b)| a + b).collect();
    let mut dz_u_c: Vec<f64> = state.dz_u.iter().zip(ddz_u.iter()).map(|(a, b)| a + b).collect();

    if dx_norm_orig > 1e-30 {
        let dx_c_norm: f64 = dx_c.iter().map(|v| v * v).sum::<f64>().sqrt();
        if dx_c_norm > 1e-30 {
            let dot: f64 = state.dx.iter().zip(dx_c.iter()).map(|(a, b)| a * b).sum::<f64>();
            let cos_angle = dot / (dx_norm_orig * dx_c_norm);
            if cos_angle < 0.7 {
                let alpha_damp = 0.3;
                for i in 0..n { dx_c[i] = (1.0 - alpha_damp) * state.dx[i] + alpha_damp * dx_c[i]; }
                for i in 0..m { dy_c[i] = (1.0 - alpha_damp) * state.dy[i] + alpha_damp * dy_c[i]; }
                for i in 0..n { dz_l_c[i] = (1.0 - alpha_damp) * state.dz_l[i] + alpha_damp * dz_l_c[i]; }
                for i in 0..n { dz_u_c[i] = (1.0 - alpha_damp) * state.dz_u[i] + alpha_damp * dz_u_c[i]; }
                log::debug!(
                    "Gondzio MCC iter {}: dampened correction (cos={:.3})",
                    iteration, cos_angle
                );
            }
        }
    }

    (dx_c, dy_c, dz_l_c, dz_u_c)
}

/// Build the Gondzio multiple-centrality-corrections (MCC)
/// corrector RHS for the augmented KKT system. For each bound,
/// projects the trial complementarity `c = z·s` at the current
/// step `alpha_mcc` onto the centrality interval
/// `[beta_min·μ_target, beta_max·μ_target]`; only complementarities
/// outside this band contribute, with sign chosen so adding the
/// correction step pulls `c` back inside.
///
/// Returns the RHS vector (zero outside the bound rows) and a
/// boolean indicating whether at least one bound contributed —
/// when no bound contributes, the caller short-circuits the loop.
fn build_mcc_corrector_rhs(
    state: &SolverState,
    kkt_dim: usize,
    alpha_mcc: f64,
    mu_target: f64,
    beta_min: f64,
    beta_max: f64,
    n: usize,
) -> (Vec<f64>, bool) {
    let mut rhs_mcc = vec![0.0_f64; kkt_dim];
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
    (rhs_mcc, needs_correction)
}

/// Maximum step size α such that the iterate stays within
/// `tau_mcc` × distance-to-boundary of the bounds for both the
/// primal variables (lower and upper) and the bound multipliers.
/// Used twice inside the Gondzio MCC loop: once on the original
/// Newton direction and once on each candidate corrected direction.
fn compute_mcc_alpha_max(
    state: &SolverState,
    dx: &[f64],
    dz_l: &[f64],
    dz_u: &[f64],
    tau_mcc: f64,
    n: usize,
) -> f64 {
    let mcc_zl = filter::fraction_to_boundary(&state.z_l, dz_l, tau_mcc);
    let mcc_zu = filter::fraction_to_boundary(&state.z_u, dz_u, tau_mcc);
    let mut alpha = mcc_zl.min(mcc_zu).min(1.0);
    for i in 0..n {
        if state.x_l[i].is_finite() && dx[i] < 0.0 {
            let s = state.x[i] - state.x_l[i];
            alpha = alpha.min(tau_mcc * s / (-dx[i]));
        }
        if state.x_u[i].is_finite() && dx[i] > 0.0 {
            let s = state.x_u[i] - state.x[i];
            alpha = alpha.min(tau_mcc * s / dx[i]);
        }
    }
    alpha.clamp(0.0, 1.0)
}

fn apply_gondzio_mcc(
    state: &mut SolverState,
    options: &SolverOptions,
    iteration: usize,
    mu_state: &MuState,
    primal_inf: f64,
    dual_inf: f64,
    compl_inf: f64,
    kkt_system_opt: &Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
) {
    if options.gondzio_mcc_max == 0 {
        return;
    }
    let Some(kkt) = kkt_system_opt.as_ref() else { return };
    let n = state.n;
    let m = state.m;

    let tau_mcc = if mu_state.mode == MuMode::Free {
        let nlp_error = primal_inf + dual_inf + compl_inf;
        (1.0 - nlp_error).max(options.tau_min)
    } else {
        (1.0 - state.mu).max(options.tau_min)
    };
    let mut alpha_mcc = compute_mcc_alpha_max(
        state, &state.dx, &state.dz_l, &state.dz_u, tau_mcc, n,
    );

    let mu_target = state.mu;
    let beta_min = 0.01_f64;
    let beta_max = 100.0_f64;

    let dx_norm_orig: f64 = state.dx.iter().map(|v| v * v).sum::<f64>().sqrt();

    for _mcc_iter in 0..options.gondzio_mcc_max {
        let (rhs_mcc, needs_correction) = build_mcc_corrector_rhs(
            state, kkt.dim, alpha_mcc, mu_target, beta_min, beta_max, n,
        );
        if !needs_correction {
            break;
        }

        match kkt::solve_with_custom_rhs_refined(&kkt.matrix, kkt.n, kkt.dim, lin_solver, &rhs_mcc) {
            Ok((ddx, ddy)) => {
                let nrm_rhs: f64 = rhs_mcc.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
                let nrm_sol: f64 = ddx.iter().chain(ddy.iter()).map(|v| v.abs()).fold(0.0_f64, f64::max);
                if nrm_sol > 1e10 * nrm_rhs.max(1.0) {
                    log::debug!(
                        "Gondzio MCC iter {}: ||sol||={:.2e}, ||rhs||={:.2e} — rejecting",
                        iteration, nrm_sol, nrm_rhs,
                    );
                    break;
                }
                let (dx_c, dy_c, dz_l_c, dz_u_c) = build_mcc_corrected_direction(
                    state, &ddx, &ddy, iteration, n, m, dx_norm_orig,
                );

                let alpha_new = compute_mcc_alpha_max(
                    state, &dx_c, &dz_l_c, &dz_u_c, tau_mcc, n,
                );

                if alpha_new >= 0.9 * alpha_mcc {
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

/// Watchdog control-flow decision returned after update.
enum WatchdogDecision {
    Proceed,
    Continue,
}

/// Maintain the watchdog state after an accepted step.
///
/// - Track consecutive shortened steps
///   (`alpha_primal < 0.99 * alpha_primal_max`).
/// - Activate the watchdog and snapshot state when the shortened-step
///   counter hits `watchdog_shortened_iter_trigger`.
/// - If the watchdog is active, check for sufficient progress versus
///   the saved point; on success deactivate, on exhaustion
///   (`watchdog_trial_iter_max`) revert to the saved state and
///   request `Continue`.
///
/// Reference: Chamberlain et al. (1982); Ipopt's
/// `IpBacktrackingLineSearch::DoBacktrackingLineSearch` watchdog.
fn update_watchdog<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    iteration: usize,
    alpha_primal_max: f64,
    filter: &mut Filter,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    consecutive_shortened: &mut usize,
    watchdog_active: &mut bool,
    watchdog_trial_count: &mut usize,
    watchdog_saved: &mut Option<WatchdogSavedState>,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> WatchdogDecision {
    if state.alpha_primal < alpha_primal_max * 0.99 {
        *consecutive_shortened += 1;
    } else {
        *consecutive_shortened = 0;
    }

    if !*watchdog_active
        && *consecutive_shortened >= options.watchdog_shortened_iter_trigger
    {
        state.diagnostics.watchdog_activations += 1;
        *watchdog_active = true;
        *watchdog_trial_count = 0;
        let wd_theta = state.constraint_violation();
        let wd_phi = state.barrier_objective(options);
        *watchdog_saved = Some(WatchdogSavedState {
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
        *consecutive_shortened = 0;
        log::debug!(
            "Watchdog activated at iteration {} (theta={:.2e}, phi={:.2e})",
            iteration, wd_theta, wd_phi
        );
    }

    if *watchdog_active {
        *watchdog_trial_count += 1;
        if let Some(ref saved) = watchdog_saved {
            let theta_now = state.constraint_violation();
            let phi_now = state.barrier_objective(options);
            let made_progress = filter.is_acceptable(theta_now, phi_now)
                && (theta_now < (1.0 - 1e-5) * saved.theta
                    || phi_now < saved.phi - 1e-5 * saved.theta);

            if made_progress {
                log::debug!(
                    "Watchdog succeeded at trial {} (theta: {:.2e} -> {:.2e})",
                    *watchdog_trial_count, saved.theta, theta_now
                );
                *watchdog_active = false;
                *watchdog_trial_count = 0;
                *watchdog_saved = None;
            } else if *watchdog_trial_count >= options.watchdog_trial_iter_max {
                log::debug!(
                    "Watchdog reverting after {} trials",
                    *watchdog_trial_count
                );
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
                let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
                update_lbfgs_hessian(lbfgs_state, state);

                *watchdog_active = false;
                *watchdog_trial_count = 0;
                *watchdog_saved = None;
                return WatchdogDecision::Continue;
            }
        }
    }

    WatchdogDecision::Proceed
}

/// Outcome of the search-direction solve.
enum DirectionSolveDecision {
    Proceed {
        dx: Vec<f64>,
        dy: Vec<f64>,
        cond_solver_for_soc: Option<DenseLdl>,
        mehrotra_aff: Option<(Vec<f64>, Vec<f64>, Vec<f64>, f64)>,
    },
    Continue,
    Return(SolveResult),
}

/// Outcome of the dense-condensed branch in `solve_for_search_direction`.
/// `Solved` carries the primal/dual step and (when the BK factorization
/// of the n×n Schur complement succeeded) the cached `DenseLdl` used by
/// SOC. `Continue` and `Return` propagate restoration / NumericalError.
enum CondensedDirectionOutcome {
    Solved {
        dx: Vec<f64>,
        dy: Vec<f64>,
        cond_solver: Option<DenseLdl>,
    },
    Continue,
    Return(SolveResult),
}

/// Outcome of `restore_after_solve_failure`. Continue = restoration
/// succeeded and the iterate was updated; Return = restoration failed
/// and the caller should bubble up `NumericalError`.
enum SolveRestoreOutcome {
    Continue,
    Return(SolveResult),
}

/// Shared helper for the "factor or solve failed → call restoration"
/// pattern that appears in `solve_dense_condensed_direction` (factor
/// failure and KKT solve failure paths) and in
/// `solve_for_search_direction`'s full-augmented error path.
///
/// On restoration success, advances `state.x`, zeros `alpha_primal`,
/// re-evaluates the linear-aware step and updates the L-BFGS Hessian
/// estimate. On failure, returns a `NumericalError` `SolveResult`.
#[allow(clippy::too_many_arguments)]
fn restore_after_solve_failure<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    n: usize,
    m: usize,
    filter: &Filter,
    restoration: &mut RestorationPhase,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    deadline: Option<Instant>,
) -> SolveRestoreOutcome {
    let (x_rest, success) = restoration.restore(
        &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
        &state.jac_rows, &state.jac_cols, n, m, options,
        &|theta, phi| filter.is_acceptable(theta, phi),
        &|x_eval, g_out| problem.constraints(x_eval, true, g_out),
        &|x_eval, jac_out| problem.jacobian_values(x_eval, true, jac_out),
        Some(&|x_eval: &[f64], obj_out: &mut f64| problem.objective(x_eval, true, obj_out)),
        deadline,
    );
    if success {
        state.x = x_rest;
        state.alpha_primal = 0.0;
        let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        return SolveRestoreOutcome::Continue;
    }
    SolveRestoreOutcome::Return(make_result(state, SolveStatus::NumericalError))
}

/// Mehrotra PC deflection revert: solve the original (μ-current)
/// RHS via iterative refinement and check the angle between the
/// PC step `dx_dir` and the original step `dx_orig`. If
/// cos(angle) < 0.7 (i.e. >~45° deflection), revert to the
/// original direction, restore the original RHS in
/// `kkt_system_opt`, and clear `mehrotra_aff`. This guards against
/// the predictor probe pushing the search direction far from a
/// reasonable Newton step in early iterations.
fn maybe_revert_mehrotra_deflection(
    dx_dir: &mut Vec<f64>,
    dy_dir: &mut Vec<f64>,
    mehrotra_aff: &mut Option<(Vec<f64>, Vec<f64>, Vec<f64>, f64)>,
    saved_rhs: &Option<Vec<f64>>,
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
) {
    let Some(orig_rhs) = saved_rhs.as_ref() else { return };
    let Some(kkt) = kkt_system_opt.as_ref() else { return };
    let Ok((dx_orig, dy_orig)) = kkt::solve_with_custom_rhs_refined(
        &kkt.matrix, kkt.n, kkt.dim, lin_solver, orig_rhs,
    ) else { return };
    let norm_orig: f64 = dx_orig.iter().map(|v| v * v).sum::<f64>().sqrt();
    let norm_pc: f64 = dx_dir.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm_orig <= 1e-30 || norm_pc <= 1e-30 {
        return;
    }
    let dot: f64 = dx_orig.iter().zip(dx_dir.iter()).map(|(a, b)| a * b).sum::<f64>();
    let cos_angle = dot / (norm_orig * norm_pc);
    if cos_angle >= 0.7 {
        return;
    }
    log::debug!(
        "Mehrotra PC deflection too large (cos={:.3}), reverting",
        cos_angle
    );
    *dx_dir = dx_orig;
    *dy_dir = dy_orig;
    kkt_system_opt.as_mut().unwrap().rhs = orig_rhs.clone();
    *mehrotra_aff = None;
}

/// Ipopt-style quality escalation around `kkt::solve_for_direction`
/// for the full augmented KKT path. After the initial solve, on a
/// `PretendSingular` error walks the escalation ladder:
///
///   1. Ruiz equilibrate (if not already scaled);
///   2. raise the linear solver's pivot tolerance;
///   3. apply δ_c constraint regularization (if `m > 0`);
///   4. apply δ_w Hessian regularization.
///
/// Each rung re-factors and re-solves; the loop exits once the
/// solve no longer reports `PretendSingular`. Returns the final
/// `dir_result` and the (possibly increased) `(ic_delta_w, ic_delta_c)`.
fn solve_with_quality_escalation(
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    inertia_params: &mut InertiaCorrectionParams,
    mut ic_delta_w: f64,
    mut ic_delta_c: f64,
    n: usize,
    m: usize,
) -> (Result<(Vec<f64>, Vec<f64>), crate::linear_solver::SolverError>, f64, f64) {
    let mut dir_result = kkt::solve_for_direction(
        kkt_system_opt.as_ref().unwrap(), lin_solver, ic_delta_w, ic_delta_c,
    );
    if matches!(dir_result, Err(crate::linear_solver::SolverError::PretendSingular)) {
        if let Some(kkt_system) = kkt_system_opt.as_mut() {
            let mut ps_resolved = false;

            if !ps_resolved {
                if !inertia_params.use_scaling && kkt_system.scale_factors.is_none() {
                    inertia_params.use_scaling = true;
                    let scale = kkt::ruiz_equilibrate(&mut kkt_system.matrix, &mut kkt_system.rhs);
                    kkt_system.scale_factors = Some(scale);
                    if lin_solver.factor(&kkt_system.matrix).is_ok() {
                        dir_result = kkt::solve_for_direction(kkt_system, lin_solver, ic_delta_w, ic_delta_c);
                        ps_resolved = !matches!(dir_result, Err(crate::linear_solver::SolverError::PretendSingular));
                    }
                }
                if !ps_resolved && lin_solver.increase_quality() {
                    if lin_solver.factor(&kkt_system.matrix).is_ok() {
                        dir_result = kkt::solve_for_direction(kkt_system, lin_solver, ic_delta_w, ic_delta_c);
                        ps_resolved = !matches!(dir_result, Err(crate::linear_solver::SolverError::PretendSingular));
                    }
                }
            }

            if !ps_resolved && m > 0 && ic_delta_c == 0.0 {
                let dc = inertia_params.delta_c_base;
                let mut perturbed = kkt_system.matrix.clone();
                perturbed.add_diagonal_range(n, n + m, -dc);
                if lin_solver.factor(&perturbed).is_ok() {
                    kkt_system.matrix = perturbed;
                    ic_delta_c = dc;
                    dir_result = kkt::solve_for_direction(kkt_system, lin_solver, ic_delta_w, ic_delta_c);
                    ps_resolved = !matches!(dir_result, Err(crate::linear_solver::SolverError::PretendSingular));
                }
            }
            if !ps_resolved && matches!(dir_result, Err(crate::linear_solver::SolverError::PretendSingular)) {
                let dw = if ic_delta_w == 0.0 { inertia_params.delta_w_init } else { ic_delta_w * inertia_params.delta_w_growth };
                let dc = if m > 0 && ic_delta_c == 0.0 { inertia_params.delta_c_base } else { ic_delta_c };
                let mut perturbed = kkt_system.matrix.clone();
                perturbed.add_diagonal_range(0, n, dw);
                if m > 0 && dc > ic_delta_c {
                    perturbed.add_diagonal_range(n, n + m, -(dc - ic_delta_c));
                }
                if lin_solver.factor(&perturbed).is_ok() {
                    kkt_system.matrix = perturbed;
                    ic_delta_w = dw;
                    ic_delta_c = dc;
                    inertia_params.delta_w_last = dw;
                    dir_result = kkt::solve_for_direction(kkt_system, lin_solver, ic_delta_w, ic_delta_c);
                }
            }
            if !inertia_params.use_scaling {
                inertia_params.use_scaling = true;
            }
        }
    }
    (dir_result, ic_delta_w, ic_delta_c)
}

/// Mehrotra predictor-corrector probe: solve the affine-predictor
/// system (rhs zeroed of barrier terms), recover bound-multiplier
/// directions via `kkt::recover_dz` at μ=0, take the
/// fraction-to-boundary step α_aff under τ=1-1e-3, derive the
/// affine complementarity μ_aff and centering parameter
/// σ = (μ_aff/μ)^3 ∈ [0,1], and rebuild the RHS at
/// μ_pc = max(σ·μ, μ_min). Returns the rebuilt RHS together with
/// the affine direction (dx_aff, dz_l_aff, dz_u_aff) and μ_pc.
/// Returns None when the affine solve fails or no bound
/// constraints contributed to μ_aff.
fn try_mehrotra_predictor(
    state: &SolverState,
    options: &SolverOptions,
    kkt: &kkt::KktSystem,
    lin_solver: &mut dyn LinearSolver,
    iteration: usize,
    n: usize,
    last_mehrotra_sigma: &mut Option<f64>,
) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, f64)> {
    let rhs_aff = kkt::affine_predictor_rhs(
        &kkt.rhs, &state.x, &state.x_l, &state.x_u, state.mu,
    );
    let (dx_aff, _) = kkt::solve_with_custom_rhs_refined(
        &kkt.matrix, kkt.n, kkt.dim, lin_solver, &rhs_aff,
    ).ok()?;
    let (dz_l_aff, dz_u_aff) = kkt::recover_dz(
        &state.x, &state.x_l, &state.x_u,
        &state.z_l, &state.z_u, &dx_aff, 0.0,
    );
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
    if nb == 0 {
        return None;
    }
    let mu_aff = mu_aff_sum / nb as f64;
    let sigma_mehr = (mu_aff / state.mu).powi(3).clamp(0.0, 1.0);
    *last_mehrotra_sigma = Some(sigma_mehr);
    let mu_pc = (sigma_mehr * state.mu).max(options.mu_min);
    log::debug!(
        "Mehrotra PC iter {}: σ={:.4} α_aff={:.4} μ: {:.2e}→{:.2e}",
        iteration, sigma_mehr, alpha_aff, state.mu, mu_pc
    );
    let new_rhs = kkt::rebuild_rhs_with_mu(
        &kkt.rhs, &state.x, &state.x_l, &state.x_u,
        state.mu, mu_pc,
    );
    Some((new_rhs, dx_aff, dz_l_aff, dz_u_aff, mu_pc))
}

/// Sparse condensed direction solve: factor the sparse Schur
/// complement S with the banded/sparse solver. On factor or
/// solve failure rebuilds the full augmented KKT, factors it with
/// inertia correction (using a fresh fallback solver), and falls
/// back to a gradient-descent step if even that fails. Always
/// returns a (dx, dy) pair — restoration is never invoked here
/// because gradient_descent_fallback acts as the floor.
#[allow(clippy::too_many_arguments)]
fn solve_sparse_condensed_direction(
    state: &SolverState,
    sc: &kkt::SparseCondensedKktSystem,
    sigma: &[f64],
    n: usize,
    m: usize,
    use_sparse: bool,
    lin_solver: &mut dyn LinearSolver,
    inertia_params: &mut InertiaCorrectionParams,
) -> (Vec<f64>, Vec<f64>) {
    let kkt_sc = KktMatrix::Sparse(sc.matrix.clone());
    let factor_ok = lin_solver.factor(&kkt_sc).is_ok();
    if factor_ok {
        match kkt::solve_sparse_condensed(sc, lin_solver) {
            Ok(d) => return d,
            Err(_) => {}
        }
    }
    let mut kkt = kkt::assemble_kkt(
        n, m,
        &state.hess_rows, &state.hess_cols, &state.hess_vals,
        &state.jac_rows, &state.jac_cols, &state.jac_vals,
        sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
        &state.y, &state.z_l, &state.z_u,
        &state.x, &state.x_l, &state.x_u, state.mu,
        use_sparse, &state.v_l, &state.v_u,
    );
    let mut fallback_solver = new_fallback_solver(use_sparse);
    if let Ok((fb_dw, fb_dc)) = kkt::factor_with_inertia_correction(
        &mut kkt, fallback_solver.as_mut(), inertia_params,
    ) {
        kkt::solve_for_direction(&kkt, fallback_solver.as_mut(), fb_dw, fb_dc)
            .unwrap_or_else(|_| gradient_descent_fallback(state)
                .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m])))
    } else {
        gradient_descent_fallback(state)
            .unwrap_or_else(|| (vec![0.0; n], vec![0.0; m]))
    }
}

/// Dense condensed direction solve: Bunch-Kaufman on the n×n Schur
/// complement `H + Σ + J^T·D_c^{-1}·J`. On BK or solve failure rebuilds
/// the full augmented KKT, factors it with inertia correction, and
/// solves; falls back to restoration on persistent failure.
#[allow(clippy::too_many_arguments)]
fn solve_dense_condensed_direction<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    n: usize,
    m: usize,
    use_sparse: bool,
    cond: &kkt::CondensedKktSystem,
    sigma: &[f64],
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    inertia_params: &mut InertiaCorrectionParams,
    filter: &Filter,
    restoration: &mut RestorationPhase,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    deadline: Option<Instant>,
) -> CondensedDirectionOutcome {
    let mut cond_solver = DenseLdl::new();
    let t_cond_bk = Instant::now();
    let cond_ok = cond_solver.bunch_kaufman_factor(&cond.matrix).is_ok();
    if options.print_level >= 5 {
        rip_log!("ripopt: Dense condensed BK factor n={}: {:.3}s (ok={})",
            n, t_cond_bk.elapsed().as_secs_f64(), cond_ok);
    }
    let cond_result = if cond_ok {
        kkt::solve_condensed(cond, &mut cond_solver).ok()
    } else {
        None
    };

    if let Some((dx, dy)) = cond_result {
        return CondensedDirectionOutcome::Solved {
            dx,
            dy,
            cond_solver: Some(cond_solver),
        };
    }

    // Condensed failed — build full KKT on demand.
    let mut kkt = kkt::assemble_kkt(
        n, m,
        &state.hess_rows, &state.hess_cols, &state.hess_vals,
        &state.jac_rows, &state.jac_cols, &state.jac_vals,
        sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
        &state.y, &state.z_l, &state.z_u,
        &state.x, &state.x_l, &state.x_u, state.mu,
        use_sparse, &state.v_l, &state.v_u,
    );
    let fb_ic = kkt::factor_with_inertia_correction(
        &mut kkt, lin_solver, inertia_params,
    );
    if fb_ic.is_err() {
        return match restore_after_solve_failure(
            state, problem, options, n, m, filter, restoration,
            lbfgs_state, lbfgs_mode, linear_constraints, deadline,
        ) {
            SolveRestoreOutcome::Continue => CondensedDirectionOutcome::Continue,
            SolveRestoreOutcome::Return(r) => CondensedDirectionOutcome::Return(r),
        };
    }
    let (fb_dw, fb_dc) = fb_ic.unwrap();
    match kkt::solve_for_direction(&kkt, lin_solver, fb_dw, fb_dc) {
        Ok((dx, dy)) => {
            *kkt_system_opt = Some(kkt);
            CondensedDirectionOutcome::Solved {
                dx,
                dy,
                cond_solver: None,
            }
        }
        Err(e) => {
            log::warn!("KKT solve failed: {}", e);
            match restore_after_solve_failure(
                state, problem, options, n, m, filter, restoration,
                lbfgs_state, lbfgs_mode, linear_constraints, deadline,
            ) {
                SolveRestoreOutcome::Continue => CondensedDirectionOutcome::Continue,
                SolveRestoreOutcome::Return(r) => CondensedDirectionOutcome::Return(r),
            }
        }
    }
}

/// Solve the KKT system for the primal-dual Newton search direction.
///
/// Dispatches over the three KKT representations:
///
///  * **Dense condensed** (`m ≥ 2n` and small `n`) — Bunch-Kaufman on the
///    `n×n` Schur complement `H + Σ + J^T·D_c^{-1}·J`. On BK failure
///    builds the full augmented KKT on demand and retries.
///  * **Sparse condensed** — sparse factor of the Schur complement.
///    On factor/solve failure falls back to the full augmented KKT
///    with `gradient_descent_fallback` as last resort.
///  * **Full augmented** (`else` branch) — optional Mehrotra
///    predictor-corrector probe to update μ, then Ipopt-style quality
///    escalation (Ruiz scaling → raise pivot tolerance → δ_c → δ_w),
///    with 30° deflection revert if the PC direction strays too far
///    from the original.
///
/// Mirrors `IpPDSearchDirCalc.cpp:81-110` / `IpPDFullSpaceSolver.cpp`
/// in Ipopt.
#[allow(clippy::too_many_arguments)]
fn solve_for_search_direction<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    iteration: usize,
    n: usize,
    m: usize,
    use_sparse: bool,
    condensed_system: &Option<kkt::CondensedKktSystem>,
    sparse_condensed_system: &Option<kkt::SparseCondensedKktSystem>,
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    sigma: &[f64],
    inertia_params: &mut InertiaCorrectionParams,
    ic_delta_w: f64,
    ic_delta_c: f64,
    filter: &Filter,
    restoration: &mut RestorationPhase,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    last_mehrotra_sigma: &mut Option<f64>,
    deadline: Option<Instant>,
) -> DirectionSolveDecision {
    let mut cond_solver_for_soc: Option<DenseLdl> = None;
    let mut mehrotra_aff: Option<(Vec<f64>, Vec<f64>, Vec<f64>, f64)> = None;

    let (dx, dy) = if let Some(cond) = condensed_system.as_ref() {
        match solve_dense_condensed_direction(
            state, problem, options, n, m, use_sparse, cond, sigma,
            kkt_system_opt, lin_solver, inertia_params, filter, restoration,
            lbfgs_state, lbfgs_mode, linear_constraints, deadline,
        ) {
            CondensedDirectionOutcome::Solved { dx, dy, cond_solver } => {
                cond_solver_for_soc = cond_solver;
                (dx, dy)
            }
            CondensedDirectionOutcome::Continue => return DirectionSolveDecision::Continue,
            CondensedDirectionOutcome::Return(r) => return DirectionSolveDecision::Return(r),
        }
    } else if let Some(sc) = sparse_condensed_system.as_ref() {
        solve_sparse_condensed_direction(
            state, sc, sigma, n, m, use_sparse, lin_solver, inertia_params,
        )
    } else {
        // Full augmented KKT with optional Mehrotra PC.
        let saved_rhs = if options.mehrotra_pc {
            kkt_system_opt.as_ref().map(|k| k.rhs.clone())
        } else {
            None
        };
        let mut mehrotra_applied = false;

        if options.mehrotra_pc {
            let has_bounds = (0..n).any(|i| state.x_l[i].is_finite() || state.x_u[i].is_finite());
            if has_bounds {
                let pc_result = try_mehrotra_predictor(
                    state, options, kkt_system_opt.as_ref().unwrap(), lin_solver,
                    iteration, n, last_mehrotra_sigma,
                );
                if let Some((new_rhs, dx_aff_v, dz_l_aff_v, dz_u_aff_v, mu_pc_used)) = pc_result {
                    kkt_system_opt.as_mut().unwrap().rhs = new_rhs;
                    mehrotra_applied = true;
                    mehrotra_aff = Some((dx_aff_v, dz_l_aff_v, dz_u_aff_v, mu_pc_used));
                }
            }
        }

        // Escalated (δ_w, δ_c) values are not read again — the helper
        // mutates `kkt_system.matrix`, `inertia_params`, and `lin_solver`
        // in place, which carries the escalation effect into subsequent
        // iterations.
        let (dir_result, _, _) = solve_with_quality_escalation(
            kkt_system_opt, lin_solver, inertia_params, ic_delta_w, ic_delta_c, n, m,
        );
        let (mut dx_dir, mut dy_dir) = match dir_result {
            Ok(d) => d,
            Err(e) => {
                log::warn!("KKT solve failed: {}", e);
                if let Some(fallback) = gradient_descent_fallback(state) {
                    fallback
                } else {
                    match restore_after_solve_failure(
                        state, problem, options, n, m, filter, restoration,
                        lbfgs_state, lbfgs_mode, linear_constraints, deadline,
                    ) {
                        SolveRestoreOutcome::Continue => return DirectionSolveDecision::Continue,
                        SolveRestoreOutcome::Return(r) => return DirectionSolveDecision::Return(r),
                    }
                }
            }
        };

        if mehrotra_applied {
            maybe_revert_mehrotra_deflection(
                &mut dx_dir, &mut dy_dir, &mut mehrotra_aff,
                &saved_rhs, kkt_system_opt, lin_solver,
            );
        }

        (dx_dir, dy_dir)
    };

    DirectionSolveDecision::Proceed {
        dx,
        dy,
        cond_solver_for_soc,
        mehrotra_aff,
    }
}

/// First-iteration bandwidth detection for the sparse condensed Schur
/// complement. On `iteration == 0` with `use_sparse_condensed`, measure
/// the bandwidth of S and either:
///   - `bw > n/2`: abandon sparse condensed, rebuild the full augmented
///     KKT via `kkt::assemble_kkt`, and keep the sparse solver;
///   - `bw*bw <= n`: swap in a `BandedLdl` solver;
///   - otherwise: keep the current sparse solver.
///
/// The dense-condensed path has no bandwidth concept, so this is a
/// no-op on that path.
#[allow(clippy::too_many_arguments)]
fn adjust_sparse_condensed_bandwidth<P: NlpProblem>(
    state: &SolverState,
    _problem: &P,
    options: &SolverOptions,
    n: usize,
    m: usize,
    use_sparse: bool,
    use_sparse_condensed: bool,
    sigma: &[f64],
    sparse_condensed_system: &mut Option<kkt::SparseCondensedKktSystem>,
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut Box<dyn LinearSolver>,
    disable_sparse_condensed: &mut bool,
) {
    if !use_sparse_condensed {
        return;
    }
    let sc_bw = sparse_condensed_system.as_ref().map(|sc| {
        BandedLdl::compute_bandwidth(&sc.matrix.triplet_rows, &sc.matrix.triplet_cols)
    });
    let Some(bw) = sc_bw else { return };

    if bw > n / 2 {
        // Condensed Schur complement is essentially dense — switch to full
        // augmented KKT with the sparse solver (rmumps/AMD/ND handles this).
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Sparse condensed S has bandwidth {} for n={}, switching to dense condensed KKT",
                bw, n
            );
        }
        *disable_sparse_condensed = true;
        *sparse_condensed_system = None;
        *kkt_system_opt = Some(kkt::assemble_kkt(
            n, m,
            &state.hess_rows, &state.hess_cols, &state.hess_vals,
            &state.jac_rows, &state.jac_cols, &state.jac_vals,
            sigma, &state.grad_f, &state.g, &state.g_l, &state.g_u,
            &state.y, &state.z_l, &state.z_u,
            &state.x, &state.x_l, &state.x_u, state.mu,
            use_sparse, &state.v_l, &state.v_u,
        ));
    } else if bw * bw <= n {
        if options.print_level >= 5 {
            rip_log!("ripopt: Sparse condensed S has bandwidth {} for n={}, using banded solver", bw, n);
        }
        *lin_solver = Box::new(BandedLdl::new());
    } else if options.print_level >= 5 {
        rip_log!("ripopt: Sparse condensed S has bandwidth {} for n={}, using sparse solver", bw, n);
    }
}

/// Post-step "acceptable" tracking. Mirrors the pre-step acceptable check
/// (see `track_consecutive_acceptable`) but is evaluated at the freshly
/// accepted iterate — catches cases where the step just taken pushes
/// the state into the acceptable region but the pre-step check at the
/// top of the loop missed it.
///
/// Increments `state.consecutive_acceptable` iff the accepted iterate
/// passes both the scaled (`options.tol`) and unscaled
/// (`constr_viol_tol` / `dual_inf_tol` / `compl_inf_tol`) Ipopt
/// acceptable thresholds. Never resets the counter here; pre-step
/// handling owns resets.
fn track_post_step_acceptable(state: &mut SolverState, options: &SolverOptions) {
    let n = state.n;
    let m = state.m;
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
    let post_du_unsc = convergence::dual_infeasibility_scaled(
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
    let post_near_scaled = post_primal <= 100.0 * options.tol
        && post_du <= 100.0 * options.tol * post_sd
        && post_compl_best <= 100.0 * options.tol * post_sd;
    let post_near_unscaled = post_primal <= 10.0 * options.constr_viol_tol
        && post_du_unsc <= 10.0 * options.dual_inf_tol
        && post_compl_best <= 10.0 * options.compl_inf_tol;
    if post_near_scaled && post_near_unscaled {
        state.consecutive_acceptable += 1;
    }
}

/// Outcome of the post-LS restoration cascade.
enum RestorationCascadeDecision {
    Continue,
    Return(SolveResult),
}

/// Run the fast Gauss–Newton restoration. Returns true if GN restoration
/// succeeded and `apply_restoration_success` was invoked (caller should
/// short-circuit the cascade with `Continue`); false otherwise (caller
/// proceeds to the recovery / NLP-restoration fallbacks).
#[allow(clippy::too_many_arguments)]
fn try_gn_restoration<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    filter: &mut Filter,
    mu_state: &mut MuState,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    restoration: &mut RestorationPhase,
    n: usize,
    m: usize,
    deadline: Option<Instant>,
) -> bool {
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
        &|x_eval, g_out| problem.constraints(x_eval, true, g_out),
        &|x_eval, jac_out| problem.jacobian_values(x_eval, true, jac_out),
        Some(&|x_eval: &[f64], obj_out: &mut f64| problem.objective(x_eval, true, obj_out)),
        deadline,
    );

    if gn_success {
        state.diagnostics.restoration_count += 1;
        apply_restoration_success(
            state, filter, mu_state, options, n, m, problem, &x_rest,
            linear_constraints, lbfgs_mode, lbfgs_state,
        );
        true
    } else {
        false
    }
}

/// Adjust mu / mu-mode and (on attempt 3+) jitter x to escape the current
/// basin after restoration failed but the cascade has not exhausted its
/// retry budget. fail_count == 1 switches Free → Fixed mode (or applies
/// the standard mu_linear/superlinear decrease in Fixed mode); subsequent
/// attempts cycle through the perturbation sequence ×10, ×0.1, ×100,
/// ×0.01, ×1000, ×0.001. The filter is reset and inertia δ_w cleared.
#[allow(clippy::too_many_arguments)]
fn apply_restoration_recovery_strategy<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    filter: &mut Filter,
    mu_state: &mut MuState,
    inertia_params: &mut InertiaCorrectionParams,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    fail_count: usize,
    n: usize,
) {
    log::debug!("Restoration failed (attempt #{}), trying recovery", fail_count);
    let mu_factors: [f64; 6] = [10.0, 0.1, 100.0, 0.01, 1000.0, 0.001];

    match fail_count {
        1 => {
            if mu_state.mode == MuMode::Free {
                state.diagnostics.mu_mode_switches += 1;
                mu_state.mode = MuMode::Fixed;
                mu_state.first_iter_in_mode = true;
                let avg_compl = compute_avg_complementarity(state);
                if avg_compl > 0.0 {
                    state.mu = (options.adaptive_mu_monotone_init_factor * avg_compl)
                        .clamp(options.mu_min, 1e5);
                } else {
                    state.mu = (options.mu_linear_decrease_factor * state.mu)
                        .max(options.mu_min);
                }
            } else {
                let new_mu = (options.mu_linear_decrease_factor * state.mu)
                    .min(state.mu.powf(options.mu_superlinear_decrease_power))
                    .max(options.mu_min);
                state.mu = new_mu;
            }
        }
        _ => {
            let factor = mu_factors[(fail_count - 2) % mu_factors.len()];
            state.mu = (state.mu * factor).max(options.mu_min).min(1e5);
        }
    }
    filter.reset();
    let theta_now = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_now);
    inertia_params.delta_w_last = 0.0;

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
        let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
    }
}

/// Classify the terminal status when the restoration cascade has exhausted its
/// retry budget. Returns LocalInfeasibility when the constraint-violation
/// gradient is near-stationary at an infeasible point, Infeasible when theta
/// has stayed > 1% of its historical minimum well into the run, otherwise
/// RestorationFailed. Caller must only invoke this once `fail_count` exceeds
/// `max_restore_attempts`.
fn classify_exhausted_restoration_attempt(
    state: &SolverState,
    options: &SolverOptions,
    iteration: usize,
    fail_count: usize,
    n: usize,
    m: usize,
    ever_feasible: bool,
    theta_history: &[f64],
    theta_history_len: usize,
) -> SolveResult {
    log::warn!(
        "Restoration failed at iteration {} (attempt #{})",
        iteration, fail_count
    );
    let current_theta = state.constraint_violation();

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
            return make_result(state, SolveStatus::LocalInfeasibility);
        }
    }

    if !ever_feasible && current_theta > 1e4 && iteration > 500
        && theta_history.len() >= theta_history_len
    {
        let min_theta = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
        if current_theta > 0.01 * min_theta {
            return make_result(state, SolveStatus::Infeasible);
        }
    }
    make_result(state, SolveStatus::RestorationFailed)
}

/// Invoke the restoration cascade after a failed line search. Runs the
/// fast Gauss–Newton restoration first; on failure, rotates through
/// (a) NLP restoration at `fail_count ∈ {2, 4}`,
/// (b) infeasibility / `RestorationFailed` / `MaxIterations` exits,
/// (c) mu-mode switching and mu perturbation,
/// (d) x-perturbation on `fail_count >= 3`.
///
/// Always ends with either `Continue` (the main loop should re-enter its
/// top with an updated iterate / recovery state) or `Return(SolveResult)`
/// (the solver is terminating).
#[allow(clippy::too_many_arguments)]
fn run_post_ls_restoration_cascade<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    filter: &mut Filter,
    mu_state: &mut MuState,
    inertia_params: &mut InertiaCorrectionParams,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    lbfgs_mode: bool,
    linear_constraints: Option<&[bool]>,
    restoration: &mut RestorationPhase,
    iteration: usize,
    n: usize,
    m: usize,
    start_time: Instant,
    deadline: Option<Instant>,
    early_timeout: f64,
    ever_feasible: bool,
    theta_current: f64,
    phi_current: f64,
    theta_history: &[f64],
    theta_history_len: usize,
) -> RestorationCascadeDecision {
    state.diagnostics.filter_rejects += 1;

    // Add current point to filter before entering restoration (Ipopt convention).
    // augment_for_restoration adds the margin entry
    // (phi - gamma_phi*theta, (1-gamma_theta)*theta) AND bumps theta_max —
    // this prevents restoration from handing back a point as bad as the entry.
    filter.add(theta_current, phi_current);
    filter.augment_for_restoration(theta_current, phi_current);

    // Phase 1: Fast GN restoration
    log::debug!("Line search failed at iteration {}, entering restoration", iteration);

    if try_gn_restoration(
        state, problem, options, filter, mu_state, lbfgs_state, lbfgs_mode,
        linear_constraints, restoration, n, m, deadline,
    ) {
        return RestorationCascadeDecision::Continue;
    }

    // GN restoration failed — recovery logic with NLP restoration as last resort.
    // Bail out of recovery cascade if wall time is nearly exhausted.
    if options.max_wall_time > 0.0 {
        let remaining = options.max_wall_time - start_time.elapsed().as_secs_f64();
        if remaining < 1.0 {
            return RestorationCascadeDecision::Return(make_result(state, SolveStatus::MaxIterations));
        }
    }

    mu_state.consecutive_restoration_failures += 1;
    let fail_count = mu_state.consecutive_restoration_failures;

    // At fail_count == 2 (or 4): try full NLP restoration. It's Ipopt's primary
    // method — most robust, but doubles the problem size, so skip for large n+m.
    let kkt_dim = n + m;
    let skip_nlp_restoration = iteration < 3
        && options.early_stall_timeout > 0.0
        && start_time.elapsed().as_secs_f64() > early_timeout * 0.5;
    if (fail_count == 2 || fail_count == 4)
        && !options.disable_nlp_restoration
        && kkt_dim <= 50000
        && !skip_nlp_restoration
    {
        state.diagnostics.nlp_restoration_count += 1;
        let (x_nlp, outcome) = attempt_nlp_restoration(
            problem, state, filter, options, theta_current, start_time,
        );
        match outcome {
            RestorationOutcome::Success => {
                apply_restoration_success(
                    state, filter, mu_state, options, n, m,
                    problem, &x_nlp,
                    linear_constraints, lbfgs_mode, lbfgs_state,
                );
                return RestorationCascadeDecision::Continue;
            }
            RestorationOutcome::LocalInfeasibility | RestorationOutcome::Failed => {
                // Fall through to continue recovery. Don't immediately return
                // LocalInfeasibility — the fail_count > 6 stationarity check below
                // is more reliable.
            }
        }
    }

    // For large problems (no NLP restoration), give up sooner.
    let max_restore_attempts = if kkt_dim > 50000 { 3 } else { 6 };
    if fail_count > max_restore_attempts {
        return RestorationCascadeDecision::Return(classify_exhausted_restoration_attempt(
            state, options, iteration, fail_count, n, m, ever_feasible,
            theta_history, theta_history_len,
        ));
    }

    apply_restoration_recovery_strategy(
        state, problem, options, filter, mu_state, inertia_params,
        lbfgs_state, lbfgs_mode, linear_constraints, fail_count, n,
    );

    RestorationCascadeDecision::Continue
}

/// Outcome of the post-step re-evaluation with α-halving recovery.
enum PostStepEvalDecision {
    Proceed,
    Continue,
    Return(SolveResult),
}

/// Re-evaluate the problem at the accepted iterate and, on evaluation
/// failure, execute the Ipopt-style α-halving → restoration recovery
/// cascade.
///
/// Ipopt's `IpBacktrackingLineSearch.cpp:776-784` treats a post-step
/// `Eval_Error` as an α backtrack rather than a fatal. Mirror that by
/// halving α and α_dual (from the accepted point back toward
/// `x_pre_step`) up to 5 times; if that fails, invoke the restoration
/// phase; if restoration also fails, return `NumericalError`.
///
/// On success (or successful recovery via α-halving or restoration)
/// returns `Proceed` / `Continue` for the main loop's control flow.
fn reevaluate_after_step<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    filter: &mut Filter,
    restoration: &mut RestorationPhase,
    timings: &mut PhaseTimings,
    x_pre_step: &[f64],
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
    deadline: Option<Instant>,
) -> PostStepEvalDecision {
    let n = state.n;
    let m = state.m;
    let t_eval = Instant::now();
    let eval_ok = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    if eval_ok {
        update_lbfgs_hessian(lbfgs_state, state);
    }
    timings.problem_eval += t_eval.elapsed();

    if eval_ok {
        return PostStepEvalDecision::Proceed;
    }

    // Post-step evaluation failure: α-halving.
    let mut retry_alpha = state.alpha_primal;
    let mut retry_alpha_dual = state.alpha_dual;
    for _ in 0..5 {
        retry_alpha *= 0.5;
        retry_alpha_dual *= 0.5;
        for i in 0..n {
            state.x[i] = x_pre_step[i] + retry_alpha * state.dx[i];
        }
        if state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode) {
            state.alpha_primal = retry_alpha;
            state.alpha_dual = retry_alpha_dual;
            update_lbfgs_hessian(lbfgs_state, state);
            return PostStepEvalDecision::Continue;
        }
    }

    // α halving exhausted: try restoration.
    let (x_rest, success) = restoration.restore(
        &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
        &state.jac_rows, &state.jac_cols, n, m, options,
        &|theta, phi| filter.is_acceptable(theta, phi),
        &|x_eval, g_out| problem.constraints(x_eval, true, g_out),
        &|x_eval, jac_out| problem.jacobian_values(x_eval, true, jac_out),
        Some(&|x_eval: &[f64], obj_out: &mut f64| problem.objective(x_eval, true, obj_out)),
        deadline,
    );
    if success {
        state.x = x_rest;
        state.alpha_primal = 0.0;
        if state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode) {
            update_lbfgs_hessian(lbfgs_state, state);
            return PostStepEvalDecision::Continue;
        }
    }

    PostStepEvalDecision::Return(make_result(state, SolveStatus::NumericalError))
}

/// Reset the slack-constraint multipliers `v_l`, `v_u` from the
/// barrier equilibrium `v = mu_ks / slack` after the post-step
/// re-evaluation.
///
/// A simple reset rather than a Newton update — the dv direction is
/// approximate (we carry no explicit slacks) and applying FTB on `v`
/// can restrict `alpha_d` too much. Only active slacks (`v > 0`) are
/// touched.
fn reset_slack_multipliers(state: &mut SolverState, mu_ks: f64) {
    let m = state.m;
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
}

/// Ipopt's TrySoftRestoStep: before invoking the expensive GN restoration,
/// try the current search direction with `alpha = min(alpha_primal_max, alpha_dual_max)`
/// and accept it if the constraint violation drops by at least `1e-4` relative
/// and the trial passes the filter. Returns `true` if the step was taken.
fn attempt_soft_restoration<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    filter: &mut Filter,
    alpha_primal_max: f64,
    alpha_dual_max: f64,
    theta_current: f64,
    phi_current: f64,
) -> bool {
    let n = state.n;
    let m = state.m;
    let soft_alpha = alpha_primal_max.min(alpha_dual_max);
    let mut x_soft = vec![0.0; n];
    for i in 0..n {
        x_soft[i] = state.x[i] + soft_alpha * state.dx[i];
        if state.x_l[i].is_finite() {
            x_soft[i] = x_soft[i].max(state.x_l[i] + 1e-14);
        }
        if state.x_u[i].is_finite() {
            x_soft[i] = x_soft[i].min(state.x_u[i] - 1e-14);
        }
    }
    let mut obj_soft = 0.0;
    let obj_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        problem.objective(&x_soft, true, &mut obj_soft)
    }));
    if !matches!(obj_result, Ok(true)) {
        return false;
    }
    let mut g_soft = vec![0.0; m];
    if !problem.constraints(&x_soft, false, &mut g_soft) {
        return false;
    }
    let theta_soft = convergence::primal_infeasibility(&g_soft, &state.g_l, &state.g_u);
    if theta_soft < (1.0 - 1e-4) * theta_current
        && filter.is_acceptable(theta_soft, obj_soft)
    {
        log::debug!(
            "Soft restoration accepted: theta {:.2e} -> {:.2e}",
            theta_current,
            theta_soft
        );
        state.x = x_soft;
        state.obj = obj_soft;
        state.g = g_soft;
        state.alpha_primal = soft_alpha;
        filter.add(theta_current, phi_current);
        true
    } else {
        false
    }
}

/// Record the current iterate as the best-feasible point seen so far
/// if it satisfies `constr_viol_tol` and strictly improves `best_obj`.
/// Used as the fallback iterate returned at `max_iter` exit.
fn track_best_feasible(
    state: &SolverState,
    options: &SolverOptions,
    best_obj: &mut f64,
    best_x: &mut Option<Vec<f64>>,
) {
    let theta_now = state.constraint_violation();
    if theta_now < options.constr_viol_tol && state.obj < *best_obj {
        *best_obj = state.obj;
        *best_x = Some(state.x.clone());
    }
}

/// Fraction-to-boundary step limits for primal and dual.
///
/// `tau` = `max(1 - NLP_error, tau_min)` in Free mode or
/// `max(1 - mu, tau_min)` in Fixed mode. Primal bound:
/// `x + alpha*dx > x_l` and `< x_u` by a factor of `tau`.
/// Dual bound: `z + alpha*dz > 0` by a factor of `tau`.
fn compute_alpha_max(
    state: &SolverState,
    options: &SolverOptions,
    mu_state: &MuState,
    primal_inf: f64,
    dual_inf: f64,
    compl_inf: f64,
) -> (f64, f64, f64) {
    let n = state.n;
    let tau = if mu_state.mode == MuMode::Free {
        let nlp_error = primal_inf + dual_inf + compl_inf;
        (1.0 - nlp_error).max(options.tau_min)
    } else {
        (1.0 - state.mu).max(options.tau_min)
    };

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

    let alpha_dual_max_l = filter::fraction_to_boundary(&state.z_l, &state.dz_l, tau);
    let alpha_dual_max_u = filter::fraction_to_boundary(&state.z_u, &state.dz_u, tau);
    let alpha_dual_max = alpha_dual_max_l.min(alpha_dual_max_u);

    (tau, alpha_primal_max, alpha_dual_max)
}

/// Ipopt-style tiny-step detection: when the relative step size is
/// below `10*eps` for two consecutive iterations and primal
/// infeasibility is small, force a monotone μ decrease and set the
/// `tiny_step` flag so downstream logic knows to accept the full step.
///
/// Resets the filter on μ change (standard Ipopt convention).
fn detect_tiny_step(
    state: &mut SolverState,
    options: &SolverOptions,
    mu_state: &mut MuState,
    filter: &mut Filter,
    consecutive_tiny_steps: &mut usize,
    alpha_primal_max: f64,
    primal_inf: f64,
) {
    let n = state.n;
    let max_rel_step: f64 = (0..n)
        .map(|i| (alpha_primal_max * state.dx[i]).abs() / (state.x[i].abs() + 1.0))
        .fold(0.0f64, f64::max);
    if max_rel_step < 1e-14 && primal_inf < 1e-4 {
        *consecutive_tiny_steps += 1;
        mu_state.tiny_step = true;
        if *consecutive_tiny_steps >= 2 {
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
            *consecutive_tiny_steps = 0;
        }
    } else {
        *consecutive_tiny_steps = 0;
        mu_state.tiny_step = false;
    }
}

/// Outcome of one factorization attempt with inertia correction.
///
/// `Continue` and `Return` correspond to the main loop's `continue`
/// and early-return paths after a recovery branch. `Proceed` passes
/// the regularization magnitudes `(delta_w, delta_c)` back to the
/// caller for use during iterative refinement.
enum FactorDecision {
    Proceed { ic_delta_w: f64, ic_delta_c: f64 },
    Continue,
    Return(SolveResult),
}

/// Factor the augmented KKT system with inertia correction, handling
/// the recovery cascade on failure.
///
/// No-op when `kkt_system_opt` is `None` (condensed paths do their own
/// factorization downstream). On success returns the inertia-correction
/// `(delta_w, delta_c)` for use in iterative refinement.
///
/// Recovery cascade on factorization failure:
/// 1. Early-iteration (< 5) perturbation sweep across scales 1e-4..1e-1
///    with post-perturbation re-factorization; success → `Continue`.
/// 2. Gradient-descent fallback with Armijo backtracking; success → `Continue`.
/// 3. Restoration phase; success → `Continue`.
/// 4. Late-iteration perturbation sweep (no re-factorization); success → `Continue`.
/// 5. All exhausted → `Return(NumericalError)`.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition.
/// Early-iteration degenerate-starting-point recovery: if KKT factorization
/// failed in the first 5 iterations, perturb x at increasing scales and try
/// to refactor a freshly assembled KKT. Returns true if a perturbation
/// recovered a successful factorization (in which case the caller should
/// `continue` the IPM loop with the perturbed point).
fn try_early_perturbation_recovery<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    iteration: usize,
    n: usize,
    m: usize,
    use_sparse: bool,
    lin_solver: &mut dyn LinearSolver,
    inertia_params: &mut InertiaCorrectionParams,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    filter: &mut Filter,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> bool {
    if iteration >= 5 {
        return false;
    }
    for &perturb_scale in &[1e-4, 1e-3, 1e-2, 5e-2, 1e-1] {
        let x_saved = state.x.clone();
        for i in 0..n {
            let mag = state.x[i].abs().max(1.0);
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
        let pert_eval_ok = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        if pert_eval_ok && !state.obj.is_nan() && !state.obj.is_infinite()
            && !state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
        {
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
                &mut kkt_p, lin_solver, inertia_params,
            ).is_ok() {
                log::debug!(
                    "Early perturbation (scale={:.0e}) recovered factorization at iter {}",
                    perturb_scale, iteration
                );
                filter.reset();
                let theta_p = state.constraint_violation();
                filter.set_theta_min_from_initial(theta_p);
                return true;
            }
        }
        state.x.copy_from_slice(&x_saved);
    }
    false
}

/// Gradient-descent fallback with Armijo backtracking after a KKT
/// factorization failure. Computes a steepest-descent direction (projected
/// to satisfy bound feasibility) and bisects alpha up to 20 times until the
/// objective decreases. On acceptance the state is re-evaluated and the
/// L-BFGS Hessian updated; returns true so the caller can `continue` the
/// IPM loop.
fn try_gradient_descent_fallback<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    n: usize,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> bool {
    let Some(fallback) = gradient_descent_fallback(state) else {
        return false;
    };
    state.dx = fallback.0;
    state.dy = fallback.1;
    state.dz_l = vec![0.0; n];
    state.dz_u = vec![0.0; n];

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
        let mut obj_trial = f64::INFINITY;
        let obj_ok = problem.objective(&x_trial, true, &mut obj_trial);
        if obj_ok && obj_trial.is_finite() && obj_trial < obj_current {
            state.x = x_trial;
            state.obj = obj_trial;
            state.alpha_primal = alpha_fb;
            fb_accepted = true;
            break;
        }
        alpha_fb *= 0.5;
    }
    if fb_accepted {
        let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        return true;
    }
    false
}

fn factor_kkt_with_recovery<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    iteration: usize,
    n: usize,
    m: usize,
    use_sparse: bool,
    kkt_system_opt: &mut Option<kkt::KktSystem>,
    lin_solver: &mut dyn LinearSolver,
    inertia_params: &mut InertiaCorrectionParams,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    filter: &mut Filter,
    restoration: &mut RestorationPhase,
    timings: &mut PhaseTimings,
    prev_ic_delta_w: &mut f64,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
    deadline: Option<Instant>,
) -> FactorDecision {
    let mut ic_delta_w = 0.0f64;
    let mut ic_delta_c = 0.0f64;
    let Some(kkt_system) = kkt_system_opt.as_mut() else {
        return FactorDecision::Proceed { ic_delta_w, ic_delta_c };
    };

    let t_fact = Instant::now();
    if options.print_level >= 5 {
        let dim = match &kkt_system.matrix {
            KktMatrix::Dense(d) => d.n,
            KktMatrix::Sparse(s) => s.n,
        };
        let nnz = match &kkt_system.matrix {
            KktMatrix::Dense(d) => d.n * (d.n + 1) / 2,
            KktMatrix::Sparse(s) => s.triplet_rows.len(),
        };
        rip_log!("ripopt: Factoring KKT dim={} nnz={}...", dim, nnz);
    }
    let inertia_result =
        kkt::factor_with_inertia_correction(kkt_system, lin_solver, inertia_params);
    if options.print_level >= 5 {
        rip_log!("ripopt: KKT factorization took {:.3}s (ok={})",
            t_fact.elapsed().as_secs_f64(), inertia_result.is_ok());
    }
    timings.factorization += t_fact.elapsed();

    if let Ok((dw, dc)) = &inertia_result {
        ic_delta_w = *dw;
        ic_delta_c = *dc;
        *prev_ic_delta_w = *dw;
    }

    if let Some(ref dump_dir) = options.kkt_dump_dir {
        if inertia_result.is_ok() {
            dump_kkt_matrix(
                dump_dir,
                &options.kkt_dump_name,
                iteration,
                kkt_system,
                Some((n, m, 0)),
                ic_delta_w,
                ic_delta_c,
            );
        }
    }

    let Err(e) = inertia_result else {
        return FactorDecision::Proceed { ic_delta_w, ic_delta_c };
    };
    log::warn!("KKT factorization failed: {}", e);

    // Early-iteration perturbation: degenerate starting point recovery.
    if try_early_perturbation_recovery(
        state, problem, iteration, n, m, use_sparse,
        lin_solver, inertia_params, lbfgs_state, filter,
        linear_constraints, lbfgs_mode,
    ) {
        return FactorDecision::Continue;
    }

    // Gradient descent fallback with Armijo backtracking.
    if try_gradient_descent_fallback(
        state, problem, n, lbfgs_state, linear_constraints, lbfgs_mode,
    ) {
        return FactorDecision::Continue;
    }

    // Restoration phase.
    let (x_rest, success) = restoration.restore(
        &state.x, &state.x_l, &state.x_u, &state.g_l, &state.g_u,
        &state.jac_rows, &state.jac_cols, n, m, options,
        &|theta, phi| filter.is_acceptable(theta, phi),
        &|x_eval, g_out| problem.constraints(x_eval, true, g_out),
        &|x_eval, jac_out| problem.jacobian_values(x_eval, true, jac_out),
        Some(&|x_eval: &[f64], obj_out: &mut f64| problem.objective(x_eval, true, obj_out)),
        deadline,
    );
    if success {
        state.x = x_rest;
        state.alpha_primal = 0.0;
        let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        return FactorDecision::Continue;
    }

    // Last resort: perturb x and retry.
    if try_last_resort_perturbation(
        state, problem, iteration, n, lbfgs_state, filter,
        linear_constraints, lbfgs_mode,
    ) {
        return FactorDecision::Continue;
    }
    FactorDecision::Return(make_result(state, SolveStatus::NumericalError))
}

/// Last-resort perturbation after restoration has also failed: cumulatively
/// perturb x at scales 1e-3, 1e-2, 1e-1 and accept the first scale at which
/// problem evaluation produces a finite objective. On success the filter is
/// reset and theta_min recomputed; returns true so the caller can continue.
fn try_last_resort_perturbation<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    iteration: usize,
    n: usize,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    filter: &mut Filter,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> bool {
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
        let pert2_ok = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        if pert2_ok && !state.obj.is_nan() && !state.obj.is_infinite() {
            filter.reset();
            let theta_new = state.constraint_violation();
            filter.set_theta_min_from_initial(theta_new);
            return true;
        }
    }
    false
}

/// Update the barrier parameter μ (interior-point centering parameter) once
/// per iteration after the step has been accepted.
///
/// Handles three cases:
/// - No variable bounds: superlinear decrease `μ^mu_superlinear_decrease_power`.
/// - `MuMode::Free` (adaptive oracle): Loqo mu oracle when enabled, with
///   monotone floor to prevent single-step collapse; falls back to `avg_compl/kappa`
///   when quality-function oracle is off; switches to Fixed mode on repeated
///   insufficient progress or dual-infeasibility stagnation.
/// - `MuMode::Fixed` (monotone decrease): decrease μ only when the barrier
///   subproblem is approximately solved (Ipopt's `IpMonotoneMuUpdate`
///   kappa_eps gate), using the min of linear and superlinear rates; may
///   switch back to Free mode in adaptive strategy.
///
/// On every μ change the filter is reset and `theta_min` is recomputed
/// (Ipopt convention, `IpFilterLSAcceptor.cpp:524-532`).
///
/// Ipopt parallel: `IpMuUpdate::UpdateBarrierParameter`, with concrete
/// subclasses `MonotoneMuUpdate` / `AdaptiveMuUpdate` inline here and the
/// `MuOracle` role played by the embedded Loqo formula.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2). Does not return a value — all effects are via
/// `&mut` on state, mu_state, filter, last_mehrotra_sigma.
/// Loqo barrier-parameter oracle (IpLoqoMuOracle.cpp).
///
/// Centrality-driven μ update: ξ = min(z_i·s_i) / avg_compl, then
/// σ = 0.1 · min(0.05·(1-ξ)/ξ, 2)³, and μ_new = σ · avg_compl,
/// floored from below by Ipopt's monotone schedule
/// `min(κ_μ·μ, μ^sldp)` (so the Loqo-proposed μ can't undershoot
/// a gradual monotone schedule on a single step) and from above
/// by 1e5. The lower clamp `mu_floor` is `mu_min` when the
/// barrier subproblem is approximately solved
/// (`barrier_err ≤ kappa_eps · μ`) and `μ/5` otherwise.
fn compute_loqo_mu(
    state: &SolverState,
    options: &SolverOptions,
    avg_compl: f64,
) -> f64 {
    let n = state.n;
    let barrier_err = compute_barrier_error(state);
    let mu_floor = if barrier_err <= options.barrier_tol_factor * state.mu {
        options.mu_min
    } else {
        (state.mu / 5.0).max(options.mu_min)
    };

    // Compute centrality measure: xi = min(z_i*s_i) / avg_compl
    let mut min_compl = f64::INFINITY;
    for i in 0..n {
        if state.x_l[i].is_finite() {
            let s = (state.x[i] - state.x_l[i]).max(1e-20);
            min_compl = min_compl.min(s * state.z_l[i]);
        }
        if state.x_u[i].is_finite() {
            let s = (state.x_u[i] - state.x[i]).max(1e-20);
            min_compl = min_compl.min(s * state.z_u[i]);
        }
    }
    let xi = if avg_compl > 0.0 && min_compl.is_finite() {
        (min_compl / avg_compl).clamp(0.0, 1.0)
    } else {
        1.0
    };

    let ratio = if xi > 1e-20 {
        (0.05 * (1.0 - xi) / xi).min(2.0)
    } else {
        2.0
    };
    let sigma = 0.1 * ratio.powi(3);
    let loqo_mu = sigma * avg_compl;
    let monotone_floor =
        (options.mu_linear_decrease_factor * state.mu)
            .min(state.mu.powf(options.mu_superlinear_decrease_power));
    let new_mu = loqo_mu
        .max(monotone_floor)
        .clamp(mu_floor, 1e5);

    if options.print_level >= 5 {
        rip_log!("ripopt: mu loqo: xi={:.4} sigma={:.4} avg_compl={:.3e} floor={:.3e} -> mu={:.3e}",
            xi, sigma, avg_compl, monotone_floor, new_mu);
    }
    new_mu
}

fn update_barrier_parameter(
    state: &mut SolverState,
    mu_state: &mut MuState,
    filter: &mut Filter,
    last_mehrotra_sigma: &mut Option<f64>,
    options: &SolverOptions,
) {
    let n = state.n;
    // When there are no variable bounds, mu serves no barrier purpose but is
    // still used for KKT regularization and the filter line search. We decrease
    // mu superlinearly (mu^1.5) rather than collapsing it instantly to mu_min,
    // which would destroy filter protection against infeasible steps. This
    // prevents the PENTAGON-type failure where mu=1e-11 at iteration 1 causes
    // the switching condition to accept a step that destroys feasibility.
    let has_bounds = (0..n).any(|i| state.x_l[i].is_finite() || state.x_u[i].is_finite());
    if !has_bounds {
        state.mu = state.mu.powf(options.mu_superlinear_decrease_power).max(options.mu_min);
        return;
    }

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

    // Track dual infeasibility for stagnation detection.
    // If du is not decreasing over 3 consecutive iterations in Free mode,
    // force switch to Fixed mode with mu = 0.8 * avg_compl.
    {
        let du_now = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        if mu_state.dual_inf_window.len() >= 3 {
            mu_state.dual_inf_window.remove(0);
        }
        mu_state.dual_inf_window.push(du_now);
    }

    // Ipopt-style barrier-subproblem stop test (IpMonotoneMuUpdate.cpp:135-194).
    // Decrease mu only when the current barrier subproblem is approximately
    // solved: barrier_err <= kappa_eps * mu (kappa_eps = barrier_tol_factor,
    // default 10). Without this gate, mu collapses every iteration regardless
    // of whether the line search is actually making progress on the current
    // subproblem — observed on cho parmest where mu went 1e-1 -> 1e-9 in 9
    // iters while inf_pr stayed pinned at 19. The check_sufficient_progress
    // gate below is a relative-history check; this is the absolute gate.
    let barrier_err_for_gate = compute_barrier_error(state);
    let barrier_subproblem_solved =
        barrier_err_for_gate <= options.barrier_tol_factor * state.mu;

    match mu_state.mode {
        MuMode::Free => {
            update_barrier_parameter_free_mode(
                state, mu_state, filter, last_mehrotra_sigma, options,
                sufficient, kkt_error, barrier_subproblem_solved,
            );
        }
        MuMode::Fixed => {
            update_barrier_parameter_fixed_mode(
                state, mu_state, filter, options, sufficient, kkt_error,
            );
        }
    }
}

/// Free-mode (adaptive) barrier-parameter update. Three branches:
/// 1) Sufficient progress + barrier subproblem solved: pick a new mu via
///    the Loqo oracle (when quality_function is on) or rate-limited Loqo
///    fallback `avg_compl/kappa`, then reset the filter.
/// 2) Insufficient progress (>=2 consecutive) or dual-infeasibility
///    stagnation: switch to Fixed mode with mu = adaptive_mu_monotone_init
///    * avg_compl, reset the filter.
/// 3) Stay in Free with conservative mu decrease only when the barrier
///    subproblem is approximately solved; otherwise mu stays put waiting
///    for the line search to make progress on the current subproblem.
#[allow(clippy::too_many_arguments)]
fn update_barrier_parameter_free_mode(
    state: &mut SolverState,
    mu_state: &mut MuState,
    filter: &mut Filter,
    last_mehrotra_sigma: &mut Option<f64>,
    options: &SolverOptions,
    sufficient: bool,
    kkt_error: f64,
    barrier_subproblem_solved: bool,
) {
    // Consume Mehrotra sigma for use as quality function candidate
    let _sigma_mu = last_mehrotra_sigma.take();
    if sufficient && !mu_state.tiny_step && barrier_subproblem_solved {
        mu_state.consecutive_insufficient = 0;
        mu_state.remember_accepted(kkt_error);
        let avg_compl = compute_avg_complementarity(state);
        if options.mu_oracle_quality_function && avg_compl > 0.0 {
            state.mu = compute_loqo_mu(state, options, avg_compl);
        } else if avg_compl > 0.0 {
            // Loqo fallback with rate limit (at most 5x decrease if subproblem not solved)
            let barrier_err = compute_barrier_error(state);
            let mu_floor = if barrier_err <= options.barrier_tol_factor * state.mu {
                options.mu_min
            } else {
                (state.mu / 5.0).max(options.mu_min)
            };
            state.mu = (avg_compl / options.kappa).clamp(mu_floor, 1e5);
        } else {
            // No complementarity products (all bounds inactive) → decrease mu
            state.mu = (options.mu_linear_decrease_factor * state.mu)
                .max(options.mu_min);
        }
        // Reset filter on mu change (Ipopt convention,
        // IpFilterLSAcceptor.cpp:524-532). Now that the barrier-
        // stop gate above ensures mu only changes when the
        // subproblem is approximately solved, this happens at
        // the same low frequency as in Ipopt — not every iter.
        filter.reset();
        let theta_new = state.constraint_violation();
        filter.set_theta_min_from_initial(theta_new);
    } else {
        // Also check dual infeasibility stagnation: if du is not improving
        // over 3 consecutive iterations, force the switch to Fixed mode
        // even if consecutive_insufficient < 2.
        let du_stagnant = if mu_state.dual_inf_window.len() >= 3 {
            let w = &mu_state.dual_inf_window;
            let recent = w[w.len()-1];
            let oldest = w[w.len()-3];
            // Stagnant if recent du is >=90% of oldest (not improving)
            // and du is still large relative to tolerance
            recent >= 0.9 * oldest && recent > options.tol * 100.0
        } else {
            false
        };
        mu_state.consecutive_insufficient += 1;
        if mu_state.consecutive_insufficient >= 2 || du_stagnant {
            // Switch to fixed mode after 2 consecutive insufficient iterations
            // or when dual infeasibility is stagnant
            mu_state.consecutive_insufficient = 0;
            log::debug!("Switching to fixed mu mode (insufficient progress or tiny step)");
            state.diagnostics.mu_mode_switches += 1;
            mu_state.mode = MuMode::Fixed;
            mu_state.first_iter_in_mode = true;
            let avg_compl = compute_avg_complementarity(state);
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
        } else if barrier_subproblem_solved {
            // Stay in Free mode with conservative mu decrease —
            // only when the barrier subproblem is approximately
            // solved. Without this gate, mu collapses unconditionally
            // even when the line search is making no progress
            // (observed on cho parmest: mu 0.1 -> 0.02 at iter 1
            // despite barrier_err=1.4e4).
            let avg_compl = compute_avg_complementarity(state);
            if avg_compl > 0.0 {
                let mu_floor = options.mu_min;
                state.mu = (avg_compl / options.kappa).clamp(mu_floor, 1e5);
            } else {
                state.mu = (options.mu_linear_decrease_factor * state.mu)
                    .max(options.mu_min);
            }
        }
        // When !barrier_subproblem_solved, mu stays put and we
        // wait for the line search to make progress on the current
        // subproblem.
    }
}

/// Fixed-mode (monotone) barrier-parameter update. Either switches back to
/// Free mode when adaptive strategy + sufficient progress is detected
/// (skipping the first iteration in Fixed mode), or — when the barrier
/// subproblem is approximately solved (barrier_err <= kappa_eps*mu) or a
/// tiny step was taken — decreases mu by the min of linear and superlinear
/// rates. Filter and theta_min are reset on every accepted decrease.
/// Mirrors Ipopt's IpMonotoneMuUpdate.cpp.
fn update_barrier_parameter_fixed_mode(
    state: &mut SolverState,
    mu_state: &mut MuState,
    filter: &mut Filter,
    options: &SolverOptions,
    sufficient: bool,
    kkt_error: f64,
) {
    if options.mu_strategy_adaptive && sufficient && !mu_state.tiny_step && !mu_state.first_iter_in_mode {
        // Switch back to free mode (only in adaptive strategy)
        log::debug!("Switching back to free mu mode (sufficient progress)");
        state.diagnostics.mu_mode_switches += 1;
        mu_state.mode = MuMode::Free;
        mu_state.remember_accepted(kkt_error);
        mu_state.first_iter_in_mode = true;
    } else {
        mu_state.first_iter_in_mode = false;
        // Check if subproblem is solved (barrier error small enough)
        let barrier_err = compute_barrier_error(state);
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

/// State threaded through the convergence-check promotion attempts.
///
/// Groups the loop-spanning history + one-shot promotion flags so the
/// `check_convergence_and_handle_promotions` helper doesn't need a 20-param
/// signature. All fields live in `solve_ipm`'s stack frame; this struct is
/// only a bundle of mutable references.
struct ConvergenceWorkspace<'a> {
    du_history: &'a mut Vec<f64>,
    iterate_history: &'a mut Vec<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)>,
    tried_iterate_averaging: &'a mut bool,
    tried_active_set: &'a mut bool,
    tried_compl_polish: &'a mut bool,
}

/// Check convergence status and, on Acceptable, try three promotion
/// strategies in order: iterate averaging (oscillation smoothing), active-set
/// reduced solve, and complementarity polishing via multiplier snap. Each
/// promotion attempt is one-shot (guarded by its `tried_*` flag). Any
/// successful promotion returns `Some(SolveResult { Optimal })`.
///
/// Returns:
/// - `Some(SolveResult)` when the solver should terminate (Converged,
///   Acceptable, promoted Optimal, or Diverging/Unbounded).
/// - `None` when iteration should continue (NotConverged).
///
/// Ipopt parallel: `IpIpoptAlg::Optimize` → `IpConvCheck::CheckConvergence`
/// plus the near-tolerance polishing heuristics (`RecalcIpoptData`,
/// multiplier snap), which live inline in Ipopt as well.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2).
#[allow(clippy::too_many_arguments)]
/// Strategy 1 (Acceptable promotion): if the dual-infeasibility history has
/// `AVG_WINDOW` entries and shows oscillation (>= AVG_WINDOW/2 sign changes
/// in consecutive differences), average the last `AVG_WINDOW` iterates and
/// re-check convergence at the averaged point. On success returns
/// `Some(Optimal)`; otherwise restores the original state and returns None.
/// One-shot (guarded by `ws.tried_iterate_averaging`).
#[allow(clippy::too_many_arguments)]
fn try_iterate_averaging_promotion<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    ws: &mut ConvergenceWorkspace,
    avg_window: usize,
    n: usize,
    m: usize,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> Option<SolveResult> {
    if *ws.tried_iterate_averaging || ws.du_history.len() != avg_window {
        return None;
    }
    let mut sign_changes = 0;
    for w in 1..ws.du_history.len() - 1 {
        let d1 = ws.du_history[w] - ws.du_history[w - 1];
        let d2 = ws.du_history[w + 1] - ws.du_history[w];
        if d1 * d2 < 0.0 {
            sign_changes += 1;
        }
    }
    if sign_changes < avg_window / 2 {
        return None;
    }
    *ws.tried_iterate_averaging = true;
    let len = ws.iterate_history.len() as f64;
    let mut avg_x = vec![0.0; n];
    let mut avg_y = vec![0.0; m];
    let mut avg_zl = vec![0.0; n];
    let mut avg_zu = vec![0.0; n];
    for (hx, hy, hzl, hzu) in ws.iterate_history.iter() {
        for i in 0..n { avg_x[i] += hx[i] / len; }
        for i in 0..m { avg_y[i] += hy[i] / len; }
        for i in 0..n { avg_zl[i] += hzl[i] / len; }
        for i in 0..n { avg_zu[i] += hzu[i] / len; }
    }
    for i in 0..n {
        avg_x[i] = avg_x[i].clamp(
            if state.x_l[i].is_finite() { state.x_l[i] + 1e-15 } else { f64::NEG_INFINITY },
            if state.x_u[i].is_finite() { state.x_u[i] - 1e-15 } else { f64::INFINITY },
        );
        avg_zl[i] = avg_zl[i].max(0.0);
        avg_zu[i] = avg_zu[i].max(0.0);
    }
    let saved_x = state.x.clone();
    let saved_y = state.y.clone();
    let saved_zl = state.z_l.clone();
    let saved_zu = state.z_u.clone();
    state.x.copy_from_slice(&avg_x);
    state.y.copy_from_slice(&avg_y);
    state.z_l.copy_from_slice(&avg_zl);
    state.z_u.copy_from_slice(&avg_zu);
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    let avg_pr = convergence::primal_infeasibility_max(&state.g, &state.g_l, &state.g_u);
    let avg_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &avg_y, &avg_zl, &avg_zu, n,
    );
    let avg_compl = convergence::complementarity_error(
        &avg_x, &state.x_l, &state.x_u, &avg_zl, &avg_zu, 0.0,
    );
    let avg_conv = ConvergenceInfo {
        primal_inf: avg_pr, dual_inf: avg_du,
        dual_inf_unscaled: avg_du,
        compl_inf: avg_compl,
        mu: state.mu, objective: state.obj,
        multiplier_sum: avg_y.iter().map(|v| v.abs()).sum::<f64>()
            + avg_zl.iter().map(|v| v.abs()).sum::<f64>()
            + avg_zu.iter().map(|v| v.abs()).sum::<f64>(),
        multiplier_count: m + 2 * n,
    };
    if let ConvergenceStatus::Converged = check_convergence(&avg_conv, options, 0) {
        if options.print_level >= 3 {
            rip_log!("ripopt: Iterate averaging promoted near-tolerance -> Optimal (du={:.2e})", avg_du);
        }
        return Some(make_result(state, SolveStatus::Optimal));
    }
    state.x.copy_from_slice(&saved_x);
    state.y.copy_from_slice(&saved_y);
    state.z_l.copy_from_slice(&saved_zl);
    state.z_u.copy_from_slice(&saved_zu);
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    None
}

/// Strategy 4 (Acceptable promotion): when complementarity is the
/// bottleneck (primal_inf and dual_inf already within 100x tol but
/// compl_inf > tol*s_d), snap bound multipliers to reduce
/// complementarity. For each variable that is clearly interior to a
/// bound (gap > 1e-6), zero out the corresponding z; keep z otherwise.
/// Re-check convergence with the snapped multipliers; on success return
/// Optimal, otherwise restore z and return None. One-shot via
/// `ws.tried_compl_polish`.
fn try_complementarity_polish_promotion(
    state: &mut SolverState,
    options: &SolverOptions,
    conv_info: &ConvergenceInfo,
    ws: &mut ConvergenceWorkspace,
    n: usize,
    m: usize,
) -> Option<SolveResult> {
    if *ws.tried_compl_polish {
        return None;
    }
    let compl_inf_now = conv_info.compl_inf;
    let s_d_now = {
        let s_max: f64 = 100.0;
        let s_d_max: f64 = 1e4;
        if conv_info.multiplier_count > 0 {
            ((s_max.max(conv_info.multiplier_sum / conv_info.multiplier_count as f64)) / s_max).min(s_d_max)
        } else { 1.0 }
    };
    let compl_tol_scaled = options.tol * s_d_now;
    if !(compl_inf_now > compl_tol_scaled
        && conv_info.primal_inf <= 100.0 * options.tol
        && conv_info.dual_inf <= 100.0 * options.tol * s_d_now)
    {
        return None;
    }
    *ws.tried_compl_polish = true;
    let saved_zl = state.z_l.clone();
    let saved_zu = state.z_u.clone();
    let gap_tol = 1e-6;
    for i in 0..n {
        let gap_l = if state.x_l[i].is_finite() { state.x[i] - state.x_l[i] } else { f64::INFINITY };
        let gap_u = if state.x_u[i].is_finite() { state.x_u[i] - state.x[i] } else { f64::INFINITY };
        if gap_l > gap_tol {
            state.z_l[i] = 0.0;
        }
        if gap_u > gap_tol {
            state.z_u[i] = 0.0;
        }
    }
    let snap_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let snap_compl = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );
    let snap_mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
    let snap_conv = ConvergenceInfo {
        primal_inf: conv_info.primal_inf,
        dual_inf: snap_du,
        dual_inf_unscaled: snap_du,
        compl_inf: snap_compl,
        mu: state.mu,
        objective: state.obj,
        multiplier_sum: snap_mult_sum,
        multiplier_count: m + 2 * n,
    };
    if let ConvergenceStatus::Converged = check_convergence(&snap_conv, options, 0) {
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Complementarity snap promoted near-tolerance -> Optimal (compl {:.2e} -> {:.2e}, du {:.2e})",
                compl_inf_now, snap_compl, snap_du
            );
        }
        return Some(make_result(state, SolveStatus::Optimal));
    }
    state.z_l.copy_from_slice(&saved_zl);
    state.z_u.copy_from_slice(&saved_zu);
    None
}

fn check_convergence_and_handle_promotions<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    primal_inf_max: f64,
    dual_inf: f64,
    dual_inf_unscaled: f64,
    compl_inf: f64,
    multiplier_sum: f64,
    multiplier_count: usize,
    ws: &mut ConvergenceWorkspace,
    timings: &PhaseTimings,
    iteration: usize,
    ipm_start: Instant,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> Option<SolveResult> {
    const AVG_WINDOW: usize = 6;
    let n = state.n;
    let m = state.m;

    let conv_info = ConvergenceInfo {
        primal_inf: primal_inf_max,
        dual_inf,
        dual_inf_unscaled,
        compl_inf,
        mu: state.mu,
        objective: state.obj,
        multiplier_sum,
        multiplier_count,
    };

    // Track iterate history for oscillation detection (Strategy 1)
    ws.du_history.push(dual_inf);
    ws.iterate_history.push((
        state.x.clone(), state.y.clone(), state.z_l.clone(), state.z_u.clone(),
    ));
    if ws.du_history.len() > AVG_WINDOW {
        ws.du_history.remove(0);
        ws.iterate_history.remove(0);
    }

    match check_convergence(&conv_info, options, state.consecutive_acceptable) {
        ConvergenceStatus::Converged => {
            if options.print_level >= 5 {
                timings.print_summary(iteration + 1, ipm_start.elapsed());
            }
            Some(make_result(state, SolveStatus::Optimal))
        }
        ConvergenceStatus::Acceptable => {
            // Strategy 1: Try iterate averaging before declaring Acceptable
            if let Some(result) = try_iterate_averaging_promotion(
                state, problem, options, ws, AVG_WINDOW, n, m,
                linear_constraints, lbfgs_mode,
            ) {
                return Some(result);
            }

            // Strategy 3: Try active set identification + reduced solve
            if !*ws.tried_active_set {
                *ws.tried_active_set = true;
                if let Some(result) = try_active_set_solve(state, problem, options, linear_constraints, lbfgs_mode) {
                    if options.print_level >= 3 {
                        rip_log!("ripopt: Active set solve promoted Acceptable -> Optimal");
                    }
                    return Some(result);
                }
            }

            // Strategy 4: Complementarity polishing via multiplier snap
            if let Some(result) = try_complementarity_polish_promotion(
                state, options, &conv_info, ws, n, m,
            ) {
                return Some(result);
            }

            if options.print_level >= 5 {
                timings.print_summary(iteration + 1, ipm_start.elapsed());
            }
            // Promoted to SolveStatus::Acceptable (matches Ipopt's
            // Solved_To_Acceptable_Level). Previously this fell through
            // to NumericalError, which caused problems that met Ipopt's
            // default acceptable-level tolerances to show as unsolved in
            // benchmarks. Benchmark reporter counts Acceptable as solved.
            Some(make_result(state, SolveStatus::Acceptable))
        }
        ConvergenceStatus::Diverging => {
            Some(make_result(state, SolveStatus::Unbounded))
        }
        ConvergenceStatus::NotConverged => None,
    }
}

/// Check wall-clock and early-stall time limits at the top of each iteration.
///
/// Returns `Some(SolveResult)` to terminate the loop if either:
/// - `max_wall_time` has been exceeded → `MaxIterations`
/// - `early_stall_timeout` was hit during the first 5 iterations
///   (scaled by problem size) → `NumericalError`
///
/// Wall-clock is polled every iteration during the first 10, then every 10
/// thereafter to keep overhead negligible on long runs.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2). Pure guard function.
fn check_time_limits(
    state: &SolverState,
    iteration: usize,
    start_time: Instant,
    early_timeout: f64,
    options: &SolverOptions,
) -> Option<SolveResult> {
    // Check wall-clock time limit (every iteration in early phase, every 10 after)
    if (iteration < 10 || iteration % 10 == 0) && options.max_wall_time > 0.0 {
        if start_time.elapsed().as_secs_f64() >= options.max_wall_time {
            return Some(make_result(state, SolveStatus::MaxIterations));
        }
    }

    // Early stall detection: bail out if stuck in early iterations.
    // `early_timeout` is pre-scaled by problem size in the caller (see
    // scaling formula in solve_ipm).
    if iteration < 5 && options.early_stall_timeout > 0.0 {
        if start_time.elapsed().as_secs_f64() > early_timeout {
            if options.print_level >= 3 {
                rip_log!(
                    "ripopt: Early stall at iteration {} ({:.1}s elapsed), terminating",
                    iteration, start_time.elapsed().as_secs_f64()
                );
            }
            return Some(make_result(state, SolveStatus::NumericalError));
        }
    }
    None
}

/// Per-iteration KKT residuals used by the log row, filter, and
/// convergence check.
struct OptimalityMeasures {
    /// 1-norm constraint violation (filter/log).
    primal_inf: f64,
    /// Max-norm primal infeasibility (convergence gate).
    primal_inf_max: f64,
    /// Iterative-z dual infeasibility `||∇f + J^T y - z_L + z_U||_∞`.
    dual_inf: f64,
    /// Component-wise scaled dual infeasibility for the unscaled gate
    /// (divides each component by `1 + |∇f_i|`).
    dual_inf_unscaled: f64,
    /// Complementarity error `max_i {(x-x_L)_i z_{L,i}, (x_U-x)_i z_{U,i}}`
    /// against target 0 (unscaled).
    compl_inf: f64,
}

/// Compute the per-iteration optimality residuals.
///
/// Ipopt parallel: `IpIpoptCalculatedQuantities::curr_dual_infeasibility`
/// / `curr_primal_infeasibility` / `curr_complementarity`.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition. Pure function of `state`.
fn compute_optimality_measures(state: &SolverState) -> OptimalityMeasures {
    let n = state.n;
    let primal_inf = state.constraint_violation();
    let primal_inf_max = convergence::primal_infeasibility_max(&state.g, &state.g_l, &state.g_u);

    // Iterative-z dual infeasibility (matches Ipopt's curr_dual_infeasibility):
    // honest KKT residual. If iterative z is inconsistent with ∇f + J^T y the
    // residual stays large and iteration continues.
    let dual_inf = convergence::dual_infeasibility(
        &state.grad_f,
        &state.jac_rows,
        &state.jac_cols,
        &state.jac_vals,
        &state.y,
        &state.z_l,
        &state.z_u,
        n,
    );
    let dual_inf_unscaled = convergence::dual_infeasibility_scaled(
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
    OptimalityMeasures {
        primal_inf,
        primal_inf_max,
        dual_inf,
        dual_inf_unscaled,
        compl_inf,
    }
}

/// Compute the multiplier-based scaling factor `s_d` used by the
/// Ipopt acceptable-tolerance gate and the overall-progress stall
/// detector, then update `state.consecutive_acceptable`.
///
/// Scaling matches `IpIpoptCalculatedQuantities::ComputeOptimalityErrorScaling`
/// with `s_max=100` and the 1e4 cap preserved for compatibility with
/// the rest of ripopt's tolerance pipeline. Acceptable thresholds
/// match `IpOptErrorConvCheck.cpp:70-121` defaults
/// (acceptable_tol=1e-6, acceptable_constr_viol_tol=1e-2,
/// acceptable_dual_inf_tol=1e10, acceptable_compl_inf_tol=1e-2).
///
/// Returns `s_d_for_acc` because it is consumed downstream by
/// `detect_and_handle_progress_stall`.
fn track_consecutive_acceptable(
    state: &mut SolverState,
    primal_inf: f64,
    dual_inf: f64,
    dual_inf_unscaled: f64,
    compl_inf: f64,
    multiplier_sum: f64,
) -> f64 {
    let n = state.n;
    let m = state.m;
    let s_d_for_acc = {
        let s_max: f64 = 100.0;
        let s_d_max: f64 = 1e4;
        if (m + 2 * n) > 0 {
            ((s_max.max(multiplier_sum / (m + 2 * n) as f64)) / s_max).min(s_d_max)
        } else {
            1.0
        }
    };
    let meets_acc_scaled = primal_inf <= 1e-6
        && dual_inf <= 1e-6 * s_d_for_acc
        && compl_inf <= 1e-6 * s_d_for_acc;
    let meets_acc_unscaled = primal_inf <= 1e-2
        && dual_inf_unscaled <= 1e10
        && compl_inf <= 1e-2;
    if meets_acc_scaled && meets_acc_unscaled {
        state.consecutive_acceptable += 1;
    } else {
        state.consecutive_acceptable = 0;
    }
    s_d_for_acc
}

/// Record the current iterate as the best-du point if its dual
/// infeasibility beats the previous best.
///
/// Used by the overall-progress stall detector to revert before a
/// `NumericalError` exit, and by the dual-stagnation detector to
/// restart from the good point.
fn update_best_du_iterate(
    state: &SolverState,
    dual_inf: f64,
    best_du_val: &mut f64,
    best_du_x: &mut Option<Vec<f64>>,
    best_du_y: &mut Option<Vec<f64>>,
    best_du_zl: &mut Option<Vec<f64>>,
    best_du_zu: &mut Option<Vec<f64>>,
) {
    if dual_inf < *best_du_val {
        *best_du_val = dual_inf;
        *best_du_x = Some(state.x.clone());
        *best_du_y = Some(state.y.clone());
        *best_du_zl = Some(state.z_l.clone());
        *best_du_zu = Some(state.z_u.clone());
    }
}

/// Detect unboundedness by requiring 10 consecutive iterations of
/// `obj < -1e20` at a feasible iterate. The counter prevents false
/// positives from transient dips.
fn detect_unbounded(
    state: &SolverState,
    options: &SolverOptions,
    primal_inf: f64,
    consecutive_unbounded: &mut usize,
) -> Option<SolveResult> {
    if state.obj < -1e20 && primal_inf < options.constr_viol_tol {
        *consecutive_unbounded += 1;
        if *consecutive_unbounded >= 10 {
            return Some(make_result(state, SolveStatus::Unbounded));
        }
    } else {
        *consecutive_unbounded = 0;
    }
    None
}

/// Primal-infeasibility divergence detector.
///
/// Tracks consecutive iterations where θ is strictly increasing (by
/// at least 1e-6 relative). After 8 such iterations *and* a cumulative
/// 20%+ growth from the start of the divergence run, sets the
/// force-restoration flag so the main loop re-enters restoration
/// rather than continuing with worsening feasibility.
///
/// Catches the pattern where after NLP restoration the filter accepts
/// steps via a slight φ decrease while θ grows steadily. Ipopt has no
/// direct analogue — this is a ripopt safety valve.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition. Returns `true` when the caller should force
/// restoration on the next step.
fn detect_primal_divergence(
    options: &SolverOptions,
    iteration: usize,
    primal_inf: f64,
    consecutive_pr_increase: &mut usize,
    pr_prev_for_divergence: &mut f64,
    pr_at_divergence_start: &mut f64,
    m: usize,
) -> bool {
    let mut force_restoration = false;
    if m > 0 && iteration > 5 && primal_inf > options.constr_viol_tol {
        if primal_inf > *pr_prev_for_divergence * (1.0 + 1e-6) {
            if *consecutive_pr_increase == 0 {
                *pr_at_divergence_start = *pr_prev_for_divergence;
            }
            *consecutive_pr_increase += 1;
        } else {
            *consecutive_pr_increase = 0;
        }
        // After 8 consecutive increases AND cumulative growth of at least
        // 20%, force restoration. The growth check prevents triggering on
        // tiny numerical oscillations.
        if *consecutive_pr_increase >= 8 && primal_inf > 1.2 * *pr_at_divergence_start {
            log::info!(
                "Primal divergence at iter {}: pr grew for {} consecutive iterations ({:.2e} -> {:.2e}), forcing restoration",
                iteration, *consecutive_pr_increase, *pr_at_divergence_start, primal_inf
            );
            if options.print_level >= 3 {
                rip_log!(
                    "ripopt: Primal divergence detected (pr grew {:.2e} -> {:.2e} over {} iters), re-entering restoration",
                    *pr_at_divergence_start, primal_inf, *consecutive_pr_increase
                );
            }
            force_restoration = true;
            *consecutive_pr_increase = 0;
        }
    } else {
        *consecutive_pr_increase = 0;
    }
    *pr_prev_for_divergence = primal_inf;
    force_restoration
}

/// Outcome of the overall-progress stall detector.
enum StallDecision {
    /// No stall (or not yet activated) — fall through to normal step.
    Proceed,
    /// Stall detected and recovered by bumping μ / resetting the filter —
    /// caller must `continue` the main loop immediately (no Newton step
    /// this iteration).
    Continue,
    /// Stall detected and no recovery possible — caller must return this
    /// result immediately.
    Return(SolveResult),
}

/// Overall-progress stall detector: terminate or recover when neither
/// primal nor dual infeasibility has improved for many iterations.
///
/// Two triggers (see block comment in body): (1) tiny steps for half
/// the `stall_iter_limit`, (2) no metric improvement for the full
/// limit. Before declaring `NumericalError`:
/// - If the current point is near-tolerance, attempt a μ boost or a
///   forced μ decrease (Fixed mode) and restart the filter.
/// - Re-evaluate with "optimal duals" recomputed from the current
///   gradient (covers cases where iterative y/z diverged even though
///   the primal is near-optimal, e.g. HS116, CONCON).
/// - Optionally revert to the best-du iterate seen so far.
///
/// No direct Ipopt parallel; this is a ripopt-specific safety net for
/// long-running stalls (see CLAUDE.md: CONCON, HS13, HS116).
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition.
#[allow(clippy::too_many_arguments)]
/// Full two-gate near-tolerance check with optimal dual multipliers. When
/// duals have diverged but the primal point is near-optimal (HS116,
/// CONCON), the simple `compl_inf`/`dual_inf` check fails because it uses
/// the current (diverged) duals. This helper recomputes optimal duals from
/// the gradient (z = max(0, ∇L) capped at kc*mu/slack), then checks both
/// the scaled (sc) and unscaled (usc) tolerance gates. Returns
/// `Some(StallDecision::Return(NumericalError))` when both gates pass —
/// in which case the caller should exit with NumericalError because the
/// iterate is good enough to not be a hard failure but the solver can't
/// drive it further. Returns `None` to let the caller continue with the
/// remaining stall-recovery logic.
fn check_stall_near_tolerance_via_optimal_duals(
    state: &SolverState,
    options: &SolverOptions,
    primal_inf: f64,
    primal_inf_max: f64,
    compl_inf: f64,
    stall_near_tol: f64,
    n: usize,
    m: usize,
) -> Option<StallDecision> {
    let mut gj = state.grad_f.clone();
    for (idx, (&row, &col)) in
        state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate()
    {
        gj[col] += state.jac_vals[idx] * state.y[row];
    }
    let mut opt_zl = vec![0.0; n];
    let mut opt_zu = vec![0.0; n];
    let kc = 1e10;
    for i in 0..n {
        if gj[i] > 0.0 && state.x_l[i].is_finite() {
            let sl = (state.x[i] - state.x_l[i]).max(1e-20);
            if gj[i] * sl <= kc * state.mu.max(1e-20) {
                opt_zl[i] = gj[i];
            }
        } else if gj[i] < 0.0 && state.x_u[i].is_finite() {
            let su = (state.x_u[i] - state.x[i]).max(1e-20);
            if (-gj[i]) * su <= kc * state.mu.max(1e-20) {
                opt_zu[i] = -gj[i];
            }
        }
    }
    let opt_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &opt_zl, &opt_zu, n,
    );
    let opt_co = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &opt_zl, &opt_zu, 0.0,
    );
    let opt_co_best = compl_inf.min(opt_co);
    let fmult: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
        + opt_zl.iter().map(|v| v.abs()).sum::<f64>()
        + opt_zu.iter().map(|v| v.abs()).sum::<f64>();
    let fsd = if (m + 2 * n) > 0 {
        ((100.0f64.max(fmult / (m + 2 * n) as f64)) / 100.0).min(1e4)
    } else {
        1.0
    };
    let stall_fdu_tol = (stall_near_tol * fsd).max(1e-2);
    let stall_fco_tol = (stall_near_tol * fsd).max(1e-2);
    let stall_fpr_tol = stall_near_tol.max(10.0 * options.constr_viol_tol);
    let sc = primal_inf_max <= stall_fpr_tol
        && opt_du <= stall_fdu_tol
        && opt_co_best <= stall_fco_tol;
    let du_u = convergence::dual_infeasibility_scaled(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let usc = primal_inf_max <= 10.0 * options.constr_viol_tol
        && du_u <= 10.0 * options.dual_inf_tol
        && opt_co_best <= 10.0 * options.compl_inf_tol;
    if sc && usc {
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Stalled but near-tolerance via optimal duals (pr={:.2e}, du_opt={:.2e}, co={:.2e}), returning NumericalError",
                primal_inf, opt_du, opt_co_best
            );
        }
        return Some(StallDecision::Return(make_result(state, SolveStatus::NumericalError)));
    }
    None
}

/// Handle a stall whose current iterate is near-tolerance (1000x tol on
/// each metric). Three sub-decisions in priority order:
/// 1. mu has outrun feasibility (mu < 0.01*pr_max while pr_max above
///    constr_viol_tol): boost mu to 0.1*pr_max, switch to Fixed mode,
///    return Continue.
/// 2. Fixed (monotone) mode + mu still above mu_min: force a mu decrease
///    by min(linear, superlinear) rate (the barrier subproblem is
///    effectively solved at the current mu), return Continue.
/// 3. Otherwise, classify as Acceptable when both unscaled (1e-2) and
///    scaled (1e-6 * s_d) tolerances pass; else return NumericalError.
#[allow(clippy::too_many_arguments)]
fn handle_near_tolerance_stall(
    state: &mut SolverState,
    options: &SolverOptions,
    primal_inf: f64,
    primal_inf_max: f64,
    dual_inf: f64,
    compl_inf: f64,
    s_d_for_acc: f64,
    filter: &mut Filter,
    mu_state: &mut MuState,
    stall_no_progress_count: &mut usize,
    stall_best_pr: &mut f64,
    stall_best_du: &mut f64,
) -> StallDecision {
    if state.mu < primal_inf_max * 0.01 && primal_inf_max > options.constr_viol_tol {
        let new_mu = (primal_inf_max * 0.1).max(1e-6);
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Near-tolerance stall: boosting mu {:.2e} -> {:.2e} (pr_max={:.2e})",
                state.mu, new_mu, primal_inf_max
            );
        }
        state.mu = new_mu;
        filter.reset();
        let theta_new = state.constraint_violation();
        filter.set_theta_min_from_initial(theta_new);
        *stall_no_progress_count = 0;
        *stall_best_pr = f64::INFINITY;
        *stall_best_du = f64::INFINITY;
        mu_state.mode = MuMode::Fixed;
        mu_state.first_iter_in_mode = true;
        return StallDecision::Continue;
    }
    if !options.mu_strategy_adaptive && state.mu > options.mu_min {
        let new_mu = (options.mu_linear_decrease_factor * state.mu)
            .min(state.mu.powf(options.mu_superlinear_decrease_power))
            .max(options.mu_min);
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Fixed mode stall near tolerance (pr={:.2e}, du={:.2e}, co={:.2e}), forcing mu {:.2e} -> {:.2e}",
                primal_inf, dual_inf, compl_inf, state.mu, new_mu
            );
        }
        state.mu = new_mu;
        filter.reset();
        let theta_new = state.constraint_violation();
        filter.set_theta_min_from_initial(theta_new);
        *stall_no_progress_count = 0;
        *stall_best_pr = f64::INFINITY;
        *stall_best_du = f64::INFINITY;
        return StallDecision::Continue;
    }
    let acc_pr_ok = primal_inf_max <= 1e-2;
    let acc_du_ok = dual_inf <= 1e10;
    let acc_co_ok = compl_inf <= 1e-2;
    let acc_scaled_ok = primal_inf_max <= 1e-6
        && dual_inf <= 1e-6 * s_d_for_acc
        && compl_inf <= 1e-6 * s_d_for_acc;
    if acc_pr_ok && acc_du_ok && acc_co_ok && acc_scaled_ok {
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Stalled at acceptable tolerance (pr={:.2e}, du={:.2e}, co={:.2e}), returning Acceptable",
                primal_inf, dual_inf, compl_inf
            );
        }
        return StallDecision::Return(make_result(state, SolveStatus::Acceptable));
    }
    if options.print_level >= 3 {
        rip_log!(
            "ripopt: Stalled but near-tolerance (pr={:.2e}, du={:.2e}, co={:.2e}), returning NumericalError",
            primal_inf, dual_inf, compl_inf
        );
    }
    StallDecision::Return(make_result(state, SolveStatus::NumericalError))
}

fn detect_and_handle_progress_stall<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    iteration: usize,
    primal_inf: f64,
    primal_inf_max: f64,
    dual_inf: f64,
    compl_inf: f64,
    s_d_for_acc: f64,
    filter: &mut Filter,
    mu_state: &mut MuState,
    stall_no_progress_count: &mut usize,
    stall_best_pr: &mut f64,
    stall_best_du: &mut f64,
    best_du_val: f64,
    best_du_x: &Option<Vec<f64>>,
    best_du_y: &Option<Vec<f64>>,
    best_du_zl: &Option<Vec<f64>>,
    best_du_zu: &Option<Vec<f64>>,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> StallDecision {
    if iteration <= 50 || options.stall_iter_limit == 0 {
        return StallDecision::Proceed;
    }
    let n = state.n;
    let m = state.m;

    let pr_improved = primal_inf_max < 0.99 * *stall_best_pr;
    let du_improved = dual_inf < 0.99 * *stall_best_du;
    if pr_improved {
        *stall_best_pr = primal_inf_max;
    }
    if du_improved {
        *stall_best_du = dual_inf;
    }
    if pr_improved || du_improved {
        *stall_no_progress_count = 0;
        return StallDecision::Proceed;
    }
    *stall_no_progress_count += 1;
    let tiny_alpha = state.alpha_primal < 1e-8 && state.alpha_dual < 1e-4;
    // Terminate after half the limit with truly negligible steps, or the full
    // limit with no metric improvement regardless of step size.
    let stall_limit = if tiny_alpha { options.stall_iter_limit / 2 } else { options.stall_iter_limit };
    if *stall_no_progress_count < stall_limit {
        return StallDecision::Proceed;
    }

    // Before declaring NumericalError, check if the current point is
    // near-tolerance (1000x tol).
    let stall_near_tol = options.tol * 1000.0;
    let stall_pr_ok = primal_inf_max <= stall_near_tol.max(10.0 * options.constr_viol_tol);
    let stall_du_ok = dual_inf <= (stall_near_tol * s_d_for_acc).max(1e-2);
    let stall_co_ok = compl_inf <= (stall_near_tol * s_d_for_acc).max(1e-2);
    if stall_pr_ok && stall_du_ok && stall_co_ok {
        return handle_near_tolerance_stall(
            state, options, primal_inf, primal_inf_max, dual_inf, compl_inf,
            s_d_for_acc, filter, mu_state,
            stall_no_progress_count, stall_best_pr, stall_best_du,
        );
    }
    if let Some(decision) = check_stall_near_tolerance_via_optimal_duals(
        state, options, primal_inf, primal_inf_max, compl_inf, stall_near_tol, n, m,
    ) {
        return decision;
    }
    // Before terminating: if primal is close but μ is very small, boost μ
    // to restore KKT regularization for dual convergence. Addresses the
    // "μ outrunning feasibility" pattern.
    if primal_inf_max < 0.1 && state.mu < primal_inf_max * 0.01 {
        let new_mu = (primal_inf_max * 0.1).max(1e-6);
        if options.print_level >= 3 {
            rip_log!(
                "ripopt: Stall recovery: boosting mu {:.2e} -> {:.2e} (pr_max={:.2e})",
                state.mu, new_mu, primal_inf_max
            );
        }
        state.mu = new_mu;
        filter.reset();
        let theta_new = state.constraint_violation();
        filter.set_theta_min_from_initial(theta_new);
        *stall_no_progress_count = 0;
        *stall_best_pr = f64::INFINITY;
        *stall_best_du = f64::INFINITY;
        mu_state.mode = MuMode::Fixed;
        mu_state.first_iter_in_mode = true;
        return StallDecision::Continue;
    }
    // Revert to the best-du iterate if it has significantly better dual
    // feasibility before declaring NumericalError. Post-stall y/z can be
    // corrupted by inertia-escalated KKT solves (CONCON: iter 48 had
    // du=7e-16, stall at iter 81 has du=1.03 at the same x).
    if let (Some(bdx), Some(bdy), Some(bdzl), Some(bdzu)) =
        (best_du_x, best_du_y, best_du_zl, best_du_zu)
    {
        if best_du_val < dual_inf * 0.1 {
            state.x.copy_from_slice(bdx);
            state.y.copy_from_slice(bdy);
            state.z_l.copy_from_slice(bdzl);
            state.z_u.copy_from_slice(bdzu);
            let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
            if options.print_level >= 3 {
                rip_log!(
                    "ripopt: Reverting to best-du iterate (du: {:.2e} -> {:.2e}) before NumericalError exit",
                    dual_inf, best_du_val
                );
            }
        }
    }
    if options.print_level >= 3 {
        rip_log!(
            "ripopt: Stalled for {} iterations without progress (alpha_p={:.2e}, pr={:.2e}, du={:.2e}), terminating",
            *stall_no_progress_count, state.alpha_primal, primal_inf, dual_inf
        );
    }
    StallDecision::Return(make_result(state, SolveStatus::NumericalError))
}

/// Track constraint-violation history and optionally short-circuit
/// with `LocalInfeasibility`.
///
/// Pushes the current `primal_inf` into `theta_history`, updates the
/// "ever feasible" flag, and (when proactive infeasibility detection
/// is enabled) looks for a stagnated θ with a near-stationary ∇θ
/// over a 100-iteration window. Returns `Some(SolveResult)` when
/// infeasibility is declared.
///
/// Ipopt uses its `RestoFilterConvCheck` for a broadly analogous
/// purpose; this path is ripopt-specific and fires before the main
/// IPM would otherwise burn iterations waiting for restoration to
/// reach the same conclusion.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition.
#[allow(clippy::too_many_arguments)]
fn track_feasibility_and_detect_infeasibility(
    state: &SolverState,
    options: &SolverOptions,
    iteration: usize,
    primal_inf: f64,
    theta_history: &mut Vec<f64>,
    theta_history_len: usize,
    theta_stall_count: &mut usize,
    ever_feasible: &mut bool,
) -> Option<SolveResult> {
    let n = state.n;
    let m = state.m;

    if theta_history.len() >= theta_history_len {
        theta_history.remove(0);
    }
    theta_history.push(primal_inf);

    if primal_inf < options.constr_viol_tol {
        *ever_feasible = true;
    }

    // Proactive infeasibility detection: if θ has stagnated over the history
    // window AND ‖∇θ‖_∞ is near zero, declare infeasibility instead of
    // waiting for restoration to reach the same conclusion.
    if options.proactive_infeasibility_detection
        && !*ever_feasible
        && m > 0
        && iteration >= 50
        && primal_inf > options.constr_viol_tol
        && theta_history.len() >= theta_history_len
    {
        let theta_min_h = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
        let theta_max_h = theta_history.iter().cloned().fold(0.0f64, f64::max);
        if theta_max_h > 0.0 && (theta_max_h - theta_min_h) < 0.01 * primal_inf {
            *theta_stall_count += 1;
        } else {
            *theta_stall_count = 0;
        }
        if *theta_stall_count >= 10 {
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
                return Some(make_result(state, SolveStatus::LocalInfeasibility));
            }
            // Stationarity not met — reset counter to check again next window.
            *theta_stall_count = 0;
        }
    } else if *ever_feasible {
        *theta_stall_count = 0;
    }
    None
}

/// Detect and recover from dual-infeasibility stagnation.
///
/// Tracks the best dual-infeasibility seen so far (`last_good_du`,
/// `last_good_iter`) and, if `du` has failed to halve for ≥500
/// iterations while a best-du iterate is available, restores that
/// iterate and resets the filter/μ/inertia state for a fresh start.
///
/// Catches restoration-cycling failure modes where the main IPM and
/// restoration NLP ping-pong and drift away from a near-converged
/// point for thousands of iterations. No Ipopt parallel — this is a
/// ripopt-specific safety net.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition. Returns `Some(SolveResult)` only when the restored
/// point already meets the near-tolerance acceptable level.
#[allow(clippy::too_many_arguments)]
fn handle_dual_stagnation<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    iteration: usize,
    filter: &mut Filter,
    mu_state: &mut MuState,
    inertia_params: &mut InertiaCorrectionParams,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    last_good_du: &mut f64,
    last_good_iter: &mut usize,
    triggered: &mut bool,
    best_x: &Option<Vec<f64>>,
    best_du_x: &Option<Vec<f64>>,
    best_du_y: &Option<Vec<f64>>,
    best_du_zl: &Option<Vec<f64>>,
    best_du_zu: &Option<Vec<f64>>,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> Option<SolveResult> {
    if iteration == 0 {
        return None;
    }
    let n = state.n;
    let m = state.m;

    let current_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    if current_du < 0.5 * *last_good_du {
        *last_good_du = current_du;
        *last_good_iter = iteration;
    }

    let stall_iters = iteration.saturating_sub(*last_good_iter);
    if stall_iters < 500
        || *triggered
        || current_du <= 100.0 * options.tol
        || best_x.is_none()
    {
        return None;
    }

    let bdx = match best_du_x {
        Some(v) => v,
        None => return None,
    };
    log::debug!(
        "Dual stagnation at iter {}: du={:.2e}, restoring best-du point (du={:.2e} at iter {})",
        iteration, current_du, *last_good_du, *last_good_iter
    );
    state.x.copy_from_slice(bdx);
    if let Some(ref bdy) = best_du_y { state.y.copy_from_slice(bdy); }
    if let Some(ref bdzl) = best_du_zl { state.z_l.copy_from_slice(bdzl); }
    if let Some(ref bdzu) = best_du_zu { state.z_u.copy_from_slice(bdzu); }
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    update_lbfgs_hessian(lbfgs_state, state);

    // Reset filter and bump mu for a fresh start from the good point.
    filter.reset();
    let theta_restart = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_restart);
    state.mu = (state.mu * 100.0).max(1e-4).min(1e-1);
    if options.mu_strategy_adaptive {
        mu_state.mode = MuMode::Free;
    }
    mu_state.first_iter_in_mode = true;
    mu_state.consecutive_restoration_failures = 0;
    inertia_params.delta_w_last = 0.0;

    // Check if the restored point meets acceptable convergence.
    let rest_pr = state.constraint_violation();
    let rest_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let rest_co = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );
    let s_max = 100.0_f64;
    let mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
    let s_d = if (m + 2 * n) > 0 {
        (s_max.max(mult_sum / (m + 2 * n) as f64) / s_max).min(1e4)
    } else { 1.0 };
    let near_tol = 100.0 * options.tol;
    let du_tol = (near_tol * s_d).max(1e-2);
    let co_tol = (near_tol * s_d).max(1e-2);
    let pr_tol = near_tol.max(10.0 * options.constr_viol_tol);
    if rest_pr <= pr_tol && rest_du <= du_tol && rest_co <= co_tol {
        log::debug!(
            "Restored best-du point passes near-tolerance (pr={:.2e}, du={:.2e}, co={:.2e})",
            rest_pr, rest_du, rest_co
        );
        return Some(make_result(state, SolveStatus::NumericalError));
    }

    *triggered = true;
    None
}

/// Emit one row of the per-iteration TSV trace (for the
/// direction-diff harness with Ipopt's IntermediateCallback).
///
/// Gated on `trace::is_enabled()` (RIP_TRACE_TSV env var). No-op in
/// production runs. Also resets the per-iteration scratch values
/// (`last_alpha_primal_max`, `last_tau_used`, `last_soc_accepted`)
/// so the next iteration starts clean.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop
/// decomposition (pre-work step 2). Pure side-effect function.
#[allow(clippy::too_many_arguments)]
fn emit_trace_row_if_enabled(
    state: &SolverState,
    iteration: usize,
    primal_inf: f64,
    dual_inf: f64,
    compl_inf: f64,
    ls_steps: usize,
    inertia_params: &InertiaCorrectionParams,
    last_mehrotra_sigma: Option<f64>,
    last_alpha_primal_max: &mut Option<f64>,
    last_tau_used: &mut Option<f64>,
    last_soc_accepted: &mut bool,
) {
    if !trace::is_enabled() {
        return;
    }
    let n = state.n;
    let dx_inf = state.dx.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let dzl_inf = state.dz_l.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let dzu_inf = state.dz_u.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    // log10(max Σ_i / min Σ_i) where Σ_i = z_l_i/s_l_i + z_u_i/s_u_i.
    // High values (~10+) signal ill-conditioning of the condensed KKT at
    // low μ and are an α=1/low-μ direction-quality suspect.
    let (sigma_min, sigma_max_val) = {
        let mut mn = f64::INFINITY;
        let mut mx = 0.0_f64;
        for i in 0..n {
            let mut s_i = 0.0_f64;
            if state.x_l[i].is_finite() {
                let s = (state.x[i] - state.x_l[i]).max(1e-20);
                s_i += state.z_l[i] / s;
            }
            if state.x_u[i].is_finite() {
                let s = (state.x_u[i] - state.x[i]).max(1e-20);
                s_i += state.z_u[i] / s;
            }
            if s_i > 0.0 {
                mn = mn.min(s_i);
                mx = mx.max(s_i);
            }
        }
        (mn, mx)
    };
    let sigma_cond = if sigma_min.is_finite() && sigma_min > 0.0 && sigma_max_val > 0.0 {
        (sigma_max_val / sigma_min).log10()
    } else {
        f64::NAN
    };
    trace::emit(&trace::TraceRow {
        iter: iteration,
        obj: state.obj / state.obj_scaling,
        inf_pr: primal_inf,
        inf_du: dual_inf,
        compl: compl_inf,
        mu: state.mu,
        alpha_pr: state.alpha_primal,
        alpha_du: state.alpha_dual,
        alpha_aff_p: f64::NAN,
        alpha_aff_d: f64::NAN,
        mu_aff: f64::NAN,
        sigma: last_mehrotra_sigma.unwrap_or(f64::NAN),
        mu_pc: f64::NAN,
        delta_w: inertia_params.delta_w_last,
        delta_c: 0.0,
        dx_inf,
        dzl_inf,
        dzu_inf,
        mcc_iters: 0,
        ls: ls_steps as u32,
        accepted: true,
        alpha_primal_max: last_alpha_primal_max.unwrap_or(f64::NAN),
        tau_used: last_tau_used.unwrap_or(f64::NAN),
        sigma_cond,
        soc_accepted: *last_soc_accepted,
    });
    *last_alpha_primal_max = None;
    *last_tau_used = None;
    *last_soc_accepted = false;
}

/// Populate the per-iteration IterateSnapshot and invoke the user
/// intermediate callback.
///
/// Sets the global "current iterate" for GetCurrentIterate / GetViolations
/// access inside the callback, then runs it. If the callback returns false
/// the solver must terminate with `UserRequestedStop`; returns
/// `Some(SolveResult)` in that case. The snapshot is cleared before return
/// regardless.
///
/// Ipopt parallel: IpIpoptAlg intermediate-callback invocation (around
/// IpIpoptAlg.cpp:560–610 and IpOptionsList intermediate_callback wiring).
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2).
fn populate_snapshot_and_invoke_callback(
    state: &SolverState,
    iteration: usize,
    primal_inf: f64,
    dual_inf: f64,
    prev_ic_delta_w: f64,
    ls_steps: usize,
    options: &SolverOptions,
) -> Option<SolveResult> {
    use crate::intermediate::IterateSnapshot;
    let n = state.n;
    let m = state.m;

    // Populate iterate snapshot for GetCurrentIterate/Violations access
    let mut x_l_viol = vec![0.0; n];
    let mut x_u_viol = vec![0.0; n];
    let mut compl_xl = vec![0.0; n];
    let mut compl_xu = vec![0.0; n];
    for i in 0..n {
        if state.x_l[i].is_finite() {
            x_l_viol[i] = (state.x_l[i] - state.x[i]).max(0.0);
            compl_xl[i] = (state.x[i] - state.x_l[i]) * state.z_l[i];
        }
        if state.x_u[i].is_finite() {
            x_u_viol[i] = (state.x[i] - state.x_u[i]).max(0.0);
            compl_xu[i] = (state.x_u[i] - state.x[i]) * state.z_u[i];
        }
    }
    // grad_lag = grad_f + J^T y - z_l + z_u
    let mut grad_lag = state.grad_f.clone();
    for (k, (&r, &c)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        grad_lag[c] += state.jac_vals[k] * state.y[r];
    }
    for i in 0..n {
        grad_lag[i] -= state.z_l[i];
        grad_lag[i] += state.z_u[i];
    }
    let mut constr_viol = vec![0.0; m];
    let mut compl_g_vec = vec![0.0; m];
    for i in 0..m {
        if state.g_l[i].is_finite() && state.g[i] < state.g_l[i] {
            constr_viol[i] = state.g_l[i] - state.g[i];
        } else if state.g_u[i].is_finite() && state.g[i] > state.g_u[i] {
            constr_viol[i] = state.g[i] - state.g_u[i];
        }
        // Complementarity: lambda_i * c_i where c_i is the active constraint slack
        if state.g_l[i].is_finite() && state.g_u[i].is_finite() {
            // Equality or range: use min slack
            compl_g_vec[i] = state.y[i] * (state.g[i] - state.g_l[i]).min(state.g_u[i] - state.g[i]);
        } else if state.g_l[i].is_finite() {
            compl_g_vec[i] = state.y[i] * (state.g[i] - state.g_l[i]);
        } else if state.g_u[i].is_finite() {
            compl_g_vec[i] = state.y[i] * (state.g_u[i] - state.g[i]);
        }
    }
    crate::intermediate::set_current_iterate(Some(IterateSnapshot {
        x: state.x.clone(),
        z_l: state.z_l.clone(),
        z_u: state.z_u.clone(),
        g: state.g.clone(),
        lambda: state.y.clone(),
        x_l_violation: x_l_viol,
        x_u_violation: x_u_viol,
        compl_x_l: compl_xl,
        compl_x_u: compl_xu,
        grad_lag_x: grad_lag,
        constraint_violation: constr_viol,
        compl_g: compl_g_vec,
    }));

    // Invoke intermediate callback (if registered)
    let d_norm = state.dx.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    let continue_ok = crate::intermediate::invoke_intermediate(
        0, // alg_mod: 0 = regular mode (not in restoration)
        iteration,
        state.obj / state.obj_scaling,
        primal_inf,
        dual_inf,
        state.mu,
        d_norm,
        prev_ic_delta_w,
        state.alpha_dual,
        state.alpha_primal,
        ls_steps,
    );
    crate::intermediate::set_current_iterate(None);
    if !continue_ok {
        if options.print_level >= 5 {
            rip_log!("ripopt: User requested stop via intermediate callback");
        }
        return Some(make_result(state, SolveStatus::UserRequestedStop));
    }
    None
}

/// Emit a single iteration-row line to the solver log.
///
/// Ipopt parallel: `IpIpoptAlg::OutputIteration` (IpIpoptAlg.cpp:520–560).
/// Header is reprinted every 25 data rows.
///
/// Extracted from `solve_ipm` as part of the v0.8 main-loop decomposition
/// (pre-work step 2). Pure side-effect function — no solver state mutation
/// beyond incrementing `log_line_count`.
fn log_iteration_row(
    iteration: usize,
    state: &SolverState,
    primal_inf: f64,
    dual_inf: f64,
    ls_steps: usize,
    log_line_count: &mut usize,
    options: &SolverOptions,
) {
    if options.print_level < 3 {
        return;
    }
    // Reprint header every 25 data rows for readability
    if *log_line_count > 0 && *log_line_count % 25 == 0 {
        rip_log!(
            "{:>4}  {:>14}  {:>10}  {:>10}  {:>7}  {:>8}  {:>8}  {:>3}",
            "iter", "objective", "inf_pr", "inf_du", "lg(mu)", "alpha_pr", "alpha_du", "ls"
        );
    }
    rip_log!(
        "{:>4}  {:>14.7e}  {:>10.2e}  {:>10.2e}  {:>10.2e}  {:>8.2e}  {:>8.2e}  {:>3}",
        iteration,
        state.obj / state.obj_scaling,
        primal_inf,
        dual_inf,
        state.mu,
        state.alpha_primal,
        state.alpha_dual,
        ls_steps,
    );
    *log_line_count += 1;
}

/// Pick the linear solver for the KKT factorization.
///
/// Dense (`DenseLdl`) for small systems (`n + m < sparse_threshold`). For
/// sparse, when the `rmumps` feature is enabled and the problem has
/// constraints (`m > 0`), use the KKT-aware `MultifrontalLdl::new_kkt(n)`
/// which enables CB pivot search for numerically stable primal-dual 2×2
/// pivots. Otherwise fall back to the user's `options.linear_solver` choice.
fn select_linear_solver(use_sparse: bool, n: usize, m: usize, options: &SolverOptions) -> Box<dyn LinearSolver> {
    if use_sparse {
        #[cfg(feature = "rmumps")]
        {
            let _ = n;
            if m > 0 {
                Box::new(MultifrontalLdl::new_kkt(n))
            } else {
                new_sparse_solver_with_choice(options.linear_solver)
            }
        }
        #[cfg(not(feature = "rmumps"))]
        {
            let _ = (n, m);
            new_sparse_solver_with_choice(options.linear_solver)
        }
    } else {
        Box::new(DenseLdl::new())
    }
}

/// Estimate whether sparse condensed KKT should be disabled in favor of the
/// full augmented system based on Jacobian structure. The Schur complement
/// `J^T·D^{-1}·J` has at most `Σ k_i*(k_i+1)/2` nonzeros (before dedup),
/// where `k_i` is the nnz in row `i`. If that exceeds `2×` the nnz of the
/// augmented KKT (`hess_nnz + jac_nnz + n`), factoring the Schur complement
/// costs more than factoring the full system, so return `true`.
fn estimate_schur_density_disable<P: NlpProblem>(
    problem: &P,
    n: usize,
    m: usize,
    use_sparse: bool,
    options: &SolverOptions,
) -> bool {
    if !(use_sparse && m > 0) {
        return false;
    }
    let (jac_rows_est, _) = problem.jacobian_structure();
    let mut row_nnz = vec![0usize; m];
    for &r in &jac_rows_est {
        row_nnz[r] += 1;
    }
    let schur_nnz_upper: usize = row_nnz.iter().map(|&k| k * (k + 1) / 2).sum();
    let (hess_rows_est, _) = problem.hessian_structure();
    let augmented_nnz = hess_rows_est.len() + jac_rows_est.len() + n;
    let disable = schur_nnz_upper > 2 * augmented_nnz;
    if disable && options.print_level >= 3 {
        rip_log!(
            "ripopt: Disabling sparse condensed KKT: Schur complement nnz estimate ({}) > 2× augmented KKT nnz ({})",
            schur_nnz_upper, augmented_nnz
        );
    }
    disable
}

/// Detect constraints that are linear in `x` (constant Jacobian, zero
/// contribution to the Hessian). When `options.detect_linear_constraints`
/// is set, run `crate::linearity::detect_linear_constraints` on the
/// original unscaled problem. Returns `Some(flags)` only when at least one
/// linear constraint is found; `None` otherwise.
fn detect_linear_constraint_flags<P: NlpProblem>(
    problem: &P,
    options: &SolverOptions,
    x0: &[f64],
    m_sc: usize,
) -> Option<Vec<bool>> {
    if !(options.detect_linear_constraints && m_sc > 0) {
        return None;
    }
    let flags = crate::linearity::detect_linear_constraints(problem, x0);
    let n_linear = flags.iter().filter(|&&f| f).count();
    if n_linear == 0 {
        return None;
    }
    if options.print_level >= 5 {
        rip_log!(
            "ripopt: Detected {}/{} linear constraints (Hessian contribution skipped)",
            n_linear, m_sc
        );
    }
    Some(flags)
}

/// Initialize constraint slack barrier multipliers `v_l`, `v_u` (Ipopt's
/// `v_L`, `v_U`). For each inequality constraint side,
/// `v = mu_init / max(slack, 1e-20)`. Equality rows (`g_l ≈ g_u`) are
/// skipped. Mirrors Ipopt's `IpDefaultIterateInitializer.cpp`.
fn initialize_constraint_slack_multipliers(state: &mut SolverState, m: usize, options: &SolverOptions) {
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
}

/// Apply user-provided warm-start multipliers (`y`, `z_l`, `z_u`) and then
/// run `WarmStartInitializer::initialize`, which pushes `x` off the bounds
/// and rescales `(z_l, z_u)` if needed so (x − x_l)·z_l ≈ μ is well-posed.
/// No-op when `options.warm_start` is false.
fn apply_warm_start(state: &mut SolverState, options: &SolverOptions) {
    if !options.warm_start {
        return;
    }
    if let Some(ref init_y) = options.warm_start_y {
        let len = init_y.len().min(state.y.len());
        state.y[..len].copy_from_slice(&init_y[..len]);
    }
    if let Some(ref init_z_l) = options.warm_start_z_l {
        let len = init_z_l.len().min(state.z_l.len());
        state.z_l[..len].copy_from_slice(&init_z_l[..len]);
    }
    if let Some(ref init_z_u) = options.warm_start_z_u {
        let len = init_z_u.len().min(state.z_u.len());
        state.z_u[..len].copy_from_slice(&init_z_u[..len]);
    }
    state.mu = WarmStartInitializer::initialize(
        &mut state.x,
        &mut state.z_l,
        &mut state.z_u,
        &state.x_l,
        &state.x_u,
        options,
    );
}

/// Evaluate the NLP at the initial point. If the evaluation fails, produces
/// NaN/Inf in `obj` or `grad_f`, try up to three bound-push perturbations
/// (1%, 10%, 50% of each bound range) and retry. If every perturbation also
/// fails, return `Err(SolveResult)` with `SolveStatus::EvaluationError`.
fn initial_evaluate_with_recovery<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    lbfgs_state: &mut Option<LbfgsIpmState>,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
    n: usize,
    options: &SolverOptions,
) -> Result<(), SolveResult> {
    let init_eval_ok = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    update_lbfgs_hessian(lbfgs_state, state);

    if init_eval_ok && !state.obj.is_nan() && !state.obj.is_infinite()
        && !state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
    {
        return Ok(());
    }

    let x_saved = state.x.clone();
    for &push_factor in &[1e-2, 1e-1, 0.5] {
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
        let perturb_ok = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        update_lbfgs_hessian(lbfgs_state, state);
        if perturb_ok && !state.obj.is_nan() && !state.obj.is_infinite()
            && !state.grad_f.iter().any(|v| v.is_nan() || v.is_infinite())
        {
            return Ok(());
        }
    }
    Err(make_result(state, SolveStatus::EvaluationError))
}

/// Compute gradient-based NLP scaling at the initial point `x0`, matching
/// Ipopt's `nlp_scaling_method = gradient-based`:
///
///   * Objective: if `||∇f(x0)||_∞ > 100`, set `obj_scaling = 100/||∇f||_∞`,
///     clamped below by `1e-2`.
///   * Constraints: row-wise on `J(x0)`. If `max_j |J_{ij}| > 100`, set
///     `g_scaling[i] = 100/max_j |J_{ij}|`, clamped below by `1e-2`.
///
/// User-provided scalings (`options.user_obj_scaling`, `options.user_g_scaling`)
/// take priority — when either is set, skips automatic scaling entirely.
///
/// Constraint scaling is skipped when the initial constraint violation exceeds
/// `1e6` (Ipopt's threshold — at highly infeasible points J has no signal).
fn compute_nlp_scaling<P: NlpProblem>(
    problem: &P,
    options: &SolverOptions,
    x0: &[f64],
    jac_rows_sc: &[usize],
) -> (f64, Vec<f64>) {
    let n_sc = problem.num_variables();
    let m_sc = problem.num_constraints();

    if options.user_obj_scaling.is_some() || options.user_g_scaling.is_some() {
        let os = options.user_obj_scaling.unwrap_or(1.0);
        let gs = options.user_g_scaling.clone().unwrap_or_else(|| vec![1.0; m_sc]);
        return (os, gs);
    }

    let nlp_scaling_max_gradient = 100.0;
    let nlp_scaling_min_value = 1e-2;
    let mut grad_f0 = vec![0.0; n_sc];
    let grad_ok = problem.gradient(x0, true, &mut grad_f0);
    let grad_max = if grad_ok {
        grad_f0.iter().map(|v| v.abs()).fold(0.0f64, f64::max)
    } else {
        0.0
    };
    let os = if grad_max > nlp_scaling_max_gradient && grad_max.is_finite() {
        (nlp_scaling_max_gradient / grad_max).max(nlp_scaling_min_value)
    } else {
        1.0
    };

    let mut gs = vec![1.0; m_sc];
    if m_sc > 0 {
        let mut g0_sc = vec![0.0; m_sc];
        let constr_ok = problem.constraints(x0, false, &mut g0_sc);
        let mut g_l_sc = vec![0.0; m_sc];
        let mut g_u_sc = vec![0.0; m_sc];
        problem.constraint_bounds(&mut g_l_sc, &mut g_u_sc);
        let init_cv = if constr_ok {
            convergence::primal_infeasibility(&g0_sc, &g_l_sc, &g_u_sc)
        } else {
            f64::INFINITY
        };

        if init_cv < 1e6 {
            let mut jac_vals0 = vec![0.0; jac_rows_sc.len()];
            if problem.jacobian_values(x0, false, &mut jac_vals0) {
                let mut row_max = vec![0.0f64; m_sc];
                for (idx, &row) in jac_rows_sc.iter().enumerate() {
                    let v = jac_vals0[idx].abs();
                    if v.is_finite() && v > row_max[row] {
                        row_max[row] = v;
                    }
                }
                for i in 0..m_sc {
                    if row_max[i] > nlp_scaling_max_gradient {
                        gs[i] = (nlp_scaling_max_gradient / row_max[i]).max(nlp_scaling_min_value);
                    }
                }
            }
        }
    }
    (os, gs)
}

/// Core IPM solver implementation.
fn solve_ipm<P: NlpProblem>(problem: &P, options: &SolverOptions) -> SolveResult {
    let n_sc = problem.num_variables();
    let m_sc = problem.num_constraints();

    let mut x0 = vec![0.0; n_sc];
    problem.initial_point(&mut x0);

    let (jac_rows_sc, _) = problem.jacobian_structure();
    let (obj_scaling, g_scaling) = compute_nlp_scaling(problem, options, &x0, &jac_rows_sc);

    if options.print_level >= 5
        && (obj_scaling != 1.0 || g_scaling.iter().any(|&s| s != 1.0))
    {
        let n_scaled_g = g_scaling.iter().filter(|&&s| s != 1.0).count();
        rip_log!(
            "ripopt: NLP scaling: obj_scaling={:.4e}, {}/{} constraints scaled",
            obj_scaling, n_scaled_g, m_sc
        );
    }

    let linear_constraints = detect_linear_constraint_flags(problem, options, &x0, m_sc);

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
            rip_log!("ripopt: Using L-BFGS Hessian approximation (limited-memory mode)");
        }
        Some(LbfgsIpmState::new(n))
    } else {
        None
    };

    apply_warm_start(&mut state, options);

    // Initialize linear solver — use sparse for large KKT systems
    let use_sparse = (n + m) >= options.sparse_threshold;
    let mut lin_solver = select_linear_solver(use_sparse, n, m, options);
    let mut inertia_params = InertiaCorrectionParams::default();
    let mut restoration = RestorationPhase::new(500);

    let mut disable_sparse_condensed = estimate_schur_density_disable(problem, n, m, use_sparse, options);
    // Flag set by bandwidth detection: when the sparse condensed Schur complement
    // is essentially dense, switch to dense condensed KKT (n×n) for all subsequent
    // iterations. This avoids the catastrophic rmumps fill-in on PDE problems.
    let _use_dense_condensed_fallback = false;

    // Initialize filter
    let mut filter = Filter::new(1e4);

    // Mehrotra centering parameter from the last iteration's predictor step.
    // Used in the Free-mode mu update: when sigma is available, mu = sigma * mu_current
    // gives a more aggressive (and adaptive) decrease than the Loqo oracle.
    let mut last_mehrotra_sigma: Option<f64> = None;
    // Per-iteration trace intermediates captured during the line-search /
    // direction-compute sub-phases, drained into the TSV at iteration-end.
    let mut last_alpha_primal_max: Option<f64> = None;
    let mut last_tau_used: Option<f64> = None;
    let mut last_soc_accepted: bool = false;

    // Free/fixed mu mode state (replaces ad-hoc stall recovery)
    let mut mu_state = MuState::new();
    // Monotone mu strategy: start in Fixed mode and never switch to Free
    if !options.mu_strategy_adaptive {
        mu_state.mode = MuMode::Fixed;
    }

    // Wall-clock time limit
    let start_time = Instant::now();
    let deadline = if options.max_wall_time > 0.0 {
        Some(start_time + Duration::from_secs_f64(options.max_wall_time))
    } else {
        None
    };

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

    // Line-search backtrack count for the previous iteration (printed in table).
    let mut ls_steps: usize = 0;
    // Hessian regularization delta from previous iteration (for intermediate callback).
    let mut prev_ic_delta_w: f64 = 0.0;

    // Primal divergence detection: track consecutive iterations where pr is growing.
    // When pr grows steadily post-restoration, re-trigger restoration rather than
    // continuing for many iterations with worsening feasibility.
    let mut pr_prev_for_divergence: f64 = f64::INFINITY;
    let mut pr_at_divergence_start: f64 = f64::INFINITY;
    let mut consecutive_pr_increase: usize = 0;

    // Consecutive iterations with obj < -1e20 for robust unbounded detection
    let mut consecutive_unbounded: usize = 0;

    // Consecutive iterations where theta (primal infeasibility) stagnated.
    // Used by proactive infeasibility detection to exit early.
    let mut theta_stall_count: usize = 0;

    // Best feasible point tracking: save the best (lowest obj) point that is feasible
    let mut best_x: Option<Vec<f64>> = None;
    let mut best_obj: f64 = f64::INFINITY;

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

    // Strategy 1: Iterate averaging for oscillation recovery
    const AVG_WINDOW: usize = 6;
    let mut du_history: Vec<f64> = Vec::with_capacity(AVG_WINDOW + 1);
    let mut iterate_history: Vec<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> = Vec::new(); // (x, y, z_l, z_u)
    let mut tried_iterate_averaging: bool = false;

    // Strategy 2: Damped multiplier updates when oscillation detected
    let mut prev_dy: Option<Vec<f64>> = None;
    let mut dy_sign_change_count: Vec<u8> = vec![0u8; m]; // per-component consecutive sign change count

    // Strategy 3: Active set reduced KKT solve
    let mut tried_active_set: bool = false;

    // Strategy 4: Complementarity polishing — force mu small when compl is bottleneck
    let mut _tried_compl_polish: bool = false;

    // Initial evaluation with NaN/Inf recovery by bound-push perturbation.
    if let Err(result) = initial_evaluate_with_recovery(
        &mut state, problem, &mut lbfgs_state, linear_constraints.as_deref(), lbfgs_mode, n, options,
    ) {
        return result;
    }

    initialize_constraint_slack_multipliers(&mut state, m, options);

    // Set filter parameters based on initial constraint violation
    let theta_init = state.constraint_violation();
    filter.set_theta_min_from_initial(theta_init);

    // Print iteration table header (shown at print_level >= 3, reprinted every 25 rows)
    let mut log_line_count: usize = 0;
    if options.print_level >= 3 {
        rip_log!(
            "{:>4}  {:>14}  {:>10}  {:>10}  {:>10}  {:>8}  {:>8}  {:>3}",
            "iter", "objective", "inf_pr", "inf_du", "mu", "alpha_pr", "alpha_du", "ls"
        );
    }

    if options.print_level >= 5 {
        rip_log!("ripopt: Starting main loop (n={}, m={})", n, m);
    }

    // Main IPM loop
    for iteration in 0..options.max_iter {
        state.iter = iteration;

        // Early-stall timeout scaled by problem size: medium-scale problems
        // (n+m > 1000) can legitimately spend 30-60s on restoration or line
        // search during early iterations.
        let early_timeout = options.early_stall_timeout * ((n + m) as f64 / 200.0).max(1.0);
        if let Some(result) = check_time_limits(&state, iteration, start_time, early_timeout, options) {
            return result;
        }

        // Dual stagnation detection (runs every iteration, including restoration).
        // Catches restoration cycling that drifts from a near-converged point.
        if let Some(result) = handle_dual_stagnation(
            &mut state,
            problem,
            options,
            iteration,
            &mut filter,
            &mut mu_state,
            &mut inertia_params,
            &mut lbfgs_state,
            &mut dual_stall_last_good_du,
            &mut dual_stall_last_good_iter,
            &mut dual_stall_triggered,
            &best_x,
            &best_du_x,
            &best_du_y,
            &best_du_zl,
            &best_du_zu,
            linear_constraints.as_deref(),
            lbfgs_mode,
        ) {
            return result;
        }

        let OptimalityMeasures {
            primal_inf,
            primal_inf_max,
            dual_inf,
            dual_inf_unscaled,
            compl_inf,
        } = compute_optimality_measures(&state);

        log_iteration_row(
            iteration,
            &state,
            primal_inf,
            dual_inf,
            ls_steps,
            &mut log_line_count,
            options,
        );

        emit_trace_row_if_enabled(
            &state,
            iteration,
            primal_inf,
            dual_inf,
            compl_inf,
            ls_steps,
            &inertia_params,
            last_mehrotra_sigma,
            &mut last_alpha_primal_max,
            &mut last_tau_used,
            &mut last_soc_accepted,
        );

        if let Some(result) = populate_snapshot_and_invoke_callback(
            &state,
            iteration,
            primal_inf,
            dual_inf,
            prev_ic_delta_w,
            ls_steps,
            options,
        ) {
            return result;
        }

        // Compute multiplier scaling (also used by the consecutive-acceptable
        // tracker further below, so kept here rather than inside the helper).
        let multiplier_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let multiplier_count = m + 2 * n;

        let mut conv_ws = ConvergenceWorkspace {
            du_history: &mut du_history,
            iterate_history: &mut iterate_history,
            tried_iterate_averaging: &mut tried_iterate_averaging,
            tried_active_set: &mut tried_active_set,
            tried_compl_polish: &mut _tried_compl_polish,
        };
        if let Some(result) = check_convergence_and_handle_promotions(
            &mut state,
            problem,
            options,
            primal_inf_max,
            dual_inf,
            dual_inf_unscaled,
            compl_inf,
            multiplier_sum,
            multiplier_count,
            &mut conv_ws,
            &timings,
            iteration,
            ipm_start,
            linear_constraints.as_deref(),
            lbfgs_mode,
        ) {
            return result;
        }

        let s_d_for_acc = track_consecutive_acceptable(
            &mut state,
            primal_inf,
            dual_inf,
            dual_inf_unscaled,
            compl_inf,
            multiplier_sum,
        );

        update_best_du_iterate(
            &state,
            dual_inf,
            &mut best_du_val,
            &mut best_du_x,
            &mut best_du_y,
            &mut best_du_zl,
            &mut best_du_zu,
        );

        if let Some(result) = track_feasibility_and_detect_infeasibility(
            &state,
            options,
            iteration,
            primal_inf,
            &mut theta_history,
            theta_history_len,
            &mut theta_stall_count,
            &mut ever_feasible,
        ) {
            return result;
        }

        if let Some(result) = detect_unbounded(
            &state,
            options,
            primal_inf,
            &mut consecutive_unbounded,
        ) {
            return result;
        }

        // Overall-progress stall detection — see helper doc comment.
        match detect_and_handle_progress_stall(
            &mut state,
            problem,
            options,
            iteration,
            primal_inf,
            primal_inf_max,
            dual_inf,
            compl_inf,
            s_d_for_acc,
            &mut filter,
            &mut mu_state,
            &mut stall_no_progress_count,
            &mut stall_best_pr,
            &mut stall_best_du,
            best_du_val,
            &best_du_x,
            &best_du_y,
            &best_du_zl,
            &best_du_zu,
            linear_constraints.as_deref(),
            lbfgs_mode,
        ) {
            StallDecision::Return(r) => return r,
            StallDecision::Continue => continue,
            StallDecision::Proceed => {}
        }
        let force_restoration = detect_primal_divergence(
            options,
            iteration,
            primal_inf,
            &mut consecutive_pr_increase,
            &mut pr_prev_for_divergence,
            &mut pr_at_divergence_start,
            m,
        );

        let t_kkt = Instant::now();
        let AssembledKkt {
            sigma,
            use_sparse_condensed,
            condensed_system,
            mut sparse_condensed_system,
            mut kkt_system_opt,
        } = assemble_kkt_systems(&state, n, m, use_sparse, disable_sparse_condensed);
        timings.kkt_assembly += t_kkt.elapsed();

        // On first iteration with sparse condensed, detect bandwidth and pick
        // the right downstream solver (full augmented / banded / sparse).
        if iteration == 0 {
            adjust_sparse_condensed_bandwidth(
                &state,
                problem,
                options,
                n,
                m,
                use_sparse,
                use_sparse_condensed,
                &sigma,
                &mut sparse_condensed_system,
                &mut kkt_system_opt,
                &mut lin_solver,
                &mut disable_sparse_condensed,
            );
        }

        let (ic_delta_w, ic_delta_c) = match factor_kkt_with_recovery(
            &mut state,
            problem,
            options,
            iteration,
            n,
            m,
            use_sparse,
            &mut kkt_system_opt,
            lin_solver.as_mut(),
            &mut inertia_params,
            &mut lbfgs_state,
            &mut filter,
            &mut restoration,
            &mut timings,
            &mut prev_ic_delta_w,
            linear_constraints.as_deref(),
            lbfgs_mode,
            deadline,
        ) {
            FactorDecision::Proceed { ic_delta_w, ic_delta_c } => (ic_delta_w, ic_delta_c),
            FactorDecision::Continue => continue,
            FactorDecision::Return(result) => return result,
        };

        // Solve for search direction via the three-way KKT dispatch.
        let t_dir = Instant::now();
        let (dx, dy);
        let mut cond_solver_for_soc: Option<DenseLdl>;
        let mehrotra_aff: Option<(Vec<f64>, Vec<f64>, Vec<f64>, f64)>;
        match solve_for_search_direction(
            &mut state,
            problem,
            options,
            iteration,
            n,
            m,
            use_sparse,
            &condensed_system,
            &sparse_condensed_system,
            &mut kkt_system_opt,
            lin_solver.as_mut(),
            &sigma,
            &mut inertia_params,
            ic_delta_w,
            ic_delta_c,
            &filter,
            &mut restoration,
            &mut lbfgs_state,
            lbfgs_mode,
            linear_constraints.as_deref(),
            &mut last_mehrotra_sigma,
            deadline,
        ) {
            DirectionSolveDecision::Proceed {
                dx: dx_val,
                dy: dy_val,
                cond_solver_for_soc: cs,
                mehrotra_aff: ma,
            } => {
                dx = dx_val;
                dy = dy_val;
                cond_solver_for_soc = cs;
                mehrotra_aff = ma;
            }
            DirectionSolveDecision::Continue => continue,
            DirectionSolveDecision::Return(r) => return r,
        }

        timings.direction_solve += t_dir.elapsed();

        // Recover bound multiplier steps. If the Mehrotra corrector was applied,
        // use the cross-term-aware recovery at μ_pc so dz stays consistent with
        // the corrector complementarity equation; otherwise use the plain formula
        // at state.mu.
        // Mirror the RHS choice above: filter-LS mode does not apply the
        // Mehrotra cross-term, so dz recovery uses the plain Fiacco formula
        // dz_L[i] = (mu - z_L*s_L)/s_L - (z_L/s_L)*dx[i] at mu_new.
        let mu_for_dz = mehrotra_aff.as_ref().map(|t| t.3).unwrap_or(state.mu);
        let (dz_l, dz_u) = kkt::recover_dz(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, &dx, mu_for_dz,
        );

        state.dx = dx;
        state.dy = dy;
        state.dz_l = dz_l;
        state.dz_u = dz_u;

        apply_gondzio_mcc(
            &mut state,
            options,
            iteration,
            &mu_state,
            primal_inf,
            dual_inf,
            compl_inf,
            &kkt_system_opt,
            lin_solver.as_mut(),
        );

        let (tau, alpha_primal_max, alpha_dual_max) = compute_alpha_max(
            &state, options, &mu_state, primal_inf, dual_inf, compl_inf,
        );
        last_alpha_primal_max = Some(alpha_primal_max);
        last_tau_used = Some(tau);

        detect_tiny_step(
            &mut state,
            options,
            &mut mu_state,
            &mut filter,
            &mut consecutive_tiny_steps,
            alpha_primal_max,
            primal_inf,
        );

        // Line search
        let t_ls = Instant::now();
        let theta_current = primal_inf;
        let phi_current = state.barrier_objective(options);
        let grad_phi_step = state.barrier_directional_derivative(options);

        let mut step_accepted;
        let min_alpha = filter.compute_alpha_min(theta_current, grad_phi_step);
        // Snapshot x before the line search so we can roll back + halve α if the
        // post-step gradient / Jacobian / Hessian eval trips the NaN/Inf guard
        // (Ipopt's line search catches Eval_Error from these with α-backtracking,
        // not a hard abort — IpBacktrackingLineSearch.cpp:776-784, 1158, 1193).
        let x_pre_step = state.x.clone();

        match run_line_search_loop(
            &mut state,
            problem,
            options,
            &mut filter,
            &condensed_system,
            &mut cond_solver_for_soc,
            &sparse_condensed_system,
            &kkt_system_opt,
            lin_solver.as_mut(),
            alpha_primal_max,
            theta_current,
            phi_current,
            grad_phi_step,
            min_alpha,
            force_restoration,
            watchdog_active,
            iteration,
            n,
            m,
            start_time,
            early_timeout,
            &mut last_soc_accepted,
            &mut ls_steps,
        ) {
            LineSearchOutcome::StepAccepted => {
                step_accepted = true;
            }
            LineSearchOutcome::Rejected => {
                step_accepted = false;
            }
            LineSearchOutcome::Return(result) => return result,
        }

        if !step_accepted {
            step_accepted = attempt_soft_restoration(
                &mut state,
                problem,
                &mut filter,
                alpha_primal_max,
                alpha_dual_max,
                theta_current,
                phi_current,
            );
        }

        if !step_accepted {
            match run_post_ls_restoration_cascade(
                &mut state,
                problem,
                options,
                &mut filter,
                &mut mu_state,
                &mut inertia_params,
                &mut lbfgs_state,
                lbfgs_mode,
                linear_constraints.as_deref(),
                &mut restoration,
                iteration,
                n,
                m,
                start_time,
                deadline,
                early_timeout,
                ever_feasible,
                theta_current,
                phi_current,
                &theta_history,
                theta_history_len,
            ) {
                RestorationCascadeDecision::Continue => continue,
                RestorationCascadeDecision::Return(result) => return result,
            }
        }

        // Step was accepted — reset consecutive restoration failure counter
        mu_state.consecutive_restoration_failures = 0;

        match update_watchdog(
            &mut state,
            problem,
            options,
            iteration,
            alpha_primal_max,
            &mut filter,
            &mut lbfgs_state,
            &mut consecutive_shortened,
            &mut watchdog_active,
            &mut watchdog_trial_count,
            &mut watchdog_saved,
            linear_constraints.as_deref(),
            lbfgs_mode,
        ) {
            WatchdogDecision::Proceed => {}
            WatchdogDecision::Continue => continue,
        }

        timings.line_search += t_ls.elapsed();

        let mu_ks = update_dual_variables(
            &mut state,
            &mu_state,
            alpha_dual_max,
            &mut prev_dy,
            &mut dy_sign_change_count,
        );

        match reevaluate_after_step(
            &mut state,
            problem,
            options,
            &mut lbfgs_state,
            &mut filter,
            &mut restoration,
            &mut timings,
            &x_pre_step,
            linear_constraints.as_deref(),
            lbfgs_mode,
            deadline,
        ) {
            PostStepEvalDecision::Proceed => {}
            PostStepEvalDecision::Continue => continue,
            PostStepEvalDecision::Return(result) => return result,
        }

        reset_slack_multipliers(&mut state, mu_ks);
        track_best_feasible(&state, options, &mut best_obj, &mut best_x);

        // --- Barrier parameter update (free/fixed mode) ---
        update_barrier_parameter(
            &mut state,
            &mut mu_state,
            &mut filter,
            &mut last_mehrotra_sigma,
            options,
        );

        track_post_step_acceptable(&mut state, options);
    }

    finalize_after_max_iter(
        &state, options, n, m, ever_feasible, &theta_history, theta_history_len,
        &timings, ipm_start,
    )
}

/// After the main IPM loop terminates at `max_iter`, decide the final
/// `SolveStatus`:
///
///   1. Log MaxIter diagnostics at print_level ≥ 5.
///   2. Infeasibility check (only if never feasible): `||∇θ|| ≈ 0`
///      → `LocalInfeasibility`; theta stagnation history → `Infeasible`.
///   3. Print phase timing summary.
///   4. Strict Optimal at final iterate (catches MGH10LS-class zero-residual
///      stalls where the counter reset).
///   5. Relaxed Acceptable at final iterate (Ipopt default thresholds).
///   6. Fallback: `MaxIterations`.
fn finalize_after_max_iter(
    state: &SolverState,
    options: &SolverOptions,
    n: usize,
    m: usize,
    ever_feasible: bool,
    theta_history: &[f64],
    theta_history_len: usize,
    timings: &PhaseTimings,
    ipm_start: Instant,
) -> SolveResult {
    if options.print_level >= 5 {
        let final_primal_inf = convergence::primal_infeasibility_max(&state.g, &state.g_l, &state.g_u);
        let final_dual_inf = convergence::dual_infeasibility(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        let final_dual_inf_unscaled = convergence::dual_infeasibility_scaled(
            &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
            &state.y, &state.z_l, &state.z_u, n,
        );
        let final_compl = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
        );
        let mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let s_max: f64 = 100.0;
        let s_d = if (m + 2 * n) > 0 {
            ((s_max.max(mult_sum / (m + 2 * n) as f64)) / s_max).min(1e4)
        } else {
            1.0
        };
        rip_log!(
            "ripopt: MaxIter diag: pr={:.2e} du={:.2e}(t={:.2e}) du_u={:.2e}(t={:.0e}) co={:.2e}(t={:.2e}) mu={:.2e} sd={:.1} ac={}",
            final_primal_inf,
            final_dual_inf, options.tol * s_d,
            final_dual_inf_unscaled, options.dual_inf_tol,
            final_compl, 10.0 * options.compl_inf_tol,
            state.mu, s_d, state.consecutive_acceptable,
        );
    }

    // Infeasibility detection (only when never feasible).
    let final_theta = state.constraint_violation();
    if !ever_feasible && final_theta > options.constr_viol_tol {
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
            return make_result(state, SolveStatus::LocalInfeasibility);
        }

        if final_theta > 1e4 && theta_history.len() >= theta_history_len {
            let min_theta = theta_history.iter().cloned().fold(f64::INFINITY, f64::min);
            if final_theta > 0.01 * min_theta {
                return make_result(state, SolveStatus::Infeasible);
            }
        }
    }

    if options.print_level >= 5 {
        timings.print_summary(options.max_iter, ipm_start.elapsed());
    }

    let primal_inf = state.constraint_violation();
    let dual_inf = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, state.n,
    );
    let compl_inf = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );
    if primal_inf <= options.constr_viol_tol
        && dual_inf <= options.dual_inf_tol
        && compl_inf <= options.compl_inf_tol
        && primal_inf <= options.tol
        && dual_inf <= options.tol
        && compl_inf <= options.tol
    {
        return make_result(state, SolveStatus::Optimal);
    }
    if primal_inf <= 1e-2
        && dual_inf <= 1e10
        && compl_inf <= 1e-2
        && primal_inf <= 1e-6
        && dual_inf <= 1e-6
        && compl_inf <= 1e-6
    {
        return make_result(state, SolveStatus::Acceptable);
    }
    make_result(state, SolveStatus::MaxIterations)
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

    let kappa_soc = 0.99;
    let tau = (1.0 - state.mu).max(options.tau_min);

    // c_soc starts at curr_c (constraint residual at current iterate) and
    // accumulates alpha_primal_soc * trial_c each iter. Matches Ipopt
    // IpFilterLSAcceptor.cpp:555-569.
    let mut c_soc = vec![0.0; m];
    let mut latest_trial_c = vec![0.0; m];
    for i in 0..m {
        let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_equality || state.g_l[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_l[i];
            latest_trial_c[i] = g_trial[i] - state.g_l[i];
        } else if state.g_u[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_u[i];
            latest_trial_c[i] = g_trial[i] - state.g_u[i];
        }
    }

    let mut alpha_primal_soc = alpha;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

    for _soc_iter in 0..options.max_soc {
        // Accumulate: c_soc += alpha_primal_soc * latest_trial_c
        for i in 0..m {
            c_soc[i] += alpha_primal_soc * latest_trial_c[i];
        }

        let mut rhs_soc = kkt.rhs.clone();
        for i in 0..m {
            rhs_soc[n + i] = -c_soc[i];
        }

        let mut sol_soc = vec![0.0; n + m];
        if solver.solve(&rhs_soc, &mut sol_soc).is_err() {
            return None;
        }
        let dx_soc = &sol_soc[..n];

        // Fresh fraction-to-boundary alpha on dx_soc (Ipopt IpFilterLSAcceptor.cpp:617-620).
        let mut alpha_new: f64 = 1.0;
        for i in 0..n {
            if state.x_l[i].is_finite() && dx_soc[i] < 0.0 {
                let slack = state.x[i] - state.x_l[i];
                alpha_new = alpha_new.min(-tau * slack / dx_soc[i]);
            }
            if state.x_u[i].is_finite() && dx_soc[i] > 0.0 {
                let slack = state.x_u[i] - state.x[i];
                alpha_new = alpha_new.min(tau * slack / dx_soc[i]);
            }
        }
        alpha_primal_soc = alpha_new.clamp(0.0, 1.0);

        #[allow(clippy::needless_range_loop)]
        let mut x_soc = vec![0.0; n];
        for i in 0..n {
            x_soc[i] = state.x[i] + alpha_primal_soc * dx_soc[i];
            if state.x_l[i].is_finite() {
                x_soc[i] = x_soc[i].max(state.x_l[i] + 1e-14);
            }
            if state.x_u[i].is_finite() {
                x_soc[i] = x_soc[i].min(state.x_u[i] - 1e-14);
            }
        }

        let mut obj_soc = f64::INFINITY;
        if !problem.objective(&x_soc, true, &mut obj_soc) || !obj_soc.is_finite() { return None; }
        let mut g_soc = vec![0.0; m];
        if !problem.constraints(&x_soc, false, &mut g_soc) { return None; }
        if g_soc.iter().any(|v| !v.is_finite()) { return None; }

        let theta_soc = convergence::primal_infeasibility(&g_soc, &state.g_l, &state.g_u);

        if theta_soc >= kappa_soc * theta_prev_soc {
            return None;
        }
        theta_prev_soc = theta_soc;

        let phi_soc = compute_barrier_phi(
            obj_soc, &x_soc, &g_soc, state, n, m, options.constraint_slack_barrier,
        );

        // Pass ORIGINAL alpha (alpha_primal_test), not alpha_primal_soc.
        // Matches Ipopt IpFilterLSAcceptor.cpp:629.
        let (acceptable, _) = filter.check_acceptability(
            theta_current,
            phi_current,
            theta_soc,
            phi_soc,
            grad_phi_step,
            alpha,
        );

        if acceptable {
            return Some((x_soc, obj_soc, g_soc, alpha_primal_soc));
        }

        // Update latest trial_c for next iter's accumulation.
        for i in 0..m {
            let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_equality || state.g_l[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_l[i];
            } else if state.g_u[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_u[i];
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

    let kappa_soc = 0.99;
    let tau = (1.0 - state.mu).max(options.tau_min);

    let mut c_soc = vec![0.0; m];
    let mut latest_trial_c = vec![0.0; m];
    for i in 0..m {
        let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_equality || state.g_l[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_l[i];
            latest_trial_c[i] = g_trial[i] - state.g_l[i];
        } else if state.g_u[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_u[i];
            latest_trial_c[i] = g_trial[i] - state.g_u[i];
        }
    }

    let mut alpha_primal_soc = alpha;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

    for _soc_iter in 0..options.max_soc {
        for i in 0..m {
            c_soc[i] += alpha_primal_soc * latest_trial_c[i];
        }

        let dx_soc = match kkt::solve_condensed_soc(condensed, solver, &c_soc) {
            Ok(d) => d,
            Err(_) => return None,
        };

        let mut alpha_new: f64 = 1.0;
        for i in 0..n {
            if state.x_l[i].is_finite() && dx_soc[i] < 0.0 {
                let slack = state.x[i] - state.x_l[i];
                alpha_new = alpha_new.min(-tau * slack / dx_soc[i]);
            }
            if state.x_u[i].is_finite() && dx_soc[i] > 0.0 {
                let slack = state.x_u[i] - state.x[i];
                alpha_new = alpha_new.min(tau * slack / dx_soc[i]);
            }
        }
        alpha_primal_soc = alpha_new.clamp(0.0, 1.0);

        #[allow(clippy::needless_range_loop)]
        let mut x_soc = vec![0.0; n];
        for i in 0..n {
            x_soc[i] = state.x[i] + alpha_primal_soc * dx_soc[i];
            if state.x_l[i].is_finite() {
                x_soc[i] = x_soc[i].max(state.x_l[i] + 1e-14);
            }
            if state.x_u[i].is_finite() {
                x_soc[i] = x_soc[i].min(state.x_u[i] - 1e-14);
            }
        }

        let mut obj_soc = f64::INFINITY;
        if !problem.objective(&x_soc, true, &mut obj_soc) || !obj_soc.is_finite() { return None; }
        let mut g_soc = vec![0.0; m];
        if !problem.constraints(&x_soc, false, &mut g_soc) { return None; }
        if g_soc.iter().any(|v| !v.is_finite()) { return None; }

        let theta_soc = convergence::primal_infeasibility(&g_soc, &state.g_l, &state.g_u);

        if theta_soc >= kappa_soc * theta_prev_soc {
            return None;
        }
        theta_prev_soc = theta_soc;

        let phi_soc = compute_barrier_phi(
            obj_soc, &x_soc, &g_soc, state, n, m, options.constraint_slack_barrier,
        );

        let (acceptable, _) = filter.check_acceptability(
            theta_current,
            phi_current,
            theta_soc,
            phi_soc,
            grad_phi_step,
            alpha,
        );

        if acceptable {
            return Some((x_soc, obj_soc, g_soc, alpha_primal_soc));
        }

        for i in 0..m {
            let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_equality || state.g_l[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_l[i];
            } else if state.g_u[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_u[i];
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

    let kappa_soc = 0.99;
    let tau = (1.0 - state.mu).max(options.tau_min);

    let mut c_soc = vec![0.0; m];
    let mut latest_trial_c = vec![0.0; m];
    for i in 0..m {
        let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
            && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
        if is_equality || state.g_l[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_l[i];
            latest_trial_c[i] = g_trial[i] - state.g_l[i];
        } else if state.g_u[i].is_finite() {
            c_soc[i] = state.g[i] - state.g_u[i];
            latest_trial_c[i] = g_trial[i] - state.g_u[i];
        }
    }

    let mut alpha_primal_soc = alpha;
    let mut theta_prev_soc = convergence::primal_infeasibility(g_trial, &state.g_l, &state.g_u);

    for _soc_iter in 0..options.max_soc {
        for i in 0..m {
            c_soc[i] += alpha_primal_soc * latest_trial_c[i];
        }

        let dx_soc = match kkt::solve_sparse_condensed_soc(condensed, solver, &c_soc) {
            Ok(d) => d,
            Err(_) => return None,
        };

        let mut alpha_new: f64 = 1.0;
        for i in 0..n {
            if state.x_l[i].is_finite() && dx_soc[i] < 0.0 {
                let slack = state.x[i] - state.x_l[i];
                alpha_new = alpha_new.min(-tau * slack / dx_soc[i]);
            }
            if state.x_u[i].is_finite() && dx_soc[i] > 0.0 {
                let slack = state.x_u[i] - state.x[i];
                alpha_new = alpha_new.min(tau * slack / dx_soc[i]);
            }
        }
        alpha_primal_soc = alpha_new.clamp(0.0, 1.0);

        let mut x_soc = vec![0.0; n];
        for i in 0..n {
            x_soc[i] = state.x[i] + alpha_primal_soc * dx_soc[i];
            if state.x_l[i].is_finite() { x_soc[i] = x_soc[i].max(state.x_l[i] + 1e-14); }
            if state.x_u[i].is_finite() { x_soc[i] = x_soc[i].min(state.x_u[i] - 1e-14); }
        }

        let mut obj_soc = f64::INFINITY;
        if !problem.objective(&x_soc, true, &mut obj_soc) || !obj_soc.is_finite() { return None; }
        let mut g_soc = vec![0.0; m];
        if !problem.constraints(&x_soc, false, &mut g_soc) { return None; }
        if g_soc.iter().any(|v| !v.is_finite()) { return None; }

        let theta_soc = convergence::primal_infeasibility(&g_soc, &state.g_l, &state.g_u);
        if theta_soc >= kappa_soc * theta_prev_soc { return None; }
        theta_prev_soc = theta_soc;

        let phi_soc = compute_barrier_phi(
            obj_soc, &x_soc, &g_soc, state, n, m, options.constraint_slack_barrier,
        );

        let (acceptable, _) = filter.check_acceptability(
            theta_current, phi_current, theta_soc, phi_soc, grad_phi_step, alpha,
        );
        if acceptable { return Some((x_soc, obj_soc, g_soc, alpha_primal_soc)); }

        for i in 0..m {
            let is_equality = state.g_l[i].is_finite() && state.g_u[i].is_finite()
                && (state.g_l[i] - state.g_u[i]).abs() < 1e-15;
            if is_equality || state.g_l[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_l[i];
            } else if state.g_u[i].is_finite() {
                latest_trial_c[i] = g_soc[i] - state.g_u[i];
            }
        }
    }

    None
}

/// Reset bound multipliers after restoration, matching Ipopt
/// IpRestoMinC_1Nrm.cpp:374-419:
///   1. Tentatively set z = mu/slack per bound (no element-wise clamp).
///   2. If max(|z|) exceeds bound_mult_reset_threshold (1e3),
///      *nuclear reset*: set ALL z_L, z_U to 1.0.
///   3. Otherwise keep z = mu/slack as-is.
/// An element-wise clamp at 1e3 leaves inf_du stuck at ~mu/slack - 1000
/// when slack is tight — the least-squares y computed after can't absorb
/// that. Returns whether the nuclear reset was triggered (for downstream
/// v_L/v_U handling).
fn reset_bound_multipliers_after_restoration(state: &mut SolverState, n: usize) -> bool {
    let bound_mult_reset_threshold = 1000.0;
    let mu_for_reset = state.mu;
    let mut z_max: f64 = 0.0;
    for i in 0..n {
        if state.x_l[i].is_finite() {
            let slack = (state.x[i] - state.x_l[i]).max(1e-12);
            state.z_l[i] = mu_for_reset / slack;
            z_max = z_max.max(state.z_l[i]);
        } else {
            state.z_l[i] = 0.0;
        }
        if state.x_u[i].is_finite() {
            let slack = (state.x_u[i] - state.x[i]).max(1e-12);
            state.z_u[i] = mu_for_reset / slack;
            z_max = z_max.max(state.z_u[i]);
        } else {
            state.z_u[i] = 0.0;
        }
    }
    let nuclear_reset = z_max > bound_mult_reset_threshold;
    if nuclear_reset {
        for i in 0..n {
            state.z_l[i] = if state.x_l[i].is_finite() { 1.0 } else { 0.0 };
            state.z_u[i] = if state.x_u[i].is_finite() { 1.0 } else { 0.0 };
        }
    }
    nuclear_reset
}

/// Recompute y at the restored point via the augmented-saddle-point
/// least-squares multiplier estimate, INCLUDING the reset z_L/z_U
/// contribution. Otherwise any deviation between z_true = mu/slack
/// (huge at tight slack) and the reset value (1.0) appears entirely
/// in inf_du, driving a bad first Newton step.
///
/// Uses the Ipopt-exact augmented saddle-point system
///   [ I   J^T ] [ r ] = [ grad_f - z_L + z_U ]
///   [ J    0  ] [ y ]   [ 0                   ]
/// (matches IpLeastSquareMults::CalculateMultipliers with W=0, δ=0).
/// This is far better conditioned than the normal equations
/// J*J^T*y = rhs when J is nearly rank-deficient (as happens on
/// AC-OPF with gauge freedom — case30_ieee hits this at a
/// post-restoration feasible iterate).
///
/// If the augmented solve fails or the result exceeds
/// `constr_mult_init_max` (matching Ipopt
/// DefaultIterateInitializer::least_square_mults), zero out y.
fn recompute_y_after_restoration(
    state: &mut SolverState,
    options: &SolverOptions,
    n: usize,
    m: usize,
) {
    if m == 0 {
        return;
    }
    let y_ls_result = compute_ls_multiplier_estimate_augmented(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.g_l, &state.g_u, n, m,
        Some(&state.z_l), Some(&state.z_u),
    );
    let y_accepted = match y_ls_result {
        Some(y_ls) => {
            let max_abs = y_ls.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if max_abs > options.constr_mult_init_max { None } else { Some(y_ls) }
        }
        None => None,
    };
    if let Some(y_ls) = y_accepted {
        state.y.copy_from_slice(&y_ls);
    } else {
        for i in 0..m {
            state.y[i] = 0.0;
        }
    }
}

/// Reset constraint-slack barrier multipliers v_L, v_U after restoration,
/// mirroring the nuclear-reset semantics of bound multipliers. Equality
/// constraints get v=0 (no slack barrier). For inequality constraints, if
/// the bound-multiplier reset triggered the all-to-1.0 path, also reset v
/// to 1.0; else use v = mu/slack uncapped.
fn reset_constraint_slack_multipliers_after_restoration(
    state: &mut SolverState,
    m: usize,
    nuclear_reset: bool,
) {
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
            state.v_l[i] = if nuclear_reset { 1.0 } else { mu_r / slack };
        } else {
            state.v_l[i] = 0.0;
        }
        if state.g_u[i].is_finite() {
            let slack = (state.g_u[i] - state.g[i]).max(1e-12);
            state.v_u[i] = if nuclear_reset { 1.0 } else { mu_r / slack };
        } else {
            state.v_u[i] = 0.0;
        }
    }
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

    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
    update_lbfgs_hessian(lbfgs_state, state);

    let nuclear_reset = reset_bound_multipliers_after_restoration(state, n);
    recompute_y_after_restoration(state, options, n, m);

    reset_constraint_slack_multipliers_after_restoration(state, m, nuclear_reset);

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

    // Reset mu_state mode (restoration is a fresh start)
    // In monotone strategy, stay in Fixed mode
    if options.mu_strategy_adaptive {
        mu_state.mode = MuMode::Free;
    }
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
/// Configure the inner SolverOptions for the restoration NLP solve. Caps
/// max_iter at restoration_max_iter (>=500), disables nested restoration to
/// prevent recursion, sets mu_init to resto_mu, disables stall detection
/// (restoration makes slow steady progress that would trip the 30-iter
/// stall limit), and relaxes tol to 1e-7 (we want feasibility, not
/// optimality). Propagates the remaining wall-time budget so the inner
/// solve can't outlive the outer fallback cascade; returns None when the
/// remaining budget is < 0.5s. Scales early_stall_timeout by restoration
/// NLP size so large restorations get the full default timeout.
fn configure_restoration_inner_options(
    options: &SolverOptions,
    resto_mu: f64,
    resto_dim: usize,
    start_time: Instant,
) -> Option<SolverOptions> {
    let mut inner_opts = options.clone();
    inner_opts.max_iter = options.restoration_max_iter.max(500);
    inner_opts.disable_nlp_restoration = true;
    inner_opts.print_level = if options.print_level >= 5 { 3 } else { 0 };
    inner_opts.mu_init = resto_mu;
    inner_opts.stall_iter_limit = 0;
    inner_opts.tol = 1e-7;

    if options.max_wall_time > 0.0 {
        let remaining = options.max_wall_time - start_time.elapsed().as_secs_f64();
        if remaining < 0.5 {
            return None;
        }
        inner_opts.max_wall_time = remaining;
    }
    inner_opts.early_stall_timeout = if options.early_stall_timeout > 0.0 {
        if resto_dim > 500 {
            options.early_stall_timeout
        } else {
            options.early_stall_timeout.min(3.0)
        }
    } else {
        3.0
    };
    Some(inner_opts)
}

fn attempt_nlp_restoration<P: NlpProblem>(
    problem: &P,
    state: &SolverState,
    filter: &Filter,
    options: &SolverOptions,
    theta_current: f64,
    start_time: Instant,
) -> (Vec<f64>, RestorationOutcome) {
    let n = state.n;
    let m = state.m;

    if options.print_level >= 5 {
        rip_log!(
            "ripopt: Entering NLP restoration (theta={:.2e}, mu={:.2e})",
            theta_current, state.mu
        );
    }

    // Ipopt's published default (IpRestoIpoptNLP.cpp:60): rho=1000 fixed.
    // The earlier adaptive rho was premature optimization — Ipopt's slack-penalty
    // analysis relies on rho being fixed and large enough that sum(p+n) is
    // monotone across inner iterations.
    let rho = 1000.0;

    // Ipopt's resto_mu = max(curr_mu, ||c(x_r)||_inf) (IpRestoIterateInitializer.cpp:58).
    // Using a mu_init consistent with the current infeasibility makes the
    // closed-form (p,n) init well-conditioned: when theta ≫ mu, the slacks
    // would otherwise be pinned near 0 with enormous bound multipliers.
    let c_inf = {
        let mut v: f64 = 0.0;
        for i in 0..m {
            if state.g_l[i].is_finite() && state.g[i] < state.g_l[i] {
                v = v.max(state.g_l[i] - state.g[i]);
            } else if state.g_u[i].is_finite() && state.g[i] > state.g_u[i] {
                v = v.max(state.g[i] - state.g_u[i]);
            }
        }
        v
    };
    let resto_mu = state.mu.max(c_inf);

    // Build restoration NLP using the same resto_mu for p/n quadratic init.
    let resto_nlp = RestorationNlp::new(problem, &state.x, resto_mu, rho, 1.0);

    let inner_opts = match configure_restoration_inner_options(
        options, resto_mu, resto_nlp.num_variables() + resto_nlp.num_constraints(), start_time,
    ) {
        Some(opts) => opts,
        None => return (state.x[..n].to_vec(), RestorationOutcome::Failed),
    };

    // Solve the restoration NLP
    let result = solve_ipm(&resto_nlp, &inner_opts);

    // Extract x_orig from the restoration solution
    let x_nlp: Vec<f64> = result.x[..n].to_vec();

    // Evaluate original constraints at the restored point
    let mut g_new = vec![0.0; m];
    if !problem.constraints(&x_nlp, true, &mut g_new)
        || g_new.iter().any(|v| !v.is_finite())
    {
        return (x_nlp, RestorationOutcome::Failed);
    }
    let theta_new = convergence::primal_infeasibility(&g_new, &state.g_l, &state.g_u);

    // Evaluate original objective at the restored point
    let mut phi_new = f64::INFINITY;
    if !problem.objective(&x_nlp, false, &mut phi_new) || !phi_new.is_finite() {
        return (x_nlp, RestorationOutcome::Failed);
    }

    if options.print_level >= 5 {
        // Log slack residuals to diagnose restoration effectiveness
        let sum_p: f64 = result.x[n..n + m].iter().sum();
        let sum_n: f64 = result.x[n + m..n + 2 * m].iter().sum();
        rip_log!(
            "ripopt: NLP restoration result: status={:?}, theta_new={:.2e} (was {:.2e}), phi_new={:.2e}, sum_p={:.2e}, sum_n={:.2e}, iters={}",
            result.status, theta_new, theta_current, phi_new, sum_p, sum_n, result.iterations
        );
    }

    let inner_converged = result.status == SolveStatus::Optimal;
    let outcome = classify_restoration_outcome(
        filter, options, theta_current, theta_new, phi_new, inner_converged,
    );
    (x_nlp, outcome)
}

/// Classify the outcome of a completed restoration solve. Decision tree:
/// 1. theta_new < constr_viol_tol → Success (achieved feasibility).
/// 2. theta_new ≤ 0.5*theta_current → Success (50% reduction, stricter than
///    Gauss-Newton's 10% to avoid marginal "success" that prevents recovery
///    mechanisms from engaging).
/// 3. theta_new < 0.9*theta_current AND filter-acceptable → Success.
/// 4. inner_converged but no feasibility improvement → LocalInfeasibility
///    (the restoration NLP itself reached a stationary point of the
///    L1-feasibility objective with positive residual).
/// 5. Otherwise → Failed.
fn classify_restoration_outcome(
    filter: &Filter,
    options: &SolverOptions,
    theta_current: f64,
    theta_new: f64,
    phi_new: f64,
    inner_converged: bool,
) -> RestorationOutcome {
    if theta_new < options.constr_viol_tol {
        return RestorationOutcome::Success;
    }
    if theta_new <= 0.5 * theta_current {
        return RestorationOutcome::Success;
    }
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
            return RestorationOutcome::Success;
        }
    }
    if inner_converged {
        return RestorationOutcome::LocalInfeasibility;
    }
    RestorationOutcome::Failed
}

/// Compute average complementarity for recomputing mu after restoration.
/// Quality function for barrier parameter selection.
///
/// Evaluates Q(mu) = dual_inf² + primal_inf² + compl_err(mu)² for log-spaced
/// candidate mu values and returns the minimizer. The first two terms are fixed
/// at the current iterate; only the complementarity term varies with mu.
///
/// This replaces the Loqo oracle (mu = avg_compl / kappa) with a global search
/// that can make aggressive mu decreases when the iterate is well-centered.
#[allow(dead_code)]
fn quality_function_mu(state: &SolverState, mu_lower: f64, mu_upper: f64, n_candidates: usize) -> f64 {
    if mu_upper <= mu_lower || n_candidates < 2 {
        return mu_upper;
    }

    let n = state.n;
    let pi = state.constraint_violation();
    let di = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let fixed_part = pi * pi + di * di;

    let log_min = mu_lower.max(1e-20).ln();
    let log_max = mu_upper.ln();

    let mut best_mu = mu_upper;
    let mut best_q = f64::INFINITY;

    for k in 0..n_candidates {
        let t = k as f64 / (n_candidates - 1) as f64;
        let mu_candidate = (log_min + t * (log_max - log_min)).exp();

        let ci = convergence::complementarity_error(
            &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, mu_candidate,
        );
        let q = fixed_part + ci * ci;

        if q < best_q {
            best_q = q;
            best_mu = mu_candidate;
        }
    }

    best_mu
}

/// Compute least-squares multiplier estimate: min ||grad_f + J^T y||^2.
/// Solves the normal equations (J J^T) y = -J grad_f.
/// Uses dense Bunch-Kaufman for small problems, sparse LDL^T for large ones.
/// Returns Some(y) if successful and all estimates are within threshold; None otherwise.
fn compute_ls_multiplier_estimate(
    grad_f: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    g_l: &[f64],
    g_u: &[f64],
    n: usize,
    m: usize,
    max_abs_threshold: f64,
) -> Option<Vec<f64>> {
    compute_ls_multiplier_estimate_with_z(
        grad_f, jac_rows, jac_cols, jac_vals, g_l, g_u, n, m, max_abs_threshold, None, None,
    )
}

/// Compute least-squares multiplier estimate via the Ipopt-exact augmented
/// saddle-point system (matches IpLeastSquareMults::CalculateMultipliers with
/// W=0, δ=0, no inertia correction).
///
/// Solves the (n+m)×(n+m) system
///   [ I    J^T ] [ r ]   [ grad_f - z_L + z_U ]
///   [ J     0  ] [ y ] = [ 0                  ]
/// and returns the y block.
///
/// This is algebraically equivalent to the normal-equations form
/// (J·J^T)·y = J·(grad_f − z_L + z_U), but is numerically far better
/// conditioned when J is nearly rank-deficient. The normal-equations form
/// fails outright (J·J^T singular) for AC-OPF problems with gauge symmetry
/// at post-restoration iterates (case30_ieee), whereas the augmented form's
/// LDL^T handles the indefinite saddle point via Bunch-Kaufman pivoting
/// without explicitly forming J·J^T.
///
/// On factorization failure, returns None. Matches Ipopt behavior at
/// IpLeastSquareMults.cpp:82-87 — no δ_c retry, no Tikhonov regularization.
/// Callers apply the `constr_mult_init_max` cap externally.
fn compute_ls_multiplier_estimate_augmented(
    grad_f: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    g_l: &[f64],
    g_u: &[f64],
    n: usize,
    m: usize,
    z_l: Option<&[f64]>,
    z_u: Option<&[f64]>,
) -> Option<Vec<f64>> {
    use crate::linear_solver::SparseSymmetricMatrix;
    if m == 0 {
        return None;
    }

    // RHS: [grad_f − z_L + z_U; 0]
    let mut rhs = vec![0.0_f64; n + m];
    for i in 0..n {
        rhs[i] = grad_f[i];
        if let Some(zl) = z_l { rhs[i] -= zl[i]; }
        if let Some(zu) = z_u { rhs[i] += zu[i]; }
    }

    // Build augmented matrix in upper-triangle triplet form (row ≤ col):
    //   (i, i) = 1         for i in 0..n           (identity (1,1) block)
    //   (col, n+row) = J_{row,col}                 (J^T off-diagonal)
    //   (n+j, n+j) = 0     for j in 0..m           (structural zero on (2,2) diag)
    let nnz_est = n + jac_rows.len() + m;
    let mut ssm = SparseSymmetricMatrix {
        n: n + m,
        triplet_rows: Vec::with_capacity(nnz_est),
        triplet_cols: Vec::with_capacity(nnz_est),
        triplet_vals: Vec::with_capacity(nnz_est),
    };
    for i in 0..n {
        ssm.triplet_rows.push(i);
        ssm.triplet_cols.push(i);
        ssm.triplet_vals.push(1.0);
    }
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        // J^T at (col, n+row). Since col < n ≤ n+row, this is always upper triangle.
        ssm.triplet_rows.push(col);
        ssm.triplet_cols.push(n + row);
        ssm.triplet_vals.push(jac_vals[idx]);
    }
    for j in 0..m {
        ssm.triplet_rows.push(n + j);
        ssm.triplet_cols.push(n + j);
        ssm.triplet_vals.push(0.0);
    }

    let matrix = KktMatrix::Sparse(ssm);
    let mut solver = new_sparse_solver();
    if solver.factor(&matrix).is_err() {
        return None;
    }
    let mut sol = vec![0.0_f64; n + m];
    if solver.solve(&rhs, &mut sol).is_err() {
        return None;
    }
    if sol.iter().any(|v| !v.is_finite()) {
        return None;
    }

    let mut y_ls: Vec<f64> = sol[n..].to_vec();
    fix_inequality_mult_signs(&mut y_ls, g_l, g_u, m);
    Some(y_ls)
}

/// LS multiplier estimate with optional bound-multiplier contributions.
///
/// Minimizes ||grad_f + J^T y - z_L + z_U||², yielding the normal equation
///   (J J^T) y = -J * (grad_f - z_L + z_U)
/// If `z_l`/`z_u` are `None`, they are treated as zero (cold-start / pre-z
/// initialization). Post-restoration callers must pass them so the reset z
/// values are absorbed into y, otherwise inf_du = ||grad_f + J^T y - z_L + z_U||
/// stays large and the next Newton step is ill-directed.
fn compute_ls_multiplier_estimate_with_z(
    grad_f: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    g_l: &[f64],
    g_u: &[f64],
    n: usize,
    m: usize,
    max_abs_threshold: f64,
    z_l: Option<&[f64]>,
    z_u: Option<&[f64]>,
) -> Option<Vec<f64>> {
    // For large problems, use sparse J*J^T factorization
    if m > 500 {
        return compute_ls_multiplier_estimate_sparse(
            grad_f, jac_rows, jac_cols, jac_vals, g_l, g_u, n, m, max_abs_threshold, None,
            z_l, z_u,
        );
    }
    if m == 0 {
        return None;
    }

    // b = -J * (grad_f - z_L + z_U)
    let mut rhs_grad = grad_f.to_vec();
    if let Some(zl) = z_l {
        for i in 0..n { rhs_grad[i] -= zl[i]; }
    }
    if let Some(zu) = z_u {
        for i in 0..n { rhs_grad[i] += zu[i]; }
    }
    let mut b = vec![0.0; m];
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        b[row] -= jac_vals[idx] * rhs_grad[col];
    }

    // Compute A = J * J^T  (m x m dense symmetric matrix)
    let mut j_dense = vec![0.0; m * n];
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        j_dense[row * n + col] = jac_vals[idx];
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

    if !solved {
        return None;
    }

    let max_abs = y_ls.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    if max_abs > max_abs_threshold {
        return None;
    }

    fix_inequality_mult_signs(&mut y_ls, g_l, g_u, m);
    Some(y_ls)
}

/// Sparse variant of LS multiplier estimate for large problems.
/// Builds sparse J*J^T in COO format and factors with the sparse solver.
/// If `solver` is Some, reuses the solver (for recalc_y caching).
fn compute_ls_multiplier_estimate_sparse(
    grad_f: &[f64],
    jac_rows: &[usize],
    jac_cols: &[usize],
    jac_vals: &[f64],
    g_l: &[f64],
    g_u: &[f64],
    n: usize,
    m: usize,
    max_abs_threshold: f64,
    solver: Option<&mut Box<dyn LinearSolver>>,
    z_l: Option<&[f64]>,
    z_u: Option<&[f64]>,
) -> Option<Vec<f64>> {
    use crate::linear_solver::SparseSymmetricMatrix;
    if m == 0 {
        return None;
    }

    // b = -J * (grad_f - z_L + z_U)
    let mut rhs_grad = grad_f.to_vec();
    if let Some(zl) = z_l {
        for i in 0..n { rhs_grad[i] -= zl[i]; }
    }
    if let Some(zu) = z_u {
        for i in 0..n { rhs_grad[i] += zu[i]; }
    }
    let mut b = vec![0.0; m];
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        b[row] -= jac_vals[idx] * rhs_grad[col];
    }

    // Build sparse J*J^T: group Jacobian entries by column, then for each
    // column accumulate outer products into COO triplets.
    // Use a HashMap to accumulate duplicates efficiently.
    let mut col_entries: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
    for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
        col_entries[col].push((row, jac_vals[idx]));
    }

    // Estimate nnz: for PDE Jacobians, each column has O(1) nonzeros
    let total_col_nnz: usize = col_entries.iter().map(|c| c.len()).sum();
    let nnz_est = total_col_nnz * 4; // rough upper bound

    use std::collections::HashMap;
    let mut triplet_map: HashMap<(usize, usize), f64> = HashMap::with_capacity(nnz_est);

    for k in 0..n {
        let entries = &col_entries[k];
        for &(i, vi) in entries.iter() {
            for &(j, vj) in entries.iter() {
                if i >= j {
                    // Upper triangle: row <= col, so store (j, i) since j <= i
                    *triplet_map.entry((j, i)).or_insert(0.0) += vi * vj;
                }
            }
        }
    }

    // Add small regularization to diagonal for numerical stability
    let reg = 1e-12;
    for i in 0..m {
        *triplet_map.entry((i, i)).or_insert(0.0) += reg;
    }

    // Build SparseSymmetricMatrix (upper triangle: row <= col)
    let nnz = triplet_map.len();
    let mut ssm = SparseSymmetricMatrix {
        n: m,
        triplet_rows: Vec::with_capacity(nnz),
        triplet_cols: Vec::with_capacity(nnz),
        triplet_vals: Vec::with_capacity(nnz),
    };
    for (&(r, c), &v) in &triplet_map {
        ssm.triplet_rows.push(r);
        ssm.triplet_cols.push(c);
        ssm.triplet_vals.push(v);
    }

    let matrix = KktMatrix::Sparse(ssm);

    // Factor and solve
    let mut y_ls = vec![0.0; m];
    let solved = if let Some(ls) = solver {
        ls.factor(&matrix).is_ok() && ls.solve(&b, &mut y_ls).is_ok()
    } else {
        let mut ls = new_sparse_solver();
        ls.factor(&matrix).is_ok() && ls.solve(&b, &mut y_ls).is_ok()
    };

    if !solved {
        return None;
    }

    let max_abs = y_ls.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    if max_abs > max_abs_threshold {
        return None;
    }

    fix_inequality_mult_signs(&mut y_ls, g_l, g_u, m);
    Some(y_ls)
}

/// Fix signs of inequality constraint multipliers from LS estimate.
/// Ipopt convention (L = f + y^T g): g >= g_l → y >= 0, g <= g_u → y <= 0.
fn fix_inequality_mult_signs(y_ls: &mut [f64], g_l: &[f64], g_u: &[f64], m: usize) {
    for i in 0..m {
        if convergence::is_equality_constraint(g_l[i], g_u[i]) {
            continue;
        }
        let has_lower = g_l[i].is_finite();
        let has_upper = g_u[i].is_finite();
        if has_lower && !has_upper && y_ls[i] < 0.0 {
            y_ls[i] = 0.0;
        } else if has_upper && !has_lower && y_ls[i] > 0.0 {
            y_ls[i] = 0.0;
        } else if !has_lower && !has_upper {
            y_ls[i] = 0.0;
        }
    }
}

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

    // Dual infeasibility safeguard: prevent the barrier subproblem from being
    // declared "solved" when the NLP dual infeasibility is still large.
    // This prevents mu from collapsing to 1e-11 while du remains huge (issue #8 Class 2).
    let unscaled_du = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let du_floor = unscaled_du * 0.1; // 10% of unscaled du as floor on barrier error

    dual_err.max(compl_err).max(primal_err).max(du_floor)
}

/// Strategy 3: Try active set identification + reduced KKT solve.
///
/// At near-optimal points, identify variables at their bounds (active set),
/// fix them, solve the reduced KKT system for free variables, and check
/// if the result meets strict convergence tolerances.
/// Identification of the working active set for `try_active_set_solve`.
struct ActiveSet {
    active_lower: Vec<bool>,
    active_upper: Vec<bool>,
    /// `free_idx[k]` is the original variable index of the k-th free variable.
    free_idx: Vec<usize>,
    /// `orig_to_free[i]` is the reduced-system index of variable `i`, or
    /// `usize::MAX` when `i` is fixed at a bound.
    orig_to_free: Vec<usize>,
    n_free: usize,
    /// Reduced KKT dimension: `n_free + m`.
    dim: usize,
}

/// Identify the working active set: a variable is "active at lower bound"
/// if x_i is close to x_l_i (relative tol 1e-6) and z_l_i > 1e-8. Returns
/// `None` when the reduced system is too large for dense solve (dim > 500),
/// when no bounds are active (dim == n + m), or when dim == 0.
fn identify_active_bounds(state: &SolverState, n: usize, m: usize) -> Option<ActiveSet> {
    let tol_bound = 1e-6;
    let mut is_free = vec![true; n];
    let mut active_lower = vec![false; n];
    let mut active_upper = vec![false; n];
    let mut n_free = 0usize;

    for i in 0..n {
        let at_lower = state.x_l[i].is_finite()
            && (state.x[i] - state.x_l[i]).abs() < tol_bound * (1.0 + state.x_l[i].abs());
        let at_upper = state.x_u[i].is_finite()
            && (state.x_u[i] - state.x[i]).abs() < tol_bound * (1.0 + state.x_u[i].abs());

        if at_lower && state.z_l[i] > 1e-8 {
            is_free[i] = false;
            active_lower[i] = true;
        } else if at_upper && state.z_u[i] > 1e-8 {
            is_free[i] = false;
            active_upper[i] = true;
        } else {
            n_free += 1;
        }
    }

    if n_free == n {
        return None;
    }
    let dim = n_free + m;
    if dim > 500 || dim == 0 {
        return None;
    }

    let mut free_idx = Vec::with_capacity(n_free);
    let mut orig_to_free = vec![usize::MAX; n];
    for i in 0..n {
        if is_free[i] {
            orig_to_free[i] = free_idx.len();
            free_idx.push(i);
        }
    }
    Some(ActiveSet { active_lower, active_upper, free_idx, orig_to_free, n_free, dim })
}

/// Build the reduced dense KKT system for the working active set:
///   [ H_ff   J_f^T ] [ dx_f ]   [ -grad_f_f       ]
///   [ J_f    0     ] [ dy   ] = [ g_target - g(x) ]
/// where H_ff is the Hessian restricted to free-free pairs, J_f is the
/// Jacobian columns for free variables, and the constraint target picks
/// either g_l or g_u depending on which side is active (or 0 for inactive
/// inequalities). Returns the dense `dim*dim` KKT (row-major) and the
/// length-`dim` RHS.
fn build_reduced_kkt_dense(
    state: &SolverState,
    orig_to_free: &[usize],
    free_idx: &[usize],
    n_free: usize,
    m: usize,
    dim: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut kkt = vec![0.0; dim * dim];
    let mut rhs = vec![0.0; dim];

    // H_ff (top-left n_free x n_free)
    for (idx, (&row, &col)) in state.hess_rows.iter().zip(state.hess_cols.iter()).enumerate() {
        let fr = orig_to_free[row];
        let fc = orig_to_free[col];
        if fr != usize::MAX && fc != usize::MAX {
            kkt[fr * dim + fc] += state.hess_vals[idx];
            if fr != fc {
                kkt[fc * dim + fr] += state.hess_vals[idx];
            }
        }
    }

    // J_f (bottom-left m x n_free) and J_f^T (top-right n_free x m)
    for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        let fc = orig_to_free[col];
        if fc != usize::MAX {
            let r = n_free + row;
            kkt[r * dim + fc] += state.jac_vals[idx];
            kkt[fc * dim + r] += state.jac_vals[idx];
        }
    }

    // RHS top: -grad_f for free variables
    for k in 0..n_free {
        rhs[k] = -state.grad_f[free_idx[k]];
    }

    // RHS bottom: g_target - g, picking the active side per constraint
    for i in 0..m {
        if (state.g_l[i] - state.g_u[i]).abs() < 1e-20 {
            rhs[n_free + i] = state.g_l[i] - state.g[i];
        } else if state.g[i] <= state.g_l[i] + 1e-10 {
            rhs[n_free + i] = state.g_l[i] - state.g[i];
        } else if state.g[i] >= state.g_u[i] - 1e-10 {
            rhs[n_free + i] = state.g_u[i] - state.g[i];
        } else {
            rhs[n_free + i] = 0.0;
        }
    }

    (kkt, rhs)
}

fn try_active_set_solve<P: NlpProblem>(
    state: &mut SolverState,
    problem: &P,
    options: &SolverOptions,
    linear_constraints: Option<&[bool]>,
    lbfgs_mode: bool,
) -> Option<SolveResult> {
    let n = state.n;
    let m = state.m;

    let active_set = match identify_active_bounds(state, n, m) {
        Some(set) => set,
        None => return None,
    };
    let ActiveSet { active_lower, active_upper, free_idx, orig_to_free, n_free, dim } = active_set;

    // Fix active variables at their bounds (save full state for restoration)
    let saved_x = state.x.clone();
    let saved_y = state.y.clone();
    let saved_zl = state.z_l.clone();
    let saved_zu = state.z_u.clone();
    for i in 0..n {
        if active_lower[i] {
            state.x[i] = state.x_l[i];
        } else if active_upper[i] {
            state.x[i] = state.x_u[i];
        }
    }

    // Re-evaluate at the snapped point
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);

    let (mut kkt, mut rhs) = build_reduced_kkt_dense(
        state, &orig_to_free, &free_idx, n_free, m, dim,
    );
    let solution = dense_symmetric_solve(dim, &mut kkt, &mut rhs);
    if solution.is_none() {
        // Singular system, restore and bail
        state.x.copy_from_slice(&saved_x);
        state.y.copy_from_slice(&saved_y);
        state.z_l.copy_from_slice(&saved_zl);
        state.z_u.copy_from_slice(&saved_zu);
        let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);
        return None;
    }
    let sol = solution.unwrap();

    // Apply the step for free variables (with damping for safety)
    let alpha = 1.0; // full Newton step
    for k in 0..n_free {
        let i = free_idx[k];
        state.x[i] += alpha * sol[k];
        // Clamp to bounds
        if state.x_l[i].is_finite() {
            state.x[i] = state.x[i].max(state.x_l[i]);
        }
        if state.x_u[i].is_finite() {
            state.x[i] = state.x[i].min(state.x_u[i]);
        }
    }

    // Update y from the solve
    for i in 0..m {
        state.y[i] = sol[n_free + i];
    }

    // Re-evaluate at the new point
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);

    // Recover z from stationarity: ∇f + J^T y = z_l - z_u
    let mut grad_jty = state.grad_f.clone();
    for (idx, (&row, &col)) in state.jac_rows.iter().zip(state.jac_cols.iter()).enumerate() {
        grad_jty[col] += state.jac_vals[idx] * state.y[row];
    }
    for i in 0..n {
        state.z_l[i] = 0.0;
        state.z_u[i] = 0.0;
        if state.x_l[i].is_finite() && grad_jty[i] > 0.0 {
            state.z_l[i] = grad_jty[i];
        } else if state.x_u[i].is_finite() && grad_jty[i] < 0.0 {
            state.z_u[i] = -grad_jty[i];
        }
    }

    // Check strict convergence (use max-norm for primal infeasibility)
    let primal_inf = convergence::primal_infeasibility_max(&state.g, &state.g_l, &state.g_u);
    let dual_inf = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    let compl_inf = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );
    let multiplier_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
        + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
    let multiplier_count = m + 2 * n;

    let conv_info = ConvergenceInfo {
        primal_inf,
        dual_inf,
        dual_inf_unscaled: dual_inf, // same z used for both
        compl_inf,
        mu: 0.0, // at the solution, mu should be zero
        objective: state.obj,
        multiplier_sum,
        multiplier_count,
    };

    if let ConvergenceStatus::Converged = check_convergence(&conv_info, options, 0) {
        return Some(make_result(state, SolveStatus::Optimal));
    }

    // Didn't converge; restore original state and re-evaluate
    state.x.copy_from_slice(&saved_x);
    state.y.copy_from_slice(&saved_y);
    state.z_l.copy_from_slice(&saved_zl);
    state.z_u.copy_from_slice(&saved_zu);
    let _ = state.evaluate_with_linear(problem, 1.0, linear_constraints, lbfgs_mode);

    None
}

/// Dense symmetric indefinite solve using diagonal pivoting (LDL^T).
/// Solves A*x = b in-place. Returns Some(x) on success, None if singular.
/// `a` is row-major dim x dim, `b` is length dim.
fn dense_symmetric_solve(dim: usize, a: &mut [f64], b: &mut [f64]) -> Option<Vec<f64>> {
    // Gaussian elimination with partial pivoting for symmetric indefinite systems.
    // Simple but sufficient for small systems (dim < 500).
    let mut piv = vec![0usize; dim];
    for i in 0..dim {
        piv[i] = i;
    }

    for k in 0..dim {
        // Find pivot: largest diagonal element in remaining submatrix
        let mut max_val = a[k * dim + k].abs();
        let mut max_idx = k;
        for i in (k + 1)..dim {
            if a[i * dim + i].abs() > max_val {
                max_val = a[i * dim + i].abs();
                max_idx = i;
            }
        }

        if max_val < 1e-15 {
            return None; // Singular
        }

        // Swap rows/cols k and max_idx
        if max_idx != k {
            piv.swap(k, max_idx);
            // Swap rows
            for j in 0..dim {
                let tmp = a[k * dim + j];
                a[k * dim + j] = a[max_idx * dim + j];
                a[max_idx * dim + j] = tmp;
            }
            // Swap cols
            for i in 0..dim {
                let tmp = a[i * dim + k];
                a[i * dim + k] = a[i * dim + max_idx];
                a[i * dim + max_idx] = tmp;
            }
            b.swap(k, max_idx);
        }

        let pivot = a[k * dim + k];
        // Eliminate below
        for i in (k + 1)..dim {
            let factor = a[i * dim + k] / pivot;
            a[i * dim + k] = factor;
            for j in (k + 1)..dim {
                a[i * dim + j] -= factor * a[k * dim + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    let mut x = b.to_vec();
    for k in (0..dim).rev() {
        for j in (k + 1)..dim {
            x[k] -= a[k * dim + j] * x[j];
        }
        x[k] /= a[k * dim + k];
    }

    Some(x)
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
    diag.final_primal_inf = convergence::primal_infeasibility_max(&state.g, &state.g_l, &state.g_u);
    diag.final_dual_inf = convergence::dual_infeasibility(
        &state.grad_f, &state.jac_rows, &state.jac_cols, &state.jac_vals,
        &state.y, &state.z_l, &state.z_u, n,
    );
    diag.final_compl = convergence::complementarity_error(
        &state.x, &state.x_l, &state.x_u, &state.z_l, &state.z_u, 0.0,
    );

    // Compute dual scaling factor s_d (same formula as check_convergence)
    {
        let s_max: f64 = 100.0;
        let s_d_max: f64 = 1e4;
        let mult_sum: f64 = state.y.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_l.iter().map(|v| v.abs()).sum::<f64>()
            + state.z_u.iter().map(|v| v.abs()).sum::<f64>();
        let mult_count = m + 2 * n;
        diag.final_s_d = if mult_count > 0 {
            ((s_max.max(mult_sum / mult_count as f64)) / s_max).min(s_d_max)
        } else {
            1.0
        };
    }

    // Use iterative z for reported bound multipliers (Ipopt semantics).
    // Unscale: z_unscaled = z_scaled / obj_scaling.
    let mut z_l_out = vec![0.0; n];
    let mut z_u_out = vec![0.0; n];
    for i in 0..n {
        z_l_out[i] = state.z_l[i] / state.obj_scaling;
        z_u_out[i] = state.z_u[i] / state.obj_scaling;
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

    // No Acceptable status anymore — pass status through directly.
    let validated_status = status;

    SolveResult {
        x: state.x.clone(),
        objective: state.obj / state.obj_scaling,
        constraint_multipliers: y_out,
        bound_multipliers_lower: z_l_out,
        bound_multipliers_upper: z_u_out,
        constraint_values: g_out,
        status: validated_status,
        iterations: state.iter,
        diagnostics: diag,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a minimal SolverState for testing private mu/complementarity helpers.
    // The caller supplies only fields the test exercises; everything else is zeroed.
    fn minimal_state(n: usize, m: usize) -> SolverState {
        SolverState {
            x: vec![0.0; n],
            y: vec![0.0; m],
            z_l: vec![0.0; n],
            z_u: vec![0.0; n],
            v_l: vec![0.0; m],
            v_u: vec![0.0; m],
            dx: vec![0.0; n],
            dy: vec![0.0; m],
            dz_l: vec![0.0; n],
            dz_u: vec![0.0; n],
            mu: 0.1,
            alpha_primal: 0.0,
            alpha_dual: 0.0,
            iter: 0,
            x_l: vec![f64::NEG_INFINITY; n],
            x_u: vec![f64::INFINITY; n],
            g_l: vec![f64::NEG_INFINITY; m],
            g_u: vec![f64::INFINITY; m],
            n,
            m,
            obj: 0.0,
            grad_f: vec![0.0; n],
            g: vec![0.0; m],
            jac_rows: Vec::new(),
            jac_cols: Vec::new(),
            jac_vals: Vec::new(),
            hess_rows: Vec::new(),
            hess_cols: Vec::new(),
            hess_vals: Vec::new(),
            consecutive_acceptable: 0,
            obj_scaling: 1.0,
            g_scaling: vec![1.0; m],
            diagnostics: SolverDiagnostics::default(),
            x_last_eval: vec![f64::NAN; n],
        }
    }

    #[test]
    fn test_avg_compl_variable_bounds_only() {
        // 2 vars, lower-bound only. x = [1.5, 2.0], x_l = [1.0, 1.0], z_l = [2.0, 3.0].
        // slacks = [0.5, 1.0]; avg_compl = (0.5*2.0 + 1.0*3.0) / 2 = 4.0 / 2 = 2.0
        let mut state = minimal_state(2, 0);
        state.x = vec![1.5, 2.0];
        state.x_l = vec![1.0, 1.0];
        state.z_l = vec![2.0, 3.0];
        let avg = compute_avg_complementarity(&state);
        assert!((avg - 2.0).abs() < 1e-12, "expected 2.0, got {}", avg);
    }

    #[test]
    fn test_avg_compl_both_bounds() {
        // 1 var, both bounds. x = 1.5, x_l = 1.0, x_u = 2.0, z_l = 2.0, z_u = 3.0.
        // avg = (0.5*2.0 + 0.5*3.0) / 2 = 2.5 / 2 = 1.25
        let mut state = minimal_state(1, 0);
        state.x = vec![1.5];
        state.x_l = vec![1.0];
        state.x_u = vec![2.0];
        state.z_l = vec![2.0];
        state.z_u = vec![3.0];
        let avg = compute_avg_complementarity(&state);
        assert!((avg - 1.25).abs() < 1e-12, "expected 1.25, got {}", avg);
    }

    #[test]
    fn test_avg_compl_inequality_fallback() {
        // No variable bounds, but an inequality constraint with v_l > 0 triggers fallback.
        // g = 2.0, g_l = 1.0, v_l = 0.5 -> slack = 1.0, contrib = 0.5; avg = 0.5 / 1 = 0.5.
        let mut state = minimal_state(1, 1);
        state.g = vec![2.0];
        state.g_l = vec![1.0];
        state.g_u = vec![f64::INFINITY];
        state.v_l = vec![0.5];
        let avg = compute_avg_complementarity(&state);
        assert!((avg - 0.5).abs() < 1e-12, "fallback path: expected 0.5, got {}", avg);
    }

    #[test]
    fn test_avg_compl_inequality_fallback_skipped_when_bounds_exist() {
        // Variable bounds present AND inequality constraint with v_l > 0: fallback is skipped.
        // Only the variable bound contributes, so avg_compl = slack * z_l / 1 = 1.0 * 1.0 = 1.0.
        let mut state = minimal_state(1, 1);
        state.x = vec![2.0];
        state.x_l = vec![1.0];
        state.z_l = vec![1.0];
        state.g = vec![2.0];
        state.g_l = vec![1.0];
        state.v_l = vec![99.0]; // Would bias avg if fallback ran
        let avg = compute_avg_complementarity(&state);
        assert!((avg - 1.0).abs() < 1e-12,
            "fallback must skip when var bounds present: got {}", avg);
    }

    #[test]
    fn test_avg_compl_no_bounds_anywhere() {
        // Unconstrained, no bounds: avg_compl = 0.0 (count stays at 0).
        let state = minimal_state(3, 0);
        let avg = compute_avg_complementarity(&state);
        assert_eq!(avg, 0.0);
    }

    #[test]
    fn test_quality_function_mu_degenerate_range() {
        // mu_upper <= mu_lower → returns mu_upper unchanged
        let state = minimal_state(1, 0);
        let mu = quality_function_mu(&state, 1.0, 1.0, 5);
        assert_eq!(mu, 1.0);
        let mu2 = quality_function_mu(&state, 2.0, 1.0, 5);
        assert_eq!(mu2, 1.0, "lower > upper still returns upper");
    }

    #[test]
    fn test_quality_function_mu_too_few_candidates() {
        // n_candidates < 2 → returns mu_upper
        let state = minimal_state(1, 0);
        let mu = quality_function_mu(&state, 1e-6, 1e-3, 1);
        assert_eq!(mu, 1e-3);
        let mu0 = quality_function_mu(&state, 1e-6, 1e-3, 0);
        assert_eq!(mu0, 1e-3);
    }

    #[test]
    fn test_quality_function_mu_picks_candidate_in_range() {
        // Well-posed state: 1 lower-bound-active variable. The quality function
        // q(mu) = pi^2 + di^2 + ci(mu)^2 where ci is complementarity_error.
        // With pi=di=0 (satisfied primal/dual) and slack*z = 0.5*0.2 = 0.1,
        // ci(mu) is minimized at mu ≈ 0.1 (interior of [1e-6, 1e-1]).
        let mut state = minimal_state(1, 0);
        state.x = vec![1.5];
        state.x_l = vec![1.0];
        state.z_l = vec![0.2];
        let mu = quality_function_mu(&state, 1e-6, 1e-1, 11);
        // Candidate grid is log-spaced over [1e-6, 1e-1]; ci(mu=0.1)=0 exactly.
        // exp(ln(1e-1)) has a ~1e-16 rounding overshoot, so allow a slack margin.
        assert!(mu >= 1e-6 * (1.0 - 1e-12) && mu <= 1e-1 * (1.0 + 1e-12),
            "mu must lie in range, got {}", mu);
        // The exact optimum 0.1 is grid point k=10 with n_candidates=11.
        assert!((mu - 0.1).abs() < 1e-10, "expected mu≈0.1, got {}", mu);
    }
}
