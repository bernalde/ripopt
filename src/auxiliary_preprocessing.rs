//! Internal auxiliary-system preprocessing utilities.
//!
//! This module is intentionally crate-private. The auxiliary preprocessor is an
//! implementation detail of `enable_preprocessing`; it must not expose a public
//! decomposition or transform API.

use crate::logging::rip_log;
use crate::options::SolverOptions;
use crate::problem::NlpProblem;
use crate::result::{SolveResult, SolveStatus};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::time::Instant;

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EqualityBlock {
    pub(crate) rows: Vec<usize>,
    pub(crate) vars: Vec<usize>,
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PresolveCandidate {
    pub(crate) blocks: Vec<EqualityBlock>,
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AuxiliarySolveOutcome {
    pub(crate) x: Vec<f64>,
    pub(crate) blocks_solved: usize,
    pub(crate) max_residual: f64,
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AuxiliarySolveError {
    InvalidBlock {
        block: EqualityBlock,
        reason: &'static str,
    },
    EvaluationFailed {
        block: EqualityBlock,
    },
    TimeBudgetExceeded {
        blocks_solved: usize,
    },
    BlockSolveFailed {
        block: EqualityBlock,
        status: SolveStatus,
        residual: f64,
    },
    RankDeficientBlock {
        block: EqualityBlock,
        rank: usize,
        expected: usize,
    },
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum BlockTriangularizationError {
    NonSquare {
        rows: usize,
        vars: usize,
    },
    ImperfectMatching {
        unmatched_rows: Vec<usize>,
        unmatched_vars: Vec<usize>,
    },
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct EqualityIncidence {
    pub(crate) n_vars: usize,
    pub(crate) m_orig: usize,
    pub(crate) row_global: Vec<usize>,
    pub(crate) row_local_for_global: Vec<Option<usize>>,
    pub(crate) row_adj_vars: Vec<Vec<usize>>,
    pub(crate) var_adj_rows: Vec<Vec<usize>>,
}

#[allow(dead_code)]
pub(crate) struct AuxiliaryBlockProblem<'a> {
    inner: &'a dyn NlpProblem,
    n_orig: usize,
    m_orig: usize,
    rows: Vec<usize>,
    vars: Vec<usize>,
    fixed_x: Vec<f64>,
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    jac_entry_map: Vec<usize>,
    inner_jac_nnz: usize,
    hess_rows: Vec<usize>,
    hess_cols: Vec<usize>,
    hess_entry_map: Vec<usize>,
    inner_hess_nnz: usize,
}

impl<'a> AuxiliaryBlockProblem<'a> {
    pub(crate) fn new(
        inner: &'a dyn NlpProblem,
        block: &EqualityBlock,
        fixed_x: &[f64],
    ) -> Result<Self, AuxiliarySolveError> {
        let n_orig = inner.num_variables();
        let m_orig = inner.num_constraints();
        validate_auxiliary_block(block, fixed_x, n_orig, m_orig)?;

        let mut row_local_for_global = vec![None; m_orig];
        for (local, &row) in block.rows.iter().enumerate() {
            row_local_for_global[row] = Some(local);
        }
        let mut var_local_for_global = vec![None; n_orig];
        for (local, &var) in block.vars.iter().enumerate() {
            var_local_for_global[var] = Some(local);
        }

        let (inner_jac_rows, inner_jac_cols) = inner.jacobian_structure();
        let mut jac_rows = Vec::new();
        let mut jac_cols = Vec::new();
        let mut jac_entry_map = Vec::new();
        for (idx, (&row, &col)) in inner_jac_rows.iter().zip(inner_jac_cols.iter()).enumerate() {
            if row >= m_orig || col >= n_orig {
                continue;
            }
            if let (Some(local_row), Some(local_col)) =
                (row_local_for_global[row], var_local_for_global[col])
            {
                jac_rows.push(local_row);
                jac_cols.push(local_col);
                jac_entry_map.push(idx);
            }
        }

        let (inner_hess_rows, inner_hess_cols) = inner.hessian_structure();
        let mut hess_rows = Vec::new();
        let mut hess_cols = Vec::new();
        let mut hess_entry_map = Vec::new();
        for (idx, (&row, &col)) in inner_hess_rows
            .iter()
            .zip(inner_hess_cols.iter())
            .enumerate()
        {
            if row >= n_orig || col >= n_orig {
                continue;
            }
            if let (Some(local_row), Some(local_col)) =
                (var_local_for_global[row], var_local_for_global[col])
            {
                if local_row >= local_col {
                    hess_rows.push(local_row);
                    hess_cols.push(local_col);
                } else {
                    hess_rows.push(local_col);
                    hess_cols.push(local_row);
                }
                hess_entry_map.push(idx);
            }
        }

        Ok(Self {
            inner,
            n_orig,
            m_orig,
            rows: block.rows.clone(),
            vars: block.vars.clone(),
            fixed_x: fixed_x.to_vec(),
            jac_rows,
            jac_cols,
            jac_entry_map,
            inner_jac_nnz: inner_jac_rows.len(),
            hess_rows,
            hess_cols,
            hess_entry_map,
            inner_hess_nnz: inner_hess_rows.len(),
        })
    }

    fn block(&self) -> EqualityBlock {
        EqualityBlock {
            rows: self.rows.clone(),
            vars: self.vars.clone(),
        }
    }

    fn expand_x(&self, x_block: &[f64]) -> Vec<f64> {
        let mut x_full = self.fixed_x.clone();
        for (local, &var) in self.vars.iter().enumerate() {
            x_full[var] = x_block[local];
        }
        x_full
    }
}

impl NlpProblem for AuxiliaryBlockProblem<'_> {
    fn num_variables(&self) -> usize {
        self.vars.len()
    }

    fn num_constraints(&self) -> usize {
        self.rows.len()
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        let mut x_l_full = vec![0.0; self.n_orig];
        let mut x_u_full = vec![0.0; self.n_orig];
        self.inner.bounds(&mut x_l_full, &mut x_u_full);
        for (local, &var) in self.vars.iter().enumerate() {
            x_l[local] = x_l_full[var];
            x_u[local] = x_u_full[var];
        }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        let mut g_l_full = vec![0.0; self.m_orig];
        let mut g_u_full = vec![0.0; self.m_orig];
        self.inner.constraint_bounds(&mut g_l_full, &mut g_u_full);
        for (local, &row) in self.rows.iter().enumerate() {
            g_l[local] = g_l_full[row];
            g_u[local] = g_u_full[row];
        }
    }

    fn initial_point(&self, x0: &mut [f64]) {
        for (local, &var) in self.vars.iter().enumerate() {
            x0[local] = self.fixed_x[var];
        }
    }

    fn objective(&self, _x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
        *obj = 0.0;
        true
    }

    fn gradient(&self, _x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
        grad.fill(0.0);
        true
    }

    fn constraints(&self, x: &[f64], new_x: bool, g: &mut [f64]) -> bool {
        let x_full = self.expand_x(x);
        let mut g_full = vec![0.0; self.m_orig];
        if !self.inner.constraints(&x_full, new_x, &mut g_full) {
            return false;
        }
        for (local, &row) in self.rows.iter().enumerate() {
            g[local] = g_full[row];
        }
        true
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.jac_rows.clone(), self.jac_cols.clone())
    }

    fn jacobian_values(&self, x: &[f64], new_x: bool, vals: &mut [f64]) -> bool {
        if self.jac_entry_map.is_empty() {
            return true;
        }
        let x_full = self.expand_x(x);
        let mut inner_vals = vec![0.0; self.inner_jac_nnz];
        if !self.inner.jacobian_values(&x_full, new_x, &mut inner_vals) {
            return false;
        }
        for (local, &inner_idx) in self.jac_entry_map.iter().enumerate() {
            vals[local] = inner_vals[inner_idx];
        }
        true
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hess_rows.clone(), self.hess_cols.clone())
    }

    fn hessian_values(
        &self,
        x: &[f64],
        new_x: bool,
        _obj_factor: f64,
        lambda: &[f64],
        vals: &mut [f64],
    ) -> bool {
        if self.hess_entry_map.is_empty() {
            return true;
        }
        let x_full = self.expand_x(x);
        let mut lambda_full = vec![0.0; self.m_orig];
        for (local, &row) in self.rows.iter().enumerate() {
            lambda_full[row] = lambda[local];
        }

        let mut inner_vals = vec![0.0; self.inner_hess_nnz];
        if !self
            .inner
            .hessian_values(&x_full, new_x, 0.0, &lambda_full, &mut inner_vals)
        {
            return false;
        }
        for (local, &inner_idx) in self.hess_entry_map.iter().enumerate() {
            vals[local] = inner_vals[inner_idx];
        }
        true
    }
}

#[allow(dead_code)]
pub(crate) struct AuxiliaryReducedProblem<'a> {
    inner: &'a dyn NlpProblem,
    n_orig: usize,
    m_orig: usize,
    fixed_x: Vec<f64>,
    var_map: Vec<usize>,
    constr_map: Vec<usize>,
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    jac_entry_map: Vec<usize>,
    inner_jac_nnz: usize,
    hess_rows: Vec<usize>,
    hess_cols: Vec<usize>,
    hess_entry_map: Vec<usize>,
    inner_hess_nnz: usize,
}

impl<'a> AuxiliaryReducedProblem<'a> {
    pub(crate) fn new(
        inner: &'a dyn NlpProblem,
        candidates: &[PresolveCandidate],
        fixed_x: Vec<f64>,
    ) -> Result<Self, AuxiliarySolveError> {
        let n_orig = inner.num_variables();
        let m_orig = inner.num_constraints();
        if fixed_x.len() != n_orig {
            return Err(AuxiliarySolveError::InvalidBlock {
                block: EqualityBlock {
                    rows: Vec::new(),
                    vars: Vec::new(),
                },
                reason: "fixed_x length does not match problem variables",
            });
        }

        let mut fixed_vars = vec![false; n_orig];
        let mut removed_constraints = vec![false; m_orig];
        for candidate in candidates {
            for block in &candidate.blocks {
                validate_auxiliary_block(block, &fixed_x, n_orig, m_orig)?;
                let rank = auxiliary_block_jacobian_rank(inner, block, &fixed_x)?;
                let expected = block.vars.len();
                if rank < expected {
                    return Err(AuxiliarySolveError::RankDeficientBlock {
                        block: block.clone(),
                        rank,
                        expected,
                    });
                }
                for &var in &block.vars {
                    fixed_vars[var] = true;
                }
                for &row in &block.rows {
                    removed_constraints[row] = true;
                }
            }
        }

        let var_map: Vec<_> = (0..n_orig).filter(|&var| !fixed_vars[var]).collect();
        let constr_map: Vec<_> = (0..m_orig)
            .filter(|&row| !removed_constraints[row])
            .collect();

        let mut orig_to_reduced_var = vec![None; n_orig];
        for (reduced, &orig) in var_map.iter().enumerate() {
            orig_to_reduced_var[orig] = Some(reduced);
        }
        let mut orig_to_reduced_constr = vec![None; m_orig];
        for (reduced, &orig) in constr_map.iter().enumerate() {
            orig_to_reduced_constr[orig] = Some(reduced);
        }

        let (inner_jac_rows, inner_jac_cols) = inner.jacobian_structure();
        let mut jac_rows = Vec::new();
        let mut jac_cols = Vec::new();
        let mut jac_entry_map = Vec::new();
        for (idx, (&row, &col)) in inner_jac_rows.iter().zip(inner_jac_cols.iter()).enumerate() {
            if row >= m_orig || col >= n_orig {
                continue;
            }
            if let (Some(reduced_row), Some(reduced_col)) =
                (orig_to_reduced_constr[row], orig_to_reduced_var[col])
            {
                jac_rows.push(reduced_row);
                jac_cols.push(reduced_col);
                jac_entry_map.push(idx);
            }
        }

        let (inner_hess_rows, inner_hess_cols) = inner.hessian_structure();
        let mut hess_rows = Vec::new();
        let mut hess_cols = Vec::new();
        let mut hess_entry_map = Vec::new();
        for (idx, (&row, &col)) in inner_hess_rows
            .iter()
            .zip(inner_hess_cols.iter())
            .enumerate()
        {
            if row >= n_orig || col >= n_orig {
                continue;
            }
            if let (Some(reduced_row), Some(reduced_col)) =
                (orig_to_reduced_var[row], orig_to_reduced_var[col])
            {
                if reduced_row >= reduced_col {
                    hess_rows.push(reduced_row);
                    hess_cols.push(reduced_col);
                } else {
                    hess_rows.push(reduced_col);
                    hess_cols.push(reduced_row);
                }
                hess_entry_map.push(idx);
            }
        }

        Ok(Self {
            inner,
            n_orig,
            m_orig,
            fixed_x,
            var_map,
            constr_map,
            jac_rows,
            jac_cols,
            jac_entry_map,
            inner_jac_nnz: inner_jac_rows.len(),
            hess_rows,
            hess_cols,
            hess_entry_map,
            inner_hess_nnz: inner_hess_rows.len(),
        })
    }

    pub(crate) fn did_reduce(&self) -> bool {
        self.var_map.len() < self.n_orig || self.constr_map.len() < self.m_orig
    }

    pub(crate) fn num_fixed(&self) -> usize {
        self.n_orig - self.var_map.len()
    }

    pub(crate) fn num_removed_constraints(&self) -> usize {
        self.m_orig - self.constr_map.len()
    }

    pub(crate) fn reduced_x_scaling(&self, scaling: &[f64]) -> Option<Vec<f64>> {
        if scaling.len() != self.n_orig {
            return None;
        }
        Some(self.var_map.iter().map(|&orig| scaling[orig]).collect())
    }

    pub(crate) fn reduced_g_scaling(&self, scaling: &[f64]) -> Option<Vec<f64>> {
        if scaling.len() != self.m_orig {
            return None;
        }
        Some(self.constr_map.iter().map(|&orig| scaling[orig]).collect())
    }

    fn expand_x(&self, x_reduced: &[f64]) -> Vec<f64> {
        let mut x_full = self.fixed_x.clone();
        for (reduced, &orig) in self.var_map.iter().enumerate() {
            x_full[orig] = x_reduced[reduced];
        }
        x_full
    }

    #[cfg(test)]
    pub(crate) fn unmap_solution(&self, reduced: &SolveResult) -> SolveResult {
        self.unmap_solution_impl(reduced, None)
    }

    pub(crate) fn unmap_solution_with_options(
        &self,
        reduced: &SolveResult,
        options: &SolverOptions,
    ) -> SolveResult {
        self.unmap_solution_impl(reduced, Some(options))
    }

    fn unmap_solution_impl(
        &self,
        reduced: &SolveResult,
        options: Option<&SolverOptions>,
    ) -> SolveResult {
        let x_full = self.expand_x(&reduced.x);

        let mut constraint_multipliers = vec![0.0; self.m_orig];
        for (reduced_idx, &orig_idx) in self.constr_map.iter().enumerate() {
            if reduced_idx < reduced.constraint_multipliers.len() {
                constraint_multipliers[orig_idx] = reduced.constraint_multipliers[reduced_idx];
            }
        }

        let mut bound_multipliers_lower = vec![0.0; self.n_orig];
        let mut bound_multipliers_upper = vec![0.0; self.n_orig];
        for (reduced_idx, &orig_idx) in self.var_map.iter().enumerate() {
            if reduced_idx < reduced.bound_multipliers_lower.len() {
                bound_multipliers_lower[orig_idx] = reduced.bound_multipliers_lower[reduced_idx];
            }
            if reduced_idx < reduced.bound_multipliers_upper.len() {
                bound_multipliers_upper[orig_idx] = reduced.bound_multipliers_upper[reduced_idx];
            }
        }

        if let Err(reason) = self.reconstruct_removed_constraint_multipliers(
            &x_full,
            &mut constraint_multipliers,
            &bound_multipliers_lower,
            &bound_multipliers_upper,
        ) {
            if let Some(options) = options {
                if options.print_level >= 5 {
                    rip_log!(
                        "ripopt: Auxiliary multiplier reconstruction skipped: {}",
                        reason
                    );
                }
            }
        }

        let mut objective = reduced.objective;
        let _ = self.inner.objective(&x_full, true, &mut objective);

        let mut constraint_values = vec![0.0; self.m_orig];
        if self.m_orig > 0 {
            let _ = self
                .inner
                .constraints(&x_full, true, &mut constraint_values);
        }

        SolveResult {
            x: x_full,
            objective,
            constraint_multipliers,
            bound_multipliers_lower,
            bound_multipliers_upper,
            constraint_values,
            status: reduced.status,
            iterations: reduced.iterations,
            diagnostics: reduced.diagnostics.clone(),
        }
    }

    fn reconstruct_removed_constraint_multipliers(
        &self,
        x_full: &[f64],
        constraint_multipliers: &mut [f64],
        bound_multipliers_lower: &[f64],
        bound_multipliers_upper: &[f64],
    ) -> Result<(), &'static str> {
        let removed_rows = removed_indices(self.m_orig, &self.constr_map);
        let removed_vars = removed_indices(self.n_orig, &self.var_map);
        if removed_rows.is_empty() && removed_vars.is_empty() {
            return Ok(());
        }
        if removed_rows.len() != removed_vars.len() {
            return Err("removed auxiliary row/variable counts are not square");
        }
        if removed_rows.is_empty() {
            return Ok(());
        }
        if self.removed_vars_have_bound_ambiguity(&removed_vars, x_full) {
            return Err("removed auxiliary variable is active at a finite bound");
        }

        let mut grad = vec![0.0; self.n_orig];
        if !self.inner.gradient(x_full, true, &mut grad) {
            return Err("full gradient evaluation failed");
        }

        let (jac_rows, jac_cols) = self.inner.jacobian_structure();
        if jac_rows.len() != jac_cols.len() {
            return Err("Jacobian structure has mismatched row/column lengths");
        }
        let mut jac_vals = vec![0.0; jac_rows.len()];
        if !self.inner.jacobian_values(x_full, true, &mut jac_vals) {
            return Err("full Jacobian evaluation failed");
        }

        let row_pos = index_positions(self.m_orig, &removed_rows);
        let var_pos = index_positions(self.n_orig, &removed_vars);
        let dim = removed_rows.len();
        let mut system = vec![0.0; dim * dim];
        let mut rhs: Vec<f64> = removed_vars
            .iter()
            .map(|&var| -(grad[var] - bound_multipliers_lower[var] + bound_multipliers_upper[var]))
            .collect();

        for (idx, (&row, &col)) in jac_rows.iter().zip(jac_cols.iter()).enumerate() {
            if row >= self.m_orig || col >= self.n_orig {
                continue;
            }
            let val = jac_vals[idx];
            if !val.is_finite() {
                return Err("Jacobian contains a non-finite value");
            }
            let Some(var_local) = var_pos[col] else {
                continue;
            };
            if let Some(row_local) = row_pos[row] {
                system[var_local * dim + row_local] += val;
            } else {
                rhs[var_local] -= val * constraint_multipliers[row];
            }
        }

        let lambda = solve_dense_square_system(system, rhs, dim)?;
        for (local, &row) in removed_rows.iter().enumerate() {
            constraint_multipliers[row] = lambda[local];
        }
        Ok(())
    }

    fn removed_vars_have_bound_ambiguity(&self, removed_vars: &[usize], x_full: &[f64]) -> bool {
        let mut x_l = vec![0.0; self.n_orig];
        let mut x_u = vec![0.0; self.n_orig];
        self.inner.bounds(&mut x_l, &mut x_u);
        removed_vars.iter().any(|&var| {
            let x = x_full[var];
            let scale_l = x.abs().max(x_l[var].abs()).max(1.0);
            let scale_u = x.abs().max(x_u[var].abs()).max(1.0);
            (x_l[var].is_finite() && x <= x_l[var] + 1e-8 * scale_l)
                || (x_u[var].is_finite() && x >= x_u[var] - 1e-8 * scale_u)
        })
    }
}

fn removed_indices(total: usize, kept: &[usize]) -> Vec<usize> {
    let mut is_kept = vec![false; total];
    for &idx in kept {
        if idx < total {
            is_kept[idx] = true;
        }
    }
    (0..total).filter(|&idx| !is_kept[idx]).collect()
}

fn index_positions(total: usize, indices: &[usize]) -> Vec<Option<usize>> {
    let mut positions = vec![None; total];
    for (pos, &idx) in indices.iter().enumerate() {
        if idx < total {
            positions[idx] = Some(pos);
        }
    }
    positions
}

fn solve_dense_square_system(
    mut matrix: Vec<f64>,
    mut rhs: Vec<f64>,
    dim: usize,
) -> Result<Vec<f64>, &'static str> {
    if dim == 0 {
        return Ok(Vec::new());
    }
    if matrix.len() != dim * dim || rhs.len() != dim {
        return Err("dense system dimensions are inconsistent");
    }
    if matrix.iter().any(|value| !value.is_finite()) {
        return Err("dense system matrix contains a non-finite value");
    }
    if rhs.iter().any(|value| !value.is_finite()) {
        return Err("dense system RHS contains a non-finite value");
    }

    let matrix_norm = matrix
        .iter()
        .fold(0.0_f64, |acc, &value| acc.max(value.abs()));
    if matrix_norm == 0.0 {
        return Err("dense system matrix is zero");
    }

    let original_matrix = matrix.clone();
    let original_rhs = rhs.clone();
    let pivot_tol = (dim as f64) * matrix_norm * 1e-10;

    for col in 0..dim {
        let pivot_row = (col..dim)
            .max_by(|&a, &b| {
                matrix[a * dim + col]
                    .abs()
                    .partial_cmp(&matrix[b * dim + col].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .ok_or("dense system has no pivot row")?;
        let pivot_abs = matrix[pivot_row * dim + col].abs();
        if pivot_abs <= pivot_tol {
            return Err("dense system is singular or ill-conditioned");
        }
        if pivot_row != col {
            for j in col..dim {
                matrix.swap(col * dim + j, pivot_row * dim + j);
            }
            rhs.swap(col, pivot_row);
        }

        let pivot = matrix[col * dim + col];
        for row in (col + 1)..dim {
            let factor = matrix[row * dim + col] / pivot;
            matrix[row * dim + col] = 0.0;
            for j in (col + 1)..dim {
                matrix[row * dim + j] -= factor * matrix[col * dim + j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    let mut solution = vec![0.0; dim];
    for i in (0..dim).rev() {
        let mut sum = rhs[i];
        for j in (i + 1)..dim {
            sum -= matrix[i * dim + j] * solution[j];
        }
        let pivot = matrix[i * dim + i];
        if pivot.abs() <= pivot_tol {
            return Err("dense system is singular or ill-conditioned");
        }
        solution[i] = sum / pivot;
    }
    if solution.iter().any(|value| !value.is_finite()) {
        return Err("dense system solution is non-finite");
    }

    let residual = dense_residual_inf(&original_matrix, &solution, &original_rhs, dim);
    let rhs_norm = original_rhs
        .iter()
        .fold(0.0_f64, |acc, &value| acc.max(value.abs()));
    if residual > 1e-8 * rhs_norm.max(1.0) {
        return Err("dense system residual is too large");
    }

    Ok(solution)
}

fn dense_residual_inf(matrix: &[f64], x: &[f64], rhs: &[f64], dim: usize) -> f64 {
    let mut residual: f64 = 0.0;
    for row in 0..dim {
        let mut ax = 0.0;
        for col in 0..dim {
            ax += matrix[row * dim + col] * x[col];
        }
        residual = residual.max((ax - rhs[row]).abs());
    }
    residual
}

impl NlpProblem for AuxiliaryReducedProblem<'_> {
    fn num_variables(&self) -> usize {
        self.var_map.len()
    }

    fn num_constraints(&self) -> usize {
        self.constr_map.len()
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        let mut x_l_full = vec![0.0; self.n_orig];
        let mut x_u_full = vec![0.0; self.n_orig];
        self.inner.bounds(&mut x_l_full, &mut x_u_full);
        for (reduced, &orig) in self.var_map.iter().enumerate() {
            x_l[reduced] = x_l_full[orig];
            x_u[reduced] = x_u_full[orig];
        }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        let mut g_l_full = vec![0.0; self.m_orig];
        let mut g_u_full = vec![0.0; self.m_orig];
        self.inner.constraint_bounds(&mut g_l_full, &mut g_u_full);
        for (reduced, &orig) in self.constr_map.iter().enumerate() {
            g_l[reduced] = g_l_full[orig];
            g_u[reduced] = g_u_full[orig];
        }
    }

    fn initial_point(&self, x0: &mut [f64]) {
        for (reduced, &orig) in self.var_map.iter().enumerate() {
            x0[reduced] = self.fixed_x[orig];
        }
    }

    fn initial_multipliers(&self, lam_g: &mut [f64], z_l: &mut [f64], z_u: &mut [f64]) -> bool {
        let mut lam_g_full = vec![0.0; self.m_orig];
        let mut z_l_full = vec![0.0; self.n_orig];
        let mut z_u_full = vec![0.0; self.n_orig];
        if !self
            .inner
            .initial_multipliers(&mut lam_g_full, &mut z_l_full, &mut z_u_full)
        {
            return false;
        }
        for (reduced, &orig) in self.constr_map.iter().enumerate() {
            lam_g[reduced] = lam_g_full[orig];
        }
        for (reduced, &orig) in self.var_map.iter().enumerate() {
            z_l[reduced] = z_l_full[orig];
            z_u[reduced] = z_u_full[orig];
        }
        true
    }

    fn objective(&self, x: &[f64], new_x: bool, obj: &mut f64) -> bool {
        let x_full = self.expand_x(x);
        self.inner.objective(&x_full, new_x, obj)
    }

    fn gradient(&self, x: &[f64], new_x: bool, grad: &mut [f64]) -> bool {
        let x_full = self.expand_x(x);
        let mut grad_full = vec![0.0; self.n_orig];
        if !self.inner.gradient(&x_full, new_x, &mut grad_full) {
            return false;
        }
        for (reduced, &orig) in self.var_map.iter().enumerate() {
            grad[reduced] = grad_full[orig];
        }
        true
    }

    fn constraints(&self, x: &[f64], new_x: bool, g: &mut [f64]) -> bool {
        let x_full = self.expand_x(x);
        let mut g_full = vec![0.0; self.m_orig];
        if !self.inner.constraints(&x_full, new_x, &mut g_full) {
            return false;
        }
        for (reduced, &orig) in self.constr_map.iter().enumerate() {
            g[reduced] = g_full[orig];
        }
        true
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.jac_rows.clone(), self.jac_cols.clone())
    }

    fn jacobian_values(&self, x: &[f64], new_x: bool, vals: &mut [f64]) -> bool {
        if self.jac_entry_map.is_empty() {
            return true;
        }
        let x_full = self.expand_x(x);
        let mut inner_vals = vec![0.0; self.inner_jac_nnz];
        if !self.inner.jacobian_values(&x_full, new_x, &mut inner_vals) {
            return false;
        }
        for (reduced, &inner_idx) in self.jac_entry_map.iter().enumerate() {
            vals[reduced] = inner_vals[inner_idx];
        }
        true
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hess_rows.clone(), self.hess_cols.clone())
    }

    fn hessian_values(
        &self,
        x: &[f64],
        new_x: bool,
        obj_factor: f64,
        lambda: &[f64],
        vals: &mut [f64],
    ) -> bool {
        if self.hess_entry_map.is_empty() {
            return true;
        }

        let x_full = self.expand_x(x);
        let mut lambda_full = vec![0.0; self.m_orig];
        for (reduced, &orig) in self.constr_map.iter().enumerate() {
            lambda_full[orig] = lambda[reduced];
        }

        let mut inner_vals = vec![0.0; self.inner_hess_nnz];
        if !self
            .inner
            .hessian_values(&x_full, new_x, obj_factor, &lambda_full, &mut inner_vals)
        {
            return false;
        }
        for (reduced, &inner_idx) in self.hess_entry_map.iter().enumerate() {
            vals[reduced] = inner_vals[inner_idx];
        }
        true
    }
}

fn auxiliary_block_jacobian_rank(
    inner: &dyn NlpProblem,
    block: &EqualityBlock,
    fixed_x: &[f64],
) -> Result<usize, AuxiliarySolveError> {
    let block_problem = AuxiliaryBlockProblem::new(inner, block, fixed_x)?;
    let rows = block_problem.rows.len();
    let cols = block_problem.vars.len();
    if rows == 0 || cols == 0 {
        return Ok(0);
    }

    let x_block: Vec<_> = block_problem.vars.iter().map(|&var| fixed_x[var]).collect();
    let mut jac_vals = vec![0.0; block_problem.jac_rows.len()];
    if !block_problem.jacobian_values(&x_block, true, &mut jac_vals) {
        return Err(AuxiliarySolveError::EvaluationFailed {
            block: block.clone(),
        });
    }

    let mut dense = vec![0.0; rows * cols];
    for (idx, (&row, &col)) in block_problem
        .jac_rows
        .iter()
        .zip(block_problem.jac_cols.iter())
        .enumerate()
    {
        dense[row * cols + col] += jac_vals[idx];
    }
    Ok(dense_numeric_rank(&mut dense, rows, cols))
}

fn dense_numeric_rank(matrix: &mut [f64], rows: usize, cols: usize) -> usize {
    let max_abs = matrix
        .iter()
        .fold(0.0_f64, |acc, &value| acc.max(value.abs()));
    if max_abs == 0.0 || !max_abs.is_finite() {
        return 0;
    }

    let tol = (rows.max(cols) as f64) * max_abs * 1e-10;
    let mut rank = 0usize;
    for col in 0..cols {
        let pivot_row = (rank..rows).max_by(|&a, &b| {
            matrix[a * cols + col]
                .abs()
                .partial_cmp(&matrix[b * cols + col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(pivot_row) = pivot_row else {
            break;
        };
        if matrix[pivot_row * cols + col].abs() <= tol {
            continue;
        }

        if pivot_row != rank {
            for j in 0..cols {
                matrix.swap(rank * cols + j, pivot_row * cols + j);
            }
        }

        let pivot = matrix[rank * cols + col];
        for row in (rank + 1)..rows {
            let factor = matrix[row * cols + col] / pivot;
            matrix[row * cols + col] = 0.0;
            for j in (col + 1)..cols {
                matrix[row * cols + j] -= factor * matrix[rank * cols + j];
            }
        }
        rank += 1;
        if rank == rows {
            break;
        }
    }
    rank
}

pub(crate) fn solve_auxiliary_blocks(
    problem: &dyn NlpProblem,
    candidates: &[PresolveCandidate],
    options: &SolverOptions,
    solve_start: Instant,
) -> Result<AuxiliarySolveOutcome, AuxiliarySolveError> {
    let mut x_full = vec![0.0; problem.num_variables()];
    problem.initial_point(&mut x_full);
    solve_auxiliary_blocks_from(problem, candidates, options, solve_start, &mut x_full)
}

pub(crate) fn solve_auxiliary_blocks_from(
    problem: &dyn NlpProblem,
    candidates: &[PresolveCandidate],
    options: &SolverOptions,
    solve_start: Instant,
    x_full: &mut [f64],
) -> Result<AuxiliarySolveOutcome, AuxiliarySolveError> {
    let mut blocks_solved = 0;
    let mut max_residual: f64 = 0.0;

    for candidate in candidates {
        for block in &candidate.blocks {
            let block_problem = AuxiliaryBlockProblem::new(problem, block, x_full)?;
            let Some(aux_options) = auxiliary_solver_options(options, solve_start) else {
                return Err(AuxiliarySolveError::TimeBudgetExceeded { blocks_solved });
            };
            let result = crate::solve(&block_problem, &aux_options);
            let residual = auxiliary_result_residual(&block_problem, &result)?;

            if !(matches!(
                result.status,
                SolveStatus::Optimal | SolveStatus::Acceptable
            ) && residual <= options.auxiliary_tol)
            {
                return Err(AuxiliarySolveError::BlockSolveFailed {
                    block: block_problem.block(),
                    status: result.status,
                    residual,
                });
            }

            for (local, &var) in block.vars.iter().enumerate() {
                x_full[var] = result.x[local];
            }
            blocks_solved += 1;
            max_residual = max_residual.max(residual);
        }
    }

    Ok(AuxiliarySolveOutcome {
        x: x_full.to_vec(),
        blocks_solved,
        max_residual,
    })
}

fn validate_auxiliary_block(
    block: &EqualityBlock,
    fixed_x: &[f64],
    n_orig: usize,
    m_orig: usize,
) -> Result<(), AuxiliarySolveError> {
    if fixed_x.len() != n_orig {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "fixed_x length does not match problem variables",
        });
    }
    if block.rows.is_empty() || block.vars.is_empty() {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "block must contain at least one row and one variable",
        });
    }
    if block.rows.iter().any(|&row| row >= m_orig) {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "block row is out of range",
        });
    }
    if block.vars.iter().any(|&var| var >= n_orig) {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "block variable is out of range",
        });
    }
    if has_duplicates(&block.rows) {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "block rows must be unique",
        });
    }
    if has_duplicates(&block.vars) {
        return Err(AuxiliarySolveError::InvalidBlock {
            block: block.clone(),
            reason: "block variables must be unique",
        });
    }
    Ok(())
}

fn has_duplicates(values: &[usize]) -> bool {
    let mut values = values.to_vec();
    values.sort_unstable();
    values.windows(2).any(|pair| pair[0] == pair[1])
}

fn auxiliary_solver_options(
    options: &SolverOptions,
    solve_start: Instant,
) -> Option<SolverOptions> {
    let mut aux_options = options.clone();
    aux_options.enable_preprocessing = false;
    aux_options.warm_start = false;
    aux_options.warm_start_y = None;
    aux_options.warm_start_z_l = None;
    aux_options.warm_start_z_u = None;
    aux_options.user_obj_scaling = None;
    aux_options.user_g_scaling = None;
    aux_options.user_x_scaling = None;
    if options.max_wall_time > 0.0 {
        let remaining = options.max_wall_time - solve_start.elapsed().as_secs_f64();
        if remaining <= 0.0 {
            return None;
        }
        aux_options.max_wall_time = remaining;
    }
    Some(aux_options)
}

fn auxiliary_result_residual(
    problem: &AuxiliaryBlockProblem<'_>,
    result: &SolveResult,
) -> Result<f64, AuxiliarySolveError> {
    let m = problem.num_constraints();
    let g = if result.constraint_values.len() == m {
        result.constraint_values.clone()
    } else {
        let mut values = vec![0.0; m];
        if !problem.constraints(&result.x, true, &mut values) {
            return Err(AuxiliarySolveError::EvaluationFailed {
                block: problem.block(),
            });
        }
        values
    };

    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);

    let mut residual: f64 = 0.0;
    for i in 0..m {
        let violation = if !g[i].is_finite() {
            f64::INFINITY
        } else if g_l[i].is_finite() && g_u[i].is_finite() && (g_u[i] - g_l[i]).abs() <= 1e-12 {
            (g[i] - 0.5 * (g_l[i] + g_u[i])).abs()
        } else {
            let lower = if g_l[i].is_finite() {
                (g_l[i] - g[i]).max(0.0)
            } else {
                0.0
            };
            let upper = if g_u[i].is_finite() {
                (g[i] - g_u[i]).max(0.0)
            } else {
                0.0
            };
            lower.max(upper)
        };
        residual = residual.max(violation);
    }
    Ok(residual)
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BipartiteMatching {
    /// Local equality-row index -> matched original variable index.
    pub(crate) row_to_var: Vec<Option<usize>>,
    /// Original variable index -> matched local equality-row index.
    pub(crate) var_to_row: Vec<Option<usize>>,
    /// Unmatched local equality-row indices.
    pub(crate) unmatched_rows: Vec<usize>,
    /// Unmatched original variable indices.
    pub(crate) unmatched_vars: Vec<usize>,
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DulmageMendelsohnPartition {
    pub(crate) matching: BipartiteMatching,
    /// Overconstrained local equality-row indices.
    pub(crate) overconstrained_rows: Vec<usize>,
    /// Overconstrained original variable indices.
    pub(crate) overconstrained_vars: Vec<usize>,
    /// Square local equality-row indices.
    pub(crate) square_rows: Vec<usize>,
    /// Square original variable indices.
    pub(crate) square_vars: Vec<usize>,
    /// Underconstrained local equality-row indices.
    pub(crate) underconstrained_rows: Vec<usize>,
    /// Underconstrained original variable indices.
    pub(crate) underconstrained_vars: Vec<usize>,
    /// Unmatched local equality-row indices.
    pub(crate) unmatched_rows: Vec<usize>,
    /// Unmatched original variable indices.
    pub(crate) unmatched_vars: Vec<usize>,
}

#[allow(dead_code)]
pub(crate) fn find_presolve_candidates(
    problem: &dyn NlpProblem,
    tol: f64,
) -> Vec<PresolveCandidate> {
    let incidence = EqualityIncidence::from_problem(problem, tol);
    if incidence.row_global.is_empty() {
        return Vec::new();
    }

    let selected_rows: Vec<_> = (0..incidence.row_adj_vars.len()).collect();
    let selected_vars: Vec<_> = incidence
        .var_adj_rows
        .iter()
        .enumerate()
        .filter_map(|(var, rows)| (!rows.is_empty()).then_some(var))
        .collect();

    incidence
        .connected_components(&selected_rows, &selected_vars)
        .into_iter()
        .filter_map(|component| presolve_candidate_from_component(&incidence, component))
        .collect()
}

fn presolve_candidate_from_component(
    incidence: &EqualityIncidence,
    component: EqualityBlock,
) -> Option<PresolveCandidate> {
    if component.rows.len() != component.vars.len() || component.rows.is_empty() {
        return None;
    }

    let local_rows: Vec<_> = component
        .rows
        .iter()
        .map(|&row| incidence.row_local_for_global[row])
        .collect::<Option<_>>()?;

    if !is_closed_equality_component(incidence, &local_rows, &component.vars) {
        return None;
    }

    incidence
        .block_triangular_decomposition(&local_rows, &component.vars)
        .ok()
        .map(|blocks| PresolveCandidate { blocks })
}

fn is_closed_equality_component(
    incidence: &EqualityIncidence,
    local_rows: &[usize],
    vars: &[usize],
) -> bool {
    let mut selected_rows = vec![false; incidence.row_adj_vars.len()];
    let mut selected_vars = vec![false; incidence.n_vars];

    for &row in local_rows {
        if row >= selected_rows.len() {
            return false;
        }
        selected_rows[row] = true;
    }
    for &var in vars {
        if var >= selected_vars.len() {
            return false;
        }
        selected_vars[var] = true;
    }

    local_rows.iter().all(|&row| {
        incidence.row_adj_vars[row]
            .iter()
            .all(|&var| selected_vars[var])
    }) && vars.iter().all(|&var| {
        incidence.var_adj_rows[var]
            .iter()
            .all(|&row| selected_rows[row])
    })
}

#[derive(Debug, Clone, Copy)]
enum BipartiteNode {
    Row(usize),
    Var(usize),
}

impl EqualityIncidence {
    #[allow(dead_code)]
    pub(crate) fn from_problem(problem: &dyn NlpProblem, tol: f64) -> Self {
        let n_vars = problem.num_variables();
        let m_orig = problem.num_constraints();

        let mut g_l = vec![0.0; m_orig];
        let mut g_u = vec![0.0; m_orig];
        if m_orig > 0 {
            problem.constraint_bounds(&mut g_l, &mut g_u);
        }

        let mut row_global = Vec::new();
        let mut row_local_for_global = vec![None; m_orig];
        for i in 0..m_orig {
            if g_l[i].is_finite() && g_u[i].is_finite() && (g_l[i] - g_u[i]).abs() <= tol {
                row_local_for_global[i] = Some(row_global.len());
                row_global.push(i);
            }
        }

        let mut row_adj_vars = vec![Vec::new(); row_global.len()];
        let mut var_adj_rows = vec![Vec::new(); n_vars];
        let (jac_rows, jac_cols) = problem.jacobian_structure();
        for (&row, &col) in jac_rows.iter().zip(jac_cols.iter()) {
            if row >= m_orig || col >= n_vars {
                continue;
            }
            if let Some(local_row) = row_local_for_global[row] {
                row_adj_vars[local_row].push(col);
                var_adj_rows[col].push(local_row);
            }
        }

        for adj in &mut row_adj_vars {
            adj.sort_unstable();
            adj.dedup();
        }
        for adj in &mut var_adj_rows {
            adj.sort_unstable();
            adj.dedup();
        }

        Self {
            n_vars,
            m_orig,
            row_global,
            row_local_for_global,
            row_adj_vars,
            var_adj_rows,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn maximum_matching(&self) -> BipartiteMatching {
        let (row_adj_vars, _) = self.deterministic_adjacency();
        hopcroft_karp(self.n_vars, &row_adj_vars)
    }

    #[allow(dead_code)]
    pub(crate) fn dulmage_mendelsohn_partition(&self) -> DulmageMendelsohnPartition {
        let (row_adj_vars, var_adj_rows) = self.deterministic_adjacency();
        let matching = hopcroft_karp(self.n_vars, &row_adj_vars);

        let n_rows = row_adj_vars.len();
        let mut over_rows = vec![false; n_rows];
        let mut over_vars = vec![false; self.n_vars];
        let mut queue = VecDeque::new();

        for &row in &matching.unmatched_rows {
            over_rows[row] = true;
            queue.push_back(BipartiteNode::Row(row));
        }

        while let Some(node) = queue.pop_front() {
            match node {
                BipartiteNode::Row(row) => {
                    for &var in &row_adj_vars[row] {
                        if matching.row_to_var[row] == Some(var) || over_vars[var] {
                            continue;
                        }
                        over_vars[var] = true;
                        queue.push_back(BipartiteNode::Var(var));
                    }
                }
                BipartiteNode::Var(var) => {
                    if let Some(row) = matching.var_to_row[var] {
                        if !over_rows[row] {
                            over_rows[row] = true;
                            queue.push_back(BipartiteNode::Row(row));
                        }
                    }
                }
            }
        }

        let mut under_rows = vec![false; n_rows];
        let mut under_vars = vec![false; self.n_vars];
        let mut queue = VecDeque::new();

        for &var in &matching.unmatched_vars {
            under_vars[var] = true;
            queue.push_back(BipartiteNode::Var(var));
        }

        while let Some(node) = queue.pop_front() {
            match node {
                BipartiteNode::Var(var) => {
                    for &row in &var_adj_rows[var] {
                        if matching.var_to_row[var] == Some(row) || under_rows[row] {
                            continue;
                        }
                        under_rows[row] = true;
                        queue.push_back(BipartiteNode::Row(row));
                    }
                }
                BipartiteNode::Row(row) => {
                    if let Some(var) = matching.row_to_var[row] {
                        if !under_vars[var] {
                            under_vars[var] = true;
                            queue.push_back(BipartiteNode::Var(var));
                        }
                    }
                }
            }
        }

        let overconstrained_rows = indices_where(&over_rows);
        let overconstrained_vars = indices_where(&over_vars);
        let underconstrained_rows = indices_where(&under_rows);
        let underconstrained_vars = indices_where(&under_vars);
        let square_rows = (0..n_rows)
            .filter(|&row| {
                !over_rows[row] && !under_rows[row] && matching.row_to_var[row].is_some()
            })
            .collect();
        let square_vars = (0..self.n_vars)
            .filter(|&var| {
                !over_vars[var] && !under_vars[var] && matching.var_to_row[var].is_some()
            })
            .collect();

        DulmageMendelsohnPartition {
            unmatched_rows: matching.unmatched_rows.clone(),
            unmatched_vars: matching.unmatched_vars.clone(),
            matching,
            overconstrained_rows,
            overconstrained_vars,
            square_rows,
            square_vars,
            underconstrained_rows,
            underconstrained_vars,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn connected_components(
        &self,
        selected_local_rows: &[usize],
        selected_vars: &[usize],
    ) -> Vec<EqualityBlock> {
        let selected_local_rows =
            sorted_unique_bounded(selected_local_rows, self.row_adj_vars.len());
        let selected_vars = sorted_unique_bounded(selected_vars, self.n_vars);
        let mut row_selected = vec![false; self.row_adj_vars.len()];
        let mut var_selected = vec![false; self.n_vars];
        for &row in &selected_local_rows {
            row_selected[row] = true;
        }
        for &var in &selected_vars {
            var_selected[var] = true;
        }

        let (row_adj_vars, var_adj_rows) = self.deterministic_adjacency();
        let mut row_seen = vec![false; self.row_adj_vars.len()];
        let mut var_seen = vec![false; self.n_vars];
        let mut components = Vec::new();

        for &start_row in &selected_local_rows {
            if row_seen[start_row] {
                continue;
            }
            components.push(self.connected_component_from(
                BipartiteNode::Row(start_row),
                &row_selected,
                &var_selected,
                &row_adj_vars,
                &var_adj_rows,
                &mut row_seen,
                &mut var_seen,
            ));
        }

        for &start_var in &selected_vars {
            if var_seen[start_var] {
                continue;
            }
            components.push(self.connected_component_from(
                BipartiteNode::Var(start_var),
                &row_selected,
                &var_selected,
                &row_adj_vars,
                &var_adj_rows,
                &mut row_seen,
                &mut var_seen,
            ));
        }

        components.sort_by_key(equality_block_sort_key);
        components
    }

    #[allow(dead_code)]
    pub(crate) fn block_triangular_decomposition(
        &self,
        selected_local_rows: &[usize],
        selected_vars: &[usize],
    ) -> Result<Vec<EqualityBlock>, BlockTriangularizationError> {
        let selected_local_rows =
            sorted_unique_bounded(selected_local_rows, self.row_adj_vars.len());
        let selected_vars = sorted_unique_bounded(selected_vars, self.n_vars);
        if selected_local_rows.len() != selected_vars.len() {
            return Err(BlockTriangularizationError::NonSquare {
                rows: selected_local_rows.len(),
                vars: selected_vars.len(),
            });
        }

        let (row_adj_vars, _) = self.deterministic_adjacency();
        let mut compact_var_for_global = vec![None; self.n_vars];
        for (compact_var, &global_var) in selected_vars.iter().enumerate() {
            compact_var_for_global[global_var] = Some(compact_var);
        }

        let restricted_row_adj_vars: Vec<Vec<usize>> = selected_local_rows
            .iter()
            .map(|&row| {
                row_adj_vars[row]
                    .iter()
                    .filter_map(|&var| compact_var_for_global[var])
                    .collect()
            })
            .collect();

        let matching = hopcroft_karp(selected_vars.len(), &restricted_row_adj_vars);
        if !matching.unmatched_rows.is_empty() || !matching.unmatched_vars.is_empty() {
            return Err(BlockTriangularizationError::ImperfectMatching {
                unmatched_rows: matching
                    .unmatched_rows
                    .iter()
                    .map(|&row| self.row_global[selected_local_rows[row]])
                    .collect(),
                unmatched_vars: matching
                    .unmatched_vars
                    .iter()
                    .map(|&var| selected_vars[var])
                    .collect(),
            });
        }

        let mut matched_row_for_compact_var = vec![None; selected_vars.len()];
        for (row, &matched_var) in matching.row_to_var.iter().enumerate() {
            if let Some(var) = matched_var {
                matched_row_for_compact_var[var] = Some(row);
            }
        }

        let mut dependencies = vec![Vec::new(); selected_local_rows.len()];
        for (row, adj_vars) in restricted_row_adj_vars.iter().enumerate() {
            for &var in adj_vars {
                if matching.row_to_var[row] == Some(var) {
                    continue;
                }
                if let Some(predecessor_row) = matched_row_for_compact_var[var] {
                    dependencies[predecessor_row].push(row);
                }
            }
        }
        for adj in &mut dependencies {
            adj.sort_unstable();
            adj.dedup();
        }

        let components = tarjan_strongly_connected_components(&dependencies);
        Ok(topologically_order_blocks(
            &components,
            &dependencies,
            &selected_local_rows,
            &selected_vars,
            &matching.row_to_var,
            &self.row_global,
        ))
    }

    fn connected_component_from(
        &self,
        start: BipartiteNode,
        row_selected: &[bool],
        var_selected: &[bool],
        row_adj_vars: &[Vec<usize>],
        var_adj_rows: &[Vec<usize>],
        row_seen: &mut [bool],
        var_seen: &mut [bool],
    ) -> EqualityBlock {
        let mut queue = VecDeque::new();
        let mut rows = Vec::new();
        let mut vars = Vec::new();

        match start {
            BipartiteNode::Row(row) => {
                row_seen[row] = true;
                queue.push_back(BipartiteNode::Row(row));
            }
            BipartiteNode::Var(var) => {
                var_seen[var] = true;
                queue.push_back(BipartiteNode::Var(var));
            }
        }

        while let Some(node) = queue.pop_front() {
            match node {
                BipartiteNode::Row(row) => {
                    rows.push(self.row_global[row]);
                    for &var in &row_adj_vars[row] {
                        if var_selected[var] && !var_seen[var] {
                            var_seen[var] = true;
                            queue.push_back(BipartiteNode::Var(var));
                        }
                    }
                }
                BipartiteNode::Var(var) => {
                    vars.push(var);
                    for &row in &var_adj_rows[var] {
                        if row_selected[row] && !row_seen[row] {
                            row_seen[row] = true;
                            queue.push_back(BipartiteNode::Row(row));
                        }
                    }
                }
            }
        }

        rows.sort_unstable();
        vars.sort_unstable();
        EqualityBlock { rows, vars }
    }

    fn deterministic_adjacency(&self) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let mut row_adj_vars = Vec::with_capacity(self.row_adj_vars.len());
        let mut var_adj_rows = vec![Vec::new(); self.n_vars];

        for (row, adj) in self.row_adj_vars.iter().enumerate() {
            let mut vars: Vec<usize> = adj
                .iter()
                .copied()
                .filter(|&var| var < self.n_vars)
                .collect();
            vars.sort_unstable();
            vars.dedup();

            for &var in &vars {
                var_adj_rows[var].push(row);
            }
            row_adj_vars.push(vars);
        }

        for adj in &mut var_adj_rows {
            adj.sort_unstable();
            adj.dedup();
        }

        (row_adj_vars, var_adj_rows)
    }
}

fn sorted_unique_bounded(indices: &[usize], upper_bound: usize) -> Vec<usize> {
    let mut indices: Vec<_> = indices
        .iter()
        .copied()
        .filter(|&idx| idx < upper_bound)
        .collect();
    indices.sort_unstable();
    indices.dedup();
    indices
}

fn equality_block_sort_key(block: &EqualityBlock) -> (usize, usize) {
    (
        block.rows.first().copied().unwrap_or(usize::MAX),
        block.vars.first().copied().unwrap_or(usize::MAX),
    )
}

fn hopcroft_karp(n_vars: usize, row_adj_vars: &[Vec<usize>]) -> BipartiteMatching {
    let n_rows = row_adj_vars.len();
    let mut row_to_var = vec![None; n_rows];
    let mut var_to_row = vec![None; n_vars];
    let mut dist = vec![usize::MAX; n_rows];

    while matching_bfs(row_adj_vars, &row_to_var, &var_to_row, &mut dist) {
        for row in 0..n_rows {
            if row_to_var[row].is_none() {
                matching_dfs(
                    row,
                    row_adj_vars,
                    &mut row_to_var,
                    &mut var_to_row,
                    &mut dist,
                );
            }
        }
    }

    let unmatched_rows = row_to_var
        .iter()
        .enumerate()
        .filter_map(|(row, var)| var.is_none().then_some(row))
        .collect();
    let unmatched_vars = var_to_row
        .iter()
        .enumerate()
        .filter_map(|(var, row)| row.is_none().then_some(var))
        .collect();

    BipartiteMatching {
        row_to_var,
        var_to_row,
        unmatched_rows,
        unmatched_vars,
    }
}

fn matching_bfs(
    row_adj_vars: &[Vec<usize>],
    row_to_var: &[Option<usize>],
    var_to_row: &[Option<usize>],
    dist: &mut [usize],
) -> bool {
    let mut queue = VecDeque::new();
    let mut found_unmatched_var = false;

    for row in 0..row_adj_vars.len() {
        if row_to_var[row].is_none() {
            dist[row] = 0;
            queue.push_back(row);
        } else {
            dist[row] = usize::MAX;
        }
    }

    while let Some(row) = queue.pop_front() {
        for &var in &row_adj_vars[row] {
            if let Some(next_row) = var_to_row[var] {
                if dist[next_row] == usize::MAX {
                    dist[next_row] = dist[row] + 1;
                    queue.push_back(next_row);
                }
            } else {
                found_unmatched_var = true;
            }
        }
    }

    found_unmatched_var
}

fn matching_dfs(
    row: usize,
    row_adj_vars: &[Vec<usize>],
    row_to_var: &mut [Option<usize>],
    var_to_row: &mut [Option<usize>],
    dist: &mut [usize],
) -> bool {
    struct Frame {
        row: usize,
        next_edge: usize,
    }

    let mut stack = vec![Frame { row, next_edge: 0 }];
    let mut path_vars = Vec::new();

    while let Some(frame) = stack.last_mut() {
        if frame.next_edge == row_adj_vars[frame.row].len() {
            dist[frame.row] = usize::MAX;
            stack.pop();
            while path_vars.len() >= stack.len() && !path_vars.is_empty() {
                path_vars.pop();
            }
            continue;
        }

        let current_row = frame.row;
        let var = row_adj_vars[current_row][frame.next_edge];
        frame.next_edge += 1;

        if let Some(next_row) = var_to_row[var] {
            if dist[next_row] == dist[current_row] + 1 {
                path_vars.push(var);
                stack.push(Frame {
                    row: next_row,
                    next_edge: 0,
                });
            }
            continue;
        }

        path_vars.push(var);
        for (frame, &path_var) in stack.iter().zip(&path_vars) {
            row_to_var[frame.row] = Some(path_var);
            var_to_row[path_var] = Some(frame.row);
        }
        return true;
    }

    false
}

fn indices_where(flags: &[bool]) -> Vec<usize> {
    flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| flag.then_some(idx))
        .collect()
}

fn tarjan_strongly_connected_components(adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    struct Frame {
        node: usize,
        next_edge: usize,
    }

    let n = adj.len();
    let mut next_index = 0usize;
    let mut index = vec![None; n];
    let mut lowlink = vec![0usize; n];
    let mut scc_stack = Vec::new();
    let mut on_stack = vec![false; n];
    let mut components = Vec::new();

    for start in 0..n {
        if index[start].is_some() {
            continue;
        }

        index[start] = Some(next_index);
        lowlink[start] = next_index;
        next_index += 1;
        scc_stack.push(start);
        on_stack[start] = true;

        let mut dfs_stack = vec![Frame {
            node: start,
            next_edge: 0,
        }];

        while !dfs_stack.is_empty() {
            let frame_idx = dfs_stack.len() - 1;
            let node = dfs_stack[frame_idx].node;

            if dfs_stack[frame_idx].next_edge < adj[node].len() {
                let next = adj[node][dfs_stack[frame_idx].next_edge];
                dfs_stack[frame_idx].next_edge += 1;

                if index[next].is_none() {
                    index[next] = Some(next_index);
                    lowlink[next] = next_index;
                    next_index += 1;
                    scc_stack.push(next);
                    on_stack[next] = true;
                    dfs_stack.push(Frame {
                        node: next,
                        next_edge: 0,
                    });
                } else if on_stack[next] {
                    lowlink[node] = lowlink[node].min(index[next].expect("visited node has index"));
                }
                continue;
            }

            dfs_stack.pop();

            if lowlink[node] == index[node].expect("visited node has index") {
                let mut component = Vec::new();
                while let Some(member) = scc_stack.pop() {
                    on_stack[member] = false;
                    component.push(member);
                    if member == node {
                        break;
                    }
                }
                component.sort_unstable();
                components.push(component);
            }

            if let Some(parent) = dfs_stack.last() {
                lowlink[parent.node] = lowlink[parent.node].min(lowlink[node]);
            }
        }
    }

    components
}

fn topologically_order_blocks(
    components: &[Vec<usize>],
    dependencies: &[Vec<usize>],
    selected_local_rows: &[usize],
    selected_vars: &[usize],
    row_to_compact_var: &[Option<usize>],
    row_global: &[usize],
) -> Vec<EqualityBlock> {
    let mut component_for_row = vec![0; selected_local_rows.len()];
    for (component, rows) in components.iter().enumerate() {
        for &row in rows {
            component_for_row[row] = component;
        }
    }

    let mut component_edges = vec![Vec::new(); components.len()];
    let mut indegree = vec![0usize; components.len()];
    for (from_row, edges) in dependencies.iter().enumerate() {
        let from_component = component_for_row[from_row];
        for &to_row in edges {
            let to_component = component_for_row[to_row];
            if from_component != to_component {
                component_edges[from_component].push(to_component);
            }
        }
    }
    for edges in &mut component_edges {
        edges.sort_unstable();
        edges.dedup();
        for &to_component in edges.iter() {
            indegree[to_component] += 1;
        }
    }

    let blocks: Vec<_> = components
        .iter()
        .map(|rows| {
            let mut block_rows: Vec<_> = rows
                .iter()
                .map(|&row| row_global[selected_local_rows[row]])
                .collect();
            let mut block_vars: Vec<_> = rows
                .iter()
                .filter_map(|&row| row_to_compact_var[row].map(|var| selected_vars[var]))
                .collect();
            block_rows.sort_unstable();
            block_vars.sort_unstable();
            EqualityBlock {
                rows: block_rows,
                vars: block_vars,
            }
        })
        .collect();

    let mut ordered_blocks = Vec::with_capacity(blocks.len());
    let mut ready = BinaryHeap::new();
    for component in 0..components.len() {
        if indegree[component] == 0 {
            ready.push(Reverse((
                equality_block_sort_key(&blocks[component]),
                component,
            )));
        }
    }

    while let Some(Reverse((_key, component))) = ready.pop() {
        ordered_blocks.push(blocks[component].clone());
        for &next in &component_edges[component] {
            indegree[next] -= 1;
            if indegree[next] == 0 {
                ready.push(Reverse((equality_block_sort_key(&blocks[next]), next)));
            }
        }
    }

    ordered_blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[derive(Clone)]
    struct GraphProblem {
        n: usize,
        gl: Vec<f64>,
        gu: Vec<f64>,
        edges: Vec<(usize, usize)>,
        objective_vars: Vec<usize>,
    }

    impl NlpProblem for GraphProblem {
        fn num_variables(&self) -> usize {
            self.n
        }

        fn num_constraints(&self) -> usize {
            self.gl.len()
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            for i in 0..self.n {
                x_l[i] = f64::NEG_INFINITY;
                x_u[i] = f64::INFINITY;
            }
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l.copy_from_slice(&self.gl);
            g_u.copy_from_slice(&self.gu);
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0.fill(0.0);
        }

        fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = self.objective_vars.iter().map(|&var| x[var]).sum();
            true
        }

        fn gradient(&self, _x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad.fill(0.0);
            for &var in &self.objective_vars {
                grad[var] = 1.0;
            }
            true
        }

        fn constraints(&self, _x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g.fill(0.0);
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            self.edges.iter().copied().unzip()
        }

        fn jacobian_values(&self, _x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals.fill(1.0);
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![], vec![])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            _obj_factor: f64,
            _lambda: &[f64],
            _vals: &mut [f64],
        ) -> bool {
            true
        }
    }

    fn graph_problem(n: usize, bounds: &[(f64, f64)], edges: &[(usize, usize)]) -> GraphProblem {
        GraphProblem {
            n,
            gl: bounds.iter().map(|b| b.0).collect(),
            gu: bounds.iter().map(|b| b.1).collect(),
            edges: edges.to_vec(),
            objective_vars: Vec::new(),
        }
    }

    fn graph_problem_with_objective(
        n: usize,
        bounds: &[(f64, f64)],
        edges: &[(usize, usize)],
        objective_vars: &[usize],
    ) -> GraphProblem {
        GraphProblem {
            n,
            gl: bounds.iter().map(|b| b.0).collect(),
            gu: bounds.iter().map(|b| b.1).collect(),
            edges: edges.to_vec(),
            objective_vars: objective_vars.to_vec(),
        }
    }

    fn equality_bounds(n_rows: usize) -> Vec<(f64, f64)> {
        vec![(0.0, 0.0); n_rows]
    }

    fn quiet_aux_options() -> SolverOptions {
        let mut options = SolverOptions::default();
        options.print_level = 0;
        options.max_iter = 200;
        options.auxiliary_tol = 1e-7;
        options
    }

    struct TriangularAuxProblem;

    impl NlpProblem for TriangularAuxProblem {
        fn num_variables(&self) -> usize {
            2
        }

        fn num_constraints(&self) -> usize {
            2
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            x_l[0] = 0.1;
            x_u[0] = 10.0;
            x_l[1] = 0.0;
            x_u[1] = 10.0;
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l[0] = 4.0;
            g_u[0] = 4.0;
            g_l[1] = 0.0;
            g_u[1] = 0.0;
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0[0] = 1.0;
            x0[1] = 0.0;
        }

        fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = x[0] + x[1];
            true
        }

        fn gradient(&self, _x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad[0] = 1.0;
            grad[1] = 1.0;
            true
        }

        fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g[0] = x[0] * x[0];
            g[1] = x[1] - x[0] - 1.0;
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0, 1, 1], vec![0, 0, 1])
        }

        fn jacobian_values(&self, x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals[0] = 2.0 * x[0];
            vals[1] = -1.0;
            vals[2] = 1.0;
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0], vec![0])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            _obj_factor: f64,
            lambda: &[f64],
            vals: &mut [f64],
        ) -> bool {
            vals[0] = 2.0 * lambda[0];
            true
        }
    }

    struct InfeasibleBoundAuxProblem;

    impl NlpProblem for InfeasibleBoundAuxProblem {
        fn num_variables(&self) -> usize {
            1
        }

        fn num_constraints(&self) -> usize {
            1
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            x_l[0] = 0.0;
            x_u[0] = 1.0;
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l[0] = 2.0;
            g_u[0] = 2.0;
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0[0] = 0.5;
        }

        fn objective(&self, _x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = 0.0;
            true
        }

        fn gradient(&self, _x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad[0] = 0.0;
            true
        }

        fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g[0] = x[0];
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0], vec![0])
        }

        fn jacobian_values(&self, _x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals[0] = 1.0;
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![], vec![])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            _obj_factor: f64,
            _lambda: &[f64],
            _vals: &mut [f64],
        ) -> bool {
            true
        }
    }

    struct ReducedMappingProblem;

    impl NlpProblem for ReducedMappingProblem {
        fn num_variables(&self) -> usize {
            4
        }

        fn num_constraints(&self) -> usize {
            4
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            x_l[0] = -10.0;
            x_u[0] = 10.0;
            x_l[1] = 2.0;
            x_u[1] = 2.0;
            x_l[2] = -20.0;
            x_u[2] = 20.0;
            x_l[3] = 5.0;
            x_u[3] = 5.0;
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l[0] = 0.0;
            g_u[0] = 0.0;
            g_l[1] = f64::NEG_INFINITY;
            g_u[1] = 100.0;
            g_l[2] = 0.0;
            g_u[2] = 0.0;
            g_l[3] = 1.0;
            g_u[3] = 1.0;
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0.copy_from_slice(&[9.0, 2.0, 8.0, 5.0]);
        }

        fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = x[0] * x[0] + 4.0 * x[0] * x[2] + 3.0 * x[2] * x[2] + 5.0 * x[1] + 7.0 * x[3];
            true
        }

        fn gradient(&self, x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad[0] = 2.0 * x[0] + 4.0 * x[2];
            grad[1] = 5.0;
            grad[2] = 4.0 * x[0] + 6.0 * x[2];
            grad[3] = 7.0;
            true
        }

        fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g[0] = x[1] - 2.0;
            g[1] = x[0] + 10.0 * x[1] + 2.0 * x[2];
            g[2] = x[3] - 5.0;
            g[3] = x[0] * x[0] + x[2] - 3.0 * x[3];
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0, 1, 1, 1, 2, 3, 3, 3], vec![1, 0, 1, 2, 3, 0, 2, 3])
        }

        fn jacobian_values(&self, x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals[0] = 1.0;
            vals[1] = 1.0;
            vals[2] = 10.0;
            vals[3] = 2.0;
            vals[4] = 1.0;
            vals[5] = 2.0 * x[0];
            vals[6] = 1.0;
            vals[7] = -3.0;
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0, 1, 2, 2, 3, 3], vec![0, 0, 0, 2, 2, 3])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            obj_factor: f64,
            lambda: &[f64],
            vals: &mut [f64],
        ) -> bool {
            vals[0] = 2.0 * obj_factor + 2.0 * lambda[3];
            vals[1] = 0.0;
            vals[2] = 4.0 * obj_factor;
            vals[3] = 6.0 * obj_factor;
            vals[4] = 0.0;
            vals[5] = 0.0;
            true
        }
    }

    struct KnownAuxMultiplierProblem {
        bound_active_aux: bool,
        objective_slope: f64,
    }

    impl NlpProblem for KnownAuxMultiplierProblem {
        fn num_variables(&self) -> usize {
            2
        }

        fn num_constraints(&self) -> usize {
            1
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            x_l[0] = if self.bound_active_aux {
                2.0
            } else {
                f64::NEG_INFINITY
            };
            x_u[0] = if self.bound_active_aux {
                2.0
            } else {
                f64::INFINITY
            };
            x_l[1] = f64::NEG_INFINITY;
            x_u[1] = f64::INFINITY;
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l[0] = 0.0;
            g_u[0] = 0.0;
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0[0] = 2.0;
            x0[1] = 5.0;
        }

        fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = self.objective_slope * x[0] + (x[1] - 5.0) * (x[1] - 5.0);
            true
        }

        fn gradient(&self, x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad[0] = self.objective_slope;
            grad[1] = 2.0 * (x[1] - 5.0);
            true
        }

        fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g[0] = x[0] - 2.0;
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0], vec![0])
        }

        fn jacobian_values(&self, _x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals[0] = 1.0;
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![1], vec![1])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            obj_factor: f64,
            _lambda: &[f64],
            vals: &mut [f64],
        ) -> bool {
            vals[0] = 2.0 * obj_factor;
            true
        }
    }

    struct RankDeficientAuxProblem;

    impl NlpProblem for RankDeficientAuxProblem {
        fn num_variables(&self) -> usize {
            2
        }

        fn num_constraints(&self) -> usize {
            2
        }

        fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
            x_l.fill(f64::NEG_INFINITY);
            x_u.fill(f64::INFINITY);
        }

        fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
            g_l[0] = 0.0;
            g_u[0] = 0.0;
            g_l[1] = f64::NEG_INFINITY;
            g_u[1] = 10.0;
        }

        fn initial_point(&self, x0: &mut [f64]) {
            x0[0] = 0.0;
            x0[1] = 1.0;
        }

        fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = x[1] * x[1];
            true
        }

        fn gradient(&self, x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad[0] = 0.0;
            grad[1] = 2.0 * x[1];
            true
        }

        fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
            g[0] = x[0] * x[0];
            g[1] = x[1];
            true
        }

        fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0, 1], vec![0, 1])
        }

        fn jacobian_values(&self, x: &[f64], _new_x: bool, vals: &mut [f64]) -> bool {
            vals[0] = 2.0 * x[0];
            vals[1] = 1.0;
            true
        }

        fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
            (vec![0, 1], vec![0, 1])
        }

        fn hessian_values(
            &self,
            _x: &[f64],
            _new_x: bool,
            obj_factor: f64,
            lambda: &[f64],
            vals: &mut [f64],
        ) -> bool {
            vals[0] = 2.0 * lambda[0];
            vals[1] = 2.0 * obj_factor;
            true
        }
    }

    fn reduced_mapping_problem<'a>(problem: &'a dyn NlpProblem) -> AuxiliaryReducedProblem<'a> {
        let candidates = vec![PresolveCandidate {
            blocks: vec![
                EqualityBlock {
                    rows: vec![0],
                    vars: vec![1],
                },
                EqualityBlock {
                    rows: vec![2],
                    vars: vec![3],
                },
            ],
        }];
        AuxiliaryReducedProblem::new(problem, &candidates, vec![9.0, 2.0, 8.0, 5.0]).unwrap()
    }

    #[test]
    fn auxiliary_reduced_problem_expands_evaluations_through_fixed_auxiliaries() {
        let problem = ReducedMappingProblem;
        let reduced = reduced_mapping_problem(&problem);

        assert!(reduced.did_reduce());
        assert_eq!(reduced.num_fixed(), 2);
        assert_eq!(reduced.num_removed_constraints(), 2);
        assert_eq!(reduced.num_variables(), 2);
        assert_eq!(reduced.num_constraints(), 2);
        assert_eq!(reduced.var_map, vec![0, 2]);
        assert_eq!(reduced.constr_map, vec![1, 3]);
        assert_eq!(
            reduced.reduced_x_scaling(&[1.0, 2.0, 3.0, 4.0]),
            Some(vec![1.0, 3.0])
        );
        assert_eq!(
            reduced.reduced_g_scaling(&[10.0, 20.0, 30.0, 40.0]),
            Some(vec![20.0, 40.0])
        );

        let mut x0 = vec![0.0; 2];
        reduced.initial_point(&mut x0);
        assert_eq!(x0, vec![9.0, 8.0]);

        let mut x_l = vec![0.0; 2];
        let mut x_u = vec![0.0; 2];
        reduced.bounds(&mut x_l, &mut x_u);
        assert_eq!(x_l, vec![-10.0, -20.0]);
        assert_eq!(x_u, vec![10.0, 20.0]);

        let mut g_l = vec![0.0; 2];
        let mut g_u = vec![0.0; 2];
        reduced.constraint_bounds(&mut g_l, &mut g_u);
        assert_eq!(g_l, vec![f64::NEG_INFINITY, 1.0]);
        assert_eq!(g_u, vec![100.0, 1.0]);

        let x_reduced = vec![11.0, 13.0];
        let mut obj = 0.0;
        assert!(reduced.objective(&x_reduced, true, &mut obj));
        assert_eq!(obj, 1245.0);

        let mut grad = vec![0.0; 2];
        assert!(reduced.gradient(&x_reduced, true, &mut grad));
        assert_eq!(grad, vec![74.0, 122.0]);

        let mut g = vec![0.0; 2];
        assert!(reduced.constraints(&x_reduced, true, &mut g));
        assert_eq!(g, vec![57.0, 119.0]);
    }

    #[test]
    fn auxiliary_reduced_problem_remaps_jacobian_entries() {
        let problem = ReducedMappingProblem;
        let reduced = reduced_mapping_problem(&problem);

        let (rows, cols) = reduced.jacobian_structure();
        assert_eq!(rows, vec![0, 0, 1, 1]);
        assert_eq!(cols, vec![0, 1, 0, 1]);

        let mut vals = vec![0.0; rows.len()];
        assert!(reduced.jacobian_values(&[11.0, 13.0], true, &mut vals));
        assert_eq!(vals, vec![1.0, 2.0, 22.0, 1.0]);
    }

    #[test]
    fn auxiliary_reduced_problem_remaps_sparse_lower_hessian_entries() {
        let problem = ReducedMappingProblem;
        let reduced = reduced_mapping_problem(&problem);

        let (rows, cols) = reduced.hessian_structure();
        assert_eq!(rows, vec![0, 1, 1]);
        assert_eq!(cols, vec![0, 0, 1]);

        let mut vals = vec![0.0; rows.len()];
        assert!(reduced.hessian_values(&[11.0, 13.0], true, 0.5, &[7.0, 11.0], &mut vals));
        assert_eq!(vals, vec![23.0, 2.0, 3.0]);
    }

    #[test]
    fn auxiliary_reduced_problem_unmaps_full_solution() {
        let problem = ReducedMappingProblem;
        let reduced = reduced_mapping_problem(&problem);
        let reduced_result = SolveResult {
            x: vec![11.0, 13.0],
            objective: -1.0,
            constraint_multipliers: vec![7.0, 11.0],
            bound_multipliers_lower: vec![0.1, 0.2],
            bound_multipliers_upper: vec![0.3, 0.4],
            constraint_values: vec![-1.0, -1.0],
            status: SolveStatus::Optimal,
            iterations: 12,
            diagnostics: Default::default(),
        };

        let full = reduced.unmap_solution(&reduced_result);

        assert_eq!(full.x, vec![11.0, 2.0, 13.0, 5.0]);
        assert_eq!(full.objective, 1245.0);
        assert_eq!(full.constraint_values, vec![0.0, 57.0, 0.0, 119.0]);
        assert_eq!(full.constraint_multipliers, vec![0.0, 7.0, 0.0, 11.0]);
        assert_eq!(full.bound_multipliers_lower, vec![0.1, 0.0, 0.2, 0.0]);
        assert_eq!(full.bound_multipliers_upper, vec![0.3, 0.0, 0.4, 0.0]);
        assert_eq!(full.status, SolveStatus::Optimal);
        assert_eq!(full.iterations, 12);
    }

    #[test]
    fn auxiliary_reduced_problem_reconstructs_removed_multiplier() {
        let problem = KnownAuxMultiplierProblem {
            bound_active_aux: false,
            objective_slope: 3.0,
        };
        let candidates = vec![PresolveCandidate {
            blocks: vec![EqualityBlock {
                rows: vec![0],
                vars: vec![0],
            }],
        }];
        let reduced = AuxiliaryReducedProblem::new(&problem, &candidates, vec![2.0, 5.0])
            .expect("auxiliary reduction");
        let reduced_result = SolveResult {
            x: vec![5.0],
            objective: 0.0,
            constraint_multipliers: vec![],
            bound_multipliers_lower: vec![0.0],
            bound_multipliers_upper: vec![0.0],
            constraint_values: vec![],
            status: SolveStatus::Optimal,
            iterations: 1,
            diagnostics: Default::default(),
        };

        let full = reduced.unmap_solution(&reduced_result);

        assert_eq!(full.x, vec![2.0, 5.0]);
        assert_eq!(full.constraint_values, vec![0.0]);
        assert!(
            (full.constraint_multipliers[0] + 3.0).abs() < 1e-10,
            "lambda_aux={}, expected -3",
            full.constraint_multipliers[0]
        );
    }

    #[test]
    fn auxiliary_reduced_problem_reconstructs_large_removed_multiplier() {
        let problem = KnownAuxMultiplierProblem {
            bound_active_aux: false,
            objective_slope: 1.0e12,
        };
        let candidates = vec![PresolveCandidate {
            blocks: vec![EqualityBlock {
                rows: vec![0],
                vars: vec![0],
            }],
        }];
        let reduced = AuxiliaryReducedProblem::new(&problem, &candidates, vec![2.0, 5.0])
            .expect("auxiliary reduction");
        let reduced_result = SolveResult {
            x: vec![5.0],
            objective: 0.0,
            constraint_multipliers: vec![],
            bound_multipliers_lower: vec![0.0],
            bound_multipliers_upper: vec![0.0],
            constraint_values: vec![],
            status: SolveStatus::Optimal,
            iterations: 1,
            diagnostics: Default::default(),
        };

        let full = reduced.unmap_solution(&reduced_result);

        let expected = -1.0e12;
        assert!(
            (full.constraint_multipliers[0] - expected).abs() <= expected.abs() * 1e-12,
            "lambda_aux={}, expected {expected}",
            full.constraint_multipliers[0]
        );
    }

    #[test]
    fn auxiliary_reduced_problem_skips_multiplier_reconstruction_for_bound_active_auxiliary() {
        let problem = KnownAuxMultiplierProblem {
            bound_active_aux: true,
            objective_slope: 3.0,
        };
        let candidates = vec![PresolveCandidate {
            blocks: vec![EqualityBlock {
                rows: vec![0],
                vars: vec![0],
            }],
        }];
        let reduced = AuxiliaryReducedProblem::new(&problem, &candidates, vec![2.0, 5.0])
            .expect("auxiliary reduction");
        let reduced_result = SolveResult {
            x: vec![5.0],
            objective: 0.0,
            constraint_multipliers: vec![],
            bound_multipliers_lower: vec![0.0],
            bound_multipliers_upper: vec![0.0],
            constraint_values: vec![],
            status: SolveStatus::Optimal,
            iterations: 1,
            diagnostics: Default::default(),
        };

        let full = reduced.unmap_solution(&reduced_result);

        assert_eq!(full.x, vec![2.0, 5.0]);
        assert_eq!(full.constraint_values, vec![0.0]);
        assert_eq!(
            full.constraint_multipliers[0], 0.0,
            "bound-active auxiliary variables should leave removed multipliers conservative"
        );
    }

    #[test]
    fn auxiliary_multiplier_reconstruction_rejects_singular_system() {
        let matrix = vec![1.0, 2.0, 2.0, 4.0];
        let rhs = vec![1.0, 2.0];

        let result = solve_dense_square_system(matrix, rhs, 2);

        assert!(
            result.is_err(),
            "singular multiplier systems should not be reconstructed"
        );
    }

    #[test]
    fn auxiliary_reduced_problem_rejects_rank_deficient_auxiliary_block() {
        let problem = RankDeficientAuxProblem;
        let candidates = vec![PresolveCandidate {
            blocks: vec![EqualityBlock {
                rows: vec![0],
                vars: vec![0],
            }],
        }];

        let err = match AuxiliaryReducedProblem::new(&problem, &candidates, vec![0.0, 1.0]) {
            Ok(_) => panic!("rank-deficient auxiliary block should not reduce"),
            Err(err) => err,
        };

        match err {
            AuxiliarySolveError::RankDeficientBlock {
                block,
                rank,
                expected,
            } => {
                assert_eq!(block.rows, vec![0]);
                assert_eq!(block.vars, vec![0]);
                assert_eq!(rank, 0);
                assert_eq!(expected, 1);
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

    #[test]
    fn auxiliary_solve_handles_one_variable_nonlinear_equality() {
        let problem = TriangularAuxProblem;
        let candidates = vec![PresolveCandidate {
            blocks: vec![EqualityBlock {
                rows: vec![0],
                vars: vec![0],
            }],
        }];
        let options = quiet_aux_options();

        let outcome =
            solve_auxiliary_blocks(&problem, &candidates, &options, std::time::Instant::now())
                .expect("auxiliary solve");

        assert_eq!(outcome.blocks_solved, 1);
        assert!(outcome.max_residual <= options.auxiliary_tol);
        assert!((outcome.x[0] - 2.0).abs() < 1e-5, "x = {:?}", outcome.x);
        assert_eq!(outcome.x[1], 0.0);
    }

    #[test]
    fn auxiliary_solve_updates_full_vector_between_triangular_blocks() {
        let problem = TriangularAuxProblem;
        let candidates = find_presolve_candidates(&problem, TOL);
        let options = quiet_aux_options();

        let outcome =
            solve_auxiliary_blocks(&problem, &candidates, &options, std::time::Instant::now())
                .expect("auxiliary solve");

        assert_eq!(outcome.blocks_solved, 2);
        assert!(outcome.max_residual <= options.auxiliary_tol);
        assert!((outcome.x[0] - 2.0).abs() < 1e-5, "x = {:?}", outcome.x);
        assert!((outcome.x[1] - 3.0).abs() < 1e-5, "x = {:?}", outcome.x);
    }

    #[test]
    fn auxiliary_solve_failure_returns_structured_failure() {
        let problem = InfeasibleBoundAuxProblem;
        let candidates = find_presolve_candidates(&problem, TOL);
        let options = quiet_aux_options();

        let err =
            solve_auxiliary_blocks(&problem, &candidates, &options, std::time::Instant::now())
                .unwrap_err();

        match err {
            AuxiliarySolveError::BlockSolveFailed {
                block,
                status: _,
                residual,
            } => {
                assert_eq!(block.rows, vec![0]);
                assert_eq!(block.vars, vec![0]);
                assert!(residual > options.auxiliary_tol || !residual.is_finite());
            }
            other => panic!("unexpected auxiliary error: {:?}", other),
        }
    }

    #[test]
    fn auxiliary_solve_stops_when_outer_wall_time_is_exhausted() {
        let problem = TriangularAuxProblem;
        let candidates = find_presolve_candidates(&problem, TOL);
        let mut options = quiet_aux_options();
        options.max_wall_time = 0.01;
        let expired_start = std::time::Instant::now() - std::time::Duration::from_secs(1);

        let err =
            solve_auxiliary_blocks(&problem, &candidates, &options, expired_start).unwrap_err();

        assert_eq!(
            err,
            AuxiliarySolveError::TimeBudgetExceeded { blocks_solved: 0 }
        );
    }

    #[test]
    fn auxiliary_block_problem_respects_original_variable_bounds() {
        let problem = InfeasibleBoundAuxProblem;
        let block = EqualityBlock {
            rows: vec![0],
            vars: vec![0],
        };
        let fixed_x = vec![0.5];
        let aux = AuxiliaryBlockProblem::new(&problem, &block, &fixed_x).unwrap();

        let mut x_l = vec![0.0; 1];
        let mut x_u = vec![0.0; 1];
        aux.bounds(&mut x_l, &mut x_u);
        assert_eq!(x_l, vec![0.0]);
        assert_eq!(x_u, vec![1.0]);

        let result = crate::solve(&aux, &quiet_aux_options());
        assert!(
            result.x[0] >= -1e-10 && result.x[0] <= 1.0 + 1e-10,
            "auxiliary solution violated bounds: {:?}",
            result.x
        );
    }

    #[test]
    fn incidence_handles_no_constraints() {
        let problem = graph_problem(2, &[], &[]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert_eq!(inc.n_vars, 2);
        assert_eq!(inc.m_orig, 0);
        assert!(inc.row_global.is_empty());
        assert!(inc.row_local_for_global.is_empty());
        assert!(inc.row_adj_vars.is_empty());
        assert_eq!(inc.var_adj_rows, vec![Vec::<usize>::new(), Vec::new()]);
    }

    #[test]
    fn incidence_ignores_inequality_rows() {
        let problem = graph_problem(
            3,
            &[(0.0, 1.0), (2.0, f64::INFINITY), (f64::NEG_INFINITY, 0.0)],
            &[(0, 0), (1, 1), (2, 2)],
        );
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert!(inc.row_global.is_empty());
        assert_eq!(inc.row_local_for_global, vec![None, None, None]);
        assert!(inc.row_adj_vars.is_empty());
        assert_eq!(
            inc.var_adj_rows,
            vec![Vec::<usize>::new(), Vec::new(), Vec::new()]
        );
    }

    #[test]
    fn incidence_detects_equality_rows_within_tolerance() {
        let problem = graph_problem(
            3,
            &[(1.0, 1.0 + TOL * 0.5), (0.0, 1.0), (2.0, 2.0)],
            &[(0, 0), (1, 1), (2, 2)],
        );
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert_eq!(inc.row_global, vec![0, 2]);
        assert_eq!(inc.row_local_for_global, vec![Some(0), None, Some(1)]);
        assert_eq!(inc.row_adj_vars, vec![vec![0], vec![2]]);
        assert_eq!(inc.var_adj_rows, vec![vec![0], Vec::new(), vec![1]]);
    }

    #[test]
    fn incidence_deduplicates_and_sorts_structural_edges() {
        let problem = graph_problem(
            4,
            &[(0.0, 0.0), (1.0, 1.0)],
            &[(0, 3), (0, 1), (0, 3), (1, 2), (1, 0), (1, 2)],
        );
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert_eq!(inc.row_global, vec![0, 1]);
        assert_eq!(inc.row_adj_vars, vec![vec![1, 3], vec![0, 2]]);
        assert_eq!(inc.var_adj_rows, vec![vec![1], vec![0], vec![1], vec![0]]);
    }

    #[test]
    fn incidence_keeps_empty_equality_rows() {
        let problem = graph_problem(2, &[(3.0, 3.0), (0.0, 1.0)], &[(1, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert_eq!(inc.row_global, vec![0]);
        assert_eq!(inc.row_adj_vars, vec![Vec::<usize>::new()]);
        assert_eq!(inc.var_adj_rows, vec![Vec::<usize>::new(), Vec::new()]);
    }

    #[test]
    fn incidence_ignores_out_of_range_structure_entries() {
        let problem = graph_problem(
            3,
            &[(0.0, 0.0), (1.0, 1.0)],
            &[(0, 0), (99, 1), (1, 99), (1, 2)],
        );
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        assert_eq!(inc.row_adj_vars, vec![vec![0], vec![2]]);
        assert_eq!(inc.var_adj_rows, vec![vec![0], Vec::new(), vec![1]]);
    }

    #[test]
    fn matching_and_dm_square_1x1() {
        let bounds = equality_bounds(1);
        let problem = graph_problem(1, &bounds, &[(0, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let matching = inc.maximum_matching();
        assert_eq!(matching.row_to_var, vec![Some(0)]);
        assert_eq!(matching.var_to_row, vec![Some(0)]);
        assert!(matching.unmatched_rows.is_empty());
        assert!(matching.unmatched_vars.is_empty());

        let dm = inc.dulmage_mendelsohn_partition();
        assert_eq!(dm.square_rows, vec![0]);
        assert_eq!(dm.square_vars, vec![0]);
        assert!(dm.overconstrained_rows.is_empty());
        assert!(dm.overconstrained_vars.is_empty());
        assert!(dm.underconstrained_rows.is_empty());
        assert!(dm.underconstrained_vars.is_empty());
        assert!(dm.unmatched_rows.is_empty());
        assert!(dm.unmatched_vars.is_empty());
    }

    #[test]
    fn matching_square_2x2_is_deterministic() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(2, &bounds, &[(0, 1), (1, 1), (0, 0), (1, 0), (0, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let matching = inc.maximum_matching();
        assert_eq!(matching.row_to_var, vec![Some(0), Some(1)]);
        assert_eq!(matching.var_to_row, vec![Some(0), Some(1)]);
        assert!(matching.unmatched_rows.is_empty());
        assert!(matching.unmatched_vars.is_empty());

        let dm = inc.dulmage_mendelsohn_partition();
        assert_eq!(dm.square_rows, vec![0, 1]);
        assert_eq!(dm.square_vars, vec![0, 1]);
        assert!(dm.overconstrained_rows.is_empty());
        assert!(dm.underconstrained_rows.is_empty());
    }

    #[test]
    fn matching_reroutes_along_augmenting_path() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(2, &bounds, &[(0, 0), (0, 1), (1, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let matching = inc.maximum_matching();
        assert_eq!(matching.row_to_var, vec![Some(1), Some(0)]);
        assert_eq!(matching.var_to_row, vec![Some(1), Some(0)]);
        assert!(matching.unmatched_rows.is_empty());
        assert!(matching.unmatched_vars.is_empty());
    }

    #[test]
    fn matching_handles_long_augmenting_path_iteratively() {
        let chain_len = 10_000;
        let bounds = equality_bounds(chain_len + 1);
        let mut edges = Vec::with_capacity(2 * chain_len + 1);
        edges.push((0, 0));
        edges.push((0, chain_len));
        for row in 1..chain_len {
            edges.push((row, row - 1));
            edges.push((row, row));
        }
        edges.push((chain_len, chain_len - 1));

        let problem = graph_problem(chain_len + 1, &bounds, &edges);
        let inc = EqualityIncidence::from_problem(&problem, TOL);
        let matching = inc.maximum_matching();

        assert!(matching.unmatched_rows.is_empty());
        assert!(matching.unmatched_vars.is_empty());
        assert_eq!(matching.row_to_var[0], Some(chain_len));
        for row in 1..chain_len {
            assert_eq!(matching.row_to_var[row], Some(row - 1));
        }
        assert_eq!(matching.row_to_var[chain_len], Some(chain_len - 1));
    }

    #[test]
    fn dm_underconstrained_1_row_2_variables() {
        let bounds = equality_bounds(1);
        let problem = graph_problem(2, &bounds, &[(0, 1), (0, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let dm = inc.dulmage_mendelsohn_partition();
        assert_eq!(dm.matching.row_to_var, vec![Some(0)]);
        assert_eq!(dm.matching.var_to_row, vec![Some(0), None]);
        assert!(dm.unmatched_rows.is_empty());
        assert_eq!(dm.unmatched_vars, vec![1]);
        assert!(dm.overconstrained_rows.is_empty());
        assert!(dm.overconstrained_vars.is_empty());
        assert!(dm.square_rows.is_empty());
        assert!(dm.square_vars.is_empty());
        assert_eq!(dm.underconstrained_rows, vec![0]);
        assert_eq!(dm.underconstrained_vars, vec![0, 1]);
    }

    #[test]
    fn dm_overconstrained_2_rows_1_variable() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(1, &bounds, &[(1, 0), (0, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let dm = inc.dulmage_mendelsohn_partition();
        assert_eq!(dm.matching.row_to_var, vec![Some(0), None]);
        assert_eq!(dm.matching.var_to_row, vec![Some(0)]);
        assert_eq!(dm.unmatched_rows, vec![1]);
        assert!(dm.unmatched_vars.is_empty());
        assert_eq!(dm.overconstrained_rows, vec![0, 1]);
        assert_eq!(dm.overconstrained_vars, vec![0]);
        assert!(dm.square_rows.is_empty());
        assert!(dm.square_vars.is_empty());
        assert!(dm.underconstrained_rows.is_empty());
        assert!(dm.underconstrained_vars.is_empty());
    }

    #[test]
    fn dm_mixed_square_over_and_under_blocks() {
        let bounds = equality_bounds(5);
        let problem = graph_problem(
            5,
            &bounds,
            &[(0, 0), (1, 1), (2, 2), (3, 2), (4, 4), (4, 3)],
        );
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let dm = inc.dulmage_mendelsohn_partition();
        assert_eq!(
            dm.matching.row_to_var,
            vec![Some(0), Some(1), Some(2), None, Some(3)]
        );
        assert_eq!(
            dm.matching.var_to_row,
            vec![Some(0), Some(1), Some(2), Some(4), None]
        );
        assert_eq!(dm.unmatched_rows, vec![3]);
        assert_eq!(dm.unmatched_vars, vec![4]);
        assert_eq!(dm.overconstrained_rows, vec![2, 3]);
        assert_eq!(dm.overconstrained_vars, vec![2]);
        assert_eq!(dm.square_rows, vec![0, 1]);
        assert_eq!(dm.square_vars, vec![0, 1]);
        assert_eq!(dm.underconstrained_rows, vec![4]);
        assert_eq!(dm.underconstrained_vars, vec![3, 4]);
    }

    #[test]
    fn connected_components_split_selected_independent_systems() {
        let bounds = equality_bounds(4);
        let problem = graph_problem(5, &bounds, &[(0, 0), (1, 2), (1, 1), (2, 3), (3, 4)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let components = inc.connected_components(&[3, 1, 0, 2], &[4, 0, 3, 2, 1]);

        assert_eq!(
            components,
            vec![
                EqualityBlock {
                    rows: vec![0],
                    vars: vec![0],
                },
                EqualityBlock {
                    rows: vec![1],
                    vars: vec![1, 2],
                },
                EqualityBlock {
                    rows: vec![2],
                    vars: vec![3],
                },
                EqualityBlock {
                    rows: vec![3],
                    vars: vec![4],
                },
            ]
        );
    }

    #[test]
    fn connected_components_return_global_equality_rows() {
        let problem = graph_problem(2, &[(0.0, 1.0), (2.0, 2.0), (3.0, 3.0)], &[(1, 0), (2, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let components = inc.connected_components(&[0, 1], &[0, 1]);

        assert_eq!(
            components,
            vec![
                EqualityBlock {
                    rows: vec![1],
                    vars: vec![0],
                },
                EqualityBlock {
                    rows: vec![2],
                    vars: vec![1],
                },
            ]
        );
    }

    #[test]
    fn btd_splits_independent_square_systems() {
        let bounds = equality_bounds(3);
        let problem = graph_problem(3, &bounds, &[(2, 2), (0, 0), (1, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let blocks = inc
            .block_triangular_decomposition(&[2, 0, 1], &[2, 0, 1])
            .unwrap();

        assert_eq!(
            blocks,
            vec![
                EqualityBlock {
                    rows: vec![0],
                    vars: vec![0],
                },
                EqualityBlock {
                    rows: vec![1],
                    vars: vec![1],
                },
                EqualityBlock {
                    rows: vec![2],
                    vars: vec![2],
                },
            ]
        );
    }

    #[test]
    fn btd_orders_triangular_system_upstream_to_downstream() {
        let bounds = equality_bounds(3);
        let problem = graph_problem(3, &bounds, &[(2, 2), (1, 1), (1, 0), (0, 0), (2, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let blocks = inc
            .block_triangular_decomposition(&[2, 1, 0], &[2, 1, 0])
            .unwrap();

        assert_eq!(
            blocks,
            vec![
                EqualityBlock {
                    rows: vec![0],
                    vars: vec![0],
                },
                EqualityBlock {
                    rows: vec![1],
                    vars: vec![1],
                },
                EqualityBlock {
                    rows: vec![2],
                    vars: vec![2],
                },
            ]
        );
    }

    #[test]
    fn btd_returns_cyclic_square_system_as_one_block() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(2, &bounds, &[(0, 1), (1, 0), (0, 0), (1, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let blocks = inc
            .block_triangular_decomposition(&[0, 1], &[0, 1])
            .unwrap();

        assert_eq!(
            blocks,
            vec![EqualityBlock {
                rows: vec![0, 1],
                vars: vec![0, 1],
            }]
        );
    }

    #[test]
    fn btd_rejects_non_square_input() {
        let bounds = equality_bounds(1);
        let problem = graph_problem(2, &bounds, &[(0, 0), (0, 1)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let err = inc
            .block_triangular_decomposition(&[0], &[0, 1])
            .unwrap_err();

        assert_eq!(
            err,
            BlockTriangularizationError::NonSquare { rows: 1, vars: 2 }
        );
    }

    #[test]
    fn btd_rejects_square_but_imperfectly_matched_input() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(2, &bounds, &[(0, 0), (1, 0)]);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let err = inc
            .block_triangular_decomposition(&[0, 1], &[0, 1])
            .unwrap_err();

        assert_eq!(
            err,
            BlockTriangularizationError::ImperfectMatching {
                unmatched_rows: vec![1],
                unmatched_vars: vec![1],
            }
        );
    }

    #[test]
    fn btd_handles_long_triangular_dependency_chain_iteratively() {
        let chain_len = 10_000;
        let bounds = equality_bounds(chain_len);
        let mut edges = Vec::with_capacity(2 * chain_len - 1);
        edges.push((0, 0));
        for row in 1..chain_len {
            edges.push((row, row - 1));
            edges.push((row, row));
        }
        let problem = graph_problem(chain_len, &bounds, &edges);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let selected: Vec<_> = (0..chain_len).collect();
        let blocks = inc
            .block_triangular_decomposition(&selected, &selected)
            .unwrap();

        assert_eq!(blocks.len(), chain_len);
        for (idx, block) in blocks.iter().enumerate() {
            assert_eq!(block.rows, vec![idx]);
            assert_eq!(block.vars, vec![idx]);
        }
    }

    #[test]
    fn btd_handles_many_independent_blocks_with_heap_ready_set() {
        let block_count = 10_000;
        let bounds = equality_bounds(block_count);
        let edges: Vec<_> = (0..block_count).map(|idx| (idx, idx)).collect();
        let problem = graph_problem(block_count, &bounds, &edges);
        let inc = EqualityIncidence::from_problem(&problem, TOL);

        let selected: Vec<_> = (0..block_count).rev().collect();
        let blocks = inc
            .block_triangular_decomposition(&selected, &selected)
            .unwrap();

        assert_eq!(blocks.len(), block_count);
        for (idx, block) in blocks.iter().enumerate() {
            assert_eq!(block.rows, vec![idx]);
            assert_eq!(block.vars, vec![idx]);
        }
    }

    #[test]
    fn find_candidates_detects_independent_auxiliary_system() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(3, &bounds, &[(0, 1), (1, 2), (1, 1)]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert_eq!(
            candidates,
            vec![PresolveCandidate {
                blocks: vec![
                    EqualityBlock {
                        rows: vec![0],
                        vars: vec![1],
                    },
                    EqualityBlock {
                        rows: vec![1],
                        vars: vec![2],
                    },
                ],
            }]
        );
    }

    #[test]
    fn find_candidates_rejects_underconstrained_equality_component() {
        let bounds = equality_bounds(1);
        let problem = graph_problem(2, &bounds, &[(0, 0), (0, 1)]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert!(candidates.is_empty());
    }

    #[test]
    fn find_candidates_rejects_overconstrained_equality_component() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(1, &bounds, &[(0, 0), (1, 0)]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert!(candidates.is_empty());
    }

    #[test]
    fn find_candidates_rejects_equality_component_coupled_to_unsolved_variable() {
        let bounds = equality_bounds(2);
        let problem = graph_problem(3, &bounds, &[(0, 0), (1, 0), (1, 1), (1, 2)]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert!(candidates.is_empty());
    }

    #[test]
    fn find_candidates_allows_candidate_variables_in_objective_terms() {
        let bounds = equality_bounds(1);
        let problem = graph_problem_with_objective(2, &bounds, &[(0, 1)], &[1]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert_eq!(
            candidates,
            vec![PresolveCandidate {
                blocks: vec![EqualityBlock {
                    rows: vec![0],
                    vars: vec![1],
                }],
            }]
        );
    }

    #[test]
    fn find_candidates_allows_candidate_variables_in_inequality_rows() {
        let problem = graph_problem(2, &[(0.0, 0.0), (0.0, 1.0)], &[(0, 1), (1, 1)]);

        let candidates = find_presolve_candidates(&problem, TOL);

        assert_eq!(
            candidates,
            vec![PresolveCandidate {
                blocks: vec![EqualityBlock {
                    rows: vec![0],
                    vars: vec![1],
                }],
            }]
        );
    }

    #[test]
    fn find_candidates_returns_btd_ordered_blocks() {
        let problem = graph_problem(
            3,
            &[(0.0, 1.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)],
            &[(1, 0), (2, 0), (2, 1), (3, 1), (3, 2)],
        );

        let candidates = find_presolve_candidates(&problem, TOL);

        assert_eq!(
            candidates,
            vec![PresolveCandidate {
                blocks: vec![
                    EqualityBlock {
                        rows: vec![1],
                        vars: vec![0],
                    },
                    EqualityBlock {
                        rows: vec![2],
                        vars: vec![1],
                    },
                    EqualityBlock {
                        rows: vec![3],
                        vars: vec![2],
                    },
                ],
            }]
        );
    }
}
