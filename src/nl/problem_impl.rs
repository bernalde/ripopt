use std::collections::HashMap;

use super::autodiff::Tape;
use super::parser::NlFileData;
use crate::NlpProblem;

/// An NLP problem parsed from an NL file, implementing the NlpProblem trait.
pub struct NlProblem {
    n: usize,
    m: usize,
    x_l: Vec<f64>,
    x_u: Vec<f64>,
    g_l: Vec<f64>,
    g_u: Vec<f64>,
    x0: Vec<f64>,
    maximize: bool,

    /// Tape for the nonlinear part of the objective.
    obj_tape: Option<Tape>,
    /// Linear coefficients for the objective: (var_idx, coeff).
    obj_linear: Vec<(usize, f64)>,

    /// Tape for the nonlinear part of each constraint.
    con_tapes: Vec<Option<Tape>>,
    /// Linear coefficients for each constraint.
    con_linear: Vec<Vec<(usize, f64)>>,

    /// Jacobian sparsity: (row_indices, col_indices).
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    /// Map from (row, col) to position in Jacobian values array.
    jac_map: HashMap<(usize, usize), usize>,

    /// Hessian sparsity (dense lower triangle).
    hess_rows: Vec<usize>,
    hess_cols: Vec<usize>,
}

impl NlProblem {
    /// Build an NlProblem from parsed NL file data.
    pub fn from_nl_data(data: NlFileData) -> Self {
        let n = data.header.n_vars;
        let m = data.header.n_constrs;

        // Build tapes for objective
        let (obj_idx, maximize, obj_expr) = data
            .obj_exprs
            .into_iter()
            .next()
            .unwrap_or((0, false, None));

        let obj_tape = obj_expr.map(|expr| Tape::build(&expr, &data.common_exprs, n));

        let obj_linear = if obj_idx < data.obj_linear.len() {
            data.obj_linear[obj_idx].clone()
        } else {
            Vec::new()
        };

        // Build tapes for constraints
        let con_tapes: Vec<Option<Tape>> = data
            .con_exprs
            .iter()
            .map(|expr| expr.as_ref().map(|e| Tape::build(e, &data.common_exprs, n)))
            .collect();

        // Build Jacobian sparsity from con_linear entries + nonlinear variables
        let mut jac_entries: Vec<(usize, usize)> = Vec::new();
        let mut jac_set: HashMap<(usize, usize), usize> = HashMap::new();

        for (i, linear) in data.con_linear.iter().enumerate() {
            for &(var_idx, _) in linear {
                let key = (i, var_idx);
                if !jac_set.contains_key(&key) {
                    let pos = jac_entries.len();
                    jac_set.insert(key, pos);
                    jac_entries.push(key);
                }
            }
        }

        // Also add entries for nonlinear variables in each constraint
        for (i, tape) in con_tapes.iter().enumerate() {
            if let Some(tape) = tape {
                for op in &tape.ops {
                    if let super::autodiff::TapeOp::Var(j) = op {
                        let key = (i, *j);
                        if !jac_set.contains_key(&key) {
                            let pos = jac_entries.len();
                            jac_set.insert(key, pos);
                            jac_entries.push(key);
                        }
                    }
                }
            }
        }

        let jac_rows: Vec<usize> = jac_entries.iter().map(|&(r, _)| r).collect();
        let jac_cols: Vec<usize> = jac_entries.iter().map(|&(_, c)| c).collect();

        // Dense lower-triangle Hessian sparsity
        let mut hess_rows = Vec::with_capacity(n * (n + 1) / 2);
        let mut hess_cols = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in 0..=i {
                hess_rows.push(i);
                hess_cols.push(j);
            }
        }

        NlProblem {
            n,
            m,
            x_l: data.x_l,
            x_u: data.x_u,
            g_l: data.g_l,
            g_u: data.g_u,
            x0: data.x0,
            maximize,
            obj_tape,
            obj_linear,
            con_tapes,
            con_linear: data.con_linear,
            jac_rows,
            jac_cols,
            jac_map: jac_set,
            hess_rows,
            hess_cols,
        }
    }

    /// Compute the gradient of the objective (nonlinear + linear).
    fn obj_gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad.iter_mut().for_each(|v| *v = 0.0);

        // Nonlinear part via reverse AD
        if let Some(tape) = &self.obj_tape {
            tape.gradient(x, grad);
        }

        // Linear part
        for &(idx, coeff) in &self.obj_linear {
            if idx < grad.len() {
                grad[idx] += coeff;
            }
        }

        // Negate if maximizing (solver minimizes)
        if self.maximize {
            for g in grad.iter_mut() {
                *g = -*g;
            }
        }
    }

    /// Compute the gradient of constraint i (nonlinear + linear).
    fn con_gradient(&self, i: usize, x: &[f64], grad: &mut [f64]) {
        grad.iter_mut().for_each(|v| *v = 0.0);

        // Nonlinear part via reverse AD
        if let Some(tape) = &self.con_tapes[i] {
            tape.gradient(x, grad);
        }

        // Linear part
        for &(idx, coeff) in &self.con_linear[i] {
            if idx < grad.len() {
                grad[idx] += coeff;
            }
        }
    }

    /// Compute the Lagrangian gradient: obj_factor * ∇f + Σ lambda[i] * ∇g_i.
    fn lagrangian_gradient(&self, x: &[f64], obj_factor: f64, lambda: &[f64], grad: &mut [f64]) {
        grad.iter_mut().for_each(|v| *v = 0.0);

        // Objective gradient
        if obj_factor != 0.0 {
            let mut obj_grad = vec![0.0; self.n];
            self.obj_gradient(x, &mut obj_grad);
            for i in 0..self.n {
                grad[i] += obj_factor * obj_grad[i];
            }
        }

        // Constraint gradients
        let mut con_grad = vec![0.0; self.n];
        for (j, &lam) in lambda.iter().enumerate() {
            if lam == 0.0 {
                continue;
            }
            con_grad.iter_mut().for_each(|v| *v = 0.0);
            self.con_gradient(j, x, &mut con_grad);
            for i in 0..self.n {
                grad[i] += lam * con_grad[i];
            }
        }
    }
}

impl NlpProblem for NlProblem {
    fn num_variables(&self) -> usize {
        self.n
    }

    fn num_constraints(&self) -> usize {
        self.m
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l.copy_from_slice(&self.x_l);
        x_u.copy_from_slice(&self.x_u);
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l.copy_from_slice(&self.g_l);
        g_u.copy_from_slice(&self.g_u);
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0.copy_from_slice(&self.x0);
    }

    fn objective(&self, x: &[f64]) -> f64 {
        let mut val = 0.0;

        // Nonlinear part
        if let Some(tape) = &self.obj_tape {
            val += tape.eval(x);
        }

        // Linear part
        for &(idx, coeff) in &self.obj_linear {
            val += coeff * x[idx];
        }

        if self.maximize {
            -val
        } else {
            val
        }
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        self.obj_gradient(x, grad);
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        for i in 0..self.m {
            let mut val = 0.0;

            // Nonlinear part
            if let Some(tape) = &self.con_tapes[i] {
                val += tape.eval(x);
            }

            // Linear part
            for &(idx, coeff) in &self.con_linear[i] {
                val += coeff * x[idx];
            }

            g[i] = val;
        }
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.jac_rows.clone(), self.jac_cols.clone())
    }

    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        vals.iter_mut().for_each(|v| *v = 0.0);

        let mut grad = vec![0.0; self.n];

        for i in 0..self.m {
            self.con_gradient(i, x, &mut grad);

            // Scatter into vals using jac_map
            for j in 0..self.n {
                if grad[j] != 0.0 {
                    if let Some(&pos) = self.jac_map.get(&(i, j)) {
                        vals[pos] = grad[j];
                    }
                }
            }
        }
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hess_rows.clone(), self.hess_cols.clone())
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]) {
        // Finite-difference Hessian of the Lagrangian.
        // H[:,j] ≈ (∇L(x + h*e_j) - ∇L(x - h*e_j)) / (2h)
        let n = self.n;
        vals.iter_mut().for_each(|v| *v = 0.0);

        let mut x_plus = x.to_vec();
        let mut x_minus = x.to_vec();
        let mut grad_plus = vec![0.0; n];
        let mut grad_minus = vec![0.0; n];

        let mut idx = 0;
        for j in 0..n {
            let h = (1e-8_f64).max(x[j].abs() * 1e-8);

            x_plus[j] = x[j] + h;
            x_minus[j] = x[j] - h;

            self.lagrangian_gradient(&x_plus, obj_factor, lambda, &mut grad_plus);
            self.lagrangian_gradient(&x_minus, obj_factor, lambda, &mut grad_minus);

            let inv_2h = 0.5 / h;

            // Fill lower triangle column j: rows j..n
            // Dense lower triangle ordering: for each row i, columns 0..=i
            // We iterate (i,j) for i >= j, which for column j means i = j..n
            // But our hess layout is row-major lower triangle:
            //   (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
            // For column j, the entries are at positions where col == j.
            // We need to find those positions.
            for i in j..n {
                let hij = (grad_plus[i] - grad_minus[i]) * inv_2h;
                // Position in lower triangle: row i, col j → index = i*(i+1)/2 + j
                let pos = i * (i + 1) / 2 + j;
                vals[pos] = hij;
            }

            x_plus[j] = x[j];
            x_minus[j] = x[j];
            let _ = idx; // suppress warning
            idx += 1;
        }
    }
}
