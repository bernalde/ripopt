//! Internal auxiliary-system preprocessing utilities.
//!
//! This module is intentionally crate-private. The auxiliary preprocessor is an
//! implementation detail of `enable_preprocessing`; it must not expose a public
//! decomposition or transform API.

use crate::problem::NlpProblem;

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

        fn objective(&self, _x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
            *obj = 0.0;
            true
        }

        fn gradient(&self, _x: &[f64], _new_x: bool, grad: &mut [f64]) -> bool {
            grad.fill(0.0);
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
        }
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
}
