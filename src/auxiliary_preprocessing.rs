//! Internal auxiliary-system preprocessing utilities.
//!
//! This module is intentionally crate-private. The auxiliary preprocessor is an
//! implementation detail of `enable_preprocessing`; it must not expose a public
//! decomposition or transform API.

use crate::problem::NlpProblem;
use std::collections::VecDeque;

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
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BipartiteMatching {
    pub(crate) row_to_var: Vec<Option<usize>>,
    pub(crate) var_to_row: Vec<Option<usize>>,
    pub(crate) unmatched_rows: Vec<usize>,
    pub(crate) unmatched_vars: Vec<usize>,
}

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DulmageMendelsohnPartition {
    pub(crate) matching: BipartiteMatching,
    pub(crate) overconstrained_rows: Vec<usize>,
    pub(crate) overconstrained_vars: Vec<usize>,
    pub(crate) square_rows: Vec<usize>,
    pub(crate) square_vars: Vec<usize>,
    pub(crate) underconstrained_rows: Vec<usize>,
    pub(crate) underconstrained_vars: Vec<usize>,
    pub(crate) unmatched_rows: Vec<usize>,
    pub(crate) unmatched_vars: Vec<usize>,
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
    for &var in &row_adj_vars[row] {
        if let Some(next_row) = var_to_row[var] {
            if dist[next_row] != dist[row] + 1 {
                continue;
            }
            if !matching_dfs(next_row, row_adj_vars, row_to_var, var_to_row, dist) {
                continue;
            }
        }

        row_to_var[row] = Some(var);
        var_to_row[var] = Some(row);
        return true;
    }

    dist[row] = usize::MAX;
    false
}

fn indices_where(flags: &[bool]) -> Vec<usize> {
    flags
        .iter()
        .enumerate()
        .filter_map(|(idx, &flag)| flag.then_some(idx))
        .collect()
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

    fn equality_bounds(n_rows: usize) -> Vec<(f64, f64)> {
        vec![(0.0, 0.0); n_rows]
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
}
