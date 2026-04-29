//! Internal auxiliary-system preprocessing utilities.
//!
//! This module is intentionally crate-private. The auxiliary preprocessor is an
//! implementation detail of `enable_preprocessing`; it must not expose a public
//! decomposition or transform API.

use crate::problem::NlpProblem;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};

#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EqualityBlock {
    pub(crate) rows: Vec<usize>,
    pub(crate) vars: Vec<usize>,
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
}
