use crate::csc::CscMatrix;
use std::collections::VecDeque;

// --- Constants ---

const ND_THRESHOLD: usize = 200;
const COARSEN_THRESHOLD: usize = 100;
const MAX_FM_PASSES: usize = 10;
const MAX_IMBALANCE: f64 = 0.4;
const MIN_COARSEN_RATIO: f64 = 0.8;

// --- Data structures ---

/// CSR-style adjacency graph (flat arrays, cache-friendly).
struct Graph {
    n: usize,
    xadj: Vec<usize>,
    adjncy: Vec<usize>,
}

/// 3-way vertex partition: left (0), right (1), separator (2).
struct Partition {
    part: Vec<u8>,
    sizes: [usize; 3],
}

/// Maps between coarsening levels.
struct CoarsenLevel {
    coarse_graph: Graph,
    fine_to_coarse: Vec<usize>,
    coarse_to_fine: Vec<Vec<usize>>,
}

// --- Graph construction ---

/// Build undirected adjacency graph from upper-triangle CSC matrix.
fn graph_from_csc(csc: &CscMatrix) -> Graph {
    let n = csc.n;
    if n == 0 {
        return Graph { n: 0, xadj: vec![0], adjncy: Vec::new() };
    }

    // Count degrees
    let mut deg = vec![0usize; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            if i != j {
                deg[i] += 1;
                deg[j] += 1;
            }
        }
    }

    // Build xadj
    let mut xadj = vec![0usize; n + 1];
    for i in 0..n {
        xadj[i + 1] = xadj[i] + deg[i];
    }
    let total = xadj[n];

    // Fill adjncy
    let mut adjncy = vec![0usize; total];
    let mut pos = vec![0usize; n];
    for j in 0..n {
        for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
            let i = csc.row_idx[idx];
            if i != j {
                adjncy[xadj[i] + pos[i]] = j;
                pos[i] += 1;
                adjncy[xadj[j] + pos[j]] = i;
                pos[j] += 1;
            }
        }
    }

    // Sort neighbor lists
    for i in 0..n {
        adjncy[xadj[i]..xadj[i + 1]].sort_unstable();
    }

    Graph { n, xadj, adjncy }
}

/// Build vertex-induced subgraph with remapped indices.
/// Returns (subgraph, old_to_new_map) where old_to_new_map[old] = new index (or usize::MAX if not included).
fn extract_subgraph(graph: &Graph, vertices: &[usize]) -> (Graph, Vec<usize>) {
    let mut old_to_new = vec![usize::MAX; graph.n];
    for (new, &old) in vertices.iter().enumerate() {
        old_to_new[old] = new;
    }

    let sub_n = vertices.len();
    let mut xadj = vec![0usize; sub_n + 1];
    let mut adjncy = Vec::new();

    for (new_v, &old_v) in vertices.iter().enumerate() {
        for k in graph.xadj[old_v]..graph.xadj[old_v + 1] {
            let nb = graph.adjncy[k];
            if old_to_new[nb] != usize::MAX {
                adjncy.push(old_to_new[nb]);
            }
        }
        xadj[new_v + 1] = adjncy.len();
    }

    (Graph { n: sub_n, xadj, adjncy }, old_to_new)
}

// --- BFS utilities ---

/// Find a pseudo-peripheral node via two BFS passes.
fn pseudo_peripheral_node(graph: &Graph) -> usize {
    if graph.n == 0 {
        return 0;
    }
    let start = 0;
    let far = bfs_farthest(graph, start);
    bfs_farthest(graph, far)
}

/// BFS from `start`, return the last node visited (farthest).
fn bfs_farthest(graph: &Graph, start: usize) -> usize {
    let mut visited = vec![false; graph.n];
    let mut queue = VecDeque::new();
    visited[start] = true;
    queue.push_back(start);
    let mut last = start;

    while let Some(v) = queue.pop_front() {
        last = v;
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let nb = graph.adjncy[k];
            if !visited[nb] {
                visited[nb] = true;
                queue.push_back(nb);
            }
        }
    }
    last
}

// --- Bisection ---

/// BFS greedy growing from pseudo-peripheral node.
/// First n/2 nodes visited → part 0, rest → part 1.
fn initial_bisection(graph: &Graph) -> Partition {
    let n = graph.n;
    let half = n / 2;

    let seed = pseudo_peripheral_node(graph);
    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut queue = VecDeque::new();

    visited[seed] = true;
    queue.push_back(seed);
    while let Some(v) = queue.pop_front() {
        order.push(v);
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let nb = graph.adjncy[k];
            if !visited[nb] {
                visited[nb] = true;
                queue.push_back(nb);
            }
        }
    }

    // Handle disconnected components: add unvisited nodes
    for v in 0..n {
        if !visited[v] {
            order.push(v);
        }
    }

    let mut part = vec![1u8; n];
    let mut sizes = [0usize; 3];
    for (i, &v) in order.iter().enumerate() {
        if i < half {
            part[v] = 0;
            sizes[0] += 1;
        } else {
            sizes[1] += 1;
        }
    }

    Partition { part, sizes }
}

// --- FM refinement ---

/// Fiduccia-Mattheyses refinement of a 2-way partition.
fn fm_refine(graph: &Graph, partition: &mut Partition, max_passes: usize) {
    let n = graph.n;
    if n <= 1 {
        return;
    }

    let mut gain = vec![0i64; n];
    let mut locked = vec![false; n];

    for pass in 0..max_passes {
        // Compute initial gains
        for v in 0..n {
            let my_side = partition.part[v];
            if my_side == 2 {
                gain[v] = 0;
                continue;
            }
            let mut g: i64 = 0;
            for k in graph.xadj[v]..graph.xadj[v + 1] {
                let nb = graph.adjncy[k];
                let nb_side = partition.part[nb];
                if nb_side == my_side {
                    g -= 1; // edge to same side (would lose this)
                } else if nb_side != 2 {
                    g += 1; // edge to other side (would gain this)
                }
            }
            gain[v] = g;
        }

        for l in &mut locked {
            *l = false;
        }

        let mut best_cut_improvement = 0i64;
        let mut cumulative_gain = 0i64;
        let mut moves: Vec<(usize, u8)> = Vec::new(); // (vertex, old_part)
        let mut best_prefix = 0usize; // number of moves to keep

        let max_moves = n;
        for _ in 0..max_moves {
            // Find best unlocked movable vertex
            let mut best_v = usize::MAX;
            let mut best_g = i64::MIN;

            for v in 0..n {
                if locked[v] || partition.part[v] == 2 {
                    continue;
                }
                // Check balance constraint
                let from = partition.part[v] as usize;
                let to = 1 - from;
                if partition.sizes[to] + 1 > ((n as f64) * (1.0 - MAX_IMBALANCE)).ceil() as usize {
                    continue;
                }
                if partition.sizes[from] <= 1 {
                    continue;
                }
                if gain[v] > best_g {
                    best_g = gain[v];
                    best_v = v;
                }
            }

            if best_v == usize::MAX {
                break;
            }

            // Move best_v
            let old_part = partition.part[best_v];
            let new_part = 1 - old_part;
            partition.part[best_v] = new_part;
            partition.sizes[old_part as usize] -= 1;
            partition.sizes[new_part as usize] += 1;
            locked[best_v] = true;
            moves.push((best_v, old_part));

            cumulative_gain += best_g;
            if cumulative_gain > best_cut_improvement {
                best_cut_improvement = cumulative_gain;
                best_prefix = moves.len();
            }

            // Update neighbor gains
            for k in graph.xadj[best_v]..graph.xadj[best_v + 1] {
                let nb = graph.adjncy[k];
                if locked[nb] || partition.part[nb] == 2 {
                    continue;
                }
                let nb_side = partition.part[nb];
                if nb_side == new_part {
                    // nb and best_v now same side: nb's gain decreases
                    gain[nb] -= 2;
                } else if nb_side != 2 {
                    // nb and best_v now different sides: nb's gain increases
                    gain[nb] += 2;
                }
            }
        }

        // Roll back to best prefix
        if best_prefix < moves.len() {
            for i in (best_prefix..moves.len()).rev() {
                let (v, old_part) = moves[i];
                let cur_part = partition.part[v] as usize;
                partition.part[v] = old_part;
                partition.sizes[cur_part] -= 1;
                partition.sizes[old_part as usize] += 1;
            }
        }

        if best_cut_improvement <= 0 || pass + 1 >= max_passes {
            break;
        }
    }
}

// --- Coarsening ---

/// Heavy-edge matching coarsening.
fn coarsen(graph: &Graph) -> CoarsenLevel {
    let n = graph.n;
    let mut matched = vec![false; n];
    let mut fine_to_coarse = vec![usize::MAX; n];
    let mut coarse_to_fine: Vec<Vec<usize>> = Vec::new();
    let mut coarse_n = 0usize;

    // Match vertices with their heaviest unmatched neighbor
    for v in 0..n {
        if matched[v] {
            continue;
        }
        // Find first unmatched neighbor
        let mut mate = usize::MAX;
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let nb = graph.adjncy[k];
            if !matched[nb] {
                mate = nb;
                break;
            }
        }

        if mate != usize::MAX {
            matched[v] = true;
            matched[mate] = true;
            fine_to_coarse[v] = coarse_n;
            fine_to_coarse[mate] = coarse_n;
            coarse_to_fine.push(vec![v, mate]);
        } else {
            matched[v] = true;
            fine_to_coarse[v] = coarse_n;
            coarse_to_fine.push(vec![v]);
        }
        coarse_n += 1;
    }

    // Build coarse graph
    let mut edge_set: Vec<Vec<usize>> = vec![Vec::new(); coarse_n];
    for v in 0..n {
        let cv = fine_to_coarse[v];
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let cu = fine_to_coarse[graph.adjncy[k]];
            if cu != cv {
                edge_set[cv].push(cu);
            }
        }
    }

    // Deduplicate and build flat adjacency
    let mut xadj = vec![0usize; coarse_n + 1];
    let mut adjncy = Vec::new();
    for cv in 0..coarse_n {
        edge_set[cv].sort_unstable();
        edge_set[cv].dedup();
        for &nb in &edge_set[cv] {
            adjncy.push(nb);
        }
        xadj[cv + 1] = adjncy.len();
    }

    CoarsenLevel {
        coarse_graph: Graph { n: coarse_n, xadj, adjncy },
        fine_to_coarse,
        coarse_to_fine,
    }
}

// --- Multilevel bisection ---

/// Multilevel bisection returning (separator, left, right) vertex lists.
fn multilevel_bisect(graph: &Graph) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    if graph.n <= 1 {
        if graph.n == 1 {
            return (vec![], vec![0], vec![]);
        }
        return (vec![], vec![], vec![]);
    }

    // Coarsen until small enough
    let mut levels: Vec<CoarsenLevel> = Vec::new();
    let mut current = graph;
    let mut owned_graphs: Vec<Graph> = Vec::new();

    loop {
        if current.n <= COARSEN_THRESHOLD {
            break;
        }
        let level = coarsen(current);
        let ratio = level.coarse_graph.n as f64 / current.n as f64;
        if ratio > MIN_COARSEN_RATIO {
            break;
        }
        owned_graphs.push(Graph {
            n: level.coarse_graph.n,
            xadj: level.coarse_graph.xadj.clone(),
            adjncy: level.coarse_graph.adjncy.clone(),
        });
        levels.push(CoarsenLevel {
            coarse_graph: Graph {
                n: owned_graphs.last().unwrap().n,
                xadj: Vec::new(), // placeholder
                adjncy: Vec::new(),
            },
            fine_to_coarse: level.fine_to_coarse,
            coarse_to_fine: level.coarse_to_fine,
        });
        current = owned_graphs.last().unwrap();
    }

    // Initial bisection on coarsest graph
    let mut partition = initial_bisection(current);
    fm_refine(current, &mut partition, MAX_FM_PASSES);

    // Uncoarsen: project partition and refine at each level
    for level_idx in (0..levels.len()).rev() {
        let fine_graph = if level_idx == 0 {
            graph
        } else {
            &owned_graphs[level_idx - 1]
        };
        let level = &levels[level_idx];
        let fine_n = fine_graph.n;

        // Project partition to finer level
        let mut fine_part = vec![0u8; fine_n];
        let mut fine_sizes = [0usize; 3];
        for fv in 0..fine_n {
            let cv = level.fine_to_coarse[fv];
            fine_part[fv] = partition.part[cv];
            fine_sizes[fine_part[fv] as usize] += 1;
        }

        partition = Partition { part: fine_part, sizes: fine_sizes };
        fm_refine(fine_graph, &mut partition, MAX_FM_PASSES);
    }

    // Convert edge partition to vertex separator
    edge_to_vertex_separator(graph, &partition)
}

/// Convert 2-way edge partition to vertex separator.
/// For each cut edge (u,v), add the lower-degree endpoint to separator.
fn edge_to_vertex_separator(graph: &Graph, partition: &Partition) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let n = graph.n;
    let mut is_sep = vec![false; n];

    // Find cut edges and add lower-degree endpoint to separator
    for v in 0..n {
        if partition.part[v] == 2 || is_sep[v] {
            continue;
        }
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let nb = graph.adjncy[k];
            if partition.part[nb] != partition.part[v] && partition.part[nb] != 2 && !is_sep[nb] && !is_sep[v] {
                let deg_v = graph.xadj[v + 1] - graph.xadj[v];
                let deg_nb = graph.xadj[nb + 1] - graph.xadj[nb];
                if deg_v <= deg_nb {
                    is_sep[v] = true;
                } else {
                    is_sep[nb] = true;
                }
            }
        }
    }

    let mut sep = Vec::new();
    let mut left = Vec::new();
    let mut right = Vec::new();

    for v in 0..n {
        if is_sep[v] || partition.part[v] == 2 {
            sep.push(v);
        } else if partition.part[v] == 0 {
            left.push(v);
        } else {
            right.push(v);
        }
    }

    // Ensure we actually have a valid separator — if both sides are empty, just split
    if left.is_empty() && right.is_empty() && sep.len() == n {
        // Degenerate case: all in separator. Split into halves with no separator.
        let half = n / 2;
        left = (0..half).collect();
        right = (half..n).collect();
        sep = Vec::new();
    }

    (sep, left, right)
}

// --- Recursive nested dissection ---

/// Recursive nested dissection core.
fn nd_recursive(graph: &Graph, global_ids: &[usize], perm: &mut Vec<usize>, threshold: usize) {
    let n = graph.n;

    if n == 0 {
        return;
    }

    // Base case: delegate to AMD for small subproblems
    if n <= threshold {
        // Build a local CSC from the subgraph for AMD
        let local_csc = graph_to_csc(graph);
        let (local_perm, _) = super::amd::amd_ordering(&local_csc);
        for &lp in &local_perm {
            perm.push(global_ids[lp]);
        }
        return;
    }

    // Multilevel bisection
    let (sep, left, right) = multilevel_bisect(graph);

    // If bisection fails to produce meaningful split, fall back to AMD
    if left.is_empty() || right.is_empty() {
        let local_csc = graph_to_csc(graph);
        let (local_perm, _) = super::amd::amd_ordering(&local_csc);
        for &lp in &local_perm {
            perm.push(global_ids[lp]);
        }
        return;
    }

    // Recurse on left partition
    {
        let left_global: Vec<usize> = left.iter().map(|&v| global_ids[v]).collect();
        let (left_graph, _) = extract_subgraph(graph, &left);
        nd_recursive(&left_graph, &left_global, perm, threshold);
    }

    // Recurse on right partition
    {
        let right_global: Vec<usize> = right.iter().map(|&v| global_ids[v]).collect();
        let (right_graph, _) = extract_subgraph(graph, &right);
        nd_recursive(&right_graph, &right_global, perm, threshold);
    }

    // Separator last (eliminated last → fill reduction)
    for &s in &sep {
        perm.push(global_ids[s]);
    }
}

/// Build a minimal CSC matrix from a Graph (for AMD fallback).
/// Creates upper-triangle only entries with 1.0 values on diagonal and off-diagonal.
fn graph_to_csc(graph: &Graph) -> CscMatrix {
    let n = graph.n;
    if n == 0 {
        return CscMatrix { n: 0, col_ptr: vec![0], row_idx: Vec::new(), vals: Vec::new() };
    }

    // Collect upper-triangle entries (including diagonal)
    let mut triplets: Vec<(usize, usize)> = Vec::new();
    for v in 0..n {
        triplets.push((v, v)); // diagonal
        for k in graph.xadj[v]..graph.xadj[v + 1] {
            let nb = graph.adjncy[k];
            if nb > v {
                triplets.push((v, nb));
            }
        }
    }

    // Sort by (col, row) for CSC
    triplets.sort_unstable_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

    // Build CSC
    let mut col_ptr = vec![0usize; n + 1];
    let mut row_idx = Vec::with_capacity(triplets.len());
    let mut vals = Vec::with_capacity(triplets.len());

    for &(r, c) in &triplets {
        col_ptr[c + 1] += 1;
        row_idx.push(r);
        vals.push(1.0);
    }
    for j in 0..n {
        col_ptr[j + 1] += col_ptr[j];
    }

    // Fix col_ptr: we accumulated counts, need cumulative
    // Actually, let's rebuild properly
    let mut col_ptr2 = vec![0usize; n + 1];
    for &(_, c) in &triplets {
        col_ptr2[c + 1] += 1;
    }
    for j in 0..n {
        col_ptr2[j + 1] += col_ptr2[j];
    }

    CscMatrix { n, col_ptr: col_ptr2, row_idx, vals }
}

// --- Public API ---

/// Nested dissection ordering for a symmetric CSC matrix.
/// Returns `(perm, perm_inv)` where `perm[new_pos] = old_pos`.
pub fn nd_ordering(csc: &CscMatrix) -> (Vec<usize>, Vec<usize>) {
    let n = csc.n;
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // For small matrices, just use AMD directly
    if n <= ND_THRESHOLD {
        return super::amd::amd_ordering(csc);
    }

    let graph = graph_from_csc(csc);
    let global_ids: Vec<usize> = (0..n).collect();
    let mut perm = Vec::with_capacity(n);

    nd_recursive(&graph, &global_ids, &mut perm, ND_THRESHOLD);

    // Safety: ensure valid permutation (should always be the case)
    debug_assert_eq!(perm.len(), n);

    let mut perm_inv = vec![0usize; n];
    for (new_pos, &old_idx) in perm.iter().enumerate() {
        perm_inv[old_idx] = new_pos;
    }

    (perm, perm_inv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coo::CooMatrix;
    use crate::csc::CscMatrix;

    fn csc_from_upper_triplets(n: usize, triplets: &[(usize, usize, f64)]) -> CscMatrix {
        let rows: Vec<usize> = triplets.iter().map(|t| t.0).collect();
        let cols: Vec<usize> = triplets.iter().map(|t| t.1).collect();
        let vals: Vec<f64> = triplets.iter().map(|t| t.2).collect();
        let coo = CooMatrix::new(n, rows, cols, vals).unwrap();
        CscMatrix::from_coo(&coo)
    }

    fn is_valid_permutation(perm: &[usize], n: usize) -> bool {
        if perm.len() != n {
            return false;
        }
        let mut seen = vec![false; n];
        for &p in perm {
            if p >= n || seen[p] {
                return false;
            }
            seen[p] = true;
        }
        true
    }

    fn count_fill(csc: &CscMatrix, perm: &[usize]) -> usize {
        let n = csc.n;
        let perm_inv: Vec<usize> = {
            let mut pi = vec![0; n];
            for (new, &old) in perm.iter().enumerate() {
                pi[old] = new;
            }
            pi
        };
        let permuted = crate::ordering::permute_symmetric_csc(csc, perm, &perm_inv);
        let sym = crate::symbolic::SymbolicFactorization::from_csc(&permuted);
        sym.l_nnz
    }

    fn build_2d_grid(nx: usize) -> CscMatrix {
        let n = nx * nx;
        let mut triplets = Vec::new();
        for iy in 0..nx {
            for ix in 0..nx {
                let idx = iy * nx + ix;
                triplets.push((idx, idx, 4.0));
                if ix + 1 < nx {
                    let right = iy * nx + ix + 1;
                    triplets.push((idx, right, -1.0));
                }
                if iy + 1 < nx {
                    let below = (iy + 1) * nx + ix;
                    triplets.push((idx, below, -1.0));
                }
            }
        }
        csc_from_upper_triplets(n, &triplets)
    }

    #[test]
    fn test_nd_valid_permutation_path_graph() {
        // Path graph: 0-1-2-3-4
        let csc = csc_from_upper_triplets(5, &[
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 1, 2.0), (1, 2, -1.0),
            (2, 2, 2.0), (2, 3, -1.0),
            (3, 3, 2.0), (3, 4, -1.0),
            (4, 4, 2.0),
        ]);
        let (perm, perm_inv) = nd_ordering(&csc);
        assert!(is_valid_permutation(&perm, 5));
        assert!(is_valid_permutation(&perm_inv, 5));
        // Check perm/perm_inv consistency
        for i in 0..5 {
            assert_eq!(perm_inv[perm[i]], i);
            assert_eq!(perm[perm_inv[i]], i);
        }
    }

    #[test]
    fn test_nd_empty() {
        let csc = CscMatrix { n: 0, col_ptr: vec![0], row_idx: vec![], vals: vec![] };
        let (perm, perm_inv) = nd_ordering(&csc);
        assert_eq!(perm.len(), 0);
        assert_eq!(perm_inv.len(), 0);
    }

    #[test]
    fn test_nd_single_node() {
        let csc = csc_from_upper_triplets(1, &[(0, 0, 1.0)]);
        let (perm, perm_inv) = nd_ordering(&csc);
        assert_eq!(perm, vec![0]);
        assert_eq!(perm_inv, vec![0]);
    }

    #[test]
    fn test_nd_tridiagonal() {
        let n = 10;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 3.0));
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);
        let (perm, perm_inv) = nd_ordering(&csc);
        assert!(is_valid_permutation(&perm, n));
        assert!(is_valid_permutation(&perm_inv, n));

        // Verify solve works with this ordering
        let permuted = crate::ordering::permute_symmetric_csc(&csc, &perm, &perm_inv);
        let sym = crate::symbolic::SymbolicFactorization::from_csc(&permuted);
        assert!(sym.l_nnz > 0);
    }

    #[test]
    fn test_nd_2d_grid_fill_vs_natural() {
        let csc = build_2d_grid(10);
        let n = 100;

        let natural_perm: Vec<usize> = (0..n).collect();
        let natural_fill = count_fill(&csc, &natural_perm);

        let (nd_perm, _) = nd_ordering(&csc);
        assert!(is_valid_permutation(&nd_perm, n));
        let nd_fill = count_fill(&csc, &nd_perm);

        assert!(
            nd_fill <= natural_fill,
            "ND fill {} should be <= natural fill {}",
            nd_fill, natural_fill
        );
    }

    #[test]
    fn test_nd_2d_grid_fill_vs_amd() {
        let csc = build_2d_grid(20);
        let n = 400;

        let (amd_perm, _) = super::super::amd::amd_ordering(&csc);
        let amd_fill = count_fill(&csc, &amd_perm);

        let (nd_perm, _) = nd_ordering(&csc);
        assert!(is_valid_permutation(&nd_perm, n));
        let nd_fill = count_fill(&csc, &nd_perm);

        // ND should be comparable to or better than AMD on structured grids
        // Allow up to 50% more fill as a loose bound (ND without METIS may not always beat AMD)
        assert!(
            nd_fill <= amd_fill * 3 / 2,
            "ND fill {} should be comparable to AMD fill {} (within 50%)",
            nd_fill, amd_fill
        );
    }

    #[test]
    fn test_nd_arrow_matrix() {
        let n = 20;
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 10.0));
            if i < n - 1 {
                triplets.push((i, n - 1, 1.0));
            }
        }
        let csc = csc_from_upper_triplets(n, &triplets);

        let natural_perm: Vec<usize> = (0..n).collect();
        let natural_fill = count_fill(&csc, &natural_perm);

        let (nd_perm, _) = nd_ordering(&csc);
        assert!(is_valid_permutation(&nd_perm, n));
        let nd_fill = count_fill(&csc, &nd_perm);

        // ND should not be catastrophically worse than natural on arrow
        assert!(
            nd_fill <= natural_fill * 2,
            "ND fill {} should be reasonable vs natural fill {} for arrow",
            nd_fill, natural_fill
        );
    }

    #[test]
    fn test_nd_solver_integration() {
        use crate::solver::{Solver, SolverOptions};
        use crate::ordering::Ordering;

        // 2D Laplacian 5x5
        let nx = 5;
        let n = nx * nx;
        let coo = {
            let mut rows = Vec::new();
            let mut cols = Vec::new();
            let mut vals = Vec::new();
            for iy in 0..nx {
                for ix in 0..nx {
                    let idx = iy * nx + ix;
                    rows.push(idx); cols.push(idx); vals.push(4.0);
                    if ix + 1 < nx {
                        let right = iy * nx + ix + 1;
                        rows.push(idx); cols.push(right); vals.push(-1.0);
                    }
                    if iy + 1 < nx {
                        let below = (iy + 1) * nx + ix;
                        rows.push(idx); cols.push(below); vals.push(-1.0);
                    }
                }
            }
            CooMatrix::new(n, rows, cols, vals).unwrap()
        };

        let mut opts = SolverOptions::default();
        opts.ordering = Ordering::NestedDissection;
        let mut solver = Solver::new(opts);
        let inertia = solver.analyze_and_factor(&coo).unwrap();
        assert_eq!(inertia.positive, n);

        // Solve Ax = b where b = A * ones
        let csc = CscMatrix::from_coo(&coo);
        let x_true: Vec<f64> = vec![1.0; n];
        let mut b = vec![0.0; n];
        csc.matvec(&x_true, &mut b);

        let mut x = vec![0.0; n];
        solver.solve(&b, &mut x).unwrap();

        // Check residual
        let mut r = vec![0.0; n];
        csc.matvec(&x, &mut r);
        let mut max_res = 0.0f64;
        for i in 0..n {
            max_res = max_res.max((r[i] - b[i]).abs());
        }
        assert!(
            max_res < 1e-10,
            "Residual {} should be < 1e-10",
            max_res
        );
    }

    #[test]
    fn test_nd_disconnected_graph() {
        use crate::solver::{Solver, SolverOptions};
        use crate::ordering::Ordering;

        // Two disconnected 3x3 blocks
        let n = 6;
        let coo = CooMatrix::new(n,
            vec![0, 0, 1, 1, 2, 3, 3, 4, 4, 5],
            vec![0, 1, 1, 2, 2, 3, 4, 4, 5, 5],
            vec![3.0, -1.0, 3.0, -1.0, 3.0, 3.0, -1.0, 3.0, -1.0, 3.0],
        ).unwrap();

        let csc = CscMatrix::from_coo(&coo);
        let (perm, perm_inv) = nd_ordering(&csc);
        assert!(is_valid_permutation(&perm, n));
        assert!(is_valid_permutation(&perm_inv, n));

        // Verify solve correctness
        let mut opts = SolverOptions::default();
        opts.ordering = Ordering::NestedDissection;
        let mut solver = Solver::new(opts);
        solver.analyze_and_factor(&coo).unwrap();

        let x_true = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut b = vec![0.0; n];
        csc.matvec(&x_true, &mut b);

        let mut x = vec![0.0; n];
        solver.solve(&b, &mut x).unwrap();

        let mut max_err = 0.0f64;
        for i in 0..n {
            max_err = max_err.max((x[i] - x_true[i]).abs());
        }
        assert!(max_err < 1e-10, "Solution error {} should be < 1e-10", max_err);
    }
}
