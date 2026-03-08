use crate::csc::CscMatrix;

/// Elimination tree for a sparse symmetric matrix.
#[derive(Debug, Clone)]
pub struct EliminationTree {
    /// `parent[i]` = parent of node i, or None if i is a root.
    pub parent: Vec<Option<usize>>,
    /// Number of nodes.
    pub n: usize,
}

impl EliminationTree {
    /// Build the elimination tree from an upper-triangle CSC matrix.
    /// Uses the row-subtree algorithm (Liu 1990) with path compression.
    pub fn from_csc(csc: &CscMatrix) -> Self {
        let n = csc.n;
        let mut parent: Vec<Option<usize>> = vec![None; n];
        let mut ancestor = vec![0usize; n]; // for path compression

        for i in 0..n {
            ancestor[i] = i;
        }

        for j in 0..n {
            // Process column j: for each row i < j in the upper triangle
            for idx in csc.col_ptr[j]..csc.col_ptr[j + 1] {
                let i = csc.row_idx[idx];
                if i >= j {
                    continue; // only process upper triangle entries where i < j
                }
                // Walk from i up to the root, setting parent = j for the root
                let mut node = i;
                loop {
                    let anc = ancestor[node];
                    if anc == j {
                        break;
                    }
                    if anc == node {
                        // node is a root — make j its parent
                        parent[node] = Some(j);
                        ancestor[node] = j;
                        break;
                    }
                    // Path compression
                    ancestor[node] = j;
                    node = anc;
                }
            }
        }

        EliminationTree { parent, n }
    }

    /// Compute the children list for each node.
    pub fn children(&self) -> Vec<Vec<usize>> {
        let mut ch = vec![Vec::new(); self.n];
        for i in 0..self.n {
            if let Some(p) = self.parent[i] {
                ch[p].push(i);
            }
        }
        ch
    }

    /// Compute a postordering of the elimination tree.
    /// Returns a permutation where leaves come before their parents.
    pub fn postorder(&self) -> Vec<usize> {
        let children = self.children();
        let mut order = Vec::with_capacity(self.n);
        let mut stack = Vec::new();

        // Find roots (nodes with no parent)
        for i in 0..self.n {
            if self.parent[i].is_none() {
                stack.push((i, false));
            }
        }

        // Iterative DFS postorder
        while let Some((node, visited)) = stack.pop() {
            if visited {
                order.push(node);
            } else {
                stack.push((node, true));
                // Push children in reverse so they're processed left-to-right
                for &child in children[node].iter().rev() {
                    stack.push((child, false));
                }
            }
        }

        order
    }
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

    #[test]
    fn test_etree_tridiagonal() {
        // Tridiagonal 4x4: nonzeros at (i,i) and (i,i+1)
        // Upper triangle: (0,0),(0,1),(1,1),(1,2),(2,2),(2,3),(3,3)
        // Elimination tree: 0->1->2->3 (chain)
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let etree = EliminationTree::from_csc(&csc);
        assert_eq!(etree.parent[0], Some(1));
        assert_eq!(etree.parent[1], Some(2));
        assert_eq!(etree.parent[2], Some(3));
        assert_eq!(etree.parent[3], None); // root
    }

    #[test]
    fn test_etree_arrow() {
        // Arrow matrix 4x4: all entries connect to last column
        // (0,0),(0,3),(1,1),(1,3),(2,2),(2,3),(3,3)
        // Etree: 0->3, 1->3, 2->3, 3=root
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 3, 1.0),
            (1, 1, 1.0), (1, 3, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let etree = EliminationTree::from_csc(&csc);
        assert_eq!(etree.parent[0], Some(3));
        assert_eq!(etree.parent[1], Some(3));
        assert_eq!(etree.parent[2], Some(3));
        assert_eq!(etree.parent[3], None);
    }

    #[test]
    fn test_etree_diagonal() {
        // Diagonal 3x3: no off-diagonal entries
        let csc = csc_from_upper_triplets(3, &[
            (0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0),
        ]);
        let etree = EliminationTree::from_csc(&csc);
        // Each node is a root (forest)
        assert_eq!(etree.parent[0], None);
        assert_eq!(etree.parent[1], None);
        assert_eq!(etree.parent[2], None);
    }

    #[test]
    fn test_postorder_chain() {
        // Chain: 0->1->2->3
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 1, 1.0),
            (1, 1, 1.0), (1, 2, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let etree = EliminationTree::from_csc(&csc);
        let order = etree.postorder();
        assert_eq!(order.len(), 4);
        // In postorder, each node must appear before its parent
        let mut pos = vec![0; 4];
        for (i, &node) in order.iter().enumerate() {
            pos[node] = i;
        }
        for i in 0..4 {
            if let Some(p) = etree.parent[i] {
                assert!(pos[i] < pos[p], "node {} should come before parent {}", i, p);
            }
        }
    }

    #[test]
    fn test_postorder_arrow() {
        // Arrow: 0,1,2 -> 3
        let csc = csc_from_upper_triplets(4, &[
            (0, 0, 1.0), (0, 3, 1.0),
            (1, 1, 1.0), (1, 3, 1.0),
            (2, 2, 1.0), (2, 3, 1.0),
            (3, 3, 1.0),
        ]);
        let etree = EliminationTree::from_csc(&csc);
        let order = etree.postorder();
        assert_eq!(order.len(), 4);
        // 3 must be last (root)
        assert_eq!(*order.last().unwrap(), 3);
        // 0, 1, 2 should come before 3
        let mut pos = vec![0; 4];
        for (i, &node) in order.iter().enumerate() {
            pos[node] = i;
        }
        for i in 0..3 {
            assert!(pos[i] < pos[3]);
        }
    }
}
