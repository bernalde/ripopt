/// Parsed NL file header containing problem dimensions.
#[derive(Debug, Clone)]
pub struct NlHeader {
    pub n_vars: usize,
    pub n_constrs: usize,
    pub n_objs: usize,
    pub n_ranges: usize,
    pub n_eqns: usize,
    pub n_nl_constrs: usize,
    pub n_nl_objs: usize,
    /// NL vars in both objectives and constraints.
    pub nlvb: usize,
    /// NL vars in constraints only.
    pub nlvc: usize,
    /// NL vars in objectives only.
    pub nlvo: usize,
    /// Common sub-expression counts: [both, constr1, obj1, constr, obj].
    pub common_expr_counts: [usize; 5],
    pub nnz_jac: usize,
    pub nnz_grad: usize,
    /// Number of AMPL imported (external) functions declared in the file.
    pub n_funcs: usize,
}

impl NlHeader {
    /// Total number of common sub-expressions (defined variables).
    pub fn n_common_exprs(&self) -> usize {
        self.common_expr_counts.iter().sum()
    }
}

/// Parse NL header from line iterator.
/// Consumes the magic line and 9 dimension lines.
pub fn parse_header<'a>(lines: &mut impl Iterator<Item = &'a str>) -> Result<NlHeader, String> {
    // Line 0: magic line "g3 ..."
    let magic = lines.next().ok_or("Empty NL file")?;
    if !magic.starts_with('g') {
        return Err(format!("Expected text NL file (starting with 'g'), got: {}", magic));
    }

    // Parse 9 dimension lines
    let mut dim_lines: Vec<Vec<usize>> = Vec::new();
    for i in 0..9 {
        let line = lines
            .next()
            .ok_or(format!("NL header truncated at line {}", i + 1))?;
        let nums: Vec<usize> = parse_nums(line);
        dim_lines.push(nums);
    }

    let get = |line: usize, field: usize| -> usize {
        dim_lines
            .get(line)
            .and_then(|l| l.get(field))
            .copied()
            .unwrap_or(0)
    };

    Ok(NlHeader {
        n_vars: get(0, 0),
        n_constrs: get(0, 1),
        n_objs: get(0, 2),
        n_ranges: get(0, 3),
        n_eqns: get(0, 4),
        n_nl_constrs: get(1, 0),
        n_nl_objs: get(1, 1),
        nlvb: get(3, 0),
        nlvc: get(3, 1),
        nlvo: get(3, 2),
        common_expr_counts: [get(8, 0), get(8, 1), get(8, 2), get(8, 3), get(8, 4)],
        nnz_jac: get(6, 0),
        nnz_grad: get(6, 1),
        // Line 4 (0-indexed dim line): `nwv nfunc arith flags`.
        n_funcs: get(4, 1),
    })
}

/// Parse unsigned integers from a line, ignoring comments after '#'.
fn parse_nums(line: &str) -> Vec<usize> {
    let line = line.split('#').next().unwrap_or("");
    line.split_whitespace()
        .filter_map(|s| s.parse::<usize>().ok())
        .collect()
}
