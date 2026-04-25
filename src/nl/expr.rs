/// Binary operation opcodes.
#[derive(Debug, Clone, Copy)]
pub enum BinaryOp {
    Add,    // o0
    Sub,    // o1
    Mul,    // o2
    Div,    // o3
    Mod,    // o4
    Pow,    // o5
    Atan2,  // o48
    Less,   // o22
    IntDiv, // o55
}

/// Unary operation opcodes.
#[derive(Debug, Clone, Copy)]
pub enum UnaryOp {
    Abs,    // o15
    Neg,    // o16
    Floor,  // o13
    Ceil,   // o14
    Tanh,   // o37
    Tan,    // o38
    Sqrt,   // o39
    Sinh,   // o40
    Sin,    // o41
    Log10,  // o42
    Log,    // o43
    Exp,    // o44
    Cosh,   // o45
    Cos,    // o46
    Atanh,  // o47
    Atan,   // o49
    Asinh,  // o50
    Asin,   // o51
    Acosh,  // o52
    Acos,   // o53
}

/// N-ary operation.
#[derive(Debug, Clone, Copy)]
pub enum NaryOp {
    Sum,    // o54 (SUMLIST)
    Min,    // o11
    Max,    // o12
}

/// Expression tree node.
#[derive(Debug, Clone)]
pub enum ExprNode {
    Const(f64),
    /// Variable reference (index into problem vars or defined vars).
    Var(usize),
    Binary(BinaryOp, Box<ExprNode>, Box<ExprNode>),
    Unary(UnaryOp, Box<ExprNode>),
    Nary(NaryOp, Vec<ExprNode>),
    /// If-then-else: condition > 0 → then_expr, else → else_expr.
    If(Box<ExprNode>, Box<ExprNode>, Box<ExprNode>),
    /// String literal (ignored in evaluation).
    StringLiteral(String),
    /// AMPL imported (external) function call. Cannot be evaluated natively —
    /// callers must surface an error before any evaluation occurs.
    Funcall { id: usize, args: Vec<ExprNode> },
}

/// Parse a prefix-notation expression tree from NL file lines.
pub fn parse_expr<'a>(lines: &mut impl Iterator<Item = &'a str>) -> Result<ExprNode, String> {
    let line = lines
        .next()
        .ok_or("Unexpected end of expression")?
        .trim();

    // Strip inline comments
    let token = line.split('#').next().unwrap_or("").trim();

    if token.starts_with('n') {
        // Numeric constant
        let val: f64 = token[1..]
            .parse()
            .map_err(|e| format!("Bad number '{}': {}", &token[1..], e))?;
        Ok(ExprNode::Const(val))
    } else if token.starts_with('v') {
        // Variable reference
        let idx: usize = token[1..]
            .parse()
            .map_err(|e| format!("Bad var index '{}': {}", &token[1..], e))?;
        Ok(ExprNode::Var(idx))
    } else if token.starts_with('h') {
        // String literal: h<len>:<string>
        let rest = &token[1..];
        if let Some(colon) = rest.find(':') {
            Ok(ExprNode::StringLiteral(rest[colon + 1..].to_string()))
        } else {
            Ok(ExprNode::StringLiteral(String::new()))
        }
    } else if token.starts_with('o') {
        let opcode: usize = token[1..]
            .parse()
            .map_err(|e| format!("Bad opcode '{}': {}", &token[1..], e))?;
        parse_op(opcode, lines)
    } else if token.starts_with('f') {
        // AMPL imported (external) function call: `f<id> <nargs>` on one line,
        // followed by nargs child expressions (each possibly multi-line).
        let rest = &token[1..];
        let mut parts = rest.split_whitespace();
        let id_str = parts.next().ok_or_else(|| format!("Missing function id in '{}'", token))?;
        let nargs_str = parts.next().ok_or_else(|| format!("Missing nargs in '{}'", token))?;
        let id: usize = id_str
            .parse()
            .map_err(|e| format!("Bad function id '{}': {}", id_str, e))?;
        let nargs: usize = nargs_str
            .parse()
            .map_err(|e| format!("Bad funcall nargs '{}': {}", nargs_str, e))?;
        let args: Result<Vec<_>, _> = (0..nargs).map(|_| parse_expr(lines)).collect();
        Ok(ExprNode::Funcall { id, args: args? })
    } else {
        Err(format!("Unknown expression token: '{}'", token))
    }
}

fn parse_op<'a>(
    opcode: usize,
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<ExprNode, String> {
    match opcode {
        // Binary ops
        0 => binary(BinaryOp::Add, lines),
        1 => binary(BinaryOp::Sub, lines),
        2 => binary(BinaryOp::Mul, lines),
        3 => binary(BinaryOp::Div, lines),
        4 => binary(BinaryOp::Mod, lines),
        5 => binary(BinaryOp::Pow, lines),
        22 => binary(BinaryOp::Less, lines),
        48 => binary(BinaryOp::Atan2, lines),
        55 => binary(BinaryOp::IntDiv, lines),

        // Unary ops
        13 => unary(UnaryOp::Floor, lines),
        14 => unary(UnaryOp::Ceil, lines),
        15 => unary(UnaryOp::Abs, lines),
        16 => unary(UnaryOp::Neg, lines),
        37 => unary(UnaryOp::Tanh, lines),
        38 => unary(UnaryOp::Tan, lines),
        39 => unary(UnaryOp::Sqrt, lines),
        40 => unary(UnaryOp::Sinh, lines),
        41 => unary(UnaryOp::Sin, lines),
        42 => unary(UnaryOp::Log10, lines),
        43 => unary(UnaryOp::Log, lines),
        44 => unary(UnaryOp::Exp, lines),
        45 => unary(UnaryOp::Cosh, lines),
        46 => unary(UnaryOp::Cos, lines),
        47 => unary(UnaryOp::Atanh, lines),
        49 => unary(UnaryOp::Atan, lines),
        50 => unary(UnaryOp::Asinh, lines),
        51 => unary(UnaryOp::Asin, lines),
        52 => unary(UnaryOp::Acosh, lines),
        53 => unary(UnaryOp::Acos, lines),

        // N-ary ops
        11 => nary(NaryOp::Min, lines),
        12 => nary(NaryOp::Max, lines),
        54 => nary(NaryOp::Sum, lines),

        // If-then-else
        35 => {
            let cond = parse_expr(lines)?;
            let then_expr = parse_expr(lines)?;
            let else_expr = parse_expr(lines)?;
            Ok(ExprNode::If(
                Box::new(cond),
                Box::new(then_expr),
                Box::new(else_expr),
            ))
        }

        // If-then (else = 0)
        65 => {
            let cond = parse_expr(lines)?;
            let then_expr = parse_expr(lines)?;
            Ok(ExprNode::If(
                Box::new(cond),
                Box::new(then_expr),
                Box::new(ExprNode::Const(0.0)),
            ))
        }

        _ => Err(format!("Unsupported opcode: o{}", opcode)),
    }
}

fn binary<'a>(
    op: BinaryOp,
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<ExprNode, String> {
    let left = parse_expr(lines)?;
    let right = parse_expr(lines)?;
    Ok(ExprNode::Binary(op, Box::new(left), Box::new(right)))
}

fn unary<'a>(op: UnaryOp, lines: &mut impl Iterator<Item = &'a str>) -> Result<ExprNode, String> {
    let arg = parse_expr(lines)?;
    Ok(ExprNode::Unary(op, Box::new(arg)))
}

fn nary<'a>(op: NaryOp, lines: &mut impl Iterator<Item = &'a str>) -> Result<ExprNode, String> {
    let count_line = lines
        .next()
        .ok_or("Expected count for n-ary op")?
        .trim();
    let count_str = count_line.split('#').next().unwrap_or("").trim();
    let count: usize = count_str
        .parse()
        .map_err(|e| format!("Bad n-ary count '{}': {}", count_str, e))?;
    let args: Result<Vec<_>, _> = (0..count).map(|_| parse_expr(lines)).collect();
    Ok(ExprNode::Nary(op, args?))
}

/// Evaluate an expression tree given variable values and common expression values.
/// `vals[0..n_vars]` = problem variables, `vals[n_vars..]` = common expression values.
pub fn eval_expr(node: &ExprNode, vals: &[f64]) -> f64 {
    match node {
        ExprNode::Const(c) => *c,
        ExprNode::Var(i) => vals[*i],
        ExprNode::Binary(op, l, r) => {
            let lv = eval_expr(l, vals);
            let rv = eval_expr(r, vals);
            match op {
                BinaryOp::Add => lv + rv,
                BinaryOp::Sub => lv - rv,
                BinaryOp::Mul => lv * rv,
                BinaryOp::Div => lv / rv,
                BinaryOp::Mod => lv % rv,
                BinaryOp::Pow => lv.powf(rv),
                BinaryOp::Atan2 => lv.atan2(rv),
                BinaryOp::Less => {
                    if lv < rv {
                        lv
                    } else {
                        rv
                    }
                }
                BinaryOp::IntDiv => (lv / rv).floor(),
            }
        }
        ExprNode::Unary(op, a) => {
            let av = eval_expr(a, vals);
            match op {
                UnaryOp::Abs => av.abs(),
                UnaryOp::Neg => -av,
                UnaryOp::Floor => av.floor(),
                UnaryOp::Ceil => av.ceil(),
                UnaryOp::Tanh => av.tanh(),
                UnaryOp::Tan => av.tan(),
                UnaryOp::Sqrt => av.sqrt(),
                UnaryOp::Sinh => av.sinh(),
                UnaryOp::Sin => av.sin(),
                UnaryOp::Log10 => av.log10(),
                UnaryOp::Log => av.ln(),
                UnaryOp::Exp => av.exp(),
                UnaryOp::Cosh => av.cosh(),
                UnaryOp::Cos => av.cos(),
                UnaryOp::Atanh => av.atanh(),
                UnaryOp::Atan => av.atan(),
                UnaryOp::Asinh => av.asinh(),
                UnaryOp::Asin => av.asin(),
                UnaryOp::Acosh => av.acosh(),
                UnaryOp::Acos => av.acos(),
            }
        }
        ExprNode::Nary(op, args) => {
            let vals_iter = args.iter().map(|a| eval_expr(a, vals));
            match op {
                NaryOp::Sum => vals_iter.sum(),
                NaryOp::Min => vals_iter.fold(f64::INFINITY, f64::min),
                NaryOp::Max => vals_iter.fold(f64::NEG_INFINITY, f64::max),
            }
        }
        ExprNode::If(cond, then_expr, else_expr) => {
            if eval_expr(cond, vals) > 0.0 {
                eval_expr(then_expr, vals)
            } else {
                eval_expr(else_expr, vals)
            }
        }
        ExprNode::StringLiteral(_) => 0.0,
        ExprNode::Funcall { .. } => unreachable!(
            "ExprNode::Funcall should have been rejected by NlProblem::from_nl_data"
        ),
    }
}
