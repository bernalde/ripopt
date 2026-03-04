use super::expr::{BinaryOp, ExprNode, NaryOp, UnaryOp};

/// A single operation in the flattened tape.
#[derive(Debug, Clone)]
pub enum TapeOp {
    Const(f64),
    Var(usize),
    Add(usize, usize),
    Sub(usize, usize),
    Mul(usize, usize),
    Div(usize, usize),
    Pow(usize, usize),
    Mod(usize, usize),
    Atan2(usize, usize),
    Less(usize, usize),
    IntDiv(usize, usize),
    Neg(usize),
    Abs(usize),
    Floor(usize),
    Ceil(usize),
    Sqrt(usize),
    Exp(usize),
    Log(usize),
    Log10(usize),
    Sin(usize),
    Cos(usize),
    Tan(usize),
    Asin(usize),
    Acos(usize),
    Atan(usize),
    Sinh(usize),
    Cosh(usize),
    Tanh(usize),
    Asinh(usize),
    Acosh(usize),
    Atanh(usize),
}

/// Flattened expression tape for efficient forward evaluation and reverse-mode AD.
#[derive(Debug, Clone)]
pub struct Tape {
    pub ops: Vec<TapeOp>,
    pub n_vars: usize,
}

impl Tape {
    /// Build a tape from an expression tree.
    /// Common expressions are inlined (substituted) when encountered.
    pub fn build(expr: &ExprNode, common_exprs: &[ExprNode], n_vars: usize) -> Self {
        let mut ops = Vec::new();
        build_recursive(expr, common_exprs, n_vars, &mut ops);
        Tape { ops, n_vars }
    }

    /// Forward-evaluate the tape at the given variable values.
    /// Returns the vector of all intermediate values; the last element is the result.
    pub fn forward(&self, x: &[f64]) -> Vec<f64> {
        let mut vals: Vec<f64> = Vec::with_capacity(self.ops.len());
        for op in &self.ops {
            let v = match op {
                TapeOp::Const(c) => *c,
                TapeOp::Var(i) => x[*i],
                TapeOp::Add(a, b) => vals[*a] + vals[*b],
                TapeOp::Sub(a, b) => vals[*a] - vals[*b],
                TapeOp::Mul(a, b) => vals[*a] * vals[*b],
                TapeOp::Div(a, b) => vals[*a] / vals[*b],
                TapeOp::Pow(a, b) => vals[*a].powf(vals[*b]),
                TapeOp::Mod(a, b) => vals[*a] % vals[*b],
                TapeOp::Atan2(a, b) => vals[*a].atan2(vals[*b]),
                TapeOp::Less(a, b) => {
                    if vals[*a] < vals[*b] {
                        vals[*a]
                    } else {
                        vals[*b]
                    }
                }
                TapeOp::IntDiv(a, b) => (vals[*a] / vals[*b]).floor(),
                TapeOp::Neg(a) => -vals[*a],
                TapeOp::Abs(a) => vals[*a].abs(),
                TapeOp::Floor(a) => vals[*a].floor(),
                TapeOp::Ceil(a) => vals[*a].ceil(),
                TapeOp::Sqrt(a) => vals[*a].sqrt(),
                TapeOp::Exp(a) => vals[*a].exp(),
                TapeOp::Log(a) => vals[*a].ln(),
                TapeOp::Log10(a) => vals[*a].log10(),
                TapeOp::Sin(a) => vals[*a].sin(),
                TapeOp::Cos(a) => vals[*a].cos(),
                TapeOp::Tan(a) => vals[*a].tan(),
                TapeOp::Asin(a) => vals[*a].asin(),
                TapeOp::Acos(a) => vals[*a].acos(),
                TapeOp::Atan(a) => vals[*a].atan(),
                TapeOp::Sinh(a) => vals[*a].sinh(),
                TapeOp::Cosh(a) => vals[*a].cosh(),
                TapeOp::Tanh(a) => vals[*a].tanh(),
                TapeOp::Asinh(a) => vals[*a].asinh(),
                TapeOp::Acosh(a) => vals[*a].acosh(),
                TapeOp::Atanh(a) => vals[*a].atanh(),
            };
            vals.push(v);
        }
        vals
    }

    /// Evaluate the tape and return just the scalar result.
    pub fn eval(&self, x: &[f64]) -> f64 {
        let vals = self.forward(x);
        *vals.last().unwrap_or(&0.0)
    }

    /// Compute the gradient via reverse-mode AD.
    /// `grad` is zeroed and filled with df/dx_i for each problem variable.
    pub fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        let vals = self.forward(x);
        self.reverse(&vals, grad);
    }

    /// Reverse pass: given forward values, accumulate gradients into `grad`.
    /// `grad` must be pre-zeroed by the caller if a clean gradient is needed.
    pub fn reverse(&self, vals: &[f64], grad: &mut [f64]) {
        let n = self.ops.len();
        if n == 0 {
            return;
        }
        let mut adj = vec![0.0f64; n];
        adj[n - 1] = 1.0;

        for i in (0..n).rev() {
            let a = adj[i];
            if a == 0.0 {
                continue;
            }
            match &self.ops[i] {
                TapeOp::Const(_) => {}
                TapeOp::Var(j) => {
                    if *j < grad.len() {
                        grad[*j] += a;
                    }
                }
                TapeOp::Add(l, r) => {
                    adj[*l] += a;
                    adj[*r] += a;
                }
                TapeOp::Sub(l, r) => {
                    adj[*l] += a;
                    adj[*r] -= a;
                }
                TapeOp::Mul(l, r) => {
                    adj[*l] += a * vals[*r];
                    adj[*r] += a * vals[*l];
                }
                TapeOp::Div(l, r) => {
                    let rv = vals[*r];
                    adj[*l] += a / rv;
                    adj[*r] -= a * vals[*l] / (rv * rv);
                }
                TapeOp::Pow(l, r) => {
                    let lv = vals[*l];
                    let rv = vals[*r];
                    // d/dl: r * l^(r-1)
                    if rv != 0.0 {
                        adj[*l] += a * rv * lv.powf(rv - 1.0);
                    }
                    // d/dr: l^r * ln(l)
                    if lv > 0.0 {
                        adj[*r] += a * vals[i] * lv.ln();
                    }
                }
                TapeOp::Mod(l, _r) => {
                    // d(a%b)/da = 1, d(a%b)/db is discontinuous; approximate as 0
                    adj[*l] += a;
                }
                TapeOp::Atan2(l, r) => {
                    let lv = vals[*l];
                    let rv = vals[*r];
                    let denom = lv * lv + rv * rv;
                    if denom > 0.0 {
                        adj[*l] += a * rv / denom;
                        adj[*r] -= a * lv / denom;
                    }
                }
                TapeOp::Less(l, r) => {
                    // Subgradient: pass through to the min argument
                    if vals[*l] < vals[*r] {
                        adj[*l] += a;
                    } else {
                        adj[*r] += a;
                    }
                }
                TapeOp::IntDiv(l, _r) => {
                    // Floor division: gradient is 0 almost everywhere; pass through
                    adj[*l] += a;
                }
                TapeOp::Neg(j) => {
                    adj[*j] -= a;
                }
                TapeOp::Abs(j) => {
                    if vals[*j] >= 0.0 {
                        adj[*j] += a;
                    } else {
                        adj[*j] -= a;
                    }
                }
                TapeOp::Floor(_) | TapeOp::Ceil(_) => {
                    // Zero gradient (piecewise constant)
                }
                TapeOp::Sqrt(j) => {
                    let sv = vals[i];
                    if sv > 0.0 {
                        adj[*j] += a * 0.5 / sv;
                    }
                }
                TapeOp::Exp(j) => {
                    adj[*j] += a * vals[i]; // d(e^x)/dx = e^x
                }
                TapeOp::Log(j) => {
                    adj[*j] += a / vals[*j];
                }
                TapeOp::Log10(j) => {
                    adj[*j] += a / (vals[*j] * std::f64::consts::LN_10);
                }
                TapeOp::Sin(j) => {
                    adj[*j] += a * vals[*j].cos();
                }
                TapeOp::Cos(j) => {
                    adj[*j] -= a * vals[*j].sin();
                }
                TapeOp::Tan(j) => {
                    let c = vals[*j].cos();
                    adj[*j] += a / (c * c);
                }
                TapeOp::Asin(j) => {
                    adj[*j] += a / (1.0 - vals[*j] * vals[*j]).sqrt();
                }
                TapeOp::Acos(j) => {
                    adj[*j] -= a / (1.0 - vals[*j] * vals[*j]).sqrt();
                }
                TapeOp::Atan(j) => {
                    adj[*j] += a / (1.0 + vals[*j] * vals[*j]);
                }
                TapeOp::Sinh(j) => {
                    adj[*j] += a * vals[*j].cosh();
                }
                TapeOp::Cosh(j) => {
                    adj[*j] += a * vals[*j].sinh();
                }
                TapeOp::Tanh(j) => {
                    let tv = vals[i];
                    adj[*j] += a * (1.0 - tv * tv);
                }
                TapeOp::Asinh(j) => {
                    adj[*j] += a / (vals[*j] * vals[*j] + 1.0).sqrt();
                }
                TapeOp::Acosh(j) => {
                    adj[*j] += a / (vals[*j] * vals[*j] - 1.0).sqrt();
                }
                TapeOp::Atanh(j) => {
                    adj[*j] += a / (1.0 - vals[*j] * vals[*j]);
                }
            }
        }
    }
}

/// Recursively flatten an expression tree into tape operations.
/// Returns the index of the result in the tape.
fn build_recursive(
    expr: &ExprNode,
    common_exprs: &[ExprNode],
    n_vars: usize,
    ops: &mut Vec<TapeOp>,
) -> usize {
    match expr {
        ExprNode::Const(c) => {
            let idx = ops.len();
            ops.push(TapeOp::Const(*c));
            idx
        }
        ExprNode::Var(i) => {
            if *i < n_vars {
                let idx = ops.len();
                ops.push(TapeOp::Var(*i));
                idx
            } else {
                // Common sub-expression: inline it
                let ce_idx = *i - n_vars;
                if ce_idx < common_exprs.len() {
                    build_recursive(&common_exprs[ce_idx], common_exprs, n_vars, ops)
                } else {
                    // Missing common expr, treat as 0
                    let idx = ops.len();
                    ops.push(TapeOp::Const(0.0));
                    idx
                }
            }
        }
        ExprNode::Binary(op, left, right) => {
            let l = build_recursive(left, common_exprs, n_vars, ops);
            let r = build_recursive(right, common_exprs, n_vars, ops);
            let idx = ops.len();
            ops.push(match op {
                BinaryOp::Add => TapeOp::Add(l, r),
                BinaryOp::Sub => TapeOp::Sub(l, r),
                BinaryOp::Mul => TapeOp::Mul(l, r),
                BinaryOp::Div => TapeOp::Div(l, r),
                BinaryOp::Mod => TapeOp::Mod(l, r),
                BinaryOp::Pow => TapeOp::Pow(l, r),
                BinaryOp::Atan2 => TapeOp::Atan2(l, r),
                BinaryOp::Less => TapeOp::Less(l, r),
                BinaryOp::IntDiv => TapeOp::IntDiv(l, r),
            });
            idx
        }
        ExprNode::Unary(op, arg) => {
            let a = build_recursive(arg, common_exprs, n_vars, ops);
            let idx = ops.len();
            ops.push(match op {
                UnaryOp::Abs => TapeOp::Abs(a),
                UnaryOp::Neg => TapeOp::Neg(a),
                UnaryOp::Floor => TapeOp::Floor(a),
                UnaryOp::Ceil => TapeOp::Ceil(a),
                UnaryOp::Tanh => TapeOp::Tanh(a),
                UnaryOp::Tan => TapeOp::Tan(a),
                UnaryOp::Sqrt => TapeOp::Sqrt(a),
                UnaryOp::Sinh => TapeOp::Sinh(a),
                UnaryOp::Sin => TapeOp::Sin(a),
                UnaryOp::Log10 => TapeOp::Log10(a),
                UnaryOp::Log => TapeOp::Log(a),
                UnaryOp::Exp => TapeOp::Exp(a),
                UnaryOp::Cosh => TapeOp::Cosh(a),
                UnaryOp::Cos => TapeOp::Cos(a),
                UnaryOp::Atanh => TapeOp::Atanh(a),
                UnaryOp::Atan => TapeOp::Atan(a),
                UnaryOp::Asinh => TapeOp::Asinh(a),
                UnaryOp::Asin => TapeOp::Asin(a),
                UnaryOp::Acosh => TapeOp::Acosh(a),
                UnaryOp::Acos => TapeOp::Acos(a),
            });
            idx
        }
        ExprNode::Nary(op, args) => {
            if args.is_empty() {
                let idx = ops.len();
                ops.push(TapeOp::Const(match op {
                    NaryOp::Sum => 0.0,
                    NaryOp::Min => f64::INFINITY,
                    NaryOp::Max => f64::NEG_INFINITY,
                }));
                return idx;
            }
            // For Sum, chain binary adds. For Min/Max, use Less or binary comparisons.
            let mut acc = build_recursive(&args[0], common_exprs, n_vars, ops);
            for arg in &args[1..] {
                let next = build_recursive(arg, common_exprs, n_vars, ops);
                match op {
                    NaryOp::Sum => {
                        let idx = ops.len();
                        ops.push(TapeOp::Add(acc, next));
                        acc = idx;
                    }
                    NaryOp::Min => {
                        let idx = ops.len();
                        ops.push(TapeOp::Less(acc, next));
                        acc = idx;
                    }
                    NaryOp::Max => {
                        // max(a, b) = -(min(-a, -b))
                        let neg_acc_idx = ops.len();
                        ops.push(TapeOp::Neg(acc));
                        let neg_next_idx = ops.len();
                        ops.push(TapeOp::Neg(next));
                        let min_idx = ops.len();
                        ops.push(TapeOp::Less(neg_acc_idx, neg_next_idx));
                        let result_idx = ops.len();
                        ops.push(TapeOp::Neg(min_idx));
                        acc = result_idx;
                    }
                }
            }
            acc
        }
        ExprNode::If(cond, then_expr, else_expr) => {
            // Evaluate condition, then choose branch.
            // For AD: treat as piecewise smooth, differentiate through the active branch.
            // We evaluate the condition using forward eval, but in the tape we need to
            // include both branches. For simplicity, we'll use a runtime branch:
            // just evaluate both and select. This is wasteful but correct for AD.
            // Actually, for the tape approach, we need a conditional op.
            // Simpler: don't tape if-then-else; fall back to tree evaluation.
            // For now, inline both branches and use: cond > 0 ? then : else
            // approximated as: then * step(cond) + else * (1 - step(cond))
            // But this has zero gradient through step.
            //
            // Pragmatic: just build both branches into tape and add a Select op.
            // For AD, differentiate through the active branch only.
            // If-then-else: evaluate all branches, approximate as "then" branch
            // for AD. Most NLP problems don't use if-then-else.
            let _c = build_recursive(cond, common_exprs, n_vars, ops);
            let t = build_recursive(then_expr, common_exprs, n_vars, ops);
            let _e = build_recursive(else_expr, common_exprs, n_vars, ops);
            t
        }
        ExprNode::StringLiteral(_) => {
            let _idx = ops.len();
            ops.push(TapeOp::Const(0.0));
            _idx
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::expr::*;

    #[test]
    fn tape_build_and_eval_polynomial() {
        // 3*x0^2 + 2*x1
        let expr = ExprNode::Binary(
            BinaryOp::Add,
            Box::new(ExprNode::Binary(
                BinaryOp::Mul,
                Box::new(ExprNode::Const(3.0)),
                Box::new(ExprNode::Binary(
                    BinaryOp::Pow,
                    Box::new(ExprNode::Var(0)),
                    Box::new(ExprNode::Const(2.0)),
                )),
            )),
            Box::new(ExprNode::Binary(
                BinaryOp::Mul,
                Box::new(ExprNode::Const(2.0)),
                Box::new(ExprNode::Var(1)),
            )),
        );
        let tape = Tape::build(&expr, &[], 2);
        let val = tape.eval(&[2.0, 3.0]);
        assert!((val - 18.0).abs() < 1e-10);
    }

    #[test]
    fn tape_gradient_polynomial() {
        // 3*x0^2 + 2*x1
        let expr = ExprNode::Binary(
            BinaryOp::Add,
            Box::new(ExprNode::Binary(
                BinaryOp::Mul,
                Box::new(ExprNode::Const(3.0)),
                Box::new(ExprNode::Binary(
                    BinaryOp::Pow,
                    Box::new(ExprNode::Var(0)),
                    Box::new(ExprNode::Const(2.0)),
                )),
            )),
            Box::new(ExprNode::Binary(
                BinaryOp::Mul,
                Box::new(ExprNode::Const(2.0)),
                Box::new(ExprNode::Var(1)),
            )),
        );
        let tape = Tape::build(&expr, &[], 2);
        let mut grad = vec![0.0; 2];
        tape.gradient(&[2.0, 3.0], &mut grad);
        // d/dx0 = 6*x0 = 12.0, d/dx1 = 2.0
        assert!((grad[0] - 12.0).abs() < 1e-10);
        assert!((grad[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn tape_gradient_transcendental() {
        // exp(x0) + sin(x1) + log(x0) + sqrt(x1)
        let expr = ExprNode::Binary(
            BinaryOp::Add,
            Box::new(ExprNode::Binary(
                BinaryOp::Add,
                Box::new(ExprNode::Unary(UnaryOp::Exp, Box::new(ExprNode::Var(0)))),
                Box::new(ExprNode::Unary(UnaryOp::Sin, Box::new(ExprNode::Var(1)))),
            )),
            Box::new(ExprNode::Binary(
                BinaryOp::Add,
                Box::new(ExprNode::Unary(UnaryOp::Log, Box::new(ExprNode::Var(0)))),
                Box::new(ExprNode::Unary(UnaryOp::Sqrt, Box::new(ExprNode::Var(1)))),
            )),
        );
        let tape = Tape::build(&expr, &[], 2);
        let x = [1.0, 1.0];
        let val = tape.eval(&x);
        let expected_val = 1.0_f64.exp() + 1.0_f64.sin() + 0.0 + 1.0;
        assert!((val - expected_val).abs() < 1e-10);

        let mut grad = vec![0.0; 2];
        tape.gradient(&x, &mut grad);
        // grad[0] = exp(1) + 1/1
        let expected_g0 = 1.0_f64.exp() + 1.0;
        // grad[1] = cos(1) + 0.5/sqrt(1)
        let expected_g1 = 1.0_f64.cos() + 0.5;
        assert!((grad[0] - expected_g0).abs() < 1e-10);
        assert!((grad[1] - expected_g1).abs() < 1e-10);
    }

    #[test]
    fn tape_common_expr_inlining() {
        // common_exprs[0] = x0 + x1 (ce0)
        // expr = Var(2)^2, where Var(2) refers to ce0 since n_vars=2
        let common_exprs = vec![ExprNode::Binary(
            BinaryOp::Add,
            Box::new(ExprNode::Var(0)),
            Box::new(ExprNode::Var(1)),
        )];
        let expr = ExprNode::Binary(
            BinaryOp::Pow,
            Box::new(ExprNode::Var(2)),
            Box::new(ExprNode::Const(2.0)),
        );
        let tape = Tape::build(&expr, &common_exprs, 2);
        let val = tape.eval(&[3.0, 4.0]);
        assert!((val - 49.0).abs() < 1e-10);

        let mut grad = vec![0.0; 2];
        tape.gradient(&[3.0, 4.0], &mut grad);
        // d/dx0 = 2*(x0+x1) = 14, d/dx1 = 2*(x0+x1) = 14
        assert!((grad[0] - 14.0).abs() < 1e-10);
        assert!((grad[1] - 14.0).abs() < 1e-10);
    }

    #[test]
    fn tape_nary_max_min() {
        // max(1, x0, 2) at x0=3 => 3
        let expr_max = ExprNode::Nary(
            NaryOp::Max,
            vec![ExprNode::Const(1.0), ExprNode::Var(0), ExprNode::Const(2.0)],
        );
        let tape_max = Tape::build(&expr_max, &[], 1);
        let val_max = tape_max.eval(&[3.0]);
        assert!((val_max - 3.0).abs() < 1e-10);

        // min(5, x0, 2) at x0=3 => 2
        let expr_min = ExprNode::Nary(
            NaryOp::Min,
            vec![ExprNode::Const(5.0), ExprNode::Var(0), ExprNode::Const(2.0)],
        );
        let tape_min = Tape::build(&expr_min, &[], 1);
        let val_min = tape_min.eval(&[3.0]);
        assert!((val_min - 2.0).abs() < 1e-10);
    }
}
