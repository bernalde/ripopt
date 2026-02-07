/// Trait defining a nonlinear programming problem.
///
/// All methods use a buffer-filling pattern to avoid allocations in the hot loop.
/// The caller provides pre-allocated slices and the implementation fills them.
pub trait NlpProblem {
    /// Number of primal variables.
    fn num_variables(&self) -> usize;

    /// Number of constraints.
    fn num_constraints(&self) -> usize;

    /// Fill variable bounds: x_l[i] <= x[i] <= x_u[i].
    /// Use f64::NEG_INFINITY / f64::INFINITY for unbounded.
    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]);

    /// Fill constraint bounds: g_l[i] <= g(x)[i] <= g_u[i].
    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]);

    /// Fill initial primal point.
    fn initial_point(&self, x0: &mut [f64]);

    /// Evaluate objective f(x).
    fn objective(&self, x: &[f64]) -> f64;

    /// Fill gradient of objective: grad[i] = df/dx_i.
    fn gradient(&self, x: &[f64], grad: &mut [f64]);

    /// Evaluate constraints: g[i] = g_i(x).
    fn constraints(&self, x: &[f64], g: &mut [f64]);

    /// Return the sparsity structure of the constraint Jacobian.
    /// Returns (row_indices, col_indices) in triplet format.
    /// Only the non-zero entries need to be specified.
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>);

    /// Fill Jacobian values at x in the same order as jacobian_structure().
    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]);

    /// Return the sparsity structure of the Lagrangian Hessian (lower triangle only).
    /// Returns (row_indices, col_indices) in triplet format.
    /// This is the Hessian of: obj_factor * f(x) + sum_i lambda[i] * g_i(x).
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>);

    /// Fill Hessian values at x with the given obj_factor and constraint multipliers lambda.
    /// Only lower triangle entries in the same order as hessian_structure().
    fn hessian_values(&self, x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]);
}
