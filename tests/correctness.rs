use ripopt::{NlpProblem, SolveStatus, SolverOptions};

// ---------------------------------------------------------------------------
// 1. Rosenbrock (unconstrained)
//    min f(x) = (1 - x1)^2 + 100*(x2 - x1^2)^2
//    x* = (1, 1), f* = 0
// ---------------------------------------------------------------------------

struct Rosenbrock;

impl NlpProblem for Rosenbrock {
    fn num_variables(&self) -> usize {
        2
    }

    fn num_constraints(&self) -> usize {
        0
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY;
        x_l[1] = f64::NEG_INFINITY;
        x_u[0] = f64::INFINITY;
        x_u[1] = f64::INFINITY;
    }

    fn constraint_bounds(&self, _g_l: &mut [f64], _g_u: &mut [f64]) {}

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = -1.2;
        x0[1] = 1.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        let a = 1.0 - x[0];
        let b = x[1] - x[0] * x[0];
        a * a + 100.0 * b * b
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        let x1 = x[0];
        let x2 = x[1];
        grad[0] = -2.0 * (1.0 - x1) - 400.0 * x1 * (x2 - x1 * x1);
        grad[1] = 200.0 * (x2 - x1 * x1);
    }

    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {}

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {}

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Lower triangle: (0,0), (1,0), (1,1)
        (vec![0, 1, 1], vec![0, 0, 1])
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        let x1 = x[0];
        let x2 = x[1];
        // H[0,0] = 2 - 400*x2 + 1200*x1^2
        vals[0] = obj_factor * (2.0 - 400.0 * x2 + 1200.0 * x1 * x1);
        // H[1,0] = -400*x1
        vals[1] = obj_factor * (-400.0 * x1);
        // H[1,1] = 200
        vals[2] = obj_factor * 200.0;
    }
}

#[test]
fn rosenbrock_unconstrained() {
    let problem = Rosenbrock;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );
    assert!(
        (result.x[0] - 1.0).abs() < 1e-4,
        "x1 should be ~1.0, got {}",
        result.x[0]
    );
    assert!(
        (result.x[1] - 1.0).abs() < 1e-4,
        "x2 should be ~1.0, got {}",
        result.x[1]
    );
    assert!(
        result.objective.abs() < 1e-3,
        "f* should be ~0.0, got {}",
        result.objective
    );
}

// ---------------------------------------------------------------------------
// 2. Simple constrained QP
//    min f(x) = 0.5*(x1^2 + x2^2)
//    s.t. x1 + x2 = 1
//    x* = (0.5, 0.5), f* = 0.25
// ---------------------------------------------------------------------------

struct SimpleQP;

impl NlpProblem for SimpleQP {
    fn num_variables(&self) -> usize {
        2
    }

    fn num_constraints(&self) -> usize {
        1
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY;
        x_l[1] = f64::NEG_INFINITY;
        x_u[0] = f64::INFINITY;
        x_u[1] = f64::INFINITY;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l[0] = 1.0;
        g_u[0] = 1.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.0;
        x0[1] = 0.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1])
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = x[0];
        grad[1] = x[1];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] + x[1];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // J = [1, 1] -> row 0 col 0, row 0 col 1
        (vec![0, 0], vec![0, 1])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        vals[0] = 1.0;
        vals[1] = 1.0;
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Diagonal: (0,0) and (1,1) — both are on the lower triangle
        (vec![0, 1], vec![0, 1])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        // Hessian of f: diag(1, 1). Constraint is linear so its Hessian is zero.
        vals[0] = obj_factor * 1.0;
        vals[1] = obj_factor * 1.0;
    }
}

#[test]
fn simple_constrained_qp() {
    let problem = SimpleQP;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );
    assert!(
        (result.x[0] - 0.5).abs() < 1e-4,
        "x1 should be ~0.5, got {}",
        result.x[0]
    );
    assert!(
        (result.x[1] - 0.5).abs() < 1e-4,
        "x2 should be ~0.5, got {}",
        result.x[1]
    );
    assert!(
        (result.objective - 0.25).abs() < 1e-3,
        "f* should be ~0.25, got {}",
        result.objective
    );
}

// ---------------------------------------------------------------------------
// 3. HS071
//    min f = x1*x4*(x1+x2+x3) + x3
//    s.t. g1 = x1*x2*x3*x4 >= 25
//         g2 = x1^2 + x2^2 + x3^2 + x4^2 = 40
//    1 <= xi <= 5, i = 1..4
//    x0 = (1, 5, 5, 1)
//    f* ≈ 17.014
// ---------------------------------------------------------------------------

struct HS071;

impl NlpProblem for HS071 {
    fn num_variables(&self) -> usize {
        4
    }

    fn num_constraints(&self) -> usize {
        2
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..4 {
            x_l[i] = 1.0;
            x_u[i] = 5.0;
        }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        // g1 >= 25 (no upper bound)
        g_l[0] = 25.0;
        g_u[0] = f64::INFINITY;
        // g2 = 40 (equality)
        g_l[1] = 40.0;
        g_u[1] = 40.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 1.0;
        x0[1] = 5.0;
        x0[2] = 5.0;
        x0[3] = 1.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        // df/dx1 = x4*(x1+x2+x3) + x1*x4 = x4*(2*x1 + x2 + x3)
        grad[0] = x[3] * (2.0 * x[0] + x[1] + x[2]);
        // df/dx2 = x1*x4
        grad[1] = x[0] * x[3];
        // df/dx3 = x1*x4 + 1
        grad[2] = x[0] * x[3] + 1.0;
        // df/dx4 = x1*(x1+x2+x3)
        grad[3] = x[0] * (x[0] + x[1] + x[2]);
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] * x[1] * x[2] * x[3];
        g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // g1 depends on all 4 vars, g2 depends on all 4 vars
        // Row 0 (g1): cols 0,1,2,3
        // Row 1 (g2): cols 0,1,2,3
        (
            vec![0, 0, 0, 0, 1, 1, 1, 1],
            vec![0, 1, 2, 3, 0, 1, 2, 3],
        )
    }

    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        // dg1/dx1 = x2*x3*x4
        vals[0] = x[1] * x[2] * x[3];
        // dg1/dx2 = x1*x3*x4
        vals[1] = x[0] * x[2] * x[3];
        // dg1/dx3 = x1*x2*x4
        vals[2] = x[0] * x[1] * x[3];
        // dg1/dx4 = x1*x2*x3
        vals[3] = x[0] * x[1] * x[2];
        // dg2/dx1 = 2*x1
        vals[4] = 2.0 * x[0];
        // dg2/dx2 = 2*x2
        vals[5] = 2.0 * x[1];
        // dg2/dx3 = 2*x3
        vals[6] = 2.0 * x[2];
        // dg2/dx4 = 2*x4
        vals[7] = 2.0 * x[3];
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Lower triangle of the 4x4 Hessian of the Lagrangian.
        // We include all lower-triangle entries that can be non-zero.
        //
        // Objective Hessian non-zeros (lower triangle):
        //   (0,0): 2*x4
        //   (1,0): x4,  (2,0): x4,  (3,0): 2*x1+x2+x3
        //   (2,1): 0,   (3,1): x1
        //   (3,2): x1
        //
        // g1 = x1*x2*x3*x4, Hessian non-zeros (lower triangle):
        //   (1,0): x3*x4, (2,0): x2*x4, (3,0): x2*x3
        //   (2,1): x1*x4, (3,1): x1*x3
        //   (3,2): x1*x2
        //
        // g2 = sum xi^2, Hessian (lower triangle):
        //   (0,0): 2, (1,1): 2, (2,2): 2, (3,3): 2
        //
        // Combined lower-triangle entries (row >= col):
        // (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3)
        (
            vec![0, 1, 1, 2, 2, 2, 3, 3, 3, 3],
            vec![0, 0, 1, 0, 1, 2, 0, 1, 2, 3],
        )
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]) {
        // Indices match hessian_structure order:
        // 0: (0,0), 1: (1,0), 2: (1,1), 3: (2,0), 4: (2,1), 5: (2,2),
        // 6: (3,0), 7: (3,1), 8: (3,2), 9: (3,3)

        // ---- Objective Hessian ----
        // d2f/dx1dx1 = 2*x4
        vals[0] = obj_factor * 2.0 * x[3];
        // d2f/dx2dx1 = x4
        vals[1] = obj_factor * x[3];
        // d2f/dx2dx2 = 0
        vals[2] = 0.0;
        // d2f/dx3dx1 = x4
        vals[3] = obj_factor * x[3];
        // d2f/dx3dx2 = 0
        vals[4] = 0.0;
        // d2f/dx3dx3 = 0
        vals[5] = 0.0;
        // d2f/dx4dx1 = 2*x1 + x2 + x3
        vals[6] = obj_factor * (2.0 * x[0] + x[1] + x[2]);
        // d2f/dx4dx2 = x1
        vals[7] = obj_factor * x[0];
        // d2f/dx4dx3 = x1
        vals[8] = obj_factor * x[0];
        // d2f/dx4dx4 = 0
        vals[9] = 0.0;

        // ---- Constraint 1 Hessian: g1 = x1*x2*x3*x4 ----
        // (0,0): 0
        // (1,0): x3*x4
        vals[1] += lambda[0] * x[2] * x[3];
        // (1,1): 0
        // (2,0): x2*x4
        vals[3] += lambda[0] * x[1] * x[3];
        // (2,1): x1*x4
        vals[4] += lambda[0] * x[0] * x[3];
        // (2,2): 0
        // (3,0): x2*x3
        vals[6] += lambda[0] * x[1] * x[2];
        // (3,1): x1*x3
        vals[7] += lambda[0] * x[0] * x[2];
        // (3,2): x1*x2
        vals[8] += lambda[0] * x[0] * x[1];
        // (3,3): 0

        // ---- Constraint 2 Hessian: g2 = sum xi^2 ----
        // Only diagonal entries: 2 each
        vals[0] += lambda[1] * 2.0;
        vals[2] += lambda[1] * 2.0;
        vals[5] += lambda[1] * 2.0;
        vals[9] += lambda[1] * 2.0;
    }
}

#[test]
fn hs071_constrained() {
    let problem = HS071;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );
    assert!(
        (result.objective - 17.014).abs() < 0.1,
        "f* should be ~17.014, got {}",
        result.objective
    );

    // Check that constraints are satisfied
    // g1 = x1*x2*x3*x4 >= 25
    let g1 = result.x[0] * result.x[1] * result.x[2] * result.x[3];
    assert!(
        g1 >= 25.0 - 1e-3,
        "g1 = {} should be >= 25",
        g1
    );
    // g2 = x1^2 + x2^2 + x3^2 + x4^2 = 40
    let g2: f64 = result.x.iter().map(|xi| xi * xi).sum();
    assert!(
        (g2 - 40.0).abs() < 1e-3,
        "g2 = {} should be ~40",
        g2
    );
    // Bounds: 1 <= xi <= 5
    for (i, &xi) in result.x.iter().enumerate() {
        assert!(
            xi >= 1.0 - 1e-4 && xi <= 5.0 + 1e-4,
            "x[{}] = {} out of bounds [1, 5]",
            i,
            xi
        );
    }
}

// ---------------------------------------------------------------------------
// 4. Bound-constrained quadratic (HS035-like)
//    min f(x) = 9 - 8x1 - 6x2 - 4x3 + 2x1^2 + 2x2^2 + x3^2 + 2x1*x2 + 2x1*x3
//    s.t. x1 + x2 + 2*x3 <= 3
//         x1 >= 0, x2 >= 0, x3 >= 0
//    x* = (4/3, 7/9, 4/9), f* = 1/9 ≈ 0.1111
//    Hessian is constant: H = [[4, 2, 2], [2, 4, 0], [2, 0, 2]]
// ---------------------------------------------------------------------------

struct BoundConstrainedQuadratic;

impl NlpProblem for BoundConstrainedQuadratic {
    fn num_variables(&self) -> usize {
        3
    }

    fn num_constraints(&self) -> usize {
        1
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..3 {
            x_l[i] = 0.0;
            x_u[i] = f64::INFINITY;
        }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        // x1 + x2 + 2*x3 <= 3
        g_l[0] = f64::NEG_INFINITY;
        g_u[0] = 3.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.5;
        x0[1] = 0.5;
        x0[2] = 0.5;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        9.0 - 8.0 * x[0] - 6.0 * x[1] - 4.0 * x[2]
            + 2.0 * x[0] * x[0]
            + 2.0 * x[1] * x[1]
            + x[2] * x[2]
            + 2.0 * x[0] * x[1]
            + 2.0 * x[0] * x[2]
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        // df/dx1 = -8 + 4*x1 + 2*x2 + 2*x3
        grad[0] = -8.0 + 4.0 * x[0] + 2.0 * x[1] + 2.0 * x[2];
        // df/dx2 = -6 + 2*x1 + 4*x2
        grad[1] = -6.0 + 2.0 * x[0] + 4.0 * x[1];
        // df/dx3 = -4 + 2*x1 + 2*x3
        grad[2] = -4.0 + 2.0 * x[0] + 2.0 * x[2];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] + x[1] + 2.0 * x[2];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // J = [1, 1, 2] -> row 0 cols 0,1,2
        (vec![0, 0, 0], vec![0, 1, 2])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        vals[0] = 1.0;
        vals[1] = 1.0;
        vals[2] = 2.0;
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Lower triangle of 3x3 constant Hessian:
        // (0,0), (1,0), (1,1), (2,0), (2,2)
        // Note: H[2,1] = 0 so we skip it
        (vec![0, 1, 1, 2, 2], vec![0, 0, 1, 0, 2])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        // Constant Hessian of objective: H = [[4, 2, 2], [2, 4, 0], [2, 0, 2]]
        // Constraint is linear, so its Hessian is zero.
        // (0,0): 4
        vals[0] = obj_factor * 4.0;
        // (1,0): 2
        vals[1] = obj_factor * 2.0;
        // (1,1): 4
        vals[2] = obj_factor * 4.0;
        // (2,0): 2
        vals[3] = obj_factor * 2.0;
        // (2,2): 2
        vals[4] = obj_factor * 2.0;
    }
}

#[test]
fn bound_constrained_quadratic() {
    let problem = BoundConstrainedQuadratic;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );

    let expected_x = [4.0 / 3.0, 7.0 / 9.0, 4.0 / 9.0];
    let expected_f = 1.0 / 9.0;

    assert!(
        (result.objective - expected_f).abs() < 1e-4,
        "f* should be ~{}, got {}",
        expected_f,
        result.objective
    );
    for i in 0..3 {
        assert!(
            (result.x[i] - expected_x[i]).abs() < 1e-3,
            "x[{}] should be ~{}, got {}",
            i,
            expected_x[i],
            result.x[i]
        );
    }
}

// ---------------------------------------------------------------------------
// 5. Pure bound-constrained (no general constraints)
//    min f(x) = (x1-1)^2 + (x2-2)^2 + (x3-3)^2 + (x4-4)^2
//    s.t. 0 <= x_i <= 3 for all i
//    x* = (1, 2, 3, 3), f* = 1.0
// ---------------------------------------------------------------------------

struct PureBoundConstrained;

impl NlpProblem for PureBoundConstrained {
    fn num_variables(&self) -> usize {
        4
    }

    fn num_constraints(&self) -> usize {
        0
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..4 {
            x_l[i] = 0.0;
            x_u[i] = 3.0;
        }
    }

    fn constraint_bounds(&self, _g_l: &mut [f64], _g_u: &mut [f64]) {}

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.0;
        x0[1] = 0.0;
        x0[2] = 0.0;
        x0[3] = 0.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        (x[0] - 1.0).powi(2)
            + (x[1] - 2.0).powi(2)
            + (x[2] - 3.0).powi(2)
            + (x[3] - 4.0).powi(2)
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = 2.0 * (x[0] - 1.0);
        grad[1] = 2.0 * (x[1] - 2.0);
        grad[2] = 2.0 * (x[2] - 3.0);
        grad[3] = 2.0 * (x[3] - 4.0);
    }

    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {}

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {}

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Diagonal: (0,0), (1,1), (2,2), (3,3)
        (vec![0, 1, 2, 3], vec![0, 1, 2, 3])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        // Hessian is diagonal with all entries = 2
        for v in vals.iter_mut() {
            *v = obj_factor * 2.0;
        }
    }
}

#[test]
fn pure_bound_constrained() {
    let problem = PureBoundConstrained;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );

    let expected_x = [1.0, 2.0, 3.0, 3.0];
    let expected_f = 1.0;

    assert!(
        (result.objective - expected_f).abs() < 1e-4,
        "f* should be ~{}, got {}",
        expected_f,
        result.objective
    );
    for i in 0..4 {
        assert!(
            (result.x[i] - expected_x[i]).abs() < 1e-3,
            "x[{}] should be ~{}, got {}",
            i,
            expected_x[i],
            result.x[i]
        );
    }
}

// ---------------------------------------------------------------------------
// 6. Multiple equality constraints
//    min f(x) = x1^2 + x2^2 + x3^2
//    s.t. x1 + x2 + x3 = 1
//         x1 - x2 = 0
//    x* = (1/3, 1/3, 1/3), f* = 1/3
// ---------------------------------------------------------------------------

struct MultipleEqualityConstraints;

impl NlpProblem for MultipleEqualityConstraints {
    fn num_variables(&self) -> usize {
        3
    }

    fn num_constraints(&self) -> usize {
        2
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..3 {
            x_l[i] = f64::NEG_INFINITY;
            x_u[i] = f64::INFINITY;
        }
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        // g1: x1 + x2 + x3 = 1
        g_l[0] = 1.0;
        g_u[0] = 1.0;
        // g2: x1 - x2 = 0
        g_l[1] = 0.0;
        g_u[1] = 0.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.0;
        x0[1] = 0.0;
        x0[2] = 0.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = 2.0 * x[0];
        grad[1] = 2.0 * x[1];
        grad[2] = 2.0 * x[2];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] + x[1] + x[2];
        g[1] = x[0] - x[1];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // g1 depends on x1, x2, x3: row 0 cols 0,1,2
        // g2 depends on x1, x2:     row 1 cols 0,1
        (vec![0, 0, 0, 1, 1], vec![0, 1, 2, 0, 1])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        // dg1/dx1 = 1
        vals[0] = 1.0;
        // dg1/dx2 = 1
        vals[1] = 1.0;
        // dg1/dx3 = 1
        vals[2] = 1.0;
        // dg2/dx1 = 1
        vals[3] = 1.0;
        // dg2/dx2 = -1
        vals[4] = -1.0;
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Hessian of objective is diagonal: diag(2, 2, 2)
        // Constraints are linear so their Hessians are zero.
        (vec![0, 1, 2], vec![0, 1, 2])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        // Hessian of objective: diag(2, 2, 2)
        vals[0] = obj_factor * 2.0;
        vals[1] = obj_factor * 2.0;
        vals[2] = obj_factor * 2.0;
    }
}

#[test]
fn multiple_equality_constraints() {
    let problem = MultipleEqualityConstraints;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);

    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "Expected Optimal or Acceptable, got {:?}",
        result.status
    );

    let expected_x = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    let expected_f = 1.0 / 3.0;

    assert!(
        (result.objective - expected_f).abs() < 1e-4,
        "f* should be ~{}, got {}",
        expected_f,
        result.objective
    );
    for i in 0..3 {
        assert!(
            (result.x[i] - expected_x[i]).abs() < 1e-3,
            "x[{}] should be ~{}, got {}",
            i,
            expected_x[i],
            result.x[i]
        );
    }
}
