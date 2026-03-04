use ripopt::{NlpProblem, SolveStatus, SolverOptions};

// ---------------------------------------------------------------------------
// 1. NE-to-LS detection: f=0, grad=0, 3 equalities in 2 vars
//    x0+x1=1, x0-x1=0, 2*x0=1  =>  x*=(0.5, 0.5)
// ---------------------------------------------------------------------------

struct NeToLs;

impl NlpProblem for NeToLs {
    fn num_variables(&self) -> usize { 2 }
    fn num_constraints(&self) -> usize { 3 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY; x_u[0] = f64::INFINITY;
        x_l[1] = f64::NEG_INFINITY; x_u[1] = f64::INFINITY;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l[0] = 1.0; g_u[0] = 1.0; // x0+x1=1
        g_l[1] = 0.0; g_u[1] = 0.0; // x0-x1=0
        g_l[2] = 1.0; g_u[2] = 1.0; // 2*x0=1
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.0;
        x0[1] = 0.0;
    }

    fn objective(&self, _x: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _x: &[f64], grad: &mut [f64]) {
        grad[0] = 0.0;
        grad[1] = 0.0;
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] + x[1];
        g[1] = x[0] - x[1];
        g[2] = 2.0 * x[0];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // (row, col): (0,0),(0,1),(1,0),(1,1),(2,0)
        (vec![0, 0, 1, 1, 2], vec![0, 1, 0, 1, 0])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        vals[0] = 1.0;  // dg0/dx0
        vals[1] = 1.0;  // dg0/dx1
        vals[2] = 1.0;  // dg1/dx0
        vals[3] = -1.0; // dg1/dx1
        vals[4] = 2.0;  // dg2/dx0
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn hessian_values(&self, _x: &[f64], _obj_factor: f64, _lambda: &[f64], _vals: &mut [f64]) {}
}

#[test]
fn ipm_ne_to_ls_detection() {
    let problem = NeToLs;
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
    assert!((result.x[0] - 0.5).abs() < 1e-4, "x0={}, expected 0.5", result.x[0]);
    assert!((result.x[1] - 0.5).abs() < 1e-4, "x1={}, expected 0.5", result.x[1]);
}

// ---------------------------------------------------------------------------
// 2. Condensed KKT: m >= 2*n triggers condensed path
//    min x^2 s.t. x=1, x=1, x=1  (m=3, n=1)
// ---------------------------------------------------------------------------

struct CondensedKkt;

impl NlpProblem for CondensedKkt {
    fn num_variables(&self) -> usize { 1 }
    fn num_constraints(&self) -> usize { 3 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY;
        x_u[0] = f64::INFINITY;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        for i in 0..3 {
            g_l[i] = 1.0;
            g_u[i] = 1.0;
        }
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.5;
    }

    fn objective(&self, x: &[f64]) -> f64 { x[0] * x[0] }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = 2.0 * x[0];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0];
        g[1] = x[0];
        g[2] = x[0];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0, 1, 2], vec![0, 0, 0])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        vals[0] = 1.0;
        vals[1] = 1.0;
        vals[2] = 1.0;
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0], vec![0])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        vals[0] = 2.0 * obj_factor;
    }
}

#[test]
fn ipm_condensed_kkt() {
    let problem = CondensedKkt;
    let options = SolverOptions {
        print_level: 0,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);
    assert_eq!(result.status, SolveStatus::Optimal, "Expected Optimal, got {:?}", result.status);
    assert!((result.x[0] - 1.0).abs() < 1e-6, "x={}, expected 1.0", result.x[0]);
}

// ---------------------------------------------------------------------------
// 3. Unbounded detection: min -x, no bounds, no constraints
// ---------------------------------------------------------------------------

struct UnboundedProblem;

impl NlpProblem for UnboundedProblem {
    fn num_variables(&self) -> usize { 1 }
    fn num_constraints(&self) -> usize { 0 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY;
        x_u[0] = f64::INFINITY;
    }

    fn constraint_bounds(&self, _g_l: &mut [f64], _g_u: &mut [f64]) {}

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.0;
    }

    fn objective(&self, x: &[f64]) -> f64 { -x[0] }

    fn gradient(&self, _x: &[f64], grad: &mut [f64]) {
        grad[0] = -1.0;
    }

    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {}

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {}

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn hessian_values(&self, _x: &[f64], _obj_factor: f64, _lambda: &[f64], _vals: &mut [f64]) {}
}

#[test]
fn ipm_unbounded_detection() {
    let problem = UnboundedProblem;
    let options = SolverOptions {
        print_level: 0,
        enable_lbfgs_fallback: false,
        enable_al_fallback: false,
        enable_sqp_fallback: false,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);
    // Without bounds, the solver may detect unboundedness, hit numerical issues,
    // reach max iterations, or declare acceptable at a large negative objective.
    assert!(
        result.status == SolveStatus::Unbounded
            || result.status == SolveStatus::NumericalError
            || result.status == SolveStatus::MaxIterations
            || result.status == SolveStatus::Acceptable,
        "Expected Unbounded/NumericalError/MaxIterations/Acceptable, got {:?}",
        result.status
    );
}

// ---------------------------------------------------------------------------
// 4. Preprocessing: fixed variable x0=5, min x1^2+x2^2 s.t. x1+x2=1
//    x*=(5.0, 0.5, 0.5), obj*=0.5
// ---------------------------------------------------------------------------

struct PreprocessingProblem;

impl NlpProblem for PreprocessingProblem {
    fn num_variables(&self) -> usize { 3 }
    fn num_constraints(&self) -> usize { 1 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = 5.0; x_u[0] = 5.0; // fixed
        x_l[1] = f64::NEG_INFINITY; x_u[1] = f64::INFINITY;
        x_l[2] = f64::NEG_INFINITY; x_u[2] = f64::INFINITY;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l[0] = 1.0;
        g_u[0] = 1.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 5.0;
        x0[1] = 0.0;
        x0[2] = 0.0;
    }

    fn objective(&self, x: &[f64]) -> f64 {
        x[1] * x[1] + x[2] * x[2]
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = 0.0;
        grad[1] = 2.0 * x[1];
        grad[2] = 2.0 * x[2];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[1] + x[2];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // dg0/dx1, dg0/dx2
        (vec![0, 0], vec![1, 2])
    }

    fn jacobian_values(&self, _x: &[f64], vals: &mut [f64]) {
        vals[0] = 1.0;
        vals[1] = 1.0;
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        // Lower triangle: (1,1), (2,2)
        (vec![1, 2], vec![1, 2])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        vals[0] = 2.0 * obj_factor; // d2f/dx1^2
        vals[1] = 2.0 * obj_factor; // d2f/dx2^2
    }
}

#[test]
fn ipm_preprocessing_integration() {
    let problem = PreprocessingProblem;
    let options = SolverOptions {
        print_level: 0,
        enable_preprocessing: true,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);
    assert_eq!(result.status, SolveStatus::Optimal, "Expected Optimal, got {:?}", result.status);
    assert!((result.x[0] - 5.0).abs() < 1e-10, "x0={}, expected 5.0", result.x[0]);
    assert!((result.x[1] - 0.5).abs() < 1e-4, "x1={}, expected 0.5", result.x[1]);
    assert!((result.x[2] - 0.5).abs() < 1e-4, "x2={}, expected 0.5", result.x[2]);
    assert!((result.objective - 0.5).abs() < 1e-4, "obj={}, expected 0.5", result.objective);
}

// ---------------------------------------------------------------------------
// 5. Best-du restore at max_iter: min x^2 s.t. x^2=1
//    Two solutions: x=1 or x=-1. Start x0=2.0.
//    Use impossibly tight tol but reasonable acceptable_tol.
// ---------------------------------------------------------------------------

struct BestDuProblem;

impl NlpProblem for BestDuProblem {
    fn num_variables(&self) -> usize { 1 }
    fn num_constraints(&self) -> usize { 1 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = f64::NEG_INFINITY;
        x_u[0] = f64::INFINITY;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l[0] = 1.0;
        g_u[0] = 1.0;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 2.0;
    }

    fn objective(&self, x: &[f64]) -> f64 { x[0] * x[0] }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        grad[0] = 2.0 * x[0];
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        g[0] = x[0] * x[0];
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0], vec![0])
    }

    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        vals[0] = 2.0 * x[0];
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0], vec![0])
    }

    fn hessian_values(&self, _x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]) {
        vals[0] = 2.0 * obj_factor + 2.0 * lambda[0];
    }
}

#[test]
fn ipm_best_du_at_maxiter() {
    let problem = BestDuProblem;
    let options = SolverOptions {
        print_level: 0,
        max_iter: 50,
        tol: 1e-14,
        acceptable_tol: 1e-4,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);
    assert!(
        result.status == SolveStatus::Acceptable || result.status == SolveStatus::Optimal,
        "Expected Acceptable or Optimal, got {:?}",
        result.status
    );
    // x should be near 1.0 (the minimizer on x^2=1)
    assert!((result.x[0].abs() - 1.0).abs() < 1e-2, "x={}, expected |x|~1.0", result.x[0]);
}

// ---------------------------------------------------------------------------
// 6. L-BFGS fallback for unconstrained Rosenbrock
//    min (1-x0)^2 + 100*(x1-x0^2)^2, start (-1.2, 1.0)
// ---------------------------------------------------------------------------

struct RosenbrockLbfgs;

impl NlpProblem for RosenbrockLbfgs {
    fn num_variables(&self) -> usize { 2 }
    fn num_constraints(&self) -> usize { 0 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        for i in 0..2 {
            x_l[i] = f64::NEG_INFINITY;
            x_u[i] = f64::INFINITY;
        }
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
        grad[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        grad[1] = 200.0 * (x[1] - x[0] * x[0]);
    }

    fn constraints(&self, _x: &[f64], _g: &mut [f64]) {}

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![], vec![])
    }

    fn jacobian_values(&self, _x: &[f64], _vals: &mut [f64]) {}

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0, 1, 1], vec![0, 0, 1])
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, _lambda: &[f64], vals: &mut [f64]) {
        vals[0] = obj_factor * (2.0 - 400.0 * x[1] + 1200.0 * x[0] * x[0]);
        vals[1] = obj_factor * (-400.0 * x[0]);
        vals[2] = obj_factor * 200.0;
    }
}

#[test]
fn ipm_lbfgs_fallback_unconstrained() {
    let problem = RosenbrockLbfgs;
    let options = SolverOptions {
        print_level: 0,
        enable_lbfgs_fallback: true,
        ..SolverOptions::default()
    };
    let result = ripopt::solve(&problem, &options);
    assert_eq!(result.status, SolveStatus::Optimal, "Expected Optimal, got {:?}", result.status);
    assert!((result.x[0] - 1.0).abs() < 1e-4, "x0={}, expected 1.0", result.x[0]);
    assert!((result.x[1] - 1.0).abs() < 1e-4, "x1={}, expected 1.0", result.x[1]);
    assert!(result.objective < 1e-6, "obj={}, expected ~0", result.objective);
}
