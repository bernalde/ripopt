#!/usr/bin/env python3
"""
Solve a comprehensive set of NLP test problems using cyipopt and output JSON
results for comparison with other solvers (e.g. a Rust Ipopt wrapper).

Problems:
  1. Rosenbrock (2 vars, unconstrained, bounds [-5, 5])
  2. HS071 (4 vars, 2 constraints: 1 inequality, 1 equality)
  3. Simple QP (2 vars, 1 equality constraint)
  4. HS035 (3 vars, 1 inequality constraint, bounds x >= 0)
  5. Bound-constrained quadratic (4 vars, no general constraints, bounds [0, 10])

Each problem is solved via the cyipopt.Problem (C-level) interface so we can
extract full multiplier information (constraint multipliers y, bound
multipliers z_L and z_U).
"""

import json
import time
import numpy as np
import cyipopt


# ---------------------------------------------------------------------------
# Problem definition helpers
# ---------------------------------------------------------------------------

def _make_problem_obj(objective, gradient, constraints, jacobian,
                      jacobianstructure, hessian, hessianstructure,
                      intermediate_cb=None):
    """
    Build a lightweight problem-object that cyipopt.Problem accepts.
    The optional *intermediate_cb* is called on every Ipopt iteration.
    """
    class ProblemDef:
        pass

    p = ProblemDef()
    p.objective = objective
    p.gradient = gradient
    p.constraints = constraints
    p.jacobian = jacobian
    p.jacobianstructure = jacobianstructure
    p.hessian = hessian
    p.hessianstructure = hessianstructure

    if intermediate_cb is not None:
        p.intermediate = intermediate_cb
    else:
        p.intermediate = lambda *a, **kw: True

    return p


# ---------------------------------------------------------------------------
# Problem 1: Rosenbrock  (2 vars, 0 constraints, bounds [-5, 5])
# min (1-x1)^2 + 100*(x2 - x1^2)^2
# ---------------------------------------------------------------------------

_ROSENBROCK_N = 2
_ROSENBROCK_M = 0
_ROSENBROCK_X0 = np.array([-1.0, 1.0])


def _rb_objective(x):
    return (1.0 - x[0])**2 + 100.0 * (x[1] - x[0]**2)**2


def _rb_gradient(x):
    return np.array([
        -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0]**2),
        200.0 * (x[1] - x[0]**2),
    ])


def _rb_constraints(x):
    return np.array([])


def _rb_jacobian(x):
    return np.array([])


def _rb_jacobianstructure():
    return (np.array([], dtype=int), np.array([], dtype=int))


def _rb_hessian(x, lagrange, obj_factor):
    # f = (1-x0)^2 + 100*(x1-x0^2)^2
    # d2f/dx0^2 = 2 - 400*x1 + 1200*x0^2
    # d2f/dx0dx1 = -400*x0
    # d2f/dx1^2 = 200
    return np.array([
        obj_factor * (2.0 - 400.0 * x[1] + 1200.0 * x[0]**2),
        obj_factor * (-400.0 * x[0]),
        obj_factor * 200.0,
    ])


def _rb_hessianstructure():
    return (np.array([0, 1, 1]), np.array([0, 0, 1]))


def rosenbrock_factory(intermediate_cb=None):
    prob = cyipopt.Problem(
        n=_ROSENBROCK_N, m=_ROSENBROCK_M,
        problem_obj=_make_problem_obj(
            _rb_objective, _rb_gradient, _rb_constraints, _rb_jacobian,
            _rb_jacobianstructure, _rb_hessian, _rb_hessianstructure,
            intermediate_cb,
        ),
        lb=np.full(_ROSENBROCK_N, -5.0),
        ub=np.full(_ROSENBROCK_N, 5.0),
        cl=np.array([]),
        cu=np.array([]),
    )
    prob.add_option("mu_strategy", "adaptive")
    prob.add_option("tol", 1e-8)
    return prob, _ROSENBROCK_X0.copy()


# ---------------------------------------------------------------------------
# Problem 2: HS071  (4 vars, 2 constraints)
#   min x1*x4*(x1+x2+x3) + x3
#   s.t. x1*x2*x3*x4 >= 25          (inequality)
#        x1^2+x2^2+x3^2+x4^2 = 40   (equality)
#        1 <= xi <= 5
#   x0 = (1, 5, 5, 1)
# ---------------------------------------------------------------------------

_HS071_N = 4
_HS071_M = 2
_HS071_X0 = np.array([1.0, 5.0, 5.0, 1.0])


def _hs071_objective(x):
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]


def _hs071_gradient(x):
    return np.array([
        x[3] * (2.0 * x[0] + x[1] + x[2]),
        x[0] * x[3],
        x[0] * x[3] + 1.0,
        x[0] * (x[0] + x[1] + x[2]),
    ])


def _hs071_constraints(x):
    return np.array([
        x[0] * x[1] * x[2] * x[3],
        x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2,
    ])


def _hs071_jacobian(x):
    return np.array([
        x[1]*x[2]*x[3], x[0]*x[2]*x[3], x[0]*x[1]*x[3], x[0]*x[1]*x[2],
        2.0*x[0],        2.0*x[1],        2.0*x[2],        2.0*x[3],
    ])


def _hs071_jacobianstructure():
    return (
        np.array([0, 0, 0, 0, 1, 1, 1, 1]),
        np.array([0, 1, 2, 3, 0, 1, 2, 3]),
    )


def _hs071_hessian(x, lagrange, obj_factor):
    h = np.zeros(10)
    # Objective Hessian (lower triangle, 4x4)
    h[0] = obj_factor * 2.0 * x[3]                       # (0,0)
    h[1] = obj_factor * x[3]                              # (1,0)
    h[2] = 0.0                                            # (1,1)
    h[3] = obj_factor * x[3]                              # (2,0)
    h[4] = 0.0                                            # (2,1)
    h[5] = 0.0                                            # (2,2)
    h[6] = obj_factor * (2.0 * x[0] + x[1] + x[2])       # (3,0)
    h[7] = obj_factor * x[0]                              # (3,1)
    h[8] = obj_factor * x[0]                              # (3,2)
    h[9] = 0.0                                            # (3,3)

    # Constraint 0: g0 = x0*x1*x2*x3
    lam0 = lagrange[0]
    h[1] += lam0 * x[2] * x[3]   # (1,0)
    h[3] += lam0 * x[1] * x[3]   # (2,0)
    h[4] += lam0 * x[0] * x[3]   # (2,1)
    h[6] += lam0 * x[1] * x[2]   # (3,0)
    h[7] += lam0 * x[0] * x[2]   # (3,1)
    h[8] += lam0 * x[0] * x[1]   # (3,2)

    # Constraint 1: g1 = x0^2+x1^2+x2^2+x3^2
    lam1 = lagrange[1]
    h[0] += lam1 * 2.0
    h[2] += lam1 * 2.0
    h[5] += lam1 * 2.0
    h[9] += lam1 * 2.0

    return h


def _hs071_hessianstructure():
    return (
        np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3]),
        np.array([0, 0, 1, 0, 1, 2, 0, 1, 2, 3]),
    )


def hs071_factory(intermediate_cb=None):
    prob = cyipopt.Problem(
        n=_HS071_N, m=_HS071_M,
        problem_obj=_make_problem_obj(
            _hs071_objective, _hs071_gradient, _hs071_constraints,
            _hs071_jacobian, _hs071_jacobianstructure,
            _hs071_hessian, _hs071_hessianstructure,
            intermediate_cb,
        ),
        lb=np.full(_HS071_N, 1.0),
        ub=np.full(_HS071_N, 5.0),
        cl=np.array([25.0, 40.0]),
        cu=np.array([2.0e19, 40.0]),
    )
    prob.add_option("mu_strategy", "adaptive")
    prob.add_option("tol", 1e-8)
    return prob, _HS071_X0.copy()


# ---------------------------------------------------------------------------
# Problem 3: Simple QP  (2 vars, 1 equality constraint)
#   min 0.5*(x1^2 + x2^2)
#   s.t. x1 + x2 = 1
#   x0 = (0, 0)
# ---------------------------------------------------------------------------

_SQP_N = 2
_SQP_M = 1
_SQP_X0 = np.array([0.0, 0.0])


def _sqp_objective(x):
    return 0.5 * (x[0]**2 + x[1]**2)


def _sqp_gradient(x):
    return np.array([x[0], x[1]])


def _sqp_constraints(x):
    return np.array([x[0] + x[1]])


def _sqp_jacobian(x):
    return np.array([1.0, 1.0])


def _sqp_jacobianstructure():
    return (np.array([0, 0]), np.array([0, 1]))


def _sqp_hessian(x, lagrange, obj_factor):
    return np.array([obj_factor * 1.0, 0.0, obj_factor * 1.0])


def _sqp_hessianstructure():
    return (np.array([0, 1, 1]), np.array([0, 0, 1]))


def simple_qp_factory(intermediate_cb=None):
    prob = cyipopt.Problem(
        n=_SQP_N, m=_SQP_M,
        problem_obj=_make_problem_obj(
            _sqp_objective, _sqp_gradient, _sqp_constraints,
            _sqp_jacobian, _sqp_jacobianstructure,
            _sqp_hessian, _sqp_hessianstructure,
            intermediate_cb,
        ),
        lb=np.full(_SQP_N, -1e19),
        ub=np.full(_SQP_N, 1e19),
        cl=np.array([1.0]),
        cu=np.array([1.0]),
    )
    prob.add_option("mu_strategy", "adaptive")
    prob.add_option("tol", 1e-8)
    return prob, _SQP_X0.copy()


# ---------------------------------------------------------------------------
# Problem 4: HS035  (3 vars, 1 inequality constraint, bounds x >= 0)
#   min 9 - 8x1 - 6x2 - 4x3 + 2x1^2 + 2x2^2 + x3^2 + 2x1*x2 + 2x1*x3
#   s.t. x1 + x2 + 2*x3 <= 3
#        x >= 0
#   x0 = (0.5, 0.5, 0.5)
#   Known optimal: f* = 1/9 at x* = (4/3, 7/9, 4/9)
# ---------------------------------------------------------------------------

_HS035_N = 3
_HS035_M = 1
_HS035_X0 = np.array([0.5, 0.5, 0.5])


def _hs035_objective(x):
    return (9.0 - 8.0*x[0] - 6.0*x[1] - 4.0*x[2]
            + 2.0*x[0]**2 + 2.0*x[1]**2 + x[2]**2
            + 2.0*x[0]*x[1] + 2.0*x[0]*x[2])


def _hs035_gradient(x):
    return np.array([
        -8.0 + 4.0*x[0] + 2.0*x[1] + 2.0*x[2],
        -6.0 + 4.0*x[1] + 2.0*x[0],
        -4.0 + 2.0*x[2] + 2.0*x[0],
    ])


def _hs035_constraints(x):
    return np.array([x[0] + x[1] + 2.0*x[2]])


def _hs035_jacobian(x):
    return np.array([1.0, 1.0, 2.0])


def _hs035_jacobianstructure():
    return (np.array([0, 0, 0]), np.array([0, 1, 2]))


def _hs035_hessian(x, lagrange, obj_factor):
    h = np.zeros(6)
    h[0] = obj_factor * 4.0    # (0,0)
    h[1] = obj_factor * 2.0    # (1,0)
    h[2] = obj_factor * 4.0    # (1,1)
    h[3] = obj_factor * 2.0    # (2,0)
    h[4] = 0.0                 # (2,1)
    h[5] = obj_factor * 2.0    # (2,2)
    # Constraint is linear => no Hessian contribution
    return h


def _hs035_hessianstructure():
    return (
        np.array([0, 1, 1, 2, 2, 2]),
        np.array([0, 0, 1, 0, 1, 2]),
    )


def hs035_factory(intermediate_cb=None):
    prob = cyipopt.Problem(
        n=_HS035_N, m=_HS035_M,
        problem_obj=_make_problem_obj(
            _hs035_objective, _hs035_gradient, _hs035_constraints,
            _hs035_jacobian, _hs035_jacobianstructure,
            _hs035_hessian, _hs035_hessianstructure,
            intermediate_cb,
        ),
        lb=np.zeros(_HS035_N),
        ub=np.full(_HS035_N, 1e19),
        cl=np.array([-1e19]),
        cu=np.array([3.0]),
    )
    prob.add_option("mu_strategy", "adaptive")
    prob.add_option("tol", 1e-8)
    return prob, _HS035_X0.copy()


# ---------------------------------------------------------------------------
# Problem 5: Bound-constrained quadratic  (4 vars, 0 general constraints)
#   min  sum_i (x_i - t_i)^2  where t = [1, 2, 3, 4]
#   bounds: 0 <= x_i <= 3
#   x0 = (0.5, 0.5, 0.5, 0.5)
#   Optimal: x* = (1, 2, 3, 3), f* = 1  (x4 clamped at upper bound)
# ---------------------------------------------------------------------------

_BQ_N = 4
_BQ_M = 0
_BQ_TARGET = np.array([1.0, 2.0, 3.0, 4.0])
_BQ_X0 = np.array([0.5, 0.5, 0.5, 0.5])


def _bq_objective(x):
    return float(np.sum((x - _BQ_TARGET)**2))


def _bq_gradient(x):
    return 2.0 * (x - _BQ_TARGET)


def _bq_constraints(x):
    return np.array([])


def _bq_jacobian(x):
    return np.array([])


def _bq_jacobianstructure():
    return (np.array([], dtype=int), np.array([], dtype=int))


def _bq_hessian(x, lagrange, obj_factor):
    h = np.zeros(10)
    h[0] = obj_factor * 2.0   # (0,0)
    h[2] = obj_factor * 2.0   # (1,1)
    h[5] = obj_factor * 2.0   # (2,2)
    h[9] = obj_factor * 2.0   # (3,3)
    return h


def _bq_hessianstructure():
    return (
        np.array([0, 1, 1, 2, 2, 2, 3, 3, 3, 3]),
        np.array([0, 0, 1, 0, 1, 2, 0, 1, 2, 3]),
    )


def bound_quad_factory(intermediate_cb=None):
    prob = cyipopt.Problem(
        n=_BQ_N, m=_BQ_M,
        problem_obj=_make_problem_obj(
            _bq_objective, _bq_gradient, _bq_constraints,
            _bq_jacobian, _bq_jacobianstructure,
            _bq_hessian, _bq_hessianstructure,
            intermediate_cb,
        ),
        lb=np.zeros(_BQ_N),
        ub=np.full(_BQ_N, 3.0),
        cl=np.array([]),
        cu=np.array([]),
    )
    prob.add_option("mu_strategy", "adaptive")
    prob.add_option("tol", 1e-8)
    return prob, _BQ_X0.copy()


# ---------------------------------------------------------------------------
# Constraint-violation and optimality-error helpers
# ---------------------------------------------------------------------------

def compute_constraint_violation(x, constraint_fn, cl, cu):
    """Max violation of cl <= g(x) <= cu."""
    if len(cl) == 0:
        return 0.0
    g = constraint_fn(x)
    viol = 0.0
    for i, gi in enumerate(g):
        if gi < cl[i] - 1e-20:
            viol = max(viol, cl[i] - gi)
        if gi > cu[i] + 1e-20:
            viol = max(viol, gi - cu[i])
    return viol


def compute_optimality_error(x, grad_fn, mult_g, jac_fn, jac_struct_fn,
                             mult_x_L, mult_x_U, n, m):
    """
    Approximate KKT optimality error using cyipopt's sign convention:
      || grad_f + J^T * lambda - z_L + z_U ||_inf

    cyipopt returns multipliers such that the Lagrangian stationarity
    condition is:  grad_f + J^T * mult_g - mult_x_L + mult_x_U = 0
    """
    g = grad_fn(x)
    kkt = g.copy()

    if m > 0:
        rows, cols = jac_struct_fn()
        jac_vals = jac_fn(x)
        lam = np.array(mult_g)
        # + J^T * lambda
        for idx in range(len(rows)):
            kkt[cols[idx]] += lam[rows[idx]] * jac_vals[idx]

    kkt -= np.array(mult_x_L)
    kkt += np.array(mult_x_U)
    return float(np.max(np.abs(kkt)))


# ---------------------------------------------------------------------------
# Solver wrapper
# ---------------------------------------------------------------------------

def solve_and_record(name, factory, n_vars, n_constraints,
                     constraint_fn, cl, cu, grad_fn, jac_fn, jac_struct_fn,
                     n_timing_runs=20):
    """
    Solve via cyipopt and collect all diagnostics.
    """
    # -- authoritative solve ------------------------------------------------
    prob, x0 = factory()
    prob.add_option("print_level", 0)
    prob.add_option("sb", "yes")

    x_opt, info = prob.solve(x0)

    obj = float(info["obj_val"])
    mult_g = info["mult_g"].tolist()
    mult_x_L = info["mult_x_L"].tolist()
    mult_x_U = info["mult_x_U"].tolist()
    status = int(info["status"])
    status_msg = info["status_msg"]
    if isinstance(status_msg, bytes):
        status_msg = status_msg.decode()
    else:
        status_msg = str(status_msg)

    # -- iteration count via intermediate callback --------------------------
    iter_counter = {"count": 0}

    def counting_intermediate(*args, **kwargs):
        iter_counter["count"] += 1
        return True

    prob2, x0_2 = factory(intermediate_cb=counting_intermediate)
    prob2.add_option("print_level", 0)
    prob2.add_option("sb", "yes")
    prob2.solve(x0_2)
    n_iters = iter_counter["count"]

    # -- timing (average over multiple runs) --------------------------------
    times = []
    for _ in range(n_timing_runs):
        p, x0_t = factory()
        p.add_option("print_level", 0)
        p.add_option("sb", "yes")
        t0 = time.perf_counter()
        p.solve(x0_t)
        t1 = time.perf_counter()
        times.append(t1 - t0)

    avg_time = float(np.mean(times))
    std_time = float(np.std(times))

    # -- constraint violation -----------------------------------------------
    cv = compute_constraint_violation(np.array(x_opt), constraint_fn, cl, cu)

    # -- optimality error (KKT) ---------------------------------------------
    oe = compute_optimality_error(
        np.array(x_opt), grad_fn, mult_g, jac_fn, jac_struct_fn,
        mult_x_L, mult_x_U, n_vars, n_constraints,
    )

    return {
        "problem": name,
        "n_vars": n_vars,
        "n_constraints": n_constraints,
        "x_opt": [float(v) for v in x_opt],
        "obj_opt": obj,
        "mult_g": mult_g,
        "mult_x_L": mult_x_L,
        "mult_x_U": mult_x_U,
        "n_iterations": n_iters,
        "solve_time_avg_s": avg_time,
        "solve_time_std_s": std_time,
        "n_timing_runs": n_timing_runs,
        "constraint_violation": cv,
        "optimality_error": oe,
        "status": status,
        "status_msg": status_msg,
    }


# ---------------------------------------------------------------------------
# Problem registry
# ---------------------------------------------------------------------------

PROBLEMS = [
    {
        "name": "rosenbrock",
        "factory": rosenbrock_factory,
        "n_vars": _ROSENBROCK_N,
        "n_constraints": _ROSENBROCK_M,
        "constraint_fn": _rb_constraints,
        "cl": np.array([]),
        "cu": np.array([]),
        "grad_fn": _rb_gradient,
        "jac_fn": _rb_jacobian,
        "jac_struct_fn": _rb_jacobianstructure,
    },
    {
        "name": "hs071",
        "factory": hs071_factory,
        "n_vars": _HS071_N,
        "n_constraints": _HS071_M,
        "constraint_fn": _hs071_constraints,
        "cl": np.array([25.0, 40.0]),
        "cu": np.array([2.0e19, 40.0]),
        "grad_fn": _hs071_gradient,
        "jac_fn": _hs071_jacobian,
        "jac_struct_fn": _hs071_jacobianstructure,
    },
    {
        "name": "simple_qp",
        "factory": simple_qp_factory,
        "n_vars": _SQP_N,
        "n_constraints": _SQP_M,
        "constraint_fn": _sqp_constraints,
        "cl": np.array([1.0]),
        "cu": np.array([1.0]),
        "grad_fn": _sqp_gradient,
        "jac_fn": _sqp_jacobian,
        "jac_struct_fn": _sqp_jacobianstructure,
    },
    {
        "name": "hs035",
        "factory": hs035_factory,
        "n_vars": _HS035_N,
        "n_constraints": _HS035_M,
        "constraint_fn": _hs035_constraints,
        "cl": np.array([-1e19]),
        "cu": np.array([3.0]),
        "grad_fn": _hs035_gradient,
        "jac_fn": _hs035_jacobian,
        "jac_struct_fn": _hs035_jacobianstructure,
    },
    {
        "name": "bound_constrained_quadratic",
        "factory": bound_quad_factory,
        "n_vars": _BQ_N,
        "n_constraints": _BQ_M,
        "constraint_fn": _bq_constraints,
        "cl": np.array([]),
        "cu": np.array([]),
        "grad_fn": _bq_gradient,
        "jac_fn": _bq_jacobian,
        "jac_struct_fn": _bq_jacobianstructure,
    },
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results = []

    for pdef in PROBLEMS:
        rec = solve_and_record(
            name=pdef["name"],
            factory=pdef["factory"],
            n_vars=pdef["n_vars"],
            n_constraints=pdef["n_constraints"],
            constraint_fn=pdef["constraint_fn"],
            cl=pdef["cl"],
            cu=pdef["cu"],
            grad_fn=pdef["grad_fn"],
            jac_fn=pdef["jac_fn"],
            jac_struct_fn=pdef["jac_struct_fn"],
            n_timing_runs=20,
        )
        results.append(rec)

    # Print JSON to stdout
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
