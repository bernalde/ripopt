"""Minimal reproducer for issue #7: call ripopt directly via discopt's PyO3 bindings.

Creates the exact same NLP that discopt's B&B would create for the exponential
activation problem with y1=1, y2=1 fixed.
"""
import os
os.environ["JAX_PLATFORMS"] = "cpu"
os.environ["JAX_ENABLE_X64"] = "1"

import numpy as np
import jax
import jax.numpy as jnp

# Problem: min exp(x1) + exp(x2) + 5*y1 + 6*y2 + 8*y3
# Constraints (discopt normalization: >= becomes negated <=):
#   g0: 3-x1-x2 <= 0
#   g1: x1-4*y1 <= 0
#   g2: x2-4*y2 <= 0
#   g3: y1+y2-1-y3 <= 0

N = 5  # x1, x2, y1, y2, y3
M = 4

def objective(x):
    return jnp.exp(x[0]) + jnp.exp(x[1]) + 5*x[2] + 6*x[3] + 8*x[4]

def constraints(x):
    return jnp.array([
        3.0 - x[0] - x[1],
        x[0] - 4*x[2],
        x[1] - 4*x[3],
        x[2] + x[3] - 1.0 - x[4],
    ])

def lagrangian(x, obj_factor, lam):
    return obj_factor * objective(x) + jnp.dot(lam, constraints(x))

grad_fn = jax.jit(jax.grad(objective))
jac_fn = jax.jit(jax.jacobian(constraints))
hess_fn = jax.jit(jax.hessian(lagrangian, argnums=0))

class Evaluator:
    def __init__(self, lb, ub):
        self._lb = lb
        self._ub = ub

    @property
    def n_variables(self):
        return N

    @property
    def n_constraints(self):
        return M

    @property
    def variable_bounds(self):
        return self._lb, self._ub

    def evaluate_objective(self, x):
        return float(objective(jnp.array(x)))

    def evaluate_gradient(self, x):
        return np.asarray(grad_fn(jnp.array(x)))

    def evaluate_constraints(self, x):
        return np.asarray(constraints(jnp.array(x)))

    def evaluate_jacobian(self, x):
        return np.asarray(jac_fn(jnp.array(x)))

    def evaluate_lagrangian_hessian(self, x, obj_factor, lam):
        return np.asarray(hess_fn(jnp.array(x), obj_factor, jnp.array(lam)))

# Test subproblems (matching discopt B&B nodes)
configs = [
    ("Root relaxation", [0,0,0,0,0], [4,4,1,1,1], [2,2,0.5,0.5,0.5]),
    ("y1=1", [0,0,1,0,0], [4,4,1,1,1], [2,2,1,0.5,0.5]),
    ("y1=0", [0,0,0,0,0], [4,4,0,1,1], [0,2,0,0.5,0.5]),
    ("y1=1,y2=1", [0,0,1,1,0], [4,4,1,1,1], [1.5,1.5,1,1,0.5]),
    ("y1=1,y2=0", [0,0,1,0,0], [4,4,1,0,1], [3,0,1,0,0.5]),
    ("y1=0,y2=1", [0,0,0,1,0], [4,4,0,1,1], [0,3,0,1,0.5]),
    ("y1=0,y2=1,y3=0", [0,0,0,1,0], [4,4,0,1,0], [0,3,0,1,0]),
]

from discopt._rust import solve_ripopt

for name, lb, ub, x0 in configs:
    lb = np.array(lb, dtype=np.float64)
    ub = np.array(ub, dtype=np.float64)
    x0 = np.array(x0, dtype=np.float64)
    g_l = np.full(M, -np.inf)
    g_u = np.zeros(M)

    ev = Evaluator(lb, ub)
    opts = {"print_level": 3, "tol": 1e-7, "max_iter": 3000}

    result = solve_ripopt(ev, x0, lb, ub, g_l, g_u, opts)
    print(f"{name:25s}: status={result['status']:15s} iters={result['iterations']:4d} obj={result['objective']:.6f}")
