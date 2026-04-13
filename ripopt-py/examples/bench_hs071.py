"""Quick wall-time benchmark for the Python interface on HS071.

Runs the solve a fixed number of times in each configuration (warm JIT,
discarding the first trace-time-inflated run) and reports mean and min
wall time.

Usage:
    python examples/bench_hs071.py
"""

import statistics
import time

import jax.numpy as jnp
import numpy as np

from ripopt import minimize


def f(x):
    return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]


def g(x):
    return jnp.array(
        [
            x[0] * x[1] * x[2] * x[3],
            x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2,
        ]
    )


X0 = [1.0, 5.0, 5.0, 1.0]
CONSTRAINTS = {"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]}
BOUNDS = (1.0, 5.0)
TARGET_F = 17.0140173


def run_once(**kwargs):
    t0 = time.perf_counter()
    res = minimize(
        f,
        x0=X0,
        bounds=BOUNDS,
        constraints=CONSTRAINTS,
        options={"tol": 1e-8, "print_level": 0},
        **kwargs,
    )
    wall = time.perf_counter() - t0
    assert res.success, f"expected Optimal, got {res.status}"
    assert abs(res.fun - TARGET_F) < 1e-4, f"f={res.fun}"
    return wall, res.wall_time_secs, res.iterations


def bench(label, n_runs=20, **kwargs):
    # One warm-up so JIT tracing doesn't pollute the first measured run.
    run_once(**kwargs)

    totals = []
    solver_only = []
    iters_list = []
    for _ in range(n_runs):
        total, solver, iters = run_once(**kwargs)
        totals.append(total)
        solver_only.append(solver)
        iters_list.append(iters)

    print(f"{label}")
    print(
        f"  total wall    (mean / min):  {statistics.mean(totals)*1e3:7.2f} ms"
        f"  /  {min(totals)*1e3:7.2f} ms"
    )
    print(
        f"  solver only   (mean / min):  {statistics.mean(solver_only)*1e3:7.2f} ms"
        f"  /  {min(solver_only)*1e3:7.2f} ms"
    )
    print(f"  iterations:                  {iters_list[0]}")
    print()


def main():
    print("HS071 — JAX autodiff + ripopt\n")
    bench("sparsity='dense' (default)", sparsity="dense")
    bench("sparsity='detect'", sparsity="detect")


if __name__ == "__main__":
    main()
