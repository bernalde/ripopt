#!/usr/bin/env python3
"""
Solve all generated HS problems using cyipopt and output JSON results
for comparison with ripopt.
"""

import json
import sys
import time
import numpy as np

# Add generated directory to path
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/generated')
from hs_cyipopt import HS_PROBLEMS, _make_problem_obj


def solve_problem(pdef, n_timing_runs=5):
    """Solve a single HS problem and collect diagnostics."""
    number = pdef['number']
    factory = pdef['factory']
    n = pdef['n']
    m = pdef['m']
    known_fopt = pdef['known_fopt']

    try:
        # Authoritative solve
        prob, x0 = factory()
        prob.add_option("print_level", 0)
        prob.add_option("sb", "yes")

        x_opt, info = prob.solve(x0)

        obj = float(info["obj_val"])
        mult_g = info["mult_g"].tolist() if m > 0 else []
        mult_x_L = info["mult_x_L"].tolist()
        mult_x_U = info["mult_x_U"].tolist()
        status = int(info["status"])
        status_msg = info["status_msg"]
        if isinstance(status_msg, bytes):
            status_msg = status_msg.decode()
        else:
            status_msg = str(status_msg)

        # Map status to comparable string
        if status == 0:
            status_str = "Optimal"
        elif status == 1:
            status_str = "Acceptable"
        elif status == 2:
            status_str = "Infeasible"
        elif status == -1:
            status_str = "MaxIterations"
        else:
            status_str = f"IpoptStatus({status})"

        # Iteration count via callback
        iter_counter = {"count": 0}
        def counting_intermediate(*args, **kwargs):
            iter_counter["count"] += 1
            return True

        prob2, x0_2 = factory(intermediate_cb=counting_intermediate)
        prob2.add_option("print_level", 0)
        prob2.add_option("sb", "yes")
        prob2.solve(x0_2)
        n_iters = iter_counter["count"]

        # Timing
        times = []
        for _ in range(n_timing_runs):
            p, x0_t = factory()
            p.add_option("print_level", 0)
            p.add_option("sb", "yes")
            t0 = time.perf_counter()
            p.solve(x0_t)
            t1 = time.perf_counter()
            times.append(t1 - t0)
        avg_time = float(np.min(times))

        # Constraint values
        try:
            # Re-evaluate constraints at solution
            prob3, _ = factory()
            constraint_values = prob3.problem_obj.constraints(np.array(x_opt)).tolist()
        except:
            constraint_values = []

        return {
            "number": number,
            "status": status_str,
            "objective": obj,
            "x": [float(v) for v in x_opt],
            "constraint_multipliers": mult_g,
            "bound_multipliers_lower": mult_x_L,
            "bound_multipliers_upper": mult_x_U,
            "constraint_values": constraint_values,
            "iterations": n_iters,
            "solve_time": avg_time,
            "known_fopt": known_fopt,
            "n": n,
            "m": m,
        }

    except Exception as e:
        return {
            "number": number,
            "status": f"Error: {e}",
            "objective": float('nan'),
            "x": [],
            "constraint_multipliers": [],
            "bound_multipliers_lower": [],
            "bound_multipliers_upper": [],
            "constraint_values": [],
            "iterations": 0,
            "solve_time": 0.0,
            "known_fopt": known_fopt,
            "n": n,
            "m": m,
        }


def main():
    print(f"Solving {len(HS_PROBLEMS)} HS problems with cyipopt...", file=sys.stderr)

    results = []
    for pdef in HS_PROBLEMS:
        rec = solve_problem(pdef)
        results.append(rec)
        status = rec['status']
        number = rec['number']
        print(f"  TP{number:03d}: {status}", file=sys.stderr)

    # Summary
    total = len(results)
    optimal = sum(1 for r in results if r['status'] == 'Optimal')
    acceptable = sum(1 for r in results if r['status'] == 'Acceptable')
    solved = optimal + acceptable
    print(f"Solved {solved}/{total} ({optimal} optimal, {acceptable} acceptable)",
          file=sys.stderr)

    # Output JSON
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
