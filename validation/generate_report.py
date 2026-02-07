#!/usr/bin/env python3
"""Generate a validation report comparing ripopt and cyipopt results."""
import json
import sys

def load_results(path):
    with open(path) as f:
        data = json.load(f)
    return {r["problem"].lower().replace(" ", "_"): r for r in data}

def main():
    ripopt = load_results("validation/ripopt_results.json")
    cyipopt = load_results("validation/cyipopt_results.json")

    # Map problem names between the two
    name_map = {
        "rosenbrock": "rosenbrock",
        "hs071": "hs071",
        "simpleqp": "simple_qp",
        "hs035": "hs035",
        "pureboundconstrained": "bound_constrained_quadratic",
    }

    print("# Ripopt vs Cyipopt Validation Report")
    print()
    print("## Summary")
    print()
    print("| Problem | Vars | Cons | ripopt Status | cyipopt Status | Obj Match | x Match |")
    print("|---------|------|------|---------------|----------------|-----------|---------|")

    all_pass = True
    for rname, cname in name_map.items():
        r = ripopt.get(rname)
        c = cyipopt.get(cname)
        if not r or not c:
            print(f"| {rname} | - | - | MISSING | MISSING | - | - |")
            all_pass = False
            continue

        obj_diff = abs(r["obj_opt"] - c["obj_opt"])
        obj_match = obj_diff < 1e-4
        x_r = r["x_opt"]
        x_c = c["x_opt"]
        x_diff = max(abs(a - b) for a, b in zip(x_r, x_c))
        x_match = x_diff < 1e-3

        obj_sym = "PASS" if obj_match else f"FAIL ({obj_diff:.2e})"
        x_sym = "PASS" if x_match else f"FAIL ({x_diff:.2e})"

        if not (obj_match and x_match):
            all_pass = False

        print(f"| {rname:20s} | {r['n_vars']} | {r['n_constraints']} | {r['status']:15s} | {c['status_msg'][:15]:15s} | {obj_sym:9s} | {x_sym:7s} |")

    print()
    print(f"**Overall: {'ALL PASS' if all_pass else 'SOME FAILURES'}**")

    print()
    print("## Detailed Comparison")
    print()

    for rname, cname in name_map.items():
        r = ripopt.get(rname)
        c = cyipopt.get(cname)
        if not r or not c:
            continue

        print(f"### {rname}")
        print()
        print(f"| Metric | ripopt | cyipopt | Difference |")
        print(f"|--------|--------|---------|------------|")

        obj_diff = abs(r["obj_opt"] - c["obj_opt"])
        print(f"| Objective | {r['obj_opt']:.10e} | {c['obj_opt']:.10e} | {obj_diff:.2e} |")

        for i in range(r["n_vars"]):
            xr = r["x_opt"][i]
            xc = c["x_opt"][i]
            print(f"| x[{i}] | {xr:.10f} | {xc:.10f} | {abs(xr-xc):.2e} |")

        # Constraint multipliers
        mult_key_r = "constraint_multipliers"
        mult_key_c = "mult_g"
        mr = r.get(mult_key_r, [])
        mc = c.get(mult_key_c, [])
        for i in range(min(len(mr), len(mc))):
            print(f"| y[{i}] | {mr[i]:.10f} | {mc[i]:.10f} | {abs(mr[i]-mc[i]):.2e} |")

        print(f"| Iterations | {r['iterations']} | {c['n_iterations']} | {abs(r['iterations']-c['n_iterations'])} |")
        print(f"| Constraint viol. | {r['constraint_violation']:.2e} | {c['constraint_violation']:.2e} | |")
        print()

    print("## Performance Comparison")
    print()
    print("| Problem | ripopt (ms) | cyipopt (ms) | Speedup |")
    print("|---------|-------------|--------------|---------|")

    for rname, cname in name_map.items():
        r = ripopt.get(rname)
        c = cyipopt.get(cname)
        if not r or not c:
            continue

        rt = r["solve_time_avg_s"] * 1000
        ct = c["solve_time_avg_s"] * 1000
        speedup = ct / rt if rt > 0 else float('inf')
        print(f"| {rname:20s} | {rt:8.3f} +/- {r['solve_time_std_s']*1000:.3f} | {ct:8.3f} +/- {c['solve_time_std_s']*1000:.3f} | {speedup:.1f}x |")

    print()
    print("*Note: ripopt is compiled in release mode. cyipopt uses MUMPS linear solver.*")
    print("*Timing includes problem setup overhead. Small problems amplify per-call overhead.*")

    print()
    print("## Test Coverage Comparison with Ipopt")
    print()
    print("Ipopt's test suite includes the following test categories:")
    print()
    print("| Category | Ipopt Tests | ripopt Tests | Coverage |")
    print("|----------|-------------|--------------|----------|")
    print("| Unconstrained NLP | HS-collection | Rosenbrock | Partial |")
    print("| Equality constrained | HS-collection | SimpleQP, MultipleEquality | Partial |")
    print("| Inequality constrained | HS-collection | HS071, HS035 | Partial |")
    print("| Bound constrained | Various | PureBoundConstrained | Partial |")
    print("| Linear solver | MA27/MA57/MUMPS | Dense LDL (Bunch-Kaufman) | Different |")
    print("| Warm start | Yes | Basic (untested) | Minimal |")
    print("| Restoration phase | Full feasibility restoration | Simplified gradient projection | Partial |")
    print("| Inertia correction | Full with multiple strategies | Basic growth factor | Partial |")
    print("| Scaling | NLP scaling | None | Missing |")
    print("| Iterative refinement | Yes | No | Missing |")
    print()
    print("### Ipopt features NOT yet implemented in ripopt:")
    print()
    print("- Sparse linear solvers (MA27, MA57, MUMPS, Pardiso)")
    print("- NLP scaling (gradient-based, user-provided)")
    print("- Iterative refinement for KKT solves")
    print("- Full restoration phase (feasibility restoration NLP)")
    print("- Adaptive barrier strategies (Mehrotra predictor-corrector)")
    print("- Quasi-Newton Hessian approximation (L-BFGS)")
    print("- HSL and SPRAL linear solver interfaces")
    print("- Pardiso and WSMP linear solver interfaces")
    print("- sIpopt (sensitivity analysis)")
    print()
    print("### ripopt unit tests (17 total):")
    print()
    print("- Dense LDL factorization: 7 tests (identity, positive definite, negative definite,")
    print("  indefinite, KKT-shaped, larger KKT, near-singular)")
    print("- Convergence: 4 tests (feasible, violated, optimal, not converged)")
    print("- Filter: 5 tests (empty filter, theta_max rejection, dominated point,")
    print("  reset, switching condition)")
    print("- Fraction-to-boundary: 1 test")
    print()
    print("### ripopt integration tests (6 total):")
    print()
    print("- Rosenbrock (unconstrained, variable bounds)")
    print("- SimpleQP (equality constraint)")
    print("- HS071 (inequality + equality constraints)")
    print("- HS035 / BoundConstrainedQuadratic (inequality + non-negative bounds)")
    print("- PureBoundConstrained (box bounds only)")
    print("- MultipleEqualityConstraints (2 equality constraints)")


if __name__ == "__main__":
    main()
