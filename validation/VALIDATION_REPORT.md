# Ripopt vs Cyipopt Validation Report

## Summary

| Problem | Vars | Cons | ripopt Status | cyipopt Status | Obj Match | x Match |
|---------|------|------|---------------|----------------|-----------|---------|
| rosenbrock           | 2 | 0 | Optimal         | Algorithm termi | PASS      | PASS    |
| hs071                | 4 | 2 | Optimal         | Algorithm termi | PASS      | PASS    |
| simpleqp             | 2 | 1 | Optimal         | Algorithm termi | PASS      | PASS    |
| hs035                | 3 | 1 | Acceptable      | Algorithm termi | PASS      | PASS    |
| pureboundconstrained | 4 | 0 | Optimal         | Algorithm termi | PASS      | PASS    |

**Overall: ALL PASS**

## Detailed Comparison

### rosenbrock

| Metric | ripopt | cyipopt | Difference |
|--------|--------|---------|------------|
| Objective | 3.7439756431e-21 | 1.5642801260e-24 | 3.74e-21 |
| x[0] | 0.9999999999 | 1.0000000000 | 5.87e-11 |
| x[1] | 0.9999999999 | 1.0000000000 | 1.19e-10 |
| Iterations | 21 | 21 | 0 |
| Constraint viol. | 0.00e+00 | 0.00e+00 | |

### hs071

| Metric | ripopt | cyipopt | Difference |
|--------|--------|---------|------------|
| Objective | 1.7014017289e+01 | 1.7014017140e+01 | 1.49e-07 |
| x[0] | 1.0000000001 | 0.9999999900 | 1.01e-08 |
| x[1] | 4.7429996368 | 4.7429996436 | 6.81e-09 |
| x[2] | 3.8211499849 | 3.8211499789 | 5.98e-09 |
| x[3] | 1.3794082929 | 1.3794082932 | 3.07e-10 |
| y[0] | -0.5522936599 | -0.5522936595 | 3.99e-10 |
| y[1] | 0.1614685664 | 0.1614685642 | 2.18e-09 |
| Iterations | 17 | 9 | 8 |
| Constraint viol. | 4.32e-10 | 2.50e-07 | |

### simpleqp

| Metric | ripopt | cyipopt | Difference |
|--------|--------|---------|------------|
| Objective | 2.5000000000e-01 | 2.5000000000e-01 | 0.00e+00 |
| x[0] | 0.5000000000 | 0.5000000000 | 0.00e+00 |
| x[1] | 0.5000000000 | 0.5000000000 | 0.00e+00 |
| y[0] | -0.5000000000 | -0.5000000000 | 0.00e+00 |
| Iterations | 1 | 2 | 1 |
| Constraint viol. | 0.00e+00 | 0.00e+00 | |

### hs035

| Metric | ripopt | cyipopt | Difference |
|--------|--------|---------|------------|
| Objective | 1.1111111111e-01 | 1.1111110445e-01 | 6.66e-09 |
| x[0] | 1.3333333287 | 1.3333333233 | 5.37e-09 |
| x[1] | 0.7777777801 | 0.7777777844 | 4.35e-09 |
| x[2] | 0.4444444456 | 0.4444444611 | 1.55e-08 |
| y[0] | 0.2222222222 | 0.2222222156 | 6.67e-09 |
| Iterations | 22 | 8 | 14 |
| Constraint viol. | 0.00e+00 | 3.00e-08 | |

### pureboundconstrained

| Metric | ripopt | cyipopt | Difference |
|--------|--------|---------|------------|
| Objective | 1.0000000014e+00 | 9.9999994336e-01 | 5.81e-08 |
| x[0] | 1.0000000001 | 1.0000000000 | 7.10e-11 |
| x[1] | 1.9999999999 | 2.0000000000 | 7.10e-11 |
| x[2] | 2.9999637141 | 2.9999421524 | 2.16e-05 |
| x[3] | 2.9999999999 | 3.0000000300 | 3.00e-08 |
| Iterations | 17 | 16 | 1 |
| Constraint viol. | 0.00e+00 | 0.00e+00 | |

## Performance Comparison

| Problem | ripopt (ms) | cyipopt (ms) | Speedup |
|---------|-------------|--------------|---------|
| rosenbrock           |    0.013 +/- 0.011 |    4.139 +/- 0.212 | 328.7x |
| hs071                |    0.017 +/- 0.006 |    1.572 +/- 0.083 | 91.9x |
| simpleqp             |    0.001 +/- 0.000 |    0.298 +/- 0.042 | 233.5x |
| hs035                |    0.036 +/- 0.003 |    1.291 +/- 0.078 | 36.2x |
| pureboundconstrained |    0.011 +/- 0.001 |    2.499 +/- 0.184 | 235.8x |

*Note: ripopt is compiled in release mode. cyipopt uses MUMPS linear solver.*
*Timing includes problem setup overhead. Small problems amplify per-call overhead.*

## Test Coverage Comparison with Ipopt

Ipopt's test suite includes the following test categories:

| Category | Ipopt Tests | ripopt Tests | Coverage |
|----------|-------------|--------------|----------|
| Unconstrained NLP | HS-collection | Rosenbrock | Partial |
| Equality constrained | HS-collection | SimpleQP, MultipleEquality | Partial |
| Inequality constrained | HS-collection | HS071, HS035 | Partial |
| Bound constrained | Various | PureBoundConstrained | Partial |
| Linear solver | MA27/MA57/MUMPS | Dense LDL (Bunch-Kaufman) | Different |
| Warm start | Yes | Basic (untested) | Minimal |
| Restoration phase | Full feasibility restoration | Simplified gradient projection | Partial |
| Inertia correction | Full with multiple strategies | Basic growth factor | Partial |
| Scaling | NLP scaling | None | Missing |
| Iterative refinement | Yes | No | Missing |

### Ipopt features NOT yet implemented in ripopt:

- Sparse linear solvers (MA27, MA57, MUMPS, Pardiso)
- NLP scaling (gradient-based, user-provided)
- Iterative refinement for KKT solves
- Full restoration phase (feasibility restoration NLP)
- Adaptive barrier strategies (Mehrotra predictor-corrector)
- Quasi-Newton Hessian approximation (L-BFGS)
- HSL and SPRAL linear solver interfaces
- Pardiso and WSMP linear solver interfaces
- sIpopt (sensitivity analysis)

### ripopt unit tests (17 total):

- Dense LDL factorization: 7 tests (identity, positive definite, negative definite,
  indefinite, KKT-shaped, larger KKT, near-singular)
- Convergence: 4 tests (feasible, violated, optimal, not converged)
- Filter: 5 tests (empty filter, theta_max rejection, dominated point,
  reset, switching condition)
- Fraction-to-boundary: 1 test

### ripopt integration tests (6 total):

- Rosenbrock (unconstrained, variable bounds)
- SimpleQP (equality constraint)
- HS071 (inequality + equality constraints)
- HS035 / BoundConstrainedQuadratic (inequality + non-negative bounds)
- PureBoundConstrained (box bounds only)
- MultipleEqualityConstraints (2 equality constraints)
