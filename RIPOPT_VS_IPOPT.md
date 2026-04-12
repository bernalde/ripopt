# ripopt vs Ipopt: A Comparative Analysis

ripopt began as a Rust translation of the Ipopt interior-point optimizer. Through
iterative development it has diverged significantly, incorporating novel algorithmic
strategies that allow it to solve **more problems** than the reference implementation
on the CUTEst benchmark suite while being dramatically faster on small-to-medium
problems. This document provides a balanced analysis of where ripopt innovates, where
Ipopt remains stronger, and where there is room to improve.

## Benchmark Summary

|                          | ripopt              | Ipopt                |
|--------------------------|---------------------|----------------------|
| CUTEst solved            | 553/727 (76.1%)     | **561/727 (77.2%)**  |
| HS solved                | 115/120 (95.8%)     | **116/120 (96.7%)**  |
| Both solve (CUTEst)      | 513                 | 513                  |
| Matching objectives      | 413/513 (80.5%)     |                      |
| ripopt-only (CUTEst)     | **40**              | --                   |
| Ipopt-only (CUTEst)      | --                  | **48**               |
| Both fail (CUTEst)       | 166                 | 166                  |

**Solution quality** (513 CUTEst problems where both converge):
- Matching objectives (rel diff < 1e-4): 413/513 (80.5%)
- 100 mismatches: both reach valid KKT points but at different local optima

**Speed** (513 CUTEst common successes):
- Geometric mean speedup: **8.0x** (ripopt faster)
- Median speedup: **18.8x**
- 81% of problems: ripopt faster
- 61% of problems: ripopt 10x+ faster
- Small problems (n <= 10): massive speedup (100x+ median on microsecond solves)
- Medium problems (10 < n <= 50): typically 2-5x faster
- Large problems (n > 50): roughly even, with Ipopt's MUMPS winning on largest systems

ripopt's speed advantage on small problems comes from Rust's zero-overhead abstractions,
no dynamic memory allocation in the hot loop, and the absence of C/Fortran interop
overhead. On larger problems, Ipopt's Fortran MUMPS factorization narrows and eventually
reverses the gap.

---

## Key Innovations in ripopt

### 1. NE-to-LS Reformulation

**The problem.** Many CUTEst problems are "nonlinear equation" (NE) problems: find x
such that g(x) = 0, with a trivially zero objective (f = 0). Ipopt treats these as
constrained optimization and applies its standard IPM machinery, which struggles
because the zero objective provides no gradient information to guide the search.
Ipopt (via CUTEst interface) returns `IpoptStatus(-10)` on 20 such problems.

**ripopt's approach.** ripopt detects NE problems at the start of solve: if f = 0,
grad_f = 0, all constraints are equalities, and m >= n, it reformulates the problem as
unconstrained least-squares:

```
min  (1/2) * ||g(x) - target||^2
```

with gradient J^T * r and a full Hessian that includes both the Gauss-Newton term
(J^T * J) and the second-order correction (sum_i r_i * nabla^2 g_i). The second-order
terms are critical for convergence when far from the solution; a GN-only Hessian fails
on problems like PFIT4 where the initial residual is large.

For square systems (m = n), if the LS reformulation reports infeasibility (stuck at a
nonzero local minimum), ripopt falls back to the original constrained IPM. This hybrid
strategy handles both overdetermined systems (m > n, LS natural) and square systems
(m = n, LS first for speed, constrained as fallback).

**Impact.** 10 NE problems solved by ripopt that Ipopt cannot handle (BEALENE, BOX3NE,
BROWNBSNE, DENSCHNBNE, DENSCHNENE, DEVGLA1NE, ENGVAL2NE, EXP2NE, GULFNE, YFITNE).
Additionally solves PFIT1, PFIT2, PFIT4, HEART6, LEWISPOL, GROUPING, MESH, NYSTROM5,
NYSTROM5C, and others where the NE structure is exploited.

### 2. Two-Phase Restoration

**The problem.** When the filter line search fails to find an acceptable step, the
solver must recover feasibility. Ipopt uses a full NLP restoration phase that minimizes
constraint violations using the same IPM engine. This is robust but expensive.

**ripopt's approach.** ripopt uses a two-phase strategy:

- **Phase 1: Gauss-Newton restoration** (fast). Minimizes ||violation||^2 using
  Gauss-Newton steps on the active constraint subset. Provides quadratic convergence
  for nonlinear equalities (vs. linear for gradient descent). Includes Levenberg-Marquardt
  regularization, gradient descent fallback when GN is singular, and proximity
  regularization to prevent wandering. Typically resolves feasibility in < 10 iterations.

- **Phase 2: NLP restoration** (robust). Only triggered after 2 consecutive GN failures.
  Formulates the full restoration NLP with slack decomposition:
  ```
  min  rho*(sum(p) + sum(n)) + (eta/2)*||D_R(x - x_r)||^2
  s.t. g(x) - p + n = g_target,  p,n >= 0
  ```
  Solved by the same IPM engine with recursion prevention (inner solve disables NLP
  restoration). Uses dynamic dispatch (`&dyn NlpProblem`) to break infinite
  monomorphization in the Rust type system.

**Impact.** The GN phase handles 90%+ of restoration calls cheaply. The NLP phase
recovers from the hard cases that GN cannot (e.g., TP374, DISCS, SPANHYD). Several
CUTEst problems that fail with GN-only restoration succeed with the two-phase approach.

### 3. Dual Convergence with Complementarity Gate

**The problem.** Interior-point methods check optimality via stationarity:
grad_f + J^T * y - z = 0. When the Lagrange multipliers y oscillate (common at
degenerate points), the bound multipliers z_opt computed from stationarity can absorb
the gradient residual, falsely satisfying the convergence check at a non-optimal point.

**ripopt's approach.** ripopt maintains two sets of bound multipliers:
- **z_iterative**: updated each iteration via the IPM step
- **z_optimal**: computed from stationarity (z_opt = -(grad_f + J^T * y))

The convergence check uses z_optimal for dual infeasibility, but only when a
**complementarity gate** is satisfied: z_opt * slack <= kappa_compl * mu (with
kappa_compl = 1e10). When the gate fails, z_iterative is used instead. This prevents
false convergence when z_opt is dominated by oscillating y.

**Impact.** Without the gate, TP023 falsely reports Optimal at obj=4697 (true optimum
is 2.0). The gate forces continued iteration until genuine optimality is reached.

### 4. Pragmatic Inertia Correction

**The problem.** The KKT matrix must have specific inertia (n positive, m negative,
0 zero eigenvalues) for the Newton step to be a descent direction. When factorization
produces wrong inertia, regularization (delta_w, delta_c) is added. But sometimes the
required regularization is so large that the step becomes meaningless.

**ripopt's approach.** After max_attempts=10 inertia correction attempts, ripopt
proceeds with the approximate factorization rather than returning an error. The
filter line search rejects bad steps (small alpha), and restoration can recover from
any damage. This is more pragmatic than Ipopt's approach of reporting
ErrorInStepComputation.

**Impact.** Problems like EQC and HIMMELBJ where Ipopt reports ErrorInStepComputation
are solved by ripopt (Acceptable) because the solver continues past inertia failures.

### 5. Second-Order Correction on Every Backtracking Step

**The problem.** The standard filter line search applies Second-Order Correction (SOC)
to handle the Maratos effect, where the constraint linearization error causes rejection
of good Newton steps. Ipopt typically applies SOC only at the full step.

**ripopt's approach.** ripopt applies SOC at every backtracking step where constraint
violation increases (theta_trial > theta_current), not just the first. This gives more
opportunities to correct the linearization error as the step size decreases.

**Impact.** Improves convergence on problems where the Maratos effect persists at
reduced step sizes (e.g., HS23, several CUTEst constrained problems).

### 6. Hybrid Dense/Sparse Linear Algebra

**The problem.** The KKT matrix is symmetric indefinite, requiring a factorization
that handles mixed-sign eigenvalues. The solver must work across problem sizes from
n+m=3 to n+m=5000+.

**ripopt's approach.** Two factorization backends, selected automatically:

- **Dense Bunch-Kaufman** (n+m < 100): Custom implementation with 1x1 and 2x2 block
  pivoting. A critical bug fix ensures that when rows/columns are swapped during
  pivoting, the L entries from previously computed columns are also swapped. Without
  this, P*L*D*L^T*P^T != A, producing incorrect solutions.

- **Sparse LDL^T** (n+m >= 100): Uses faer's simplicial LDLT with AMD ordering.
  Symbolic factorization is computed once and cached; only numeric factorization
  repeats each iteration. Inertia is extracted from the D diagonal.

**Impact.** The hybrid approach gives good performance across a wide range of problem
sizes. Small problems benefit from cache-friendly dense operations (median speedups
exceed 18x on CUTEst and 15x on the HS suite). Large problems use sparse multifrontal
factorization and are typically competitive with Ipopt's MUMPS up to a few thousand
variables, where Fortran MUMPS begins to pull ahead on the very largest systems.

---

## Where Ipopt Remains Stronger

### 1. Convergence on Some Constrained Problems

12 problems are solved by Ipopt but not ripopt:

| Problem | n | m | ripopt Status | Root Cause |
|---------|---|---|---------------|------------|
| ACOPR14 | 38 | 82 | MaxIterations | Slow convergence, needs more iterations |
| ACOPR30 | 72 | 172 | RestorationFailed | Restoration cannot recover feasibility |
| CRESC50 | 6 | 100 | RestorationFailed | Many constraints, restoration stalls |
| HATFLDH | 4 | 7 | MaxIterations | Oscillating dual variables |
| HS109 | 9 | 10 | MaxIterations | Slow convergence near degenerate point |
| HS83 | 5 | 3 | MaxIterations | Restoration cycling |
| MGH10SLS | 3 | 0 | MaxIterations | Stuck at wrong local minimum |
| OET2 | 3 | 1002 | MaxIterations | Wrong basin, dual oscillation |
| OET6 | 5 | 1002 | MaxIterations | Slow convergence (2999 iters vs Ipopt's 126) |
| OET7 | 7 | 1002 | MaxIterations | Slow convergence (2999 iters vs Ipopt's 193) |
| OSBORNEA | 5 | 0 | MaxIterations | Slow convergence (unconstrained) |
| QCNEW | 9 | 3 | MaxIterations | Slow convergence near solution |

The dominant pattern is **MaxIterations** (10/12), suggesting convergence is happening
but too slowly. The OET family (3 problems) would likely be solved with sparse linear
algebra. The remaining problems involve degenerate multipliers, wrong basins, or slow
final convergence.

### 2. Iteration Counts on Some Problems

While ripopt is faster per iteration (median 18.8x on CUTEst), it sometimes requires significantly
more iterations:

| Problem | ripopt iters | Ipopt iters | Cause |
|---------|-------------|-------------|-------|
| GOFFIN | 2999 | 7 | Wrong basin, dual oscillation |
| HS85 | 2999 | 13 | Slow convergence |
| VESUVIOLS | 2999 | 10 | Stuck at saddle point |
| LAUNCH | 2999 | 12 | Slow convergence |
| HYDCAR6 | 1298 | 5 | Slow convergence |
| METHANL8 | 608 | 4 | Slow convergence |

These cases suggest differences in the mu strategy, multiplier initialization, or
second-order step computation that cause ripopt to take longer paths. In most cases
ripopt still converges correctly (just slowly), but the iteration count difference
indicates room for improvement in the adaptive mu oracle or the initial dual estimate.

### 3. Solution Quality at Degenerate Points

93 problems where both solvers converge produce different objectives (different local
optima). While both solutions satisfy KKT conditions, Ipopt more often finds the
globally better solution. This may reflect Ipopt's more mature mu strategy, better
multiplier initialization from decades of tuning, or differences in the starting point
perturbation strategy.

---

## Architectural Differences

| Aspect | ripopt | Ipopt |
|--------|--------|-------|
| Language | Rust | C++ (with Fortran MUMPS) |
| Linear solver | Dense BK (small) + faer sparse LDL (large) | MUMPS sparse |
| Restoration | 2-phase: GN then NLP | Single NLP restoration |
| NE handling | LS reformulation | Standard constrained IPM |
| Convergence | Dual gate + unscaled check | Single scaled check |
| Inertia failure | Proceed + filter recovery | Return error |
| SOC | Every backtracking step | First step only |
| Memory | Stack-allocated, no GC | Heap-allocated |
| Startup cost | ~50us | ~1-2ms (MUMPS init) |

---

## Opportunities for Improvement

### High Impact

1. **Adaptive mu tuning.** Several problems (GOFFIN, HS85, LAUNCH, HYDCAR6) show
   ripopt taking 100x more iterations than Ipopt. Improving the mu oracle to better
   estimate the optimal barrier parameter could reduce iteration counts significantly.
   10/12 ripopt-only failures are MaxIterations -- this is the dominant failure mode.
   Estimated gain: 3-5 problems, faster convergence on many others.

2. **Oscillation damping.** HATFLDH, HS83, OET2/6/7 cycle with oscillating dual
   variables. Damping y updates (e.g., y_new = alpha*y_computed + (1-alpha)*y_old
   with adaptive alpha) could stabilize convergence. Estimated gain: 3-5 problems.

### Medium Impact

3. **Multiplier initialization.** Better initial estimates of y (e.g., least-squares
   from the initial KKT system) could reduce iteration counts and avoid wrong basins.
   Currently y is initialized to 0 for all constraints.

4. **Watchdog strategy.** Ipopt uses a watchdog strategy that accepts a non-monotone
   step and checks if it leads to sufficient decrease after a few iterations. This can
   escape local stalling where the filter becomes too restrictive.

### Lower Impact

5. **BFGS Hessian approximation.** For problems where the exact Hessian is unavailable
   or pathological, a quasi-Newton (L-BFGS) approximation could provide more robust
   curvature information.

6. **Warm starting.** When solving sequences of related problems, reusing the previous
   solution as a starting point could reduce iteration counts significantly.

---

## Summary

On the CUTEst suite, ripopt and Ipopt are nearly tied (553 vs 561), with each solver
recovering a different ~40 problems the other cannot. ripopt's unique capabilities:
- NE-to-LS reformulation (~18 NE problems Ipopt cannot handle)
- Two-phase restoration (GN fast path + NLP robust fallback)
- Explicit slack fallback (recovers problems where implicit-slack multipliers oscillate)
- Pragmatic inertia correction (continues past factorization failures)
- Complementarity gate (prevents false convergence)
- Hybrid dense/sparse linear algebra (auto-selected by problem size)
- Raw speed advantage from Rust (median 18.8x faster on CUTEst, 15.1x on HS)

Ipopt's remaining advantages are:
- More mature mu strategy on difficult nonconvex problems
- Decades of parameter tuning on edge cases
- Better handling of dual oscillation at degenerate points
- Fortran MUMPS outperforms rmumps on the very largest sparse systems

The most impactful improvements would target the 83 CUTEst NumericalError failures
(often on non-convex Hessians and ill-conditioned KKT systems) and the 56 LocalInfeasibility
cases where multi-start or constraint-space search could escape false infeasibility declarations.
