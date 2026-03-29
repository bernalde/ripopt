# Algorithm

ripopt implements the primal-dual interior point method (IPM) described in [Wächter & Biegler (2006)](https://link.springer.com/article/10.1007/s10107-004-0559-y), with several extensions for robustness and performance.

## Problem form

```
min   f(x)
s.t.  g_l ≤ g(x) ≤ g_u
      x_l ≤ x    ≤ x_u
```

## Barrier formulation

Inequality constraints and variable bounds are handled via a logarithmic barrier:

```
φ(x, μ) = f(x) − μ Σ log(xᵢ − xₗᵢ) − μ Σ log(xᵤᵢ − xᵢ) − μ Σ log(gᵤⱼ − gⱼ) − μ Σ log(gⱼ − gₗⱼ)
```

As μ → 0, the solution of min φ(x, μ) converges to the solution of the original NLP.

## Perturbed KKT system

Stationarity of the barrier problem gives the perturbed KKT conditions:

- **Stationarity:** `∇f(x) + J(x)ᵀy − z_l + z_u = 0`
- **Feasibility:** `g_l ≤ g(x) ≤ g_u`
- **Complementarity:** `z_lᵢ · (xᵢ − xₗᵢ) = μ` and `z_uᵢ · (xᵤᵢ − xᵢ) = μ`

where `y` are constraint multipliers and `z_l, z_u` are bound multipliers.

## Newton step

Each iteration solves the augmented KKT system:

```
[H + Σ + δ_w I,  Jᵀ      ] [dx]   [r_d]
[J,              −δ_c I  ] [dy] = [r_p]
```

- **H** = ∇²_xx L (Lagrangian Hessian, lower triangle)
- **Σ** = diag(z_l/s_l + z_u/s_u) (barrier diagonal)
- **δ_w, δ_c** = inertia correction perturbations
- **r_d** = −(∇f + Jᵀy − z_l + z_u) (dual residual)
- **r_p** = −g(x) + slack (primal residual)

### Linear solvers

| Problem size | Solver | Notes |
|---|---|---|
| n + m < 110 | Dense Bunch-Kaufman LDL^T | Stack allocation, exact inertia |
| m ≥ 2n, n ≤ 100 | Dense condensed KKT (Schur complement) | 100-800x speedup for m >> n |
| n + m ≥ 110 | rmumps multifrontal LDL^T | SuiteSparse AMD reordering |

### Inertia correction

After factorization, the inertia (sign counts of D in LDL^T) is checked. The correct inertia is (n positive, m negative, 0 zero). If wrong, `δ_w` and `δ_c` are increased and the matrix is re-factored.

## Mehrotra predictor-corrector

Instead of a fixed centering parameter σ = 0.1, ripopt uses Mehrotra's adaptive rule:

1. **Predictor step** (σ = 0): compute affine-scaling direction `dx_aff`
2. **Adaptive centering**: `σ = (μ_aff/μ)³` where `μ_aff` is the complementarity after the affine step
3. **Corrector step**: solve the same KKT system with a modified RHS using σ·μ and affine cross-terms

This reuses the same LDL^T factorization (one extra triangular solve) and typically reduces iteration counts by 20–40%. **Gondzio centrality corrections** (up to 3 per iteration) further drive outlier complementarity pairs back to the central path.

## Filter line search

The step size α is chosen by the filter method (Fletcher & Leyffer, 2002):

- A **filter** maintains pairs (θ, φ) where θ = constraint violation and φ = barrier objective
- A trial point is **acceptable** if it is not dominated by any filter entry: θ_trial < (1−γ_θ)θ or φ_trial < φ − γ_φ·θ
- **Switching condition**: when sufficiently feasible (θ < θ_min), switch from filter to Armijo criterion
- **Second-order corrections (SOC)**: when a step is rejected, correct the constraint residual to account for nonlinearity (up to 4 SOC per iteration)

## Step size rules

- **Fraction-to-boundary**: α ≤ τ · max{α : x + α·dx > x_l, z + α·dz > 0} with τ = 0.99
- Separate primal and dual step sizes, with Ipopt's `alpha_y = alpha_d` convention

## Barrier parameter update

**Free mode** (default): oracle-based μ selection from complementarity, with filter reset each iteration.

**Fixed mode**: monotone decrease μ_new = factor · μ until barrier subproblem converges to `barrier_tol_factor · μ`. Triggered automatically when the free-mode oracle stalls.

## Restoration phase

When the filter line search fails (all backtracking steps rejected), restoration recovers feasibility:

1. **Gauss-Newton restoration** (fast): minimize `(1/2)||g(x)||²` using GN steps. Quadratic convergence for nonlinear equalities; falls back to gradient descent when Jacobian is rank-deficient.

2. **NLP restoration subproblem** (robust, triggered after 5 GN failures): solve an auxiliary NLP:
   ```
   min ρ·(Σpᵢ + Σnᵢ) + (η/2)·||D_R(x − x_r)||²
   s.t. g(x) − p + n = g_target
   ```
   This minimizes constraint violation with a proximity term to the reference point x_r.

## Fallback cascade

When the primary IPM fails, ripopt automatically tries:

1. **L-BFGS Hessian approximation**: replaces exact Hessian with L-BFGS curvature pairs (no second derivatives needed)
2. **Augmented Lagrangian**: for equality-only problems; L-BFGS inner solver with multiplier updates
3. **SQP**: sequential quadratic programming for small constrained problems
4. **Slack reformulation**: converts inequality constraints to equality + bounds on slacks (`g(x) − s = 0`)

Each fallback uses the same interface and inherits the time budget from the primary solve.

## Convergence criteria

Three scaled KKT residuals must be below tolerance:

| Residual | Formula | Default tol |
|---|---|---|
| Primal infeasibility | `max |gᵢ(x) − proj(gᵢ, [gₗᵢ, gᵤᵢ])|` | 1e-4 |
| Dual infeasibility | `||∇f + Jᵀy − z_l + z_u|| / (1 + ||y,z||)` | 1e-8 |
| Complementarity | `max |zₗᵢ·sₗᵢ − μ|, |zᵤᵢ·sᵤᵢ − μ| / (1 + ||z||)` | 1e-4 |

## Key papers

- Wächter & Biegler (2006). *On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming.* Mathematical Programming.
- Mehrotra (1992). *On the implementation of a primal-dual interior point method.* SIAM J. Optimization.
- Gondzio (1996). *Multiple centrality corrections in a primal-dual method for linear programming.* Computational Optimization and Applications.
- Fletcher & Leyffer (2002). *Nonlinear programming without a penalty function.* Mathematical Programming.
