# Reference Gap Roadmap

A citation-grounded comparison of ripopt and rmumps against the canonical
solvers they were derived from (Ipopt 3.14.x C++ and MUMPS 5.8.2 Fortran),
written as a development roadmap.

Every claim below points at `file:line` in either the ripopt tree or the
vendored reference trees at `ref/Ipopt/` and `ref/mumps/`. The two sections
were drafted by specialist agents reading the actual source of both sides;
the roadmap summary at the top is distilled from those sections and is
where a reader short on time should start.

---

## Roadmap summary (ranked, cross-cutting)

### Correctness-first (may produce wrong answers today)

1. ~~**ripopt: separate `s_d` / `s_c` scaling denominators and remove the
   1e4 cap.**~~ **DONE (commit 698a79f).** `src/convergence.rs` now
   computes separate denominators with `s_max=100` and no upper cap,
   matching `IpIpoptCalculatedQuantities.cpp:3663-3700`.
2. **ripopt: post-hoc unscaled complementarity gate (HS13 root cause).**
   *Re-scoped*: the expert pass on `IpIpoptAlg.cpp:1055-1135` showed
   that Ipopt's actual mechanism is `correct_bound_multiplier` (a
   `kappa_sigma=1e10` projection on `z` after each accepted step), not
   a post-hoc recompute. **Already implemented** in
   `apply_kappa_sigma_bound_multiplier_reset` at `src/ipm.rs:2665-2696`,
   called unconditionally from `update_dual_variables`. The convergence
   gate at `src/convergence.rs:73` therefore reads a `z` that has
   already been clamped, which is what Ipopt's gate sees too. No
   additional code change needed for this item.
3. **rmumps: wire up CNTL(4)-equivalent static-pivot threshold.**
   `frontal.rs:745` hardcodes `seuil = 0.0`, so the "must-eliminate"
   branch at `frontal.rs:864-881` accepts a tiny pivot at its literal
   value. A 1e-18 pivot produces a factorization whose inertia is noise;
   the IPM then steps along a meaningless direction. Biggest
   correctness risk for the entire pipeline.
4. **rmumps: add ICNTL(24)-style null-pivot detection.** Zero columns
   are currently accepted silently (`pivot.rs:217-219`,
   `frontal.rs:772-773`). MUMPS fixes them to `±CNTL(5) * ||A||` and
   reports the count via `INFOG(28)`. Without this, rank-deficient
   iterates are invisible to the IPM.
5. ~~**ripopt: drop the `±1` inertia acceptance heuristic.**~~
   **DONE (commit 53f4275).** All `approx_ok` branches in
   `factor_with_inertia_correction` removed; only exact inertia is
   accepted now.
6. ~~**ripopt: mu-dependent `delta_c` regularization.**~~
   **DONE (commit 83fcbc0).** `assemble_kkt` and
   `factor_with_inertia_correction` now scale `delta_c_base` by
   `mu^0.25` (matching `IpPDPerturbationHandler.cpp:82-94`).
   `InertiaCorrectionParams::default().delta_c_base` is now `1e-8`.
7. **rmumps: track tiny/static pivots and growth factor.** Without an
   `NBTINYW`/`RINFOG` equivalent, the IPM cannot tell a clean
   factorization from one that papered over 40 tiny pivots.
8. **ripopt: preprocessing redundancy detection false-positive.**
   `src/preprocessing.rs:216-334` uses a two-point probe; can drop a
   genuinely-independent constraint when the two points happen to
   coincide.
9. **ripopt: `user_x_scaling` option is declared but not applied**
   (`src/options.rs:188` vs `src/ipm.rs:2094-2170`). Silent no-op.
10. **rmumps: unify the two pivot-search paths.** The classic-BK code
    path runs when `pivot_threshold = 0.0` and cannot emit delayed
    pivots at all. Callers who disable thresholding think they are
    getting a more aggressive BK factorization; they are silently
    losing the entire delayed-pivot machinery.

### Convergence and robustness

11. **rmumps: real MC64 scaling.** The current `compute_mc64_kkt_scaling`
    is a bespoke heuristic; a true MC64 would let the pivot threshold
    drop from 0.01 toward Ipopt's 1e-6 without losing stability.
12. **ripopt: soft-restoration phase.** When filter rejects, Ipopt first
    tries a primal-dual soft restoration before restarting a nested
    IPM. ripopt jumps directly to full restoration or failure.
13. **ripopt: port `RestoFilterConvCheck`.** Restoration currently exits
    on feasibility alone; Ipopt additionally requires filter-acceptance
    in the original problem. Ping-pong between main IPM and restoration
    is a real failure mode for several CUTEst problems.
14. **ripopt: quality-function mu oracle must include centrality term.**
    Production Free-mode oracle is `compute_loqo_mu` at `src/ipm.rs`,
    which already incorporates centrality via the Loqo σ formula
    `0.1·min(0.05·(1-ξ)/ξ, 2)³`. The standalone reference QF oracle
    `quality_function_mu` (still `#[allow(dead_code)]`) now matches
    Ipopt's `IpQualityFunctionMuOracle.cpp:622-646` formula structure
    (1-norm averages summed, plus `compl_inf / xi` when
    `quality_function_centrality=true`; default false matches Ipopt's
    `centrality=none`). Wiring it as a selectable production oracle
    is deferred — Loqo with embedded centrality already covers the
    documented "drops mu too aggressively off-center" failure mode.
15. **rmumps: port analysis-time front sizing.** `expand_for_delayed`
    (`frontal.rs:161-183`) reallocates and copies dense fronts at
    runtime; MUMPS sizes fronts symbolically to include delayed slots.
    Predictability at scale.
16. **rmumps: backward-error-based refinement stop.** Current stop at
    `solver.rs:218` is a geometric-decrease heuristic; MUMPS uses
    Arioli–Demmel–Duff omega1+omega2. CONCON-style stalls probably
    benefit.

### Ecosystem / parity gaps

17. **rmumps: METIS fallback via optional feature flag.** AMD on a 50k
    KKT is measurably worse than METIS for both fill and flop count.
    `ana_set_ordering.F` is the reference auto-fallback.
18. **ripopt: `warm_start_target_mu` + per-component bound_push.**
    Parametric / MPC loops cannot resume at a user-specified mu today.
19. **ripopt: full `PDPerturbationHandler` port including `delta_s`
    and `delta_d`.** Required for correct regularization trajectory
    on inequality-only rank deficiency.
20. **rmumps: static 2×2 pivots.** `frontal.rs:864-881` only does 1×1
    static fallback; MUMPS keeps trying 2×2 under `STATICMODE`.

Items 1–10 are the correctness bar for a v1.0. Items 11–16 are the
delta from "solves" to "solves competitively with Ipopt on large and
degenerate problems." Items 17–20 are parity polish.

---

## ripopt vs Ipopt (interior-point core)

### 1. Algorithmic parity audit

The following table is grounded in a direct read of both sources.
"Faithful" means the math, constants, and escalation logic match Ipopt
within reasonable tolerance; "partial" means the mechanism exists but a
sub-component is missing or simplified; "absent" means ripopt has no
equivalent.

| Ipopt mechanism | ripopt status | Where (ripopt) | Where (Ipopt) |
|---|---|---|---|
| Full-space KKT with dense/sparse LDL^T | Faithful | `src/kkt.rs:52-261` (assemble), `src/kkt.rs:690-842` (solve + iter-refinement) | `IpPDFullSpaceSolver.cpp` |
| Sigma = Z / (x - x_L) + Z / (x_U - x), bound-mult elimination | Faithful | `src/kkt.rs:263-295` | `IpPDFullSpaceSolver.cpp` |
| `PDPerturbationHandler`: delta_w, delta_c, delta_d escalation with history | **Partial / simplified** | `src/kkt.rs:333-651` (`factor_with_inertia_correction` + `InertiaCorrectionParams`) | `IpPDPerturbationHandler.cpp:17-356` |
| delta_s (separate slack perturbation) | **Absent** — ripopt lumps slacks into x | n/a | `IpPDPerturbationHandler.cpp:130-181` |
| `perturb_dec_fact` (kappa_w^-) warm-shrink on success | Partial (`delta_w_last / delta_w_growth`, `src/kkt.rs:571`) | | `IpPDPerturbationHandler.cpp:73` |
| `perturb_inc_fact_first` distinct from `perturb_inc_fact` | Absent (single growth factor 4.0, `src/kkt.rs:358`) | | `IpPDPerturbationHandler.cpp:50-65` |
| delta_c tied to mu^kappa_c (`jacobian_regularization_exponent=0.25`) | Absent (`delta_c_base = 1e-4`, fixed) | `src/kkt.rs:357` | `IpPDPerturbationHandler.cpp:82-94` |
| Filter: (theta, phi) entries with gamma_theta/gamma_phi envelope | Faithful | `src/filter.rs:64-79,171-178` | `IpFilter.cpp`, `IpFilterLSAcceptor.cpp:311-437` |
| Switching condition, f/h-type classification | Faithful | `src/filter.rs:86-99,133-168` | `IpFilterLSAcceptor.cpp` |
| Armijo on phi | Faithful | `src/filter.rs:102-110` | `IpFilterLSAcceptor.cpp` |
| Second-order correction (SOC) | Faithful (full + condensed + sparse variants, up to `max_soc`) | `src/ipm.rs:5554-5995` | `IpFilterLSAcceptor.cpp`, `IpBacktrackingLineSearch.cpp` |
| Watchdog (max `watchdog_trial_iter_max` steps, rollback to `watchdog_saved`) | Faithful | `src/ipm.rs:4895-4977` | `IpBacktrackingLineSearch.cpp:376-...` |
| Soft-restoration phase (`SoftRestoLSAcceptor`, pderror_reduction) | **Absent** | n/a | `IpBacktrackingLineSearch.cpp:173-220` |
| Fraction-to-boundary with tau_min | Faithful | `src/ipm.rs:4217,4381,4408,5576` (tau = max(1-mu, tau_min)) | `IpFilterLSAcceptor.cpp` |
| Alpha_min floor (filter-problem-dependent) | Faithful | `src/filter.rs:186-204` | `IpFilterLSAcceptor.cpp:450-469` |
| Filter augmentation at restoration entry | Faithful | `src/filter.rs:211-216` | `IpFilterLSAcceptor.cpp:898-901` |
| Monotone mu (Fiacco-McCormick, kappa_mu, theta_mu) | Partial (barrier-subproblem stop test exists, `src/ipm.rs:5187`) | | `IpMonotoneMuUpdate.cpp:135-194` |
| Quality-function mu oracle | Partial (grid search over log(mu), fixed-part + compl only) | `src/ipm.rs:6331-6366` | `IpQualityFunctionMuOracle.cpp` |
| LOQO mu oracle | Absent | n/a | `IpLoqoMuOracle.cpp` |
| Probing (Mehrotra affine) mu oracle | Partial (affine predictor + Mehrotra corrector exists at `src/ipm.rs:3961-4340`, not a full oracle) | | `IpProbingMuOracle.cpp` |
| Adaptive / free-fixed mode switching globalization | Partial (`adaptive_mu_monotone_init_factor` hook) | `src/ipm.rs:4829,5313` | `IpAdaptiveMuUpdate.cpp` |
| Restoration NLP: rho*(sum p + sum n) + (eta/2)*||D_R(x-x_r)||^2 | Faithful math, DR_i = 1/max(1, |x_r[i]|) | `src/restoration_nlp.rs:7-208` | `IpRestoIpoptNLP.cpp:470-478,763` |
| eta = resto_proximity_weight * mu^kappa_eta (kappa_eta=0.5) | Faithful (`eta = eta_f * sqrt(mu_entry)`) | `src/restoration_nlp.rs:56` | `IpRestoIpoptNLP.cpp:34,763` |
| p/n closed-form init from `p_i * n_i = mu`, `p_i - n_i = c_i` | Faithful | `src/restoration_nlp.rs:108-138` | `IpRestoIterateInitializer.cpp:79-97` |
| Restoration convergence check (`RestoFilterConvCheck`) | **Partial / ad hoc** | `src/restoration.rs:42-...`, post-hoc hand-off in `src/ipm.rs:5996-6157` | `IpRestoFilterConvCheck.cpp` |
| User scaling (`user_obj_scaling`, `user_g_scaling`, `user_x_scaling`) | obj and g supported; **x_scaling not applied in IPM** | `src/ipm.rs:2105-2108` | `IpNLPScaling.cpp`, `IpUserScaling.cpp` |
| Gradient-based scaling (`nlp_scaling_max_gradient=100`, `min_value=1e-2`) | Faithful constants | `src/ipm.rs:2111-2158` | `IpGradientScaling.cpp` |
| Equilibration / MC19 scaling | Absent (only Ruiz on the KKT matrix, `src/kkt.rs:297-351`) | | `IpEquilibrationScaling.cpp` |
| Scaled-KKT error via s_d, s_c denominators | **Divergent formula** | `src/convergence.rs:51-64` (single `s_d`, no `s_c`) | `IpIpoptCalculatedQuantities.cpp:3078-3099, 3663-3700` |
| `dual_inf_tol`, `compl_inf_tol`, `constr_viol_tol` gates | Wired, but see §2 | `src/convergence.rs:68-76` | `IpOptErrorConvCheck.cpp:149-227` |
| TNLPAdapter: `fixed_variable_treatment` (make_constraint / make_parameter / relax_bounds) | **Only "make_parameter" equivalent** | `src/preprocessing.rs:56-74` | `IpTNLPAdapter.cpp:101-108,240-501` |
| `DetermineDependentConstraints` (MA28/Mumps-based dep detector) | Replaced with pattern+value-match heuristic at two probe points | `src/preprocessing.rs:216-334` | `IpTNLPAdapter.cpp:3396-3566` |
| One-sided g(x) ≤ u / l ≤ g(x) kept explicit (no artificial slack doubling) | Handled by ripopt's single slack-less formulation | `src/ipm.rs`/`src/convergence.rs` | `IpTNLPAdapter.cpp` |
| Warm start: `warm_start_init_point`, `warm_start_bound_push`, `warm_start_mult_bound_push`, `warm_start_bound_frac` | Faithful for primal and bound multipliers | `src/warmstart.rs:25-100` | `IpWarmStartIterateInitializer.cpp` |
| Warm start of y (constraint multipliers) | Present (passes through y from user) | `src/warmstart.rs` | |
| Warm start target-mu (`warm_start_target_mu`, `mu_init` override) | Absent | | `IpWarmStartIterateInitializer.cpp` |
| sIPOPT parametric sensitivity | One mode only (one-shot KKT back-solve at optimum) | `src/sensitivity.rs:67-285` | `contrib/sIPOPT/src/SensAlgorithm.cpp` |
| sIPOPT approximation modes (1=sens_only, 2=sens+update_bounds, 3=...) | Not exposed; ripopt does sens_only | | `contrib/sIPOPT/src/SensAlgorithm.hpp` |

### 2. Known deficiencies in ripopt (correctness risks)

Ranked by the severity of the failure mode — wrong-answer first, then
silently-weakened tolerance, then missing feature.

**D1. `s_d` / `s_c` conflation and arbitrary cap (WRONG-ANSWER risk).**
`src/convergence.rs:51-64` computes a single scaling `s_d` from the
average over *all* multipliers (`info.multiplier_sum / info.multiplier_count`)
and uses the same `s_d` to scale both the dual tolerance and the
complementarity tolerance. Ipopt computes two separate denominators at
`IpIpoptCalculatedQuantities.cpp:3677-3699`:

- `s_c = (||z_L||_1 + ||z_U||_1 + ||v_L||_1 + ||v_U||_1) / n_bounds`,
  then `s_c = max(s_max, s_c) / s_max` — **bound multipliers only**;
- `s_d = (||y_c||_1 + ||y_d||_1 + ||z_L||_1 + ||z_U||_1 + ||v_L||_1 + ||v_U||_1) / (m+n_bounds)`,
  then `s_d = max(s_max, s_d) / s_max`.

Additionally ripopt caps `s_d` at `1e4` (`src/convergence.rs:55`); Ipopt
has no such cap. Consequence: on a problem where `y` is O(1) but `z_L`
blows up (or vice-versa), ripopt's single denominator is driven by the
large group and the tolerance on the small group is loosened by the same
factor. A point where `||grad L||_∞ ≈ 1e-4`, `y` large, but `compl ≈ 1e-4`
can pass ripopt's scaled check while Ipopt would reject compl because
`s_c = 1`. The 1e4 cap means ripopt *also* loosens dual_tol by up to
10000x silently — a point with raw `||grad L||_∞ = 1e-4` can be declared
converged when `tol=1e-8`. The unscaled gate at `src/convergence.rs:71-73`
partially guards this only because `dual_inf_tol=1.0` by default, but any
user who tightens `dual_inf_tol` or `compl_inf_tol` below 1e-4 loses the
safety net.

**D2. Unscaled gate uses the scaled iterative `compl_inf` (HS13-style
termination).** `src/convergence.rs:73` tests `info.compl_inf <= options.compl_inf_tol`
where `compl_inf` is the same iterative-z complementarity used on the
scaled side (and itself is evaluated against `mu`, not the target 0 —
see `complementarity_error` at `src/convergence.rs:235-255`). Ipopt's
unscaled gate uses `unscaled_curr_complementarity(mu_target_, NORM_MAX)`
at `IpOptErrorConvCheck.cpp:211`. For HS13-style problems where the
solver quits at a KKT-infeasible point (x near a spurious stationary
point), iterative z is recomputed at the quit iterate and made to look
near-mu by construction. Ipopt's unscaled check recomputes `z_L * (x - x_L)`
against the current iterate without least-squares patching and so fails.
Root cause: ripopt's `z_l`/`z_u` are evolved by the IPM (not re-derived
post-hoc), so by the time ripopt reaches an HS13-style quit,
`(x-x_l) * z_l - mu` may be O(mu) by construction even when the true KKT
residual is large. **Fix:** add a post-hoc unscaled `z` (clip and
re-derive from x alone, or use primal stationarity residual without
z_l, z_u) for the final gate.

**D3. Missing delta_s perturbation (incorrect inertia for implicit-slack
formulation).** Ipopt's augmented system has four perturbation scalars
(delta_x, delta_s, delta_c, delta_d). ripopt uses *two* (`src/kkt.rs:335-360, 414-651`).
Because ripopt folds slacks into x via `SlackFormulation`
(`src/slack_formulation.rs`), the missing delta_s is partially absorbed —
but the distinction between delta_c (equality) and delta_d (inequality,
gets the Sigma_s contribution) matters when the constraint Jacobian is
rank-deficient only on the inequality subset. On such problems (common
in network flow with active bounds), ripopt perturbs either all
constraints or none, and the inertia-correction loop can overshoot
delta_c into the equality block when only the d block was singular.

**D4. Inertia acceptance with `±1` tolerance at n+m ≥ 100 (WRONG-ANSWER
risk for degenerate problems).** `src/kkt.rs:438-444` accepts the
factorization when inertia counts are off by ±1 at ≥100 rows. Ipopt
never accepts off-by-one inertia; `IpPDPerturbationHandler.cpp`
re-regularizes until exact match or skips to restoration. A system with
one zero eigenvalue can be accepted as "one positive, n-1 negatives"
(wrong direction), producing a step that climbs rather than descends.
Filter line search usually rejects it, but not on the very first
iteration where the filter is effectively empty.

**D5. Constant `delta_c_base` (1e-4 or 1e-8, `src/kkt.rs:228,357`)
versus mu-dependent delta_c.** Ipopt's delta_c scales with mu via
`delta_cd = jacobian_regularization_value * mu^jacobian_regularization_exponent`
(`IpPDPerturbationHandler.cpp:82-94`, default `1e-8 * mu^0.25`).
ripopt's fixed 1e-4 becomes overwhelming as mu → 1e-10: a 1e-4
constraint-block regularization perturbs the converged system far more
than Ipopt does and biases the final multipliers `y` by
`O(delta_c / ||J||)`. Users targeting `tol=1e-10` will see the final
`dual_inf_unscaled` stuck near `delta_c`.

**D6. Redundant-constraint detection can drop a genuinely-independent
constraint.** `src/preprocessing.rs:216-334` decides two constraints
are redundant if they have the same Jacobian pattern at two probe
points AND same g(x) value AND same bounds. On quadratic constraints
that happen to coincide at the two probes, the probe test fails to
separate them. Ipopt uses MA28/Mumps-based dependency detection on the
full symbolic Jacobian, which does not have this false-positive mode.
Preprocessing is enabled by default and is invisible.

**D7. `user_x_scaling` accepted but not applied.**
`src/options.rs:188` declares the option but `src/ipm.rs:2094-2170`
only uses `user_obj_scaling` and `user_g_scaling`. Users supplying
x-scaling silently get no scaling.

**D8. No `soft_resto` intermediate phase.** When a step fails the
filter, ripopt jumps directly to full restoration NLP
(`src/ipm.rs:6158-6330`) or declares failure. Ipopt's
`BacktrackingLineSearch.cpp:376-...` first tries a soft-restoration
mode that reuses the primal-dual system with a pderror-reduction
criterion before restarting a nested IPM.

**D9. Quality-function oracle ignores the centrality measure.**
`quality_function_mu` at `src/ipm.rs:6331-6366` minimizes
`(primal_inf)^2 + (dual_inf)^2 + (compl_inf(mu_candidate))^2`. Ipopt
also includes the centrality term `||X*Z*e - mu*e|| / mu`. Effect:
ripopt's adaptive mu can pick small mu too aggressively when the
iterate is off-center, producing overly-aggressive steps that get
rejected.

**D10. `perturb_dec_fact` behavior.** Ipopt defaults `kappa_w^- = 1/3`
(geometric shrink) with a floor at `delta_xs_min = 1e-20`. ripopt
divides `delta_w_last / delta_w_growth` (= /4.0) and floors at
`delta_w_init` (1e-4). The floor is ~16 orders of magnitude larger.
Once ripopt has ever needed regularization, it cannot return to the
unperturbed regime for the remainder of the solve — causing the
late-mu bias noted in D5.

### 3. Genuine advantages of ripopt

Not "faster because Rust" — these are structural/algorithmic
differences where ripopt's choice is defensible or superior.

**A1. Structural degeneracy escalation ladder (`src/kkt.rs:429-491,565`).**
ripopt maintains `degeneracy_count` and a `structurally_degenerate`
flag so that after 3 consecutive iterations that needed perturbation,
it *skips* the wasted unperturbed factorization. Ipopt always tries
the unregularized factor first.

**A2. Augmented LS-y at `src/ipm.rs:6388-...` instead of normal-equations.**
Ipopt's least-square multiplier estimator
(`IpLeastSquareMults.cpp:82-87`) solves the saddle-point system but
falls back to normal equations on failure. ripopt *always* uses the
saddle-point form with Bunch-Kaufman, avoiding the `J * J^T`
singularity mode on gauge-symmetric problems.

**A3. NE-problem detection and LS reformulation (`src/ipm.rs:1077-...`).**
`detect_ne_problem` identifies nonlinear equation systems (m = n, no
objective, no bounds) and dispatches to a least-squares formulation.
Ipopt has no such heuristic and stalls because `grad f ≡ 0` makes the
Armijo condition vacuous.

**A4. Explicit backward-error probe on accepted factorizations
(`src/kkt.rs:450,516,604`).** ripopt does a test solve
(`check_factorization_backward_error`) even when inertia looks
correct. Ipopt trusts the inertia counts; on AC-OPF gauge cases the
inertia can be (n, m, 0) while the system is rank-deficient.

**A5. Closed-form p/n initialization in restoration NLP that enforces
complementarity at iteration 0 (`src/restoration_nlp.rs:108-138`).**
Parity with Ipopt (`IpRestoIterateInitializer.cpp:79-97`), but ripopt
additionally clamps to `max(safe_mu/rho)` and documents the
50-iteration savings explicitly.

**A6. Two-point linearity verification in preprocessing.**
`src/preprocessing.rs:101-135` evaluates the Jacobian at two points
before claiming a constraint is linear for bound tightening. Ipopt's
bound tightening is off by default and, when on, assumes the user
flagged linear constraints correctly.

**A7. Fallback solvers (`augmented_lagrangian.rs`, `sqp.rs`, `lbfgs.rs`).**
Ipopt has zero fallback; a failed restoration produces
`Restoration_Failed`. ripopt can attempt an augmented Lagrangian, SQP,
or L-BFGS pass on the failed problem. Out of scope for the IPM core
but materially improves HS and adversary benchmark pass rates.

### 4. Roadmap recommendations (ripopt-side)

Prioritized for HS/CUTEst pass-rate improvement and correctness-first.
Effort: S < 1 day, M 1-5 days, L > 1 week.

1. **Separate `s_d` and `s_c`, remove the 1e4 cap.** Port
   `ComputeOptimalityErrorScaling` from
   `IpIpoptCalculatedQuantities.cpp:3663-3700` verbatim. Fixes D1.
   **Effort: S.** Test: `check_convergence` with y large / z small and
   vice-versa.

2. **Fix the unscaled complementarity gate (HS13 root cause).** Add a
   post-hoc unscaled `z` recomputation in `check_convergence`, or
   separately evaluate `unscaled_curr_complementarity` using
   `mu_target = 0` rather than the current barrier mu.
   Where: `src/convergence.rs:71-73` and new helper equivalent to
   `IpIpoptCalculatedQuantities.cpp:1720-1873`. Fixes D2, HS13, likely
   several CUTEst `NumericalError` misclassifications. **Effort: M.**

3. **Port `PDPerturbationHandler` faithfully, including mu-dependent
   `delta_c = delta_cd_val * mu^kappa_c`.** Where: `src/kkt.rs:333-651`
   → replace with a port of `IpPDPerturbationHandler.cpp:17-356`. Add
   delta_s for the slack block, delta_d for the d-block. Addresses D3,
   D5, D10. **Effort: L.**

4. **Remove `±1` inertia acceptance; keep only the `delta_c-only`
   selective perturbation path.** Where: `src/kkt.rs:438-444,505-511,590-594`.
   Replace with exact match; when exact match fails, escalate delta_c
   first (the existing selective path at lines 524-563 is correct),
   then delta_w. Fixes D4. **Effort: S.**

5. **Implement `soft_resto` acceptor.** Port
   `IpBacktrackingLineSearch.cpp:173-220` and `IpSoftRestoLSAcceptor`.
   Fixes D8. **Effort: M.**

6. **Fix `user_x_scaling` application.** `src/options.rs:188` is wired
   in, but `src/ipm.rs:2094-2170` and the evaluations must apply it.
   Fixes D7. **Effort: S.**

7. **Guard preprocessing redundancy detection with a third probe point
   OR disable by default.** Where: `src/preprocessing.rs:216-334`.
   Fixes D6. **Effort: S.**

8. **Port the full quality-function oracle including centrality term.**
   Where: `src/ipm.rs:6331-6366`. Fixes D9. Expected outcome: ~10-20%
   iteration reduction on off-center problems. **Effort: M.**

9. **Port restoration convergence check (`RestoFilterConvCheck`).**
   Where: `src/restoration.rs:42-...`. Several CUTEst failures where
   ripopt re-enters main IPM with an inferior iterate that re-enters
   restoration (ping-pong) would be resolved. **Effort: M.**

10. **Expose `warm_start_target_mu` and per-component
    `warm_start_bound_push` paths.** Where: `src/warmstart.rs:25-100`.
    Critical for sequential solves (parametric NLP, MPC). **Effort: S.**

**Not recommended** (acceptable simplifications for ripopt's scope):
full MA28-style dependency detection (recommendation 7 is sufficient);
LOQO mu oracle (quality function + Mehrotra predictor already cover
the space); `fixed_variable_treatment = relax_bounds`; porting the
full sIPOPT Schur-complement driver (only needed for many-parameter
batch sensitivities, not a stated use case).

---

## rmumps vs MUMPS (sparse linear algebra)

The target configuration in both cases is KKT factorization: `SYM=2`
(general symmetric indefinite) with threshold partial pivoting.

### 1. Algorithmic parity audit

| Mechanism | MUMPS reference | rmumps status |
|---|---|---|
| AMD ordering | `ana_orderings.F` (`MUMPS_ANA_H`, `MUMPS_HAMF4`, `MUMPS_QAMD`) | Present via the external `amd` crate (SuiteSparse AMD) — `rmumps/src/ordering/amd.rs:60-73` |
| METIS / SCOTCH / PORD | `mumps_metis.c`, `mumps_scotch.c`, `PORD/`, selected by `ICNTL(7)` with automatic fallback in `ana_set_ordering.F` | Not present. `Ordering::NestedDissection` is an in-tree implementation at `rmumps/src/ordering/nested_dissection.rs`, *not* METIS. No SCOTCH, no PORD. |
| Maximum-transversal preprocessing (MC64) | `dana_mtrans.F` (1196 lines), activated by `ICNTL(6)` (default 7 = auto) | Absent in that form. rmumps has a *structural* bipartite matching in `ordering/kkt_matching.rs:24-53` that matches duals to primals by `|J_{ij}|`, then compresses the pair into AMD. Not MC64. |
| Supernodal / multifrontal | Multifrontal with elimination tree, fundamental + relaxed supernodes, assembly tree (`STEP`, `FILS`, `DAD_STEPS`) | Multifrontal with Liu-1990 etree (`etree.rs:15-52`), fundamental supernodes (`symbolic.rs:117-181`), relaxed amalgamation (`symbolic.rs:196-317`), NEMIN-style merging gated on `n < 10000` (`symbolic.rs:325-390`). Structurally equivalent. |
| Frontal assembly / extend-add | `dfac_asm.F`, `dfac_process_contrib_type1.F` — uses IW index lists off `IOLDPS+XSIZE` | `frontal.rs:149-235` (`extend_add`) with dynamic front expansion for delayed indices. The explicit `expand_for_delayed` path (`frontal.rs:185-235`) copies the whole dense matrix on each expansion — MUMPS instead sizes fronts symbolically to include delayed slots. |
| Threshold partial pivoting (CNTL(1)) | `UU` threshold applied in `DMUMPS_FAC_I_LDLT` at `dfac_front_aux.F:1231`; test at `:1781-1792` | `find_pivot_ldlt` at `frontal.rs:739-862`. Same first-acceptable pivot logic, same AMAX/RMAX separation, same 2×2 modified-Bunch-Kaufman test. Default `CNTL(1)=0.01` is correctly replicated at `pivot.rs:118`. |
| 2×2 pivots for indefinite | `DMUMPS_FAC_I_LDLT` computes `DETPIV`, checks both `(|d_jj|*RMAX + AMAX*TMAX)*UU <= |DETPIV|` and the symmetric test | Same formulas at `frontal.rs:824-830`. |
| Delayed pivots (failure → parent) | MUMPS sets `IW(IOLDPS+1+XSIZE)=NPIV`, decrements NASS, writes the failed column into the contribution block; parent receives it via the normal assembly path | Implemented at two levels. Leaf/type-1 path: columns that fail return `INOPV=2` and are left in the trailing block (`frontal.rs:864-881`). Parent side, `numeric.rs:273-298` detects child FS indices in the CB and calls `front.promote_cb_to_fs(&delayed_cols)` (`frontal.rs:69-117`). Same mechanism, re-implemented in user code rather than falling out of the symbolic structure. |
| Static pivoting / perturbation (CNTL(4), CNTL(5)) | `SEUIL = CNTL(4)`, `FIXA = DKEEP(2) = CNTL(5)*||A||`, tiny pivots replaced by `±CSEUIL`, counted as `NBTINYW`. `KEEP(97)` gates `STATICMODE` | Only partially present. rmumps has a "must-eliminate" path (`frontal.rs:452-456`, `:864-881`) that force-accepts a failed pivot for promoted-from-child columns and at the root supernode, but: (a) `seuil = 0.0` hard-coded at `frontal.rs:745` (no `CNTL(4)`), and (b) no `CNTL(5)`-style null-pivot fixation. A tiny pivot is not replaced by a bounded perturbation — it is silently accepted. |
| Scaling (ICNTL(8)) | `dfac_scalings.F` supports row/col equilibration, MC64-based scaling, ICNTL(8)=77 (auto) with SimScale schedule | Ruiz and diagonal equilibration at `scaling.rs:95-199`. SimScale schedule correctly replicated (`scaling.rs:126-144`). MC64-based scaling absent; a hand-rolled `compute_mc64_kkt_scaling` is referenced from `solver.rs:117` but only applied when `n_primal.is_some()`. |
| Iterative refinement (ICNTL(10)) | `dsol_driver.F:1908-5610` — Wilkinson-style refinement with adaptive stopping | Present in `solver.rs:186-242`. Refines in the *original* (unscaled, unpermuted) space, which matches MUMPS's default. Adaptive stop at `solver.rs:218` (`res_norm > 0.9 * prev_res_norm`). Default 10 steps. |
| Error analysis / condition estimation (ICNTL(11)) | `dsol_driver.F` computes `omega1`, `omega2`, forward and backward error estimates | Absent. No `ICNTL(11)` equivalent; no backward-error output to the caller. |
| Null-space detection (ICNTL(24)) | `DKEEP(1)=PIVNUL` threshold, null pivots added to `PIVNUL_LIST_STRUCT`; `INFOG(28)` returns count | Absent. Zero columns return `PivotResult::Delayed` (`pivot.rs:217-219`); at the root the "must-eliminate" branch accepts them as-is with no tracking. |
| Parallelism | MPI (type-2 nodes) + OpenMP (within-front BLAS, panel threading). Controlled by `KEEP(400)`, `ICNTL(16)`. | Sequential + rayon level-set parallelism over the assembly tree (`numeric.rs:140-165`), gated at `total_front_size > 4096`. No MPI, no shared-memory BLAS threading within a single front. Dense kernels optionally use `faer` (`frontal.rs:1087-1178`) or NEON SIMD (`pivot.rs:5-70`). |
| Memory management | `KEEP8`, `IS`, `S`, growth via `ICNTL(14)` relaxation; OOC; save/restore | Straight `Vec<f64>` allocations per front, no relaxation parameter, no OOC, no save/restore. Re-factor with new values (`solver.rs:274`) supported. |

### 2. Known deficiencies in rmumps (correctness risks)

Ordered by severity.

**(a) Static-pivoting value is zero, not CNTL(4).** `find_pivot_ldlt`
at `frontal.rs:745` hard-codes `let seuil: f64 = 0.0`, and the
"must-eliminate" fallback at `frontal.rs:864-881` accepts the failed
pivot at its current value without replacing it with `±max(|SEUIL|, eps)`
the way MUMPS does at `dfac_front_aux.F:1251-1258`. For a root
supernode or a promoted column whose current diagonal is `1e-18`,
rmumps will happily divide by it in `update_within_block_ldlt`
(`frontal.rs:923-924`) — the `d.abs() <= 1e-30` guard just skips the
update, leaving the tiny pivot in `d_diag` unchanged. Downstream,
`solve.rs:144` raises `SingularMatrix` rather than perturbing, and the
inertia count is whatever sign the tiny pivot happened to have.
**This is the single biggest correctness risk for the ipopt pipeline**:
a factorization can succeed with a 1e-20 pivot whose sign is noise,
the IPM reads the inertia as `(n, m, 0)`, and the step direction is
meaningless.

**(b) Zero-column / null pivot handling is silent.**
`find_pivot_threshold` at `pivot.rs:217-219` and `find_pivot_ldlt` at
`frontal.rs:772-773` both treat a fully zero column as "delay". At a
root or must-eliminate position this becomes an accepted pivot with
value 0 and the subsequent `1/d` in the within-block update is skipped
but `d_diag[npiv]=0` is stored. There is no `ICNTL(24)`/`CNTL(3)`
null-pivot list, no `INFOG(28)` count, no bound-preserving fixation
(`CNTL(5)*||A||`). The caller has no way to know a row was
structurally singular.

**(c) "Must-eliminate" static pivoting does not count tiny pivots.**
In MUMPS, every tiny pivot accepted through the static path increments
`NBTINYW`, and any diagonal sign flip increments `NNEGW`. rmumps's
`frontal.rs:877-879` counts only negative pivots into `nneg` and
exposes no "how many pivots were accepted only because we had nowhere
to delay them" metric. An inertia answer built on many forced pivots
is not the inertia of the original matrix.

**(d) Threshold-pivoting and classic-BK code paths diverge.** rmumps
ships two disjoint pivot searches: classic Bunch-Kaufman in
`pivot.rs:135-182` (used by `partial_factor` with no delays) and
MUMPS-style in `frontal.rs:739-892` (used by
`partial_factor_threshold_inner`). `Solver` only exercises the
MUMPS-style one when `pivot_threshold > 0.0` (`numeric.rs:306-310`).
With `pivot_threshold = 0.0` the solver uses the non-MUMPS path and
cannot produce delayed pivots at all — `numeric.rs:291-298` elides the
`promote_cb_to_fs` step when `pivot_threshold == 0.0`. The default is
0.01, so this is the common case, but any caller that sets threshold
to zero silently loses the entire delayed-pivot machinery.

**(e) The symbolic-vs-numeric contract is fragile under dynamic front
expansion.** `extend_add` calls `expand_for_delayed`
(`frontal.rs:161-183`) when a child's contribution references global
indices the parent's symbolic structure did not predict. Each
expansion reallocates and copies the whole dense matrix. Functionally
correct but it means the symbolic phase does not bound memory. MUMPS
sizes fronts at analysis time to accept any delayed pivot from
children (`NE_STEPS`, `ND_STEPS`).

**(f) Scaling fallback path is asymmetric with ordering.**
`solver.rs:116-120` routes to `compute_mc64_kkt_scaling` when
`n_primal.is_some()`, otherwise to generic `compute_scaling`. But the
permutation always comes from `compute_ordering_with_kkt`, and scaling
factors are permuted *after* being computed (`solver.rs:123-130`). If
a caller sets `n_primal` but picks `Ordering::Amd`, the KKT-specific
scaling is applied but the ordering is not KKT-aware — the pairings
that the scaling was built for are scattered across the tree. No code
path rejects this combination.

**(g) Iterative refinement stops too aggressively.**
`solver.rs:218` exits if `res_norm > 0.9 * prev_res_norm`. For KKT
systems with modest stagnation this trips on the second iteration.
MUMPS at `dsol_driver.F` uses the Arioli-Demmel-Duff stopping
criterion based on componentwise backward error, not a geometric
decrease ratio.

**(h) Growth-factor / backward-stability claims.** The module doc at
`lib.rs:1-18` and `MUMPS_LDLT_ALGORITHM.md` describe the factorization
as "Bunch-Kaufman" with "inertia detection". Classic backward-stability
theorems for BK assume pivots satisfying the full alpha test — but
rmumps accepts any pivot that passes the *threshold* test with
`uu=0.01`, same as MUMPS, which relaxes the classical guarantee.
There is no growth-factor monitoring (`DKEEP(2)`, `RINFOG` in MUMPS).

**(i) 2×2 pivot sign accounting in the "must-eliminate" branch.**
`find_pivot_ldlt` at `frontal.rs:864-881` only handles the 1×1 static
fallback. If two consecutive promoted columns both fail and a 2×2
would still be numerically adequate (`|DETPIV| > 0` but
threshold-rejected), rmumps does not attempt a static 2×2 — it just
accepts the first one as 1×1. MUMPS keeps trying the 2×2 under
`STATICMODE`.

### 3. Genuine advantages of rmumps

- **Lock-free level-set factorization.** `numeric.rs:119-165` computes
  topological levels on the supernodal tree and dispatches same-level
  nodes to rayon. MUMPS gets the analogous parallelism only through
  the MPI type-2 node split, which requires cluster-sized problems to
  amortize. The `SyncCell` pattern (`numeric.rs:11-17`) makes the
  invariant explicit: same level ⇒ disjoint writes.
- **Caller-visible inertia, including 2×2 block eigenvalues.**
  `NumericFactorization::min_diagonal` (`numeric.rs:38-66`) returns
  the minimum eigenvalue across all D blocks (smaller root of the 2×2
  when active). MUMPS exposes `RINFOG(12)` (minimum pivot), but to get
  a 2×2 eigenvalue you have to reconstruct it from `DKEEP`. For IPM
  inertia correction this is a real ergonomic win.
- **KKT-matching + compressed AMD as a single orderings pass.**
  `ordering/kkt_matching.rs:61-210` implements a maximum-cardinality
  matching between primals and duals, then feeds the compressed graph
  through AMD. Close to MUMPS's `ICNTL(12)=2` but lives in one file,
  with no external dependency.
- **Explicit, auditable data layout.** `FrontalMatrix`
  (`frontal.rs:17-25`) and `PartialFactorResult` (`frontal.rs:28-42`)
  are plain structs with documented field meanings. MUMPS encodes the
  same information in integer offsets into `IW` at `IOLDPS+XSIZE`,
  `IOLDPS+2+XSIZE`, etc.
- **Re-factor-with-same-pattern fast path.** `Solver::factor_csc`
  (`solver.rs:274-304`) lets the IPM refactorize at each iteration
  without re-running AMD or the symbolic phase. MUMPS has the same
  capability via `JOB=2`, but exposing it cleanly through a Rust API
  removes a whole class of misuse.
- **One-shot build, no Fortran toolchain.** The whole solver builds
  with a Rust compiler plus optional `faer`. MUMPS requires a Fortran
  compiler, BLAS, LAPACK, an MPI stack (even for sequential use), and
  one of {METIS, SCOTCH, PORD}.

### 4. Roadmap recommendations (rmumps-side)

Effort: S = <1 day, M = 1–3 days, L = 1–2 weeks.

1. **Wire up CNTL(4)-equivalent static-pivot threshold.** Replace
   `let seuil: f64 = 0.0` at `frontal.rs:745` with a parameter plumbed
   from `SolverOptions`. In the "must-eliminate" branch at
   `frontal.rs:864-881`, replace any pivot with `|d| < seuil` by
   `sign(d) * max(seuil, eps * ||A||)` and count it into a
   `TinyPivotCount` field on `NumericFactorization`. **Effort: S.**
   Addresses (a), (c). Likely fixes IPM silent-failure modes on
   CONCON-style near-singular KKTs and HS13 suboptimal termination.

2. **Add ICNTL(24)-equivalent null-pivot detection.** Define
   `CNTL(3)` (structural zero threshold) and `CNTL(5)` (fixation
   magnitude) in `SolverOptions`. In `find_pivot_ldlt`, when
   `col_max <= CNTL(3)`, add the column to a `Vec<usize>` of null
   pivots on `NumericFactorization` and fix the diagonal to
   `CNTL(5) * ||A||`. **Effort: S–M.** Addresses (b).

3. **Report NBTINY / growth-factor info.** Add counters (`nb_tiny`,
   `nb_static`, `max_growth`) to `NumericFactorization`, incremented
   inside `find_pivot_ldlt` and `update_within_block_ldlt`. MUMPS
   keeps these in `RINFOG` / `INFOG(25)`. Without them the IPM cannot
   distinguish "factorization was clean" from "we papered over 40
   tiny pivots". **Effort: S.** Addresses (c), (h).

4. **Unify the two pivot-search paths.** Delete
   `dense_ldlt_bunch_kaufman` / `partial_factor`
   (`pivot.rs:334-443`, `frontal.rs:243-437`) as the default code
   path, or make them explicitly for tests only. Always go through
   `partial_factor_threshold_inner`. **Effort: M.** Addresses (d).

5. **Port MUMPS's analysis-time front sizing.** In `symbolic.rs`,
   compute upper bounds on front sizes that include potential delayed
   pivots from children (MUMPS's `ND_STEPS = NE_STEPS + delayed_budget`).
   Eliminates `expand_for_delayed` at runtime. Reference:
   `ana_driver.F` and `ana_reordertree.F`. **Effort: L.** Addresses
   (e); improves predictability for large CUTEst problems.

6. **Real MC64 scaling.** `dana_mtrans.F` (1196 lines) is the
   authoritative implementation; for KKT-only use, an auction-algorithm
   MC64 is ~300 lines (Duff–Koster 1999). The current
   `compute_mc64_kkt_scaling` in `scaling.rs` is a bespoke heuristic,
   not MC64. Real MC64 would let the pivot threshold drop from 0.01
   toward Ipopt's 1e-6 without loss of stability. **Effort: L.**
   Addresses (f) indirectly.

7. **Backward-error-based iterative refinement stop.** Replace the
   `0.9 * prev_res_norm` heuristic at `solver.rs:218` with an
   Arioli–Demmel–Duff criterion: stop when `omega1 + omega2 < eps` or
   when they fail to decrease by more than a factor. MUMPS at
   `dsol_driver.F:5477-5610` is the reference. **Effort: S.**
   Addresses (g); may help CONCON stalls.

8. **Allow METIS and return to natural when unavailable.** Add an
   optional `metis-sys` dependency with a clean compile-time feature
   flag; `Ordering::Metis` should fall back to AMD with a warning
   rather than panicking. MUMPS does this automatically in
   `ana_set_ordering.F`. **Effort: M.** Addresses large-problem
   slowdown: AMD on a 50k-variable CUTEst problem is measurably worse
   than METIS for both fill and flop count.

9. **Static 2×2 pivots under `AVOID_DELAYED`.** Extend the
   must-eliminate branch in `find_pivot_ldlt` to try a 2×2 before
   falling back to 1×1 static. Follow `dfac_front_aux.F:1231-1240`
   where `STATICMODE=.TRUE.` adjusts the threshold but preserves 2×2
   search. **Effort: M.** Addresses (i).

10. **Expose an inertia-correction API that distinguishes numerical
    from structural failure.** If static pivoting were tracked
    (items 1, 3), the IPM could bump `δ` *before* re-factoring once
    the tiny-pivot count exceeds a threshold, saving a factorization.
    **Effort: S once items 1–3 land.**

Items 1–4 are the prerequisites for trusting rmumps on non-trivial
KKTs. Items 5–7 address performance and convergence at scale. Items
8–10 bring the API surface closer to MUMPS's operational model.

---

## Key files referenced

### ripopt / Ipopt
- ripopt: `src/convergence.rs`, `src/kkt.rs`, `src/filter.rs`, `src/ipm.rs`,
  `src/preprocessing.rs`, `src/restoration_nlp.rs`, `src/restoration.rs`,
  `src/sensitivity.rs`, `src/warmstart.rs`, `src/options.rs`,
  `src/slack_formulation.rs`
- Ipopt: `ref/Ipopt/src/Algorithm/IpIpoptCalculatedQuantities.cpp`,
  `IpOptErrorConvCheck.cpp`, `IpPDPerturbationHandler.cpp`,
  `IpFilterLSAcceptor.cpp`, `IpBacktrackingLineSearch.cpp`,
  `IpRestoIpoptNLP.cpp`, `IpRestoIterateInitializer.cpp`,
  `IpTNLPAdapter.cpp`, `contrib/sIPOPT/src/SensAlgorithm.cpp`

### rmumps / MUMPS
- rmumps: `rmumps/src/pivot.rs`, `frontal.rs`, `numeric.rs`, `solver.rs`,
  `scaling.rs`, `symbolic.rs`, `ordering/kkt_matching.rs`,
  `ordering/amd.rs`, `solve.rs`
- MUMPS: `ref/mumps/src/dfac_front_LDLT_type1.F`, `dfac_front_aux.F`
  (`DMUMPS_FAC_I_LDLT`, `DMUMPS_FAC_MQ_LDLT`, `DMUMPS_FAC_SQ_LDLT`,
  `DMUMPS_SWAP_LDLT`), `dfac_scalings.F`, `dana_mtrans.F`,
  `dsol_driver.F`, `dini_defaults.F`
