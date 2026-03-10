# ripopt Autonomous Optimization Program

You are an autonomous optimization agent. Your goal is to solve a nonlinear
optimization problem using ripopt by iteratively adjusting solver options and
starting points, guided by structured diagnostics.

## Setup

The problem to solve is specified by the `PROBLEM` variable below. Set it
before starting.

```
PROBLEM=debug_tp374        # Rust example name (cargo run --example $PROBLEM)
MAX_ATTEMPTS=10             # Maximum experiment attempts
MAX_WALL_TIME=30            # Max seconds per solve attempt
```

## Experiment Log

All results are recorded in `experiments.jsonl` — one JSON object per line.
This is your memory across attempts. **Never delete this file.** Append only.

Each entry has this shape:

```json
{
  "attempt": 1,
  "options": {"mu_init": 0.1, "kappa": 10.0},
  "x0_change": null,
  "status": "MaxIterations",
  "objective": 0.877,
  "iterations": 2999,
  "diagnostics": {
    "filter_rejects": 47,
    "restoration_count": 8,
    "nlp_restoration_count": 1,
    "mu_mode_switches": 22,
    "final_mu": 1.2e-2,
    "final_primal_inf": 0.12,
    "final_dual_inf": 3.4e-3,
    "final_compl": 1.1e-2,
    "wall_time_secs": 4.2,
    "fallback_used": null
  },
  "notes": "Baseline run with default options"
}
```

## Protocol

### Step 0: Read history

If `experiments.jsonl` exists, read it. Identify:
- How many attempts have been made
- The best result so far (lowest objective among Optimal/Acceptable entries)
- What options and starting points have already been tried
- Any patterns in the diagnostics across attempts

If MAX_ATTEMPTS have been reached, stop and report the best result.

### Step 1: Decide what to try

Based on the experiment history (or lack thereof), pick the next experiment.
Use the decision rules below. **Do not ask the user — just decide.**

**If this is attempt 1:** Run with default SolverOptions (baseline).

**If the previous attempt failed, read its diagnostics and apply the FIRST
matching rule:**

| Priority | Diagnostic pattern | Action |
|---|---|---|
| 1 | `filter_rejects` > 10 | Increase `mu_init` to 1.0 (or 10x previous) |
| 2 | `restoration_count` > 5 | Enable `enable_slack_fallback: true` |
| 3 | `mu_mode_switches` > 15 | Set `mu_strategy_adaptive: false` |
| 4 | `final_mu` > 1e-3 | Reduce `kappa` to 3.0, reduce `mu_linear_decrease_factor` to 0.1 |
| 5 | `final_primal_inf` > 0.1 | Try new starting point (see Starting Point Rules) |
| 6 | `watchdog_activations` > 0 | Set `hessian_approximation_lbfgs: true` |
| 7 | `fallback_used` is set | Disable that fallback, try `mehrotra_pc: true` |
| 8 | None of the above | Combine: `mu_init: 1.0, kappa: 3.0, mehrotra_pc: true, gondzio_mcc_max: 3` |

**If a rule was already tried and didn't help:** Skip it, try the next rule.
The experiment log tells you what's been tried.

**If all rules have been tried:** Try starting point changes (see below).

### Starting Point Rules

Apply these in order when options alone aren't working:

1. **Constraint-informed**: Analyze the constraint functions. Pick `x0` values
   that approximately satisfy the most binding constraints. Write reasoning in
   the `notes` field.

2. **Perturbed from best**: Take `result.x` from the best attempt so far.
   Perturb each component by ±10%. Use `warm_start: true`.

3. **Scaled default**: If default `x0` is uniform (e.g., all 0.1), try
   non-uniform values spanning the variable bounds.

4. **Zero-centered**: Try `x0 = 0` (or midpoint of bounds if bounded).

### Step 2: Apply changes

Edit the example source file to apply the chosen SolverOptions and/or starting
point. The changes should be minimal — only modify the `SolverOptions` struct
and/or `initial_point()` body.

Set `print_level: 5` and `max_wall_time` to the configured limit.

### Step 3: Run

```bash
cargo run --example $PROBLEM 2>&1
```

### Step 4: Record

Parse the `--- ripopt diagnostics ---` block from stderr. Extract all fields.
Append a JSON line to `experiments.jsonl` with:

- `attempt`: sequential number
- `options`: the SolverOptions that differ from defaults (as a dict)
- `x0_change`: description of starting point change, or null
- `status`, `objective`, `iterations`: from the result
- `diagnostics`: all fields from the diagnostics block
- `notes`: your reasoning for this attempt (1-2 sentences)

### Step 5: Evaluate and loop

Compare this result to the best so far:

- **Optimal**: You're done. Report success.
- **Acceptable**: Record as best so far, but keep trying for Optimal.
- **Better objective than previous best** (even if not converged): Note improvement.
- **Worse or same**: Note that this approach didn't help.

If attempts < MAX_ATTEMPTS, go to Step 0.

If attempts >= MAX_ATTEMPTS, stop and report:
- Best result found (status, objective, iterations, which attempt)
- Summary table of all attempts
- What worked and what didn't

## Output Format

At the end, produce a summary like:

```
## Optimization Results: $PROBLEM

Best: Attempt 3, status=Acceptable, obj=0.2341, iters=847
Known optimal: 0.233264

| # | Status | Objective | Iters | Key change | Notes |
|---|--------|-----------|-------|------------|-------|
| 1 | MaxIter | 0.877 | 2999 | baseline | High filter_rejects |
| 2 | MaxIter | 0.544 | 2999 | mu_init=1.0 | Reduced filter_rejects |
| 3 | Accept | 0.234 | 847 | +slack+x0 | Near optimal |
```

## Important Rules

1. **Never ask the user questions.** Use the decision rules above.
2. **Never delete experiments.jsonl.** Append only.
3. **Always record every attempt**, even failures and errors.
4. **Don't try the same configuration twice.** Read the log.
5. **Revert source changes between attempts** — each attempt starts from the
   original source, applies its own changes cleanly.
6. **Keep wall time bounded.** Use `max_wall_time` to prevent runaway solves.
7. **Be specific in notes.** Future-you reads these to avoid repeating mistakes.
