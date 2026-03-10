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
  "hypothesis": "Baseline run to establish diagnostic profile",
  "outcome": "MaxIterations — high filter_rejects (47) and restoration_count (8) indicate the filter line search is struggling. mu stuck at 1.2e-2."
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

### Step 1b: Explain your reasoning

Before running each attempt, print a reasoning block to the user:

```
## Attempt N: [short title, e.g. "Increase mu_init"]

**What:** [1 sentence describing the concrete change being made]
**Why:** [1-2 sentences explaining the hypothesis — what diagnostic pattern
triggered this choice, what the change does mechanically in the solver, and
why that should improve convergence]
**Risk:** [1 sentence on what could go wrong or why this might not help]
```

This is mandatory for every attempt including the baseline. For the baseline:
```
## Attempt 1: Baseline

**What:** Running with default SolverOptions to establish a diagnostic baseline.
**Why:** We need to see how the solver behaves on this problem before making
any adjustments. The diagnostics will reveal the primary bottleneck.
**Risk:** The problem may be hard enough that defaults fail entirely, but we
need this data to guide subsequent attempts.
```

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
- `hypothesis`: why this change should help (written before running)
- `outcome`: what actually happened and what the diagnostics reveal (written after running)

### Step 5: Analyze and explain

After each attempt, print a post-run analysis:

```
**Result:** [status] — objective [value] in [N] iterations
**Diagnosis:** [1-2 sentences interpreting the key diagnostics — what they
reveal about solver behavior this time. Reference specific numbers.]
**Verdict:** [Improved/No change/Worse] compared to best so far — [why,
citing specific metric changes, e.g. "filter_rejects dropped from 47 to 3
but restoration_count increased from 8 to 12"]
```

Then compare this result to the best so far:

- **Optimal**: You're done. Report success (see Output Format).
- **Acceptable**: Record as best so far, but keep trying for Optimal.
- **Better objective than previous best** (even if not converged): Note improvement.
- **Worse or same**: Note that this approach didn't help.

If attempts < MAX_ATTEMPTS, go to Step 0.

If attempts >= MAX_ATTEMPTS, stop and report (see Output Format).

## Output Format

### On success (Optimal found)

```
## Optimization Results: $PROBLEM — SOLVED

Solved at attempt N with status=Optimal, obj=[value], iters=[N].

### How we got here
[2-3 sentences tracing the key turning points across attempts. What was the
main bottleneck and what change overcame it?]

| # | Status | Objective | Iters | Key change | Verdict |
|---|--------|-----------|-------|------------|---------|
| 1 | MaxIter | 0.877 | 2999 | baseline | High filter_rejects |
| 2 | MaxIter | 0.544 | 2999 | mu_init=1.0 | Improved but still stuck |
| 3 | Optimal | 0.233 | 847 | +slack+x0 | Solved |
```

### On failure (MAX_ATTEMPTS exhausted)

```
## Optimization Results: $PROBLEM — NOT SOLVED

Best: Attempt N, status=[status], obj=[value], iters=[N]

| # | Status | Objective | Iters | Key change | Verdict |
|---|--------|-----------|-------|------------|---------|
| ... | ... | ... | ... | ... | ... |

### Why this problem is hard
[Explain the structural difficulty based on what you observed across all
attempts. Reference specific diagnostic patterns. Examples of explanations:
- "Non-convex with multiple local minima — the solver converges to different
  stationary points depending on the starting point (attempts 3,5,7 found
  three distinct local solutions)."
- "Severely ill-conditioned constraints — restoration_count remained high
  (>10) across all attempts regardless of mu_init or slack reformulation,
  suggesting near-degenerate constraint qualification."
- "Barrier parameter unable to decrease — final_mu never dropped below 1e-3
  in any attempt, indicating the central path is disrupted, possibly by
  complementarity degeneracy."]

### What we learned
[Summarize which approaches showed improvement and which didn't, with
specific metrics. Identify the primary bottleneck that blocked convergence.
Example: "Increasing mu_init reduced filter_rejects from 47 to 3 (attempt 2)
but restoration_count remained at 12. Slack reformulation (attempt 4) reduced
restoration_count to 2 but objective worsened. The core issue appears to be
a narrow feasible region where the barrier subproblems become ill-conditioned
as mu decreases."]

### Suggestions for further investigation
[2-3 concrete, specific next steps. Not generic advice — reference the
actual diagnostic patterns observed. Examples:
- "Try bound_push=1e-2 to keep iterates further from bounds (final_compl
  was 1e-1 in all attempts, suggesting bound proximity issues)."
- "The problem may benefit from a penalty method instead of barrier —
  consider an external SQP solver."
- "Attempt 6 found obj=0.234 at a local minimum. A multistart with 20
  random initial points covering the full bound range may find the basin
  of attraction for the global optimum."]
```

## Variant: .nl File Problems

When the problem is an `.nl` file (from AMPL, CUTEst, or Pyomo) instead of a
Rust example, the protocol is the same but the mechanics differ:

### Setup

```
NL_FILE=problem.nl             # Path to the .nl file
MAX_ATTEMPTS=10
MAX_WALL_TIME=30
```

### Running (replaces Step 3)

Use the `ripopt_ampl` binary with key=value options on the command line:

```bash
cargo run --bin ripopt_ampl -- $NL_FILE print_level=5 max_wall_time=30 2>&1
```

To change solver options, pass them as key=value arguments:

```bash
cargo run --bin ripopt_ampl -- $NL_FILE \
  print_level=5 max_wall_time=30 \
  mu_init=1.0 kappa=3.0 mu_strategy=adaptive \
  slack_fallback=yes mehrotra_pc=yes 2>&1
```

Available CLI options: `tol`, `max_iter`, `acceptable_tol`, `mu_init`,
`print_level`, `max_wall_time`, `bound_push`, `kappa`,
`mu_linear_decrease_factor`, `mu_strategy` (adaptive/monotone),
`warm_start_init_point` (yes/no), `max_soc`, `slack_fallback` (yes/no),
`al_fallback` (yes/no), `sqp_fallback` (yes/no), `mehrotra_pc` (yes/no),
`gondzio_mcc_max`, `proactive_infeasibility_detection` (yes/no),
`hessian_approximation` (limited-memory/exact).

### Applying changes (replaces Step 2)

**Options:** No source editing needed — pass options as CLI arguments.

**Starting points:** The `.nl` file contains the initial point in its `x`
segment. To change the starting point:

1. **Edit the `.nl` file directly.** The `x` segment starts with a line `x<N>`
   (where N is the number of variables) followed by lines of the form
   `<index> <value>`. Modify the values. Keep a backup of the original.

2. **Warm-start from a previous solution.** After a solve, ripopt writes a
   `.sol` file. To warm-start from it, pass `warm_start_init_point=yes`. The
   AMPL driver reads the `.sol` file's primal values as the new starting point.

3. **Generate a new `.nl` file** from the modeling language (AMPL/Pyomo) with
   different initial values, if the source model is available.

### Recording (same as Step 4)

The diagnostics block format is identical. Parse `--- ripopt diagnostics ---`
from stderr. The experiment log entry should include:

```json
{
  "attempt": 2,
  "mode": "nl",
  "nl_file": "problem.nl",
  "cli_options": "mu_init=1.0 kappa=3.0",
  "x0_change": "edited x segment: x0=[1.0, 0.5, 0.5, 1.0]",
  "status": "MaxIterations",
  "objective": 0.544,
  "iterations": 2999,
  "diagnostics": { ... },
  "hypothesis": "Higher mu_init should reduce filter rejects by starting further from the boundary",
  "outcome": "filter_rejects dropped from 47 to 3, but restoration_count increased to 12"
}
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
