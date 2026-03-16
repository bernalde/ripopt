# ripopt Development Guide

## Running Tests
- `cargo test` — run all tests (230 tests, ~3 seconds)

## Benchmarks
- `make benchmark` — full benchmark: HS + CUTEst + large-scale + report
- `make benchmark-report` — regenerate report from existing results
- `make hs-run` — HS suite only (ripopt + ipopt, ~2 min)
- `cargo run --release --features hs --bin hs_suite` — HS suite (ripopt only)
- Individual CUTEst problems: `cargo run --bin cutest_suite --features cutest,ipopt-native --release -- PROBLEM1 PROBLEM2`
- Full CUTEst suite: `RESULTS_FILE=cutest_suite/results.json cargo run --bin cutest_suite --features cutest,ipopt-native --release`

## Code Coverage
- Install: `cargo install cargo-llvm-cov` (requires llvm-tools-preview: `rustup component add llvm-tools-preview`)
- Summary: `cargo llvm-cov test` (runs tests with instrumentation and prints coverage)
- Line-by-line: `cargo llvm-cov test --text`
- HTML: `cargo llvm-cov test --html`
- After adding tests, update the coverage table in README.md with current numbers from `cargo llvm-cov test`

## Test Guidelines
- Every test must exercise a specific code branch (not just return true)
- Tests must verify correctness: check status codes, objective values, solution quality
- Prefer small hand-crafted problems over large benchmark problems
- Keep test execution under 1 second per test
- Unit tests (`#[cfg(test)] mod tests`) for module-internal functions
- Integration tests (`tests/`) for cross-module behavior and solver paths

## Honesty in Benchmarks and Tests
**No misleading benchmarks or problem-specific hacks.** The following are explicitly prohibited:
- Counting `NumericalError`, `MaxIterations`, or any non-`Optimal` status as a "solve" in benchmark summaries
- Writing tests that accept failure statuses (e.g., `|| NumericalError`) just to make the pass rate look better
- Tuning solver parameters specifically for individual benchmark problems to inflate scores
- Adding special-case code paths triggered only by specific problem structures seen in benchmarks
- Hiding known failures behind lenient statuses (`Acceptable` was removed for this reason)

**Tests must be honest:** If the solver cannot solve a problem to `Optimal`, the test should either fail (exposing the real limitation), be marked `#[ignore]` with a clear explanation, or be removed. A failing test that documents a real limitation is more valuable than a passing test that hides one.
