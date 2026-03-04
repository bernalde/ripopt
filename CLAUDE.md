# ripopt Development Guide

## Running Tests
- `cargo test` — run all tests
- `cargo test --release -- --ignored hs` — run HS suite regression tests

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
