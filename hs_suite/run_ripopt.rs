// Benchmark binary: solve all HS problems with ripopt and output JSON results.
// Runs multiple timing passes to get stable timing measurements.

#[path = "generated/hs_problems.rs"]
mod hs_problems;

use hs_problems::solve_all;
use ripopt::SolverOptions;

fn main() {
    let n_timing_runs: usize = std::env::var("RIPOPT_TIMING_RUNS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    let options = SolverOptions {
        tol: 1e-8,
        max_iter: 3000,
        print_level: 0,
        mu_strategy_adaptive: true,
        ..SolverOptions::default()
    };

    eprintln!("Solving all HS problems with ripopt ({} timing runs)...", n_timing_runs);

    // First run: get correctness results
    let mut results = solve_all(&options);

    // Additional timing runs: keep minimum solve_time per problem
    for run in 1..n_timing_runs {
        eprintln!("  Timing run {}/{}...", run + 1, n_timing_runs);
        let timing_results = solve_all(&options);
        for (r, t) in results.iter_mut().zip(timing_results.iter()) {
            if t.solve_time < r.solve_time {
                r.solve_time = t.solve_time;
            }
        }
    }

    // Summary to stderr
    let total = results.len();
    let optimal = results
        .iter()
        .filter(|r| r.status == "Optimal")
        .count();
    let acceptable = results
        .iter()
        .filter(|r| r.status == "Acceptable")
        .count();
    let solved = optimal + acceptable;
    eprintln!(
        "Solved {}/{} ({} optimal, {} acceptable)",
        solved, total, optimal, acceptable
    );

    // JSON to stdout
    let json = serde_json::to_string_pretty(&results).unwrap();
    println!("{}", json);
}
