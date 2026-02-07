// Benchmark binary: solve all HS problems with ripopt and output JSON results.

#[path = "generated/hs_problems.rs"]
mod hs_problems;

use hs_problems::solve_all;
use ripopt::SolverOptions;

fn main() {
    let options = SolverOptions {
        tol: 1e-8,
        max_iter: 3000,
        print_level: 0,
        mu_strategy_adaptive: true,
        ..SolverOptions::default()
    };

    eprintln!("Solving all HS problems with ripopt...");
    let results = solve_all(&options);

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
