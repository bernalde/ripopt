// Quick regression diagnostic
// Run: cargo run --release --example regression_test

#[path = "../hs_suite/generated/hs_problems.rs"]
mod hs_problems;

use ripopt::SolverOptions;

macro_rules! test_configs {
    ($pname:expr, $problem:expr, $configs:expr) => {
        for (cname, opts) in $configs.iter() {
            let result = ripopt::solve(&$problem, opts);
            let status = format!("{:?}", result.status);
            let marker = if status == "Optimal" || status == "Acceptable" { "" } else { " <-- FAIL" };
            println!("{:<6} {:<30} {:<22} {:>10.4e} {:>6}{}",
                $pname, cname, status, result.objective, result.iterations, marker);
        }
        println!();
    };
}

fn main() {
    let base = SolverOptions {
        tol: 1e-8,
        max_iter: 3000,
        print_level: 0,
        mu_strategy_adaptive: true,
        least_squares_mult_init: false,
        constraint_slack_barrier: false,
        mu_allow_increase: false,
        ..SolverOptions::default()
    };

    println!("{:<6} {:<30} {:<22} {:>10} {:>6}", "TP#", "Options", "Status", "Objective", "Iters");
    println!("{}", "-".repeat(80));

    // TP116: test all pairwise combinations
    let configs: Vec<(&str, SolverOptions)> = vec![
        ("old defaults (all OFF)", SolverOptions {
            least_squares_mult_init: false,
            constraint_slack_barrier: false,
            mu_allow_increase: false,
            ..base.clone()
        }),
        ("all new features ON", SolverOptions {
            least_squares_mult_init: true,
            constraint_slack_barrier: true,
            mu_allow_increase: true,
            ..base.clone()
        }),
        ("ls_mult + slack_barrier", SolverOptions {
            least_squares_mult_init: true,
            constraint_slack_barrier: true,
            mu_allow_increase: false,
            ..base.clone()
        }),
        ("ls_mult + mu_increase", SolverOptions {
            least_squares_mult_init: true,
            constraint_slack_barrier: false,
            mu_allow_increase: true,
            ..base.clone()
        }),
        ("slack_barrier + mu_increase", SolverOptions {
            least_squares_mult_init: false,
            constraint_slack_barrier: true,
            mu_allow_increase: true,
            ..base.clone()
        }),
        ("only ls_mult_init", SolverOptions {
            least_squares_mult_init: true,
            constraint_slack_barrier: false,
            mu_allow_increase: false,
            ..base.clone()
        }),
        ("only constr_slack_barrier", SolverOptions {
            least_squares_mult_init: false,
            constraint_slack_barrier: true,
            mu_allow_increase: false,
            ..base.clone()
        }),
        ("only mu_allow_increase", SolverOptions {
            least_squares_mult_init: false,
            constraint_slack_barrier: false,
            mu_allow_increase: true,
            ..base.clone()
        }),
    ];

    test_configs!("TP116", hs_problems::HsTp116, configs);
}
