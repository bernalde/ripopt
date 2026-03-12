use ripopt::nl::{parse_nl_file, write_sol, NlProblem};
use ripopt::SolverOptions;
use std::fs;
use std::io::BufWriter;

fn main() {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();

    // Handle -v / --version flag for AMPL solver protocol
    if args.len() >= 2 && (args[1] == "-v" || args[1] == "--version" || args[1] == "-AMPL") {
        if args.len() == 2 && (args[1] == "-v" || args[1] == "--version") {
            println!("ripopt {}", env!("CARGO_PKG_VERSION"));
            return;
        }
    }

    if args.len() < 2 {
        eprintln!("Usage: ripopt <problem.nl> [-AMPL] [key=value ...]");
        std::process::exit(1);
    }

    let nl_path = &args[1];
    let mut options = SolverOptions::default();

    // Parse key=value options from command line
    for arg in &args[2..] {
        if arg == "-AMPL" || arg == "--AMPL" {
            continue; // AMPL mode flag, acknowledged
        }
        if let Some((key, value)) = arg.split_once('=') {
            apply_option(&mut options, key.trim(), value.trim());
        }
    }

    // Read and parse NL file
    let content = match fs::read_to_string(nl_path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error reading {}: {}", nl_path, e);
            std::process::exit(1);
        }
    };

    let nl_data = match parse_nl_file(&content) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error parsing NL file: {}", e);
            std::process::exit(1);
        }
    };

    let n_vars = nl_data.header.n_vars;
    let n_constrs = nl_data.header.n_constrs;

    let problem = NlProblem::from_nl_data(nl_data);

    // Solve
    let result = ripopt::solve(&problem, &options);

    // Print summary to stdout
    println!(
        "ripopt {}: {} after {} iterations",
        env!("CARGO_PKG_VERSION"),
        match result.status {
            ripopt::SolveStatus::Optimal => "Optimal",
            ripopt::SolveStatus::Acceptable => "Acceptable",
            ripopt::SolveStatus::Infeasible => "Infeasible",
            ripopt::SolveStatus::LocalInfeasibility => "LocalInfeasibility",
            ripopt::SolveStatus::MaxIterations => "MaxIterations",
            ripopt::SolveStatus::NumericalError => "NumericalError",
            ripopt::SolveStatus::Unbounded => "Unbounded",
            ripopt::SolveStatus::RestorationFailed => "RestorationFailed",
            ripopt::SolveStatus::InternalError => "InternalError",
        },
        result.iterations
    );
    println!("Objective: {:.15e}", result.objective);

    // Print diagnostics to stderr
    result.diagnostics.print_summary(result.status, result.iterations);

    // Write SOL file (replace .nl extension with .sol)
    let sol_path = if nl_path.ends_with(".nl") {
        format!("{}sol", &nl_path[..nl_path.len() - 2])
    } else {
        format!("{}.sol", nl_path)
    };

    let sol_file = match fs::File::create(&sol_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating {}: {}", sol_path, e);
            std::process::exit(1);
        }
    };

    let mut writer = BufWriter::new(sol_file);
    if let Err(e) = write_sol(&mut writer, &result, n_vars, n_constrs) {
        eprintln!("Error writing SOL file: {}", e);
        std::process::exit(1);
    }
}

/// Apply a key=value option to SolverOptions.
fn apply_option(opts: &mut SolverOptions, key: &str, value: &str) {
    match key {
        "tol" => {
            if let Ok(v) = value.parse() {
                opts.tol = v;
            }
        }
        "max_iter" => {
            if let Ok(v) = value.parse() {
                opts.max_iter = v;
            }
        }
        "acceptable_tol" => {
            if let Ok(v) = value.parse() {
                opts.acceptable_tol = v;
            }
        }
        "acceptable_iter" => {
            if let Ok(v) = value.parse() {
                opts.acceptable_iter = v;
            }
        }
        "acceptable_constr_viol_tol" => {
            if let Ok(v) = value.parse() {
                opts.acceptable_constr_viol_tol = v;
            }
        }
        "acceptable_dual_inf_tol" => {
            if let Ok(v) = value.parse() {
                opts.acceptable_dual_inf_tol = v;
            }
        }
        "acceptable_compl_inf_tol" => {
            if let Ok(v) = value.parse() {
                opts.acceptable_compl_inf_tol = v;
            }
        }
        "mu_init" => {
            if let Ok(v) = value.parse() {
                opts.mu_init = v;
            }
        }
        "print_level" => {
            if let Ok(v) = value.parse() {
                opts.print_level = v;
            }
        }
        "max_wall_time" => {
            if let Ok(v) = value.parse() {
                opts.max_wall_time = v;
            }
        }
        "bound_push" => {
            if let Ok(v) = value.parse() {
                opts.bound_push = v;
            }
        }
        "bound_frac" => {
            if let Ok(v) = value.parse() {
                opts.bound_frac = v;
            }
        }
        "slack_bound_push" => {
            if let Ok(v) = value.parse() {
                opts.slack_bound_push = v;
            }
        }
        "slack_bound_frac" => {
            if let Ok(v) = value.parse() {
                opts.slack_bound_frac = v;
            }
        }
        "constr_viol_tol" => {
            if let Ok(v) = value.parse() {
                opts.constr_viol_tol = v;
            }
        }
        "dual_inf_tol" => {
            if let Ok(v) = value.parse() {
                opts.dual_inf_tol = v;
            }
        }
        "compl_inf_tol" => {
            if let Ok(v) = value.parse() {
                opts.compl_inf_tol = v;
            }
        }
        "kappa" => {
            if let Ok(v) = value.parse() {
                opts.kappa = v;
            }
        }
        "mu_linear_decrease_factor" => {
            if let Ok(v) = value.parse() {
                opts.mu_linear_decrease_factor = v;
            }
        }
        "mu_superlinear_decrease_power" => {
            if let Ok(v) = value.parse() {
                opts.mu_superlinear_decrease_power = v;
            }
        }
        "mu_min" => {
            if let Ok(v) = value.parse() {
                opts.mu_min = v;
            }
        }
        "tau_min" => {
            if let Ok(v) = value.parse() {
                opts.tau_min = v;
            }
        }
        "mu_allow_increase" => {
            opts.mu_allow_increase = value == "yes" || value == "true" || value == "1";
        }
        "adaptive_mu_monotone_init_factor" => {
            if let Ok(v) = value.parse() {
                opts.adaptive_mu_monotone_init_factor = v;
            }
        }
        "barrier_tol_factor" => {
            if let Ok(v) = value.parse() {
                opts.barrier_tol_factor = v;
            }
        }
        "mu_strategy" => {
            opts.mu_strategy_adaptive = value == "adaptive";
        }
        "least_squares_mult_init" => {
            opts.least_squares_mult_init = value == "yes" || value == "true" || value == "1";
        }
        "constr_mult_init_max" => {
            if let Ok(v) = value.parse() {
                opts.constr_mult_init_max = v;
            }
        }
        "constraint_slack_barrier" => {
            opts.constraint_slack_barrier = value == "yes" || value == "true" || value == "1";
        }
        "sparse_threshold" => {
            if let Ok(v) = value.parse() {
                opts.sparse_threshold = v;
            }
        }
        "warm_start_init_point" => {
            opts.warm_start = value == "yes" || value == "true" || value == "1";
        }
        "warm_start_bound_push" => {
            if let Ok(v) = value.parse() {
                opts.warm_start_bound_push = v;
            }
        }
        "warm_start_bound_frac" => {
            if let Ok(v) = value.parse() {
                opts.warm_start_bound_frac = v;
            }
        }
        "warm_start_mult_bound_push" => {
            if let Ok(v) = value.parse() {
                opts.warm_start_mult_bound_push = v;
            }
        }
        "nlp_lower_bound_inf" => {
            if let Ok(v) = value.parse() {
                opts.nlp_lower_bound_inf = v;
            }
        }
        "nlp_upper_bound_inf" => {
            if let Ok(v) = value.parse() {
                opts.nlp_upper_bound_inf = v;
            }
        }
        "max_soc" => {
            if let Ok(v) = value.parse() {
                opts.max_soc = v;
            }
        }
        "watchdog_shortened_iter_trigger" => {
            if let Ok(v) = value.parse() {
                opts.watchdog_shortened_iter_trigger = v;
            }
        }
        "watchdog_trial_iter_max" => {
            if let Ok(v) = value.parse() {
                opts.watchdog_trial_iter_max = v;
            }
        }
        "restoration_max_iter" => {
            if let Ok(v) = value.parse() {
                opts.restoration_max_iter = v;
            }
        }
        "disable_nlp_restoration" => {
            opts.disable_nlp_restoration = value == "yes" || value == "true" || value == "1";
        }
        "slack_fallback" | "enable_slack_fallback" => {
            opts.enable_slack_fallback = value == "yes" || value == "true" || value == "1";
        }
        "lbfgs_fallback" | "enable_lbfgs_fallback" => {
            opts.enable_lbfgs_fallback = value == "yes" || value == "true" || value == "1";
        }
        "al_fallback" | "enable_al_fallback" => {
            opts.enable_al_fallback = value == "yes" || value == "true" || value == "1";
        }
        "sqp_fallback" | "enable_sqp_fallback" => {
            opts.enable_sqp_fallback = value == "yes" || value == "true" || value == "1";
        }
        "lbfgs_hessian_fallback" | "enable_lbfgs_hessian_fallback" => {
            opts.enable_lbfgs_hessian_fallback = value == "yes" || value == "true" || value == "1";
        }
        "enable_preprocessing" => {
            opts.enable_preprocessing = value == "yes" || value == "true" || value == "1";
        }
        "detect_linear_constraints" => {
            opts.detect_linear_constraints = value == "yes" || value == "true" || value == "1";
        }
        "mehrotra_pc" => {
            opts.mehrotra_pc = value == "yes" || value == "true" || value == "1";
        }
        "gondzio_mcc_max" => {
            if let Ok(v) = value.parse() {
                opts.gondzio_mcc_max = v;
            }
        }
        "proactive_infeasibility_detection" => {
            opts.proactive_infeasibility_detection = value == "yes" || value == "true" || value == "1";
        }
        "hessian_approximation" => {
            match value {
                "limited-memory" => opts.hessian_approximation_lbfgs = true,
                "exact" => opts.hessian_approximation_lbfgs = false,
                _ => eprintln!("Warning: unknown hessian_approximation '{}'", value),
            }
        }
        "linear_solver" => {
            match value {
                "direct" => opts.linear_solver = ripopt::LinearSolverChoice::Direct,
                "iterative" | "minres" => opts.linear_solver = ripopt::LinearSolverChoice::Iterative,
                "hybrid" | "auto" => opts.linear_solver = ripopt::LinearSolverChoice::Hybrid,
                _ => eprintln!("Warning: unknown linear_solver '{}' (use 'direct', 'iterative', or 'hybrid')", value),
            }
        }
        "stall_iter_limit" => {
            if let Ok(v) = value.parse() {
                opts.stall_iter_limit = v;
            }
        }
        _ => {
            eprintln!("Warning: unknown option '{}'", key);
        }
    }
}
