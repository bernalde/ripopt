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
        "mu_strategy" => {
            opts.mu_strategy_adaptive = value == "adaptive";
        }
        "warm_start_init_point" => {
            opts.warm_start = value == "yes";
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
        _ => {
            eprintln!("Warning: unknown option '{}'", key);
        }
    }
}
