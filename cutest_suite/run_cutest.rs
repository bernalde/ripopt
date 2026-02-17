// Benchmark binary: solve CUTEst problems with both ripopt and native Ipopt,
// output JSON results for comparison.
//
// Usage:
//   cargo run --bin cutest_suite --features cutest,ipopt-native --release -- ROSENBR HS35 HS71
//   cargo run --bin cutest_suite --features cutest,ipopt-native --release  # reads problem_list.txt

mod cutest_ffi;
mod cutest_problem;

use cutest_problem::CutestProblem;
use ripopt::{NlpProblem, SolverOptions, SolveStatus};
use serde::Serialize;
use std::ffi::CString;
use std::os::raw::c_void;
use std::path::Path;
use std::time::Instant;

// ---- Result type ----

#[derive(Serialize, serde::Deserialize)]
struct CutestResult {
    name: String,
    solver: String,
    n: usize,
    m: usize,
    status: String,
    objective: f64,
    x: Vec<f64>,
    constraint_violation: f64,
    iterations: usize,
    solve_time: f64,
}

// ---- Ipopt C API FFI (copied from hs_suite/run_ipopt_native.rs) ----

type IpoptProblem = *mut c_void;

extern "C" {
    fn CreateIpoptProblem(
        n: i32, x_l: *mut f64, x_u: *mut f64,
        m: i32, g_l: *mut f64, g_u: *mut f64,
        nele_jac: i32, nele_hess: i32, index_style: i32,
        eval_f: EvalFCB, eval_g: EvalGCB, eval_grad_f: EvalGradFCB,
        eval_jac_g: EvalJacGCB, eval_h: EvalHCB,
    ) -> IpoptProblem;
    fn FreeIpoptProblem(problem: IpoptProblem);
    fn AddIpoptStrOption(problem: IpoptProblem, keyword: *const i8, val: *const i8) -> bool;
    fn AddIpoptNumOption(problem: IpoptProblem, keyword: *const i8, val: f64) -> bool;
    fn AddIpoptIntOption(problem: IpoptProblem, keyword: *const i8, val: i32) -> bool;
    fn SetIntermediateCallback(problem: IpoptProblem, cb: IntermediateCB) -> bool;
    fn IpoptSolve(
        problem: IpoptProblem,
        x: *mut f64, g: *mut f64, obj_val: *mut f64,
        mult_g: *mut f64, mult_x_l: *mut f64, mult_x_u: *mut f64,
        user_data: *mut c_void,
    ) -> i32;
}

type EvalFCB = extern "C" fn(i32, *const f64, bool, *mut f64, *mut c_void) -> bool;
type EvalGradFCB = extern "C" fn(i32, *const f64, bool, *mut f64, *mut c_void) -> bool;
type EvalGCB = extern "C" fn(i32, *const f64, bool, i32, *mut f64, *mut c_void) -> bool;
type EvalJacGCB = extern "C" fn(i32, *const f64, bool, i32, i32, *mut i32, *mut i32, *mut f64, *mut c_void) -> bool;
type EvalHCB = extern "C" fn(i32, *const f64, bool, f64, i32, *const f64, bool, i32, *mut i32, *mut i32, *mut f64, *mut c_void) -> bool;
type IntermediateCB = extern "C" fn(i32, i32, f64, f64, f64, f64, f64, f64, f64, f64, i32, *mut c_void) -> bool;

struct ProblemWrapper<'a> {
    problem: &'a dyn NlpProblem,
    jac_rows: Vec<i32>,
    jac_cols: Vec<i32>,
    hess_rows: Vec<i32>,
    hess_cols: Vec<i32>,
    iterations: i32,
}

extern "C" fn eval_f_cb(
    n: i32, x: *const f64, _new_x: bool,
    obj_value: *mut f64, user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &*(user_data as *const ProblemWrapper);
        let x_slice = std::slice::from_raw_parts(x, n as usize);
        *obj_value = wrapper.problem.objective(x_slice);
        true
    }
}

extern "C" fn eval_grad_f_cb(
    n: i32, x: *const f64, _new_x: bool,
    grad_f: *mut f64, user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &*(user_data as *const ProblemWrapper);
        let x_slice = std::slice::from_raw_parts(x, n as usize);
        let grad_slice = std::slice::from_raw_parts_mut(grad_f, n as usize);
        wrapper.problem.gradient(x_slice, grad_slice);
        true
    }
}

extern "C" fn eval_g_cb(
    n: i32, x: *const f64, _new_x: bool,
    _m: i32, g: *mut f64, user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &*(user_data as *const ProblemWrapper);
        let x_slice = std::slice::from_raw_parts(x, n as usize);
        let m = wrapper.problem.num_constraints();
        if m > 0 {
            let g_slice = std::slice::from_raw_parts_mut(g, m);
            wrapper.problem.constraints(x_slice, g_slice);
        }
        true
    }
}

extern "C" fn eval_jac_g_cb(
    n: i32, x: *const f64, _new_x: bool,
    _m: i32, _nele_jac: i32,
    i_row: *mut i32, j_col: *mut i32, values: *mut f64,
    user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &*(user_data as *const ProblemWrapper);
        if values.is_null() {
            let nele = wrapper.jac_rows.len();
            let rows = std::slice::from_raw_parts_mut(i_row, nele);
            let cols = std::slice::from_raw_parts_mut(j_col, nele);
            for k in 0..nele {
                rows[k] = wrapper.jac_rows[k];
                cols[k] = wrapper.jac_cols[k];
            }
        } else {
            let x_slice = std::slice::from_raw_parts(x, n as usize);
            let nele = wrapper.jac_rows.len();
            let vals = std::slice::from_raw_parts_mut(values, nele);
            wrapper.problem.jacobian_values(x_slice, vals);
        }
        true
    }
}

extern "C" fn eval_h_cb(
    n: i32, x: *const f64, _new_x: bool,
    obj_factor: f64, _m: i32, lambda: *const f64, _new_lambda: bool,
    _nele_hess: i32,
    i_row: *mut i32, j_col: *mut i32, values: *mut f64,
    user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &*(user_data as *const ProblemWrapper);
        if values.is_null() {
            let nele = wrapper.hess_rows.len();
            let rows = std::slice::from_raw_parts_mut(i_row, nele);
            let cols = std::slice::from_raw_parts_mut(j_col, nele);
            for k in 0..nele {
                rows[k] = wrapper.hess_rows[k];
                cols[k] = wrapper.hess_cols[k];
            }
        } else {
            let x_slice = std::slice::from_raw_parts(x, n as usize);
            let m = wrapper.problem.num_constraints();
            let lambda_slice = if m > 0 {
                std::slice::from_raw_parts(lambda, m)
            } else {
                &[]
            };
            let nele = wrapper.hess_rows.len();
            let vals = std::slice::from_raw_parts_mut(values, nele);
            wrapper.problem.hessian_values(x_slice, obj_factor, lambda_slice, vals);
        }
        true
    }
}

extern "C" fn intermediate_cb(
    _alg_mod: i32, _iter_count: i32, _obj_value: f64,
    _inf_pr: f64, _inf_du: f64, _mu: f64,
    _d_norm: f64, _regularization_size: f64,
    _alpha_du: f64, _alpha_pr: f64, _ls_trials: i32,
    user_data: *mut c_void,
) -> bool {
    unsafe {
        let wrapper = &mut *(user_data as *mut ProblemWrapper);
        wrapper.iterations = _iter_count;
        true
    }
}

fn set_str_option(problem: IpoptProblem, key: &str, val: &str) {
    let k = CString::new(key).unwrap();
    let v = CString::new(val).unwrap();
    unsafe { AddIpoptStrOption(problem, k.as_ptr(), v.as_ptr()); }
}

fn set_num_option(problem: IpoptProblem, key: &str, val: f64) {
    let k = CString::new(key).unwrap();
    unsafe { AddIpoptNumOption(problem, k.as_ptr(), val); }
}

fn set_int_option(problem: IpoptProblem, key: &str, val: i32) {
    let k = CString::new(key).unwrap();
    unsafe { AddIpoptIntOption(problem, k.as_ptr(), val); }
}

// ---- Solve with Ipopt ----

struct IpoptResult {
    status: i32,
    objective: f64,
    x: Vec<f64>,
    constraint_violation: f64,
    iterations: i32,
    solve_time: f64, // Time for IpoptSolve only (excludes setup/teardown)
}

fn solve_with_ipopt(problem: &dyn NlpProblem) -> IpoptResult {
    let n = problem.num_variables();
    let m = problem.num_constraints();

    let mut x_l = vec![0.0; n];
    let mut x_u = vec![0.0; n];
    problem.bounds(&mut x_l, &mut x_u);

    let mut g_l = vec![0.0; m.max(1)];
    let mut g_u = vec![0.0; m.max(1)];
    if m > 0 {
        problem.constraint_bounds(&mut g_l, &mut g_u);
    }

    let (jac_rows_usize, jac_cols_usize) = problem.jacobian_structure();
    let nele_jac = jac_rows_usize.len();
    let jac_rows: Vec<i32> = jac_rows_usize.iter().map(|&r| r as i32).collect();
    let jac_cols: Vec<i32> = jac_cols_usize.iter().map(|&c| c as i32).collect();

    let (hess_rows_usize, hess_cols_usize) = problem.hessian_structure();
    let nele_hess = hess_rows_usize.len();
    let hess_rows: Vec<i32> = hess_rows_usize.iter().map(|&r| r as i32).collect();
    let hess_cols: Vec<i32> = hess_cols_usize.iter().map(|&c| c as i32).collect();

    let mut wrapper = ProblemWrapper {
        problem,
        jac_rows,
        jac_cols,
        hess_rows,
        hess_cols,
        iterations: 0,
    };

    unsafe {
        let ipopt_problem = CreateIpoptProblem(
            n as i32, x_l.as_mut_ptr(), x_u.as_mut_ptr(),
            m as i32, g_l.as_mut_ptr(), g_u.as_mut_ptr(),
            nele_jac as i32, nele_hess as i32,
            0, // C-style indexing
            eval_f_cb, eval_g_cb, eval_grad_f_cb,
            eval_jac_g_cb, eval_h_cb,
        );

        if ipopt_problem.is_null() {
            return IpoptResult {
                status: -199, objective: f64::NAN, x: vec![],
                constraint_violation: f64::NAN, iterations: 0,
                solve_time: 0.0,
            };
        }

        set_str_option(ipopt_problem, "sb", "yes");
        set_str_option(ipopt_problem, "mu_strategy", "adaptive");
        set_num_option(ipopt_problem, "tol", 1e-8);
        set_int_option(ipopt_problem, "max_iter", 3000);
        set_int_option(ipopt_problem, "print_level", 0);

        SetIntermediateCallback(ipopt_problem, intermediate_cb);

        let mut x = vec![0.0; n];
        problem.initial_point(&mut x);
        let mut g = vec![0.0; m.max(1)];
        let mut obj_val = 0.0;
        let mut mult_g = vec![0.0; m.max(1)];
        let mut mult_x_l = vec![0.0; n];
        let mut mult_x_u = vec![0.0; n];

        let user_data = &mut wrapper as *mut ProblemWrapper as *mut c_void;

        // Time only the IpoptSolve call (excludes setup/teardown)
        let t0 = Instant::now();
        let status = IpoptSolve(
            ipopt_problem,
            x.as_mut_ptr(), g.as_mut_ptr(), &mut obj_val,
            mult_g.as_mut_ptr(), mult_x_l.as_mut_ptr(), mult_x_u.as_mut_ptr(),
            user_data,
        );
        let solve_time = t0.elapsed().as_secs_f64();

        let iterations = wrapper.iterations;
        FreeIpoptProblem(ipopt_problem);

        // Compute constraint violation
        let cv = if m > 0 {
            let mut c = vec![0.0; m];
            problem.constraints(&x, &mut c);
            let mut g_l2 = vec![0.0; m];
            let mut g_u2 = vec![0.0; m];
            problem.constraint_bounds(&mut g_l2, &mut g_u2);
            compute_constraint_violation(&c, &g_l2, &g_u2)
        } else {
            0.0
        };

        IpoptResult {
            status,
            objective: obj_val,
            x,
            constraint_violation: cv,
            iterations,
            solve_time,
        }
    }
}

fn ipopt_status_to_string(status: i32) -> String {
    match status {
        0 => "Optimal".to_string(),
        1 => "Acceptable".to_string(),
        2 => "Infeasible".to_string(),
        -1 => "MaxIterations".to_string(),
        -2 => "RestorationFailed".to_string(),
        -3 => "ErrorInStepComputation".to_string(),
        -13 => "InvalidNumberDetected".to_string(),
        other => format!("IpoptStatus({})", other),
    }
}

fn ripopt_status_to_string(status: SolveStatus) -> String {
    match status {
        SolveStatus::Optimal => "Optimal".to_string(),
        SolveStatus::Acceptable => "Acceptable".to_string(),
        SolveStatus::Infeasible => "Infeasible".to_string(),
        SolveStatus::MaxIterations => "MaxIterations".to_string(),
        SolveStatus::NumericalError => "NumericalError".to_string(),
        SolveStatus::Unbounded => "Unbounded".to_string(),
        SolveStatus::RestorationFailed => "RestorationFailed".to_string(),
        SolveStatus::InternalError => "InternalError".to_string(),
        SolveStatus::LocalInfeasibility => "LocalInfeasibility".to_string(),
    }
}

fn compute_constraint_violation(c: &[f64], g_l: &[f64], g_u: &[f64]) -> f64 {
    let mut max_viol = 0.0f64;
    for i in 0..c.len() {
        if c[i] < g_l[i] {
            max_viol = max_viol.max(g_l[i] - c[i]);
        }
        if c[i] > g_u[i] {
            max_viol = max_viol.max(c[i] - g_u[i]);
        }
    }
    max_viol
}

// ---- Main ----

fn get_problem_list_from_args_or_file(suite_dir: &Path) -> Vec<String> {
    // Skip args[0] (binary), args[1..] are problem names (unless --single mode handled above)
    let args: Vec<String> = std::env::args().skip(1).collect();
    if !args.is_empty() {
        return args;
    }

    // Read from problem_list.txt
    let list_path = suite_dir.join("problem_list.txt");
    if list_path.exists() {
        let contents = std::fs::read_to_string(&list_path)
            .expect("Failed to read problem_list.txt");
        contents
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect()
    } else {
        eprintln!("No problems specified and problem_list.txt not found.");
        eprintln!("Usage: cutest_suite PROBLEM1 PROBLEM2 ...");
        std::process::exit(1);
    }
}

/// Solve a single problem in subprocess mode. Outputs JSON lines to stdout.
fn run_single_problem(name: &str) {
    let n_timing_runs: usize = std::env::var("RIPOPT_TIMING_RUNS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);

    let suite_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("cutest_suite");
    let problems_dir = suite_dir.join("problems");

    let lib_path = problems_dir.join(format!("lib{}.dylib", name));
    let outsdif_path = problems_dir.join(format!("{}_OUTSDIF.d", name));

    let lib_str = lib_path.to_str().unwrap();
    let outsdif_str = outsdif_path.to_str().unwrap();

    let problem = match CutestProblem::load(name, lib_str, outsdif_str) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  SKIP {} (load failed: {})", name, e);
            std::process::exit(1);
        }
    };

    let print_level: u8 = std::env::var("RIPOPT_PRINT_LEVEL")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let options = SolverOptions {
        tol: 1e-8,
        max_iter: 3000,
        print_level,
        mu_strategy_adaptive: true,
        max_wall_time: 30.0,
        ..SolverOptions::default()
    };

    // Solve with ripopt
    let mut best_ripopt_time = f64::MAX;
    let mut ripopt_result = None;
    for run in 0..n_timing_runs {
        let t0 = Instant::now();
        let r = ripopt::solve(&problem, &options);
        let elapsed = t0.elapsed().as_secs_f64();
        if elapsed < best_ripopt_time {
            best_ripopt_time = elapsed;
        }
        if run == 0 {
            ripopt_result = Some(r);
        }
    }
    let ripopt_result = ripopt_result.unwrap();
    let ripopt_cv = if problem.m > 0 {
        let mut g_l = vec![0.0; problem.m];
        let mut g_u = vec![0.0; problem.m];
        problem.constraint_bounds(&mut g_l, &mut g_u);
        compute_constraint_violation(&ripopt_result.constraint_values, &g_l, &g_u)
    } else {
        0.0
    };

    let r1 = CutestResult {
        name: name.to_string(),
        solver: "ripopt".to_string(),
        n: problem.n,
        m: problem.m,
        status: ripopt_status_to_string(ripopt_result.status),
        objective: if ripopt_result.objective.is_finite() { ripopt_result.objective } else { 0.0 },
        x: ripopt_result.x.iter().map(|v| if v.is_finite() { *v } else { 0.0 }).collect(),
        constraint_violation: if ripopt_cv.is_finite() { ripopt_cv } else { 0.0 },
        iterations: ripopt_result.iterations,
        solve_time: best_ripopt_time,
    };
    println!("{}", serde_json::to_string(&r1).unwrap());

    // Solve with Ipopt — use internal solve_time (IpoptSolve only, excludes setup/teardown)
    let mut best_ipopt_time = f64::MAX;
    let mut ipopt_result = None;
    for run in 0..n_timing_runs {
        let r = solve_with_ipopt(&problem);
        if r.solve_time < best_ipopt_time {
            best_ipopt_time = r.solve_time;
        }
        if run == 0 {
            ipopt_result = Some(r);
        }
    }
    let ipopt_result = ipopt_result.unwrap();

    let r2 = CutestResult {
        name: name.to_string(),
        solver: "ipopt".to_string(),
        n: problem.n,
        m: problem.m,
        status: ipopt_status_to_string(ipopt_result.status),
        objective: ipopt_result.objective,
        x: ipopt_result.x.clone(),
        constraint_violation: ipopt_result.constraint_violation,
        iterations: ipopt_result.iterations as usize,
        solve_time: best_ipopt_time,
    };
    println!("{}", serde_json::to_string(&r2).unwrap());

    eprintln!(
        "    ripopt: {} (obj={:.6e}, {:.1}ms)  ipopt: {} (obj={:.6e}, {:.1}ms)",
        r1.status, r1.objective, best_ripopt_time * 1000.0,
        r2.status, r2.objective, best_ipopt_time * 1000.0,
    );

    problem.cleanup();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    // Subprocess mode: --single PROBLEM
    if args.len() == 3 && args[1] == "--single" {
        run_single_problem(&args[2]);
        return;
    }

    let n_timing_runs: usize = std::env::var("RIPOPT_TIMING_RUNS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);

    let max_n: usize = std::env::var("CUTEST_MAX_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(100);

    let timeout_secs: u64 = std::env::var("CUTEST_TIMEOUT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(60);

    // Find the cutest_suite directory (next to the binary source)
    let suite_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("cutest_suite");
    let problems_dir = suite_dir.join("problems");

    let problem_names = get_problem_list_from_args_or_file(&suite_dir);
    eprintln!(
        "CUTEst benchmark: {} problems, {} timing runs, max_n={}, timeout={}s",
        problem_names.len(),
        n_timing_runs,
        max_n,
        timeout_secs,
    );

    let self_exe = std::env::current_exe().expect("cannot find self executable");

    let mut all_results: Vec<CutestResult> = Vec::new();

    for name in &problem_names {
        let lib_path = problems_dir.join(format!("lib{}.dylib", name));
        let outsdif_path = problems_dir.join(format!("{}_OUTSDIF.d", name));

        if !lib_path.exists() || !outsdif_path.exists() {
            eprintln!("  SKIP {} (not prepared — run prepare.sh first)", name);
            continue;
        }

        // Check dimensions by loading briefly
        let lib_str = lib_path.to_str().unwrap();
        let outsdif_str = outsdif_path.to_str().unwrap();
        let problem = match CutestProblem::load(name, lib_str, outsdif_str) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("  SKIP {} (load failed: {})", name, e);
                continue;
            }
        };
        let n = problem.n;
        let m = problem.m;
        problem.cleanup();

        if n > max_n {
            eprintln!("  SKIP {} (n={} > max_n={})", name, n, max_n);
            continue;
        }

        eprint!("  {} (n={}, m={}) ... ", name, n, m);

        // Run as subprocess with timeout
        let output = std::process::Command::new("timeout")
            .arg(format!("{}s", timeout_secs))
            .arg(&self_exe)
            .arg("--single")
            .arg(name)
            .env("RIPOPT_TIMING_RUNS", n_timing_runs.to_string())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .output();

        match output {
            Ok(out) => {
                let stdout = String::from_utf8_lossy(&out.stdout);
                let stderr = String::from_utf8_lossy(&out.stderr);

                if !out.status.success() && stdout.is_empty() {
                    // Timeout or crash
                    eprintln!("TIMEOUT/CRASH (exit={})", out.status);
                    all_results.push(CutestResult {
                        name: name.clone(), solver: "ripopt".to_string(),
                        n, m, status: "Timeout".to_string(),
                        objective: f64::NAN, x: vec![],
                        constraint_violation: f64::NAN, iterations: 0,
                        solve_time: timeout_secs as f64,
                    });
                    all_results.push(CutestResult {
                        name: name.clone(), solver: "ipopt".to_string(),
                        n, m, status: "Timeout".to_string(),
                        objective: f64::NAN, x: vec![],
                        constraint_violation: f64::NAN, iterations: 0,
                        solve_time: timeout_secs as f64,
                    });
                    continue;
                }

                // Parse JSON lines from stdout
                let mut parsed_any = false;
                for line in stdout.lines() {
                    let line = line.trim();
                    if line.is_empty() {
                        continue;
                    }
                    if let Ok(result) = serde_json::from_str::<CutestResult>(line) {
                        all_results.push(result);
                        parsed_any = true;
                    }
                }

                // Print the solver's stderr summary line
                for line in stderr.lines() {
                    if line.trim().starts_with("ripopt:") || line.trim().starts_with("ipopt:") {
                        eprintln!("{}", line.trim());
                    }
                }

                if !parsed_any {
                    eprintln!("PARSE ERROR");
                }
            }
            Err(e) => {
                eprintln!("SPAWN ERROR: {}", e);
            }
        }
    }

    // Write JSON to results.json (always saved) and also to stdout
    let json = serde_json::to_string_pretty(&all_results).unwrap();
    let results_path = suite_dir.join("results.json");
    if let Err(e) = std::fs::write(&results_path, &json) {
        eprintln!("WARNING: Failed to write {}: {}", results_path.display(), e);
    } else {
        eprintln!("Results written to {}", results_path.display());
    }
    println!("{}", json);

    // Summary to stderr
    let ripopt_results: Vec<_> = all_results.iter().filter(|r| r.solver == "ripopt").collect();
    let ipopt_results: Vec<_> = all_results.iter().filter(|r| r.solver == "ipopt").collect();
    let n_problems = ripopt_results.len();
    let ripopt_solved = ripopt_results
        .iter()
        .filter(|r| r.status == "Optimal" || r.status == "Acceptable")
        .count();
    let ipopt_solved = ipopt_results
        .iter()
        .filter(|r| r.status == "Optimal" || r.status == "Acceptable")
        .count();
    eprintln!("\nSummary: {} problems", n_problems);
    eprintln!("  ripopt solved: {}/{}", ripopt_solved, n_problems);
    eprintln!("  ipopt  solved: {}/{}", ipopt_solved, n_problems);
}
