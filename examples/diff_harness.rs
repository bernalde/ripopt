// Direction-diff harness: runs the same 3-variable eval-failure-prone problem
// through both ripopt and Ipopt, emitting per-iteration TSV files so the two
// trajectories can be joined on `iter` and diffed column-by-column.
//
// The target problem is a synthesis of what @prehner's feos-campd CAMD test
// exhibits: small variable count, nonlinear non-convex constraints, a
// log-domain region where the objective/constraints can produce NaN if the
// line search steps outside the feasible shell. Starts from an infeasible
// interior point so both solvers have to make real progress, and both
// trajectories can be compared iteration-by-iteration.
//
// Run with:
//
//   RIP_TRACE_TSV=/tmp/ripopt_trace.tsv \
//     cargo run --release --features ipopt-native --example diff_harness
//
// Produces:
//   - /tmp/ripopt_trace.tsv  (ripopt per-iter TSV, via src/trace.rs)
//   - /tmp/ipopt_trace.tsv   (Ipopt per-iter TSV, via IntermediateCallback)
// and prints a summary header + a "column difference" sketch based on the
// first 20 iterations of each.

use ripopt::{NlpProblem, SolverOptions};
use std::cell::RefCell;
use std::ffi::CString;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::os::raw::c_void;
use std::time::Instant;

// ---------------------------------------------------------------------------
// Ipopt C API FFI — identical to compare_large_scale.rs. Duplicated rather
// than factored out because examples can't share modules.
// ---------------------------------------------------------------------------

type IpoptProblemPtr = *mut c_void;

extern "C" {
    fn CreateIpoptProblem(
        n: i32, x_l: *mut f64, x_u: *mut f64,
        m: i32, g_l: *mut f64, g_u: *mut f64,
        nele_jac: i32, nele_hess: i32, index_style: i32,
        eval_f: EvalFCB, eval_g: EvalGCB, eval_grad_f: EvalGradFCB,
        eval_jac_g: EvalJacGCB, eval_h: EvalHCB,
    ) -> IpoptProblemPtr;

    fn FreeIpoptProblem(problem: IpoptProblemPtr);
    fn AddIpoptStrOption(problem: IpoptProblemPtr, keyword: *const i8, val: *const i8) -> bool;
    fn AddIpoptNumOption(problem: IpoptProblemPtr, keyword: *const i8, val: f64) -> bool;
    fn AddIpoptIntOption(problem: IpoptProblemPtr, keyword: *const i8, val: i32) -> bool;
    fn SetIntermediateCallback(problem: IpoptProblemPtr, cb: IntermediateCB) -> bool;
    fn IpoptSolve(
        problem: IpoptProblemPtr,
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

// ---------------------------------------------------------------------------
// The problem.
//
//   min  (x1 - 2)^2 + (x2 - 1)^2 + (x3 + 0.5)^2
//   s.t. log(x1) + log(x2) - x3 >= 0         (c0: eval-failure if x1<=0 or x2<=0)
//        x1^2 + x2^2 <= 2                     (c1)
//        x1 + x2 >= 0.5                       (c2)
//        0.01 <= x1, x2 <= 10
//        -5   <= x3  <= 5
//
// Quadratic objective bounds the problem. Start from x = (0.5, 0.5, -0.5):
// constraint c0 is infeasible (log(0.5)+log(0.5) = -1.386 < -0.5), c1 and
// c2 are satisfied. Both solvers must move to reach the optimum. The
// log-domain eval-failure on x_i <= 0 mirrors the phase-equilibria pattern
// in the feos-campd CAMD reproducer.
// ---------------------------------------------------------------------------

struct LogBand;

impl NlpProblem for LogBand {
    fn num_variables(&self) -> usize { 3 }
    fn num_constraints(&self) -> usize { 3 }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l[0] = 0.01; x_u[0] = 10.0;
        x_l[1] = 0.01; x_u[1] = 10.0;
        x_l[2] = -5.0; x_u[2] = 5.0;
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        // c0: log(x1) + log(x2) - x3 >= 0   => g_l=0, g_u=+inf
        g_l[0] = 0.0; g_u[0] = f64::INFINITY;
        // c1: x1^2 + x2^2 <= 2             => g_l=-inf, g_u=2
        g_l[1] = f64::NEG_INFINITY; g_u[1] = 2.0;
        // c2: x1 + x2 >= 0.5               => g_l=0.5, g_u=+inf
        g_l[2] = 0.5; g_u[2] = f64::INFINITY;
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0[0] = 0.5; x0[1] = 0.5; x0[2] = -0.5;
    }

    fn objective(&self, x: &[f64], _new_x: bool, obj: &mut f64) -> bool {
        if x[0] <= 0.0 || x[1] <= 0.0 { return false; }
        let a = x[0] - 2.0;
        let b = x[1] - 1.0;
        let c = x[2] + 0.5;
        *obj = a*a + b*b + c*c;
        true
    }

    fn gradient(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
        if x[0] <= 0.0 || x[1] <= 0.0 { return false; }
        g[0] = 2.0 * (x[0] - 2.0);
        g[1] = 2.0 * (x[1] - 1.0);
        g[2] = 2.0 * (x[2] + 0.5);
        true
    }

    fn constraints(&self, x: &[f64], _new_x: bool, g: &mut [f64]) -> bool {
        if x[0] <= 0.0 || x[1] <= 0.0 { return false; }
        g[0] = x[0].ln() + x[1].ln() - x[2];
        g[1] = x[0]*x[0] + x[1]*x[1];
        g[2] = x[0] + x[1];
        true
    }

    // Jacobian: c0 = log(x1)+log(x2)-x3  -> [1/x1, 1/x2, -1]
    //           c1 = x1^2+x2^2           -> [2x1, 2x2, 0]
    //           c2 = x1+x2                -> [1, 1, 0]
    // Sparsity: c0 all 3 vars; c1 first 2; c2 first 2. Total nnz = 7.
    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (
            vec![0, 0, 0, 1, 1, 2, 2],
            vec![0, 1, 2, 0, 1, 0, 1],
        )
    }

    fn jacobian_values(&self, x: &[f64], _new_x: bool, v: &mut [f64]) -> bool {
        if x[0] <= 0.0 || x[1] <= 0.0 { return false; }
        v[0] = 1.0 / x[0];
        v[1] = 1.0 / x[1];
        v[2] = -1.0;
        v[3] = 2.0 * x[0];
        v[4] = 2.0 * x[1];
        v[5] = 1.0;
        v[6] = 1.0;
        true
    }

    // Hessian of Lagrangian = σ·∇²f + Σ λ_i ∇²g_i.
    // ∇²f = 0. ∇²g0 = diag(-1/x1^2, -1/x2^2, 0). ∇²g1 = diag(2, 2, 0).
    // ∇²g2 = 0. Lower triangle → entries (0,0), (1,1) only. Total nnz = 2.
    // ∇²f = diag(2, 2, 2). ∇²g0 = diag(-1/x1², -1/x2², 0). ∇²g1 = diag(2,2,0).
    // ∇²g2 = 0. Lower triangle → (0,0), (1,1), (2,2).
    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (vec![0, 1, 2], vec![0, 1, 2])
    }

    fn hessian_values(
        &self, x: &[f64], _new_x: bool, obj_factor: f64,
        lambda: &[f64], v: &mut [f64],
    ) -> bool {
        if x[0] <= 0.0 || x[1] <= 0.0 { return false; }
        v[0] = 2.0 * obj_factor - lambda[0] / (x[0]*x[0]) + 2.0 * lambda[1];
        v[1] = 2.0 * obj_factor - lambda[0] / (x[1]*x[1]) + 2.0 * lambda[1];
        v[2] = 2.0 * obj_factor;
        true
    }
}

// ---------------------------------------------------------------------------
// Ipopt-side TSV writer. Mirrors the schema in src/trace.rs. Some columns
// aren't exposed by IntermediateCallback and stay NaN (sigma, mu_aff, etc.).
// ---------------------------------------------------------------------------

thread_local! {
    static IPOPT_TSV: RefCell<Option<File>> = const { RefCell::new(None) };
}

fn ipopt_tsv_open(path: &str) {
    let mut f = OpenOptions::new()
        .create(true).truncate(true).write(true).open(path).unwrap();
    // Same schema as src/trace.rs
    writeln!(f, "iter\tobj\tinf_pr\tinf_du\tcompl\tmu\talpha_pr\talpha_du\talpha_affP\talpha_affD\tmu_aff\tsigma\tmu_pc\tdelta_w\tdelta_c\tdx_inf\tdzl_inf\tdzu_inf\tmcc_iters\tls\taccepted").unwrap();
    IPOPT_TSV.with(|c| *c.borrow_mut() = Some(f));
}

struct IpoptWrapper<'a> {
    problem: &'a dyn NlpProblem,
    jac_rows: Vec<i32>,
    jac_cols: Vec<i32>,
    hess_rows: Vec<i32>,
    hess_cols: Vec<i32>,
}

extern "C" fn eval_f_cb(n: i32, x: *const f64, _: bool, obj: *mut f64, ud: *mut c_void) -> bool {
    unsafe {
        let w = &*(ud as *const IpoptWrapper);
        let xs = std::slice::from_raw_parts(x, n as usize);
        w.problem.objective(xs, true, &mut *obj)
    }
}
extern "C" fn eval_grad_f_cb(n: i32, x: *const f64, _: bool, g: *mut f64, ud: *mut c_void) -> bool {
    unsafe {
        let w = &*(ud as *const IpoptWrapper);
        let xs = std::slice::from_raw_parts(x, n as usize);
        let gs = std::slice::from_raw_parts_mut(g, n as usize);
        w.problem.gradient(xs, true, gs)
    }
}
extern "C" fn eval_g_cb(n: i32, x: *const f64, _: bool, _m: i32, g: *mut f64, ud: *mut c_void) -> bool {
    unsafe {
        let w = &*(ud as *const IpoptWrapper);
        let m = w.problem.num_constraints();
        let xs = std::slice::from_raw_parts(x, n as usize);
        let gs = std::slice::from_raw_parts_mut(g, m);
        w.problem.constraints(xs, true, gs)
    }
}
extern "C" fn eval_jac_g_cb(
    n: i32, x: *const f64, _: bool, _m: i32, _nj: i32,
    ir: *mut i32, jc: *mut i32, vals: *mut f64, ud: *mut c_void,
) -> bool {
    unsafe {
        let w = &*(ud as *const IpoptWrapper);
        if vals.is_null() {
            let nele = w.jac_rows.len();
            std::slice::from_raw_parts_mut(ir, nele).copy_from_slice(&w.jac_rows);
            std::slice::from_raw_parts_mut(jc, nele).copy_from_slice(&w.jac_cols);
        } else {
            let xs = std::slice::from_raw_parts(x, n as usize);
            let vs = std::slice::from_raw_parts_mut(vals, w.jac_rows.len());
            if !w.problem.jacobian_values(xs, true, vs) { return false; }
        }
        true
    }
}
extern "C" fn eval_h_cb(
    n: i32, x: *const f64, _: bool, obj_factor: f64,
    _m: i32, lambda: *const f64, _: bool, _nh: i32,
    ir: *mut i32, jc: *mut i32, vals: *mut f64, ud: *mut c_void,
) -> bool {
    unsafe {
        let w = &*(ud as *const IpoptWrapper);
        if vals.is_null() {
            let nele = w.hess_rows.len();
            std::slice::from_raw_parts_mut(ir, nele).copy_from_slice(&w.hess_rows);
            std::slice::from_raw_parts_mut(jc, nele).copy_from_slice(&w.hess_cols);
        } else {
            let xs = std::slice::from_raw_parts(x, n as usize);
            let m = w.problem.num_constraints();
            let ls = std::slice::from_raw_parts(lambda, m);
            let vs = std::slice::from_raw_parts_mut(vals, w.hess_rows.len());
            if !w.problem.hessian_values(xs, true, obj_factor, ls, vs) { return false; }
        }
        true
    }
}

extern "C" fn intermediate_cb(
    _alg_mod: i32, iter: i32, obj: f64, inf_pr: f64, inf_du: f64,
    mu: f64, d_norm: f64, reg_size: f64, alpha_du: f64, alpha_pr: f64,
    ls_trials: i32, _ud: *mut c_void,
) -> bool {
    IPOPT_TSV.with(|c| {
        if let Some(ref mut f) = *c.borrow_mut() {
            // Columns that map from IntermediateCB:
            //   iter, obj, inf_pr, inf_du, mu, alpha_pr, alpha_du, delta_w (reg_size), dx_inf (d_norm), ls
            // compl, alpha_affP/D, mu_aff, sigma, mu_pc, delta_c, dzl/dzu_inf, mcc_iters: NaN
            let nan = f64::NAN;
            writeln!(f,
                "{}\t{:.17e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{}\t{}\t{}",
                iter, obj, inf_pr, inf_du, nan, mu, alpha_pr, alpha_du,
                nan, nan, nan, nan, nan, reg_size, nan, d_norm, nan, nan,
                0, ls_trials, 1
            ).ok();
            let _ = f.flush();
        }
    });
    true
}

fn set_str(p: IpoptProblemPtr, k: &str, v: &str) {
    let ks = CString::new(k).unwrap();
    let vs = CString::new(v).unwrap();
    unsafe { AddIpoptStrOption(p, ks.as_ptr(), vs.as_ptr()); }
}
fn set_num(p: IpoptProblemPtr, k: &str, v: f64) {
    let ks = CString::new(k).unwrap();
    unsafe { AddIpoptNumOption(p, ks.as_ptr(), v); }
}
fn set_int(p: IpoptProblemPtr, k: &str, v: i32) {
    let ks = CString::new(k).unwrap();
    unsafe { AddIpoptIntOption(p, ks.as_ptr(), v); }
}

fn solve_ipopt(problem: &dyn NlpProblem) -> (String, f64, i32, f64) {
    let n = problem.num_variables();
    let m = problem.num_constraints();
    let mut x_l = vec![0.0; n];
    let mut x_u = vec![0.0; n];
    problem.bounds(&mut x_l, &mut x_u);
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    let (jr, jc) = problem.jacobian_structure();
    let (hr, hc) = problem.hessian_structure();
    let mut w = IpoptWrapper {
        problem,
        jac_rows: jr.iter().map(|&r| r as i32).collect(),
        jac_cols: jc.iter().map(|&c| c as i32).collect(),
        hess_rows: hr.iter().map(|&r| r as i32).collect(),
        hess_cols: hc.iter().map(|&c| c as i32).collect(),
    };
    unsafe {
        let ip = CreateIpoptProblem(
            n as i32, x_l.as_mut_ptr(), x_u.as_mut_ptr(),
            m as i32, g_l.as_mut_ptr(), g_u.as_mut_ptr(),
            w.jac_rows.len() as i32, w.hess_rows.len() as i32,
            0, eval_f_cb, eval_g_cb, eval_grad_f_cb, eval_jac_g_cb, eval_h_cb,
        );
        set_str(ip, "sb", "yes");
        set_str(ip, "mu_strategy", "adaptive");
        set_num(ip, "tol", 1e-8);
        set_int(ip, "max_iter", 500);
        set_int(ip, "print_level", 0);
        SetIntermediateCallback(ip, intermediate_cb);
        let mut x = vec![0.0; n];
        problem.initial_point(&mut x);
        let mut g = vec![0.0; m];
        let mut obj = 0.0;
        let mut mg = vec![0.0; m];
        let mut ml = vec![0.0; n];
        let mut mu = vec![0.0; n];
        let ud = &mut w as *mut IpoptWrapper as *mut c_void;
        let t0 = Instant::now();
        let status = IpoptSolve(ip, x.as_mut_ptr(), g.as_mut_ptr(), &mut obj,
            mg.as_mut_ptr(), ml.as_mut_ptr(), mu.as_mut_ptr(), ud);
        let elapsed = t0.elapsed().as_secs_f64();
        FreeIpoptProblem(ip);
        let status_str = match status {
            0 => "Optimal".to_string(),
            1 => "Acceptable".to_string(),
            2 => "Infeasible".to_string(),
            -1 => "MaxIter".to_string(),
            -2 => "RestorationFailed".to_string(),
            other => format!("Status({})", other),
        };
        // Use line count (minus header) in the TSV as iteration count
        let iters = count_tsv_rows("/tmp/ipopt_trace.tsv");
        (status_str, obj, iters, elapsed)
    }
}

fn count_tsv_rows(path: &str) -> i32 {
    let Ok(f) = File::open(path) else { return 0 };
    let r = BufReader::new(f);
    (r.lines().count() as i32 - 1).max(0)
}

fn solve_ripopt(problem: &LogBand) -> (String, f64, i32, f64) {
    let opts = SolverOptions {
        tol: 1e-8,
        max_iter: 500,
        max_wall_time: 60.0,
        print_level: 0,
        ..SolverOptions::default()
    };
    let t0 = Instant::now();
    let res = ripopt::solve(problem, &opts);
    let t = t0.elapsed().as_secs_f64();
    (format!("{:?}", res.status), res.objective, res.iterations as i32, t)
}

fn main() {
    // Ensure the Ipopt-side trace is captured to a deterministic path.
    let ipopt_tsv = std::env::var("IPOPT_TRACE_TSV").unwrap_or_else(|_| "/tmp/ipopt_trace.tsv".into());
    let ripopt_tsv = std::env::var("RIP_TRACE_TSV").unwrap_or_else(|_| {
        // Set it for the ripopt solve below; src/trace.rs reads the env at
        // first-call time, so setting here is effective for this process.
        std::env::set_var("RIP_TRACE_TSV", "/tmp/ripopt_trace.tsv");
        "/tmp/ripopt_trace.tsv".into()
    });

    ipopt_tsv_open(&ipopt_tsv);

    let problem = LogBand;

    println!("Direction-diff harness: 3-var log-domain problem");
    println!("  ripopt TSV: {}", ripopt_tsv);
    println!("  ipopt  TSV: {}", ipopt_tsv);
    println!();

    let (rstat, robj, riters, rtime) = solve_ripopt(&problem);
    let (istat, iobj, iiters, itime) = solve_ipopt(&problem);

    println!("{:<10} {:<20} {:>16} {:>10} {:>10}", "solver", "status", "obj", "iter", "time_s");
    println!("{:<10} {:<20} {:>16.10e} {:>10} {:>10.3}", "ripopt", rstat, robj, riters, rtime);
    println!("{:<10} {:<20} {:>16.10e} {:>10} {:>10.3}", "ipopt",  istat, iobj, iiters, itime);
    println!();
    println!("Inspect the TSVs side-by-side:");
    println!("  paste <(cut -f1-8 {}) <(cut -f1-8 {}) | head -30", ripopt_tsv, ipopt_tsv);
}
