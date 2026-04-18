//! Per-iteration TSV tracing for direction-diff harness against Ipopt.
//!
//! Activate by setting the `RIP_TRACE_TSV` environment variable to a file
//! path. The IPM main loop then writes one TSV row per iteration with
//! algorithmic intermediates (mu, step norms, Mehrotra sigma, inertia,
//! etc.). The schema is chosen to mirror what Ipopt's IntermediateCallback
//! + print_level=8 trace provide so the two can be joined on `iter` for a
//! side-by-side comparison.
//!
//! Doing nothing when the env var is absent keeps this zero-cost for normal
//! benchmark runs.

use std::fs::OpenOptions;
use std::io::Write;
use std::sync::Mutex;
use std::sync::OnceLock;

pub const COLUMNS: &[&str] = &[
    "iter",
    "obj",
    "inf_pr",
    "inf_du",
    "compl",
    "mu",
    "alpha_pr",
    "alpha_du",
    "alpha_affP",
    "alpha_affD",
    "mu_aff",
    "sigma",
    "mu_pc",
    "delta_w",
    "delta_c",
    "dx_inf",
    "dzl_inf",
    "dzu_inf",
    "mcc_iters",
    "ls",
    "accepted",
];

struct TraceSink {
    file: std::fs::File,
}

static SINK: OnceLock<Option<Mutex<TraceSink>>> = OnceLock::new();

fn sink() -> Option<&'static Mutex<TraceSink>> {
    SINK.get_or_init(|| {
        let path = std::env::var("RIP_TRACE_TSV").ok()?;
        let mut file = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&path)
            .ok()?;
        writeln!(file, "{}", COLUMNS.join("\t")).ok()?;
        Some(Mutex::new(TraceSink { file }))
    })
    .as_ref()
}

/// Per-iteration row; every field is f64 to keep the schema flat. Use NaN
/// when a column doesn't apply (e.g. Mehrotra σ when PC is skipped).
#[derive(Clone, Copy)]
pub struct TraceRow {
    pub iter: usize,
    pub obj: f64,
    pub inf_pr: f64,
    pub inf_du: f64,
    pub compl: f64,
    pub mu: f64,
    pub alpha_pr: f64,
    pub alpha_du: f64,
    pub alpha_aff_p: f64,
    pub alpha_aff_d: f64,
    pub mu_aff: f64,
    pub sigma: f64,
    pub mu_pc: f64,
    pub delta_w: f64,
    pub delta_c: f64,
    pub dx_inf: f64,
    pub dzl_inf: f64,
    pub dzu_inf: f64,
    pub mcc_iters: u32,
    pub ls: u32,
    pub accepted: bool,
}

impl TraceRow {
    pub fn blank(iter: usize) -> Self {
        let nan = f64::NAN;
        Self {
            iter,
            obj: nan, inf_pr: nan, inf_du: nan, compl: nan, mu: nan,
            alpha_pr: nan, alpha_du: nan,
            alpha_aff_p: nan, alpha_aff_d: nan, mu_aff: nan, sigma: nan, mu_pc: nan,
            delta_w: nan, delta_c: nan,
            dx_inf: nan, dzl_inf: nan, dzu_inf: nan,
            mcc_iters: 0, ls: 0, accepted: false,
        }
    }
}

/// Emit one TSV row if `RIP_TRACE_TSV` is set; no-op otherwise.
pub fn emit(row: &TraceRow) {
    let Some(sink) = sink() else { return };
    let Ok(mut guard) = sink.lock() else { return };
    let _ = writeln!(
        guard.file,
        "{}\t{:.17e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{}\t{}\t{}",
        row.iter,
        row.obj, row.inf_pr, row.inf_du, row.compl, row.mu,
        row.alpha_pr, row.alpha_du,
        row.alpha_aff_p, row.alpha_aff_d, row.mu_aff, row.sigma, row.mu_pc,
        row.delta_w, row.delta_c,
        row.dx_inf, row.dzl_inf, row.dzu_inf,
        row.mcc_iters, row.ls, if row.accepted { 1 } else { 0 },
    );
    let _ = guard.file.flush();
}
