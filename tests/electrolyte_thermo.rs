//! Electrolyte thermodynamics test suite for ripopt.
//!
//! Run with: cargo test electrolyte -- --nocapture

#[path = "common/electrolyte_problems.rs"]
mod problems;
use problems::*;

use ripopt::{NlpProblem, SolveStatus, SolverOptions};
use std::time::Instant;

fn default_options() -> SolverOptions {
    SolverOptions {
        tol: 1e-6,
        max_iter: 3000,
        print_level: 0,
        ..SolverOptions::default()
    }
}

/// Compute max constraint violation.
fn max_cv(problem: &dyn NlpProblem, g: &[f64]) -> f64 {
    let m = problem.num_constraints();
    if m == 0 {
        return 0.0;
    }
    let mut g_l = vec![0.0; m];
    let mut g_u = vec![0.0; m];
    problem.constraint_bounds(&mut g_l, &mut g_u);
    let mut cv = 0.0_f64;
    for i in 0..m {
        cv = cv.max((g_l[i] - g[i]).max(0.0)).max((g[i] - g_u[i]).max(0.0));
    }
    cv
}

macro_rules! electrolyte_test {
    ($name:ident, $problem:expr, $check:expr) => {
        #[test]
        fn $name() {
            let problem = $problem;
            let options = default_options();
            let start = Instant::now();
            let result = ripopt::solve(&problem, &options);
            let elapsed = start.elapsed();
            let cv = max_cv(&problem, &result.constraint_values);
            eprintln!(
                "{}: status={:?}, obj={:.6e}, cv={:.2e}, iters={}, time={:.3}s",
                stringify!($name), result.status, result.objective, cv,
                result.iterations, elapsed.as_secs_f64()
            );
            assert!(
                result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
                "Expected Optimal/Acceptable, got {:?}", result.status
            );
            if problem.num_constraints() > 0 {
                assert!(cv < 1e-4, "Constraint violation {:.2e} too large", cv);
            }
            // Problem-specific checks
            #[allow(clippy::redundant_closure_call)]
            ($check)(&result);
        }
    };
}

// ---------------------------------------------------------------------------
// Category 1: Speciation / Chemical Equilibrium
// ---------------------------------------------------------------------------

electrolyte_test!(electrolyte_01_water_autoionization, WaterAutoionization, |result: &ripopt::SolveResult| {
    // Solution: m_H ~ 1.005e-7, x ~ ln(1.005e-7) ~ -16.1
    let m = result.x[0].exp();
    assert!(m > 5e-8 && m < 5e-7, "m_H={:.3e} out of range", m);
});

electrolyte_test!(electrolyte_02_co2_water, Co2WaterSpeciation, |result: &ripopt::SolveResult| {
    // pH ~ 5.65, dominant species H2CO3 ~ 9.95e-4
    let m_h = result.x[3].exp();
    let ph = -(m_h.log10());
    assert!(ph > 4.0 && ph < 7.0, "pH={:.2} out of range", ph);
    let m_h2co3 = result.x[0].exp();
    assert!(m_h2co3 > 1e-4 && m_h2co3 < 2e-3, "m_H2CO3={:.3e} unexpected", m_h2co3);
});

electrolyte_test!(electrolyte_03_nacl_speciation, NaClSpeciation, |result: &ripopt::SolveResult| {
    // Na+ = Cl- ~ 0.1, pH ~ 7
    let m_na = result.x[0].exp();
    let m_cl = result.x[1].exp();
    assert!((m_na - 0.1).abs() < 1e-3, "m_Na={:.4e}", m_na);
    assert!((m_cl - 0.1).abs() < 1e-3, "m_Cl={:.4e}", m_cl);
    let m_h = result.x[2].exp();
    let ph = -(m_h.log10());
    assert!(ph > 5.0 && ph < 9.0, "pH={:.2}", ph);
});

electrolyte_test!(electrolyte_04_cacl2_nacl_mixed, CaCl2NaClMixed, |result: &ripopt::SolveResult| {
    let m_na = result.x[1].exp();
    assert!((m_na - 0.1).abs() < 1e-3, "m_Na={:.4e}", m_na);
    let m_cl = result.x[2].exp();
    assert!((m_cl - 0.2).abs() < 1e-3, "m_Cl={:.4e}", m_cl);
});

electrolyte_test!(electrolyte_05_phosphoric_acid, PhosphoricAcid, |result: &ripopt::SolveResult| {
    // pH ~ 2.25, dominant H3PO4 and H2PO4-
    let m_h = result.x[4].exp();
    let ph = -(m_h.log10());
    assert!(ph > 1.0 && ph < 4.0, "pH={:.2}", ph);
    // PO4^3- should be very small (trace species)
    let m_po4 = result.x[3].exp();
    assert!(m_po4 < 1e-6, "m_PO4={:.3e} too large", m_po4);
});

// ---------------------------------------------------------------------------
// Category 2: Phase Equilibrium
// ---------------------------------------------------------------------------

electrolyte_test!(electrolyte_06_hcl_activity, HclMeanActivity, |result: &ripopt::SolveResult| {
    // Solution: m ~ 1.0
    assert!((result.x[0] - 1.0).abs() < 0.05, "m={:.4}", result.x[0]);
    assert!(result.objective < 1e-6, "f={:.3e}", result.objective);
});

electrolyte_test!(electrolyte_07_nacl_solubility, NaClSolubility, |result: &ripopt::SolveResult| {
    // Solution: m ~ 6.14
    assert!(result.x[0] > 5.0 && result.x[0] < 8.0,
        "m_sat={:.3}", result.x[0]);
    assert!(result.objective < 1e-4, "f={:.3e}", result.objective);
});

electrolyte_test!(electrolyte_08_butanol_lle, ButanolWaterLle, |result: &ripopt::SolveResult| {
    // x_BuOH_aq should be small — salting out reduces aqueous solubility
    assert!(result.x[0] > 1e-4 && result.x[0] < 0.02,
        "x_BuOH_aq={:.4}", result.x[0]);
    assert!(result.x[1] > 0.3 && result.x[1] < 0.95,
        "x_BuOH_org={:.4}", result.x[1]);
});

electrolyte_test!(electrolyte_09_saturated_brine, SaturatedBrine, |result: &ripopt::SolveResult| {
    // m_NaCl ~ 6.14, a_w ~ 0.75, p_w ~ 2.4 kPa
    assert!(result.x[0] > 4.0 && result.x[0] < 9.0,
        "m_NaCl={:.3}", result.x[0]);
    assert!(result.x[1] > 0.6 && result.x[1] < 0.9,
        "a_w={:.3}", result.x[1]);
});

// ---------------------------------------------------------------------------
// Category 3: Parameter Fitting
// ---------------------------------------------------------------------------

electrolyte_test!(electrolyte_10_pitzer_fit, PitzerNaClFit, |result: &ripopt::SolveResult| {
    // True: [0.0765, 0.2664, 0.00127], f* ~ 0
    assert!(result.objective < 1e-6, "f={:.3e}", result.objective);
    assert!((result.x[0] - 0.0765).abs() < 0.02, "beta0={:.4}", result.x[0]);
    assert!((result.x[1] - 0.2664).abs() < 0.05, "beta1={:.4}", result.x[1]);
});

electrolyte_test!(electrolyte_11_multi_salt_dh_fit, MultiSaltDhFit, |result: &ripopt::SolveResult| {
    // f* should be near zero (synthetic data from true params)
    assert!(result.objective < 1e-4, "f={:.3e}", result.objective);
});

electrolyte_test!(electrolyte_12_enrtl_fit, EnrtlTempFit, |result: &ripopt::SolveResult| {
    // Multi-minima landscape; accept any reasonable local minimum
    assert!(result.objective < 1e-2, "f={:.3e}", result.objective);
});

// ---------------------------------------------------------------------------
// Category 4: Scale-Up
// ---------------------------------------------------------------------------

#[test]
#[ignore] // May be slow
fn electrolyte_13_seawater() {
    let problem = SeawaterSpeciation;
    let options = SolverOptions {
        tol: 1e-5,
        max_iter: 5000,
        print_level: 0,
        ..SolverOptions::default()
    };
    let start = Instant::now();
    let result = ripopt::solve(&problem, &options);
    let elapsed = start.elapsed();
    let cv = max_cv(&problem, &result.constraint_values);
    eprintln!(
        "seawater: status={:?}, obj={:.6e}, cv={:.2e}, iters={}, time={:.3}s",
        result.status, result.objective, cv, result.iterations, elapsed.as_secs_f64()
    );
    assert!(
        result.status == SolveStatus::Optimal || result.status == SolveStatus::Acceptable,
        "got {:?}", result.status
    );
    assert!(cv < 1e-3, "cv={:.2e}", cv);
    // pH should be around 8
    let m_h = result.x[8].exp();
    let ph = -(m_h.log10());
    eprintln!("  pH={:.2}, [Na+]={:.4e}, [Cl-]={:.4e}, [MgSO4]={:.4e}",
        ph, result.x[0].exp(), result.x[4].exp(), result.x[10].exp());
    assert!(ph > 6.0 && ph < 10.0, "pH={:.2} out of range", ph);
}
