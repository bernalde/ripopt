//! PyO3 bindings for ripopt.
//!
//! Exposes a single `_solve` function that accepts Python callables for the
//! objective, gradient, constraints, Jacobian, and Lagrangian Hessian. The
//! high-level `minimize` API in the Python package builds those callables
//! from user-supplied `f` and `g` via JAX autodiff.

use std::cell::RefCell;

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::wrap_pyfunction;

use ripopt::{solve as ripopt_solve, NlpProblem, SolveStatus, SolverOptions};

/// Adapter that implements `NlpProblem` by calling Python callables.
///
/// Because the ripopt callback methods cannot return errors, any `PyErr`
/// raised inside a callback is stashed in `err` and re-raised after the
/// solve completes. Subsequent callbacks short-circuit once an error is set.
struct PyProblem {
    n: usize,
    m: usize,
    x0: Vec<f64>,
    x_l: Vec<f64>,
    x_u: Vec<f64>,
    g_l: Vec<f64>,
    g_u: Vec<f64>,
    f: PyObject,
    grad_f: PyObject,
    g_fn: Option<PyObject>,
    jac_g: Option<PyObject>,
    hess_l: Option<PyObject>,
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    hes_rows: Vec<usize>,
    hes_cols: Vec<usize>,
    /// Warm-start initial constraint multipliers. When Some, the IPM uses
    /// these instead of its least-squares estimate (requires warm_start=true).
    init_lam_g: Option<Vec<f64>>,
    /// Warm-start initial lower-bound multipliers.
    init_z_l: Option<Vec<f64>>,
    /// Warm-start initial upper-bound multipliers.
    init_z_u: Option<Vec<f64>>,
    err: RefCell<Option<PyErr>>,
}

impl PyProblem {
    fn stash<E: Into<PyErr>>(&self, e: E) {
        let mut slot = self.err.borrow_mut();
        if slot.is_none() {
            *slot = Some(e.into());
        }
    }

    fn has_err(&self) -> bool {
        self.err.borrow().is_some()
    }

    /// Poll Python for pending signals (e.g. KeyboardInterrupt). Returns
    /// `false` if a signal error was raised, in which case it has been
    /// stashed and the caller should bail out of the current callback.
    fn check_signals(&self, py: Python<'_>) -> bool {
        match py.check_signals() {
            Ok(()) => true,
            Err(e) => {
                self.stash(e);
                false
            }
        }
    }
}

impl NlpProblem for PyProblem {
    fn num_variables(&self) -> usize {
        self.n
    }

    fn num_constraints(&self) -> usize {
        self.m
    }

    fn bounds(&self, x_l: &mut [f64], x_u: &mut [f64]) {
        x_l.copy_from_slice(&self.x_l);
        x_u.copy_from_slice(&self.x_u);
    }

    fn constraint_bounds(&self, g_l: &mut [f64], g_u: &mut [f64]) {
        g_l.copy_from_slice(&self.g_l);
        g_u.copy_from_slice(&self.g_u);
    }

    fn initial_point(&self, x0: &mut [f64]) {
        x0.copy_from_slice(&self.x0);
    }

    fn initial_multipliers(
        &self,
        lam_g: &mut [f64],
        z_l: &mut [f64],
        z_u: &mut [f64],
    ) -> bool {
        match (
            self.init_lam_g.as_ref(),
            self.init_z_l.as_ref(),
            self.init_z_u.as_ref(),
        ) {
            (Some(lg), Some(zl), Some(zu)) => {
                if lg.len() == lam_g.len() && zl.len() == z_l.len() && zu.len() == z_u.len() {
                    lam_g.copy_from_slice(lg);
                    z_l.copy_from_slice(zl);
                    z_u.copy_from_slice(zu);
                    true
                } else {
                    false
                }
            }
            _ => false,
        }
    }

    fn objective(&self, x: &[f64]) -> f64 {
        if self.has_err() {
            return 0.0;
        }
        Python::with_gil(|py| {
            if !self.check_signals(py) {
                return 0.0;
            }
            let arr = PyArray1::from_slice_bound(py, x);
            match self.f.call1(py, (arr,)) {
                Ok(r) => r.extract::<f64>(py).unwrap_or_else(|e| {
                    self.stash(e);
                    0.0
                }),
                Err(e) => {
                    self.stash(e);
                    0.0
                }
            }
        })
    }

    fn gradient(&self, x: &[f64], grad: &mut [f64]) {
        if self.has_err() {
            grad.fill(0.0);
            return;
        }
        Python::with_gil(|py| {
            if !self.check_signals(py) {
                grad.fill(0.0);
                return;
            }
            let arr = PyArray1::from_slice_bound(py, x);
            match self.grad_f.call1(py, (arr,)) {
                Ok(r) => match r.extract::<PyReadonlyArray1<f64>>(py) {
                    Ok(v) => match v.as_slice() {
                        Ok(s) if s.len() == grad.len() => grad.copy_from_slice(s),
                        Ok(s) => {
                            self.stash(PyValueError::new_err(format!(
                                "gradient returned length {} (expected {})",
                                s.len(),
                                grad.len()
                            )));
                            grad.fill(0.0);
                        }
                        Err(e) => {
                            self.stash(e);
                            grad.fill(0.0);
                        }
                    },
                    Err(e) => {
                        self.stash(e);
                        grad.fill(0.0);
                    }
                },
                Err(e) => {
                    self.stash(e);
                    grad.fill(0.0);
                }
            }
        });
    }

    fn constraints(&self, x: &[f64], g: &mut [f64]) {
        if self.m == 0 {
            return;
        }
        if self.has_err() {
            g.fill(0.0);
            return;
        }
        let gfn = self.g_fn.as_ref().expect("constraints present => g_fn set");
        Python::with_gil(|py| {
            if !self.check_signals(py) {
                g.fill(0.0);
                return;
            }
            let arr = PyArray1::from_slice_bound(py, x);
            match gfn.call1(py, (arr,)) {
                Ok(r) => match r.extract::<PyReadonlyArray1<f64>>(py) {
                    Ok(v) => match v.as_slice() {
                        Ok(s) if s.len() == g.len() => g.copy_from_slice(s),
                        Ok(s) => {
                            self.stash(PyValueError::new_err(format!(
                                "constraints returned length {} (expected {})",
                                s.len(),
                                g.len()
                            )));
                            g.fill(0.0);
                        }
                        Err(e) => {
                            self.stash(e);
                            g.fill(0.0);
                        }
                    },
                    Err(e) => {
                        self.stash(e);
                        g.fill(0.0);
                    }
                },
                Err(e) => {
                    self.stash(e);
                    g.fill(0.0);
                }
            }
        });
    }

    fn jacobian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.jac_rows.clone(), self.jac_cols.clone())
    }

    fn jacobian_values(&self, x: &[f64], vals: &mut [f64]) {
        if self.m == 0 {
            return;
        }
        if self.has_err() {
            vals.fill(0.0);
            return;
        }
        let jfn = self.jac_g.as_ref().expect("constraints present => jac_g set");
        Python::with_gil(|py| {
            if !self.check_signals(py) {
                vals.fill(0.0);
                return;
            }
            let arr = PyArray1::from_slice_bound(py, x);
            let res = match jfn.call1(py, (arr,)) {
                Ok(r) => r,
                Err(e) => {
                    self.stash(e);
                    vals.fill(0.0);
                    return;
                }
            };
            let v = match res.extract::<PyReadonlyArray2<f64>>(py) {
                Ok(v) => v,
                Err(e) => {
                    self.stash(e);
                    vals.fill(0.0);
                    return;
                }
            };
            let n = self.n;
            if let Ok(slice) = v.as_slice() {
                for (k, (&i, &j)) in self.jac_rows.iter().zip(&self.jac_cols).enumerate() {
                    vals[k] = slice[i * n + j];
                }
            } else {
                let a = v.as_array();
                for (k, (&i, &j)) in self.jac_rows.iter().zip(&self.jac_cols).enumerate() {
                    vals[k] = a[[i, j]];
                }
            }
        });
    }

    fn hessian_structure(&self) -> (Vec<usize>, Vec<usize>) {
        (self.hes_rows.clone(), self.hes_cols.clone())
    }

    fn hessian_values(&self, x: &[f64], obj_factor: f64, lambda: &[f64], vals: &mut [f64]) {
        if self.has_err() {
            vals.fill(0.0);
            return;
        }
        let hfn = match self.hess_l.as_ref() {
            Some(h) => h,
            None => {
                vals.fill(0.0);
                return;
            }
        };
        Python::with_gil(|py| {
            if !self.check_signals(py) {
                vals.fill(0.0);
                return;
            }
            let xa = PyArray1::from_slice_bound(py, x);
            let la = PyArray1::from_slice_bound(py, lambda);
            let res = match hfn.call1(py, (xa, obj_factor, la)) {
                Ok(r) => r,
                Err(e) => {
                    self.stash(e);
                    vals.fill(0.0);
                    return;
                }
            };
            let v = match res.extract::<PyReadonlyArray2<f64>>(py) {
                Ok(v) => v,
                Err(e) => {
                    self.stash(e);
                    vals.fill(0.0);
                    return;
                }
            };
            let n = self.n;
            if let Ok(slice) = v.as_slice() {
                for (k, (&i, &j)) in self.hes_rows.iter().zip(&self.hes_cols).enumerate() {
                    vals[k] = slice[i * n + j];
                }
            } else {
                let a = v.as_array();
                for (k, (&i, &j)) in self.hes_rows.iter().zip(&self.hes_cols).enumerate() {
                    vals[k] = a[[i, j]];
                }
            }
        });
    }
}

/// Translate a Python options dict into ripopt's SolverOptions.
///
/// Known keys are whitelisted; unknown keys raise `ValueError`.
fn build_options(options: &Bound<'_, PyDict>) -> PyResult<SolverOptions> {
    let mut opts = SolverOptions::default();
    for (k, v) in options.iter() {
        let key: String = k.extract()?;
        match key.as_str() {
            "tol" => opts.tol = v.extract()?,
            "max_iter" => opts.max_iter = v.extract()?,
            "print_level" => opts.print_level = v.extract()?,
            "constr_viol_tol" => opts.constr_viol_tol = v.extract()?,
            "dual_inf_tol" => opts.dual_inf_tol = v.extract()?,
            "compl_inf_tol" => opts.compl_inf_tol = v.extract()?,
            "max_wall_time" => opts.max_wall_time = v.extract()?,
            "mu_init" => opts.mu_init = v.extract()?,
            "warm_start" => opts.warm_start = v.extract()?,
            "warm_start_bound_push" => opts.warm_start_bound_push = v.extract()?,
            "warm_start_bound_frac" => opts.warm_start_bound_frac = v.extract()?,
            "warm_start_mult_bound_push" => {
                opts.warm_start_mult_bound_push = v.extract()?
            }
            "hessian_approximation" => {
                let s: String = v.extract()?;
                match s.as_str() {
                    "limited-memory" => opts.hessian_approximation_lbfgs = true,
                    "exact" => opts.hessian_approximation_lbfgs = false,
                    other => {
                        return Err(PyValueError::new_err(format!(
                            "hessian_approximation: unknown value {:?} (expected 'exact' or 'limited-memory')",
                            other
                        )))
                    }
                }
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "unknown option: {}",
                    other
                )))
            }
        }
    }
    Ok(opts)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (
    n, m, x0, x_l, x_u, g_l, g_u,
    f, grad_f, g_fn, jac_g, hess_l,
    jac_rows, jac_cols, hes_rows, hes_cols,
    options,
    init_lam_g=None, init_z_l=None, init_z_u=None,
))]
fn _solve<'py>(
    py: Python<'py>,
    n: usize,
    m: usize,
    x0: PyReadonlyArray1<'py, f64>,
    x_l: PyReadonlyArray1<'py, f64>,
    x_u: PyReadonlyArray1<'py, f64>,
    g_l: PyReadonlyArray1<'py, f64>,
    g_u: PyReadonlyArray1<'py, f64>,
    f: PyObject,
    grad_f: PyObject,
    g_fn: Option<PyObject>,
    jac_g: Option<PyObject>,
    hess_l: Option<PyObject>,
    jac_rows: Vec<usize>,
    jac_cols: Vec<usize>,
    hes_rows: Vec<usize>,
    hes_cols: Vec<usize>,
    options: &Bound<'py, PyDict>,
    init_lam_g: Option<PyReadonlyArray1<'py, f64>>,
    init_z_l: Option<PyReadonlyArray1<'py, f64>>,
    init_z_u: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Py<PyDict>> {
    if x0.len()? != n || x_l.len()? != n || x_u.len()? != n {
        return Err(PyValueError::new_err("x0, x_l, x_u must have length n"));
    }
    if g_l.len()? != m || g_u.len()? != m {
        return Err(PyValueError::new_err("g_l, g_u must have length m"));
    }
    if m > 0 && (g_fn.is_none() || jac_g.is_none()) {
        return Err(PyValueError::new_err(
            "constraints present but g_fn or jac_g missing",
        ));
    }

    // Validate warm-start multiplier arrays if supplied. All three must be
    // present (or all absent); any partial spec is an error.
    let init_lam_g_vec = match &init_lam_g {
        Some(a) => {
            if a.len()? != m {
                return Err(PyValueError::new_err(format!(
                    "init_lam_g has length {}, expected m={}",
                    a.len()?,
                    m
                )));
            }
            Some(a.as_slice()?.to_vec())
        }
        None => None,
    };
    let init_z_l_vec = match &init_z_l {
        Some(a) => {
            if a.len()? != n {
                return Err(PyValueError::new_err(format!(
                    "init_z_l has length {}, expected n={}",
                    a.len()?,
                    n
                )));
            }
            Some(a.as_slice()?.to_vec())
        }
        None => None,
    };
    let init_z_u_vec = match &init_z_u {
        Some(a) => {
            if a.len()? != n {
                return Err(PyValueError::new_err(format!(
                    "init_z_u has length {}, expected n={}",
                    a.len()?,
                    n
                )));
            }
            Some(a.as_slice()?.to_vec())
        }
        None => None,
    };
    let any_init = init_lam_g_vec.is_some()
        || init_z_l_vec.is_some()
        || init_z_u_vec.is_some();
    let all_init = init_lam_g_vec.is_some()
        && init_z_l_vec.is_some()
        && init_z_u_vec.is_some();
    if any_init && !all_init {
        return Err(PyValueError::new_err(
            "warm-start multipliers must be supplied as a complete triple (init_lam_g, init_z_l, init_z_u) or not at all",
        ));
    }

    let prob = PyProblem {
        n,
        m,
        x0: x0.as_slice()?.to_vec(),
        x_l: x_l.as_slice()?.to_vec(),
        x_u: x_u.as_slice()?.to_vec(),
        g_l: g_l.as_slice()?.to_vec(),
        g_u: g_u.as_slice()?.to_vec(),
        f,
        grad_f,
        g_fn,
        jac_g,
        hess_l,
        jac_rows,
        jac_cols,
        hes_rows,
        hes_cols,
        init_lam_g: init_lam_g_vec,
        init_z_l: init_z_l_vec,
        init_z_u: init_z_u_vec,
        err: RefCell::new(None),
    };

    let mut opts = build_options(options)?;
    // If the user supplied warm-start multipliers, enable warm_start so that
    // SolverState::new actually calls our initial_multipliers hook.
    if all_init {
        opts.warm_start = true;
    }
    let result = ripopt_solve(&prob, &opts);

    if let Some(e) = prob.err.into_inner() {
        return Err(e);
    }

    let d = PyDict::new_bound(py);
    d.set_item("x", PyArray1::from_slice_bound(py, &result.x))?;
    d.set_item("fun", result.objective)?;
    d.set_item("status", format!("{:?}", result.status))?;
    d.set_item("success", matches!(result.status, SolveStatus::Optimal))?;
    d.set_item("iterations", result.iterations)?;
    d.set_item(
        "constraint_multipliers",
        PyArray1::from_slice_bound(py, &result.constraint_multipliers),
    )?;
    d.set_item(
        "bound_multipliers_lower",
        PyArray1::from_slice_bound(py, &result.bound_multipliers_lower),
    )?;
    d.set_item(
        "bound_multipliers_upper",
        PyArray1::from_slice_bound(py, &result.bound_multipliers_upper),
    )?;
    d.set_item(
        "constraint_values",
        PyArray1::from_slice_bound(py, &result.constraint_values),
    )?;
    d.set_item("wall_time_secs", result.diagnostics.wall_time_secs)?;
    d.set_item("final_primal_inf", result.diagnostics.final_primal_inf)?;
    d.set_item("final_dual_inf", result.diagnostics.final_dual_inf)?;
    d.set_item("final_compl", result.diagnostics.final_compl)?;
    Ok(d.into())
}

#[pymodule]
fn _ripopt(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_solve, m)?)?;
    Ok(())
}
