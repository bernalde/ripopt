use super::{Inertia, KktMatrix, LinearSolver, SolverError};

/// Hybrid linear solver that starts with a direct solver and automatically
/// switches to an iterative solver (preconditioned MINRES) when the direct
/// solver fails or becomes too expensive.
///
/// Switching logic:
/// - Starts in direct mode.
/// - If `factor()` fails or takes longer than `time_threshold`, switches to iterative.
/// - If the iterative solver's `solve()` fails (MINRES doesn't converge),
///   switches back to direct for the next call.
/// - Tracks consecutive failures to avoid oscillating.
pub struct HybridSolver {
    direct: super::multifrontal::MultifrontalLdl,
    iterative: super::iterative::IterativeMinres,
    /// Which solver is currently active.
    mode: HybridMode,
    /// Wall-clock time of the last direct factorization (seconds).
    last_factor_time: f64,
    /// Threshold (seconds) above which we switch to iterative. Default: 1.0s.
    time_threshold: f64,
    /// Number of consecutive iterative solve failures.
    iterative_failures: usize,
    /// Maximum consecutive iterative failures before switching back to direct.
    max_iterative_failures: usize,
    /// Number of consecutive direct factorization failures.
    direct_failures: usize,
    /// Track the last matrix for re-factoring on mode switch.
    needs_refactor: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HybridMode {
    Direct,
    Iterative,
}

impl Default for HybridSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl HybridSolver {
    pub fn new() -> Self {
        Self {
            direct: super::multifrontal::MultifrontalLdl::new(),
            iterative: super::iterative::IterativeMinres::new(),
            mode: HybridMode::Direct,
            last_factor_time: 0.0,
            time_threshold: 1.0,
            iterative_failures: 0,
            max_iterative_failures: 3,
            direct_failures: 0,
            needs_refactor: false,
        }
    }

    /// Create with a custom time threshold for switching to iterative.
    pub fn with_time_threshold(mut self, seconds: f64) -> Self {
        self.time_threshold = seconds;
        self
    }
}

impl LinearSolver for HybridSolver {
    fn factor(&mut self, matrix: &KktMatrix) -> Result<Option<Inertia>, SolverError> {
        match self.mode {
            HybridMode::Direct => {
                let start = std::time::Instant::now();
                let result = self.direct.factor(matrix);
                self.last_factor_time = start.elapsed().as_secs_f64();

                match result {
                    Ok(inertia) => {
                        self.direct_failures = 0;

                        // Check if factorization was too slow
                        if self.last_factor_time > self.time_threshold {
                            log::info!(
                                "Hybrid solver: direct factor took {:.2}s (> {:.2}s threshold), switching to iterative",
                                self.last_factor_time, self.time_threshold
                            );
                            self.mode = HybridMode::Iterative;
                            self.iterative_failures = 0;
                            // Factor with iterative too so it's ready for solve
                            let iter_result = self.iterative.factor(matrix);
                            if let Ok(iter_inertia) = iter_result {
                                return Ok(iter_inertia);
                            }
                            // If iterative factor fails, stay with the direct result
                            self.mode = HybridMode::Direct;
                        }
                        Ok(inertia)
                    }
                    Err(e) => {
                        self.direct_failures += 1;
                        log::info!(
                            "Hybrid solver: direct factor failed ({}), switching to iterative",
                            e
                        );
                        // Try iterative as fallback
                        self.mode = HybridMode::Iterative;
                        self.iterative_failures = 0;
                        self.iterative.factor(matrix)
                    }
                }
            }
            HybridMode::Iterative => {
                let result = self.iterative.factor(matrix);
                match result {
                    Ok(inertia) => Ok(inertia),
                    Err(e) => {
                        // Iterative factor failed (shouldn't normally happen),
                        // switch back to direct
                        log::info!(
                            "Hybrid solver: iterative factor failed ({}), switching to direct",
                            e
                        );
                        self.mode = HybridMode::Direct;
                        self.direct.factor(matrix)
                    }
                }
            }
        }
    }

    fn solve(&mut self, rhs: &[f64], solution: &mut [f64]) -> Result<(), SolverError> {
        match self.mode {
            HybridMode::Direct => {
                self.direct.solve(rhs, solution)
            }
            HybridMode::Iterative => {
                let result = self.iterative.solve(rhs, solution);
                match result {
                    Ok(()) => {
                        self.iterative_failures = 0;
                        Ok(())
                    }
                    Err(e) => {
                        self.iterative_failures += 1;
                        if self.iterative_failures >= self.max_iterative_failures {
                            log::info!(
                                "Hybrid solver: {} consecutive MINRES failures, switching back to direct",
                                self.iterative_failures
                            );
                            self.mode = HybridMode::Direct;
                            self.iterative_failures = 0;
                            self.needs_refactor = true;
                            // Can't solve with direct right now — it hasn't been factored
                            // with the current matrix. Return the iterative error.
                        }
                        Err(e)
                    }
                }
            }
        }
    }

    fn provides_inertia(&self) -> bool {
        // Both solvers provide (approximate) inertia
        true
    }

    fn min_diagonal(&self) -> Option<f64> {
        match self.mode {
            HybridMode::Direct => self.direct.min_diagonal(),
            HybridMode::Iterative => self.iterative.min_diagonal(),
        }
    }
}
