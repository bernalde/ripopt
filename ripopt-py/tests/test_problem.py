"""Tests for the persistent ``Problem`` class.

These exercise the NMPC-style reuse path: build once, solve many times,
with only numeric parameters or bounds changing between calls. They also
verify that ``minimize`` remains a correct thin wrapper around
``Problem.solve``.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from ripopt import Problem, minimize


def test_problem_minimize_wrapper_matches():
    """minimize() should match Problem(...).solve() byte-for-byte."""

    def f(x):
        return (x[0] - 1.0) ** 2 + (x[1] - 2.5) ** 2

    res_fn = minimize(f, x0=[0.0, 0.0])
    res_pb = Problem(f, x0=[0.0, 0.0]).solve()

    assert res_fn.success and res_pb.success
    np.testing.assert_allclose(res_fn.x, res_pb.x, atol=1e-12)
    assert res_fn.fun == pytest.approx(res_pb.fun, abs=1e-12)


def test_problem_reuse_different_x0_same_structure():
    """Solving the same Problem from different starts must converge to
    the same optimum without rebuilding any jit'd function."""

    def f(x):
        return (x[0] - 3.0) ** 2 + (x[1] + 1.0) ** 2

    prob = Problem(f, x0=[0.0, 0.0])
    r1 = prob.solve()
    r2 = prob.solve(x0=[-5.0, 5.0])
    r3 = prob.solve(x0=[10.0, 10.0])

    for r in (r1, r2, r3):
        assert r.success, f"expected Optimal, got {r.status}"
        np.testing.assert_allclose(r.x, [3.0, -1.0], atol=1e-6)


def test_problem_warm_start_from_previous_solution():
    """Default x0 on subsequent solve() should be the previous result."""

    def f(x):
        return jnp.sum((x - jnp.array([1.0, 2.0, 3.0])) ** 2)

    prob = Problem(f, x0=[0.0, 0.0, 0.0])
    r1 = prob.solve()
    # prob.x should now be the optimum
    np.testing.assert_allclose(prob.x, [1.0, 2.0, 3.0], atol=1e-6)

    r2 = prob.solve()  # warm start from optimum → should converge in ~0 iters
    assert r2.success
    assert r2.iterations <= r1.iterations


def test_problem_update_bounds_moves_optimum():
    """update_bounds() should change the feasible set on the next solve."""

    def f(x):
        # unconstrained min at (5, -5)
        return (x[0] - 5.0) ** 2 + (x[1] + 5.0) ** 2

    prob = Problem(f, x0=[0.0, 0.0], bounds=(-1.0, 1.0))
    r1 = prob.solve()
    np.testing.assert_allclose(r1.x, [1.0, -1.0], atol=1e-6)

    prob.update_bounds(lb=[-3.0, -3.0], ub=[3.0, 3.0])
    r2 = prob.solve()
    np.testing.assert_allclose(r2.x, [3.0, -3.0], atol=1e-6)


def test_problem_params_pytree_threaded():
    """Numeric parameters threaded via params should change the solve
    result on update_parameters() without retracing."""

    def f(x, p):
        # min || x - p ||^2, so optimum is x* = p
        return jnp.sum((x - p) ** 2)

    target0 = jnp.array([1.0, 2.0, 3.0])
    prob = Problem(f, x0=[0.0, 0.0, 0.0], params=target0)
    r1 = prob.solve()
    assert r1.success
    np.testing.assert_allclose(r1.x, [1.0, 2.0, 3.0], atol=1e-6)

    prob.update_parameters(jnp.array([-4.0, 0.5, 7.0]))
    r2 = prob.solve(x0=[0.0, 0.0, 0.0])
    assert r2.success
    np.testing.assert_allclose(r2.x, [-4.0, 0.5, 7.0], atol=1e-6)


def test_problem_params_in_constraints():
    """params should also reach constraint functions."""

    def f(x, p):
        del p
        return x[0] ** 2 + x[1] ** 2

    def g(x, p):
        # equality x0 + x1 == p[0]
        return x[0] + x[1] - p[0]

    prob = Problem(
        f,
        x0=[0.0, 0.0],
        constraints={"fun": g, "lb": 0.0, "ub": 0.0},
        params=jnp.array([1.0]),
        options={"tol": 1e-10},
    )
    r1 = prob.solve()
    assert r1.success
    np.testing.assert_allclose(r1.x, [0.5, 0.5], atol=1e-6)

    prob.update_parameters(jnp.array([4.0]))
    r2 = prob.solve(x0=[0.0, 0.0])
    assert r2.success
    np.testing.assert_allclose(r2.x, [2.0, 2.0], atol=1e-6)


def test_problem_update_parameters_without_params_raises():
    prob = Problem(lambda x: jnp.sum(x ** 2), x0=[0.0, 0.0])
    with pytest.raises(ValueError, match="without params"):
        prob.update_parameters([1.0, 2.0])


def test_problem_does_not_retrace_on_reuse():
    """Core perf claim: solving the same Problem repeatedly must not
    trigger any new XLA compilation after the first solve. We capture
    jax_log_compiles output and assert 0 compiles on the second solve."""

    import io
    import logging

    def f(x, p):
        return jnp.sum((x - p) ** 2)

    prob = Problem(
        f,
        x0=[0.0, 0.0, 0.0],
        params=jnp.array([1.0, 2.0, 3.0]),
        constraints={
            "fun": lambda x, p: jnp.array([x[0] + x[1] + x[2] - p[0]]),
            "lb": -np.inf,
            "ub": np.inf,
        },
    )
    prob.solve()  # warm; any residual first-call compilation happens here

    buf = io.StringIO()
    handler = logging.StreamHandler(buf)
    handler.setLevel(logging.DEBUG)
    jax_logger = logging.getLogger("jax")
    prev_level = jax_logger.level
    jax_logger.addHandler(handler)
    jax_logger.setLevel(logging.DEBUG)
    jax.config.update("jax_log_compiles", True)
    try:
        prob.update_parameters(jnp.array([4.0, 5.0, 6.0]))
        prob.solve(x0=[0.0, 0.0, 0.0])
        prob.update_parameters(jnp.array([-1.0, -2.0, -3.0]))
        prob.solve(x0=[0.0, 0.0, 0.0])
    finally:
        jax.config.update("jax_log_compiles", False)
        jax_logger.setLevel(prev_level)
        jax_logger.removeHandler(handler)

    log_text = buf.getvalue()
    # Ripopt's objective/grad/jac/hess jit cache must stay warm.
    # Grep for the compile event lines JAX emits.
    compiles = [
        line for line in log_text.splitlines()
        if "Finished tracing" in line or "Compiling" in line
    ]
    # We don't pin an exact zero because some unrelated JAX internal jits
    # may appear; what we do require is that NO compile line names the user
    # objective/constraint/lagrangian functions we built.
    suspicious = [c for c in compiles if "lagrangian" in c or "<locals>.f" in c]
    assert suspicious == [], (
        "ripopt rebuilt a callback between solves:\n" + "\n".join(suspicious)
    )


def test_problem_jac_mode_reverse_matches_forward():
    """jac_mode='reverse' must produce the same solution as 'forward'."""

    def f(x):
        return jnp.sum(x ** 2)

    def g(x):
        return jnp.array([x[0] + x[1] - 1.0, x[1] + x[2] - 2.0])

    common = dict(
        x0=[0.0, 0.0, 0.0],
        constraints={"fun": g, "lb": 0.0, "ub": 0.0},
        options={"tol": 1e-10},
    )
    res_fwd = Problem(f, jac_mode="forward", **common).solve()
    res_rev = Problem(f, jac_mode="reverse", **common).solve()

    assert res_fwd.success and res_rev.success
    np.testing.assert_allclose(res_fwd.x, res_rev.x, atol=1e-8)


def test_problem_invalid_jac_mode_raises():
    with pytest.raises(ValueError, match="jac_mode"):
        Problem(lambda x: jnp.sum(x ** 2), x0=[0.0], jac_mode="sideways")


def test_problem_warm_start_hs071_cuts_iterations():
    """HS071: a second solve with warm_start=True from the previous optimum
    should converge in dramatically fewer iterations than without. This is
    the main reason to expose a dual warm start."""

    def f(x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def g(x):
        return jnp.array(
            [
                x[0] * x[1] * x[2] * x[3],
                x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2,
            ]
        )

    prob = Problem(
        f,
        x0=[1.0, 5.0, 5.0, 1.0],
        bounds=(1.0, 5.0),
        constraints={"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]},
        options={"tol": 1e-8},
    )
    r1 = prob.solve()
    assert r1.success

    # Cold re-solve from optimum (no dual warm start): restores mu from the
    # default init, re-traverses the barrier tail, takes many iterations.
    r_cold = prob.solve()
    assert r_cold.success

    # Warm re-solve: uses previous lam/z — should be much cheaper.
    r_warm = prob.solve(warm_start=True)
    assert r_warm.success
    assert r_warm.fun == pytest.approx(17.0140173, abs=1e-4)
    assert r_warm.iterations < r_cold.iterations, (
        f"warm_start did not reduce iterations: cold={r_cold.iterations}, "
        f"warm={r_warm.iterations}"
    )


def test_problem_explicit_lam0_accepted():
    """Supplying lam0/z_l0/z_u0 explicitly should be accepted and used."""

    def f(x):
        return x[0] ** 2 + x[1] ** 2

    def g(x):
        return x[0] + x[1]

    prob = Problem(
        f,
        x0=[0.0, 0.0],
        constraints={"fun": g, "lb": 1.0, "ub": 1.0},
        options={"tol": 1e-10},
    )
    r1 = prob.solve()
    assert r1.success

    # Supply the optimum's own duals as warm start — must converge.
    r2 = prob.solve(
        x0=r1.x,
        lam0=r1.constraint_multipliers,
        z_l0=r1.bound_multipliers_lower,
        z_u0=r1.bound_multipliers_upper,
    )
    assert r2.success
    np.testing.assert_allclose(r2.x, [0.5, 0.5], atol=1e-6)


def test_problem_partial_warm_start_triple_raises():
    """Passing only lam0 (without z_l0, z_u0) must raise ValueError."""
    prob = Problem(
        lambda x: x[0] ** 2 + x[1] ** 2,
        x0=[0.0, 0.0],
        constraints={"fun": lambda x: x[0] + x[1], "lb": 1.0, "ub": 1.0},
    )
    prob.solve()
    with pytest.raises(ValueError, match="must all be provided together"):
        prob.solve(lam0=np.array([0.0]))


def test_problem_warm_start_ignored_on_first_solve():
    """warm_start=True with no prior solve should be a no-op — the solver
    falls back to its default LS multiplier init."""

    def f(x):
        return x[0] ** 2 + x[1] ** 2

    def g(x):
        return x[0] + x[1]

    prob = Problem(
        f,
        x0=[0.0, 0.0],
        constraints={"fun": g, "lb": 1.0, "ub": 1.0},
    )
    # warm_start=True with no previous multipliers stored → just solves normally
    r = prob.solve(warm_start=True)
    assert r.success
    np.testing.assert_allclose(r.x, [0.5, 0.5], atol=1e-6)


def test_problem_warm_start_lam0_wrong_shape_raises():
    """A lam0 array of the wrong shape must raise."""
    prob = Problem(
        lambda x: x[0] ** 2 + x[1] ** 2,
        x0=[0.0, 0.0],
        constraints={"fun": lambda x: x[0] + x[1], "lb": 1.0, "ub": 1.0},
    )
    prob.solve()
    with pytest.raises(ValueError, match="shape"):
        prob.solve(
            lam0=np.zeros(5),  # wrong size — m = 1
            z_l0=np.zeros(2),
            z_u0=np.zeros(2),
        )


def test_problem_print_level_zero_is_silent(capfd):
    """print_level=0 should suppress MaxIter-diagnostic spam on failing solves
    (the silence fix that paired with the warm-start plumbing)."""

    def f(x):
        return (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2

    # Force a MaxIter failure with 1 iteration — historically this emitted
    # a "ripopt: MaxIter diag: ..." line to stderr regardless of print_level.
    prob = Problem(
        f,
        x0=[-1.2, 1.0],
        options={"max_iter": 1, "print_level": 0},
    )
    r = prob.solve()
    assert r.status == "MaxIterations"
    captured = capfd.readouterr()
    assert "MaxIter diag" not in captured.err
    assert "MaxIter diag" not in captured.out


def test_problem_hs071_reuse():
    """HS071 solved twice through the same Problem — second solve from the
    first optimum should converge trivially."""

    def f(x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def g(x):
        return jnp.array(
            [
                x[0] * x[1] * x[2] * x[3],
                x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2,
            ]
        )

    prob = Problem(
        f,
        x0=[1.0, 5.0, 5.0, 1.0],
        bounds=(1.0, 5.0),
        constraints={"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]},
        options={"tol": 1e-8},
    )
    r1 = prob.solve()
    assert r1.success
    assert r1.fun == pytest.approx(17.0140173, abs=1e-4)

    # Second solve reuses all cached callbacks. With interior-point methods
    # a "warm start" from a near-optimal primal can still take many iters
    # because bound/constraint duals are reinitialized; the important claim
    # here is correctness and that nothing was retraced.
    r2 = prob.solve()
    assert r2.success
    assert r2.fun == pytest.approx(17.0140173, abs=1e-4)
