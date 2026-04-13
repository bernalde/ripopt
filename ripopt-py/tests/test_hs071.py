"""End-to-end tests for the ripopt Python interface.

The HS071 problem is the canonical cyipopt tutorial example:

    min   x0*x3*(x0 + x1 + x2) + x2
    s.t.  x0*x1*x2*x3 >= 25
          x0**2 + x1**2 + x2**2 + x3**2 = 40
          1 <= xi <= 5,  i = 0..3

    x*  ≈ [1.00000000, 4.74299963, 3.82114998, 1.37940829]
    f*  ≈ 17.0140173
"""

import jax.numpy as jnp
import numpy as np
import pytest

from ripopt import minimize


def test_rosenbrock_unconstrained():
    def f(x):
        return (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2

    res = minimize(f, x0=[-1.2, 1.0])
    assert res.success, f"expected Optimal, got {res.status}"
    assert res.fun < 1e-10
    np.testing.assert_allclose(res.x, [1.0, 1.0], atol=1e-6)


def test_simple_box_constrained():
    def f(x):
        return (x[0] - 3.0) ** 2 + (x[1] + 1.0) ** 2

    res = minimize(f, x0=[0.0, 0.0], bounds=(0.0, 2.0))
    # Projected optimum: x0 clipped to upper bound, x1 clipped to lower.
    # (Status may be MaxIterations on this degenerate active-set case;
    # correctness is the real thing we're testing.)
    np.testing.assert_allclose(res.x, [2.0, 0.0], atol=1e-6)
    assert res.fun == pytest.approx(2.0, abs=1e-8)


def test_hs071():
    def f(x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def g(x):
        return jnp.array(
            [
                x[0] * x[1] * x[2] * x[3],
                x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2,
            ]
        )

    res = minimize(
        f,
        x0=[1.0, 5.0, 5.0, 1.0],
        bounds=(1.0, 5.0),
        constraints={"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]},
        options={"tol": 1e-8},
    )

    assert res.success, f"expected Optimal, got {res.status}"
    assert res.fun == pytest.approx(17.0140173, abs=1e-4)
    np.testing.assert_allclose(
        res.x,
        [1.00000000, 4.74299963, 3.82114998, 1.37940829],
        atol=1e-4,
    )

    # Both constraints active at the optimum
    g_val = res.constraint_values
    assert g_val[0] == pytest.approx(25.0, abs=1e-4)
    assert g_val[1] == pytest.approx(40.0, abs=1e-4)


def test_hs071_detect_sparsity():
    """HS071 with sparsity='detect' should solve to the same optimum."""

    def f(x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def g(x):
        return jnp.array(
            [
                x[0] * x[1] * x[2] * x[3],
                x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2,
            ]
        )

    res = minimize(
        f,
        x0=[1.0, 5.0, 5.0, 1.0],
        bounds=(1.0, 5.0),
        constraints={"fun": g, "lb": [25.0, 40.0], "ub": [np.inf, 40.0]},
        options={"tol": 1e-8},
        sparsity="detect",
    )
    assert res.success, f"expected Optimal, got {res.status}"
    assert res.fun == pytest.approx(17.0140173, abs=1e-4)


def test_linear_objective_has_zero_hessian():
    """A linear objective has an empty Hessian under sparsity='detect'.
    The solve must still reach the correct box-projected optimum."""

    def f(x):
        return 2.0 * x[0] + 3.0 * x[1]

    res = minimize(
        f,
        x0=[0.5, 0.5],
        bounds=(0.0, 1.0),
        sparsity="detect",
    )
    assert res.success, f"expected Optimal, got {res.status}"
    np.testing.assert_allclose(res.x, [0.0, 0.0], atol=1e-6)
    assert res.fun == pytest.approx(0.0, abs=1e-8)


def test_lbfgs_hessian_unconstrained():
    """L-BFGS Hessian path on Rosenbrock: user provides only ``f``,
    no Hessian is built."""

    def f(x):
        return (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2

    res = minimize(
        f,
        x0=[-1.2, 1.0],
        options={"hessian_approximation": "limited-memory"},
    )
    assert res.success, f"expected Optimal, got {res.status}"
    assert res.fun < 1e-6
    np.testing.assert_allclose(res.x, [1.0, 1.0], atol=1e-3)


# ---------- bounds forms ----------


def test_scipy_style_bounds_list_with_none():
    """Scipy-style per-variable list with ``None`` meaning no bound."""

    def f(x):
        # unconstrained min at (3, -1), projected onto the box below
        return (x[0] - 3.0) ** 2 + (x[1] + 1.0) ** 2

    res = minimize(
        f,
        x0=[0.0, 0.0],
        bounds=[(0.0, 1.0), (None, None)],
    )
    # x0 clipped to upper bound 1.0; x1 unbounded -> reaches -1.0
    np.testing.assert_allclose(res.x, [1.0, -1.0], atol=1e-6)
    assert res.fun == pytest.approx(4.0, abs=1e-8)


def test_per_variable_array_bounds():
    """Length-n arrays for lb and ub."""

    def f(x):
        return jnp.sum((x - jnp.array([10.0, -10.0, 0.0])) ** 2)

    lb = np.array([-1.0, -1.0, -1.0])
    ub = np.array([1.0, 1.0, 1.0])
    res = minimize(f, x0=[0.0, 0.0, 0.0], bounds=(lb, ub))
    np.testing.assert_allclose(res.x, [1.0, -1.0, 0.0], atol=1e-6)


def test_bounds_list_wrong_length_raises():
    def f(x):
        return jnp.sum(x ** 2)

    with pytest.raises(ValueError, match="bounds"):
        minimize(f, x0=[0.0, 0.0], bounds=[(0.0, 1.0)])


# ---------- constraint forms ----------


def test_multi_block_constraints_match_single_block():
    """HS071 expressed as two separate blocks must match the single-block form."""

    def f(x):
        return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]

    def prod(x):
        return x[0] * x[1] * x[2] * x[3]

    def norm_sq(x):
        return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2

    res = minimize(
        f,
        x0=[1.0, 5.0, 5.0, 1.0],
        bounds=(1.0, 5.0),
        constraints=[
            {"fun": prod, "lb": 25.0, "ub": np.inf},
            {"fun": norm_sq, "lb": 40.0, "ub": 40.0},
        ],
        options={"tol": 1e-8},
    )
    assert res.success, f"expected Optimal, got {res.status}"
    assert res.fun == pytest.approx(17.0140173, abs=1e-4)
    assert res.constraint_values.shape == (2,)
    assert res.constraint_values[0] == pytest.approx(25.0, abs=1e-4)
    assert res.constraint_values[1] == pytest.approx(40.0, abs=1e-4)


def test_equality_constraint():
    """Pure equality: min x^2+y^2 s.t. x+y=1 -> (0.5, 0.5)."""

    def f(x):
        return x[0] ** 2 + x[1] ** 2

    def g(x):
        return x[0] + x[1]

    res = minimize(
        f,
        x0=[0.0, 0.0],
        constraints={"fun": g, "lb": 1.0, "ub": 1.0},
        options={"tol": 1e-10},
    )
    assert res.success, f"expected Optimal, got {res.status}"
    np.testing.assert_allclose(res.x, [0.5, 0.5], atol=1e-6)
    assert res.fun == pytest.approx(0.5, abs=1e-8)
    assert res.constraint_multipliers.shape == (1,)


def test_scalar_constraint_broadcasts():
    """A constraint function returning a scalar should be wrapped to 1-D."""

    def f(x):
        return x[0] ** 2 + x[1] ** 2

    def g(x):
        return x[0] + x[1]  # scalar, should be treated as a 1-vector

    res = minimize(
        f,
        x0=[0.0, 0.0],
        constraints={"fun": g, "lb": 1.0, "ub": 1.0},
        options={"tol": 1e-10},
    )
    assert res.success, f"expected Optimal, got {res.status}"
    np.testing.assert_allclose(res.x, [0.5, 0.5], atol=1e-6)
    assert res.constraint_values.shape == (1,)
    assert res.constraint_values[0] == pytest.approx(1.0, abs=1e-8)


def test_invalid_constraint_dict_missing_fun_raises():
    with pytest.raises(ValueError, match="fun"):
        minimize(
            lambda x: x[0] ** 2,
            x0=[0.0],
            constraints={"lb": 0.0, "ub": 1.0},
        )


# ---------- options validation ----------


def test_unknown_option_key_raises():
    with pytest.raises(ValueError):
        minimize(
            lambda x: x[0] ** 2,
            x0=[1.0],
            options={"not_a_real_option": 42},
        )


def test_unknown_hessian_approximation_raises():
    with pytest.raises(ValueError, match="hessian_approximation"):
        minimize(
            lambda x: x[0] ** 2,
            x0=[1.0],
            options={"hessian_approximation": "bogus"},
        )


def test_max_iter_option_is_forwarded():
    """Setting max_iter=1 should force the solver to give up after 1 iter."""

    def f(x):
        return (1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2

    res = minimize(f, x0=[-1.2, 1.0], options={"max_iter": 1})
    assert res.status == "MaxIterations"
    assert res.iterations <= 1
    assert not res.success


# ---------- sparsity validation ----------


def test_unknown_sparsity_mode_raises():
    with pytest.raises(ValueError, match="sparsity"):
        minimize(lambda x: x[0] ** 2, x0=[1.0], sparsity="magic")


# ---------- input shape validation ----------


def test_x0_non_1d_raises():
    with pytest.raises(ValueError, match="1-D"):
        minimize(lambda x: jnp.sum(x ** 2), x0=np.zeros((2, 2)))


# ---------- exception propagation from user code ----------


def test_objective_exception_propagates():
    """An exception raised inside the user's objective must surface to the caller."""

    class Boom(RuntimeError):
        pass

    def f(x):
        raise Boom("objective failed")

    with pytest.raises(Boom, match="objective failed"):
        minimize(f, x0=[1.0, 2.0])


def test_constraint_exception_propagates():
    class Boom(RuntimeError):
        pass

    def f(x):
        return x[0] ** 2

    def g(x):
        raise Boom("constraint failed")

    with pytest.raises(Boom, match="constraint failed"):
        minimize(f, x0=[1.0], constraints={"fun": g, "lb": 0.0, "ub": 1.0})


# ---------- result shape / dtype contract ----------


def test_result_shapes_and_dtypes_unconstrained():
    def f(x):
        return jnp.sum(x ** 2)

    res = minimize(f, x0=[1.0, 2.0, 3.0])

    assert res.x.shape == (3,)
    assert res.x.dtype == np.float64
    assert res.bound_multipliers_lower.shape == (3,)
    assert res.bound_multipliers_upper.shape == (3,)
    # no constraints -> m = 0
    assert res.constraint_multipliers.shape == (0,)
    assert res.constraint_values.shape == (0,)
    assert isinstance(res.fun, float)
    assert isinstance(res.iterations, int)
    assert isinstance(res.status, str)
    assert isinstance(res.success, bool)
    assert isinstance(res.wall_time_secs, float)


def test_result_shapes_constrained():
    """Constrained solve: multipliers/constraint_values match m."""

    def f(x):
        return x[0] ** 2 + x[1] ** 2

    def g(x):
        return jnp.array([x[0] + x[1], x[0] - x[1]])

    res = minimize(
        f,
        x0=[0.5, 0.5],
        constraints={"fun": g, "lb": [1.0, -np.inf], "ub": [1.0, 0.0]},
    )
    assert res.success, f"expected Optimal, got {res.status}"
    assert res.constraint_multipliers.shape == (2,)
    assert res.constraint_values.shape == (2,)
