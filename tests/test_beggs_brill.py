import numpy as np
import pytest

from gaslift.beggs_brill import (
    DISTRIBUTED,
    INTERMITTENT,
    SEGREGATED,
    TRANSITION,
    churchill_friction_factor,
    closure_relations,
    flow_pattern,
)
from gaslift.fluids import FluidProperties


@pytest.mark.parametrize(
    ("lambda_l", "froude", "expected"),
    [
        (0.001, 1.0, SEGREGATED),
        (0.5, 1.0, INTERMITTENT),
        (0.5, 1000.0, DISTRIBUTED),
        (0.05, 3.0, TRANSITION),
        (0.9, 0.2, INTERMITTENT),
        (0.9, 0.0005, SEGREGATED),
    ],
)
def test_flow_pattern_regions(lambda_l, froude, expected):
    weights, _ = flow_pattern(np.array([lambda_l]), np.array([froude]))
    assert weights.index[0] == pytest.approx(expected, abs=1e-6)


def test_churchill_laminar_limit():
    re = np.array([100.0, 500.0])
    f = churchill_friction_factor(re, np.zeros(2))
    assert np.allclose(f, 64.0 / re, rtol=0.02)


def test_churchill_turbulent_smooth_pipe():
    f = churchill_friction_factor(np.array([1e5]), np.array([0.0]))
    assert 0.017 < f[0] < 0.019  # Blasius ~0.0178


def test_holdup_bounds_and_slip():
    fluid = FluidProperties()
    p = np.full(5, 150e5)
    v_sl = np.array([0.1, 0.5, 1.0, 2.0, 3.0])
    v_sg = np.array([2.0, 1.0, 1.0, 0.5, 0.1])
    res = closure_relations(
        p, v_sl, v_sg, np.full(5, 0.1), np.full(5, 1e-4), np.full(5, np.pi / 2), fluid
    )
    assert np.all(res.liquid_holdup >= res.no_slip_holdup - 1e-9)
    assert np.all(res.liquid_holdup <= 1.0)
    assert np.all(res.liquid_holdup + res.gas_holdup == pytest.approx(1.0))
    assert np.all(res.gas_velocity >= res.liquid_velocity - 1e-6)


def test_direction_follows_liquid_flow_sign():
    fluid = FluidProperties()
    fwd = closure_relations(150e5, 0.5, 1.0, 0.1, 1e-4, 0.0, fluid)
    bwd = closure_relations(150e5, -0.5, -1.0, 0.1, 1e-4, 0.0, fluid)
    assert fwd.shear_stress > 0
    assert bwd.shear_stress == pytest.approx(-fwd.shear_stress)
    assert bwd.liquid_holdup == pytest.approx(fwd.liquid_holdup)


def test_gas_density_polytropic():
    fluid = FluidProperties()
    assert fluid.gas_density(1e7) == pytest.approx(80.0)
    assert fluid.gas_density(2e7) == pytest.approx(80.0 * 2 ** (1 / 1.25))
