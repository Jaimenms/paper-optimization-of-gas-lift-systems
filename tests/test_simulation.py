import numpy as np
import pytest

from gaslift import (
    IntegrationOptions,
    ModelOptions,
    NetworkModel,
    SimulationOptions,
    case1_single_well,
    case2_multi_well,
)
from gaslift.linear_model import linear_initial_estimate
from gaslift.units import BAR


def _options(rtol: float, profiles: bool = False) -> SimulationOptions:
    return SimulationOptions(
        model=ModelOptions(integration=IntegrationOptions(rtol=rtol)), compute_profiles=profiles
    )


def test_linear_model_satisfies_mass_balances_and_specifications():
    case = case1_single_well(q_oil=5.0, q_gas=1.0)
    model = NetworkModel(case.network, case.specifications)
    x0 = linear_initial_estimate(model)
    v = model.split(x0)
    assert np.allclose(model.M @ v.qL + v.wL, 0.0, atol=1e-12)
    assert np.allclose(model.M @ v.qG + v.wG, 0.0, atol=1e-12)
    assert np.allclose(model.E @ x0, model.e, atol=1e-12)
    assert np.all(v.PW > 0)


@pytest.mark.parametrize("point", [(0, 0), (10, 3), (2, 1), (5, 5)])
def test_case1_matches_matlab_state(case1_grid, point):
    i, j = point
    case = case1_single_well(q_oil=case1_grid["qoil"][i], q_gas=case1_grid["qgas"][j])
    sol = case.simulator(_options(1e-3)).solve()
    assert sol.converged
    x_ref = np.array(case1_grid["state_vectors"][i][j], dtype=float)
    assert np.max(np.abs(sol.x - x_ref)) < 1e-7
    assert sol.pressure_bar("well_1:downhole") == pytest.approx(
        case1_grid["p_downhole_bar"][i][j][0], abs=1e-5
    )


def test_case1_three_inch_gas_lines(case1_grid_3in):
    i, j = 10, 3
    case = case1_single_well(
        q_oil=case1_grid_3in["qoil"][i], q_gas=case1_grid_3in["qgas"][j], gas_line_diameter_in=3.0
    )
    sol = case.simulator(_options(1e-3)).solve()
    x_ref = np.array(case1_grid_3in["state_vectors"][i][j], dtype=float)
    assert np.max(np.abs(sol.x - x_ref)) < 1e-7


@pytest.mark.parametrize("point", [(5, 5), (8, 2)])
def test_case2_matches_matlab_state(case2_grid, point):
    i, j = point
    case = case2_multi_well(q_oil=case2_grid["qoil"][i], q_gas=case2_grid["qgas"][j])
    sol = case.simulator(_options(1e-2)).solve()
    assert sol.converged
    x_ref = np.array(case2_grid["state_vectors"][i][j], dtype=float)
    assert np.max(np.abs(sol.x - x_ref)) < 1e-7


def test_solution_tables_and_profiles():
    case = case1_single_well(q_oil=10.0, q_gas=2.0)
    sol = case.simulator(_options(1e-2, profiles=True)).solve()
    nodes = sol.nodes_table()
    pipes = sol.pipes_table()
    assert nodes.loc["riser_top", "pressure_bar"] == pytest.approx(15.0)
    assert pipes.loc["riser", "liquid_flow_kg_s"] == pytest.approx(10.0)
    assert pipes.loc["riser", "gas_flow_kg_s"] == pytest.approx(2.0 + 1e-6, rel=1e-6)
    prof = sol.profile("riser")
    assert prof["pressure_bar"].iloc[0] == pytest.approx(sol.pipe_inlet_pressure[2] / BAR)
    assert prof["pressure_bar"].iloc[-1] == pytest.approx(15.0, abs=0.05)
    assert prof["pressure_bar"].is_monotonic_decreasing
    wells = sol.wells_table()
    assert wells.loc["well_1", "oil_rate_kg_s"] == pytest.approx(10.0)


def test_gas_side_pressures_are_slaves_of_production_side():
    """Changing the gas-line diameter must not change the downhole pressure."""
    a = case1_single_well(q_oil=10.0, q_gas=2.0).simulator(_options(1e-2)).solve()
    b = (
        case1_single_well(q_oil=10.0, q_gas=2.0, gas_line_diameter_in=3.0)
        .simulator(_options(1e-2))
        .solve()
    )
    assert a.pressure("well_1:downhole") == pytest.approx(b.pressure("well_1:downhole"), rel=1e-9)
    assert a.pressure("well_1:gas_junction") > b.pressure("well_1:gas_junction")


def test_warm_start_and_fallback_on_hard_point():
    case = case1_single_well(q_oil=2.0, q_gas=1.0)
    sim = case.simulator(_options(1e-3))
    assert sim.solve().converged
    case.set_well_rates(20.0, 4.0)
    sol = sim.solve()
    assert sol.converged
    assert sol.pressure_bar("well_1:downhole") == pytest.approx(314.2916, abs=1e-3)
