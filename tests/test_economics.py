import numpy as np
import pytest

from gaslift import SimulationOptions, case1_single_well, case2_multi_well
from gaslift.economics import (
    CompressorParameters,
    CostModel,
    compressor_capital_cost,
    compressor_power,
    evaluate_economics,
    pipeline_capital_cost,
    turbine_capital_cost,
)
from gaslift.units import INCH


def test_compressor_power_matches_table_13():
    # Table 13, scenario "Max": QG = 1.812 kg/s, PG = 279.68 bar -> 1108.1 kW
    power = compressor_power(279.68, 1.8118, CompressorParameters())
    assert power == pytest.approx(1108.1, abs=0.5)


def test_capital_cost_correlations_match_table_18():
    power = 1797.5
    assert compressor_capital_cost(power) == pytest.approx(1.79, abs=0.01)
    assert turbine_capital_cost(power) == pytest.approx(3.95, abs=0.01)
    d = [1.77, 1.89, 1.75, 1.62]
    length = [6020.0, 6020.0, 4520.0, 4520.0]
    assert pipeline_capital_cost(d, length, "matlab") == pytest.approx(1.05, abs=0.01)
    assert pipeline_capital_cost(d, length, "mccoy") == pytest.approx(1.12, abs=0.01)


def test_case1_max_profit_economics_match_matlab(matlab_optima):
    r = matlab_optima["case1_maxprofit_oil12.79"]
    case = case1_single_well(
        q_oil=float(np.ravel(r["q_oil"])[0]), q_gas=float(np.ravel(r["q_gas"])[0])
    )
    sol = case.simulator(SimulationOptions(compute_profiles=False)).solve()
    cost = CostModel(oil_price=r["oil_price"], gas_price=r["gas_price"], include_fuel_gas_cost=True)
    ec = evaluate_economics(sol, case.wells, cost)
    assert ec.profit == pytest.approx(-r["fob"], rel=1e-5)
    assert ec.compressor_power == pytest.approx(r["out_CompressorPower"], rel=1e-5)
    assert ec.discharge_pressure == pytest.approx(r["out_DischargePressure"], rel=1e-5)
    # Table 13, scenario 1
    assert ec.profit == pytest.approx(26.96, abs=0.01)
    assert ec.compressor_power == pytest.approx(1002.7, abs=0.1)


def test_case2_design_economics_match_matlab_and_tables_18_19(matlab_optima):
    r = matlab_optima["case2_design_220bar"]
    d_in = np.array(r["xopt"][8:12])
    case = case2_multi_well(q_oil=r["q_oil"], q_gas=r["q_gas"], gas_line_diameter_in=d_in)
    sol = case.simulator(SimulationOptions(compute_profiles=False)).solve()
    assert sol.converged
    cost = CostModel(include_fuel_gas_cost=True, include_capital_cost=True)
    ec = evaluate_economics(sol, case.wells, cost)
    assert ec.profit == pytest.approx(-r["fob"], rel=1e-5)
    assert ec.compressor_cost == pytest.approx(r["out_CompressorCapitalCost"], rel=1e-5)
    assert ec.turbine_cost == pytest.approx(r["out_TurbineCapitalCost"], rel=1e-5)
    assert ec.pipeline_cost == pytest.approx(r["out_PipesCapitalCost"], rel=1e-5)
    assert np.allclose(
        ec.pipeline_lengths_m, [6020.0126, 6020.0126, 4520.0126, 4520.0126], atol=1e-3
    )
    # Tables 18 and 19, scenario 1
    assert ec.compressor_power == pytest.approx(1797.5, abs=0.1)
    assert ec.isbl == pytest.approx(6.79, abs=0.01)
    assert ec.osbl == pytest.approx(20.36, abs=0.01)
    assert ec.revenue == pytest.approx(220.47, abs=0.01)
    assert ec.profit == pytest.approx(210.82, abs=0.01)
    assert ec.injected_gas_volume.sum() == pytest.approx(81.39, abs=0.01)
    assert ec.fuel_gas_volume == pytest.approx(15.19, abs=0.01)
    assert case.network.pipe(case.wells[0].gas_line_1).geometry.inlet_diameter == pytest.approx(
        d_in[0] * INCH
    )
