import numpy as np
import pytest

from gaslift import case1_single_well, case2_multi_well
from gaslift.optimization import GasLiftProblem, Objective, Reservoir


def test_problem_vector_layout_and_bounds():
    case = case2_multi_well()
    prob = GasLiftProblem(case, [Reservoir(220.0, 1.0)] * 4, Objective.OPTIMAL_DESIGN)
    x0 = prob.initial_vector([4, 5, 5, 4], [0.7, 0.8, 0.9, 0.8], 2.0)
    assert x0.shape == (12,)
    _q_oil, _q_gas, d = prob.split(x0)
    assert np.allclose(d, 2.0)
    assert len(prob.bounds()) == 12
    assert prob.cost.include_capital_cost and prob.cost.include_fuel_gas_cost


def test_constraint_uses_reservoir_inflow():
    case = case1_single_well()
    prob = GasLiftProblem(case, [Reservoir(220.0, 1.0)], Objective.MAX_PRODUCTION)
    x = prob.initial_vector(9.6, 1.81)
    c = prob.constraint_value(x)
    _sol, econ = prob.simulate(x)
    assert econ is not None
    assert c[0] == pytest.approx(220.0 - econ.downhole_pressure[0] - 9.6 / 1.0)
    assert econ.fuel_gas_cost == 0.0  # max production ignores costs


@pytest.mark.slow
def test_case1_max_production_reproduces_matlab(matlab_optima):
    r = matlab_optima["case1_maxprod_pi1_220bar"]
    case = case1_single_well()
    prob = GasLiftProblem(
        case,
        [Reservoir(220.0, 1.0)],
        Objective.MAX_PRODUCTION,
        oil_bounds=(5, 20),
        gas_bounds=(1, 4),
    )
    res = prob.solve(prob.initial_vector(10.0, 2.0))
    assert res.oil_rates[0] == pytest.approx(r["xopt"][0], abs=5e-3)
    assert res.gas_rates[0] == pytest.approx(r["xopt"][1], abs=2e-2)
    assert res.economics.profit == pytest.approx(-r["fob"], rel=1e-3)


@pytest.mark.slow
def test_case1_max_profit_reproduces_table_13(matlab_optima):
    r = matlab_optima["case1_maxprofit_oil12.79"]
    case = case1_single_well()
    from gaslift.economics import CostModel

    prob = GasLiftProblem(
        case,
        [Reservoir(220.0, 1.0)],
        Objective.MAX_PROFIT,
        cost=CostModel(oil_price=12.7875),
        oil_bounds=(5, 20),
        gas_bounds=(1, 4),
    )
    res = prob.solve(prob.initial_vector(10.0, 2.0))
    assert res.economics.profit == pytest.approx(26.96, abs=0.02)
    assert res.gas_rates[0] == pytest.approx(r["xopt"][1], abs=2e-2)
