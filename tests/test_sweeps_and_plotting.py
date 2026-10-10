import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from gaslift import SimulationOptions, case1_single_well
from gaslift.plotting import (
    plot_case_topology,
    plot_digraph,
    plot_flow_pattern_map,
    plot_pressure_path,
    plot_rate_contours,
)
from gaslift.sweeps import sweep_rates


def test_sweep_small_grid_and_contour_plot():
    case = case1_single_well()
    sweep = sweep_rates(case, np.array([1.0, 2.0]), np.array([5.0, 10.0]))
    assert sweep.converged.all()
    assert sweep.downhole_pressure.shape == (2, 2, 1)
    # more oil at the same gas rate needs a higher downhole pressure
    assert sweep.downhole(0)[1, 0] > sweep.downhole(0)[0, 0]
    # the platform node is below the gas junction node by the gas column of line 2
    assert np.all(sweep.gas_inlet(0) < sweep.gas_junction(0))
    fig, ax = plt.subplots()
    plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), ax, levels=5)
    plt.close(fig)


def test_topology_digraph_and_profile_plots():
    case = case1_single_well(q_oil=10.0, q_gas=2.0)
    fig = plot_case_topology(case.network)
    plt.close(fig)
    fig, ax = plt.subplots()
    plot_digraph(
        case.network, ax, node_numbers=case.paper_node_numbers, pipe_numbers=case.paper_pipe_numbers
    )
    plt.close(fig)
    fig, ax = plt.subplots()
    plot_flow_pattern_map(ax, n=40)
    plt.close(fig)
    sol = case.simulator(SimulationOptions()).solve()
    fig, ax = plt.subplots()
    plot_pressure_path(
        sol, ["well_1:wellbore_lower", "well_1:wellbore_upper", "well_1:flowline"], ax
    )
    assert len(ax.lines) == 1
    plt.close(fig)
