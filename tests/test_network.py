import numpy as np
import pytest

from gaslift import (
    FluidProperties,
    Network,
    NodeRole,
    Specifications,
    WellDesign,
    case1_single_well,
    case2_multi_well,
)
from gaslift.units import INCH


def test_case1_dimensions_and_equation_count():
    case = case1_single_well()
    lay = case.network.layout
    assert (lay.n_nodes, lay.n_pipes) == (10, 9)
    assert lay.size == 66  # paper: 66 equations
    assert case.specifications.count == lay.n_specifications == 19


def test_case2_dimensions_and_equation_count():
    case = case2_multi_well()
    lay = case.network.layout
    assert (lay.n_nodes, lay.n_pipes) == (28, 27)
    # The paper states 232 equations for this case, but 3N + 4S = 192 (see docs/errata.md, B1).
    assert lay.size == 192
    assert case.specifications.count == lay.n_specifications == 55


def test_incidence_matrix_matches_matlab(case1_grid, case2_grid):
    assert np.array_equal(
        case1_single_well().network.incidence_matrix(), np.array(case1_grid["param_M"])
    )
    assert np.array_equal(
        case2_multi_well().network.incidence_matrix(), np.array(case2_grid["param_M"])
    )


def test_paper_numbering_of_case1_nodes():
    case = case1_single_well()
    nums = case.paper_node_numbers
    assert nums["riser_top"] == 4
    assert nums["well_1:wellhead"] == 5
    assert nums["well_1:injection_point"] == 6
    assert nums["well_1:downhole"] == 7
    assert nums["well_1:annulus_top"] == 8
    assert nums["well_1:gas_junction"] == 9
    assert nums["well_1:gas_inlet"] == 10
    assert case.paper_pipe_numbers["well_1:gas_line_2"] == 9


def test_specification_matrix_shape_and_scaling():
    case = case1_single_well(q_oil=7.5, q_gas=1.5)
    e_mat, e_vec = case.specifications.matrix(pressure_scale=1e7)
    lay = case.network.layout
    assert e_mat.shape == (lay.n_specifications, lay.size)
    assert np.all(e_mat.sum(axis=1) == 1)
    i_top = case.network.node_index("riser_top")
    row = np.flatnonzero(e_mat[:, lay.index_PW(i_top)])[0]
    assert e_vec[row] == pytest.approx(0.15)
    i_dh = case.network.node_index("well_1:downhole")
    row = np.flatnonzero(e_mat[:, lay.index_wL(i_dh)])[0]
    assert e_vec[row] == pytest.approx(7.5)


def test_specifications_validation():
    net = Network(FluidProperties())
    net.add_node("a", 0, 0, 0)
    net.add_node("b", 0, 100, 0)
    net.add_pipe("p", "a", "b", diameter=0.1)
    specs = Specifications(net)
    with pytest.raises(ValueError):
        specs.matrix()
    specs.set_liquid_flow("a", 1.0).set_pressure("b", 1e6)
    with pytest.raises(ValueError):  # two specifications, three required
        specs.matrix()
    specs.set_gas_flow("a", 0.0)
    e_mat, _ = specs.matrix()
    assert e_mat.shape[0] == 3


def test_validate_flags_matlab_annulus_orientation():
    faithful = case1_single_well()
    issues = faithful.network.validate()
    assert any("annulus" in i for i in issues)
    physical = case1_single_well(annulus_orientation="physical")
    assert physical.network.validate() == []
    assert physical.network.pipe("well_1:annulus").geometry.segment_angles[0] == pytest.approx(
        -np.pi / 2
    )


def test_well_template_and_roundtrip():
    net = Network(FluidProperties(), "custom")
    net.add_node("manifold", 0, 0, -1000)
    net.add_node("outlet", 0, 1000, 0, NodeRole.OUTLET)
    net.add_pipe("riser", "manifold", "outlet", diameter=6 * INCH)
    well = net.add_well(
        WellDesign(
            name="w",
            attach_to="manifold",
            wellhead=(100, 0, -1000),
            wellbore_diameter=4 * INCH,
            wellbore_length=2000,
            annulus_length=1500,
            annulus_outer_diameter=10 * INCH,
            gas_line_diameter=2 * INCH,
            gas_line_1_path=[(0, 1000, -1000)],
            gas_line_2_path=[(0, 1000, 0), (0, 1000, -1000)],
        )
    )
    assert net.n_nodes == 8 and net.n_pipes == 7
    assert net.node(well.downhole).z == pytest.approx(-3000)
    assert net.pipe(well.gas_line_1).length == pytest.approx(np.hypot(100, 1000))
    net2 = Network.from_dict(net.to_dict())
    assert net2.node_names == net.node_names
    assert net2.pipe_names == net.pipe_names
    assert np.array_equal(net2.incidence_matrix(), net.incidence_matrix())
    assert net2.wells[0].downhole == well.downhole


def test_downstream_path():
    case = case1_single_well()
    path = case.network.downstream_path("well_1:wellbore_lower")
    assert path == [
        "well_1:wellbore_lower",
        "well_1:wellbore_upper",
        "well_1:flowline",
        "production_line_1",
        "production_line_2",
        "riser",
    ]
