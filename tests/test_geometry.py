import numpy as np
import pytest

from gaslift import PipeGeometry, annulus_equivalent_diameter, case1_single_well
from gaslift.units import INCH


def test_straight_pipe_length_and_angle():
    g = PipeGeometry.straight((0, 0, 0), (0, 0, -100), diameter=0.1)
    assert g.length == pytest.approx(100.0)
    assert g.segment_angles[0] == pytest.approx(-np.pi / 2)
    assert g.elevation_change == pytest.approx(-100.0)


def test_knots_follow_matlab_scheme():
    g = PipeGeometry.from_points([(0, 0, 0), (0, 1500, 0), (0, 2505, -100)], diameter=0.15)
    s, angle, _d, _r = g.knots
    assert s[0] == 0.0
    assert s[1] == pytest.approx(1499.0)
    assert s[2] == pytest.approx(1501.0)
    assert s[-1] == pytest.approx(g.length)
    assert angle[0] == angle[1] == 0.0
    assert angle[2] == angle[3] == pytest.approx(np.arcsin(-100 / np.hypot(1005, 100)))
    a, dd, _rr = g.properties_at([0.0, 1500.0, 2000.0])
    assert a[0] == 0.0
    assert a[1] == pytest.approx(0.5 * (angle[1] + angle[2]))
    assert a[2] == pytest.approx(angle[3])
    assert np.all(dd == 0.15)


def test_case1_pipe_lengths_match_paper_tables():
    net = case1_single_well().network
    expected = {
        "production_line_1": 1500.0,
        "production_line_2": 2504.98756,
        "riser": 1500.0,
        "well_1:flowline": 500.0,
        "well_1:wellbore_upper": 1710.0,
        "well_1:wellbore_lower": 338.0,
        "well_1:annulus": 1710.0,
        "well_1:gas_line_1": 4505.01256,
        "well_1:gas_line_2": 1515.0,
    }
    for name, length in expected.items():
        assert net.pipe(name).length == pytest.approx(length, abs=1e-4)


def test_annulus_equivalent_diameter_models():
    do, di = 10 * INCH, 4 * INCH
    assert annulus_equivalent_diameter(do, di, "hydraulic") == pytest.approx(6 * INCH)
    assert annulus_equivalent_diameter(do, di, "matlab") == pytest.approx(
        (do**2 + di**2) / (do + di)
    )
    with pytest.raises(ValueError):
        annulus_equivalent_diameter(do, di, "other")


def test_uniform_diameter_copy_and_serialisation_roundtrip():
    g = PipeGeometry.from_points([(0, 0, 0), (10, 0, 0)], diameter=[0.1, 0.2], roughness=1e-5)
    g2 = g.with_uniform_diameter(0.3)
    assert np.all(g2.diameters == 0.3)
    g3 = PipeGeometry.from_dict(g.to_dict())
    assert np.allclose(g3.points, g.points)
    assert np.allclose(g3.diameters, g.diameters)
