import json
from pathlib import Path

import numpy as np
import pytest

from gaslift.cases import CASES_DIR, PAPER_CASE_FILES, case1_single_well, load_case, save_case
from gaslift.cli import main
from gaslift.study import default_study

ROOT = Path(__file__).resolve().parents[1]


def test_case_files_match_python_builders():
    """cases/*.json are generated from the builders; they must not drift apart."""
    from gaslift.cases import case2_multi_well

    for name, builder in (("case1", case1_single_well), ("case2", case2_multi_well)):
        case = builder()
        case.study = default_study(name)
        with open(CASES_DIR / PAPER_CASE_FILES[name]) as f:
            on_disk = json.load(f)
        assert on_disk == json.loads(json.dumps(case.to_dict()))


def test_case_round_trip(tmp_path):
    case = case1_single_well(q_oil=7.0, q_gas=1.5)
    case.study = default_study("case1")
    path = save_case(case, tmp_path / "c.json")
    again = load_case(path)
    assert again.network.node_names == case.network.node_names
    assert again.network.pipe_names == case.network.pipe_names
    assert np.array_equal(again.network.incidence_matrix(), case.network.incidence_matrix())
    assert again.specifications.to_dict() == case.specifications.to_dict()
    assert again.paper_node_numbers == case.paper_node_numbers
    assert again.study == case.study
    assert [w.name for w in again.wells] == ["well_1"]
    assert again.network.pipe("well_1:annulus").length == pytest.approx(1710.0)


def test_load_case_by_name_and_missing():
    assert load_case("case2").n_wells == 4
    with pytest.raises(FileNotFoundError):
        load_case("case9")


def test_cli_info(capsys):
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert "10.1016/j.petrol.2010.03.028" in out
    assert "case1" in out and "case2" in out


def test_cli_validate(capsys):
    assert main(["validate", "case1"]) == 0
    out = capsys.readouterr().out
    assert "3N + 4S = 66" in out
    assert "19 of N + S = 19" in out
    assert main(["validate", str(CASES_DIR / "case2_multi_well.json")]) == 0
    assert "3N + 4S = 192" in capsys.readouterr().out


def test_cli_run_no_optimize_and_eval(tmp_path, capsys):
    results = tmp_path / "results"
    assert main(["run", "case1", "--no-optimize", "--out", str(results)]) == 0
    points = json.loads((results / "case1" / "operating_points.json").read_text())
    assert all(p["converged"] for p in points)
    ref = json.loads((results / "case1" / "reference_optima.json").read_text())
    by_name = {r["name"]: r for r in ref}
    # Table 13, scenario 1 evaluated at the published decision variables
    ec = by_name["table13_scenario1"]["economics"]
    assert ec["profit_1e6_usd_yr"] == pytest.approx(26.96, abs=0.01)
    assert ec["compressor_power_kW"] == pytest.approx(1002.7, abs=0.2)
    assert by_name["table13_scenario1"]["max_abs_state_difference"] < 1e-7
    assert not (results / "case1" / "optimizations.json").exists()

    out = tmp_path / "eval"
    assert main(["eval", "--results", str(results), "--out", str(out)]) == 0
    report = (out / "report.md").read_text()
    assert "Table 13 evaluated at the published decision variables" in report
    assert "reproduced" in report
    assert (out / "simulation_vs_matlab.csv").exists()
    assert (out / "table13_at_published_optimum.csv").exists()
    text = capsys.readouterr().out
    assert "Evaluation of the reproduction" in text
