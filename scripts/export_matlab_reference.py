"""Export the MATLAB results that back the paper into JSON reference files.

The original ``.mat`` files live in ``docs/raw/Programa Atual/bbnet10/CASES`` (not
versioned, 500 MB of raw material). This script converts the parts needed for
validation into ``data/reference/*.json``:

* ``case1_grid_2in.json`` / ``case1_grid_3in.json`` - rate sweeps of Case 1 with
  2-inch and 3-inch gas lines (``results_Subsea_1_2.mat``, ``results_Subsea_2.mat``);
* ``case2_grid.json`` - rate sweep of Case 2 (``results_Subsea2_1.mat``);
* ``matlab_optimization_results.json`` - optima of all ``result_opt_*.mat`` files
  used in the paper (Tables 13-19 and Figs. 13-16).

Run with ``uv run --with scipy python scripts/export_matlab_reference.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "docs" / "raw" / "Programa Atual" / "bbnet10" / "CASES"
OUT = ROOT / "data" / "reference"

GRIDS = {
    "case1_grid_2in.json": CASES / "Offshore Single Well" / "results_Subsea_1_2.mat",
    "case1_grid_3in.json": CASES / "Offshore Single Well" / "results_Subsea_2.mat",
    "case2_grid.json": CASES / "Offhore Multi Well" / "results_Subsea2_1.mat",
}

OPTIMA = {
    "case1_maxprod_pi1_220bar": "Offshore Single Well/result_opt_subsea_v01.mat",
    "case1_maxprod_pi2_220bar": "Offshore Single Well/result_opt_subsea_v01_2.mat",
    "case1_maxprod_pi5_220bar": "Offshore Single Well/result_opt_subsea_v01_3.mat",
    "case1_maxprod_pi10_220bar": "Offshore Single Well/result_opt_subsea_v01_4.mat",
    "case1_maxprod_pi0.5_220bar": "Offshore Single Well/result_opt_subsea_v01_5.mat",
    "case1_maxprod_piinf_220bar": "Offshore Single Well/result_opt_subsea_vx1.mat",
    "case1_maxprod_piinf_180bar": "Offshore Single Well/result_opt_subsea_vx2.mat",
    "case1_maxprod_piinf_200bar": "Offshore Single Well/result_opt_subsea_vx3.mat",
    "case1_maxprod_piinf_240bar": "Offshore Single Well/result_opt_subsea_vx4.mat",
    "case1_maxprod_piinf_260bar": "Offshore Single Well/result_opt_subsea_vx5.mat",
    "case1_maxprofit_oil51.15": "Offshore Single Well/result_opt_subsea_v21.mat",
    "case1_maxprofit_oil25.58": "Offshore Single Well/result_opt_subsea_v31.mat",
    "case1_maxprofit_oil12.79": "Offshore Single Well/result_opt_subsea_v41.mat",
    "case1_maxprofit_oil102.3": "Offshore Single Well/result_opt_subsea_v51.mat",
    "case2_maxprod_240bar": "Offhore Multi Well/result_opt_subsea2_v01.mat",
    "case2_maxprod_220bar": "Offhore Multi Well/result_opt_subsea2_v02.mat",
    "case2_design_220bar": "Offhore Multi Well/result_opt_subsea2_v101.mat",
    "case2_design_240bar": "Offhore Multi Well/result_opt_subsea2_v102.mat",
}


def _tolist(a):
    a = np.asarray(a).squeeze()
    return a.tolist() if a.ndim else float(a)


def _clean(values):
    return [None if np.isnan(v) else float(v) for v in values]


def export_grid(path: Path, out: Path) -> None:
    d = loadmat(path, squeeze_me=False)
    p = d["parameters"]
    pr = float(p["PR"][0, 0].squeeze())
    states, q_gas, q_oil = d["X"], d["vetor_qgas"].ravel(), d["vetor_qoleo"].ravel()
    nx, ny = states.shape
    ix_pw = p["ix_Pwell"][0, 0].astype(int).ravel() - 1
    ix_pg = p["ix_Pgaslift"][0, 0].astype(int).ravel() - 1
    res = {
        "source": str(path.relative_to(ROOT)),
        "qgas": q_gas.tolist(),
        "qoil": q_oil.tolist(),
        "shape": [nx, ny],
        "p_downhole_bar": [],
        "p_gasinlet_bar": [],
        "state_vectors": [],
        "ode_rtol": float(p["OdeRelTol"][0, 0].squeeze()),
    }
    for i in range(nx):
        rows = ([], [], [])
        for j in range(ny):
            x = np.asarray(states[i, j], dtype=float).ravel()
            if x.size == 0 or np.isnan(x).any():
                for r in rows:
                    r.append(None)
                continue
            rows[0].append(_clean(x[ix_pw] * pr / 1e5))
            rows[1].append(_clean(x[ix_pg] * pr / 1e5))
            rows[2].append(_clean(x))
        res["p_downhole_bar"].append(rows[0])
        res["p_gasinlet_bar"].append(rows[1])
        res["state_vectors"].append(rows[2])
    for fld in ["N", "S", "pipe_length", "HL0", "e", "x", "y", "z"]:
        res["param_" + fld] = _tolist(p[fld][0, 0])
    res["param_M"] = np.asarray(p["M"][0, 0]).tolist()
    res["param_E"] = np.asarray(p["E"][0, 0]).tolist()
    for fld in ["pipe_L", "pipe_teta", "pipe_D", "pipe_rough"]:
        res["param_" + fld] = [np.asarray(c).ravel().tolist() for c in p[fld][0, 0].ravel()]
    out.write_text(json.dumps(res))
    print(f"{out.name}: {nx} x {ny} grid")


def export_optima(out: Path) -> None:
    result = {}
    for key, rel in OPTIMA.items():
        d = loadmat(CASES / rel, squeeze_me=False)
        po, so, p = d["parametersOpt"][0, 0], d["saidasOpt"][0, 0], d["parameters"][0, 0]
        pr = float(p["PR"].squeeze())
        x = d["x"]
        x = x[:, -1] if x.ndim == 2 and x.shape[1] > 1 else x.ravel()
        rec = {
            "source": rel,
            "xopt": _tolist(d["xopt"]),
            "xopt0": _tolist(d["xopt0"]),
            "fob": float(d["fob"].squeeze()),
            "LB": _tolist(d["LB"]),
            "UB": _tolist(d["UB"]),
            "pwell_bar": _tolist(po["pwell"]),
            "IP": _tolist(po["IP"]),
            "oil_price": float(po["oilPrice"].squeeze()),
            "gas_price": float(po["gasPrice"].squeeze()),
            "fixedCost": int(po["fixedCost"].squeeze()),
            "fuelGasCost": int(po["fuelGasCost"].squeeze()),
            "p_downhole_bar": _tolist(x[p["ix_Pwell"].astype(int).ravel() - 1] * pr / 1e5),
            "p_gasinlet_bar": _tolist(x[p["ix_Pgaslift"].astype(int).ravel() - 1] * pr / 1e5),
            "q_oil": _tolist(x[p["ix_qLwell"].astype(int).ravel() - 1]),
            "q_gas": _tolist(x[p["ix_qGgaslift"].astype(int).ravel() - 1]),
            "state_vector": _tolist(x),
        }
        for fld in so.dtype.names:
            v = so[fld]
            if v.dtype.kind in "fiu" and v.size > 0:
                rec["out_" + fld] = _tolist(v)
        result[key] = rec
    out.write_text(json.dumps(result, indent=1))
    print(f"{out.name}: {len(result)} optima")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, path in GRIDS.items():
        export_grid(path, OUT / name)
    export_optima(OUT / "matlab_optimization_results.json")


if __name__ == "__main__":
    main()
