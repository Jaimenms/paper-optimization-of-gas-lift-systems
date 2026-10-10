"""Compare the results of :func:`gaslift.study.run_study` with the reference data.

Two references are used:

* ``data/reference/*.json`` - results of the original MATLAB program (state vectors
  of the rate sweeps and the optima of ``result_opt_*.mat``);
* ``data/reference/published_tables.json`` - the values printed in Tables 13-19 of
  the paper.

:func:`evaluate` reads ``results/<case>/*.json`` and writes a markdown report plus
one CSV file per comparison into ``results/eval/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .study import REFERENCE_DIR, load_reference

CASE_DIRS = ("case1", "case2")

TABLE13_ROWS = {
    "produced_oil_kg_s": ("economics", "total_oil_rate_kg_s"),
    "injected_gas_kg_s": ("economics", "total_gas_injection_kg_s"),
    "downhole_pressure_bar": ("downhole_pressure_bar", 0),
    "compressor_power_kW": ("economics", "compressor_power_kW"),
    "discharge_pressure_bar": ("economics", "discharge_pressure_bar"),
    "profit_1e6_usd_yr": ("economics", "profit_1e6_usd_yr"),
}
TABLE18_ROWS = {
    "compressor_power_kW": "compressor_power_kW",
    "discharge_pressure_bar": "discharge_pressure_bar",
    "compressor_cost_1e6_usd": "compressor_cost_1e6_usd",
    "turbine_cost_1e6_usd": "turbine_cost_1e6_usd",
    "pipeline_cost_1e6_usd": "pipeline_cost_1e6_usd",
    "fuel_gas_1e6_m3_yr": "fuel_gas_1e6_m3_yr",
}
TABLE19_ROWS = {
    "total_injected_gas_1e6_m3_yr": "total_injected_gas_1e6_m3_yr",
    "total_oil_volume_1e6_bbl_yr": "total_oil_volume_1e6_bbl_yr",
    "fuel_gas_cost_1e6_usd_yr": "fuel_gas_cost_1e6_usd_yr",
    "annualized_capital_cost_1e6_usd_yr": "annualized_capital_cost_1e6_usd_yr",
    "ISBL_1e6_usd": "ISBL_1e6_usd",
    "OSBL_1e6_usd": "OSBL_1e6_usd",
    "revenue_1e6_usd_yr": "revenue_1e6_usd_yr",
    "profit_1e6_usd_yr": "profit_1e6_usd_yr",
}


# --------------------------------------------------------------------------- loading
def _load(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as f:
        return json.load(f)


def load_results(results_dir: str | Path) -> dict[str, dict[str, list[dict]]]:
    """Results per case: operating points, reference optima, optimisations and summary."""
    root = Path(results_dir)
    out = {}
    for case in CASE_DIRS:
        d = root / case
        if not d.exists():
            continue
        out[case] = {
            "operating_points": _load(d / "operating_points.json"),
            "reference_optima": _load(d / "reference_optima.json"),
            "optimizations": _load(d / "optimizations.json"),
            "summary": _load(d / "summary.json") or {},
        }
    return out


def _get(rec: dict, path) -> float:
    if isinstance(path, tuple):
        a, b = path
        v = rec.get(a)
        if v is None:
            return np.nan
        return float(v[b]) if v is not None else np.nan
    return float(rec.get(path, np.nan))


def _grid_value(grid: dict, key: str, i: int, j: int, well: int = 0) -> float:
    cell = grid[key][i][j]
    if cell is None:
        return np.nan
    return float(np.atleast_1d(cell)[well])


# --------------------------------------------------------------------------- comparisons
def simulation_vs_matlab(results: dict) -> pd.DataFrame:
    """Operating points with a MATLAB grid reference: state-vector and pressure deviations."""
    rows = []
    grids: dict[str, dict] = {}
    for case, res in results.items():
        for rec in res["operating_points"]:
            ref = rec.get("reference")
            if not ref:
                continue
            grid = grids.setdefault(ref["grid"], load_reference(ref["grid"]))
            i, j = ref["index"]
            x_ref = grid["state_vectors"][i][j]
            if x_ref is None or not rec["converged"]:
                rows.append(
                    {
                        "case": case,
                        "point": rec["name"],
                        "q_oil": rec["q_oil"],
                        "q_gas": rec["q_gas"],
                        "converged": rec["converged"],
                        "matlab_converged": x_ref is not None,
                    }
                )
                continue
            x = np.asarray(rec["state_vector"], dtype=float)
            x_ref = np.asarray(x_ref, dtype=float)
            p_ref = _grid_value(grid, "p_downhole_bar", i, j)
            pg_ref = _grid_value(grid, "p_gasinlet_bar", i, j)
            rows.append(
                {
                    "case": case,
                    "point": rec["name"],
                    "q_oil": rec["q_oil"],
                    "q_gas": rec["q_gas"],
                    "converged": True,
                    "matlab_converged": True,
                    "downhole_pressure_bar": rec["downhole_pressure_bar"][0],
                    "downhole_pressure_matlab_bar": p_ref,
                    "gas_junction_pressure_bar": rec["gas_junction_pressure_bar"][0],
                    "gas_junction_pressure_matlab_bar": pg_ref,
                    "max_abs_state_difference": float(np.max(np.abs(x - x_ref))),
                    "max_rel_state_difference": float(
                        np.max(np.abs(x - x_ref) / np.maximum(np.abs(x_ref), 1e-12))
                    ),
                }
            )
    return pd.DataFrame(rows)


def _long_rows(table: str, item: str, column: str, published, matlab, computed, unit=""):
    dev = np.nan if computed is None or published is None else float(computed) - float(published)
    rel = np.nan if published in (None, 0) or np.isnan(dev) else dev / float(published)
    return {
        "table": table,
        "item": item,
        "column": column,
        "published": published,
        "matlab": matlab,
        "computed": computed,
        "deviation": dev,
        "relative_deviation": rel,
    }


def _records_by_name(recs: list[dict]) -> dict[str, dict]:
    return {r["name"]: r for r in recs}


def table13(results: dict, published: dict, source: str) -> pd.DataFrame:
    """Table 13 from ``source`` records (``"optimizations"`` or ``"reference_optima"``)."""
    res = results.get("case1")
    if res is None:
        return pd.DataFrame()
    recs = _records_by_name(res[source])
    pub = published["table13"]
    names = {
        "1": "table13_scenario1",
        "2": "table13_scenario2",
        "3": "table13_scenario3",
        "Max": "fig15_pi1",
    }
    rows = []
    for k, col in enumerate(pub["columns"]):
        rec = recs.get(names[col])
        if rec is None or not rec.get("converged", rec.get("success", True)):
            continue
        for item, path in TABLE13_ROWS.items():
            matlab = None
            m = rec.get("matlab", {})
            if source == "optimizations" and m:
                matlab = {
                    "produced_oil_kg_s": float(np.sum(m["q_oil"])),
                    "injected_gas_kg_s": float(np.sum(m["q_gas"])),
                    "profit_1e6_usd_yr": m["profit"],
                    "downhole_pressure_bar": m["downhole_pressure_bar"][0],
                }.get(item)
            elif m:
                matlab = {
                    "produced_oil_kg_s": float(np.sum(np.atleast_1d(m["xopt"])[:1])),
                    "injected_gas_kg_s": float(np.atleast_1d(m["xopt"])[1]),
                    "profit_1e6_usd_yr": -m["fob"],
                    "downhole_pressure_bar": float(np.atleast_1d(m["p_downhole_bar"])[0]),
                    "compressor_power_kW": m.get("out_CompressorPower"),
                    "discharge_pressure_bar": m.get(
                        "out_DischargePressure", float(np.max(np.atleast_1d(m["p_gasinlet_bar"])))
                    ),
                }.get(item)
            rows.append(
                _long_rows("Table 13", item, col, pub["rows"][item][k], matlab, _get(rec, path))
            )
    return pd.DataFrame(rows)


def table14_15(results: dict, published: dict) -> pd.DataFrame:
    res = results.get("case2")
    if res is None:
        return pd.DataFrame()
    recs = _records_by_name(res["optimizations"])
    rows = []
    for table, name in (
        ("table14", "table14_maxprod_220bar"),
        ("table15", "table15_maxprod_240bar"),
    ):
        rec = recs.get(name)
        if rec is None:
            continue
        pub = published[table]
        m = rec.get("matlab", {})
        for w in range(4):
            col = pub["columns"][w]
            rows.append(
                _long_rows(
                    f"Table {table[-2:]}",
                    "downhole_pressure_bar",
                    col,
                    pub["rows"]["downhole_pressure_bar"][w],
                    m.get("downhole_pressure_bar", [None] * 4)[w],
                    rec["downhole_pressure_bar"][w],
                )
            )
            rows.append(
                _long_rows(
                    f"Table {table[-2:]}",
                    "produced_oil_kg_s",
                    col,
                    pub["rows"]["produced_oil_kg_s"][w],
                    m.get("q_oil", [None] * 4)[w],
                    rec["q_oil"][w],
                )
            )
            rows.append(
                _long_rows(
                    f"Table {table[-2:]}",
                    "injected_gas_kg_s",
                    col,
                    pub["rows"]["injected_gas_kg_s"][w],
                    m.get("q_gas", [None] * 4)[w],
                    rec["q_gas"][w],
                )
            )
    return pd.DataFrame(rows)


def table17(results: dict, published: dict) -> pd.DataFrame:
    res = results.get("case2")
    if res is None:
        return pd.DataFrame()
    recs = _records_by_name(res["optimizations"])
    pub = published["table17"]
    rows = []
    for s, name in enumerate(("table17_19_design_220bar", "table17_19_design_240bar")):
        rec = recs.get(name)
        if rec is None:
            continue
        m = rec.get("matlab", {})
        col = pub["columns"][s]
        for w in range(4):
            for item, key, mkey in (
                ("injected_gas_kg_s", "q_gas", "q_gas"),
                ("produced_oil_kg_s", "q_oil", "q_oil"),
                ("gas_line_diameter_in", "diameters_in", "diameters_in"),
            ):
                rows.append(
                    _long_rows(
                        "Table 17",
                        f"{item} well {w + 1}",
                        col,
                        pub["rows"][item][s][w],
                        (m.get(mkey) or [None] * 4)[w],
                        rec[key][w],
                    )
                )
    return pd.DataFrame(rows)


def table18_19(results: dict, published: dict, source: str) -> pd.DataFrame:
    res = results.get("case2")
    if res is None:
        return pd.DataFrame()
    recs = _records_by_name(res[source])
    rows = []
    for s, name in enumerate(("table17_19_design_220bar", "table17_19_design_240bar")):
        rec = recs.get(name)
        if rec is None or "economics" not in rec:
            continue
        ec = rec["economics"]
        m = rec.get("matlab", {})
        matlab_map = {
            "compressor_power_kW": m.get("out_CompressorPower"),
            "discharge_pressure_bar": m.get("out_DischargePressure"),
            "compressor_cost_1e6_usd": m.get("out_CompressorCapitalCost"),
            "turbine_cost_1e6_usd": m.get("out_TurbineCapitalCost"),
            "pipeline_cost_1e6_usd": m.get("out_PipesCapitalCost"),
            "fuel_gas_1e6_m3_yr": m.get("out_FuelGasFlowRate"),
            "total_injected_gas_1e6_m3_yr": (
                float(np.sum(m["out_GasLiftFlowRate"])) if "out_GasLiftFlowRate" in m else None
            ),
            "total_oil_volume_1e6_bbl_yr": (
                float(np.sum(m["out_OilFlowRate"])) if "out_OilFlowRate" in m else None
            ),
            "fuel_gas_cost_1e6_usd_yr": m.get("out_AnnualFuelGasCost"),
            "annualized_capital_cost_1e6_usd_yr": m.get("out_AnnualTotalCapitalCost"),
            "ISBL_1e6_usd": m.get("out_ISBL"),
            "OSBL_1e6_usd": m.get("out_OSBL"),
            "revenue_1e6_usd_yr": m.get("out_InValue"),
            "profit_1e6_usd_yr": m.get("out_Profit", m.get("profit")),
        }
        for table, mapping in (("table18", TABLE18_ROWS), ("table19", TABLE19_ROWS)):
            pub = published[table]
            col = pub["columns"][s]
            for item, key in mapping.items():
                rows.append(
                    _long_rows(
                        f"Table {table[-2:]}",
                        item,
                        col,
                        pub["rows"][item][s],
                        matlab_map.get(item),
                        ec.get(key),
                    )
                )
    return pd.DataFrame(rows)


def optimization_vs_matlab(results: dict) -> pd.DataFrame:
    """Every optimisation with a MATLAB reference: decision variables and objective."""
    rows = []
    for case, res in results.items():
        for rec in res["optimizations"]:
            m = rec.get("matlab")
            if not m:
                continue
            row = {
                "case": case,
                "problem": rec["name"],
                "paper": ", ".join(rec.get("paper", [])),
                "objective": rec["objective"],
                "success": rec["success"],
                "n_evaluations": rec["n_evaluations"],
                "elapsed_s": rec["elapsed_s"],
                "oil_kg_s": float(np.sum(rec["q_oil"])),
                "oil_matlab_kg_s": float(np.sum(m["q_oil"])),
                "gas_kg_s": float(np.sum(rec["q_gas"])),
                "gas_matlab_kg_s": float(np.sum(m["q_gas"])),
                "objective_1e6_usd_yr": -rec["objective_value"],
                "objective_matlab_1e6_usd_yr": m["profit"],
                "max_abs_oil_difference": float(
                    np.max(np.abs(np.subtract(rec["q_oil"], m["q_oil"])))
                ),
                "max_abs_gas_difference": float(
                    np.max(np.abs(np.subtract(rec["q_gas"], m["q_gas"])))
                ),
            }
            if rec.get("diameters_in") is not None and m.get("diameters_in") is not None:
                row["max_abs_diameter_difference_in"] = float(
                    np.max(np.abs(np.subtract(rec["diameters_in"], m["diameters_in"])))
                )
            row["objective_relative_difference"] = (row["objective_1e6_usd_yr"] - m["profit"]) / m[
                "profit"
            ]
            rows.append(row)
    return pd.DataFrame(rows)


def reference_optima_vs_matlab(results: dict) -> pd.DataFrame:
    """Simulation at the MATLAB decision variables: state vectors and economics vs MATLAB."""
    rows = []
    for case, res in results.items():
        for rec in res["reference_optima"]:
            m = rec.get("matlab", {})
            row = {
                "case": case,
                "problem": rec["name"],
                "reference": rec["reference"],
                "converged": rec["converged"],
                "max_abs_state_difference": rec.get("max_abs_state_difference"),
            }
            if "economics" in rec:
                ec = rec["economics"]
                row.update(
                    {
                        "profit_1e6_usd_yr": ec["profit_1e6_usd_yr"],
                        "profit_matlab_1e6_usd_yr": -m["fob"],
                        "compressor_power_kW": ec["compressor_power_kW"],
                        "compressor_power_matlab_kW": m.get("out_CompressorPower"),
                        "discharge_pressure_bar": ec["discharge_pressure_bar"],
                        "discharge_pressure_matlab_bar": m.get("out_DischargePressure"),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- report
def _md_table(df: pd.DataFrame, floatfmt: str = "{:.4g}") -> str:
    if df.empty:
        return "_no data_\n"
    df = df.copy()
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines) + "\n"


def _wide(df: pd.DataFrame) -> pd.DataFrame:
    """Long comparison table -> one row per item and column with published/matlab/computed."""
    if df.empty:
        return df
    out = df[["table", "item", "column", "published", "matlab", "computed", "deviation"]].copy()
    return out


def _status(df: pd.DataFrame, abs_tol: float, rel_tol: float) -> str:
    if df.empty:
        return "not computed"
    dev = df["deviation"].abs()
    rel = df["relative_deviation"].abs()
    ok = (dev <= abs_tol) | (rel <= rel_tol)
    n_bad = int((~ok).sum())
    worst = float(rel.max()) if rel.notna().any() else float("nan")
    if n_bad == 0:
        return f"reproduced ({len(df)} values, max relative deviation {worst:.2e})"
    return f"{n_bad} of {len(df)} values outside tolerance (max relative deviation {worst:.2e})"


def evaluate(results_dir: str | Path, out_dir: str | Path, progress=None) -> Path:
    """Compare ``results_dir`` with the reference data; write CSVs and ``report.md``."""
    results = load_results(results_dir)
    if not results:
        raise FileNotFoundError(f"no results found in {results_dir}; run `gaslift run` first")
    published = load_reference("published_tables.json")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    sim = simulation_vs_matlab(results)
    ref_opt = reference_optima_vs_matlab(results)
    t13_ref = table13(results, published, "reference_optima")
    t18_ref = table18_19(results, published, "reference_optima")
    opt = optimization_vs_matlab(results)
    t13 = table13(results, published, "optimizations")
    t1415 = table14_15(results, published)
    t17 = table17(results, published)
    t1819 = table18_19(results, published, "optimizations")

    files = {
        "simulation_vs_matlab.csv": sim,
        "reference_optima_vs_matlab.csv": ref_opt,
        "table13_at_published_optimum.csv": t13_ref,
        "table18_19_at_published_optimum.csv": t18_ref,
        "optimization_vs_matlab.csv": opt,
        "table13.csv": t13,
        "table14_15.csv": t1415,
        "table17.csv": t17,
        "table18_19.csv": t1819,
    }
    for name, df in files.items():
        if not df.empty:
            df.to_csv(out / name, index=False)

    has_opt = any(res["optimizations"] for res in results.values())
    lines = ["# Evaluation of the reproduction\n"]
    lines.append(
        "Results of `gaslift run` compared with the original MATLAB results "
        "(`data/reference/*.json`) and with the tables printed in the paper "
        "(`data/reference/published_tables.json`).\n"
    )
    for case, res in results.items():
        s = res.get("summary", {})
        lines.append(
            f"* `{case}`: {s.get('n_operating_points', 0)} operating points, "
            f"{s.get('n_reference_optima', 0)} published optima re-evaluated, "
            f"{s.get('n_optimizations', 0)} optimisations run "
            f"({s.get('elapsed_s', 0):.0f} s, gaslift {s.get('gaslift_version', '?')}, "
            f"{s.get('timestamp', '')})."
        )
    lines.append("")

    lines.append("## Summary\n")
    summary = pd.DataFrame(
        [
            {
                "item": "Steady-state simulation vs MATLAB state vectors (Figs. 6-8, 11-12 points)",
                "status": (
                    f"max |dx| = {sim['max_abs_state_difference'].max():.2e} over "
                    f"{int(sim['converged'].sum())} points"
                    if not sim.empty and "max_abs_state_difference" in sim
                    else "not computed"
                ),
            },
            {
                "item": "Simulation at the published optima vs MATLAB (Tables 13-19 points)",
                "status": (
                    f"max |dx| = {ref_opt['max_abs_state_difference'].max():.2e} over "
                    f"{len(ref_opt)} optima"
                    if not ref_opt.empty
                    else "not computed"
                ),
            },
            {
                "item": "Table 13 evaluated at the published decision variables",
                "status": _status(t13_ref, 0.05, 1e-3),
            },
            {
                "item": "Tables 18-19 evaluated at the published decision variables",
                "status": _status(t18_ref, 0.05, 1e-3),
            },
            {
                "item": "Optimisation (SLSQP) vs MATLAB fmincon optima (Figs. 13-16, Tables 13-19)",
                "status": (
                    f"max relative objective difference "
                    f"{opt['objective_relative_difference'].abs().max():.2e}; "
                    f"max |dQ_oil| {opt['max_abs_oil_difference'].max():.3f} kg/s; "
                    f"max |dQ_gas| {opt['max_abs_gas_difference'].max():.3f} kg/s"
                    if not opt.empty
                    else "not run (`--no-optimize`)"
                ),
            },
            {"item": "Table 13 (optimised)", "status": _status(t13, 0.05, 5e-3)},
            {"item": "Tables 14-15 (optimised)", "status": _status(t1415, 0.05, 5e-3)},
            {"item": "Table 17 (optimised)", "status": _status(t17, 0.03, 2e-2)},
            {"item": "Tables 18-19 (optimised)", "status": _status(t1819, 0.05, 5e-3)},
        ]
    )
    lines.append(_md_table(summary))

    lines.append("\n## Steady-state simulation vs MATLAB\n")
    lines.append(
        "Selected points of the rate sweeps of Figs. 6-8 (Case 1, `ode45` rtol 1e-3) and "
        "Figs. 11-12 (Case 2, rtol 1e-2). `max_abs_state_difference` is the largest absolute "
        "difference of the scaled state vector `[qL, qG, PI, PO, PW, wL, wG]` "
        "(pressures in units of 100 bar).\n"
    )
    lines.append(
        _md_table(
            sim.drop(columns=[c for c in ("matlab_converged",) if c in sim.columns]),
            "{:.6g}",
        )
    )

    lines.append("\n## Simulation and economics at the published decision variables\n")
    lines.append(
        "The MATLAB optima (`matlab_optimization_results.json`) are re-simulated and the "
        "economic model is evaluated; this reproduces Tables 13, 18 and 19 without running "
        "the optimiser.\n"
    )
    lines.append(_md_table(ref_opt, "{:.6g}"))
    lines.append("\n### Table 13 at the published decision variables\n")
    lines.append(_md_table(_wide(t13_ref), "{:.6g}"))
    lines.append("\n### Tables 18 and 19 at the published decision variables\n")
    lines.append(_md_table(_wide(t18_ref), "{:.6g}"))

    lines.append("\n## Optimisation results\n")
    if not has_opt:
        lines.append("_The run was made with `--no-optimize`; no optimisation results._\n")
    else:
        lines.append(_md_table(opt, "{:.5g}"))
        lines.append("\n### Table 13\n")
        lines.append(_md_table(_wide(t13), "{:.6g}"))
        lines.append("\n### Tables 14 and 15\n")
        lines.append(_md_table(_wide(t1415), "{:.6g}"))
        lines.append("\n### Table 17\n")
        lines.append(_md_table(_wide(t17), "{:.6g}"))
        lines.append("\n### Tables 18 and 19\n")
        lines.append(_md_table(_wide(t1819), "{:.6g}"))
    lines.append(
        "\nDiscrepancies between the paper text and the original code (well geometry, annulus "
        "orientation, discharge-pressure node, cost correlations) are documented in "
        "`docs/reproduction.md` and `docs/errata.md`; the package reproduces the code by "
        "default, so the published numbers are the target of this evaluation.\n"
    )
    report = out / "report.md"
    report.write_text("\n".join(lines))
    if progress:
        progress(f"evaluation written to {out}")
    return report


__all__ = ["REFERENCE_DIR", "evaluate", "load_results"]
