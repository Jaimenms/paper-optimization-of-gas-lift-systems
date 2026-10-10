"""Export the figures and tables of the paper (``gaslift figures``).

Figures 2-12 are regenerated from the model (flow map, topologies, digraphs and rate
sweeps); Figs. 13-16 and Tables 13-19 use the optimisation results written by
``gaslift run`` when they exist, otherwise the published decision variables are
re-simulated (fast) and the figures are labelled accordingly.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import tables
from .cases import CaseSetup, load_case
from .plotting import (
    plot_case_topology,
    plot_digraph,
    plot_flow_pattern_map,
    plot_rate_contours,
)
from .study import load_reference, run_reference_optima, simulation_options
from .sweeps import RateSweep, sweep_rates
from .units import INCH

Progress = Callable[[str], None]


def _save(fig, path: Path, written: list[Path], progress: Progress | None) -> None:
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    written.append(path)
    if progress:
        progress(f"  wrote {path.name}")


def _csv(df: pd.DataFrame, path: Path, written: list[Path], progress: Progress | None) -> None:
    df.to_csv(path)
    written.append(path)
    if progress:
        progress(f"  wrote {path.name}")


def _sweep_frames(sweep: RateSweep, well: int = 0) -> dict[str, pd.DataFrame]:
    idx = pd.Index(sweep.q_oil, name="q_oil_kg_s")
    cols = pd.Index(sweep.q_gas, name="q_gas_kg_s")
    return {
        "downhole": pd.DataFrame(sweep.downhole(well), index=idx, columns=cols),
        "gas_junction": pd.DataFrame(sweep.gas_junction(well), index=idx, columns=cols),
        "gas_inlet": pd.DataFrame(sweep.gas_inlet(well), index=idx, columns=cols),
    }


def _grid(spec: dict, quick: bool) -> tuple[np.ndarray, np.ndarray]:
    q_gas = np.asarray(spec["q_gas"], dtype=float)
    q_oil = np.asarray(spec["q_oil"], dtype=float)
    if quick:
        q_gas = q_gas[::2]
        q_oil = q_oil[::2]
    return q_gas, q_oil


def _case_tables(case: CaseSetup, out: Path, numbers: dict, written, progress) -> None:
    _csv(tables.pipes_table(case), out / numbers["pipes"], written, progress)
    _csv(tables.nodes_table(case), out / numbers["nodes"], written, progress)
    _csv(tables.specifications_table(case), out / numbers["specs"], written, progress)
    for pipe, name in numbers["profiles"].items():
        _csv(tables.profile_table(case, pipe), out / name, written, progress)


def _optimization_records(results_dir: Path | None, case: CaseSetup, progress) -> tuple[dict, str]:
    """Records of the optimisation problems of ``case`` and a label of their origin."""
    if results_dir is not None:
        path = results_dir / case.name / "optimizations.json"
        if path.exists():
            with open(path) as f:
                recs = json.load(f)
            if recs:
                return {r["name"]: r for r in recs}, "SLSQP optimum"
    if progress:
        progress(f"  no optimisation results for {case.name}: re-evaluating the MATLAB optima")
    recs = run_reference_optima(case)
    for r in recs:
        r["matlab"] = {
            "q_oil": r["q_oil"],
            "q_gas": r["q_gas"],
            "diameters_in": r["diameters_in"],
            "profit": -float(r["matlab"]["fob"]),
            "downhole_pressure_bar": np.atleast_1d(r["matlab"]["p_downhole_bar"]).tolist(),
        }
    return {r["name"]: r for r in recs}, "MATLAB optimum (re-evaluated)"


def export_figures(
    out_dir: str | Path,
    *,
    results_dir: str | Path | None = "results",
    quick: bool = False,
    progress: Progress | None = None,
) -> list[Path]:
    """Write the figures (PNG) and tables (CSV) of the paper into ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results = Path(results_dir) if results_dir is not None else None
    published = load_reference("published_tables.json")
    written: list[Path] = []
    t0 = time.perf_counter()
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3})

    case1 = load_case("case1")
    case2 = load_case("case2")

    # --- Table 1 and Fig. 2
    _csv(tables.fluid_table(case1.network.fluid), out / "table01_fluids.csv", written, progress)
    fig, ax = plt.subplots(figsize=(7, 5))
    plot_flow_pattern_map(ax)
    ax.set_title("Fig. 2 - Beggs & Brill horizontal flow pattern map")
    _save(fig, out / "fig02_flow_pattern_map.png", written, progress)

    # --- Case 1: Tables 2-6, Figs. 3-5
    _case_tables(
        case1,
        out,
        {
            "pipes": "table02_case1_lines.csv",
            "nodes": "table05_case1_nodes.csv",
            "specs": "table06_case1_specifications.csv",
            "profiles": {
                "production_line_2": "table03_case1_line2_profile.csv",
                "well_1:gas_line_1": "table04_case1_line8_profile.csv",
            },
        },
        written,
        progress,
    )
    fig = plot_case_topology(case1.network)
    fig.suptitle("Fig. 3 - Case 1, simplified isometric projection", y=1.02)
    _save(fig, out / "fig03_case1_topology.png", written, progress)
    fig, ax = plt.subplots(figsize=(12, 7))
    plot_digraph(
        case1.network,
        ax,
        node_numbers=case1.paper_node_numbers,
        pipe_numbers=case1.paper_pipe_numbers,
    )
    ax.set_title("Figs. 4 and 5 - Case 1 digraph (paper numbering)")
    _save(fig, out / "fig04_05_case1_digraph.png", written, progress)

    # --- Case 1: Figs. 6-8 (rate sweeps)
    rtol1 = float(case1.study.get("integration_rtol", 1e-3))
    sweeps1: dict[str, RateSweep] = {}
    for spec in case1.study["sweeps"]:
        q_gas, q_oil = _grid(spec, quick)
        c = case1.copy()
        if spec.get("gas_line_diameter_in") is not None:
            for w in c.wells:
                c.set_gas_line_diameter(w, float(spec["gas_line_diameter_in"]) * INCH)
        t1 = time.perf_counter()
        sw = sweep_rates(c, q_gas, q_oil, simulation_options(float(spec.get("rtol", rtol1))))
        if progress:
            progress(
                f"  sweep {spec['name']}: {sw.converged.sum()}/{sw.converged.size} points "
                f"converged in {time.perf_counter() - t1:.0f} s"
            )
        sweeps1[spec["name"]] = sw
    sw2in = sweeps1["fig06_fig07_2in"]
    frames = _sweep_frames(sw2in)
    _csv(frames["downhole"], out / "fig06_case1_downhole_pressure.csv", written, progress)
    _csv(frames["gas_junction"], out / "fig07_case1_gas_inlet_pressure_2in.csv", written, progress)
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_rate_contours(
        sw2in.q_gas,
        sw2in.q_oil,
        sw2in.downhole(0),
        ax,
        np.arange(140, 320, 20),
        title="Fig. 6 - Case 1 downhole pressure, node 7 (bar)",
    )
    _save(fig, out / "fig06_case1_downhole_pressure.png", written, progress)
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_rate_contours(
        sw2in.q_gas,
        sw2in.q_oil,
        sw2in.gas_junction(0),
        ax,
        np.arange(200, 560, 50),
        title="Fig. 7 - Case 1 gas lift inlet pressure, node 9 (bar), 2 in gas lines",
    )
    _save(fig, out / "fig07_case1_gas_inlet_pressure_2in.png", written, progress)
    sw3in = sweeps1["fig08_3in"]
    _csv(
        _sweep_frames(sw3in)["gas_junction"],
        out / "fig08_case1_gas_inlet_pressure_3in.csv",
        written,
        progress,
    )
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_rate_contours(
        sw3in.q_gas,
        sw3in.q_oil,
        sw3in.gas_junction(0),
        ax,
        np.arange(180, 320, 20),
        title="Fig. 8 - Case 1 gas lift inlet pressure, node 9 (bar), 3 in gas lines",
    )
    _save(fig, out / "fig08_case1_gas_inlet_pressure_3in.png", written, progress)

    # --- Case 2: Tables 7-12, Figs. 9-12
    _case_tables(
        case2,
        out,
        {
            "pipes": "table07_case2_lines.csv",
            "nodes": "table11_case2_nodes.csv",
            "specs": "table12_case2_specifications.csv",
            "profiles": {
                "production_line_2": "table08_case2_line2_profile.csv",
                "well_1:gas_line_1": "table09_case2_lines_8_14_profile.csv",
                "well_3:gas_line_1": "table10_case2_lines_20_26_profile.csv",
            },
        },
        written,
        progress,
    )
    fig = plot_case_topology(case2.network)
    fig.suptitle("Fig. 9 - Case 2, simplified isometric projection", y=1.02)
    _save(fig, out / "fig09_case2_topology.png", written, progress)
    fig, ax = plt.subplots(figsize=(16, 8))
    plot_digraph(
        case2.network,
        ax,
        node_numbers=case2.paper_node_numbers,
        pipe_numbers=case2.paper_pipe_numbers,
    )
    ax.set_title("Fig. 10 - Case 2 digraph (paper numbering)")
    _save(fig, out / "fig10_case2_digraph.png", written, progress)
    rtol2 = float(case2.study.get("integration_rtol", 1e-2))
    spec = case2.study["sweeps"][0]
    q_gas, q_oil = _grid(spec, quick)
    t1 = time.perf_counter()
    sw = sweep_rates(case2.copy(), q_gas, q_oil, simulation_options(rtol2))
    if progress:
        progress(
            f"  sweep {spec['name']}: {sw.converged.sum()}/{sw.converged.size} points "
            f"converged in {time.perf_counter() - t1:.0f} s"
        )
    levels = np.arange(140, 380, 20)
    for well, fig_no, nodes in ((0, 11, "7 and 13"), (2, 12, "19 and 25")):
        _csv(
            _sweep_frames(sw, well)["downhole"],
            out / f"fig{fig_no}_case2_downhole_pressure_nodes_{nodes.replace(' and ', '_')}.csv",
            written,
            progress,
        )
        fig, ax = plt.subplots(figsize=(7, 6))
        plot_rate_contours(
            sw.q_gas,
            sw.q_oil,
            sw.downhole(well),
            ax,
            levels,
            title=f"Fig. {fig_no} - Case 2 downhole pressure, nodes {nodes} (bar)",
        )
        _save(
            fig,
            out / f"fig{fig_no}_case2_downhole_pressure_nodes_{nodes.replace(' and ', '_')}.png",
            written,
            progress,
        )

    # --- Figs. 13-16 and Tables 13-19
    recs1, label1 = _optimization_records(results, case1, progress)
    recs2, label2 = _optimization_records(results, case2, progress)

    f13 = tables.fig13_table(recs1)
    if not f13.empty:
        _csv(f13, out / "fig13_case1_max_production_vs_pressure.csv", written, progress)
        fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        axes[0].plot(f13.index, f13["oil (kg/s)"], "ko-", label=f"gaslift ({label1})")
        axes[0].plot(f13.index, f13["oil MATLAB (kg/s)"], "r+", ms=10, label="MATLAB (paper)")
        axes[0].set_ylabel("Produced oil (kg/s)")
        axes[0].legend()
        axes[0].set_title("Fig. 13 - Case 1 maximum production for different well pressures")
        axes[1].plot(f13.index, f13["gas (kg/s)"], "ko-")
        axes[1].plot(f13.index, f13["gas MATLAB (kg/s)"], "r+", ms=10)
        axes[1].set_ylabel("Injected gas (kg/s)")
        axes[1].set_xlabel("Downhole pressure (bar)")
        _save(fig, out / "fig13_case1_max_production_vs_pressure.png", written, progress)

        fig, ax = plt.subplots(figsize=(7, 6))
        plot_rate_contours(
            sw2in.q_gas,
            sw2in.q_oil,
            sw2in.downhole(0),
            ax,
            np.arange(140, 320, 20),
            title="Fig. 14 - Case 1 maximum production locus",
        )
        ax.plot(f13["gas (kg/s)"], f13["oil (kg/s)"], "ro--", label=f"optimum locus ({label1})")
        ax.legend(loc="upper right")
        _save(fig, out / "fig14_case1_max_production_locus.png", written, progress)

    f15 = tables.fig15_table(recs1)
    if not f15.empty:
        _csv(f15, out / "fig15_case1_max_production_vs_pi.csv", written, progress)
        fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
        for ax, col, ref, lab in zip(
            axes,
            ["oil (kg/s)", "gas (kg/s)", "P_downhole (bar)"],
            ["oil MATLAB (kg/s)", "gas MATLAB (kg/s)", "P_downhole MATLAB (bar)"],
            ["Produced oil (kg/s)", "Injected gas (kg/s)", "Well pressure (bar)"],
            strict=True,
        ):
            ax.plot(f15.index, f15[col], "ko-", label=f"gaslift ({label1})")
            ax.plot(f15.index, f15[ref], "r+", ms=10, label="MATLAB (paper)")
            ax.set_ylabel(lab)
        axes[0].legend()
        axes[0].set_title("Fig. 15 - Case 1 maximum production vs productivity index (220 bar)")
        axes[-1].set_xlabel("Productivity index ((kg/s)/bar)")
        _save(fig, out / "fig15_case1_max_production_vs_pi.png", written, progress)

    t13 = tables.table13(recs1, published)
    _csv(t13, out / "table13_case1_max_profit.csv", written, progress)
    cols = [c for c in t13.columns if "(paper)" not in c]
    ratio = t13.loc["Cost ratio (10³ (USD/Nm³)/(USD/bbl))", cols].astype(float) / 1e3
    if ratio.notna().any():
        fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
        for ax, row, lab in zip(
            axes,
            ["Produced oil (kg/s)", "Injected gas (kg/s)", "Profit (10⁶ USD/yr)"],
            ["Produced oil (kg/s)", "Injected gas (kg/s)", "Profit (10⁶ USD/yr)"],
            strict=True,
        ):
            ax.plot(ratio, t13.loc[row, cols].astype(float), "ko-", label=f"gaslift ({label1})")
            pcols = [c + " (paper)" for c in cols]
            ax.plot(ratio, t13.loc[row, pcols].astype(float), "r+", ms=10, label="paper (Table 13)")
            ax.set_ylabel(lab)
        axes[0].legend()
        axes[0].set_title("Fig. 16 - Case 1 maximum profit for different price scenarios")
        axes[-1].set_xlabel("Cost ratio ((USD/Nm³)/(USD/barrel))")
        _save(fig, out / "fig16_case1_max_profit.png", written, progress)

    pi = [1.0, 2.5, 1.0, 0.5]
    _csv(
        tables.table14_15(recs2.get("table14_maxprod_220bar"), published["table14"], pi),
        out / "table14_case2_max_production_220bar.csv",
        written,
        progress,
    )
    _csv(
        tables.table14_15(recs2.get("table15_maxprod_240bar"), published["table15"], pi),
        out / "table15_case2_max_production_240bar.csv",
        written,
        progress,
    )
    d1, d2 = recs2.get("table17_19_design_220bar"), recs2.get("table17_19_design_240bar")
    _csv(
        tables.table17(d1, d2, published),
        out / "table17_case2_design_variables.csv",
        written,
        progress,
    )
    _csv(
        tables.table18(d1, d2, published),
        out / "table18_case2_design_equipment.csv",
        written,
        progress,
    )
    _csv(
        tables.table19(d1, d2, published),
        out / "table19_case2_design_results.csv",
        written,
        progress,
    )
    (out / "README.txt").write_text(
        "Figures and tables of de Souza et al. (2010) exported by `gaslift figures`.\n"
        f"Figs. 13-16 and Tables 13-19 of Case 1 use: {label1}; Case 2: {label2}.\n"
        f"quick = {quick}; generated in {time.perf_counter() - t0:.0f} s.\n"
    )
    if progress:
        progress(f"{len(written)} files written to {out} in {time.perf_counter() - t0:.0f} s")
    return written


__all__ = ["export_figures"]
