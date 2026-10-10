"""Paper studies attached to the cases: operating points, rate sweeps and optimisations.

A *study* is a JSON-serialisable dictionary stored in the ``"study"`` block of a
case file (``cases/*.json``). :func:`default_study` builds the studies of the two
paper cases; :func:`run_study` executes one and writes its results as JSON files
(``operating_points.json``, ``reference_optima.json``, ``optimizations.json``),
which :mod:`gaslift.evaluation` compares with the MATLAB reference results and
with the published tables.

Study blocks::

    {
      "integration_rtol": 1e-3,          # ode45 relative tolerance of the MATLAB runs
      "operating_points": [ {name, q_oil, q_gas, gas_line_diameter_in?, reference?} ],
      "sweeps":           [ {name, paper, q_oil, q_gas, gas_line_diameter_in?, rtol?, reference?} ],
      "optimizations":    [ {name, paper, objective, reservoirs, cost?, oil_bounds, gas_bounds,
                             diameter_bounds?, x0, reference?, published?} ]
    }

``reference`` names an entry of ``data/reference`` (a grid file plus a grid index for
operating points and sweeps, a key of ``matlab_optimization_results.json`` for
optimisations). ``x0`` is either ``{"q_oil", "q_gas", "diameters_in"?}`` or
``{"from": "<optimization name>", "diameters_in": ...}`` to start from the optimum
of a previous problem (as the paper does for the design problem).
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from . import __version__
from .cases import CaseSetup
from .economics import CompressorParameters, CostModel, EconomicEvaluation, evaluate_economics
from .model import ModelOptions
from .optimization import GasLiftProblem, Objective, OptimizationResult, Reservoir
from .pipe_model import IntegrationOptions
from .results import NetworkSolution
from .simulation import SimulationOptions
from .units import BAR, INCH

REFERENCE_DIR = Path(__file__).resolve().parents[2] / "data" / "reference"

Progress = Callable[[str], None]


# --------------------------------------------------------------------------- study definitions
def _grid_values(start: float, stop: float, step: float) -> list[float]:
    return [float(v) for v in np.arange(start, stop + 0.5 * step, step)]


def _case1_study() -> dict:
    q_gas = _grid_values(0.5, 4.0, 0.5)
    q_oil = [0.5, *_grid_values(1.0, 20.0, 1.0)]
    points = [
        {
            "name": "text_section_4.1",
            "paper": "Section 4.1 text (reading of Fig. 6)",
            "q_oil": 10.0,
            "q_gas": 2.0,
        }
    ]
    # The MATLAB sweeps of Figs. 6-8 were run with ode45 rtol = 1e-3; every other result
    # of the paper (optimisations, Case 2) with rtol = 1e-2.
    for i, j in [(0, 0), (2, 1), (5, 5), (10, 3), (15, 7), (20, 2)]:
        points.append(
            {
                "name": f"matlab_grid_2in_{i}_{j}",
                "paper": "Fig. 6 / Fig. 7 sweep point",
                "q_oil": q_oil[i],
                "q_gas": q_gas[j],
                "rtol": 1e-3,
                "reference": {"grid": "case1_grid_2in.json", "index": [i, j]},
            }
        )
    for i, j in [(10, 3), (4, 6)]:
        points.append(
            {
                "name": f"matlab_grid_3in_{i}_{j}",
                "paper": "Fig. 8 sweep point",
                "q_oil": q_oil[i],
                "q_gas": q_gas[j],
                "gas_line_diameter_in": 3.0,
                "rtol": 1e-3,
                "reference": {"grid": "case1_grid_3in.json", "index": [i, j]},
            }
        )
    sweeps = [
        {
            "name": "fig06_fig07_2in",
            "paper": ["Fig. 6", "Fig. 7"],
            "q_gas": q_gas,
            "q_oil": q_oil,
            "rtol": 1e-3,
            "reference": "case1_grid_2in.json",
        },
        {
            "name": "fig08_3in",
            "paper": ["Fig. 8"],
            "q_gas": q_gas,
            "q_oil": q_oil,
            "gas_line_diameter_in": 3.0,
            "rtol": 1e-3,
            "reference": "case1_grid_3in.json",
        },
    ]
    bounds = {"oil_bounds": [5.0, 20.0], "gas_bounds": [1.0, 4.0]}
    x0 = {"q_oil": 10.0, "q_gas": 2.0}
    opts: list[dict] = []
    for p in (180, 200, 220, 240, 260):
        opts.append(
            {
                "name": f"fig13_{p}bar",
                "paper": ["Fig. 13", "Fig. 14"],
                "objective": "max_production",
                "reservoirs": [{"pressure": float(p), "productivity_index": None}],
                **bounds,
                "x0": x0,
                "reference": f"case1_maxprod_piinf_{p}bar",
            }
        )
    for pi in (0.5, 1.0, 2.0, 5.0, 10.0):
        opts.append(
            {
                "name": f"fig15_pi{pi:g}",
                "paper": ["Fig. 15"] + (["Table 13 (Max)"] if pi == 1.0 else []),
                "objective": "max_production",
                "reservoirs": [{"pressure": 220.0, "productivity_index": pi}],
                **bounds,
                "x0": x0,
                "reference": f"case1_maxprod_pi{pi:g}_220bar",
                **({"published": {"table": "table13", "column": "Max"}} if pi == 1.0 else {}),
            }
        )
    for scenario, price, key in (
        ("1", 51.15 / 4, "case1_maxprofit_oil12.79"),
        ("2", 51.15 / 2, "case1_maxprofit_oil25.58"),
        ("3", 51.15, "case1_maxprofit_oil51.15"),
    ):
        opts.append(
            {
                "name": f"table13_scenario{scenario}",
                "paper": ["Table 13", "Fig. 16"],
                "objective": "max_profit",
                "cost": {"oil_price": price},
                "reservoirs": [{"pressure": 220.0, "productivity_index": 1.0}],
                **bounds,
                "x0": x0,
                "reference": key,
                "published": {"table": "table13", "column": scenario},
            }
        )
    return {
        "integration_rtol": 1e-2,
        "operating_points": points,
        "sweeps": sweeps,
        "optimizations": opts,
    }


def _case2_study() -> dict:
    q_gas = _grid_values(0.25, 4.0, 0.25)
    q_oil = _grid_values(0.5, 8.0, 0.5)
    points = [
        {
            "name": "notebook_point",
            "paper": "Section 4.2 (3 kg/s of oil and 1 kg/s of gas per well)",
            "q_oil": 3.0,
            "q_gas": 1.0,
        }
    ]
    for i, j in [(0, 3), (5, 5), (8, 2), (12, 6), (6, 14), (15, 3)]:
        points.append(
            {
                "name": f"matlab_grid_{i}_{j}",
                "paper": "Fig. 11 / Fig. 12 sweep point",
                "q_oil": q_oil[i],
                "q_gas": q_gas[j],
                "reference": {"grid": "case2_grid.json", "index": [i, j]},
            }
        )
    sweeps = [
        {
            "name": "fig11_fig12",
            "paper": ["Fig. 11", "Fig. 12"],
            "q_gas": q_gas,
            "q_oil": q_oil,
            "reference": "case2_grid.json",
        }
    ]
    pi = [1.0, 2.5, 1.0, 0.5]
    opts = []
    for p, table in ((220.0, "table14"), (240.0, "table15")):
        opts.append(
            {
                "name": f"{table}_maxprod_{p:g}bar",
                "paper": [f"Table {table[-2:]}"],
                "objective": "max_production",
                "reservoirs": [{"pressure": p, "productivity_index": v} for v in pi],
                "oil_bounds": [0.2, 20.0],
                "gas_bounds": [0.2, 10.0],
                "x0": {"q_oil": 3.0, "q_gas": 1.0},
                "reference": f"case2_maxprod_{p:g}bar",
                "published": {"table": table},
            }
        )
    for scenario, p in (("1", 220.0), ("2", 240.0)):
        opts.append(
            {
                "name": f"table17_19_design_{p:g}bar",
                "paper": ["Table 17", "Table 18", "Table 19"],
                "objective": "optimal_design",
                "reservoirs": [{"pressure": p, "productivity_index": v} for v in pi],
                "oil_bounds": [0.5, 20.0],
                "gas_bounds": [0.1, 4.0],
                "diameter_bounds": [1.0, 4.0],
                "x0": {
                    "from": f"table1{4 if p == 220 else 5}_maxprod_{p:g}bar",
                    "diameters_in": 2.0,
                },
                "reference": f"case2_design_{p:g}bar",
                "published": {"table": "table17", "column": scenario},
            }
        )
    return {
        "integration_rtol": 1e-2,
        "operating_points": points,
        "sweeps": sweeps,
        "optimizations": opts,
    }


def default_study(case_name: str) -> dict:
    """Study block of a paper case (``"case1"`` or ``"case2"``)."""
    if case_name == "case1":
        return _case1_study()
    if case_name == "case2":
        return _case2_study()
    raise KeyError(f"no default study for case {case_name!r}")


# --------------------------------------------------------------------------- helpers
def simulation_options(rtol: float, profiles: bool = False) -> SimulationOptions:
    return SimulationOptions(
        model=ModelOptions(integration=IntegrationOptions(rtol=rtol)), compute_profiles=profiles
    )


def load_reference(name: str) -> dict:
    with open(REFERENCE_DIR / name) as f:
        return json.load(f)


def _reservoirs(spec: dict) -> list[Reservoir]:
    out = []
    for r in spec["reservoirs"]:
        pi = r.get("productivity_index")
        out.append(Reservoir(float(r["pressure"]), np.inf if pi is None else float(pi)))
    return out


def _cost(spec: dict) -> CostModel:
    return CostModel(**spec.get("cost", {}))


def _compressor(spec: dict) -> CompressorParameters:
    return CompressorParameters(**spec.get("compressor", {}))


def build_problem(case: CaseSetup, spec: dict) -> GasLiftProblem:
    """Instantiate the :class:`GasLiftProblem` described by an optimisation spec."""
    kwargs: dict[str, Any] = {}
    for key in ("oil_bounds", "gas_bounds", "diameter_bounds"):
        if key in spec:
            kwargs[key] = tuple(float(v) for v in spec[key])
    rtol = float(spec.get("rtol", case.study.get("integration_rtol", 1e-2)))
    return GasLiftProblem(
        case.copy(),
        _reservoirs(spec),
        Objective(spec["objective"]),
        cost=_cost(spec),
        compressor=_compressor(spec),
        simulation=simulation_options(rtol),
        **kwargs,
    )


def _case_with_diameters(case: CaseSetup, diameters_in) -> CaseSetup:
    c = case.copy()
    if diameters_in is not None:
        d = np.broadcast_to(np.asarray(diameters_in, dtype=float), (c.n_wells,))
        for w, di in zip(c.wells, d, strict=True):
            c.set_gas_line_diameter(w, float(di) * INCH)
    return c


def _well_pressures(sol: NetworkSolution, case: CaseSetup) -> dict:
    return {
        "downhole_pressure_bar": (sol.well_downhole_pressures(case.wells) / BAR).tolist(),
        "injection_point_pressure_bar": [sol.pressure_bar(w.injection_point) for w in case.wells],
        "wellhead_pressure_bar": [sol.pressure_bar(w.wellhead) for w in case.wells],
        "gas_junction_pressure_bar": [sol.pressure_bar(w.gas_junction) for w in case.wells],
        "gas_inlet_pressure_bar": [sol.pressure_bar(w.gas_inlet) for w in case.wells],
    }


def _economics_dict(ec: EconomicEvaluation) -> dict:
    d = {k: float(v) for k, v in ec.to_series().items()}
    d["oil_rate_kg_s"] = ec.oil_rate.tolist()
    d["gas_rate_kg_s"] = ec.gas_rate.tolist()
    d["oil_volume_1e6_bbl_yr"] = ec.oil_volume.tolist()
    d["injected_gas_volume_1e6_m3_yr"] = ec.injected_gas_volume.tolist()
    d["pipeline_diameters_in"] = ec.pipeline_diameters_in.tolist()
    d["pipeline_lengths_m"] = ec.pipeline_lengths_m.tolist()
    return d


def _cost_for_objective(spec: dict) -> CostModel:
    """Cost switches used by the objective of ``spec`` (same rule as GasLiftProblem)."""
    from dataclasses import replace

    cost = _cost(spec)
    obj = Objective(spec["objective"])
    if obj is Objective.MAX_PRODUCTION:
        return replace(
            cost,
            include_fuel_gas_cost=False,
            include_capital_cost=False,
            include_injected_gas_cost=False,
        )
    if obj is Objective.MAX_PROFIT:
        return replace(cost, include_fuel_gas_cost=True, include_capital_cost=False)
    return replace(cost, include_fuel_gas_cost=True, include_capital_cost=True)


# --------------------------------------------------------------------------- runners
def solve_point(case: CaseSetup, q_oil, q_gas, options: SimulationOptions) -> NetworkSolution:
    """Solve ``case`` at the given rates; on failure, continue from lower rates.

    The linear initial estimate is not always good enough far from the paper's
    operating region (very high gas rates); the rates are then ramped up in a few
    steps, warm-starting each solve from the previous one (as a rate sweep does).
    """
    sim = case.simulator(options)
    case.set_well_rates(q_oil, q_gas)
    sol = sim.solve()
    if sol.converged:
        return sol
    q_oil = np.asarray(q_oil, dtype=float)
    q_gas = np.asarray(q_gas, dtype=float)
    sim.reset()
    for f in (0.25, 0.5, 0.75, 1.0):
        case.set_well_rates(np.maximum(f * q_oil, 0.5), np.maximum(f * q_gas, 0.25))
        sol = sim.solve()
        if not sol.converged:
            break
    case.set_well_rates(q_oil, q_gas)
    if sol.converged and f == 1.0:
        return sol
    return sim.solve(warm_start=sol.converged)


def run_operating_points(case: CaseSetup, progress: Progress | None = None) -> list[dict]:
    """Steady-state simulations at the operating points of the study."""
    rtol = float(case.study.get("integration_rtol", 1e-2))
    records = []
    for spec in case.study.get("operating_points", []):
        c = _case_with_diameters(case, spec.get("gas_line_diameter_in"))
        t0 = time.perf_counter()
        sol = solve_point(
            c, spec["q_oil"], spec["q_gas"], simulation_options(float(spec.get("rtol", rtol)))
        )
        rec = {
            "name": spec["name"],
            "paper": spec.get("paper", ""),
            "q_oil": spec["q_oil"],
            "q_gas": spec["q_gas"],
            "gas_line_diameter_in": spec.get("gas_line_diameter_in"),
            "reference": spec.get("reference"),
            "converged": bool(sol.converged),
            "iterations": int(sol.status.iterations),
            "residual_norm": float(sol.status.residual_norm),
            "elapsed_s": time.perf_counter() - t0,
            "state_vector": sol.x.tolist(),
            "node_pressure_bar": (sol.node_pressure / BAR).tolist(),
        }
        if sol.converged:
            rec.update(_well_pressures(sol, c))
        records.append(rec)
        if progress:
            p = rec.get("downhole_pressure_bar", [np.nan])[0]
            progress(
                f"  {spec['name']}: {sol.status}; downhole {p:.2f} bar"
                if sol.converged
                else f"  {spec['name']}: {sol.status}"
            )
    return records


def run_reference_optima(case: CaseSetup, progress: Progress | None = None) -> list[dict]:
    """Simulate and evaluate the economics at the published (MATLAB) decision variables.

    This reproduces Tables 13, 18 and 19 without running the optimiser.
    """
    matlab = load_reference("matlab_optimization_results.json")
    rtol = float(case.study.get("integration_rtol", 1e-2))
    records = []
    for spec in case.study.get("optimizations", []):
        key = spec.get("reference")
        if key is None or key not in matlab:
            continue
        ref = matlab[key]
        n = case.n_wells
        xopt = np.atleast_1d(np.asarray(ref["xopt"], dtype=float))
        q_oil, q_gas = xopt[:n], xopt[n : 2 * n]
        d_in = xopt[2 * n : 3 * n] if xopt.size >= 3 * n else None
        c = _case_with_diameters(case, d_in)
        c.set_well_rates(q_oil, q_gas)
        t0 = time.perf_counter()
        sol = c.simulator(simulation_options(float(spec.get("rtol", rtol)))).solve()
        rec = {
            "name": spec["name"],
            "paper": spec.get("paper", []),
            "objective": spec["objective"],
            "reference": key,
            "published": spec.get("published"),
            "q_oil": q_oil.tolist(),
            "q_gas": q_gas.tolist(),
            "diameters_in": None if d_in is None else d_in.tolist(),
            "converged": bool(sol.converged),
            "elapsed_s": time.perf_counter() - t0,
            "matlab": {k: v for k, v in ref.items() if k not in ("state_vector",)},
            "state_vector": sol.x.tolist(),
        }
        if sol.converged:
            ec = evaluate_economics(sol, c.wells, _cost_for_objective(spec), _compressor(spec))
            rec.update(_well_pressures(sol, c))
            rec["economics"] = _economics_dict(ec)
            x_ref = np.asarray(ref["state_vector"], dtype=float)
            rec["max_abs_state_difference"] = float(np.max(np.abs(sol.x - x_ref)))
        records.append(rec)
        if progress:
            msg = f"  {spec['name']} @ MATLAB optimum: {sol.status}"
            if sol.converged:
                msg += (
                    f"; profit {rec['economics']['profit_1e6_usd_yr']:.3f} "
                    f"(MATLAB {-ref['fob']:.3f}), max|dx| {rec['max_abs_state_difference']:.1e}"
                )
            progress(msg)
    return records


def optimization_record(
    name: str,
    spec: dict,
    res: OptimizationResult,
    elapsed: float | None = None,
    matlab_ref: dict | None = None,
) -> dict:
    """JSON-serialisable record of an :class:`OptimizationResult` (``run`` and notebooks)."""
    n = res.problem.n_wells
    rec = {
        "name": name,
        "paper": spec.get("paper", []),
        "objective": spec["objective"],
        "reference": spec.get("reference"),
        "published": spec.get("published"),
        "reservoirs": spec["reservoirs"],
        "success": bool(res.success),
        "message": res.message,
        "n_evaluations": int(res.n_evaluations),
        "n_iterations": int(res.n_iterations),
        "elapsed_s": elapsed,
        "q_oil": res.oil_rates.tolist(),
        "q_gas": res.gas_rates.tolist(),
        "diameters_in": None if res.diameters_in is None else res.diameters_in.tolist(),
        "objective_value": float(res.objective),
        "economics": _economics_dict(res.economics),
        **_well_pressures(res.solution, res.problem.case),
    }
    if matlab_ref is not None:
        xopt = np.atleast_1d(np.asarray(matlab_ref["xopt"], dtype=float))
        rec["matlab"] = {
            "q_oil": xopt[:n].tolist(),
            "q_gas": xopt[n : 2 * n].tolist(),
            "diameters_in": xopt[2 * n : 3 * n].tolist() if xopt.size >= 3 * n else None,
            "profit": -float(matlab_ref["fob"]),
            "downhole_pressure_bar": np.atleast_1d(matlab_ref["p_downhole_bar"]).tolist(),
        }
    return rec


def initial_vector(prob: GasLiftProblem, spec: dict, previous: dict[str, dict]) -> np.ndarray:
    """Initial decision vector of ``spec`` (explicit, or taken from a previous optimum)."""
    x0_spec = spec["x0"]
    if "from" in x0_spec:
        start = previous.get(x0_spec["from"])
        if start is None:
            raise ValueError(f"{spec['name']}: initial point {x0_spec['from']!r} was not computed")
        return prob.initial_vector(start["q_oil"], start["q_gas"], x0_spec.get("diameters_in"))
    return prob.initial_vector(x0_spec["q_oil"], x0_spec["q_gas"], x0_spec.get("diameters_in"))


def run_optimizations(
    case: CaseSetup,
    progress: Progress | None = None,
    only: list[str] | None = None,
) -> list[dict]:
    """Solve the optimisation problems of the study (Figs. 13-16, Tables 13-19)."""
    matlab = load_reference("matlab_optimization_results.json")
    records: list[dict] = []
    by_name: dict[str, dict] = {}
    for spec in case.study.get("optimizations", []):
        name = spec["name"]
        if only and not any(name.startswith(p) for p in only):
            continue
        prob = build_problem(case, spec)
        x0 = initial_vector(prob, spec, by_name)
        t0 = time.perf_counter()
        res = prob.solve(x0)
        elapsed = time.perf_counter() - t0
        ref = matlab.get(spec.get("reference", ""), None)
        rec = optimization_record(name, spec, res, elapsed, ref)
        records.append(rec)
        by_name[name] = rec
        if progress:
            ec = res.economics
            msg = (
                f"  {name}: {res.n_evaluations} simulations in {elapsed:.0f} s; "
                f"oil {np.sum(res.oil_rates):.3f} kg/s, gas {np.sum(res.gas_rates):.3f} kg/s, "
                f"profit {ec.profit:.3f} 10^6 USD/yr"
            )
            if ref is not None:
                msg += f" (MATLAB {-ref['fob']:.3f})"
            progress(msg)
    return records


def run_study(
    case: CaseSetup,
    out_dir: str | Path,
    *,
    optimize: bool = True,
    only: list[str] | None = None,
    progress: Progress | None = None,
) -> dict:
    """Run the study of ``case`` and write the JSON results into ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    if progress:
        progress(f"{case.name}: operating points")
    points = run_operating_points(case, progress)
    (out / "operating_points.json").write_text(json.dumps(points, indent=1))
    if progress:
        progress(f"{case.name}: economics at the published decision variables")
    reference = run_reference_optima(case, progress)
    (out / "reference_optima.json").write_text(json.dumps(reference, indent=1))
    optimizations: list[dict] = []
    if optimize:
        if progress:
            progress(f"{case.name}: optimisation problems")
        optimizations = run_optimizations(case, progress, only)
        (out / "optimizations.json").write_text(json.dumps(optimizations, indent=1))
    summary = {
        "case": case.name,
        "description": case.description,
        "gaslift_version": __version__,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": time.perf_counter() - t0,
        "optimize": optimize,
        "n_operating_points": len(points),
        "n_reference_optima": len(reference),
        "n_optimizations": len(optimizations),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary
