"""pandas tables in the layout of the paper (Tables 1-19).

The input tables (1-12) are built from a :class:`~gaslift.cases.CaseSetup`; the
result tables (13-19) are built from the JSON-like records produced by
:mod:`gaslift.study` (one record per optimisation problem) and show the published
values side by side with the recomputed ones.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from .cases import CaseSetup
from .fluids import FluidProperties
from .units import INCH


# --------------------------------------------------------------------------- input tables
def fluid_table(fluid: FluidProperties) -> pd.DataFrame:
    """Table 1 - properties of the fluids."""
    return pd.DataFrame(
        {
            "Oil": {
                "Density (kg/m³)": fluid.liquid_density,
                "Viscosity (Pa s)": fluid.liquid_viscosity,
                "Interfacial tension (N/m)": fluid.interfacial_tension,
            },
            "Gas": {
                "Reference density (kg/m³)": fluid.gas_reference_density,
                "Normal density (kg/m³)": fluid.gas_normal_density,
                "Reference pressure (bar)": fluid.reference_pressure / 1e5,
                "Polytropic coefficient (-)": fluid.gas_polytropic_coefficient,
                "Viscosity (Pa s)": fluid.gas_viscosity,
            },
        }
    )


def pipes_table(case: CaseSetup) -> pd.DataFrame:
    """Tables 2 and 7 - lines of the network (numbered as in the paper)."""
    rows = []
    for p in case.network.pipes:
        g = p.geometry
        angles = np.degrees(g.segment_angles)
        incl = f"{angles[0]:.0f}°" if np.allclose(angles, angles[0]) else "profile"
        rows.append(
            {
                "Line": case.paper_pipe_numbers.get(p.name, ""),
                "Name": p.name,
                "Type": p.category.value,
                "From": p.start,
                "To": p.end,
                "Diameter (in)": g.inlet_diameter / INCH,
                "Length (m)": g.length,
                "Inclination": incl,
                "Roughness (µm)": g.roughness[0] * 1e6,
            }
        )
    return pd.DataFrame(rows).set_index("Line")


def nodes_table(case: CaseSetup) -> pd.DataFrame:
    """Tables 5 and 11 - nodal description."""
    rows = [
        {
            "Node": case.paper_node_numbers.get(n.name, ""),
            "Name": n.name,
            "Role": n.role.value,
            "X (m)": n.x,
            "Y (m)": n.y,
            "Height (m)": n.z,
        }
        for n in case.network.nodes
    ]
    return pd.DataFrame(rows).set_index("Node")


def profile_table(case: CaseSetup, pipe: str) -> pd.DataFrame:
    """Tables 3, 4, 8, 9 and 10 - vertices of a line (cumulative length, height, x, y)."""
    g = case.network.pipe(pipe).geometry
    return pd.DataFrame(
        {
            "Length (m)": g.cumulative_lengths,
            "Height (m)": g.points[:, 2],
            "x (m)": g.points[:, 0],
            "y (m)": g.points[:, 1],
        }
    ).T


def specifications_table(case: CaseSetup) -> pd.DataFrame:
    """Tables 6 and 12 - nodal specifications."""
    rows = []
    for name, s in case.specifications.items():
        parts = []
        if s.pressure is not None:
            parts.append(f"P = {s.pressure / 1e5:g} bar")
        if s.liquid_flow is not None:
            parts.append(f"wL = {s.liquid_flow:g} kg/s")
        if s.gas_flow is not None:
            parts.append(f"wG = {s.gas_flow:g} kg/s")
        rows.append(
            {
                "Node": case.paper_node_numbers.get(name, ""),
                "Name": name,
                "Specifications": "; ".join(parts),
            }
        )
    return pd.DataFrame(rows).set_index("Node")


# --------------------------------------------------------------------------- result tables
Record = Mapping[str, object]

TABLE13_NAMES = {
    "1": "table13_scenario1",
    "2": "table13_scenario2",
    "3": "table13_scenario3",
    "Max": "fig15_pi1",
}


def _ec(rec: Record, key: str) -> float:
    ec = rec.get("economics") or {}
    return float(ec.get(key, np.nan))


def table13(records: Mapping[str, Record], published: dict) -> pd.DataFrame:
    """Table 13 - Case 1 maximum profit: computed and published columns side by side.

    ``records`` maps the problem names of the study (``table13_scenario1..3`` and
    ``fig15_pi1`` for the "Max" column) to result records.
    """
    pub = published["table13"]
    items = [
        ("Cost ratio (10³ (USD/Nm³)/(USD/bbl))", "cost_ratio_1e3"),
        ("Fuel gas cost (USD/Nm³)", "fuel_gas_cost_usd_m3"),
        ("Oil price (USD/barrel)", "oil_price_usd_bbl"),
        ("Produced oil (kg/s)", "produced_oil_kg_s"),
        ("Injected gas (kg/s)", "injected_gas_kg_s"),
        ("Downhole pressure (bar)", "downhole_pressure_bar"),
        ("Compressor power (kW)", "compressor_power_kW"),
        ("Discharge pressure P_G (bar)", "discharge_pressure_bar"),
        ("Profit (10⁶ USD/yr)", "profit_1e6_usd_yr"),
    ]
    data = {}
    for k, col in enumerate(pub["columns"]):
        rec = records.get(TABLE13_NAMES[col])
        computed = {}
        if rec is not None:
            oil_price = pub["rows"]["oil_price_usd_bbl"][k]
            fuel = _ec(rec, "fuel_gas_cost_1e6_usd_yr")
            vol = _ec(rec, "fuel_gas_1e6_m3_yr")
            gas_price = fuel / vol if vol > 0 else 0.0
            computed = {
                "cost_ratio_1e3": 1e3 * gas_price / oil_price,
                "fuel_gas_cost_usd_m3": gas_price,
                "oil_price_usd_bbl": oil_price,
                "produced_oil_kg_s": float(np.sum(rec["q_oil"])),
                "injected_gas_kg_s": float(np.sum(rec["q_gas"])),
                "downhole_pressure_bar": float(rec["downhole_pressure_bar"][0]),
                "compressor_power_kW": _ec(rec, "compressor_power_kW"),
                "discharge_pressure_bar": _ec(rec, "discharge_pressure_bar"),
                "profit_1e6_usd_yr": _ec(rec, "profit_1e6_usd_yr"),
            }
        data[f"Scenario {col}"] = {label: computed.get(key, np.nan) for label, key in items}
        data[f"Scenario {col} (paper)"] = {label: pub["rows"][key][k] for label, key in items}
    return pd.DataFrame(data)


def table14_15(rec: Record | None, published_table: dict, productivity_index) -> pd.DataFrame:
    """Tables 14 and 15 - Case 2 maximum production (computed and published rows)."""
    pub = published_table["rows"]
    df = pd.DataFrame(
        {
            "Well": [1, 2, 3, 4],
            "Node": pub["node"],
            "Productivity index ((kg/s)/bar)": list(productivity_index),
            "Downhole pressure (bar)": list(rec["downhole_pressure_bar"]) if rec else [np.nan] * 4,
            "Downhole pressure (bar) paper": pub["downhole_pressure_bar"],
            "Produced oil (kg/s)": list(rec["q_oil"]) if rec else [np.nan] * 4,
            "Produced oil (kg/s) paper": pub["produced_oil_kg_s"],
            "Injected gas (kg/s)": list(rec["q_gas"]) if rec else [np.nan] * 4,
            "Injected gas (kg/s) paper": pub["injected_gas_kg_s"],
        }
    ).set_index("Well")
    return df.T


def table17(rec1: Record | None, rec2: Record | None, published: dict) -> pd.DataFrame:
    """Table 17 - decision variables of the optimal design (two scenarios)."""
    pub = published["table17"]["rows"]
    rows = {}
    for i in range(4):
        rows[f"Injected gas QG{i + 1} (kg/s)"] = (
            rec1["q_gas"][i] if rec1 else np.nan,
            pub["injected_gas_kg_s"][0][i],
            rec2["q_gas"][i] if rec2 else np.nan,
            pub["injected_gas_kg_s"][1][i],
        )
    for i in range(4):
        rows[f"Oil flow rate QL{i + 1} (kg/s)"] = (
            rec1["q_oil"][i] if rec1 else np.nan,
            pub["produced_oil_kg_s"][0][i],
            rec2["q_oil"][i] if rec2 else np.nan,
            pub["produced_oil_kg_s"][1][i],
        )
    for i in range(4):
        rows[f"Pipe diameter DG{i + 1} (in)"] = (
            rec1["diameters_in"][i] if rec1 and rec1.get("diameters_in") else np.nan,
            pub["gas_line_diameter_in"][0][i],
            rec2["diameters_in"][i] if rec2 and rec2.get("diameters_in") else np.nan,
            pub["gas_line_diameter_in"][1][i],
        )
    return pd.DataFrame(
        rows,
        index=["Scenario 1", "Scenario 1 (paper)", "Scenario 2", "Scenario 2 (paper)"],
    ).T


def _table18_19(rec1, rec2, published, table, items):
    pub = published[table]["rows"]
    data = {}
    for s, rec in ((0, rec1), (1, rec2)):
        data[f"Scenario {s + 1}"] = {
            label: (_ec(rec, key) if rec else np.nan) for label, key, _u in items
        }
        data[f"Scenario {s + 1} (paper)"] = {label: pub[pkey][s] for label, _k, pkey in items}
    df = pd.DataFrame(data)
    df["Unit"] = [u for _l, _k, u in _units(items)]
    return df


def _units(items):
    return [(label, key, unit) for label, key, unit in items]


def table18(rec1: Record | None, rec2: Record | None, published: dict) -> pd.DataFrame:
    """Table 18 - pipelines, compressor and turbine of the optimal design."""
    items = [
        ("Compressor power wC", "compressor_power_kW", "compressor_power_kW"),
        ("Compressor discharge pressure", "discharge_pressure_bar", "discharge_pressure_bar"),
        ("Compressor capital cost CC", "compressor_cost_1e6_usd", "compressor_cost_1e6_usd"),
        ("Turbine capital cost CT", "turbine_cost_1e6_usd", "turbine_cost_1e6_usd"),
        ("Gas pipelines capital cost CP", "pipeline_cost_1e6_usd", "pipeline_cost_1e6_usd"),
        ("Fuel gas flow rate", "fuel_gas_1e6_m3_yr", "fuel_gas_1e6_m3_yr"),
    ]
    units = ["kW", "bar", "10⁶ USD", "10⁶ USD", "10⁶ USD", "10⁶ m³/yr"]
    df = _table18_19(rec1, rec2, published, "table18", items)
    df["Unit"] = units
    return df


def table19(rec1: Record | None, rec2: Record | None, published: dict) -> pd.DataFrame:
    """Table 19 - optimal design results."""
    items = [
        (
            "Total injected gas flow rate",
            "total_injected_gas_1e6_m3_yr",
            "total_injected_gas_1e6_m3_yr",
        ),
        ("Total oil flow rate", "total_oil_volume_1e6_bbl_yr", "total_oil_volume_1e6_bbl_yr"),
        ("Fuel gas cost AGC", "fuel_gas_cost_1e6_usd_yr", "fuel_gas_cost_1e6_usd_yr"),
        (
            "Annualized capital cost",
            "annualized_capital_cost_1e6_usd_yr",
            "annualized_capital_cost_1e6_usd_yr",
        ),
        ("ISBL", "ISBL_1e6_usd", "ISBL_1e6_usd"),
        ("OSBL", "OSBL_1e6_usd", "OSBL_1e6_usd"),
        ("Revenue R", "revenue_1e6_usd_yr", "revenue_1e6_usd_yr"),
        ("Profit", "profit_1e6_usd_yr", "profit_1e6_usd_yr"),
    ]
    units = [
        "10⁶ m³/yr",
        "10⁶ barrel/yr",
        "10⁶ USD/yr",
        "10⁶ USD/yr",
        "10⁶ USD",
        "10⁶ USD",
        "10⁶ USD/yr",
        "10⁶ USD/yr",
    ]
    df = _table18_19(rec1, rec2, published, "table19", items)
    df["Unit"] = units
    return df


def fig13_table(records: Mapping[str, Record]) -> pd.DataFrame:
    """Fig. 13 data - maximum production versus well pressure (PI = infinity)."""
    rows = []
    for p in (180, 200, 220, 240, 260):
        rec = records.get(f"fig13_{p}bar")
        if rec is None:
            continue
        m = rec.get("matlab") or {}
        rows.append(
            {
                "P_well (bar)": p,
                "oil (kg/s)": float(np.sum(rec["q_oil"])),
                "gas (kg/s)": float(np.sum(rec["q_gas"])),
                "oil MATLAB (kg/s)": float(np.sum(m.get("q_oil", np.nan))),
                "gas MATLAB (kg/s)": float(np.sum(m.get("q_gas", np.nan))),
                "evaluations": rec.get("n_evaluations"),
                "time (s)": rec.get("elapsed_s"),
            }
        )
    return pd.DataFrame(rows).set_index("P_well (bar)")


def fig15_table(records: Mapping[str, Record]) -> pd.DataFrame:
    """Fig. 15 data - maximum production versus productivity index (220 bar)."""
    rows = []
    for pi in (0.5, 1.0, 2.0, 5.0, 10.0):
        rec = records.get(f"fig15_pi{pi:g}")
        if rec is None:
            continue
        m = rec.get("matlab") or {}
        rows.append(
            {
                "PI ((kg/s)/bar)": pi,
                "oil (kg/s)": float(np.sum(rec["q_oil"])),
                "gas (kg/s)": float(np.sum(rec["q_gas"])),
                "P_downhole (bar)": float(rec["downhole_pressure_bar"][0]),
                "oil MATLAB (kg/s)": float(np.sum(m.get("q_oil", np.nan))),
                "gas MATLAB (kg/s)": float(np.sum(m.get("q_gas", np.nan))),
                "P_downhole MATLAB (bar)": float(
                    np.atleast_1d(m.get("downhole_pressure_bar", np.nan))[0]
                ),
                "time (s)": rec.get("elapsed_s"),
            }
        )
    return pd.DataFrame(rows).set_index("PI ((kg/s)/bar)")


__all__ = [
    "TABLE13_NAMES",
    "fig13_table",
    "fig15_table",
    "fluid_table",
    "nodes_table",
    "pipes_table",
    "profile_table",
    "specifications_table",
    "table13",
    "table14_15",
    "table17",
    "table18",
    "table19",
]
