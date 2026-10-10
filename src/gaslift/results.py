"""Containers for simulation results with convenient, unit-aware accessors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .beggs_brill import PATTERN_NAMES
from .network import Network, Specifications, StateLayout, Well
from .solver import SolverStatus
from .units import BAR


@dataclass
class PipeProfiles:
    """Flow state at equally spaced points along every pipe (arrays of shape ``(S, n)``)."""

    positions: np.ndarray
    """Normalised axial positions in ``[0, 1]``."""
    axial_position: np.ndarray
    pressure: np.ndarray
    liquid_holdup: np.ndarray
    gas_holdup: np.ndarray
    no_slip_holdup: np.ndarray
    liquid_velocity: np.ndarray
    gas_velocity: np.ndarray
    superficial_liquid_velocity: np.ndarray
    superficial_gas_velocity: np.ndarray
    gas_density: np.ndarray
    mixture_density: np.ndarray
    pattern_index: np.ndarray
    angle: np.ndarray
    diameter: np.ndarray
    froude: np.ndarray
    pressure_gradient: np.ndarray
    """``dP/dx`` with respect to the physical axial coordinate, Pa/m."""


class NetworkSolution:
    """Solution of the two-phase network model for one set of specifications."""

    def __init__(
        self,
        network: Network,
        specifications: Specifications,
        x: np.ndarray,
        status: SolverStatus,
        pressure_scale: float,
        profiles: PipeProfiles | None = None,
    ) -> None:
        self.network = network
        self.specifications = specifications
        self.x = np.asarray(x, dtype=float)
        self.status = status
        self.pressure_scale = pressure_scale
        self.profiles = profiles
        self.layout: StateLayout = network.layout

    # ------------------------------------------------------------------ raw vectors
    @property
    def converged(self) -> bool:
        return self.status.converged

    @property
    def pipe_liquid_flow(self) -> np.ndarray:
        """Internal liquid mass flow rate of each pipe, kg/s."""
        return self.x[self.layout.qL]

    @property
    def pipe_gas_flow(self) -> np.ndarray:
        return self.x[self.layout.qG]

    @property
    def pipe_inlet_pressure(self) -> np.ndarray:
        """Pressure at the inlet of each pipe, Pa."""
        return self.x[self.layout.PI] * self.pressure_scale

    @property
    def pipe_outlet_pressure(self) -> np.ndarray:
        return self.x[self.layout.PO] * self.pressure_scale

    @property
    def node_pressure(self) -> np.ndarray:
        """Stagnant pressure at each node, Pa."""
        return self.x[self.layout.PW] * self.pressure_scale

    @property
    def node_liquid_flow(self) -> np.ndarray:
        """External liquid mass flow rate at each node, kg/s (positive entering)."""
        return self.x[self.layout.wL]

    @property
    def node_gas_flow(self) -> np.ndarray:
        return self.x[self.layout.wG]

    # ------------------------------------------------------------------ by name
    def pressure(self, node: str) -> float:
        """Node pressure, Pa."""
        return float(self.node_pressure[self.network.node_index(node)])

    def pressure_bar(self, node: str) -> float:
        return self.pressure(node) / BAR

    def liquid_inflow(self, node: str) -> float:
        return float(self.node_liquid_flow[self.network.node_index(node)])

    def gas_inflow(self, node: str) -> float:
        return float(self.node_gas_flow[self.network.node_index(node)])

    def pipe_flows(self, pipe: str) -> tuple[float, float]:
        j = self.network.pipe_index(pipe)
        return float(self.pipe_liquid_flow[j]), float(self.pipe_gas_flow[j])

    # ------------------------------------------------------------------ wells
    def well_downhole_pressures(
        self, wells: list[Well] | tuple[Well, ...] | None = None
    ) -> np.ndarray:
        """Downhole pressure of each well, Pa."""
        wells = self.network.wells if wells is None else wells
        return np.array([self.pressure(w.downhole) for w in wells])

    def well_gas_inlet_pressures(self, wells=None) -> np.ndarray:
        """Pressure at the gas injection point on the platform of each well, Pa."""
        wells = self.network.wells if wells is None else wells
        return np.array([self.pressure(w.gas_inlet) for w in wells])

    def well_oil_rates(self, wells=None) -> np.ndarray:
        wells = self.network.wells if wells is None else wells
        return np.array([self.liquid_inflow(w.downhole) for w in wells])

    def well_gas_rates(self, wells=None) -> np.ndarray:
        wells = self.network.wells if wells is None else wells
        return np.array([self.gas_inflow(w.gas_inlet) for w in wells])

    # ------------------------------------------------------------------ tables
    def nodes_table(self) -> pd.DataFrame:
        net = self.network
        return pd.DataFrame(
            {
                "node": net.node_names,
                "role": [n.role.value for n in net.nodes],
                "x_m": [n.x for n in net.nodes],
                "y_m": [n.y for n in net.nodes],
                "z_m": [n.z for n in net.nodes],
                "pressure_bar": self.node_pressure / BAR,
                "liquid_inflow_kg_s": self.node_liquid_flow,
                "gas_inflow_kg_s": self.node_gas_flow,
            }
        ).set_index("node")

    def pipes_table(self) -> pd.DataFrame:
        net = self.network
        data = {
            "pipe": net.pipe_names,
            "category": [p.category.value for p in net.pipes],
            "from": [p.start for p in net.pipes],
            "to": [p.end for p in net.pipes],
            "length_m": net.pipe_lengths,
            "diameter_in": [p.geometry.inlet_diameter / 0.0254 for p in net.pipes],
            "liquid_flow_kg_s": self.pipe_liquid_flow,
            "gas_flow_kg_s": self.pipe_gas_flow,
            "inlet_pressure_bar": self.pipe_inlet_pressure / BAR,
            "outlet_pressure_bar": self.pipe_outlet_pressure / BAR,
        }
        df = pd.DataFrame(data).set_index("pipe")
        df["pressure_drop_bar"] = df["inlet_pressure_bar"] - df["outlet_pressure_bar"]
        if self.profiles is not None:
            df["mean_liquid_holdup"] = self.profiles.liquid_holdup.mean(axis=1)
            idx = np.rint(np.median(self.profiles.pattern_index, axis=1)).astype(int)
            df["dominant_pattern"] = [PATTERN_NAMES.get(int(i), "mixed") for i in idx]
        return df

    def wells_table(self, wells=None) -> pd.DataFrame:
        wells = self.network.wells if wells is None else wells
        rows = []
        for w in wells:
            rows.append(
                {
                    "well": w.name,
                    "oil_rate_kg_s": self.liquid_inflow(w.downhole),
                    "gas_injection_kg_s": self.gas_inflow(w.gas_inlet),
                    "downhole_pressure_bar": self.pressure_bar(w.downhole),
                    "injection_point_pressure_bar": self.pressure_bar(w.injection_point),
                    "wellhead_pressure_bar": self.pressure_bar(w.wellhead),
                    "gas_inlet_pressure_bar": self.pressure_bar(w.gas_inlet),
                }
            )
        return pd.DataFrame(rows).set_index("well")

    def profile(self, pipe: str) -> pd.DataFrame:
        """Flow state along one pipe as a DataFrame (requires profiles)."""
        if self.profiles is None:
            raise ValueError("this solution was computed without pipe profiles")
        j = self.network.pipe_index(pipe)
        p = self.profiles
        return pd.DataFrame(
            {
                "s": p.positions,
                "axial_position_m": p.axial_position[j],
                "pressure_bar": p.pressure[j] / BAR,
                "liquid_holdup": p.liquid_holdup[j],
                "no_slip_holdup": p.no_slip_holdup[j],
                "liquid_velocity_m_s": p.liquid_velocity[j],
                "gas_velocity_m_s": p.gas_velocity[j],
                "superficial_liquid_velocity_m_s": p.superficial_liquid_velocity[j],
                "superficial_gas_velocity_m_s": p.superficial_gas_velocity[j],
                "gas_density_kg_m3": p.gas_density[j],
                "mixture_density_kg_m3": p.mixture_density[j],
                "pattern_index": p.pattern_index[j],
                "froude": p.froude[j],
                "pressure_gradient_pa_m": p.pressure_gradient[j],
            }
        )

    def __repr__(self) -> str:
        return f"NetworkSolution({self.network!r}, {self.status})"
