"""Fluid property model: incompressible liquid and polytropic gas (paper, Section 3.2)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FluidProperties:
    """Properties of the liquid and gas phases.

    The liquid is incompressible. The gas density follows the polytropic relation
    ``rho_G = rho_G0 * (P / P0) ** (1 / gamma)`` where ``gamma`` is the polytropic
    coefficient (``gamma = 1`` recovers isothermal ideal-gas behaviour).

    Default values are those of Table 1 of de Souza et al. (2010).
    """

    liquid_density: float = 850.0
    """Liquid density, kg/m³."""
    liquid_viscosity: float = 0.012
    """Liquid dynamic viscosity, Pa·s."""
    gas_viscosity: float = 8.0e-5
    """Gas dynamic viscosity, Pa·s."""
    gas_reference_density: float = 80.0
    """Gas density at ``reference_pressure``, kg/m³."""
    reference_pressure: float = 1.0e7
    """Reference pressure for the gas density, Pa (100 bar)."""
    gas_polytropic_coefficient: float = 1.25
    """Polytropic coefficient of the gas (dimensionless)."""
    interfacial_tension: float = 0.071
    """Gas-liquid interfacial tension, N/m."""
    gas_normal_density: float = 1.2
    """Gas density at normal conditions, kg/Nm³ (used to report volumetric gas rates)."""
    gas_molecular_weight: float = 25.0
    """Gas molecular weight, kg/kmol (used by the compressor model)."""

    def gas_density(self, pressure):
        """Gas density (kg/m³) at the given absolute pressure (Pa)."""
        pressure = np.asarray(pressure, dtype=float)
        return self.gas_reference_density * (pressure / self.reference_pressure) ** (
            1.0 / self.gas_polytropic_coefficient
        )

    def gas_density_derivative(self, pressure):
        """d(rho_G)/dP (kg/m³/Pa) at the given absolute pressure (Pa)."""
        pressure = np.asarray(pressure, dtype=float)
        return self.gas_density(pressure) / (self.gas_polytropic_coefficient * pressure)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}

    @classmethod
    def from_dict(cls, data: dict) -> FluidProperties:
        return cls(**data)
