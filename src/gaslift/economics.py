"""Economic evaluation of a gas-lift system (paper Section 5: Eqs. 25-35).

Revenue from produced oil, fuel-gas cost of the compressor driver, and the
annualised capital cost of compressor, gas turbine and gas pipelines. Defaults
reproduce the numbers of the paper; the ``*_formulation`` switches expose the
textbook forms where the original code deviated from them (see
``docs/reproduction.md`` and ``docs/errata.md``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .fluids import FluidProperties
from .network import Well
from .results import NetworkSolution
from .units import BAR, BARREL, INCH, KW_TO_HP, R_GAS, SECONDS_PER_YEAR


@dataclass(frozen=True)
class CompressorParameters:
    """Multi-stage compressor with complete intercooling (Eq. 31, Table 16)."""

    stages: int = 3
    polytropic_exponent: float = 1.4
    """Exponent ``gamma`` of Eq. 31. The paper's Table 16 lists 1.25 but the published
    numbers were produced with 1.4."""
    suction_temperature: float = 300.0
    """K"""
    suction_pressure: float = 5.0
    """bar"""
    efficiency: float = 0.8
    molecular_weight: float = 25.0
    """kg/kmol"""


@dataclass(frozen=True)
class CostModel:
    """Prices, cost correlations and which cost terms enter the objective function."""

    oil_price: float = 51.15
    """USD per barrel."""
    gas_price: float = 5.34 * 0.035315
    """USD per Nm³ (5.34 USD/kcf at the US wellhead in 2007)."""
    capital_charge_factor: float = 1.0 / 3.0
    """CCF of Eq. 29, 1/yr (15 % interest, 12 years)."""
    osbl_factor: float = 3.0
    """OSBL = osbl_factor * ISBL (Douglas, 1988)."""
    fuel_gas_per_kw: float = 8.45e-3
    """Fuel gas consumption of the turbine driver, 1e6 Nm³/yr per kW (Eq. 34)."""
    include_fuel_gas_cost: bool = False
    include_capital_cost: bool = False
    include_injected_gas_cost: bool = False
    """Charge the injected gas itself at ``gas_price`` (never used in the paper)."""
    compressor_cost_power_unit: str = "hp"
    """Unit of the power fed to the Couper et al. correlations: ``"hp"`` (as the original
    code) or ``"kW"`` (as printed in Eqs. 30 and 32)."""
    pipeline_cost_formulation: str = "matlab"
    """``"matlab"``: diameter multiplied by the exponent (original code); ``"mccoy"``:
    diameter raised to the exponent (McCoy, 2008; Eq. 33)."""
    discharge_pressure_node: str = "gas_junction"
    """Which well node is taken as compressor discharge: ``"gas_junction"`` (node n+4,
    used by the original code) or ``"gas_inlet"`` (node n+5, the platform)."""

    @property
    def price_ratio(self) -> float:
        """Fuel gas price over oil price, (USD/Nm³)/(USD/barrel)."""
        return self.gas_price / self.oil_price


def compressor_power(
    discharge_pressure_bar: float, gas_mass_rate: float, comp: CompressorParameters
) -> float:
    """Compression power, kW (Eq. 31)."""
    n, g = comp.stages, comp.polytropic_exponent
    ratio = max(discharge_pressure_bar, comp.suction_pressure) / comp.suction_pressure
    work = (
        n
        * g
        * R_GAS
        * comp.suction_temperature
        / (comp.efficiency * (g - 1.0))
        * (ratio ** ((g - 1.0) / (n * g)) - 1.0)
    )  # kJ/kmol
    return work * gas_mass_rate / comp.molecular_weight


def _cost_power(power_kw: float, unit: str) -> float:
    if unit == "hp":
        return power_kw * KW_TO_HP
    if unit == "kW":
        return power_kw
    raise ValueError(f"unknown power unit {unit!r}")


def compressor_capital_cost(power_kw: float, unit: str = "hp") -> float:
    """Installed cost of a centrifugal compressor, 1e6 USD (Couper et al., 2005; Eq. 30)."""
    return 14.3e-3 * _cost_power(power_kw, unit) ** 0.62


def turbine_capital_cost(power_kw: float, unit: str = "hp") -> float:
    """Installed cost of a gas turbine driver, 1e6 USD (Couper et al., 2005; Eq. 32)."""
    return 7.19e-3 * _cost_power(power_kw, unit) ** 0.81


_MCCOY = (
    # (a0 + a4, length exponent a6, diameter exponent a7) for materials, labour, misc.
    (3.112, 0.901, 1.590),
    (4.487 - 0.216, 0.820, 0.940),
    (4.390, 0.783, 0.791),
)


def pipeline_capital_cost(diameters_in, lengths_m, formulation: str = "matlab") -> float:
    """Installed cost of gas pipelines, 1e6 USD (McCoy, 2008, US South-West, no right of way).

    ``diameters_in`` in inches and ``lengths_m`` in metres, one entry per pipeline.
    """
    d = np.atleast_1d(np.asarray(diameters_in, dtype=float))
    length_km = np.atleast_1d(np.asarray(lengths_m, dtype=float)) / 1000.0
    total = 0.0
    for a04, a6, a7 in _MCCOY:
        if formulation == "matlab":
            d_term = d * a7
        elif formulation == "mccoy":
            d_term = d**a7
        else:
            raise ValueError(f"unknown pipeline cost formulation {formulation!r}")
        total += float(np.sum(10.0**a04 * length_km**a6 * d_term))
    return total / 1.0e6


def fuel_gas_volume(power_kw: float, cost: CostModel) -> float:
    """Fuel gas burnt by the turbine driver, 1e6 Nm³/yr."""
    return cost.fuel_gas_per_kw * power_kw


def oil_volume_rate(mass_rate, fluid: FluidProperties) -> np.ndarray:
    """Oil production, 1e6 barrel/yr, from mass rate(s) in kg/s."""
    return (
        np.asarray(mass_rate, dtype=float) / fluid.liquid_density / BARREL * SECONDS_PER_YEAR / 1e6
    )


def gas_volume_rate(mass_rate, fluid: FluidProperties) -> np.ndarray:
    """Gas rate, 1e6 Nm³/yr, from mass rate(s) in kg/s."""
    return np.asarray(mass_rate, dtype=float) / fluid.gas_normal_density * SECONDS_PER_YEAR / 1e6


@dataclass
class EconomicEvaluation:
    """All economic quantities of one operating point (money in 1e6 USD or 1e6 USD/yr)."""

    oil_rate: np.ndarray
    """Produced oil per well, kg/s."""
    gas_rate: np.ndarray
    """Injected gas per well, kg/s."""
    downhole_pressure: np.ndarray
    """bar"""
    discharge_pressure: float
    """Compressor discharge pressure, bar."""
    compressor_power: float
    """kW"""
    oil_volume: np.ndarray
    """1e6 barrel/yr per well."""
    injected_gas_volume: np.ndarray
    """1e6 Nm³/yr per well."""
    fuel_gas_volume: float
    """1e6 Nm³/yr."""
    revenue: float
    fuel_gas_cost: float
    injected_gas_cost: float
    compressor_cost: float
    turbine_cost: float
    pipeline_cost: float
    isbl: float
    osbl: float
    annualized_capital_cost: float
    profit: float
    pipeline_diameters_in: np.ndarray = field(default_factory=lambda: np.zeros(0))
    pipeline_lengths_m: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def objective(self) -> float:
        """Value to minimise (``-profit``)."""
        return -self.profit

    @property
    def operating_cost(self) -> float:
        return self.fuel_gas_cost + self.injected_gas_cost

    def to_series(self) -> pd.Series:
        data = {
            "total_oil_rate_kg_s": float(np.sum(self.oil_rate)),
            "total_gas_injection_kg_s": float(np.sum(self.gas_rate)),
            "total_oil_volume_1e6_bbl_yr": float(np.sum(self.oil_volume)),
            "total_injected_gas_1e6_m3_yr": float(np.sum(self.injected_gas_volume)),
            "fuel_gas_1e6_m3_yr": self.fuel_gas_volume,
            "discharge_pressure_bar": self.discharge_pressure,
            "compressor_power_kW": self.compressor_power,
            "revenue_1e6_usd_yr": self.revenue,
            "fuel_gas_cost_1e6_usd_yr": self.fuel_gas_cost,
            "compressor_cost_1e6_usd": self.compressor_cost,
            "turbine_cost_1e6_usd": self.turbine_cost,
            "pipeline_cost_1e6_usd": self.pipeline_cost,
            "ISBL_1e6_usd": self.isbl,
            "OSBL_1e6_usd": self.osbl,
            "annualized_capital_cost_1e6_usd_yr": self.annualized_capital_cost,
            "profit_1e6_usd_yr": self.profit,
        }
        return pd.Series(data)


def evaluate_economics(
    solution: NetworkSolution,
    wells: list[Well] | tuple[Well, ...],
    cost: CostModel,
    compressor: CompressorParameters | None = None,
    pipeline_diameters_in=None,
    pipeline_lengths_m=None,
) -> EconomicEvaluation:
    """Evaluate revenue, costs and profit for a converged network solution.

    Pipeline capital cost uses one gas pipeline per well with the given diameters
    (inches) and lengths (m); by default they are read from the network geometry.
    """
    compressor = compressor or CompressorParameters()
    fluid = solution.network.fluid
    oil = solution.well_oil_rates(wells)
    gas = solution.well_gas_rates(wells)
    p_down = solution.well_downhole_pressures(wells) / BAR
    if cost.discharge_pressure_node == "gas_junction":
        p_dis = max(solution.pressure(w.gas_junction) for w in wells) / BAR
    elif cost.discharge_pressure_node == "gas_inlet":
        p_dis = max(solution.pressure(w.gas_inlet) for w in wells) / BAR
    else:
        raise ValueError(f"unknown discharge pressure node {cost.discharge_pressure_node!r}")

    power = compressor_power(p_dis, float(np.sum(gas)), compressor)
    oil_vol = oil_volume_rate(oil, fluid)
    gas_vol = gas_volume_rate(gas, fluid)
    fuel = fuel_gas_volume(power, cost) if cost.include_fuel_gas_cost else 0.0

    revenue = float(np.sum(oil_vol) * cost.oil_price)
    fuel_cost = fuel * cost.gas_price
    inj_cost = float(np.sum(gas_vol) * cost.gas_price) if cost.include_injected_gas_cost else 0.0

    if cost.include_capital_cost:
        if pipeline_diameters_in is None:
            pipeline_diameters_in = [
                solution.network.pipe(w.gas_line_1).geometry.inlet_diameter / INCH for w in wells
            ]
        if pipeline_lengths_m is None:
            pipeline_lengths_m = [
                sum(solution.network.pipe(p).length for p in w.gas_lines) for w in wells
            ]
        c_comp = compressor_capital_cost(power, cost.compressor_cost_power_unit)
        c_turb = turbine_capital_cost(power, cost.compressor_cost_power_unit)
        c_pipe = pipeline_capital_cost(
            pipeline_diameters_in, pipeline_lengths_m, cost.pipeline_cost_formulation
        )
        isbl = c_comp + c_turb + c_pipe
        osbl = cost.osbl_factor * isbl
        atcc = cost.capital_charge_factor * osbl
        d_in, l_m = np.asarray(pipeline_diameters_in, float), np.asarray(pipeline_lengths_m, float)
    else:
        c_comp = c_turb = c_pipe = isbl = osbl = atcc = 0.0
        d_in, l_m = np.zeros(0), np.zeros(0)

    profit = revenue - fuel_cost - inj_cost - atcc
    return EconomicEvaluation(
        oil_rate=oil,
        gas_rate=gas,
        downhole_pressure=p_down,
        discharge_pressure=p_dis,
        compressor_power=power,
        oil_volume=oil_vol,
        injected_gas_volume=gas_vol,
        fuel_gas_volume=fuel,
        revenue=revenue,
        fuel_gas_cost=fuel_cost,
        injected_gas_cost=inj_cost,
        compressor_cost=c_comp,
        turbine_cost=c_turb,
        pipeline_cost=c_pipe,
        isbl=isbl,
        osbl=osbl,
        annualized_capital_cost=atcc,
        profit=profit,
        pipeline_diameters_in=d_in,
        pipeline_lengths_m=l_m,
    )
