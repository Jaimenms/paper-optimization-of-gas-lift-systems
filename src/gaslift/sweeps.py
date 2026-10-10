"""Parametric sweeps of injected gas and produced oil (paper Figs. 6-8, 11-12)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from .cases import CaseSetup
from .results import NetworkSolution
from .simulation import SimulationOptions, Simulator
from .units import BAR


@dataclass
class RateSweep:
    """Pressures over a grid of injected gas (columns) and produced oil (rows) rates."""

    q_gas: np.ndarray
    """Injected gas rates, kg/s (per well, same for all wells)."""
    q_oil: np.ndarray
    """Produced oil rates, kg/s (per well, same for all wells)."""
    downhole_pressure: np.ndarray
    """bar, shape ``(len(q_oil), len(q_gas), n_wells)``; NaN where not converged."""
    gas_junction_pressure: np.ndarray
    """bar, pressure at the gas-lines interconnection (node n+4) of each well."""
    gas_inlet_pressure: np.ndarray
    """bar, pressure at the platform gas inlet (node n+5) of each well."""
    converged: np.ndarray
    iterations: np.ndarray

    @property
    def grid(self) -> tuple[np.ndarray, np.ndarray]:
        """``(Q_GAS, Q_OIL)`` meshgrid arrays matching the pressure arrays."""
        return np.meshgrid(self.q_gas, self.q_oil)

    def downhole(self, well: int = 0) -> np.ndarray:
        return self.downhole_pressure[:, :, well]

    def gas_junction(self, well: int = 0) -> np.ndarray:
        return self.gas_junction_pressure[:, :, well]

    def gas_inlet(self, well: int = 0) -> np.ndarray:
        return self.gas_inlet_pressure[:, :, well]


def sweep_rates(
    case: CaseSetup,
    q_gas: np.ndarray,
    q_oil: np.ndarray,
    options: SimulationOptions | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> RateSweep:
    """Simulate the case for every combination of gas and oil rates.

    The grid is traversed column by column in a serpentine order, warm-starting
    each point from its predecessor; failed points fall back to the linear initial
    estimate (same strategy as the original scripts).
    """
    q_gas = np.asarray(q_gas, dtype=float)
    q_oil = np.asarray(q_oil, dtype=float)
    nx, ny, nw = len(q_oil), len(q_gas), case.n_wells
    opts = options or SimulationOptions(compute_profiles=False)
    sim = Simulator(case.network, case.specifications, opts)
    p_down = np.full((nx, ny, nw), np.nan)
    p_junc = np.full((nx, ny, nw), np.nan)
    p_inlet = np.full((nx, ny, nw), np.nan)
    conv = np.zeros((nx, ny), dtype=bool)
    iters = np.zeros((nx, ny), dtype=int)
    total, count = nx * ny, 0
    for jj in range(ny):
        rows = range(nx - 1, -1, -1) if jj % 2 == 1 else range(nx)
        for ii in rows:
            case.set_well_rates(q_oil[ii], q_gas[jj])
            sol: NetworkSolution = sim.solve()
            if not sol.converged:
                sol = sim.solve(warm_start=False)
            conv[ii, jj] = sol.converged
            iters[ii, jj] = sol.status.iterations
            if sol.converged:
                p_down[ii, jj] = sol.well_downhole_pressures(case.wells) / BAR
                p_junc[ii, jj] = np.array([sol.pressure(w.gas_junction) for w in case.wells]) / BAR
                p_inlet[ii, jj] = sol.well_gas_inlet_pressures(case.wells) / BAR
            count += 1
            if progress is not None:
                progress(count, total)
    return RateSweep(q_gas, q_oil, p_down, p_junc, p_inlet, conv, iters)
