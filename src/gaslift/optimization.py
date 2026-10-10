"""Gas-lift optimisation coupled to the network simulator (paper Section 5, Eq. 24).

Decision variables are the produced oil rate and injected gas rate of each well
and, for the design problem, the diameter of each well's gas supply lines. The
reservoir inflow performance ``Q_L = PI (P_res - P_downhole)`` enters as an
equality constraint, so that the simulator (which specifies rates) and the
reservoir model agree at the optimum. The problem is solved with SciPy's SLSQP,
a sequential quadratic programming method like the one used in the paper.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .cases import CaseSetup
from .economics import CompressorParameters, CostModel, EconomicEvaluation, evaluate_economics
from .results import NetworkSolution
from .simulation import SimulationOptions, Simulator
from .units import BAR, INCH

FAILED_OBJECTIVE = 1.0e6
"""Objective value returned when the simulation does not converge (penalty)."""


class Objective(StrEnum):
    MAX_PRODUCTION = "max_production"
    """Omega_1 = R (Eq. 25)."""
    MAX_PROFIT = "max_profit"
    """Omega_2 = R - AGC (Eq. 26)."""
    OPTIMAL_DESIGN = "optimal_design"
    """Omega_3 = R - AGC - ATCC (Eq. 27), with gas-line diameters as decision variables."""


@dataclass(frozen=True)
class Reservoir:
    """Reservoir pressure (bar) and productivity index ((kg/s)/bar) of one well."""

    pressure: float
    productivity_index: float = np.inf

    @property
    def inverse_productivity(self) -> float:
        return 0.0 if np.isinf(self.productivity_index) else 1.0 / self.productivity_index


@dataclass
class EvaluationRecord:
    x: np.ndarray
    converged: bool
    objective: float
    constraint: np.ndarray
    economics: EconomicEvaluation | None


@dataclass
class OptimizationResult:
    x: np.ndarray
    objective: float
    success: bool
    message: str
    n_evaluations: int
    n_iterations: int
    solution: NetworkSolution
    economics: EconomicEvaluation
    history: list[EvaluationRecord]
    problem: GasLiftProblem

    @property
    def oil_rates(self) -> np.ndarray:
        return self.problem.split(self.x)[0]

    @property
    def gas_rates(self) -> np.ndarray:
        return self.problem.split(self.x)[1]

    @property
    def diameters_in(self) -> np.ndarray | None:
        return self.problem.split(self.x)[2]

    def summary(self) -> pd.DataFrame:
        """Per-well summary of the optimum (rates, pressures, diameters)."""
        q_oil, q_gas, d = self.problem.split(self.x)
        rows = {
            "productivity_index": [r.productivity_index for r in self.problem.reservoirs],
            "reservoir_pressure_bar": [r.pressure for r in self.problem.reservoirs],
            "downhole_pressure_bar": self.economics.downhole_pressure,
            "oil_rate_kg_s": q_oil,
            "gas_injection_kg_s": q_gas,
        }
        if d is not None:
            rows["gas_line_diameter_in"] = d
        return pd.DataFrame(rows, index=[w.name for w in self.problem.case.wells])

    def history_frame(self) -> pd.DataFrame:
        recs = []
        for k, r in enumerate(self.history):
            q_oil, q_gas, d = self.problem.split(r.x)
            rec = {
                "evaluation": k,
                "converged": r.converged,
                "objective": r.objective,
                "max_abs_constraint": float(np.max(np.abs(r.constraint)))
                if r.constraint.size
                else 0.0,
            }
            for i, (qo, qg) in enumerate(zip(q_oil, q_gas, strict=True)):
                rec[f"q_oil_{i + 1}"] = qo
                rec[f"q_gas_{i + 1}"] = qg
            if d is not None:
                for i, di in enumerate(d):
                    rec[f"d_gas_{i + 1}_in"] = di
            if r.economics is not None:
                rec["profit"] = r.economics.profit
                rec["compressor_power_kW"] = r.economics.compressor_power
            recs.append(rec)
        return pd.DataFrame(recs).set_index("evaluation")


@dataclass
class GasLiftProblem:
    """Optimisation problem definition (Eq. 24)."""

    case: CaseSetup
    reservoirs: Sequence[Reservoir]
    objective: Objective = Objective.MAX_PRODUCTION
    cost: CostModel = field(default_factory=CostModel)
    compressor: CompressorParameters = field(default_factory=CompressorParameters)
    oil_bounds: tuple[float, float] = (0.0, 20.0)
    """Bounds of the produced oil rate per well, kg/s."""
    gas_bounds: tuple[float, float] = (0.0, 4.0)
    """Bounds of the injected gas rate per well, kg/s."""
    diameter_bounds: tuple[float, float] = (1.0, 4.0)
    """Bounds of the gas-line diameters, inches (design problem only)."""
    simulation: SimulationOptions = field(
        default_factory=lambda: SimulationOptions(compute_profiles=False)
    )
    history: list[EvaluationRecord] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.reservoirs) != self.case.n_wells:
            raise ValueError("one Reservoir per well is required")
        self.cost = self._cost_for_objective(self.cost)
        self._simulator = Simulator(self.case.network, self.case.specifications, self.simulation)
        self._cache: dict[bytes, tuple[NetworkSolution, EconomicEvaluation | None]] = {}

    def _cost_for_objective(self, cost: CostModel) -> CostModel:
        from dataclasses import replace

        if self.objective is Objective.MAX_PRODUCTION:
            return replace(
                cost,
                include_fuel_gas_cost=False,
                include_capital_cost=False,
                include_injected_gas_cost=False,
            )
        if self.objective is Objective.MAX_PROFIT:
            return replace(cost, include_fuel_gas_cost=True, include_capital_cost=False)
        return replace(cost, include_fuel_gas_cost=True, include_capital_cost=True)

    # ------------------------------------------------------------------ decision vector
    @property
    def n_wells(self) -> int:
        return self.case.n_wells

    @property
    def has_design_variables(self) -> bool:
        return self.objective is Objective.OPTIMAL_DESIGN

    @property
    def n_variables(self) -> int:
        return 3 * self.n_wells if self.has_design_variables else 2 * self.n_wells

    def split(self, x) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        x = np.asarray(x, dtype=float)
        n = self.n_wells
        d = x[2 * n : 3 * n] if self.has_design_variables else None
        return x[:n], x[n : 2 * n], d

    def initial_vector(self, q_oil, q_gas, diameters_in=None) -> np.ndarray:
        n = self.n_wells
        parts = [
            np.broadcast_to(np.asarray(q_oil, float), (n,)),
            np.broadcast_to(np.asarray(q_gas, float), (n,)),
        ]
        if self.has_design_variables:
            if diameters_in is None:
                diameters_in = [
                    self.case.network.pipe(w.gas_line_1).geometry.inlet_diameter / INCH
                    for w in self.case.wells
                ]
            parts.append(np.broadcast_to(np.asarray(diameters_in, float), (n,)))
        return np.concatenate(parts)

    def bounds(self) -> list[tuple[float, float]]:
        n = self.n_wells
        b = [self.oil_bounds] * n + [self.gas_bounds] * n
        if self.has_design_variables:
            b += [self.diameter_bounds] * n
        return b

    # ------------------------------------------------------------------ simulation
    def simulate(self, x) -> tuple[NetworkSolution, EconomicEvaluation | None]:
        """Simulate the network at decision vector ``x`` (cached by value)."""
        x = np.asarray(x, dtype=float)
        key = x.tobytes()
        if key in self._cache:
            return self._cache[key]
        q_oil, q_gas, d = self.split(x)
        self.case.set_well_rates(q_oil, q_gas)
        if d is not None:
            for w, di in zip(self.case.wells, d, strict=True):
                self.case.set_gas_line_diameter(w, float(di) * INCH)
        sol = self._simulator.solve()
        econ = None
        if sol.converged:
            econ = evaluate_economics(
                sol,
                self.case.wells,
                self.cost,
                self.compressor,
                pipeline_diameters_in=d if d is not None else None,
            )
        if len(self._cache) > 64:
            self._cache.clear()
        self._cache[key] = (sol, econ)
        return sol, econ

    def objective_value(self, x) -> float:
        sol, econ = self.simulate(x)
        value = FAILED_OBJECTIVE if econ is None else econ.objective
        self.history.append(
            EvaluationRecord(
                np.array(x, float), sol.converged, value, self.constraint_value(x), econ
            )
        )
        return value

    def constraint_value(self, x) -> np.ndarray:
        """Reservoir inflow residuals ``P_res - P_downhole - Q_L / PI`` per well, bar."""
        sol, _ = self.simulate(x)
        q_oil, _, _ = self.split(x)
        if not sol.converged:
            return np.zeros(self.n_wells)
        p_down = sol.well_downhole_pressures(self.case.wells) / BAR
        return np.array(
            [
                r.pressure - pd - qo * r.inverse_productivity
                for r, pd, qo in zip(self.reservoirs, p_down, q_oil, strict=True)
            ]
        )

    # ------------------------------------------------------------------ solve
    def solve(
        self,
        x0,
        *,
        max_iterations: int = 100,
        ftol: float = 1e-6,
        finite_difference_step: float = 1e-5,
        verbose: bool = False,
    ) -> OptimizationResult:
        """Run SLSQP from ``x0`` and return the best evaluated point."""
        x0 = np.asarray(x0, dtype=float)
        if x0.size != self.n_variables:
            raise ValueError(f"x0 must have {self.n_variables} entries")
        self.history.clear()
        self._simulator.reset()
        n_iter = [0]

        def cb(_x):
            n_iter[0] += 1
            if verbose:
                best = min(
                    (r for r in self.history if r.converged),
                    key=lambda r: r.objective,
                    default=None,
                )
                value = best.objective if best else float("nan")
                print(
                    f"iteration {n_iter[0]}: evaluations={len(self.history)} "
                    f"best objective={value:.4f}"
                )

        res = minimize(
            self.objective_value,
            x0,
            method="SLSQP",
            bounds=self.bounds(),
            constraints=[{"type": "eq", "fun": self.constraint_value}],
            options={
                "maxiter": max_iterations,
                "ftol": ftol,
                "eps": finite_difference_step,
                "disp": False,
            },
            callback=cb,
        )
        # Best feasible evaluation (SLSQP's last iterate is not always the best one).
        feasible = [r for r in self.history if r.converged and np.max(np.abs(r.constraint)) < 1e-3]
        best = min(feasible, key=lambda r: r.objective) if feasible else None
        x_best = best.x if best is not None else np.asarray(res.x, float)
        sol, econ = self.simulate(x_best)
        if econ is None:
            raise RuntimeError(
                "optimisation finished at a point where the simulation does not converge"
            )
        return OptimizationResult(
            x=x_best,
            objective=econ.objective,
            success=bool(res.success),
            message=str(res.message),
            n_evaluations=len(self.history),
            n_iterations=n_iter[0],
            solution=sol,
            economics=econ,
            history=list(self.history),
            problem=self,
        )
