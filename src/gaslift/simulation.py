"""High-level simulation facade: initial estimate, Newton solution and fall-backs."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from .linear_model import linear_initial_estimate
from .model import ModelOptions, NetworkModel
from .network import Network, Specifications
from .pipe_model import profile_states
from .results import NetworkSolution, PipeProfiles
from .solver import NewtonOptions, SolverStatus, newton_solve


@dataclass(frozen=True)
class SimulationOptions:
    model: ModelOptions = field(default_factory=ModelOptions)
    newton: NewtonOptions = field(default_factory=NewtonOptions)
    fallback_nonrigorous: bool = True
    """On failure, retry with the non-rigorous model and then refine rigorously."""
    compute_profiles: bool = True


class Simulator:
    """Solves the stationary two-phase network model for given specifications.

    Typical use::

        sim = Simulator(network, specifications)
        solution = sim.solve()
        solution.wells_table()

    The simulator keeps the last converged state and uses it as the initial estimate
    of the next call, which is what makes parametric sweeps and optimisations fast.
    """

    def __init__(
        self,
        network: Network,
        specifications: Specifications,
        options: SimulationOptions | None = None,
    ) -> None:
        self.network = network
        self.specifications = specifications
        self.options = options or SimulationOptions()
        self.last_state: np.ndarray | None = None

    # ------------------------------------------------------------------ helpers
    def model(
        self, specifications: Specifications | None = None, rigorous: bool | None = None
    ) -> NetworkModel:
        specs = specifications or self.specifications
        opts = (
            self.options.model
            if rigorous is None
            else replace(self.options.model, rigorous=rigorous)
        )
        return NetworkModel(self.network, specs, opts)

    def initial_estimate(self, specifications: Specifications | None = None) -> np.ndarray:
        """Initial estimate from the linear laminar model."""
        return linear_initial_estimate(self.model(specifications))

    def reset(self) -> None:
        self.last_state = None

    # ------------------------------------------------------------------ solve
    def solve(
        self,
        specifications: Specifications | None = None,
        x0: np.ndarray | None = None,
        warm_start: bool = True,
    ) -> NetworkSolution:
        """Solve the network for ``specifications`` (defaults to the simulator's).

        Strategy (same as the original code): Newton on the rigorous model from the
        initial estimate; if that fails, Newton on the non-rigorous model followed by
        a rigorous refinement; if everything fails the linear estimate is retried as
        starting point. The initial estimate is ``x0`` if given, otherwise the last
        converged state (when ``warm_start``) or the linear model.
        """
        specs = specifications or self.specifications
        model = self.model(specs)
        starts: list[np.ndarray] = []
        if x0 is not None:
            starts.append(np.asarray(x0, dtype=float))
        elif warm_start and self.last_state is not None:
            starts.append(self.last_state)
        linear = None
        if not starts:
            linear = linear_initial_estimate(model)
            starts.append(linear)

        x, status = None, None
        for start in starts:
            x, status = self._solve_from(model, start)
            if status.converged:
                break
        if (x is None or not status.converged) and linear is None:
            linear = linear_initial_estimate(model)
            x, status = self._solve_from(model, linear)

        profiles = None
        if status.converged:
            self.last_state = x.copy()
            if self.options.compute_profiles:
                profiles = self._profiles(model, x)
        return NetworkSolution(self.network, specs, x, status, model.pressure_scale, profiles)

    def _solve_from(
        self, model: NetworkModel, start: np.ndarray
    ) -> tuple[np.ndarray, SolverStatus]:
        rigorous = model.options.rigorous
        x, status = newton_solve(model, start, self.options.newton, rigorous=rigorous)
        if status.converged or not rigorous or not self.options.fallback_nonrigorous:
            return x, status
        x1, st1 = newton_solve(model, start, self.options.newton, rigorous=False)
        if not st1.converged:
            return x, status
        x2, st2 = newton_solve(model, x1, self.options.newton, rigorous=True)
        return (x2, st2) if st2.converged else (x, status)

    def _profiles(self, model: NetworkModel, x: np.ndarray) -> PipeProfiles:
        v = model.split(x)
        _, integration = model.evaluate(x)
        n_pts = model.options.profile_points
        s, state = profile_states(
            model.batch,
            model.fluid,
            integration,
            v.PI * model.pressure_scale,
            v.qL,
            v.qG,
            n_pts,
            model.options.min_pressure,
        )
        n_pipes = model.layout.n_pipes
        shape = (n_pipes, n_pts)
        lengths = model.batch.lengths[:, None]
        return PipeProfiles(
            positions=s,
            axial_position=s[None, :] * lengths,
            pressure=state.pressure.reshape(shape),
            liquid_holdup=state.liquid_holdup.reshape(shape),
            gas_holdup=state.gas_holdup.reshape(shape),
            no_slip_holdup=state.no_slip_holdup.reshape(shape),
            liquid_velocity=state.liquid_velocity.reshape(shape),
            gas_velocity=state.gas_velocity.reshape(shape),
            superficial_liquid_velocity=state.superficial_liquid_velocity.reshape(shape),
            superficial_gas_velocity=state.superficial_gas_velocity.reshape(shape),
            gas_density=state.gas_density.reshape(shape),
            mixture_density=state.mixture_density.reshape(shape),
            pattern_index=state.pattern_index.reshape(shape),
            angle=state.angle.reshape(shape),
            diameter=state.diameter.reshape(shape),
            froude=state.froude.reshape(shape),
            pressure_gradient=state.dpds.reshape(shape) / lengths,
        )


def simulate(
    network: Network,
    specifications: Specifications,
    options: SimulationOptions | None = None,
    x0: np.ndarray | None = None,
) -> NetworkSolution:
    """One-shot convenience wrapper around :class:`Simulator`."""
    return Simulator(network, specifications, options).solve(x0=x0)
