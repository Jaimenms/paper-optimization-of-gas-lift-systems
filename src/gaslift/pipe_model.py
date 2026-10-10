"""Pressure gradient along a pipe (paper Eqs. 1-4) and its integration (Eqs. 8-10).

The momentum balance of the mixture is integrated along the normalised axial
coordinate ``s`` in ``[0, 1]`` of every pipe, in both directions (from the inlet
pressure forwards and from the outlet pressure backwards). All pipes of a
network, and several perturbed copies of them, are integrated together by the
batched Runge-Kutta solver of :mod:`gaslift.integrator`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .beggs_brill import closure_relations
from .fluids import FluidProperties
from .integrator import BatchSolution, integrate_batch
from .network import Network
from .units import G


@dataclass(frozen=True)
class PipeBatch:
    """Geometry of all pipes of a network in padded arrays for vectorised evaluation."""

    lengths: np.ndarray
    equivalent_length_factor: np.ndarray
    knot_s: np.ndarray
    knot_angle: np.ndarray
    knot_diameter: np.ndarray
    knot_roughness: np.ndarray
    inlet_area: np.ndarray
    outlet_area: np.ndarray

    @classmethod
    def from_network(cls, network: Network) -> PipeBatch:
        geoms = [p.geometry for p in network.pipes]
        n = len(geoms)
        k_max = max(len(g.knots[0]) for g in geoms)
        knot_s = np.zeros((n, k_max))
        knot_angle = np.zeros((n, k_max))
        knot_d = np.zeros((n, k_max))
        knot_r = np.zeros((n, k_max))
        for i, g in enumerate(geoms):
            s, a, d, r = g.knots
            k = len(s)
            knot_s[i, :k], knot_angle[i, :k], knot_d[i, :k], knot_r[i, :k] = s, a, d, r
            knot_s[i, k:], knot_angle[i, k:], knot_d[i, k:], knot_r[i, k:] = (
                s[-1],
                a[-1],
                d[-1],
                r[-1],
            )
        lengths = np.array([g.length for g in geoms])
        leq = np.array([g.equivalent_length for g in geoms])
        return cls(
            lengths=lengths,
            equivalent_length_factor=1.0 + leq / lengths,
            knot_s=knot_s,
            knot_angle=knot_angle,
            knot_diameter=knot_d,
            knot_roughness=knot_r,
            inlet_area=np.pi * np.array([g.inlet_diameter for g in geoms]) ** 2 / 4.0,
            outlet_area=np.pi * np.array([g.outlet_diameter for g in geoms]) ** 2 / 4.0,
        )

    @property
    def n_pipes(self) -> int:
        return len(self.lengths)

    def properties_at(self, pipe_idx: np.ndarray, s: np.ndarray):
        """Inclination, diameter and roughness of ``pipe_idx`` at axial positions ``s`` (m)."""
        xs = self.knot_s[pipe_idx]
        s = np.clip(s, xs[:, 0], xs[:, -1])
        k = np.sum(xs[:, 1:] <= s[:, None], axis=1)
        k = np.clip(k, 0, xs.shape[1] - 2)
        rows = np.arange(len(pipe_idx))
        x0, x1 = xs[rows, k], xs[rows, k + 1]
        width = x1 - x0
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(width > 0, (s - x0) / width, 0.0)
        out = []
        for table in (self.knot_angle, self.knot_diameter, self.knot_roughness):
            t = table[pipe_idx]
            f0, f1 = t[rows, k], t[rows, k + 1]
            out.append(f0 + w * (f1 - f0))
        return tuple(out)


@dataclass(frozen=True)
class GradientState:
    """Local two-phase flow state along pipes (arrays over the evaluated points)."""

    pressure: np.ndarray
    dpds: np.ndarray
    """Pressure derivative with respect to the normalised axial coordinate, Pa."""
    angle: np.ndarray
    diameter: np.ndarray
    area: np.ndarray
    liquid_holdup: np.ndarray
    gas_holdup: np.ndarray
    no_slip_holdup: np.ndarray
    liquid_velocity: np.ndarray
    gas_velocity: np.ndarray
    superficial_liquid_velocity: np.ndarray
    superficial_gas_velocity: np.ndarray
    gas_density: np.ndarray
    mixture_density: np.ndarray
    shear_stress: np.ndarray
    froude: np.ndarray
    pattern_index: np.ndarray
    two_phase_friction_factor: np.ndarray


def pressure_gradient(
    batch: PipeBatch,
    fluid: FluidProperties,
    pipe_idx: np.ndarray,
    position: np.ndarray,
    dp: np.ndarray,
    p_edge: np.ndarray,
    q_l: np.ndarray,
    q_g: np.ndarray,
    direction: np.ndarray,
    min_pressure: float = 1.0e5,
) -> GradientState:
    """Evaluate ``dP/ds`` (Eq. 2 times the pipe length) at points along pipes.

    Parameters
    ----------
    pipe_idx : int array
        Pipe of each evaluated point.
    position : array
        Integration coordinate ``t`` in ``[0, 1]`` (forward: measured from the inlet;
        backward: measured from the outlet).
    dp : array
        Accumulated pressure change since the start of the integration, Pa.
    p_edge : array
        Pressure at the start of the integration (inlet for forward, outlet for backward), Pa.
    q_l, q_g : array
        Liquid and gas mass flow rates along the digraph orientation, kg/s.
    direction : array
        ``+1`` for forward integration, ``-1`` for backward integration.
    """
    pipe_idx = np.asarray(pipe_idx)
    position = np.asarray(position, dtype=float)
    direction = np.asarray(direction, dtype=float)
    s_norm = np.where(direction > 0, position, 1.0 - position)
    pressure = np.maximum(p_edge + direction * dp, min_pressure)

    length = batch.lengths[pipe_idx]
    angle, diameter, roughness = batch.properties_at(pipe_idx, s_norm * length)
    factor = batch.equivalent_length_factor[pipe_idx]
    area = np.pi * diameter**2 / 4.0

    rho_l = fluid.liquid_density
    rho_g = fluid.gas_density(pressure)
    drho_dp = rho_g / (fluid.gas_polytropic_coefficient * pressure)
    v_sl = q_l / (area * rho_l)
    v_sg = q_g / (area * rho_g)

    cl = closure_relations(pressure, v_sl, v_sg, diameter, roughness, angle, fluid)
    v_sg_p = v_sg * (1.0 + 1e-6) + 1e-9
    cl_p = closure_relations(pressure, v_sl, v_sg_p, diameter, roughness, angle, fluid)
    dhl_dvsg = (cl_p.liquid_holdup - cl.liquid_holdup) / (v_sg_p - v_sg)

    comp = drho_dp * v_sg / rho_g
    a_l = comp * (np.abs(cl.liquid_velocity) * dhl_dvsg / (cl.liquid_holdup + 1e-6))
    a_g = -comp * (np.abs(cl.gas_velocity) * dhl_dvsg / (cl.gas_holdup + 1e-6) + 1.0)

    dpdl = -(
        area * cl.mixture_density * G * np.sin(angle) + np.pi * diameter * cl.shear_stress * factor
    ) / (area + a_l * q_l + a_g * q_g)
    dpds = length * dpdl
    dpds = np.where((pressure <= min_pressure) & (dpds < 0.0), 0.0, dpds)

    return GradientState(
        pressure=pressure,
        dpds=dpds,
        angle=angle,
        diameter=diameter,
        area=area,
        liquid_holdup=cl.liquid_holdup,
        gas_holdup=cl.gas_holdup,
        no_slip_holdup=cl.no_slip_holdup,
        liquid_velocity=cl.liquid_velocity,
        gas_velocity=cl.gas_velocity,
        superficial_liquid_velocity=v_sl,
        superficial_gas_velocity=v_sg,
        gas_density=rho_g,
        mixture_density=cl.mixture_density,
        shear_stress=cl.shear_stress,
        froude=cl.froude,
        pattern_index=cl.pattern.index,
        two_phase_friction_factor=cl.two_phase_friction_factor,
    )


@dataclass(frozen=True)
class IntegrationOptions:
    """Tolerances of the pressure-drop integration (MATLAB ``odeset`` values of the paper)."""

    rtol: float = 1e-2
    atol: float = 1e-6
    first_step: float = 1e-2
    max_step: float = 0.5
    min_pressure: float = 1.0e5
    max_steps: int = 200
    """Problems needing more accepted steps are flagged as failed (NaN pressure drop)."""


@dataclass
class PipeIntegration:
    """Pressure drops and end states of a batch of pipe integrations."""

    dp_forward: np.ndarray
    """``P_out - P_in`` integrated from the inlet, Pa (Eq. 9)."""
    dp_backward: np.ndarray
    """``P_out - P_in`` integrated from the outlet, Pa (Eq. 10)."""
    inlet: GradientState
    outlet: GradientState
    failed: np.ndarray
    forward_solution: BatchSolution | None = None

    @property
    def dp_mean(self) -> np.ndarray:
        return 0.5 * (self.dp_forward + self.dp_backward)


def integrate_pressure_drop(
    batch: PipeBatch,
    fluid: FluidProperties,
    pipe_idx: np.ndarray,
    p_in: np.ndarray,
    p_out: np.ndarray,
    q_l: np.ndarray,
    q_g: np.ndarray,
    options: IntegrationOptions | None = None,
    keep_trajectories: bool = False,
) -> PipeIntegration:
    """Integrate the momentum balance of each ``(pipe, p_in, p_out, q_l, q_g)`` problem.

    Forward and backward integrations of all problems run in a single batch.
    Pressures are in Pa.
    """
    options = options or IntegrationOptions()
    n = len(pipe_idx)
    idx2 = np.concatenate([pipe_idx, pipe_idx])
    p_edge = np.concatenate([p_in, p_out])
    ql2 = np.concatenate([q_l, q_l])
    qg2 = np.concatenate([q_g, q_g])
    direction = np.concatenate([np.ones(n), -np.ones(n)])

    def rhs(t, y, active):
        return pressure_gradient(
            batch,
            fluid,
            idx2[active],
            t,
            y,
            p_edge[active],
            ql2[active],
            qg2[active],
            direction[active],
            options.min_pressure,
        ).dpds

    sol = integrate_batch(
        rhs,
        np.zeros(2 * n),
        rtol=options.rtol,
        atol=options.atol,
        first_step=options.first_step,
        max_step=options.max_step,
        max_steps=options.max_steps,
        keep_trajectories=keep_trajectories,
    )
    dp_f = np.maximum(p_in + sol.y_end[:n], options.min_pressure) - p_in
    dp_b = p_out - np.maximum(p_out - sol.y_end[n:], options.min_pressure)
    failed = sol.failed[:n] | sol.failed[n:]
    dp_f = np.where(failed, np.nan, dp_f)
    dp_b = np.where(failed, np.nan, dp_b)
    zeros = np.zeros(n)
    inlet = pressure_gradient(
        batch, fluid, pipe_idx, zeros, zeros, p_in, q_l, q_g, np.ones(n), options.min_pressure
    )
    outlet = pressure_gradient(
        batch, fluid, pipe_idx, zeros, zeros, p_out, q_l, q_g, -np.ones(n), options.min_pressure
    )
    forward = None
    if keep_trajectories:
        forward = BatchSolution(
            y_end=sol.y_end[:n],
            trajectories=sol.trajectories[:n],
            steps=sol.steps[:n],
            failed=sol.failed[:n],
        )
    return PipeIntegration(dp_f, dp_b, inlet, outlet, failed, forward)


def midpoint_pressure_drop(
    batch: PipeBatch,
    fluid: FluidProperties,
    pipe_idx: np.ndarray,
    p_in: np.ndarray,
    p_out: np.ndarray,
    q_l: np.ndarray,
    q_g: np.ndarray,
    min_pressure: float = 1.0e5,
) -> tuple[np.ndarray, GradientState]:
    """Non-rigorous pressure drop: gradient at mid-length and mean pressure times the length."""
    n = len(pipe_idx)
    state = pressure_gradient(
        batch,
        fluid,
        pipe_idx,
        np.full(n, 0.5),
        np.zeros(n),
        0.5 * (p_in + p_out),
        q_l,
        q_g,
        np.ones(n),
        min_pressure,
    )
    return state.dpds, state


def profile_states(
    batch: PipeBatch,
    fluid: FluidProperties,
    integration: PipeIntegration,
    p_in: np.ndarray,
    q_l: np.ndarray,
    q_g: np.ndarray,
    n_points: int = 101,
    min_pressure: float = 1.0e5,
) -> tuple[np.ndarray, GradientState]:
    """Flow state at ``n_points`` equally spaced positions of each pipe (from the forward pass)."""
    if integration.forward_solution is None:
        raise ValueError("integration must be run with keep_trajectories=True")
    n = len(p_in)
    s = np.linspace(0.0, 1.0, n_points)
    dp = integration.forward_solution.interpolate(s)  # (n, n_points)
    pipe_idx = np.repeat(np.arange(n), n_points)
    state = pressure_gradient(
        batch,
        fluid,
        pipe_idx,
        np.tile(s, n),
        dp.ravel(),
        np.repeat(p_in, n_points),
        np.repeat(q_l, n_points),
        np.repeat(q_g, n_points),
        np.ones(n * n_points),
        min_pressure,
    )
    return s, state
