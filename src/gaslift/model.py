"""Non-linear two-phase pipe-network model (paper Section 3.3).

Residual vector ``F(y) = 0`` with, in order: nodal liquid mass balances (Eq. 6),
nodal gas mass balances (Eq. 7), pipe momentum balances (Eq. 8), nodal momentum
balances at pipe outlets (Eq. 12) and inlets (Eq. 11), and the specification
equations (Eq. 14). Pressures in the state vector are scaled by
``ModelOptions.pressure_scale`` (10 MPa by default) so that all unknowns are O(1).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .network import Network, Specifications, StateLayout
from .pipe_model import (
    IntegrationOptions,
    PipeBatch,
    PipeIntegration,
    integrate_pressure_drop,
    midpoint_pressure_drop,
)


@dataclass(frozen=True)
class ModelOptions:
    """Numerical settings of the network model."""

    pressure_scale: float = 1.0e7
    """Pressure unit of the scaled state vector, Pa."""
    min_pressure: float = 1.0e5
    """Lowest admissible pressure, Pa."""
    integration: IntegrationOptions = field(default_factory=IntegrationOptions)
    jacobian_relative_step: float = 1e-3
    """Relative perturbation of the pipe variables in the finite-difference Jacobian blocks
    (value of the original code; robust far from the solution)."""
    jacobian_absolute_step: float = 1e-6
    jacobian_fine_relative_step: float = 1e-5
    """Smaller perturbation used once the residual norm drops below
    ``jacobian_fine_residual``. It keeps the perturbed integrations on the same adaptive
    step sequence as the base one, which restores quadratic convergence near the solution."""
    jacobian_fine_absolute_step: float = 1e-8
    jacobian_fine_residual: float = 1e-4
    rigorous: bool = True
    """Integrate the pressure gradient along the pipes (True) or evaluate it at mid-length."""
    profile_points: int = 101


@dataclass
class StateView:
    """Named views of a (scaled) state vector."""

    qL: np.ndarray
    qG: np.ndarray
    PI: np.ndarray
    PO: np.ndarray
    PW: np.ndarray
    wL: np.ndarray
    wG: np.ndarray


class NetworkModel:
    """Assembles residuals and Jacobian of the two-phase network equations."""

    def __init__(
        self,
        network: Network,
        specifications: Specifications,
        options: ModelOptions | None = None,
    ) -> None:
        self.network = network
        self.specifications = specifications
        self.options = options or ModelOptions()
        self.fluid = network.fluid
        self.layout: StateLayout = network.layout
        self.M = network.incidence_matrix()
        self.M1 = (self.M < 0).astype(float)
        self.M2 = (self.M > 0).astype(float)
        self.E, self.e = specifications.matrix(self.options.pressure_scale)
        self.batch = PipeBatch.from_network(network)
        self.pipe_index = np.arange(self.layout.n_pipes)
        self.lower_bounds, self.upper_bounds = self._state_bounds()
        self._integration = IntegrationOptions(
            rtol=self.options.integration.rtol,
            atol=self.options.integration.atol,
            first_step=self.options.integration.first_step,
            max_step=self.options.integration.max_step,
            min_pressure=self.options.min_pressure,
            max_steps=self.options.integration.max_steps,
        )

    # ------------------------------------------------------------------ helpers
    @property
    def size(self) -> int:
        return self.layout.size

    @property
    def pressure_scale(self) -> float:
        return self.options.pressure_scale

    def split(self, x: np.ndarray) -> StateView:
        lay = self.layout
        return StateView(
            x[lay.qL], x[lay.qG], x[lay.PI], x[lay.PO], x[lay.PW], x[lay.wL], x[lay.wG]
        )

    def _state_bounds(self) -> tuple[np.ndarray, np.ndarray]:
        lay = self.layout
        lb = np.full(lay.size, -np.inf)
        ub = np.full(lay.size, np.inf)
        p_min = self.options.min_pressure / self.options.pressure_scale
        for sl in (lay.PI, lay.PO, lay.PW):
            lb[sl] = p_min
        lb[lay.qL] = 0.0
        lb[lay.qG] = 0.0
        return lb, ub

    # ------------------------------------------------------------------ momentum balances
    def _momentum(
        self,
        qL: np.ndarray,
        qG: np.ndarray,
        PI: np.ndarray,
        PO: np.ndarray,
        PW: np.ndarray,
        rigorous: bool,
        keep_trajectories: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, PipeIntegration | None]:
        """Momentum residuals for ``K`` variants of the pipe variables (arrays ``(K, S)``)."""
        pr = self.pressure_scale
        k, s = qL.shape
        idx = np.tile(self.pipe_index, k)
        ql, qg = qL.ravel(), qG.ravel()
        p_in, p_out = PI.ravel() * pr, PO.ravel() * pr
        integration = None
        if rigorous:
            integration = integrate_pressure_drop(
                self.batch,
                self.fluid,
                idx,
                p_in,
                p_out,
                ql,
                qg,
                self._integration,
                keep_trajectories,
            )
            dp = integration.dp_mean
            inlet, outlet = integration.inlet, integration.outlet
        else:
            dp, state = midpoint_pressure_drop(
                self.batch, self.fluid, idx, p_in, p_out, ql, qg, self.options.min_pressure
            )
            inlet = outlet = state
        rho_l = self.fluid.liquid_density
        a_in = np.tile(self.batch.inlet_area, k)
        a_out = np.tile(self.batch.outlet_area, k)
        kin_out = (
            np.abs(outlet.liquid_velocity) * a_out * ql / (2.0 * rho_l**2)
            + np.abs(outlet.gas_velocity) * a_out * qg / (2.0 * outlet.gas_density**2)
        ) / pr
        kin_in = (
            np.abs(inlet.liquid_velocity) * a_in * ql / (2.0 * rho_l**2)
            + np.abs(inlet.gas_velocity) * a_in * qg / (2.0 * inlet.gas_density**2)
        ) / pr
        f3 = PO - PI - (dp / pr).reshape(k, s)
        f41 = (self.M2.T @ PW)[None, :] - PO - kin_out.reshape(k, s)
        f42 = (self.M1.T @ PW)[None, :] - PI - kin_in.reshape(k, s)
        return f3, f41, f42, integration

    # ------------------------------------------------------------------ public API
    def residuals(self, x: np.ndarray, rigorous: bool | None = None) -> np.ndarray:
        """Residual vector ``F(x)``."""
        rig = self.options.rigorous if rigorous is None else rigorous
        v = self.split(x)
        f1 = self.M @ v.qL + v.wL
        f2 = self.M @ v.qG + v.wG
        f3, f41, f42, _ = self._momentum(v.qL[None], v.qG[None], v.PI[None], v.PO[None], v.PW, rig)
        f5 = self.E @ x - self.e
        return np.concatenate([f1, f2, f3[0], f41[0], f42[0], f5])

    def residuals_and_jacobian(
        self, x: np.ndarray, rigorous: bool | None = None, fine: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        """Residual vector and Jacobian matrix.

        The mass-balance, nodal-momentum-vs-pressure and specification blocks are
        analytical; the derivatives of the momentum balances with respect to the pipe
        variables are diagonal blocks obtained by forward finite differences, with all
        perturbed problems integrated in one batch.
        """
        rig = self.options.rigorous if rigorous is None else rigorous
        lay = self.layout
        n, s = lay.n_nodes, lay.n_pipes
        v = self.split(x)
        if fine:
            rel = self.options.jacobian_fine_relative_step
            absd = self.options.jacobian_fine_absolute_step
        else:
            rel, absd = self.options.jacobian_relative_step, self.options.jacobian_absolute_step

        def perturb(a):
            return a * (1.0 + rel) + absd

        ql_p, qg_p, pi_p, po_p = perturb(v.qL), perturb(v.qG), perturb(v.PI), perturb(v.PO)
        qL = np.vstack([v.qL, ql_p, v.qL, v.qL, v.qL])
        qG = np.vstack([v.qG, v.qG, qg_p, v.qG, v.qG])
        PI = np.vstack([v.PI, v.PI, v.PI, pi_p, v.PI])
        PO = np.vstack([v.PO, v.PO, v.PO, v.PO, po_p])
        f3, f41, f42, _ = self._momentum(qL, qG, PI, PO, v.PW, rig)

        f1 = self.M @ v.qL + v.wL
        f2 = self.M @ v.qG + v.wG
        f5 = self.E @ x - self.e
        F = np.concatenate([f1, f2, f3[0], f41[0], f42[0], f5])

        inv = [1.0 / (ql_p - v.qL), 1.0 / (qg_p - v.qG), 1.0 / (pi_p - v.PI), 1.0 / (po_p - v.PO)]
        J = np.zeros((lay.size, lay.size))
        r1, r2 = slice(0, n), slice(n, 2 * n)
        r3, r41, r42 = (
            slice(2 * n, 2 * n + s),
            slice(2 * n + s, 2 * n + 2 * s),
            slice(2 * n + 2 * s, 2 * n + 3 * s),
        )
        r5 = slice(2 * n + 3 * s, lay.size)
        cols = (lay.qL, lay.qG, lay.PI, lay.PO)
        J[r1, lay.qL] = self.M
        J[r1, lay.wL] = np.eye(n)
        J[r2, lay.qG] = self.M
        J[r2, lay.wG] = np.eye(n)
        diag = np.arange(s)
        for rows, f in ((r3, f3), (r41, f41), (r42, f42)):
            for k, col in enumerate(cols):
                block = np.zeros((s, s))
                block[diag, diag] = (f[k + 1] - f[0]) * inv[k]
                J[rows, col] = block
        J[r41, lay.PW] = self.M2.T
        J[r42, lay.PW] = self.M1.T
        J[r5, :] = self.E
        return F, J

    def evaluate(self, x: np.ndarray) -> tuple[np.ndarray, PipeIntegration]:
        """Residuals plus the rigorous pipe integration with trajectories (for profiles)."""
        v = self.split(x)
        f1 = self.M @ v.qL + v.wL
        f2 = self.M @ v.qG + v.wG
        f3, f41, f42, integration = self._momentum(
            v.qL[None], v.qG[None], v.PI[None], v.PO[None], v.PW, True, keep_trajectories=True
        )
        f5 = self.E @ x - self.e
        return np.concatenate([f1, f2, f3[0], f41[0], f42[0], f5]), integration

    def clip_to_bounds(self, x: np.ndarray) -> np.ndarray:
        return np.clip(x, self.lower_bounds, self.upper_bounds)
