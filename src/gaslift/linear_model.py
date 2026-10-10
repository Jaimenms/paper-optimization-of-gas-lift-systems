"""Linear laminar model used to generate the initial estimate (paper Section 3.4).

Each phase is assumed to flow in laminar regime without interaction with the
other phase, occupying the fraction of the cross-section given by the pipe's
``initial_liquid_holdup``. The hydrostatic head uses the same hold-up. The
resulting linear system ``A y = c`` (Eqs. 21-23) is solved directly.
"""

from __future__ import annotations

import numpy as np

from .model import NetworkModel
from .units import G


def linear_initial_estimate(model: NetworkModel) -> np.ndarray:
    """Solve the linear network model and return a (scaled) state vector."""
    net, fluid, lay = model.network, model.fluid, model.layout
    n, s = lay.n_nodes, lay.n_pipes
    pr = model.pressure_scale
    m, e_mat, e_vec = model.M, model.E, model.e

    length = net.pipe_lengths
    diameter = np.array([p.geometry.mean_knot_diameter for p in net.pipes])
    hl0 = np.array([p.initial_liquid_holdup for p in net.pipes])
    rho_l, rho_g = fluid.liquid_density, fluid.gas_reference_density

    f_l = (
        np.diag(128.0 * hl0 * fluid.liquid_viscosity * length / (diameter**4 * np.pi)) / rho_l / pr
    )
    f_g = (
        np.diag(128.0 * (1.0 - hl0) * fluid.gas_viscosity * length / (diameter**4 * np.pi))
        / rho_g
        / pr
    )

    z_ns, z_nn, z_ss, z_sn = np.zeros((n, s)), np.zeros((n, n)), np.zeros((s, s)), np.zeros((s, n))
    i_n, i_s = np.eye(n), np.eye(s)
    m1, m2 = (m < 0).astype(float), (m > 0).astype(float)
    a = np.vstack(
        [
            np.hstack([m, z_ns, z_ns, z_ns, z_nn, i_n, z_nn]),
            np.hstack([z_ns, m, z_ns, z_ns, z_nn, z_nn, i_n]),
            np.hstack([-f_l, -f_g, i_s, -i_s, z_sn, z_sn, z_sn]),
            np.hstack([z_ss, z_ss, z_ss, -i_s, m2.T, z_sn, z_sn]),
            np.hstack([z_ss, z_ss, -i_s, z_ss, m1.T, z_sn, z_sn]),
            e_mat,
        ]
    )
    head = (rho_g * G / pr * (1.0 - hl0) + rho_l * G / pr * hl0) * (m.T @ net.elevations)
    c = np.concatenate([np.zeros(n), np.zeros(n), head, np.zeros(2 * s), e_vec])
    return np.linalg.solve(a, c)
