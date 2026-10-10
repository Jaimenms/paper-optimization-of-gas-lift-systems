"""Beggs & Brill (1973) two-phase flow correlation.

Vectorised implementation of the closure relations used by the pipe momentum
balance: flow-pattern identification, liquid hold-up for any inclination,
Churchill single-phase friction factor and the Beggs & Brill two-phase friction
multiplier. All ``if`` tests of the original correlation are replaced by very
sharp logistic functions, exactly as in the MATLAB code that produced the paper
results, so that the model is continuous in its arguments.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from .fluids import FluidProperties
from .units import G

#: Sharpness of the logistic "soft if" used for flow pattern selection.
SIGMOID_SHARPNESS = 1.0e6
#: Minimum hold-up of either phase (numerical guard).
HOLDUP_FLOOR = 1.0e-6
#: Lowest pressure used inside the correlation, Pa.
MIN_PRESSURE = 1.0e5

SEGREGATED, TRANSITION, INTERMITTENT, DISTRIBUTED = 1, 2, 3, 4
PATTERN_NAMES = {
    SEGREGATED: "segregated",
    TRANSITION: "transition",
    INTERMITTENT: "intermittent",
    DISTRIBUTED: "distributed",
}


@dataclass(frozen=True)
class FlowPatternWeights:
    """Continuous (0..1) membership of each Beggs & Brill flow pattern."""

    segregated: np.ndarray
    transition: np.ndarray
    intermittent: np.ndarray
    distributed: np.ndarray

    @property
    def index(self) -> np.ndarray:
        """Weighted pattern index: 1 segregated, 2 transition, 3 intermittent, 4 distributed."""
        return (
            SEGREGATED * self.segregated
            + TRANSITION * self.transition
            + INTERMITTENT * self.intermittent
            + DISTRIBUTED * self.distributed
        )


@dataclass(frozen=True)
class ClosureResult:
    """Outputs of :func:`closure_relations` (all arrays share the input broadcast shape)."""

    pressure: np.ndarray
    shear_stress: np.ndarray
    """Wall shear stress of the mixture, Pa (signed with the flow direction)."""
    mixture_density: np.ndarray
    """Slip mixture density ``HL rho_L + HG rho_G``, kg/m³."""
    liquid_velocity: np.ndarray
    gas_velocity: np.ndarray
    liquid_holdup: np.ndarray
    gas_holdup: np.ndarray
    no_slip_holdup: np.ndarray
    froude: np.ndarray
    liquid_velocity_number: np.ndarray
    reynolds: np.ndarray
    friction_factor: np.ndarray
    """Single-phase Darcy friction factor (Churchill)."""
    two_phase_friction_factor: np.ndarray
    gas_density: np.ndarray
    pattern: FlowPatternWeights


def soft_greater(x1, x2, sharpness: float = SIGMOID_SHARPNESS):
    """Logistic approximation of ``x1 > x2`` (0.5 when equal)."""
    return expit(sharpness * (np.asarray(x1, float) - np.asarray(x2, float)))


def pattern_boundaries(lambda_l):
    """Froude-number boundaries L1..L4 of the Beggs & Brill horizontal flow map."""
    lambda_l = np.asarray(lambda_l, dtype=float)
    with np.errstate(divide="ignore", over="ignore"):
        l1 = 316.0 * lambda_l**0.3020
        l2 = 0.0009252 * lambda_l ** (-2.4684)
        l3 = 0.10 * lambda_l ** (-1.4516)
        l4 = 0.5 * lambda_l ** (-6.7380)
    return l1, l2, l3, l4


def flow_pattern(lambda_l, froude) -> tuple[FlowPatternWeights, tuple]:
    """Flow pattern weights for the given no-slip hold-up and Froude number."""
    l1, l2, l3, l4 = pattern_boundaries(lambda_l)
    a = soft_greater(lambda_l, 0.01)  # lambda > 0.01
    b = soft_greater(lambda_l, 0.4)  # lambda > 0.4
    c = soft_greater(l1, froude)  # Fr < L1
    d = soft_greater(l2, froude)  # Fr < L2
    e = soft_greater(l3, froude)  # Fr < L3
    f = soft_greater(l4, froude)  # Fr < L4

    segregated = (1 - a) * c + a * d
    transition = a * e * (1 - d)
    intermittent = (1 - b) * a * c * (1 - e) + b * f * (1 - e)
    distributed = (1 - b) * (1 - c) + b * (1 - f)
    clip = lambda w: np.clip(w, 0.0, 1.0)  # noqa: E731
    return (
        FlowPatternWeights(
            clip(segregated), clip(transition), clip(intermittent), clip(distributed)
        ),
        (l1, l2, l3, l4),
    )


def horizontal_holdup(lambda_l, froude):
    """Horizontal liquid hold-up for the segregated, intermittent and distributed patterns."""
    out = []
    with np.errstate(divide="ignore", invalid="ignore"):
        for a, b, c in ((0.98, 0.4846, 0.0868), (0.845, 0.5351, 0.0173), (1.065, 0.5824, 0.0609)):
            hl0 = a * lambda_l**b / froude**c
            hl0 = np.maximum(hl0, lambda_l)
            out.append(np.minimum(hl0, 1.0))
    return tuple(out)


def _inclination_correction(lambda_l, froude, n_lv, angle, d, e, f, g):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        c = (1.0 - lambda_l) * np.log(d * lambda_l**e * n_lv**f * froude**g)
    c = np.where(np.isfinite(c), c, 0.0)
    c = np.maximum(c, 0.0)
    s = np.sin(1.8 * angle)
    return 1.0 + c * (s - 0.333 * s**3)


def inclination_factors(lambda_l, froude, n_lv, angle, direction):
    """Beggs & Brill inclination correction for each pattern (uphill/downhill aware)."""
    uphill = 0.5 * (1.0 + np.sign(direction) * np.sign(angle))
    downhill = 1.0 - uphill
    psi_down = _inclination_correction(
        lambda_l, froude, n_lv, angle, 4.70, -0.3692, 0.1244, -0.5056
    )
    psi_seg = _inclination_correction(lambda_l, froude, n_lv, angle, 0.011, -3.768, 3.539, -1.614)
    psi_int = _inclination_correction(lambda_l, froude, n_lv, angle, 2.96, 0.305, -0.4473, 0.0978)
    return (
        psi_seg * uphill + psi_down * downhill,
        psi_int * uphill + psi_down * downhill,
        1.0 * uphill + psi_down * downhill,
    )


def liquid_holdup(lambda_l, froude, n_lv, angle, direction):
    """Liquid hold-up for any inclination (pattern-weighted), clipped to ``[0, 1]``."""
    weights, (_, l2, l3, _) = flow_pattern(lambda_l, froude)
    hl0_seg, hl0_int, hl0_dis = horizontal_holdup(lambda_l, froude)
    psi_seg, psi_int, psi_dis = inclination_factors(lambda_l, froude, n_lv, angle, direction)
    hl_seg = hl0_seg * psi_seg
    hl_int = hl0_int * psi_int
    hl_dis = hl0_dis * psi_dis
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = (l3 - froude) / (l3 - l2)
    hl_tra = phi * hl_seg + (1.0 - phi) * hl_int
    hl = (
        hl_seg * weights.segregated
        + hl_int * weights.intermittent
        + hl_dis * weights.distributed
        + hl_tra * weights.transition
    )
    return np.clip(hl, 0.0, 1.0), weights


def churchill_friction_factor(reynolds, relative_roughness):
    """Darcy friction factor from Churchill (1977), valid for all flow regimes."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        a = (2.457 * np.log(1.0 / (0.27 * relative_roughness + (7.0 / reynolds) ** 0.9))) ** 16
        b = (37530.0 / reynolds) ** 16
        fanning = 2.0 * ((8.0 / reynolds) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)
    return 4.0 * fanning


def two_phase_friction_multiplier(lambda_l, holdup, extreme):
    """``exp(S)`` of Beggs & Brill, with ``S = 0`` at the single-phase limits."""
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.log(lambda_l / holdup**2)
        s = np.log(2.2 * np.exp(k) - 1.2)
        s_poly = 1.0 / (-0.0523 / k + 3.1820 - 0.8725 * k + 0.01853 * k**3)
    use_poly = (k <= 0.0) | (k >= 0.1823)
    s = np.where(use_poly, s_poly, s)
    s = np.where(extreme, 0.0, s)
    return np.exp(s)


def closure_relations(
    pressure,
    v_sl,
    v_sg,
    diameter,
    roughness,
    angle,
    fluid: FluidProperties,
) -> ClosureResult:
    """Evaluate the Beggs & Brill closure relations.

    Parameters
    ----------
    pressure : array_like
        Absolute pressure, Pa.
    v_sl, v_sg : array_like
        Superficial velocities of liquid and gas, m/s (signed along the pipe axis).
    diameter, roughness : array_like
        Pipe diameter and absolute roughness, m.
    angle : array_like
        Pipe inclination with respect to the horizontal, rad (positive = rising).
    fluid : FluidProperties
    """
    pressure = np.maximum(np.asarray(pressure, dtype=float), MIN_PRESSURE)
    v_sl = np.asarray(v_sl, dtype=float)
    v_sg = np.asarray(v_sg, dtype=float)
    diameter = np.asarray(diameter, dtype=float)
    roughness = np.asarray(roughness, dtype=float)
    angle = np.asarray(angle, dtype=float)

    direction = np.where(v_sl > 0.0, 1.0, -1.0)
    v_sl = np.abs(v_sl)
    v_sg = np.abs(v_sg)

    rho_l = fluid.liquid_density
    rho_g = fluid.gas_density(pressure)
    v_m = v_sl + v_sg

    lambda_l = v_sl / (v_m + 1e-9)
    lambda_g = v_sg / (v_m + 1e-9)
    over = lambda_l > 1.0
    lambda_l = np.where(over, 1.0, lambda_l)
    lambda_g = np.where(over, 0.0, lambda_g)
    over = lambda_g > 1.0
    lambda_g = np.where(over, 1.0, lambda_g)
    lambda_l = np.where(over, 0.0, lambda_l)

    rho_ns = rho_l * lambda_l + rho_g * lambda_g
    mu_ns = fluid.liquid_viscosity * lambda_l + fluid.gas_viscosity * lambda_g

    froude = v_m**2 / (G * diameter)
    n_lv = v_sl * (rho_l / fluid.interfacial_tension / G) ** 0.25
    reynolds = rho_ns * v_m * diameter / mu_ns

    hl, weights = liquid_holdup(lambda_l, froude, n_lv, angle, direction)
    hg = 1.0 - hl
    gas_limit = hg < HOLDUP_FLOOR
    hg = np.where(gas_limit, HOLDUP_FLOOR, hg)
    hl = np.where(gas_limit, 1.0 - HOLDUP_FLOOR, hl)
    liquid_limit = hl < HOLDUP_FLOOR
    hl = np.where(liquid_limit, HOLDUP_FLOOR, hl)
    hg = np.where(liquid_limit, 1.0 - HOLDUP_FLOOR, hg)

    v_l = direction * v_sl / hl
    v_g = direction * v_sg / hg
    rho_s = hl * rho_l + hg * rho_g

    f_sp = churchill_friction_factor(reynolds, roughness / diameter)
    f_tp = f_sp * two_phase_friction_multiplier(lambda_l, hl, gas_limit | liquid_limit)
    tau = direction * f_tp * rho_s * v_m**2 / 8.0

    return ClosureResult(
        pressure=pressure,
        shear_stress=tau,
        mixture_density=rho_s,
        liquid_velocity=v_l,
        gas_velocity=v_g,
        liquid_holdup=hl,
        gas_holdup=hg,
        no_slip_holdup=lambda_l,
        froude=froude,
        liquid_velocity_number=n_lv,
        reynolds=reynolds,
        friction_factor=f_sp,
        two_phase_friction_factor=f_tp,
        gas_density=rho_g,
        pattern=weights,
    )
