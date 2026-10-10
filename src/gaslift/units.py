"""Physical constants and unit conversion factors (SI base units everywhere)."""

from __future__ import annotations

#: Standard gravity used throughout the model, m/s².
G = 9.80665

#: Universal gas constant, kJ/(kmol K).
R_GAS = 8.314

#: One inch in metres.
INCH = 0.0254

#: One bar in pascal.
BAR = 1.0e5

#: One oil barrel in cubic metres.
BARREL = 0.1589873

#: Seconds in a (non-leap) year.
SECONDS_PER_YEAR = 365.0 * 24.0 * 3600.0

#: Kilowatt to brake horsepower.
KW_TO_HP = 1.34102209


def inch_to_m(value: float) -> float:
    """Convert inches to metres."""
    return value * INCH


def m_to_inch(value: float) -> float:
    """Convert metres to inches."""
    return value / INCH


def pa_to_bar(value):
    """Convert pascal to bar (works on scalars and arrays)."""
    return value / BAR


def bar_to_pa(value):
    """Convert bar to pascal (works on scalars and arrays)."""
    return value * BAR
