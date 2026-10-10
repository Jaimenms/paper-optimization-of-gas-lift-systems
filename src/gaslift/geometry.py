"""Pipe geometry: 3-D polylines with axially varying diameter, roughness and inclination.

A pipe is described by a sequence of vertices. Each straight segment between two
consecutive vertices has a constant inclination angle, diameter and roughness.
To keep the pressure-gradient integrand continuous, the segment properties are
joined by short linear transitions of half-width ``transition_length`` around
each interior vertex (the ``tinyl`` trick of the original MATLAB code).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

Point = tuple[float, float, float]


def _as_array(values, n: int, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        arr = np.full(n, float(arr))
    if arr.shape != (n,):
        raise ValueError(f"{name} must be a scalar or have one value per vertex ({n})")
    return arr


@dataclass(frozen=True)
class PipeGeometry:
    """Geometry of a single pipe (SI units)."""

    points: np.ndarray
    """Vertices of the pipe axis, shape ``(n, 3)``, ordered from inlet to outlet."""
    diameters: np.ndarray
    """Internal diameter at each vertex, m."""
    roughness: np.ndarray
    """Absolute roughness at each vertex, m."""
    equivalent_length: float = 0.0
    """Equivalent length of fittings and accidents, m."""
    transition_length: float = 1.0
    """Half-width of the property transition around interior vertices, m."""

    def __post_init__(self) -> None:
        points = np.atleast_2d(np.asarray(self.points, dtype=float))
        if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 2:
            raise ValueError("points must be an (n >= 2, 3) array of x, y, z coordinates")
        n = points.shape[0]
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "diameters", _as_array(self.diameters, n, "diameters"))
        object.__setattr__(self, "roughness", _as_array(self.roughness, n, "roughness"))
        if np.any(self.diameters <= 0):
            raise ValueError("diameters must be positive")
        if np.any(self.segment_lengths <= 0):
            raise ValueError("consecutive pipe vertices must not coincide")

    # ------------------------------------------------------------------ constructors
    @classmethod
    def from_points(
        cls,
        points,
        diameter,
        roughness=0.0,
        equivalent_length: float = 0.0,
        transition_length: float = 1.0,
    ) -> PipeGeometry:
        """Build a geometry from vertices and scalar-or-per-vertex properties."""
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        n = pts.shape[0]
        return cls(
            points=pts,
            diameters=_as_array(diameter, n, "diameter"),
            roughness=_as_array(roughness, n, "roughness"),
            equivalent_length=float(equivalent_length),
            transition_length=float(transition_length),
        )

    @classmethod
    def straight(
        cls,
        start: Point,
        end: Point,
        diameter: float,
        roughness: float = 0.0,
        equivalent_length: float = 0.0,
    ) -> PipeGeometry:
        """Single straight segment between two points."""
        return cls.from_points([start, end], diameter, roughness, equivalent_length)

    # ------------------------------------------------------------------ derived data
    @cached_property
    def segment_vectors(self) -> np.ndarray:
        return np.diff(self.points, axis=0)

    @cached_property
    def segment_lengths(self) -> np.ndarray:
        return np.linalg.norm(self.segment_vectors, axis=1)

    @cached_property
    def length(self) -> float:
        """Total developed length of the pipe, m."""
        return float(self.segment_lengths.sum())

    @cached_property
    def cumulative_lengths(self) -> np.ndarray:
        """Axial position of each vertex, m (first vertex at 0)."""
        return np.concatenate([[0.0], np.cumsum(self.segment_lengths)])

    @cached_property
    def segment_angles(self) -> np.ndarray:
        """Inclination of each segment with respect to the horizontal plane, rad.

        Positive values mean the segment rises along the flow (digraph) direction.
        """
        dz = self.segment_vectors[:, 2]
        return np.arcsin(np.clip(dz / self.segment_lengths, -1.0, 1.0))

    @cached_property
    def segment_diameters(self) -> np.ndarray:
        return 0.5 * (self.diameters[1:] + self.diameters[:-1])

    @cached_property
    def segment_roughness(self) -> np.ndarray:
        return 0.5 * (self.roughness[1:] + self.roughness[:-1])

    @cached_property
    def elevation_change(self) -> float:
        """Outlet elevation minus inlet elevation, m."""
        return float(self.points[-1, 2] - self.points[0, 2])

    @cached_property
    def knots(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Breakpoints ``(s, angle, diameter, roughness)`` for piecewise-linear interpolation.

        Each segment contributes two knots located ``transition_length`` inside its
        ends (except at the pipe extremities), so segment properties are constant
        along the segment and vary linearly across interior vertices.
        """
        n_seg = len(self.segment_lengths)
        cum = self.cumulative_lengths
        t = self.transition_length
        s = np.empty(2 * n_seg)
        s[0::2] = cum[:-1] + t
        s[1::2] = cum[1:] - t
        s[0] = cum[0]
        s[-1] = cum[-1]
        return (
            s,
            np.repeat(self.segment_angles, 2),
            np.repeat(self.segment_diameters, 2),
            np.repeat(self.segment_roughness, 2),
        )

    @property
    def inlet_diameter(self) -> float:
        return float(self.knots[2][0])

    @property
    def outlet_diameter(self) -> float:
        return float(self.knots[2][-1])

    @property
    def mean_knot_diameter(self) -> float:
        """Average of the knot diameters (used by the linear initial-estimate model)."""
        return float(self.knots[2].mean())

    @property
    def start(self) -> np.ndarray:
        return self.points[0]

    @property
    def end(self) -> np.ndarray:
        return self.points[-1]

    # ------------------------------------------------------------------ evaluation
    def properties_at(self, s) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Inclination (rad), diameter (m) and roughness (m) at axial positions ``s`` (m)."""
        s = np.asarray(s, dtype=float)
        ks, ka, kd, kr = self.knots
        return np.interp(s, ks, ka), np.interp(s, ks, kd), np.interp(s, ks, kr)

    def with_uniform_diameter(self, diameter: float) -> PipeGeometry:
        """Copy of this geometry with a single diameter along the whole pipe."""
        return PipeGeometry(
            points=self.points.copy(),
            diameters=np.full(len(self.points), float(diameter)),
            roughness=self.roughness.copy(),
            equivalent_length=self.equivalent_length,
            transition_length=self.transition_length,
        )

    def reversed(self) -> PipeGeometry:
        """Same pipe traversed from outlet to inlet."""
        return PipeGeometry(
            points=self.points[::-1].copy(),
            diameters=self.diameters[::-1].copy(),
            roughness=self.roughness[::-1].copy(),
            equivalent_length=self.equivalent_length,
            transition_length=self.transition_length,
        )

    # ------------------------------------------------------------------ serialisation
    def to_dict(self) -> dict:
        return {
            "points": self.points.tolist(),
            "diameters": self.diameters.tolist(),
            "roughness": self.roughness.tolist(),
            "equivalent_length": self.equivalent_length,
            "transition_length": self.transition_length,
        }

    @classmethod
    def from_dict(cls, data: dict) -> PipeGeometry:
        return cls(
            points=np.asarray(data["points"], dtype=float),
            diameters=np.asarray(data["diameters"], dtype=float),
            roughness=np.asarray(data["roughness"], dtype=float),
            equivalent_length=float(data.get("equivalent_length", 0.0)),
            transition_length=float(data.get("transition_length", 1.0)),
        )


def annulus_equivalent_diameter(
    outer_diameter: float, inner_diameter: float, model: str = "matlab"
) -> float:
    """Equivalent diameter used to treat a casing-tubing annulus as a circular pipe.

    ``model="hydraulic"`` returns the textbook hydraulic diameter
    ``D_outer - D_inner``. ``model="matlab"`` reproduces the expression of the
    original bbnet10 code, ``4 (A_outer + A_inner) / (P_outer + P_inner)``, which
    equals ``(D_outer² + D_inner²) / (D_outer + D_inner)``.
    """
    if model == "hydraulic":
        return outer_diameter - inner_diameter
    if model == "matlab":
        return (outer_diameter**2 + inner_diameter**2) / (outer_diameter + inner_diameter)
    raise ValueError(f"unknown annulus diameter model {model!r}")
