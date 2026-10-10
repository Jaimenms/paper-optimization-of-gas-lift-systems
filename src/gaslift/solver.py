"""Newton-Raphson solver for the network equations (paper Section 3.5)."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from .model import NetworkModel


@dataclass(frozen=True)
class NewtonOptions:
    tolerance: float = 1e-9
    """Convergence criterion on the Euclidean norm of the residual vector."""
    max_iterations: int = 15
    step_tolerance: float = 0.0
    """Stop when the Newton step norm falls below this value (0 disables the test)."""
    verbose: bool = False


@dataclass
class SolverStatus:
    converged: bool
    iterations: int
    residual_norm: float
    reason: str
    elapsed: float = 0.0
    rigorous: bool = True

    def __str__(self) -> str:
        state = "converged" if self.converged else "FAILED"
        return (
            f"{state} ({self.reason}) after {self.iterations} iterations, "
            f"|F| = {self.residual_norm:.3e}, {self.elapsed:.2f} s"
        )


def newton_solve(
    model: NetworkModel,
    x0: np.ndarray,
    options: NewtonOptions | None = None,
    rigorous: bool | None = None,
) -> tuple[np.ndarray, SolverStatus]:
    """Solve ``F(x) = 0`` by Newton-Raphson with bound clipping of the iterates.

    Returns the final iterate and a :class:`SolverStatus`. Failures (non-finite
    values, singular Jacobian, iteration limit) are reported in the status rather
    than raised.
    """
    opts = options or NewtonOptions()
    rig = model.options.rigorous if rigorous is None else rigorous
    t0 = time.perf_counter()
    x = model.clip_to_bounds(np.array(x0, dtype=float))
    res = np.inf
    for it in range(opts.max_iterations + 1):
        try:
            fine = res < model.options.jacobian_fine_residual
            f, jac = model.residuals_and_jacobian(x, rig, fine=fine)
        except FloatingPointError as exc:  # pragma: no cover - defensive
            return x, SolverStatus(
                False, it, res, f"evaluation error: {exc}", time.perf_counter() - t0, rig
            )
        if not np.all(np.isfinite(f)):
            return x, SolverStatus(
                False, it, np.inf, "non-finite residuals", time.perf_counter() - t0, rig
            )
        res = float(np.linalg.norm(f))
        if opts.verbose:
            print(f"  it {it:2d}  |F| = {res:.3e}")
        if res <= opts.tolerance:
            return x, SolverStatus(
                True, it, res, "residual below tolerance", time.perf_counter() - t0, rig
            )
        if it == opts.max_iterations:
            return x, SolverStatus(
                False, it, res, "iteration limit reached", time.perf_counter() - t0, rig
            )
        try:
            dx = np.linalg.solve(jac, f)
        except np.linalg.LinAlgError:
            return x, SolverStatus(
                False, it, res, "singular Jacobian", time.perf_counter() - t0, rig
            )
        if not np.all(np.isfinite(dx)):
            return x, SolverStatus(
                False, it, res, "non-finite Newton step", time.perf_counter() - t0, rig
            )
        x = model.clip_to_bounds(x - dx)
        if opts.step_tolerance > 0 and np.linalg.norm(dx) < opts.step_tolerance:
            f = model.residuals(x, rig)
            res = float(np.linalg.norm(f))
            return x, SolverStatus(
                res <= opts.tolerance,
                it + 1,
                res,
                "step below tolerance",
                time.perf_counter() - t0,
                rig,
            )
    return x, SolverStatus(
        False, opts.max_iterations, res, "iteration limit reached", time.perf_counter() - t0, rig
    )
