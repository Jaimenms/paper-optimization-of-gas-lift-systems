"""Batched explicit Runge-Kutta integrator (Dormand-Prince 5(4)).

Integrates many *independent* scalar initial-value problems simultaneously,
each with its own adaptive step size. The step-size control mirrors MATLAB's
``ode45`` (error estimate relative to ``max(|y|, atol/rtol)``, factors 0.8/1.25,
rejection factors 0.1 and 0.5, last step stretched by up to 10 %), which is the
solver the original gas-lift model was built with.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

# Dormand-Prince tableau (MATLAB ode45 ordering).
_C = np.array([0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
_A = np.array(
    [
        [0, 0, 0, 0, 0, 0],
        [1 / 5, 0, 0, 0, 0, 0],
        [3 / 40, 9 / 40, 0, 0, 0, 0],
        [44 / 45, -56 / 15, 32 / 9, 0, 0, 0],
        [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0],
        [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0],
        [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
    ]
)
_E = np.array([71 / 57600, 0.0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40])

RightHandSide = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]
"""``f(t, y, idx) -> dy/dt`` evaluated for the problems listed in ``idx``."""


@dataclass
class BatchSolution:
    """Result of :func:`integrate_batch`."""

    y_end: np.ndarray
    """State at ``t_end`` for each problem."""
    trajectories: list[tuple[np.ndarray, np.ndarray]]
    """Accepted ``(t, y)`` points of each problem (only if requested)."""
    steps: np.ndarray
    """Number of accepted steps per problem."""
    failed: np.ndarray
    """Problems that could not be integrated (step size underflow or NaN)."""

    def interpolate(self, t: np.ndarray) -> np.ndarray:
        """Linear interpolation of every trajectory at the common abscissae ``t`` (n x len(t))."""
        return np.vstack([np.interp(t, tt, yy) for tt, yy in self.trajectories])


def integrate_batch(
    fun: RightHandSide,
    y0: np.ndarray,
    t_end: float = 1.0,
    *,
    rtol: float = 1e-2,
    atol: float = 1e-6,
    first_step: float = 1e-2,
    max_step: float = 0.5,
    max_steps: int = 500,
    max_iterations: int = 5_000,
    keep_trajectories: bool = False,
    trace: list | None = None,
) -> BatchSolution:
    """Integrate ``dy/dt = f(t, y)`` from ``t = 0`` to ``t_end`` for a batch of problems.

    A problem that needs more than ``max_steps`` accepted steps, or whose step size
    underflows, is marked as failed (its ``y_end`` is then meaningless). When
    ``trace`` is a list, a tuple ``(problem, t, h, err, accepted)`` is appended for
    every attempted step (debugging aid).
    """
    y0 = np.asarray(y0, dtype=float)
    n = y0.size
    t = np.zeros(n)
    y = y0.copy()
    h = np.full(n, min(max_step, first_step))
    done = np.zeros(n, dtype=bool)
    failed = np.zeros(n, dtype=bool)
    nofailed = np.ones(n, dtype=bool)
    steps = np.zeros(n, dtype=int)
    threshold = atol / rtol
    eps = np.finfo(float).eps
    traj_t: list[list[float]] = [[0.0] for _ in range(n)] if keep_trajectories else []
    traj_y: list[list[float]] = [[float(v)] for v in y0] if keep_trajectories else []

    f_cur = fun(t, y, np.arange(n))
    stages = np.zeros((7, n))

    for _ in range(max_iterations):
        active = np.flatnonzero(~done & ~failed)
        if active.size == 0:
            break
        ta, ya, ha = t[active], y[active], h[active]
        hmin = np.maximum(16.0 * eps * np.abs(ta), 1e-300)
        ha = np.minimum(max_step, np.maximum(hmin, ha))
        remaining = t_end - ta
        last = 1.1 * ha >= remaining
        ha = np.where(last, remaining, ha)

        stages[0, active] = f_cur[active]
        for s in range(1, 7):
            incr = np.zeros(active.size)
            for k in range(s):
                incr += _A[s, k] * stages[k, active]
            stages[s, active] = fun(ta + _C[s] * ha, ya + ha * incr, active)
        y_new = ya + ha * (
            _A[6, 0] * stages[0, active]
            + _A[6, 2] * stages[2, active]
            + _A[6, 3] * stages[3, active]
            + _A[6, 4] * stages[4, active]
            + _A[6, 5] * stages[5, active]
        )
        err_vec = np.zeros(active.size)
        for k in range(7):
            err_vec += _E[k] * stages[k, active]
        scale = np.maximum(np.maximum(np.abs(ya), np.abs(y_new)), threshold)
        err = ha * np.abs(err_vec) / scale
        bad = ~np.isfinite(err) | ~np.isfinite(y_new)
        err = np.where(bad, np.inf, err)

        accept = err <= rtol
        if trace is not None:
            for p, tp, hp, ep, ap in zip(active, ta, ha, err, accept, strict=True):
                trace.append((int(p), float(tp), float(hp), float(ep), bool(ap)))
        # --- rejected steps
        rej = active[~accept]
        if rej.size:
            e_r = err[~accept]
            h_r = ha[~accept]
            first = nofailed[rej]
            with np.errstate(divide="ignore", invalid="ignore"):
                shrink = np.where(
                    first,
                    np.maximum(0.1, 0.8 * (rtol / e_r) ** 0.2),
                    0.5,
                )
            shrink = np.where(np.isfinite(shrink), shrink, 0.1)
            h_new = np.maximum(hmin[~accept], h_r * shrink)
            underflow = h_r <= hmin[~accept]
            failed[rej[underflow]] = True
            h[rej] = h_new
            nofailed[rej] = False
        # --- accepted steps
        acc = active[accept]
        if acc.size:
            e_a = err[accept]
            h_a = ha[accept]
            grow = np.ones(acc.size)
            nf = nofailed[acc]
            with np.errstate(divide="ignore", invalid="ignore"):
                temp = 1.25 * (e_a / rtol) ** 0.2
            temp = np.where(np.isfinite(temp), temp, 0.0)
            grow = np.where(temp > 0.2, 1.0 / np.maximum(temp, 1e-300), 5.0)
            h[acc] = np.where(nf, h_a * grow, h_a)
            t[acc] = ta[accept] + h_a
            y[acc] = y_new[accept]
            f_cur[acc] = stages[6, acc]
            steps[acc] += 1
            nofailed[acc] = True
            done[acc[last[accept]]] = True
            failed[acc[(steps[acc] >= max_steps) & ~done[acc]]] = True
            if keep_trajectories:
                for i, ti, yi in zip(acc, t[acc], y[acc], strict=True):
                    traj_t[i].append(float(ti))
                    traj_y[i].append(float(yi))
    else:
        failed[~done] = True

    trajectories = (
        [(np.asarray(tt), np.asarray(yy)) for tt, yy in zip(traj_t, traj_y, strict=True)]
        if keep_trajectories
        else []
    )
    return BatchSolution(y_end=y, trajectories=trajectories, steps=steps, failed=failed)
