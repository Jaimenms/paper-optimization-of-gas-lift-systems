# gaslift — notes for coding agents

Python port of the MATLAB model behind de Souza et al. (2010), *J. Pet. Sci. Eng.*
Goal: reproduce every table and figure of the paper in `notebooks/`.

## Commands
- `uv sync --all-extras` — create the environment (Python 3.12).
- `uv run pytest` — fast tests (~10 s). `uv run pytest -m slow` — optimisation tests.
- `uv run ruff format src tests scripts && uv run ruff check src tests scripts`.
- `uv run python scripts/build_notebooks.py --execute` — regenerate and execute notebooks.
- `uv run --with scipy python scripts/export_matlab_reference.py` — rebuild
  `data/reference/*.json` from the raw `.mat` files (needs `docs/raw`, not versioned).

## Layout
- `src/gaslift/` — package: `beggs_brill` (closure), `pipe_model` + `integrator`
  (batched Dormand–Prince, ode45-compatible), `model` (residuals/Jacobian), `linear_model`
  (initial estimate), `solver` (Newton), `simulation` (facade), `cases` (paper cases),
  `economics`, `optimization` (SLSQP), `sweeps`, `plotting`.
- `legacy/matlab/bbnet10/` — original MATLAB source (read-only reference).
- `data/reference/` — MATLAB results used by tests; `docs/paper_vs_code.md` lists where
  the paper text and the code disagree. Defaults reproduce the code, not the text.

## Conventions
- SI units internally (Pa, m, kg/s); bar/inch only at the API edges and in tables.
- State vector `[qL, qG, PI, PO, PW, wL, wG]`, pressures scaled by 1e7 Pa.
- Keep the faithful behaviour as default; add new "corrected" behaviours as options.
