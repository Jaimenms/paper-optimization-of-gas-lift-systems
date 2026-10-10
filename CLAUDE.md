# gaslift — notes for coding agents

Python port of the MATLAB model behind de Souza, de Medeiros, Costa & Nunes (2010),
*J. Pet. Sci. Eng.* 72, 277–289. Goal: reproduce every table and figure of the paper.
This repository follows the conventions of a three-repository series: read
`docs/conventions.md` first (layout, Makefile targets, CLI verbs, README and docs set).

## Commands
- `make setup` (`uv sync --all-extras`) — environment (Python 3.12).
- `make test` (`uv run pytest`, ~20 s) · `make test-all` (adds the `slow` optimisation tests).
- `make lint` / `make format` — ruff, line length 100.
- `make run` → `results/` (both cases, ~2 min; `uv run gaslift run case1 --no-optimize` for seconds).
- `make eval` → `results/eval/report.md` + CSV (comparison with MATLAB and the published tables).
- `make figures` → `figures/` (Figs. 2–16, Tables 1–19; `--quick` for coarser sweeps).
- `make notebooks` — regenerate and execute `notebooks/` (~4 min); never edit `.ipynb` by hand,
  edit `scripts/build_notebooks.py`.
- `uv run python scripts/export_cases.py` — rewrite `cases/*.json` after changing
  `gaslift.cases` or `gaslift.study.default_study` (a test checks they agree).
- `uv run --with scipy python scripts/export_matlab_reference.py` — rebuild
  `data/reference/*.json` from the `.mat` files kept outside the repository
  (`~/code/personal/paper-raw-material/paper-optimization-of-gas-lift-systems/raw/`).

## Layout
- `src/gaslift/`: `beggs_brill` (closure), `pipe_model` + `integrator` (batched
  Dormand–Prince, ode45-compatible), `model` (residuals/Jacobian), `linear_model` (initial
  estimate), `solver` (Newton), `simulation` (facade), `network` (digraph, wells,
  specifications), `cases` (paper cases, JSON I/O), `economics`, `optimization` (SLSQP),
  `study` (paper studies as data + runners), `evaluation` (eval), `figures` (figures),
  `tables`, `sweeps`, `plotting`, `cli`. Map in `docs/architecture.md`.
- `cases/` JSON cases (network + specifications + study); `data/reference/` MATLAB results
  and published tables; `legacy/matlab/bbnet10/` original MATLAB source (read-only).
- `docs/`: `formulation.md` (equations ↔ code), `reproduction.md` (status and paper ↔
  program differences), `errata.md`, `pt/paper_vs_code.md` (Portuguese notes).
- Paper PDF, figures and raw material are **not** in the repository (nor in its history).

## Conventions
- SI units internally (Pa, m, kg/s); bar/inch only at the API edges and in tables.
- State vector `[qL, qG, PI, PO, PW, wL, wG]`, pressures scaled by 1e7 Pa.
- Defaults reproduce the original program, not the paper text; add "corrected" behaviours as
  options and document them in `docs/reproduction.md` (and `docs/errata.md` if a paper error).
- Keep tests green and ruff clean; do not commit `results/` or `figures/`.
