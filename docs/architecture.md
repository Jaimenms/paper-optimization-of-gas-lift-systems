# Architecture

`gaslift` is a Python port of the MATLAB program `bbnet10` behind de Souza et al. (2010),
*J. Pet. Sci. Eng.* 72, 277–289. The goal is to reproduce every table and figure of the paper
with a package that is also usable for other networks. The repository follows the layout and
commands of the series ([`conventions.md`](conventions.md)).

## Package layout (`src/gaslift/`)

| Module | Content |
|---|---|
| `units.py` | constants and conversions (g, R, inch, bar, barrel, seconds per year, kW → hp) |
| `fluids.py` | `FluidProperties`: incompressible liquid, polytropic gas (Table 1 defaults) |
| `geometry.py` | `PipeGeometry`: 3-D polylines with per-vertex diameter and roughness, knots with short property transitions (the `tinyl` trick of the original code); `annulus_equivalent_diameter` |
| `beggs_brill.py` | vectorised Beggs & Brill (1973) closure: flow map with logistic "soft if" tests (as in the MATLAB code, so the model is continuous), hold-up, Churchill friction factor, two-phase multiplier |
| `pipe_model.py` | pressure gradient along a pipe (Eqs. 1–4), batched forward/backward integration (Eqs. 8–10), mid-point (non-rigorous) variant, profiles along the pipes |
| `integrator.py` | batched Dormand–Prince 5(4) with the `ode45` step-size control; every pipe and every Jacobian perturbation is an independent problem with its own adaptive step |
| `network.py` | `Network` (digraph of named nodes and pipes), `WellDesign`/`Well` templates (6 nodes + 6 pipes per gas-lifted well, Fig. 5), `StateLayout` (index bookkeeping of `[qL, qG, PI, PO, PW, wL, wG]`), `Specifications` (matrix E, vector e), serialisation |
| `model.py` | `NetworkModel`: residuals (Eqs. 6–14) and Jacobian; `ModelOptions` (pressure scale 10⁷ Pa, tolerances, finite-difference steps) |
| `linear_model.py` | laminar linear model for the initial estimate (Eqs. 15–23) |
| `solver.py` | Newton–Raphson with bound clipping and status reporting |
| `simulation.py` | `Simulator` facade: initial estimate, rigorous/non-rigorous fall-backs, warm start, profiles; `simulate()` one-shot helper |
| `results.py` | `NetworkSolution` (unit-aware accessors, pandas tables), `PipeProfiles` |
| `cases.py` | `CaseSetup` (network + specifications + wells + paper numbering + study block), builders `case1_single_well` / `case2_multi_well` (faithful to `case_03.m` / `case_04.m`, with keyword switches for the paper's variants), JSON `load_case` / `save_case` |
| `sweeps.py` | rate sweeps over injected gas × produced oil (Figs. 6–8, 11–12) with serpentine warm start |
| `economics.py` | revenue, fuel gas, compressor/turbine/pipeline capital costs (Eqs. 28–35), `CostModel` and `CompressorParameters` with the program's defaults and switches for the printed forms |
| `optimization.py` | `GasLiftProblem` (Eq. 24) solved with SciPy SLSQP; reservoir inflow as an equality constraint; evaluation cache and history |
| `study.py` | the paper studies as data (`default_study`): operating points, sweeps and optimisation problems; runners that write JSON results (`gaslift run`) |
| `tables.py` | pandas tables in the layout of Tables 1–19, with the published values side by side |
| `evaluation.py` | comparison of `results/` with the MATLAB results and the published tables; markdown report + CSV (`gaslift eval`) |
| `figures.py` | export of Figs. 2–16 and Tables 1–19 to `figures/` (`gaslift figures`) |
| `plotting.py` | matplotlib helpers: contour plots, 3-D topology, schematic digraph, flow map, pressure paths, optimisation history |
| `cli.py` | `gaslift info | validate | run | eval | figures` (argparse) |

Data flow: `cases/*.json` → `load_case` → `CaseSetup` → `Simulator` (→ `NetworkModel` →
`pipe_model` → `beggs_brill`) → `NetworkSolution` → `evaluate_economics` → `GasLiftProblem`.
`study.run_study` drives the whole chain for a case and writes `results/<case>/*.json`;
`evaluation.evaluate` and `figures.export_figures` only read those files and the
reference data.

## Design decisions

* **Faithful by default, corrected as options.** The package reproduces the MATLAB program
  that produced the published numbers (state vectors agree to ~10⁻¹⁰), which is what makes
  the port verifiable. Every place where the paper text differs from the program is exposed as
  a parameter (`WellDesign.annulus_orientation`, `annulus_diameter_model`,
  `CostModel.discharge_pressure_node`, `pipeline_cost_formulation`,
  `compressor_cost_power_unit`, `CompressorParameters.polytropic_exponent`, the geometry
  keywords of `case1_single_well`/`case2_multi_well`). See [`reproduction.md`](reproduction.md).
* **Batched integration.** The MATLAB code integrated each pipe with `ode45` inside a numerical
  Jacobian loop. Here all pipes, both integration directions and all four perturbed copies of
  the pipe variables are integrated together by one vectorised Dormand–Prince solver with
  per-problem adaptive steps that mirror `ode45`'s control. This keeps the results identical
  while making a 192-unknown simulation take ~50 ms.
* **Two finite-difference steps for the Jacobian.** A coarse relative perturbation (10⁻³, the
  value of the original code) is robust far from the solution; a fine one (10⁻⁵) is used once
  the residual is small so that the perturbed integrations follow the same step sequence as
  the base one, restoring quadratic convergence.
* **Scaled state.** Pressures are divided by 10⁷ Pa so that all unknowns are O(1), as in the
  program (`PR`); the reference state vectors in `data/reference` are stored in this scaling.
* **Specifications are flow rates.** The paper specifies produced oil and injected gas and
  computes pressures; the reservoir inflow relation enters only as an optimisation constraint
  (Eq. 24). Consequently the gas side does not feed back on the production side.
* **Cases are data.** The paper cases live in `cases/*.json` (network, specifications and the
  study to run); the Python builders generate those files and the tests check that the two
  agree. New networks can be written in JSON or assembled with `Network.add_node` /
  `add_pipe` / `add_well` and saved with `save_case`.
* **Notebooks are generated.** `scripts/build_notebooks.py` writes and executes the
  notebooks so that their content is reviewable as code and always matches the package.
* **Reference data in the repository, raw material outside.** `data/reference/*.json` holds
  the MATLAB results exported from the original `.mat` files
  (`scripts/export_matlab_reference.py`); the `.mat` files, the manuscript and the paper PDF
  stay outside the repository (`~/code/personal/paper-raw-material/`).

## Conventions

* SI units internally; `bar`, `in`, `10⁶ USD` only in tables, CLI output and the economic
  model's public parameters.
* ruff (line length 100), dataclasses for the data model, pytest with the `slow` marker for
  the optimisation regressions.
* Keep the faithful behaviour as the default; add a new "corrected" behaviour as an option and
  document it in `reproduction.md` (and in `errata.md` if it is an error of the paper).
