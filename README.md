# gaslift — two-phase pipe-network simulation and continuous gas-lift optimization

Python reimplementation of the model, the case studies and the optimization problems of

> de Souza, J. N. M., de Medeiros, J. L., Costa, A. L. H., Nunes, G. C. (2010).
> **Modeling, simulation and optimization of continuous gas lift systems for deepwater
> offshore petroleum production.** *Journal of Petroleum Science and Engineering*, 72 (3–4),
> 277–289. [doi:10.1016/j.petrol.2010.03.028](https://doi.org/10.1016/j.petrol.2010.03.028)

The original MATLAB program (`bbnet10`) that produced the published numbers is preserved in
`legacy/matlab/`; the `gaslift` package reproduces its results to machine precision (state
vectors identical to ~10⁻¹⁰) and regenerates every table, figure and digraph of the paper.
This repository is part of a three-repository series with common conventions
([`docs/conventions.md`](docs/conventions.md)).

## What the package does

* **Two-phase flow** in pipes of any inclination with the Beggs & Brill (1973) correlation
  (flow patterns, hold-up, Churchill friction factor, polytropic gas) — `gaslift.beggs_brill`,
  `gaslift.pipe_model`.
* **Pipe networks** as digraphs (N nodes, S lines): nodal mass balances, momentum balance
  integrated along each pipe in both directions, nodal momentum balances and N + S nodal
  specifications, solved by Newton–Raphson from a laminar linear estimate —
  `gaslift.network`, `gaslift.model`, `gaslift.linear_model`, `gaslift.solver`,
  `gaslift.simulation`.
* **Batched Runge–Kutta integrator** (Dormand–Prince 5(4) with the `ode45` step control) that
  integrates all pipes and all Jacobian perturbations at once — `gaslift.integrator`.
* **Gas-lifted wells** as a reusable template (`WellDesign`: wellbore, annulus, gas lines,
  injection point, downhole; 6 nodes and 6 lines per well).
* **Economics and optimization**: oil revenue, fuel gas, capital cost of compressor, turbine
  and gas lines (Couper et al., McCoy) and the three problems of the paper (maximum
  production, maximum profit, optimal design) solved with SQP (SciPy SLSQP) coupled to the
  simulator — `gaslift.economics`, `gaslift.optimization`.
* **Reproduction tooling**: the paper cases as editable JSON (`cases/`), a command line that
  runs the studies, compares them with the MATLAB results and the published tables and exports
  the figures, and executed notebooks.

### Reproduction status

| Paper item | Status |
|---|---|
| Table 1, Fig. 2 (fluids, Beggs & Brill map) | reproduced |
| Tables 2–12 (case geometries) | reproduced with the values actually used by the original program (4 in wellbore, 1710 m annulus, …); the printed geometry differs — see errata A1 |
| Figs. 3–5, 9–10 (topologies, digraphs) | reproduced (schematic layout) |
| Figs. 6–8, 11–12 (pressure contours) | reproduced; state vectors equal to the MATLAB sweeps to ~10⁻¹⁰ |
| Figs. 13–15 (maximum production) | reproduced; optima within 10⁻³ kg/s of the MATLAB ones (at 260 bar SLSQP finds a slightly better point) |
| Table 13, Fig. 16 (maximum profit) | all printed digits reproduced at the published decision variables; re-optimised within 3·10⁻⁴ |
| Tables 14–15 (Case 2 maximum production) | reproduced (rates within 6·10⁻³ kg/s) |
| Tables 16–19 (optimal design) | all printed digits of Tables 18–19 reproduced at the published decision variables; Table 17 within 0.012 in / 5·10⁻³ kg/s |

Details per item, the differences between the paper text and the original program (well
geometry, annulus orientation, discharge-pressure node, cost correlations, equation counts)
and the notes on the MATLAB code are in [`docs/reproduction.md`](docs/reproduction.md); the
formal errata in [`docs/errata.md`](docs/errata.md). The package reproduces the program by
default and offers the paper's variants as options.

## Quick start

Requires Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
make setup        # uv sync --all-extras
make test         # fast tests (~20 s); make test-all adds the slow optimisation regressions
make run          # both paper cases: operating points + the 17 optimisation problems -> results/ (~2 min)
make eval         # compare results/ with the MATLAB results and the published tables -> results/eval/
make figures      # export Figs. 2-16 and Tables 1-19 -> figures/ (~2 min; `gaslift figures --quick` for coarser sweeps)
make notebooks    # regenerate and execute notebooks/ (~4 min)
make lint         # ruff
```

Command line (`uv run gaslift ...`):

```bash
gaslift info                               # paper reference, cases, version
gaslift validate case1                     # or a JSON file: topology, specifications, equation count (--simulate solves it)
gaslift run case2 --out results            # one case; `gaslift run` runs both; --no-optimize for a quick run; --only fig13
gaslift eval --results results --out results/eval
gaslift figures --out figures --quick
```

Python:

```python
from gaslift.cases import load_case
from gaslift.optimization import GasLiftProblem, Objective, Reservoir

case = load_case("case1")                               # cases/case1_single_well.json (Tables 2-6)
case.set_well_rates(q_oil=10.0, q_gas=2.0)
solution = case.simulator().solve()                     # steady state (66 unknowns)
print(solution.wells_table())                           # downhole, injection and platform pressures

problem = GasLiftProblem(case, [Reservoir(pressure=220.0, productivity_index=1.0)],
                         Objective.MAX_PRODUCTION, oil_bounds=(5, 20), gas_bounds=(1, 4))
result = problem.solve(problem.initial_vector(10.0, 2.0))
print(result.summary())                                 # 9.598 kg/s of oil with 1.812 kg/s of gas (Table 13, "Max")
```

## Notebooks

Generated and executed by `scripts/build_notebooks.py` (`make notebooks`); narrative in
Portuguese. Every table and figure cell names the paper item it reproduces and shows the
published values next to the recomputed ones.

| Notebook | Paper items |
|---|---|
| `01_two_phase_flow_model` | Table 1, Fig. 2, momentum balance along one pipe (Eqs. 1–4, 8–10) |
| `02_case1_single_well_simulation` | Tables 2–6, Figs. 3–8, validation against the MATLAB sweep |
| `03_case1_optimization` | Figs. 13–16, Table 13 |
| `04_case2_multiwell_simulation` | Tables 7–12, Figs. 9–12, validation against the MATLAB sweep |
| `05_case2_optimization` | Tables 14–19 |
| `06_paper_vs_code` | effect of the differences between the paper text and the original program |

## Configuring a case

A case is a JSON file (`cases/*.json`, documented in [`cases/README.md`](cases/README.md))
with the network (nodes with coordinates, pipes as 3-D polylines with diameter and roughness
per vertex, wells), the nodal specifications (pressure and/or external liquid and gas flows,
N + S values) and a `study` block (operating points, rate sweeps and optimisation problems
with their reservoirs, prices, bounds and starting points). Copy a file, edit it and run
`gaslift validate my_case.json` / `gaslift run my_case.json`.

The same objects are available in Python: `Network.add_node`, `Network.add_pipe`,
`Network.add_well(WellDesign(...))`, `Specifications`, `CaseSetup`, `save_case`. The
builders `case1_single_well()` and `case2_multi_well()` accept keywords for the paper's
geometry variants (`wellbore_diameter_in`, `annulus_length`, `annulus_orientation`,
`annulus_diameter_model`, `gas_line_diameter_in`, …), and `CostModel` /
`CompressorParameters` expose the printed forms of the cost correlations.

## Layout

```
README.md  CITATION.cff  LICENSE  CLAUDE.md  Makefile  pyproject.toml  uv.lock  .python-version  .gitignore
src/gaslift/      the package (see docs/architecture.md for the module map)
cases/            case1_single_well.json, case2_multi_well.json (+ README.md)
data/reference/   MATLAB results (sweeps, optima) and the published Tables 13-19
notebooks/        executed notebooks 01-06, generated by scripts/build_notebooks.py
scripts/          build_notebooks.py, export_cases.py, export_matlab_reference.py
tests/            pytest suite (fast by default; `slow` marker for the optimisation regressions)
docs/             formulation.md, architecture.md, reproduction.md, errata.md, conventions.md, pt/
legacy/matlab/    original MATLAB program bbnet10 and the case/optimisation scripts (read-only)
results/          outputs of `make run` and `make eval` (not versioned)
figures/          figures and tables exported by `make figures` (not versioned)
```

The paper PDF and the raw material (manuscripts, `.mat` result files) are not kept in the
repository; they live in `~/code/personal/paper-raw-material/paper-optimization-of-gas-lift-systems/`.

## Documents

* [`docs/formulation.md`](docs/formulation.md) — equation ↔ code map (Sections 2–6, Eqs. 1–35).
* [`docs/architecture.md`](docs/architecture.md) — package layout and design decisions.
* [`docs/reproduction.md`](docs/reproduction.md) — status per table/figure, paper ↔ program
  differences, notes on the MATLAB code and on `data/reference`.
* [`docs/errata.md`](docs/errata.md) — formal errata of the paper (location, published text,
  correction, impact).
* [`docs/conventions.md`](docs/conventions.md) — conventions of the series.
* [`docs/pt/paper_vs_code.md`](docs/pt/paper_vs_code.md) — working notes in Portuguese.

## Resumo em português

`gaslift` reimplementa em Python o modelo de escoamento bifásico em redes de dutos
(Beggs & Brill, balanços nodais e de momento integrados, Newton–Raphson) e os problemas de
otimização de gas lift contínuo (produção máxima, lucro máximo e projeto ótimo) do artigo de
de Souza, de Medeiros, Costa e Nunes (*J. Pet. Sci. Eng.* 72, 2010). O código MATLAB original
está em `legacy/matlab/` e é reproduzido com precisão de máquina; os casos do artigo estão em
`cases/*.json`, editáveis, e a linha de comando `gaslift run | eval | figures` regenera os
resultados, compara-os com as tabelas publicadas e exporta as figuras. Os notebooks
(`notebooks/`, em português) reproduzem as Tabelas 1–19 e as Figs. 2–16 e discutem as
divergências entre o texto do artigo e o programa original (`docs/pt/paper_vs_code.md`,
errata em `docs/errata.md`).

## License

MIT for the code (see `LICENSE`). The paper and its figures belong to the publisher and are
not distributed with this repository.
