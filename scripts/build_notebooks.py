"""Generate the Jupyter notebooks of ``notebooks/`` from Python sources.

The notebooks are kept as generated artefacts so that their content is
reviewable as code. Run::

    uv run python scripts/build_notebooks.py            # write the .ipynb files
    uv run python scripts/build_notebooks.py --execute  # also execute them in place (~15 min)
    uv run python scripts/build_notebooks.py --execute --only 02 05

Every table or figure cell names the item of the paper it reproduces; the original
figures of the paper are not embedded (the paper is not part of the repository), so
each reproduced figure carries a "compare with Fig. N of the paper" note. The
narrative is in Portuguese.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

ROOT = Path(__file__).resolve().parents[1]
NB_DIR = ROOT / "notebooks"

NAMES = {
    "01": "01_two_phase_flow_model.ipynb",
    "02": "02_case1_single_well_simulation.ipynb",
    "03": "03_case1_optimization.ipynb",
    "04": "04_case2_multiwell_simulation.ipynb",
    "05": "05_case2_optimization.ipynb",
    "06": "06_paper_vs_code.ipynb",
}

REFERENCE = """Referência: de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010).
*Modeling, simulation and optimization of continuous gas lift systems for deepwater
offshore petroleum production*. Journal of Petroleum Science and Engineering 72 (3–4),
277–289, doi:10.1016/j.petrol.2010.03.028."""

SETUP = """import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / "pyproject.toml").exists())
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
REF = ROOT / "data" / "reference"
published = json.load(open(REF / "published_tables.json"))   # Tables 13-19 as printed in the paper
pd.set_option("display.precision", 4)
pd.set_option("display.width", 140)
plt.rcParams.update({"figure.dpi": 100, "figure.figsize": (7, 5), "axes.grid": True, "grid.alpha": 0.3})
"""


def md(text: str):
    return new_markdown_cell(text.strip("\n"))


def code(text: str):
    return new_code_cell(text.strip("\n"))


def compare(fig: str, extra: str = "") -> str:
    """Markdown note placed after a reproduced figure."""
    text = (
        f"> Compare com a {fig} do artigo (a figura original não é distribuída com o repositório)."
    )
    if extra:
        text += " " + extra
    return text


def write(name: str, cells: list) -> Path:
    nb = new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    }
    nb.metadata["language_info"] = {"name": "python"}
    path = NB_DIR / name
    nbformat.write(nb, path)
    return path


# =============================================================================
# 01 - two-phase flow model
# =============================================================================
def nb01() -> Path:
    cells = [
        md(
            f"""
# 01 · Modelo de escoamento bifásico: Beggs & Brill e balanço de momento

Este notebook apresenta os blocos de base do simulador: as propriedades dos fluidos
(Tabela 1 do artigo), o mapa de padrões de escoamento da correlação de Beggs & Brill
(1973) (Fig. 2) e a integração do balanço de momento ao longo de um duto (Eqs. 1–4 e 8–10).

{REFERENCE}

O mapa equação ↔ código está em `docs/formulation.md`.
"""
        ),
        code(SETUP),
        md("## Tabela 1 · Propriedades dos fluidos (reproduz a Tabela 1 do artigo)"),
        code(
            """
from gaslift import paper_fluid
from gaslift.tables import fluid_table

fluid = paper_fluid()
fluid_table(fluid)
"""
        ),
        md(
            """
A massa específica do gás segue a relação politrópica $\\rho_G = \\rho_{G0}\\,(P/P_0)^{1/\\gamma}$;
o líquido é incompressível.
"""
        ),
        code(
            """
p = np.linspace(1e5, 400e5, 200)
plt.plot(p / 1e5, fluid.gas_density(p))
plt.xlabel("Pressure (bar)"); plt.ylabel("Gas density (kg/m³)"); plt.title("Polytropic gas density (Table 1)");
"""
        ),
        md(
            """
## Fig. 2 · Mapa de padrões de escoamento horizontal (Beggs & Brill)

Os testes lógicos da correlação original são substituídos por funções logísticas muito
abruptas (fator $10^6$), como no código MATLAB que gerou o artigo, para que o modelo seja
contínuo. O mapa abaixo é construído exatamente com essa função (`gaslift.beggs_brill.flow_pattern`).
"""
        ),
        code(
            """
from gaslift.plotting import plot_flow_pattern_map

fig, ax = plt.subplots(figsize=(7, 5.5))
plot_flow_pattern_map(ax)
ax.set_title("Fig. 2 (reproduced) — Beggs & Brill flow map");
"""
        ),
        md(compare("Fig. 2")),
        md(
            """
## Hold-up e gradiente de pressão em um duto vertical

Para um duto de 4 in com 850 kg/m³ de óleo, o gradiente de pressão passa por um mínimo
quando a vazão de gás aumenta: o gás reduz a parcela hidrostática mas aumenta o atrito.
Esse é o mecanismo que produz a vazão ótima de injeção do gas lift.
"""
        ),
        code(
            """
from gaslift.beggs_brill import closure_relations, PATTERN_NAMES
from gaslift.units import G, INCH

D = 4 * INCH
area = np.pi * D**2 / 4
pressure = 150e5
q_oil = 10.0  # kg/s
q_gas = np.linspace(0.02, 6.0, 300)
v_sl = np.full_like(q_gas, q_oil / (area * fluid.liquid_density))
v_sg = q_gas / (area * fluid.gas_density(pressure))
res_v = closure_relations(pressure, v_sl, v_sg, D, 100e-6, np.pi / 2, fluid)
res_h = closure_relations(pressure, v_sl, v_sg, D, 100e-6, 0.0, fluid)
grad_v = (res_v.mixture_density * G + 4 * res_v.shear_stress / D) / 1e5 * 100  # bar per 100 m
grad_h = (4 * res_h.shear_stress / D) / 1e5 * 100

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
axes[0].plot(q_gas, res_v.liquid_holdup, label="vertical"); axes[0].plot(q_gas, res_h.liquid_holdup, label="horizontal")
axes[0].plot(q_gas, res_v.no_slip_holdup, "k--", label="no-slip λ")
axes[0].set_xlabel("Gas mass rate (kg/s)"); axes[0].set_ylabel("Liquid hold-up"); axes[0].legend()
axes[1].plot(q_gas, grad_v); axes[1].set_xlabel("Gas mass rate (kg/s)"); axes[1].set_ylabel("Pressure gradient, vertical (bar/100 m)")
axes[1].axvline(q_gas[np.argmin(grad_v)], color="r", ls=":", label=f"min at {q_gas[np.argmin(grad_v)]:.2f} kg/s"); axes[1].legend()
axes[2].plot(q_gas, res_v.pattern.index); axes[2].set_yticks([1, 2, 3, 4]); axes[2].set_yticklabels([PATTERN_NAMES[i] for i in (1, 2, 3, 4)])
axes[2].set_xlabel("Gas mass rate (kg/s)"); axes[2].set_title("Flow pattern, vertical pipe")
plt.tight_layout()
"""
        ),
        md(
            """
## Integração do balanço de momento em um duto (Eqs. 8–10)

Um poço vertical de 2000 m e 4 in produzindo 10 kg/s de óleo com 2 kg/s de gás, com
15 bar no topo. O modelo rigoroso integra $dP/dx$ nas duas direções (Dormand–Prince 5(4)
com o mesmo controle de passo do `ode45`) e usa a média das duas estimativas; o modelo
não rigoroso avalia o gradiente no ponto médio com a pressão média.
"""
        ),
        code(
            """
from gaslift import Network, NodeRole, Specifications, Simulator, SimulationOptions, ModelOptions

net = Network(fluid, "single vertical pipe")
net.add_node("bottom", 0, 0, -2000)
net.add_node("top", 0, 0, 0, NodeRole.OUTLET)
net.add_pipe("wellbore", "bottom", "top", diameter=4 * INCH, roughness=100e-6)
specs = Specifications(net).set_liquid_source("bottom", 10.0, 2.0).set_pressure("top", 15e5)

rigorous = Simulator(net, specs).solve()
nonrig = Simulator(net, specs, SimulationOptions(model=ModelOptions(rigorous=False), compute_profiles=False)).solve()
print(rigorous.status)
print(nonrig.status)
prof = rigorous.profile("wellbore")
plt.plot(prof["pressure_bar"], -2000 + prof["axial_position_m"], label="rigorous (integrated)")
plt.axhline(-2000, color="k", lw=0.5)
plt.plot([nonrig.pressure_bar("bottom"), nonrig.pressure_bar("top")], [-2000, 0], "o--", label="non-rigorous (mid-point)")
plt.xlabel("Pressure (bar)"); plt.ylabel("Elevation (m)"); plt.legend(); plt.title("Pressure along a gas-lifted wellbore")
print(f"bottom pressure: rigorous {rigorous.pressure_bar('bottom'):.2f} bar, non-rigorous {nonrig.pressure_bar('bottom'):.2f} bar")
"""
        ),
        code(
            """
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
axes[0].plot(prof["axial_position_m"], prof["liquid_holdup"], label="H_L"); axes[0].plot(prof["axial_position_m"], prof["no_slip_holdup"], "--", label="λ_L")
axes[0].set_xlabel("Distance from bottom (m)"); axes[0].set_ylabel("Hold-up"); axes[0].legend()
axes[1].plot(prof["axial_position_m"], prof["gas_density_kg_m3"]); axes[1].set_xlabel("Distance from bottom (m)"); axes[1].set_ylabel("Gas density (kg/m³)")
plt.tight_layout()
"""
        ),
    ]
    return write(NAMES["01"], cells)


# =============================================================================
# 02 - case 1 simulation
# =============================================================================
def nb02() -> Path:
    cells = [
        md(
            f"""
# 02 · Caso 1 — poço offshore único: rede, digrafo e simulação

Reprodução da Seção 4.1 do artigo: Tabelas 2–6, Figs. 3–5 (topologia e digrafo) e
Figs. 6–8 (curvas de nível da pressão de fundo e da pressão de injeção de gás).

A rede é lida do arquivo `cases/case1_single_well.json` (o mesmo usado pela linha de
comando `gaslift`); ela também pode ser construída com `gaslift.case1_single_well()`:
nós e linhas da rede de coleta mais um modelo de poço (`WellDesign`) que se expande em
6 nós e 6 linhas (Fig. 5).

{REFERENCE}
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import SimulationOptions, ModelOptions, IntegrationOptions
from gaslift.cases import load_case
from gaslift.plotting import plot_case_topology, plot_digraph, plot_rate_contours, plot_pressure_path
from gaslift.sweeps import sweep_rates
from gaslift.tables import pipes_table, nodes_table, profile_table, specifications_table

case = load_case("case1")
net = case.network
print(net)
print(f"Unknowns / equations: 3N + 4S = {net.layout.size}  (N = {net.n_nodes}, S = {net.n_pipes});"
      f" specifications: {case.specifications.count}")
"""
        ),
        md("## Tabela 2 · Linhas da rede (reproduz a Tabela 2 do artigo)"),
        code("pipes_table(case)"),
        md(
            """
> Os valores seguem os scripts MATLAB que geraram os resultados publicados (`case_03.m`):
> coluna de produção de 4 in, anular de 1710 m, trecho inferior de 338 m e rugosidade nula
> nas linhas de gás. A Tabela 2 impressa no artigo lista 5 in, 1500 m, 300 m e 100 µm —
> ver `docs/reproduction.md` (resumo em português em `docs/pt/paper_vs_code.md`).
"""
        ),
        md("## Tabelas 3 e 4 · Perfis das linhas 2 e 8"),
        code(
            """
display(profile_table(case, "production_line_2").rename_axis("Line 2 (Table 3)"))
display(profile_table(case, "well_1:gas_line_1").rename_axis("Line 8 (Table 4)"))
"""
        ),
        md("## Tabelas 5 e 6 · Nós e especificações"),
        code(
            'display(nodes_table(case).rename_axis("Node (Table 5)"))\n'
            'display(specifications_table(case).rename_axis("Node (Table 6)"))'
        ),
        md("## Fig. 3 · Projeção isométrica"),
        code(
            """
fig = plot_case_topology(net)
fig.suptitle("Fig. 3 (reproduced) — Case 1, simplified isometric projection", y=1.02);
"""
        ),
        md(compare("Fig. 3")),
        md("## Figs. 4 e 5 · Digrafo da rede e do poço"),
        code(
            """
fig, ax = plt.subplots(figsize=(12, 7))
plot_digraph(net, ax, node_numbers=case.paper_node_numbers, pipe_numbers=case.paper_pipe_numbers)
ax.set_title("Figs. 4 and 5 (reproduced) — Case 1 digraph (node numbers as in the paper; well nodes highlighted)");
"""
        ),
        md(compare("Fig. 4 e a Fig. 5")),
        code(
            """
M = net.incidence_matrix()
pd.DataFrame(M.astype(int), index=[f"node {i+1}" for i in range(net.n_nodes)], columns=[f"line {j+1}" for j in range(net.n_pipes)])
"""
        ),
        md(
            """
## Uma simulação: 10 kg/s de óleo e 2 kg/s de gás

O estado estacionário é obtido por Newton–Raphson a partir da estimativa do modelo linear
laminar (Seção 3.4). Abaixo, as tabelas de nós, linhas e poço, e os perfis de pressão ao
longo do caminho de produção (fundo do poço → riser) e do caminho do gás (plataforma →
ponto de injeção).
"""
        ),
        code(
            """
sim = case.simulator()
sol = sim.solve()
print(sol.status)
display(sol.wells_table().T)
display(sol.nodes_table())
"""
        ),
        code("sol.pipes_table()"),
        code(
            """
fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
plot_pressure_path(sol, ["well_1:wellbore_lower", "well_1:wellbore_upper", "well_1:flowline", "production_line_1", "production_line_2", "riser"], axes[0])
axes[0].set_title("Production path: downhole → riser outlet")
plot_pressure_path(sol, ["well_1:gas_line_2", "well_1:gas_line_1", "well_1:annulus"], axes[1])
axes[1].set_title("Gas path: platform → gas line → annulus")
plt.tight_layout()
"""
        ),
        md(
            """
## Figs. 6 e 7 · Varredura de vazões (linhas de gás de 2 in)

O artigo varre a vazão de gás injetado (0,5 a 4 kg/s) e a vazão de óleo (0,5 a 20 kg/s).
Os resultados originais desta varredura foram gerados com tolerância relativa `1e-3` na
integração (`results_Subsea_1_2.mat`), por isso usamos o mesmo valor aqui (a grade e a
tolerância vêm do bloco `study.sweeps` do arquivo do caso). A varredura percorre a grade
em serpentina com partida quente, como o script original.
"""
        ),
        code(
            """
sweep_spec = next(s for s in case.study["sweeps"] if s["name"] == "fig06_fig07_2in")
q_gas, q_oil = np.array(sweep_spec["q_gas"]), np.array(sweep_spec["q_oil"])
opts = SimulationOptions(model=ModelOptions(integration=IntegrationOptions(rtol=sweep_spec["rtol"])), compute_profiles=False)
t0 = time.perf_counter()
sweep = sweep_rates(case.copy(), q_gas, q_oil, opts)
print(f"{sweep.converged.sum()} of {sweep.converged.size} points converged in {time.perf_counter() - t0:.1f} s")
"""
        ),
        code(
            """
fig, ax = plt.subplots(figsize=(7.5, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), ax, np.arange(140, 320, 20), title="Fig. 6 (reproduced) — Downhole pressure, node 7 (bar)");
"""
        ),
        md(compare("Fig. 6")),
        md(
            """
A "pressão de entrada do gas lift" dos gráficos do artigo é lida no nó 9 (interconexão das
linhas de gás na base do riser, −1495 m), que é o nó apontado pelo índice `ix_Pgaslift` do
código original — e não no nó 10 (plataforma), como descreve o texto. A diferença entre os
dois nós é a coluna hidrostática de gás da linha vertical de 1515 m.
"""
        ),
        code(
            """
fig, ax = plt.subplots(figsize=(7.5, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.gas_junction(0), ax, np.arange(200, 560, 50), title="Fig. 7 (reproduced) — Gas lift inlet pressure, node 9 (bar), 2 in gas line");
"""
        ),
        md(compare("Fig. 7")),
        code(
            """
fig, ax = plt.subplots(figsize=(7, 6))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.gas_inlet(0), ax, np.arange(150, 560, 50), title="Platform gas inlet pressure, node 10 (bar), 2 in gas line (not in the paper)");
"""
        ),
        md("## Fig. 8 · Linhas de gás de 3 in"),
        code(
            """
from gaslift.units import INCH
case3 = case.copy()
case3.set_gas_line_diameter("well_1", 3.0 * INCH)
sweep3 = sweep_rates(case3, q_gas, q_oil, opts)
fig, ax = plt.subplots(figsize=(7.5, 6.5))
plot_rate_contours(sweep3.q_gas, sweep3.q_oil, sweep3.gas_junction(0), ax, np.arange(180, 320, 20), title="Fig. 8 (reproduced) — Gas lift inlet pressure, node 9 (bar), 3 in gas line");
"""
        ),
        md(compare("Fig. 8")),
        md(
            """
## Leitura da Fig. 6 (texto da Seção 4.1 do artigo)

> "considering a downhole pressure of 220 bar, to produce 10.0 kg/s of petroleum it is
> necessary to inject 1.0 kg/s of gas or a higher value of gas flow rate, equal to 2.7 kg/s
> [...] at downhole pressure of 220 bar, the maximum oil production is 10.7 kg/s with gas
> injection of 1.8 kg/s."
"""
        ),
        code(
            """
from scipy.interpolate import RegularGridInterpolator

text = published["text"]["section_4.1_fig6"]
interp = RegularGridInterpolator((sweep.q_oil, sweep.q_gas), sweep.downhole(0), bounds_error=False)
qg_fine = np.linspace(0.5, 4.0, 701)
p_at_10 = interp(np.column_stack([np.full_like(qg_fine, 10.0), qg_fine]))
roots = qg_fine[np.flatnonzero(np.diff(np.sign(p_at_10 - 220.0)))]
print("gas rates giving 220 bar at 10 kg/s of oil:", np.round(roots, 2), "kg/s  (paper:", text["gas_kg_s"], ")")

qo_fine = np.linspace(0.5, 20.0, 781)
QO, QG = np.meshgrid(qo_fine, qg_fine, indexing="ij")
P = interp(np.column_stack([QO.ravel(), QG.ravel()])).reshape(QO.shape)
feasible = np.where(P <= 220.0, QO, np.nan)
i, j = np.unravel_index(np.nanargmax(feasible), feasible.shape)
print(f"maximum oil at 220 bar: {qo_fine[i]:.2f} kg/s with {qg_fine[j]:.2f} kg/s of gas  (paper: {text['max_oil_kg_s']} kg/s with {text['gas_at_max_oil_kg_s']} kg/s)")
"""
        ),
        md(
            "## Validação contra os resultados MATLAB originais (`data/reference/case1_grid_2in.json`)"
        ),
        code(
            """
ref = json.load(open(REF / "case1_grid_2in.json"))
ref_p = np.array([[c[0] if c else np.nan for c in row] for row in ref["p_downhole_bar"]], dtype=float)
ref_pg = np.array([[c[0] if c else np.nan for c in row] for row in ref["p_gasinlet_bar"]], dtype=float)
diff_p = np.abs(sweep.downhole(0) - ref_p)
diff_pg = np.abs(sweep.gas_junction(0) - ref_pg)
print(f"MATLAB reference: {np.isfinite(ref_p).sum()} converged points; Python: {sweep.converged.sum()}")
print(f"max |ΔP_downhole| = {np.nanmax(diff_p):.2e} bar, max |ΔP_gas| = {np.nanmax(diff_pg):.2e} bar over common points")
pd.DataFrame(diff_p, index=sweep.q_oil, columns=sweep.q_gas).rename_axis(index="q_oil", columns="q_gas").style.format("{:.1e}")
"""
        ),
    ]
    return write(NAMES["02"], cells)


# =============================================================================
# 03 - case 1 optimisation
# =============================================================================
def nb03() -> Path:
    cells = [
        md(
            f"""
# 03 · Caso 1 — otimização: produção máxima e lucro máximo

Reprodução da Seção 6.1 do artigo: Figs. 13–16 e Tabela 13. O problema de otimização
(Eq. 24) acopla o simulador de rede (que especifica vazões) ao modelo de reservatório
$Q_L = PI\\,(P_{{res}} - P_{{downhole}})$ por meio de uma restrição de igualdade, resolvido com
programação quadrática sequencial (SLSQP do SciPy, equivalente ao `fmincon` usado no artigo).

Os problemas são os do bloco `study.optimizations` de `cases/case1_single_well.json`
(os mesmos que `gaslift run case1` resolve); cada resultado é comparado com o ótimo do
MATLAB (`data/reference/matlab_optimization_results.json`) e, quando impresso, com o valor
do artigo (`data/reference/published_tables.json`).

{REFERENCE}
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import tables
from gaslift.cases import load_case
from gaslift.plotting import plot_rate_contours
from gaslift.study import build_problem, initial_vector, load_reference, optimization_record, simulation_options
from gaslift.sweeps import sweep_rates

case = load_case("case1")
matlab = load_reference("matlab_optimization_results.json")
specs = {s["name"]: s for s in case.study["optimizations"]}
records = {}

def optimise(name):
    spec = specs[name]
    prob = build_problem(case, spec)
    t0 = time.perf_counter()
    res = prob.solve(initial_vector(prob, spec, records))
    records[name] = optimization_record(name, spec, res, time.perf_counter() - t0, matlab.get(spec.get("reference")))
    return res

pd.DataFrame([{"problem": s["name"], "paper": ", ".join(s["paper"]), "objective": s["objective"],
               "reservoir": s["reservoirs"][0], "bounds": (s["oil_bounds"], s["gas_bounds"])} for s in specs.values()]).set_index("problem")
"""
        ),
        md(
            """
## 6.1.1 · Produção máxima para diferentes pressões de fundo (Figs. 13 e 14)

Índice de produtividade infinito: a pressão de fundo é igual à pressão do reservatório.
A superfície da produção é muito plana na direção da vazão de gás perto do ótimo, de modo que
otimizadores diferentes podem parar em vazões de gás ligeiramente distintas com a mesma
produção (compare as colunas Python e MATLAB; a 260 bar o SLSQP encontra produção um pouco
maior que a publicada).
"""
        ),
        code(
            """
for p in (180, 200, 220, 240, 260):
    optimise(f"fig13_{p}bar")
fig13 = tables.fig13_table(records)
fig13
"""
        ),
        code(
            """
fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
axes[0].plot(fig13.index, fig13["oil (kg/s)"], "ko-", label="Python (SLSQP)"); axes[0].plot(fig13.index, fig13["oil MATLAB (kg/s)"], "r+", ms=10, label="MATLAB (paper)")
axes[0].set_ylabel("Produced oil (kg/s)"); axes[0].legend(); axes[0].set_title("Fig. 13 (reproduced) — Maximum production operational conditions")
axes[1].plot(fig13.index, fig13["gas (kg/s)"], "ko-"); axes[1].plot(fig13.index, fig13["gas MATLAB (kg/s)"], "r+", ms=10)
axes[1].set_ylabel("Injected gas (kg/s)"); axes[1].set_xlabel("Downhole pressure (bar)")
plt.tight_layout()
"""
        ),
        md(
            compare(
                "Fig. 13",
                "O texto do artigo (Seção 6.1.1) afirma que a vazão ótima de gás varia de 1,55 a 1,85 kg/s.",
            )
        ),
        code(
            """
sweep_spec = next(s for s in case.study["sweeps"] if s["name"] == "fig06_fig07_2in")
sweep = sweep_rates(case.copy(), np.array(sweep_spec["q_gas"]), np.array(sweep_spec["q_oil"]), simulation_options(sweep_spec["rtol"]))
fig, ax = plt.subplots(figsize=(7.5, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), ax, np.arange(140, 320, 20), title="Fig. 14 (reproduced) — Maximum production locus")
ax.plot(fig13["gas (kg/s)"], fig13["oil (kg/s)"], "ro--", label="Optimum locus"); ax.legend(loc="upper right");
"""
        ),
        md(compare("Fig. 14")),
        md("## 6.1.2 · Produção máxima para diferentes índices de produtividade (Fig. 15)"),
        code(
            """
for pi in (0.5, 1.0, 2.0, 5.0, 10.0):
    optimise(f"fig15_pi{pi:g}")
fig15 = tables.fig15_table(records)
fig15
"""
        ),
        code(
            """
fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
for ax, col, ref, lab in zip(axes, ["oil (kg/s)", "gas (kg/s)", "P_downhole (bar)"], ["oil MATLAB (kg/s)", "gas MATLAB (kg/s)", "P_downhole MATLAB (bar)"],
                             ["Produced oil (kg/s)", "Injected gas (kg/s)", "Well pressure (bar)"]):
    ax.plot(fig15.index, fig15[col], "ko-", label="Python"); ax.plot(fig15.index, fig15[ref], "r+", ms=10, label="MATLAB"); ax.set_ylabel(lab)
axes[0].legend(); axes[0].set_title("Fig. 15 (reproduced) — Maximum production vs productivity index (220 bar)"); axes[-1].set_xlabel("Productivity index ((kg/s)/bar)")
plt.tight_layout()
"""
        ),
        md(compare("Fig. 15")),
        md(
            """
## 6.1.3 · Lucro máximo (Tabela 13 e Fig. 16)

Função objetivo $\\Omega_2 = R - AGC$: receita do óleo menos custo do gás combustível da
turbina do compressor (Eqs. 26, 31, 34 e 35). Os cenários variam o preço do óleo com o preço
do gás fixo em 0,1886 USD/Nm³; o cenário "Max" corresponde à produção máxima (custo nulo),
já calculado acima (`fig15_pi1`). A tabela mostra, lado a lado, os valores recalculados e os
impressos na Tabela 13 do artigo.
"""
        ),
        code(
            """
for scenario in ("1", "2", "3"):
    optimise(f"table13_scenario{scenario}")
table13 = tables.table13(records, published)
table13.rename_axis("Table 13 — Case 1 maximum profit results")
"""
        ),
        code(
            """
cols = [c for c in table13.columns if "(paper)" not in c]
pcols = [c + " (paper)" for c in cols]
ratio = table13.loc["Cost ratio (10³ (USD/Nm³)/(USD/bbl))", cols].astype(float) / 1e3
fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
for ax, row in zip(axes, ["Produced oil (kg/s)", "Injected gas (kg/s)", "Profit (10⁶ USD/yr)"]):
    ax.plot(ratio, table13.loc[row, cols].astype(float), "ko-", label="Python (SLSQP)")
    ax.plot(ratio, table13.loc[row, pcols].astype(float), "r+", ms=10, label="paper (Table 13)")
    ax.set_ylabel(row)
axes[0].legend(); axes[0].set_title("Fig. 16 (reproduced) — Maximum profit operational conditions"); axes[2].set_xlabel("Cost ratio ((USD/Nm³)/(USD/barrel))")
plt.tight_layout()
"""
        ),
        md(compare("Fig. 16")),
    ]
    return write(NAMES["03"], cells)


# =============================================================================
# 04 - case 2 simulation
# =============================================================================
def nb04() -> Path:
    cells = [
        md(
            f"""
# 04 · Caso 2 — múltiplos poços offshore: rede e simulação

Reprodução da Seção 4.2 do artigo: Tabelas 7–12, Figs. 9–10 (topologia e digrafo) e
Figs. 11–12 (curvas de nível da pressão de fundo dos poços 1/2 e 3/4). Quatro poços com
gas lift compartilham a linha de produção e o riser: os poços 1 e 2 entram no nó 1 e os
poços 3 e 4 no nó 2. A rede é lida de `cases/case2_multi_well.json`.

{REFERENCE}
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import SimulationOptions, ModelOptions, IntegrationOptions
from gaslift.cases import load_case
from gaslift.plotting import plot_case_topology, plot_digraph, plot_rate_contours
from gaslift.sweeps import sweep_rates
from gaslift.tables import pipes_table, nodes_table, profile_table, specifications_table

case = load_case("case2")
net = case.network
print(net)
print(f"Unknowns / equations: 3N + 4S = {net.layout.size}  (N = {net.n_nodes}, S = {net.n_pipes});"
      f" specifications: {case.specifications.count}")
"""
        ),
        md(
            """
> O artigo afirma que o modelo deste caso tem 232 equações; com N = 28 nós e S = 27 linhas,
> 3N + 4S = 192. O valor 232 (e a decomposição "27 + 27 + 27 + 56 + 95") não corresponde à
> formulação da Seção 3.3 (`docs/errata.md`, item B1).
"""
        ),
        md("## Tabela 7 · Linhas da rede (reproduz a Tabela 7 do artigo)"),
        code("pipes_table(case)"),
        md("## Tabelas 8, 9 e 10 · Perfis das linhas 2, 8/14 e 20/26"),
        code(
            """
display(profile_table(case, "production_line_2").rename_axis("Line 2 (Table 8)"))
display(profile_table(case, "well_1:gas_line_1").rename_axis("Lines 8 and 14 (Table 9)"))
display(profile_table(case, "well_3:gas_line_1").rename_axis("Lines 20 and 26 (Table 10)"))
"""
        ),
        md(
            """
> A Tabela 10 do artigo lista comprimentos acumulados 0, 1005, 2505 e 4005 m para as linhas
> 20 e 26; a geometria real (e a Tabela 16, que usa $L_{G3} = L_{G4} = 4520$ m) corresponde a
> 3005 m + 1515 m (`docs/errata.md`, item B3).
"""
        ),
        md("## Tabelas 11 e 12 · Nós e especificações"),
        code(
            'display(nodes_table(case).rename_axis("Node (Table 11)"))\n'
            'display(specifications_table(case).rename_axis("Node (Table 12)"))'
        ),
        md("## Fig. 9 · Projeção isométrica"),
        code(
            """
fig = plot_case_topology(net)
fig.suptitle("Fig. 9 (reproduced) — Case 2, simplified isometric projection", y=1.02);
"""
        ),
        md(compare("Fig. 9")),
        md("## Fig. 10 · Digrafo"),
        code(
            """
fig, ax = plt.subplots(figsize=(16, 8))
plot_digraph(net, ax, node_numbers=case.paper_node_numbers, pipe_numbers=case.paper_pipe_numbers)
ax.set_title("Fig. 10 (reproduced) — Case 2 digraph (node numbers as in the paper)");
"""
        ),
        md(compare("Fig. 10")),
        md("## Simulação de um ponto de operação: 3 kg/s de óleo e 1 kg/s de gás por poço"),
        code(
            """
sol = case.simulator().solve()
print(sol.status)
display(sol.wells_table())
sol.pipes_table()[["category", "length_m", "liquid_flow_kg_s", "gas_flow_kg_s", "inlet_pressure_bar", "outlet_pressure_bar", "dominant_pattern"]]
"""
        ),
        md(
            """
## Figs. 11 e 12 · Varredura de vazões

Gás injetado de 0,25 a 4 kg/s e óleo de 0,5 a 8 kg/s, iguais em todos os poços
(256 simulações de 192 incógnitas; tolerância de integração `1e-2`, como no original;
grade definida em `study.sweeps` do arquivo do caso).
"""
        ),
        code(
            """
sweep_spec = case.study["sweeps"][0]
q_gas, q_oil = np.array(sweep_spec["q_gas"]), np.array(sweep_spec["q_oil"])
opts = SimulationOptions(model=ModelOptions(integration=IntegrationOptions(rtol=case.study["integration_rtol"])), compute_profiles=False)
t0 = time.perf_counter()
sweep = sweep_rates(case.copy(), q_gas, q_oil, opts)
print(f"{sweep.converged.sum()} of {sweep.converged.size} points converged in {time.perf_counter() - t0:.1f} s")
"""
        ),
        code(
            """
levels = np.arange(140, 380, 20)
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), axes[0], levels, title="Fig. 11 (reproduced) — Pressure at nodes 7 and 13 (bar)")
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(2), axes[1], levels, title="Fig. 12 (reproduced) — Pressure at nodes 19 and 25 (bar)")
plt.tight_layout()
"""
        ),
        md(compare("Fig. 11 e a Fig. 12")),
        md("## Validação contra os resultados MATLAB originais (`data/reference/case2_grid.json`)"),
        code(
            """
ref = json.load(open(REF / "case2_grid.json"))
ref_p = np.array([[[np.nan] * 4 if c is None else c for c in row] for row in ref["p_downhole_bar"]], dtype=float)
diff = np.abs(sweep.downhole_pressure - ref_p)
both = np.isfinite(diff[:, :, 0])
print(f"MATLAB reference: {np.isfinite(ref_p[:, :, 0]).sum()} converged points; Python: {sweep.converged.sum()}; common: {both.sum()}")
print(f"max |ΔP_downhole| over common points: {np.nanmax(diff):.3f} bar; 99th percentile: {np.nanpercentile(diff, 99):.2e} bar")
pd.DataFrame(diff[:, :, 0], index=sweep.q_oil, columns=sweep.q_gas).rename_axis(index="q_oil", columns="q_gas").style.format("{:.1e}")
"""
        ),
    ]
    return write(NAMES["04"], cells)


# =============================================================================
# 05 - case 2 optimisation
# =============================================================================
def nb05() -> Path:
    cells = [
        md(
            f"""
# 05 · Caso 2 — otimização: produção máxima e projeto ótimo do sistema de gas lift

Reprodução da Seção 6.2 do artigo: Tabelas 14 e 15 (produção máxima para reservatórios a
220 e 240 bar) e Tabelas 16–19 (projeto ótimo com custos anualizados de compressor, turbina,
linhas de gás e gás combustível). As variáveis de decisão do projeto são as vazões de óleo e
gás de cada poço e o diâmetro das linhas de gás de cada poço (12 variáveis).

Os problemas são os do bloco `study.optimizations` de `cases/case2_multi_well.json`; as
tabelas mostram os valores recalculados ao lado dos impressos no artigo.

{REFERENCE}
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import tables
from gaslift.cases import load_case
from gaslift.economics import CostModel, CompressorParameters
from gaslift.study import build_problem, initial_vector, load_reference, optimization_record

case = load_case("case2")
matlab = load_reference("matlab_optimization_results.json")
specs = {s["name"]: s for s in case.study["optimizations"]}
records, results = {}, {}
PI = [r["productivity_index"] for r in specs["table14_maxprod_220bar"]["reservoirs"]]

def optimise(name):
    spec = specs[name]
    prob = build_problem(case, spec)
    t0 = time.perf_counter()
    res = prob.solve(initial_vector(prob, spec, records))
    records[name] = optimization_record(name, spec, res, time.perf_counter() - t0, matlab.get(spec.get("reference")))
    results[name] = res
    print(f"{name}: {res.n_evaluations} simulations in {time.perf_counter() - t0:.0f} s; total oil {res.oil_rates.sum():.3f} kg/s, profit {res.economics.profit:.3f} 10⁶ USD/yr")
    return res
"""
        ),
        md("## 6.2.1 · Produção máxima (Tabelas 14 e 15)"),
        code(
            """
optimise("table14_maxprod_220bar")
optimise("table15_maxprod_240bar")
"""
        ),
        code(
            'tables.table14_15(records["table14_maxprod_220bar"], published["table14"], PI).rename_axis("Table 14 — 220 bar")'
        ),
        code(
            'tables.table14_15(records["table15_maxprod_240bar"], published["table15"], PI).rename_axis("Table 15 — 240 bar")'
        ),
        md(
            """
Observações do artigo que podem ser verificadas acima: (1) para poços de mesmo índice de
produtividade e geometria semelhante (1 e 3), o poço mais próximo da plataforma recebe mais
gás e produz mais; (2) para distâncias semelhantes (3 e 4), mais gás vai ao poço de maior
índice; (3) a 240 bar produz-se mais óleo com menos gás injetado.
"""
        ),
        md(
            """
## 6.2.2 · Projeto ótimo (Tabelas 16–19)

Função objetivo $\\Omega_3 = R - AGC - ATCC$. As estimativas iniciais são as soluções de
produção máxima acima, com linhas de gás de 2 in. Limites: diâmetros entre 1 e 4 in.
A Tabela 16 reúne os parâmetros de projeto usados (o coeficiente politrópico do compressor é
1,4, e não 1,25 como impresso — `docs/errata.md`, item A6).
"""
        ),
        code(
            """
comp = CompressorParameters()
cost = CostModel()
pub16 = published["table16"]["rows"]
table16 = pd.DataFrame({"gaslift": {
    "Number of compressor stages (n)": comp.stages, "Overall compressor efficiency (η)": comp.efficiency,
    "Polytropic coefficient used in Eq. 31 (γ)": comp.polytropic_exponent, "Suction temperature (K)": comp.suction_temperature,
    "Suction pressure (bar)": comp.suction_pressure, "Molecular weight (kg/kmol)": comp.molecular_weight,
    "Fuel gas cost (USD/Nm³)": cost.gas_price, "Oil price (USD/barrel)": cost.oil_price,
    "Pipe length LG1, LG2 (m)": case.gas_line_length("well_1"), "Pipe length LG3, LG4 (m)": case.gas_line_length("well_3"),
    "Diameter bounds (in)": str(specs["table17_19_design_220bar"]["diameter_bounds"]), "Productivity indexes": str(PI), "Reservoir pressure scenarios (bar)": "220, 240"},
    "paper (Table 16)": {
    "Number of compressor stages (n)": pub16["compressor_stages"], "Overall compressor efficiency (η)": pub16["compressor_efficiency"],
    "Polytropic coefficient used in Eq. 31 (γ)": pub16["polytropic_coefficient_printed"], "Suction temperature (K)": pub16["suction_temperature_K"],
    "Suction pressure (bar)": pub16["suction_pressure_bar"], "Molecular weight (kg/kmol)": pub16["molecular_weight_kg_kmol"],
    "Fuel gas cost (USD/Nm³)": pub16["fuel_gas_cost_usd_m3"], "Oil price (USD/barrel)": pub16["oil_price_usd_bbl"],
    "Pipe length LG1, LG2 (m)": pub16["gas_line_length_m"][0], "Pipe length LG3, LG4 (m)": pub16["gas_line_length_m"][2],
    "Diameter bounds (in)": str(pub16["diameter_bounds_in"]), "Productivity indexes": str(pub16["productivity_index"]), "Reservoir pressure scenarios (bar)": "220, 240"}})
table16.rename_axis("Table 16 — design parameters")
"""
        ),
        code(
            """
des220 = optimise("table17_19_design_220bar")
des240 = optimise("table17_19_design_240bar")
"""
        ),
        code(
            """
r1, r2 = records["table17_19_design_220bar"], records["table17_19_design_240bar"]
tables.table17(r1, r2, published).rename_axis("Table 17 — decision variables")
"""
        ),
        code(
            'tables.table18(r1, r2, published).rename_axis("Table 18 — pipelines, compressor and turbine")'
        ),
        code('tables.table19(r1, r2, published).rename_axis("Table 19 — optimal design results")'),
        md(
            """
As pequenas diferenças em relação ao MATLAB (terceira casa decimal) decorrem do otimizador:
a superfície do lucro é muito plana em torno do ótimo, e `fmincon` (2009) e SLSQP param em
pontos ligeiramente distintos com valores de função praticamente iguais. Os valores da
Tabela 18–19 avaliados exatamente nas variáveis de decisão publicadas reproduzem todos os
dígitos impressos (ver `gaslift eval`).
"""
        ),
        code(
            """
from gaslift.plotting import plot_optimization_history
fig, axes = plt.subplots(1, 2, figsize=(14, 4))
plot_optimization_history(des220.history, axes[0]); axes[0].set_title("Optimal design, 220 bar")
plot_optimization_history(des240.history, axes[1]); axes[1].set_title("Optimal design, 240 bar")
plt.tight_layout()
"""
        ),
    ]
    return write(NAMES["05"], cells)


# =============================================================================
# 06 - paper vs code
# =============================================================================
def nb06() -> Path:
    cells = [
        md(
            f"""
# 06 · Artigo × código: divergências encontradas na migração e seu efeito nos resultados

Ao portar o código MATLAB (`bbnet10`) que gerou os números publicados, algumas diferenças
entre o texto do artigo e o código foram identificadas (detalhes em `docs/reproduction.md`
e `docs/errata.md`; resumo em português em `docs/pt/paper_vs_code.md`).
O pacote reproduz o código por padrão e expõe interruptores para as formas "corrigidas".
Este notebook quantifica o efeito de cada uma delas.

| # | Item | Código (padrão, reproduz o artigo) | Forma descrita no artigo / física |
|---|------|------|------|
| 1 | Pressão de descarga $P_G$ | nó n+4 (interconexão das linhas de gás, −1495 m) | nó n+5 (plataforma, +20 m) |
| 2 | Orientação do anular | geometria escrita de baixo para cima (termo hidrostático de escoamento ascendente) | gás desce: −90° |
| 3 | Diâmetro equivalente do anular | $(D_o^2 + D_i^2)/(D_o + D_i)$ | $D_o - D_i$ (hidráulico) |
| 4 | Custo das linhas (McCoy) | $D \\cdot a_7$ | $D^{{a_7}}$ (Eq. 33) |
| 5 | Expoente politrópico no compressor | 1,4 | 1,25 (Tabela 16) |
| 6 | Potência nas correlações de custo | hp | kW (Eqs. 30 e 32) |
| 7 | Geometria do poço | 4 in; 1710 + 338 m; rugosidade 0 nas linhas de gás | 5 in; 1500 + 300 m; 100 µm (Tabela 2) |

{REFERENCE}
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import case1_single_well, case2_multi_well, SimulationOptions
from gaslift.economics import CostModel, CompressorParameters, evaluate_economics, pipeline_capital_cost, compressor_capital_cost, turbine_capital_cost
matlab = json.load(open(REF / "matlab_optimization_results.json"))
sim_opts = SimulationOptions(compute_profiles=False)
"""
        ),
        md(
            """
## Itens 1–3 · Lado do gás: onde se lê a pressão de descarga e como o anular é modelado

O lado do gás não realimenta o lado da produção (as especificações são vazões), portanto
a pressão de fundo, a produção máxima e as Tabelas 14–15 não mudam. O que muda é a pressão
de descarga do compressor e tudo que dela deriva (potência, custos, lucro). Avaliamos o ponto
"Max" da Tabela 13 (9,598 kg/s de óleo, 1,812 kg/s de gás).
"""
        ),
        code(
            """
ref = matlab["case1_maxprod_pi1_220bar"]
q_oil, q_gas = float(np.ravel(ref["q_oil"])[0]), float(np.ravel(ref["q_gas"])[0])
variants = {
    "code (paper numbers)": dict(annulus_orientation="matlab", annulus_diameter_model="matlab", node="gas_junction"),
    "1: discharge at platform node": dict(annulus_orientation="matlab", annulus_diameter_model="matlab", node="gas_inlet"),
    "2: physical annulus orientation": dict(annulus_orientation="physical", annulus_diameter_model="matlab", node="gas_junction"),
    "3: hydraulic annulus diameter": dict(annulus_orientation="matlab", annulus_diameter_model="hydraulic", node="gas_junction"),
    "1+2+3 combined": dict(annulus_orientation="physical", annulus_diameter_model="hydraulic", node="gas_inlet"),
}
rows = {}
for name, v in variants.items():
    case = case1_single_well(q_oil=q_oil, q_gas=q_gas, annulus_orientation=v["annulus_orientation"], annulus_diameter_model=v["annulus_diameter_model"])
    sol = case.simulator(sim_opts).solve()
    ec = evaluate_economics(sol, case.wells, CostModel(include_fuel_gas_cost=True, discharge_pressure_node=v["node"]))
    rows[name] = {"Downhole pressure (bar)": ec.downhole_pressure[0], "Injection point pressure (bar)": sol.pressure_bar("well_1:injection_point"),
                  "Annulus top pressure (bar)": sol.pressure_bar("well_1:annulus_top"), "Gas junction, node 9 (bar)": sol.pressure_bar("well_1:gas_junction"),
                  "Platform gas inlet, node 10 (bar)": sol.pressure_bar("well_1:gas_inlet"), "Discharge pressure used (bar)": ec.discharge_pressure,
                  "Compressor power (kW)": ec.compressor_power, "Fuel gas cost (10⁶ USD/yr)": ec.fuel_gas_cost}
t = pd.DataFrame(rows)
t["Table 13 (Max, paper)"] = pd.Series({"Downhole pressure (bar)": published["table13"]["rows"]["downhole_pressure_bar"][3],
                                        "Discharge pressure used (bar)": published["table13"]["rows"]["discharge_pressure_bar"][3],
                                        "Compressor power (kW)": published["table13"]["rows"]["compressor_power_kW"][3]})
t
"""
        ),
        md(
            """
Com a orientação física do anular, a pressão no topo do anular cai cerca de 2× a coluna de
gás de 1710 m (≈ 43 bar) e a pressão no nó 9 cai ≈ 31 bar; lendo a descarga na plataforma,
subtrai-se ainda a coluna da linha vertical de 1515 m (≈ 18 bar). Como o trabalho de compressão
depende fracamente da razão de pressões (expoente $(\\gamma-1)/(n\\gamma) \\approx 0{,}1$), a
potência correspondente ao ponto da Tabela 13 cai apenas de ≈ 1108 kW para ≈ 1069 kW (anular
físico) ou ≈ 1049 kW (as três correções juntas, −5 %). O custo com gás combustível e os
custos de capital do compressor são superestimados nessa mesma proporção; as conclusões
qualitativas do artigo não mudam.
"""
        ),
        md("## Itens 4–6 · Correlações de custo (Tabela 18, cenário 1)"),
        code(
            """
k = matlab["case2_design_220bar"]
d_in, lengths = np.array(k["xopt"][8:12]), [6020.0126, 6020.0126, 4520.0126, 4520.0126]
power = k["out_CompressorPower"]
q_gas_total = sum(k["q_gas"])
pg = k["out_DischargePressure"]
pub18 = published["table18"]["rows"]
from gaslift.economics import compressor_power
rows = {
    "paper (Table 18, scenario 1)": {"Compressor power (kW)": pub18["compressor_power_kW"][0], "CC (10⁶ USD)": pub18["compressor_cost_1e6_usd"][0], "CT (10⁶ USD)": pub18["turbine_cost_1e6_usd"][0], "CP (10⁶ USD)": pub18["pipeline_cost_1e6_usd"][0]},
    "code (paper numbers)": {"Compressor power (kW)": power, "CC (10⁶ USD)": compressor_capital_cost(power, "hp"), "CT (10⁶ USD)": turbine_capital_cost(power, "hp"), "CP (10⁶ USD)": pipeline_capital_cost(d_in, lengths, "matlab")},
    "4: McCoy with D^a7": {"Compressor power (kW)": power, "CC (10⁶ USD)": compressor_capital_cost(power, "hp"), "CT (10⁶ USD)": turbine_capital_cost(power, "hp"), "CP (10⁶ USD)": pipeline_capital_cost(d_in, lengths, "mccoy")},
    "5: γ = 1.25": {"Compressor power (kW)": compressor_power(pg, q_gas_total, CompressorParameters(polytropic_exponent=1.25)), "CC (10⁶ USD)": np.nan, "CT (10⁶ USD)": np.nan, "CP (10⁶ USD)": np.nan},
    "6: Couper correlations in kW": {"Compressor power (kW)": power, "CC (10⁶ USD)": compressor_capital_cost(power, "kW"), "CT (10⁶ USD)": turbine_capital_cost(power, "kW"), "CP (10⁶ USD)": pipeline_capital_cost(d_in, lengths, "matlab")},
}
t = pd.DataFrame(rows).T
t["ISBL (10⁶ USD)"] = t[["CC (10⁶ USD)", "CT (10⁶ USD)", "CP (10⁶ USD)"]].sum(axis=1, min_count=3)
t
"""
        ),
        md(
            """
## Item 7 · Geometria impressa nas Tabelas 2 e 7

Simulando o Caso 1 com a geometria impressa no artigo (coluna de produção de 5 in, anular de
1500 m, trecho inferior de 300 m, rugosidade de 100 µm nas linhas de gás) os valores da
Fig. 6 mudam substancialmente — o que confirma que os números publicados vieram do
`case_03.m` e não da Tabela 2.
"""
        ),
        code(
            """
points = [(10.0, 2.0), (5.0, 1.0), (15.0, 3.0)]
rows = {}
for qo, qg in points:
    a = case1_single_well(q_oil=qo, q_gas=qg).simulator(sim_opts).solve()
    b = case1_single_well(q_oil=qo, q_gas=qg, wellbore_diameter_in=5.0, annulus_length=1500.0, wellbore_length=1800.0, gas_line_roughness=100e-6).simulator(sim_opts).solve()
    rows[f"oil {qo:g}, gas {qg:g} kg/s"] = {"Downhole, code geometry (bar)": a.pressure_bar("well_1:downhole"), "Downhole, Table 2 geometry (bar)": b.pressure_bar("well_1:downhole"),
                                             "Node 9, code geometry (bar)": a.pressure_bar("well_1:gas_junction"), "Node 9, Table 2 geometry (bar)": b.pressure_bar("well_1:gas_junction")}
pd.DataFrame(rows).T
"""
        ),
        md(
            """
## Contagem de equações

Para o Caso 2 o artigo informa 232 equações; a formulação da Seção 3.3 produz
$3N + 4S = 3 \\cdot 28 + 4 \\cdot 27 = 192$ (o Caso 1 tem $3 \\cdot 10 + 4 \\cdot 9 = 66$, como publicado).
"""
        ),
        code(
            """
for case in (case1_single_well(), case2_multi_well()):
    lay = case.network.layout
    print(f"{case.network.name}: N = {lay.n_nodes}, S = {lay.n_pipes}, 3N + 4S = {lay.size}, specifications N + S = {lay.n_specifications}")
"""
        ),
    ]
    return write(NAMES["06"], cells)


BUILDERS = [nb01, nb02, nb03, nb04, nb05, nb06]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execute", action="store_true", help="execute the notebooks in place after writing them"
    )
    parser.add_argument("--only", nargs="*", help="prefixes of the notebooks to build (e.g. 02 05)")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    NB_DIR.mkdir(exist_ok=True)
    paths = []
    for build in BUILDERS:
        name = build.__name__[2:]
        if args.only and not any(name.startswith(p) for p in args.only):
            continue
        paths.append(build())
        print("wrote", paths[-1].relative_to(ROOT))
    if args.execute:
        for path in paths:
            cmd = [
                sys.executable,
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--inplace",
                f"--ExecutePreprocessor.timeout={args.timeout}",
                str(path),
            ]
            print("executing", path.name, flush=True)
            subprocess.run(cmd, check=True, cwd=NB_DIR)


if __name__ == "__main__":
    main()
