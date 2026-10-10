"""Generate the Jupyter notebooks of ``notebooks/`` from Python sources.

The notebooks are kept as generated artefacts so that their content is
reviewable as code. Run::

    uv run python scripts/build_notebooks.py          # write the .ipynb files
    uv run python scripts/build_notebooks.py --execute  # also execute them in place
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

SETUP = """import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Image, display

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / "pyproject.toml").exists())
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
REF = ROOT / "data" / "reference"
FIG = ROOT / "docs" / "paper" / "figures"
pd.set_option("display.precision", 4)
pd.set_option("display.width", 140)
plt.rcParams.update({"figure.dpi": 100, "figure.figsize": (7, 5), "axes.grid": True, "grid.alpha": 0.3})
"""


def md(text: str):
    return new_markdown_cell(text.strip("\n"))


def code(text: str):
    return new_code_cell(text.strip("\n"))


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
            """
# 01 · Modelo de escoamento bifásico: Beggs & Brill e balanço de momento

Este notebook apresenta os blocos de base do simulador: as propriedades dos fluidos
(Tabela 1 do artigo), o mapa de padrões de escoamento da correlação de Beggs & Brill
(1973) (Fig. 2) e a integração do balanço de momento ao longo de um duto (Eqs. 1–4 e 8–10).

Referência: de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010).
*Modeling, simulation and optimization of continuous gas lift systems for deepwater
offshore petroleum production*. J. Pet. Sci. Eng., doi:10.1016/j.petrol.2010.03.028.
"""
        ),
        code(SETUP),
        md("## Tabela 1 · Propriedades dos fluidos"),
        code(
            """
from gaslift import FluidProperties, paper_fluid

fluid = paper_fluid()
table1 = pd.DataFrame(
    {
        "Oil": {"Density (kg/m³)": fluid.liquid_density, "Viscosity (Pa s)": fluid.liquid_viscosity,
                "Interfacial tension (N/m)": fluid.interfacial_tension},
        "Gas": {"Reference density (kg/m³)": fluid.gas_reference_density,
                "Normal density (kg/m³)": fluid.gas_normal_density,
                "Reference pressure (bar)": fluid.reference_pressure / 1e5,
                "Polytropic coefficient (-)": fluid.gas_polytropic_coefficient,
                "Viscosity (Pa s)": fluid.gas_viscosity},
    }
)
table1
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
plt.xlabel("Pressure (bar)"); plt.ylabel("Gas density (kg/m³)"); plt.title("Polytropic gas density");
"""
        ),
        md(
            """
## Fig. 2 · Mapa de padrões de escoamento horizontal (Beggs & Brill)

Os testes lógicos da correlação original são substituídos por funções logísticas muito
abruptas (fator $10^6$), como no código MATLAB que gerou o artigo, para que o modelo seja
contínuo. O mapa abaixo é construído exatamente com essa função.
"""
        ),
        code(
            """
from gaslift.plotting import plot_flow_pattern_map

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
plot_flow_pattern_map(ax1)
ax1.set_title("Fig. 2 (reproduced)")
ax2.imshow(plt.imread(FIG / "FIG_02.png")); ax2.axis("off"); ax2.set_title("Fig. 2 (paper)")
plt.tight_layout()
"""
        ),
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
    return write("01_modelo_escoamento_bifasico.ipynb", cells)


# =============================================================================
# 02 - case 1 simulation
# =============================================================================
TABLES_CASE_CODE = """
from gaslift.units import INCH

def pipes_table(case):
    net = case.network
    rows = []
    for p in net.pipes:
        g = p.geometry
        angles = np.degrees(g.segment_angles)
        incl = f"{angles[0]:.0f}°" if np.allclose(angles, angles[0]) else "profile"
        rows.append({"Line": case.paper_pipe_numbers[p.name], "Name": p.name, "Type": p.category.value,
                     "Diameter (in)": g.inlet_diameter / INCH, "Length (m)": g.length,
                     "Inclination": incl, "Roughness (µm)": g.roughness[0] * 1e6})
    return pd.DataFrame(rows).set_index("Line")

def nodes_table(case):
    net = case.network
    rows = [{"Node": case.paper_node_numbers[n.name], "Name": n.name, "Role": n.role.value,
             "X (m)": n.x, "Y (m)": n.y, "Height (m)": n.z} for n in net.nodes]
    return pd.DataFrame(rows).set_index("Node")

def profile_table(case, pipe):
    g = case.network.pipe(pipe).geometry
    return pd.DataFrame({"Length (m)": g.cumulative_lengths, "Height (m)": g.points[:, 2],
                         "x (m)": g.points[:, 0], "y (m)": g.points[:, 1]}).T

def specs_table(case):
    rows = []
    for name, s in case.specifications.items():
        parts = []
        if s.pressure is not None: parts.append(f"P = {s.pressure / 1e5:g} bar")
        if s.liquid_flow is not None: parts.append(f"wL = {s.liquid_flow:g} kg/s")
        if s.gas_flow is not None: parts.append(f"wG = {s.gas_flow:g} kg/s")
        rows.append({"Node": case.paper_node_numbers[name], "Name": name, "Specifications": "; ".join(parts)})
    return pd.DataFrame(rows).set_index("Node")
"""


def nb02() -> Path:
    cells = [
        md(
            """
# 02 · Caso 1 — poço offshore único: rede, digrafo e simulação

Reprodução da Seção 4.1 do artigo: Tabelas 2–6, Figs. 3–5 (topologia e digrafo) e
Figs. 6–8 (curvas de nível da pressão de fundo e da pressão de injeção de gás).

A rede é construída com a API do pacote `gaslift`: nós e linhas da rede de coleta mais um
modelo de poço (`WellDesign`) que se expande em 6 nós e 6 linhas (Fig. 5).
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import case1_single_well, SimulationOptions, ModelOptions, IntegrationOptions
from gaslift.plotting import plot_case_topology, plot_digraph, plot_rate_contours, plot_pressure_path
from gaslift.sweeps import sweep_rates

case = case1_single_well(q_oil=10.0, q_gas=2.0)
net = case.network
print(net)
print(f"Unknowns / equations: 3N + 4S = {net.layout.size}  (N = {net.n_nodes}, S = {net.n_pipes});"
      f" specifications: {case.specifications.count}")
"""
        ),
        md("## Tabelas 2, 3 e 4 · Linhas da rede e perfis das linhas 2 e 8"),
        code(TABLES_CASE_CODE + "\npipes_table(case)"),
        md(
            """
> Os valores seguem os scripts MATLAB que geraram os resultados publicados (`case_03.m`):
> coluna de produção de 4 in, anular de 1710 m, trecho inferior de 338 m e rugosidade nula
> nas linhas de gás. A Tabela 2 impressa no artigo lista 5 in, 1500 m, 300 m e 100 µm —
> ver `docs/paper_vs_code.md`.
"""
        ),
        code(
            """
display(profile_table(case, "production_line_2").rename_axis("Line 2 (Table 3)"))
display(profile_table(case, "well_1:gas_line_1").rename_axis("Line 8 (Table 4)"))
"""
        ),
        md("## Tabelas 5 e 6 · Nós e especificações"),
        code("display(nodes_table(case))\ndisplay(specs_table(case))"),
        md("## Fig. 3 · Projeção isométrica"),
        code(
            """
fig = plot_case_topology(net)
fig.suptitle("Fig. 3 — Case 1, simplified isometric projection", y=1.02);
"""
        ),
        code('display(Image(filename=str(FIG / "FIG_03.png"), width=900))'),
        md("## Figs. 4 e 5 · Digrafo da rede e do poço"),
        code(
            """
fig, ax = plt.subplots(figsize=(12, 7))
plot_digraph(net, ax, node_numbers=case.paper_node_numbers, pipe_numbers=case.paper_pipe_numbers)
ax.set_title("Fig. 4/5 — Case 1 digraph (node numbers as in the paper; well nodes highlighted)");
"""
        ),
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
integração (`results_Subsea_1_2.mat`), por isso usamos o mesmo valor aqui. A varredura
percorre a grade em serpentina com partida quente, como o script original.
"""
        ),
        code(
            """
q_gas = np.arange(0.5, 4.01, 0.5)
q_oil = np.concatenate([[0.5], np.arange(1.0, 20.01, 1.0)])
opts = SimulationOptions(model=ModelOptions(integration=IntegrationOptions(rtol=1e-3)), compute_profiles=False)
import time
t0 = time.perf_counter()
sweep = sweep_rates(case1_single_well(), q_gas, q_oil, opts)
print(f"{sweep.converged.sum()} of {sweep.converged.size} points converged in {time.perf_counter() - t0:.1f} s")
"""
        ),
        code(
            """
levels = np.arange(140, 320, 20)
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), axes[0], levels, title="Fig. 6 — Downhole pressure, node 7 (bar)")
axes[1].imshow(plt.imread(FIG / "FIG_06.png")); axes[1].axis("off"); axes[1].set_title("Fig. 6 (paper)")
plt.tight_layout()
"""
        ),
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
levels = np.arange(200, 560, 50)
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.gas_junction(0), axes[0], levels, title="Fig. 7 — Gas lift inlet pressure, node 9 (bar), 2 in gas line")
axes[1].imshow(plt.imread(FIG / "FIG_07.png")); axes[1].axis("off"); axes[1].set_title("Fig. 7 (paper)")
plt.tight_layout()
"""
        ),
        code(
            """
fig, ax = plt.subplots(figsize=(7, 6))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.gas_inlet(0), ax, np.arange(150, 560, 50), title="Platform gas inlet pressure, node 10 (bar), 2 in gas line");
"""
        ),
        md("## Fig. 8 · Linhas de gás de 3 in"),
        code(
            """
sweep3 = sweep_rates(case1_single_well(gas_line_diameter_in=3.0), q_gas, q_oil, opts)
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
plot_rate_contours(sweep3.q_gas, sweep3.q_oil, sweep3.gas_junction(0), axes[0], np.arange(180, 320, 20), title="Fig. 8 — Gas lift inlet pressure, node 9 (bar), 3 in gas line")
axes[1].imshow(plt.imread(FIG / "FIG_08.png")); axes[1].axis("off"); axes[1].set_title("Fig. 8 (paper)")
plt.tight_layout()
"""
        ),
        md(
            """
## Leitura da Fig. 6 (texto do artigo)

> "considering a downhole pressure of 220 bar, to produce 10.0 kg/s of petroleum it is
> necessary to inject 1.0 kg/s of gas or a higher value of gas flow rate, equal to 2.7 kg/s
> [...] at downhole pressure of 220 bar, the maximum oil production is 10.7 kg/s with gas
> injection of 1.8 kg/s."
"""
        ),
        code(
            """
from scipy.interpolate import RegularGridInterpolator

interp = RegularGridInterpolator((sweep.q_oil, sweep.q_gas), sweep.downhole(0), bounds_error=False)
qg_fine = np.linspace(0.5, 4.0, 701)
p_at_10 = interp(np.column_stack([np.full_like(qg_fine, 10.0), qg_fine]))
roots = qg_fine[np.flatnonzero(np.diff(np.sign(p_at_10 - 220.0)))]
print("gas rates giving 220 bar at 10 kg/s of oil:", np.round(roots, 2), "kg/s")

qo_fine = np.linspace(0.5, 20.0, 781)
QO, QG = np.meshgrid(qo_fine, qg_fine, indexing="ij")
P = interp(np.column_stack([QO.ravel(), QG.ravel()])).reshape(QO.shape)
feasible = np.where(P <= 220.0, QO, np.nan)
i, j = np.unravel_index(np.nanargmax(feasible), feasible.shape)
print(f"maximum oil at 220 bar: {qo_fine[i]:.2f} kg/s with {qg_fine[j]:.2f} kg/s of gas")
"""
        ),
        md("## Validação contra os resultados MATLAB originais"),
        code(
            """
import json
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
    return write("02_caso1_poco_unico_simulacao.ipynb", cells)


# =============================================================================
# 03 - case 1 optimisation
# =============================================================================
def nb03() -> Path:
    cells = [
        md(
            """
# 03 · Caso 1 — otimização: produção máxima e lucro máximo

Reprodução da Seção 6.1 do artigo: Figs. 13–16 e Tabela 13. O problema de otimização
(Eq. 24) acopla o simulador de rede (que especifica vazões) ao modelo de reservatório
$Q_L = PI\\,(P_{res} - P_{downhole})$ por meio de uma restrição de igualdade, resolvido com
programação quadrática sequencial (SLSQP do SciPy, equivalente ao `fmincon` usado no artigo).
"""
        ),
        code(SETUP),
        code(
            """
import json, time
from gaslift import case1_single_well, SimulationOptions, ModelOptions, IntegrationOptions
from gaslift.optimization import GasLiftProblem, Objective, Reservoir
from gaslift.economics import CostModel
from gaslift.plotting import plot_rate_contours
from gaslift.sweeps import sweep_rates

matlab = json.load(open(REF / "matlab_optimization_results.json"))

def optimise(objective, reservoir, cost=None, x0=(10.0, 2.0)):
    case = case1_single_well()
    prob = GasLiftProblem(case, [reservoir], objective, cost=cost or CostModel(), oil_bounds=(5.0, 20.0), gas_bounds=(1.0, 4.0))
    t0 = time.perf_counter()
    res = prob.solve(prob.initial_vector(*x0))
    res.elapsed = time.perf_counter() - t0
    return res
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
pressures = [180, 200, 220, 240, 260]
rows = []
for p in pressures:
    res = optimise(Objective.MAX_PRODUCTION, Reservoir(p, np.inf))
    ref = matlab[f"case1_maxprod_piinf_{p}bar"]
    rows.append({"P_well (bar)": p, "oil (kg/s)": res.oil_rates[0], "gas (kg/s)": res.gas_rates[0],
                 "oil MATLAB": ref["xopt"][0], "gas MATLAB": ref["xopt"][1],
                 "evaluations": res.n_evaluations, "time (s)": res.elapsed})
fig13 = pd.DataFrame(rows).set_index("P_well (bar)")
fig13
"""
        ),
        code(
            """
fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
axes[0].plot(fig13.index, fig13["oil (kg/s)"], "ko-", label="Python (SLSQP)"); axes[0].plot(fig13.index, fig13["oil MATLAB"], "r+", ms=10, label="MATLAB (paper)")
axes[0].set_ylabel("Produced oil (kg/s)"); axes[0].legend(); axes[0].set_title("Fig. 13 — Maximum production operational conditions")
axes[1].plot(fig13.index, fig13["gas (kg/s)"], "ko-"); axes[1].plot(fig13.index, fig13["gas MATLAB"], "r+", ms=10)
axes[1].set_ylabel("Injected gas (kg/s)"); axes[1].set_xlabel("Downhole pressure (bar)")
plt.tight_layout()
"""
        ),
        code(
            """
q_gas = np.arange(0.5, 4.01, 0.5)
q_oil = np.concatenate([[0.5], np.arange(1.0, 20.01, 1.0)])
opts = SimulationOptions(model=ModelOptions(integration=IntegrationOptions(rtol=1e-3)), compute_profiles=False)
sweep = sweep_rates(case1_single_well(), q_gas, q_oil, opts)
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), axes[0], np.arange(140, 320, 20), title="Fig. 14 — Maximum production locus")
axes[0].plot(fig13["gas (kg/s)"], fig13["oil (kg/s)"], "ro--", label="Optimum locus"); axes[0].legend(loc="upper right")
axes[1].imshow(plt.imread(FIG / "FIG_14.png")); axes[1].axis("off"); axes[1].set_title("Fig. 14 (paper)")
plt.tight_layout()
"""
        ),
        md("## 6.1.2 · Produção máxima para diferentes índices de produtividade (Fig. 15)"),
        code(
            """
pis = [0.5, 1.0, 2.0, 5.0, 10.0]
rows = []
for pi in pis:
    res = optimise(Objective.MAX_PRODUCTION, Reservoir(220.0, pi))
    ref = matlab[f"case1_maxprod_pi{pi:g}_220bar"]
    rows.append({"PI ((kg/s)/bar)": pi, "oil (kg/s)": res.oil_rates[0], "gas (kg/s)": res.gas_rates[0],
                 "P_downhole (bar)": res.economics.downhole_pressure[0],
                 "oil MATLAB": ref["xopt"][0], "gas MATLAB": ref["xopt"][1], "P_downhole MATLAB": ref["p_downhole_bar"],
                 "time (s)": res.elapsed})
fig15 = pd.DataFrame(rows).set_index("PI ((kg/s)/bar)")
fig15
"""
        ),
        code(
            """
fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
for ax, col, ref, lab in zip(axes, ["oil (kg/s)", "gas (kg/s)", "P_downhole (bar)"], ["oil MATLAB", "gas MATLAB", "P_downhole MATLAB"],
                             ["Produced oil (kg/s)", "Injected gas (kg/s)", "Well pressure (bar)"]):
    ax.plot(fig15.index, fig15[col], "ko-", label="Python"); ax.plot(fig15.index, fig15[ref], "r+", ms=10, label="MATLAB"); ax.set_ylabel(lab)
axes[0].legend(); axes[0].set_title("Fig. 15 — Maximum production vs productivity index (220 bar)"); axes[-1].set_xlabel("Productivity index ((kg/s)/bar)")
plt.tight_layout()
"""
        ),
        md(
            """
## 6.1.3 · Lucro máximo (Tabela 13 e Fig. 16)

Função objetivo $\\Omega_2 = R - AGC$: receita do óleo menos custo do gás combustível da
turbina do compressor (Eqs. 26, 31, 34 e 35). Os cenários variam o preço do óleo com o preço
do gás fixo em 0,1886 USD/Nm³; o cenário "Max" corresponde à produção máxima (custo nulo).
"""
        ),
        code(
            """
scenarios = [("1", 51.15 / 4, Objective.MAX_PROFIT), ("2", 51.15 / 2, Objective.MAX_PROFIT), ("3", 51.15, Objective.MAX_PROFIT), ("Max", 51.15, Objective.MAX_PRODUCTION)]
keys = {"1": "case1_maxprofit_oil12.79", "2": "case1_maxprofit_oil25.58", "3": "case1_maxprofit_oil51.15", "Max": "case1_maxprod_pi1_220bar"}
rows = {}
for name, oil_price, obj in scenarios:
    res = optimise(obj, Reservoir(220.0, 1.0), CostModel(oil_price=oil_price))
    ec = res.economics
    gas_price = 0.0 if obj is Objective.MAX_PRODUCTION else ec.fuel_gas_cost / max(ec.fuel_gas_volume, 1e-12)
    ref = matlab[keys[name]]
    rows[name] = {"Cost ratio ×10³ ((USD/Nm³)/(USD/bbl))": 1e3 * gas_price / oil_price, "Fuel gas cost (USD/Nm³)": gas_price,
                  "Oil price (USD/bbl)": oil_price, "Produced oil (kg/s)": ec.oil_rate[0], "Injected gas (kg/s)": ec.gas_rate[0],
                  "Downhole pressure (bar)": ec.downhole_pressure[0], "Compressor power (kW)": ec.compressor_power,
                  "Discharge pressure PG (bar)": ec.discharge_pressure, "Profit (10⁶ USD/yr)": ec.profit,
                  "Profit MATLAB": -ref["fob"], "Injected gas MATLAB": ref["xopt"][1]}
table13 = pd.DataFrame(rows)
table13
"""
        ),
        code(
            """
ratio = table13.loc["Cost ratio ×10³ ((USD/Nm³)/(USD/bbl))"] / 1e3
fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
axes[0].plot(ratio, table13.loc["Produced oil (kg/s)"], "ko-"); axes[0].set_ylabel("Produced oil (kg/s)"); axes[0].set_title("Fig. 16 — Maximum profit operational conditions")
axes[1].plot(ratio, table13.loc["Injected gas (kg/s)"], "ko-"); axes[1].set_ylabel("Injected gas (kg/s)")
axes[2].plot(ratio, table13.loc["Profit (10⁶ USD/yr)"], "ko-"); axes[2].set_ylabel("Profit (10⁶ USD/yr)"); axes[2].set_xlabel("Cost ratio ((USD/Nm³)/(USD/barrel))")
plt.tight_layout()
"""
        ),
        code('display(Image(filename=str(FIG / "FIG_16.png"), width=500))'),
    ]
    return write("03_caso1_otimizacao.ipynb", cells)


# =============================================================================
# 04 - case 2 simulation
# =============================================================================
def nb04() -> Path:
    cells = [
        md(
            """
# 04 · Caso 2 — múltiplos poços offshore: rede e simulação

Reprodução da Seção 4.2 do artigo: Tabelas 7–12, Figs. 9–10 (topologia e digrafo) e
Figs. 11–12 (curvas de nível da pressão de fundo dos poços 1/2 e 3/4). Quatro poços com
gas lift compartilham a linha de produção e o riser: os poços 1 e 2 entram no nó 1 e os
poços 3 e 4 no nó 2.
"""
        ),
        code(SETUP),
        code(
            """
from gaslift import case2_multi_well, SimulationOptions, ModelOptions, IntegrationOptions
from gaslift.plotting import plot_case_topology, plot_digraph, plot_rate_contours
from gaslift.sweeps import sweep_rates

case = case2_multi_well()
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
> formulação da Seção 3.3.
"""
        ),
        md("## Tabelas 7–12 · Linhas, perfis, nós e especificações"),
        code(TABLES_CASE_CODE + "\npipes_table(case)"),
        code(
            """
display(profile_table(case, "well_1:gas_line_1").rename_axis("Lines 8 and 14 (Table 9)"))
display(profile_table(case, "well_3:gas_line_1").rename_axis("Lines 20 and 26 (Table 10)"))
"""
        ),
        md(
            """
> A Tabela 10 do artigo lista comprimentos acumulados 0, 1005, 2505 e 4005 m para as linhas
> 20 e 26; a geometria real (e a Tabela 16, que usa $L_{G3} = L_{G4} = 4520$ m) corresponde a
> 3005 m + 1515 m.
"""
        ),
        code("display(nodes_table(case))\ndisplay(specs_table(case))"),
        md("## Fig. 9 · Projeção isométrica"),
        code(
            """
fig = plot_case_topology(net)
fig.suptitle("Fig. 9 — Case 2, simplified isometric projection", y=1.02);
display(Image(filename=str(FIG / "FIG_09.png"), width=900))
"""
        ),
        md("## Fig. 10 · Digrafo"),
        code(
            """
fig, ax = plt.subplots(figsize=(16, 8))
plot_digraph(net, ax, node_numbers=case.paper_node_numbers, pipe_numbers=case.paper_pipe_numbers)
ax.set_title("Fig. 10 — Case 2 digraph (node numbers as in the paper)");
"""
        ),
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
(256 simulações de 192 incógnitas; tolerância de integração `1e-2`, como no original).
"""
        ),
        code(
            """
import time
q_gas = np.arange(0.25, 4.01, 0.25)
q_oil = np.arange(0.5, 8.01, 0.5)
opts = SimulationOptions(model=ModelOptions(integration=IntegrationOptions(rtol=1e-2)), compute_profiles=False)
t0 = time.perf_counter()
sweep = sweep_rates(case2_multi_well(), q_gas, q_oil, opts)
print(f"{sweep.converged.sum()} of {sweep.converged.size} points converged in {time.perf_counter() - t0:.1f} s")
"""
        ),
        code(
            """
levels = np.arange(140, 380, 20)
fig, axes = plt.subplots(2, 2, figsize=(15, 13))
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(0), axes[0, 0], levels, title="Fig. 11 — Pressure at nodes 7 and 13 (bar)")
axes[0, 1].imshow(plt.imread(FIG / "FIG_11.png")); axes[0, 1].axis("off"); axes[0, 1].set_title("Fig. 11 (paper)")
plot_rate_contours(sweep.q_gas, sweep.q_oil, sweep.downhole(2), axes[1, 0], levels, title="Fig. 12 — Pressure at nodes 19 and 25 (bar)")
axes[1, 1].imshow(plt.imread(FIG / "FIG_12.png")); axes[1, 1].axis("off"); axes[1, 1].set_title("Fig. 12 (paper)")
plt.tight_layout()
"""
        ),
        md("## Validação contra os resultados MATLAB originais"),
        code(
            """
import json
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
    return write("04_caso2_multipocos_simulacao.ipynb", cells)


# =============================================================================
# 05 - case 2 optimisation
# =============================================================================
def nb05() -> Path:
    cells = [
        md(
            """
# 05 · Caso 2 — otimização: produção máxima e projeto ótimo do sistema de gas lift

Reprodução da Seção 6.2 do artigo: Tabelas 14 e 15 (produção máxima para reservatórios a
220 e 240 bar) e Tabelas 16–19 (projeto ótimo com custos anualizados de compressor, turbina,
linhas de gás e gás combustível). As variáveis de decisão do projeto são as vazões de óleo e
gás de cada poço e o diâmetro das linhas de gás de cada poço (12 variáveis).
"""
        ),
        code(SETUP),
        code(
            """
import json, time
from gaslift import case2_multi_well
from gaslift.optimization import GasLiftProblem, Objective, Reservoir
from gaslift.economics import CostModel, CompressorParameters

matlab = json.load(open(REF / "matlab_optimization_results.json"))
PI = [1.0, 2.5, 1.0, 0.5]
"""
        ),
        md("## 6.2.1 · Produção máxima (Tabelas 14 e 15)"),
        code(
            """
def max_production(p_res):
    case = case2_multi_well()
    prob = GasLiftProblem(case, [Reservoir(p_res, pi) for pi in PI], Objective.MAX_PRODUCTION, oil_bounds=(0.2, 20.0), gas_bounds=(0.2, 10.0))
    t0 = time.perf_counter()
    res = prob.solve(prob.initial_vector(3.0, 1.0))
    print(f"P_res = {p_res} bar: {res.n_evaluations} simulations in {time.perf_counter() - t0:.0f} s; total oil {res.oil_rates.sum():.3f} kg/s")
    return res

def production_table(res, key):
    ref = matlab[key]
    df = pd.DataFrame({"Well": [1, 2, 3, 4], "Node": [7, 13, 19, 25], "Productivity index ((kg/s)/bar)": PI,
                       "Downhole pressure (bar)": res.economics.downhole_pressure, "Produced oil (kg/s)": res.oil_rates,
                       "Injected gas (kg/s)": res.gas_rates, "Produced oil MATLAB": ref["q_oil"], "Injected gas MATLAB": ref["q_gas"]}).set_index("Well")
    return df.T

res220 = max_production(220.0)
res240 = max_production(240.0)
"""
        ),
        code(
            'display(production_table(res220, "case2_maxprod_220bar").rename_axis("Table 14 — 220 bar"))'
        ),
        code(
            'display(production_table(res240, "case2_maxprod_240bar").rename_axis("Table 15 — 240 bar"))'
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
"""
        ),
        code(
            """
comp = CompressorParameters()
cost = CostModel()
table16 = pd.Series({
    "Number of compressor stages (n)": comp.stages, "Overall compressor efficiency (η)": comp.efficiency,
    "Polytropic coefficient used in Eq. 31 (γ)": comp.polytropic_exponent, "Suction temperature (K)": comp.suction_temperature,
    "Suction pressure (bar)": comp.suction_pressure, "Molecular weight (kg/kmol)": comp.molecular_weight,
    "Fuel gas cost (USD/Nm³)": cost.gas_price, "Oil price (USD/barrel)": cost.oil_price,
    "Pipe length LG1, LG2 (m)": case2_multi_well().gas_line_length("well_1"), "Pipe length LG3, LG4 (m)": case2_multi_well().gas_line_length("well_3"),
    "Diameter bounds (in)": "1 – 4", "Productivity indexes": str(PI), "Reservoir pressure scenarios (bar)": "220, 240"}, name="Table 16")
table16.to_frame()
"""
        ),
        code(
            """
def optimal_design(p_res, start):
    case = case2_multi_well()
    prob = GasLiftProblem(case, [Reservoir(p_res, pi) for pi in PI], Objective.OPTIMAL_DESIGN,
                          oil_bounds=(0.5, 20.0), gas_bounds=(0.1, 4.0), diameter_bounds=(1.0, 4.0))
    x0 = prob.initial_vector(start.oil_rates, start.gas_rates, 2.0)
    t0 = time.perf_counter()
    res = prob.solve(x0)
    print(f"P_res = {p_res} bar: {res.n_evaluations} simulations in {time.perf_counter() - t0:.0f} s; profit {res.economics.profit:.2f} 10⁶ USD/yr")
    return res

des220 = optimal_design(220.0, res220)
des240 = optimal_design(240.0, res240)
"""
        ),
        code(
            """
def design_tables(r1, r2):
    k1, k2 = matlab["case2_design_220bar"], matlab["case2_design_240bar"]
    dec = {}
    for i in range(4):
        dec[f"Injected gas QG{i+1} (kg/s)"] = (r1.gas_rates[i], r2.gas_rates[i], k1["xopt"][4 + i], k2["xopt"][4 + i])
    for i in range(4):
        dec[f"Oil flow rate QL{i+1} (kg/s)"] = (r1.oil_rates[i], r2.oil_rates[i], k1["xopt"][i], k2["xopt"][i])
    for i in range(4):
        dec[f"Pipe diameter DG{i+1} (in)"] = (r1.diameters_in[i], r2.diameters_in[i], k1["xopt"][8 + i], k2["xopt"][8 + i])
    t17 = pd.DataFrame(dec, index=["Scenario 1 (220 bar)", "Scenario 2 (240 bar)", "Scenario 1 MATLAB", "Scenario 2 MATLAB"]).T
    e1, e2 = r1.economics, r2.economics
    t18 = pd.DataFrame({
        "Scenario 1": [e1.compressor_power, e1.discharge_pressure, e1.compressor_cost, e1.turbine_cost, e1.pipeline_cost, e1.fuel_gas_volume],
        "Scenario 2": [e2.compressor_power, e2.discharge_pressure, e2.compressor_cost, e2.turbine_cost, e2.pipeline_cost, e2.fuel_gas_volume],
        "Scenario 1 MATLAB": [k1["out_CompressorPower"], k1["out_DischargePressure"], k1["out_CompressorCapitalCost"], k1["out_TurbineCapitalCost"], k1["out_PipesCapitalCost"], k1["out_FuelGasFlowRate"]],
        "Scenario 2 MATLAB": [k2["out_CompressorPower"], k2["out_DischargePressure"], k2["out_CompressorCapitalCost"], k2["out_TurbineCapitalCost"], k2["out_PipesCapitalCost"], k2["out_FuelGasFlowRate"]],
        "Unit": ["kW", "bar", "10⁶ USD", "10⁶ USD", "10⁶ USD", "10⁶ m³/yr"]},
        index=["Compressor power (wC)", "Compressor discharge pressure", "Compressor capital cost (CC)", "Turbine capital cost (CT)", "Gas pipelines capital cost (CP)", "Fuel gas flow rate"])
    t19 = pd.DataFrame({
        "Scenario 1": [e1.injected_gas_volume.sum(), e1.oil_volume.sum(), e1.fuel_gas_cost, e1.annualized_capital_cost, e1.isbl, e1.osbl, e1.revenue, e1.profit],
        "Scenario 2": [e2.injected_gas_volume.sum(), e2.oil_volume.sum(), e2.fuel_gas_cost, e2.annualized_capital_cost, e2.isbl, e2.osbl, e2.revenue, e2.profit],
        "Scenario 1 MATLAB": [sum(k1["out_GasLiftFlowRate"]), sum(k1["out_OilFlowRate"]), k1["out_AnnualFuelGasCost"], k1["out_AnnualTotalCapitalCost"], k1["out_ISBL"], k1["out_OSBL"], k1["out_InValue"], k1["out_Profit"]],
        "Scenario 2 MATLAB": [sum(k2["out_GasLiftFlowRate"]), sum(k2["out_OilFlowRate"]), k2["out_AnnualFuelGasCost"], k2["out_AnnualTotalCapitalCost"], k2["out_ISBL"], k2["out_OSBL"], k2["out_InValue"], k2["out_Profit"]],
        "Unit": ["10⁶ m³/yr", "10⁶ barrel/yr", "10⁶ USD/yr", "10⁶ USD/yr", "10⁶ USD", "10⁶ USD", "10⁶ USD/yr", "10⁶ USD/yr"]},
        index=["Total injected gas flow rate", "Total oil flow rate", "Fuel gas cost (AGC)", "Annualized capital cost", "ISBL", "OSBL", "Revenue (R)", "Profit"])
    return t17, t18, t19

table17, table18, table19 = design_tables(des220, des240)
display(table17.rename_axis("Table 17 — decision variables"))
"""
        ),
        code('display(table18.rename_axis("Table 18 — pipelines, compressor and turbine"))'),
        code('display(table19.rename_axis("Table 19 — optimal design results"))'),
        md(
            """
As pequenas diferenças em relação ao MATLAB (terceira casa decimal) decorrem do otimizador:
a superfície do lucro é muito plana em torno do ótimo, e `fmincon` (2009) e SLSQP param em
pontos ligeiramente distintos com valores de função praticamente iguais.
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
    return write("05_caso2_otimizacao.ipynb", cells)


# =============================================================================
# 06 - paper vs code
# =============================================================================
def nb06() -> Path:
    cells = [
        md(
            """
# 06 · Artigo × código: divergências encontradas na migração e seu efeito nos resultados

Ao portar o código MATLAB (`bbnet10`) que gerou os números publicados, algumas diferenças
entre o texto do artigo e o código foram identificadas (detalhes em `docs/paper_vs_code.md`).
O pacote reproduz o código por padrão e expõe interruptores para as formas "corrigidas".
Este notebook quantifica o efeito de cada uma delas.

| # | Item | Código (padrão, reproduz o artigo) | Forma descrita no artigo / física |
|---|------|------|------|
| 1 | Pressão de descarga $P_G$ | nó n+4 (interconexão das linhas de gás, −1495 m) | nó n+5 (plataforma, +20 m) |
| 2 | Orientação do anular | geometria escrita de baixo para cima (termo hidrostático de escoamento ascendente) | gás desce: −90° |
| 3 | Diâmetro equivalente do anular | $(D_o^2 + D_i^2)/(D_o + D_i)$ | $D_o - D_i$ (hidráulico) |
| 4 | Custo das linhas (McCoy) | $D \\cdot a_7$ | $D^{a_7}$ (Eq. 33) |
| 5 | Expoente politrópico no compressor | 1,4 | 1,25 (Tabela 16) |
| 6 | Potência nas correlações de custo | hp | kW (Eqs. 30 e 32) |
| 7 | Geometria do poço | 4 in; 1710 + 338 m; rugosidade 0 nas linhas de gás | 5 in; 1500 + 300 m; 100 µm (Tabela 2) |
"""
        ),
        code(SETUP),
        code(
            """
import json
from dataclasses import replace
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
pd.DataFrame(rows)
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
from gaslift.economics import compressor_power
rows = {
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
    return write("06_artigo_vs_codigo.ipynb", cells)


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
