# gaslift — simulação de redes bifásicas e otimização de gas lift contínuo

Reimplementação em Python do modelo e dos estudos de caso de

> de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010). **Modeling,
> simulation and optimization of continuous gas lift systems for deepwater offshore
> petroleum production.** *Journal of Petroleum Science and Engineering*,
> doi:[10.1016/j.petrol.2010.03.028](https://doi.org/10.1016/j.petrol.2010.03.028).

O código original (MATLAB, `bbnet10`) está preservado em `legacy/matlab/`. O pacote
`gaslift` reproduz os seus resultados com precisão de máquina (vetores de estado idênticos a
~1e-10) e os notebooks regeneram todas as tabelas, figuras e digrafos do artigo.

## O que há no modelo

* **Escoamento bifásico** em dutos com qualquer inclinação pela correlação de Beggs & Brill
  (1973), com padrões de escoamento, hold-up, fator de atrito de Churchill e compressibilidade
  do gás (politrópica) — `gaslift.beggs_brill`, `gaslift.pipe_model`.
* **Rede de dutos** como digrafo (N nós, S linhas): balanços de massa nodais, balanço de
  momento integrado em cada duto (nas duas direções) e balanços de momento nodais, resolvidos
  por Newton–Raphson a partir de um modelo linear laminar — `gaslift.network`,
  `gaslift.model`, `gaslift.linear_model`, `gaslift.solver`, `gaslift.simulation`.
* **Integrador Runge–Kutta em lote** (Dormand–Prince 5(4) com o controle de passo do
  `ode45`) que integra todos os dutos e todas as perturbações do Jacobiano de uma só vez —
  `gaslift.integrator`.
* **Poços com gas lift** como modelo reutilizável (`WellDesign`): coluna de produção, anular,
  linhas de gás, ponto de injeção e fundo do poço — 6 nós e 6 linhas por poço.
* **Economia e otimização**: receita do óleo, gás combustível, custos de capital de
  compressor, turbina e linhas de gás (Couper et al., McCoy), e os três problemas do artigo
  (produção máxima, lucro máximo e projeto ótimo) com SQP (SciPy/SLSQP) acoplado ao
  simulador — `gaslift.economics`, `gaslift.optimization`.

## Instalação e uso

```bash
uv sync --all-extras            # Python 3.12, numpy/scipy/pandas/matplotlib/networkx/jupyter
uv run pytest                   # testes rápidos (~10 s); `-m slow` para as otimizações
uv run jupyter lab notebooks/   # abrir os notebooks
```

```python
from gaslift import case1_single_well
from gaslift.optimization import GasLiftProblem, Objective, Reservoir

case = case1_single_well(q_oil=10.0, q_gas=2.0)        # rede do Caso 1 (Tabelas 2–6)
solution = case.simulator().solve()                     # estado estacionário
print(solution.wells_table())                           # pressões de fundo, injeção, plataforma

problem = GasLiftProblem(case, [Reservoir(pressure=220.0, productivity_index=1.0)],
                         Objective.MAX_PRODUCTION, oil_bounds=(5, 20), gas_bounds=(1, 4))
result = problem.solve(problem.initial_vector(10.0, 2.0))
print(result.summary())                                 # 9.598 kg/s de óleo com 1.812 kg/s de gás
```

Redes arbitrárias são montadas com `Network.add_node`, `Network.add_pipe` (polilinhas 3-D
com diâmetro e rugosidade por vértice) e `Network.add_well`, e especificadas com
`Specifications` (pressão e/ou vazões externas em cada nó; o modelo exige N + S
especificações). `Network.to_dict()/from_dict()` serializa a rede em JSON.

## Notebooks (`notebooks/`)

| Notebook | Conteúdo do artigo |
|---|---|
| `01_modelo_escoamento_bifasico` | Tabela 1, Fig. 2, integração do balanço de momento em um duto |
| `02_caso1_poco_unico_simulacao` | Tabelas 2–6, Figs. 3–8, validação contra o MATLAB |
| `03_caso1_otimizacao` | Figs. 13–16, Tabela 13 |
| `04_caso2_multipocos_simulacao` | Tabelas 7–12, Figs. 9–12, validação contra o MATLAB |
| `05_caso2_otimizacao` | Tabelas 14–19 |
| `06_artigo_vs_codigo` | Efeito das divergências entre o texto do artigo e o código original |

Os notebooks são gerados a partir de `scripts/build_notebooks.py`
(`uv run python scripts/build_notebooks.py --execute` reconstrói e executa todos; ~15 min).

## Validação

`data/reference/` contém os resultados MATLAB originais (varreduras das Figs. 6–8 e 11–12 e
os ótimos das Tabelas 13–19), exportados dos `.mat` por `scripts/export_matlab_reference.py`.
Os testes em `tests/` verificam que:

* os vetores de estado completos coincidem com o MATLAB (diferença máxima ~1e-10) em pontos
  das varreduras dos Casos 1 e 2 e em todos os ótimos publicados;
* potência de compressão, custos de capital, gás combustível, receita e lucro reproduzem as
  Tabelas 13, 18 e 19 em todos os dígitos impressos;
* a otimização SLSQP reencontra os ótimos do `fmincon` (produção máxima e lucro máximo).

## Divergências entre o artigo e o código original

Durante a migração foram encontradas diferenças entre o texto publicado e o código que gerou
os números (geometria dos poços, orientação do anular, nó usado como descarga do compressor,
forma da correlação de custo de dutos, expoente politrópico, contagem de equações).
O pacote reproduz o código por padrão e oferece opções para as formas do artigo;
ver `docs/paper_vs_code.md` e o notebook 06.

## Estrutura

```
src/gaslift/        pacote (ver CLAUDE.md para o mapa dos módulos)
tests/              pytest (rápidos + marcados `slow`)
notebooks/          notebooks executados
scripts/            geração dos notebooks e exportação das referências MATLAB
data/reference/     resultados MATLAB usados na validação
legacy/matlab/      código MATLAB original (bbnet10) e scripts dos casos
docs/paper/         artigo (prova tipográfica) e figuras originais
docs/raw/           material bruto do projeto (não versionado)
```
