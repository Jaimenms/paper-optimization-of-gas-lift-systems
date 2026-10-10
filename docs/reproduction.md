# Reproduction status

Paper: de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010).
**Modeling, simulation and optimization of continuous gas lift systems for deepwater
offshore petroleum production.** *Journal of Petroleum Science and Engineering* 72 (3–4),
277–289. doi:10.1016/j.petrol.2010.03.028.

The package reproduces the MATLAB program (`bbnet10`, `legacy/matlab/`) that generated every
number of the paper. `make run` + `make eval` regenerate the comparison below
(`results/eval/report.md`); the notebooks show the same items with commentary (Portuguese).
The formal errata of the paper are in [`errata.md`](errata.md); the Portuguese working notes
on the paper ↔ code differences are in [`pt/paper_vs_code.md`](pt/paper_vs_code.md).

## Status per table and figure

| Item | Content | Status | Where |
|---|---|---|---|
| Table 1 | fluid properties | input, reproduced | `FluidProperties`, notebook 01 |
| Fig. 1 | control volume | sketch, not reproduced | — |
| Fig. 2 | Beggs & Brill flow map | reproduced from the implemented pattern function | `plot_flow_pattern_map`, notebook 01, `figures/fig02_*` |
| Tables 2–6 | Case 1 lines, profiles, nodes, specifications | reproduced **with the values of the MATLAB script `case_03.m`** (diameter 4 in, annulus 1710 m, lower wellbore 338 m, gas-line roughness 0; the printed 5 in / 1500 m / 300 m / 100 µm were not what the program used — errata A1) | `cases/case1_single_well.json`, notebook 02 |
| Figs. 3–5 | Case 1 isometric view and digraphs | reproduced (schematic layout) | `plot_case_topology`, `plot_digraph` |
| Fig. 6 | Case 1 downhole pressure contours | reproduced; state vectors equal to the MATLAB sweep (`results_Subsea_1_2.mat`) to ~10⁻¹⁰ at every converged point (2 of the 168 grid points did not converge in MATLAB, 1 in Python); the text readings (1.0 / 2.7 kg/s of gas for 10 kg/s at 220 bar; maximum 10.7 kg/s at 1.8 kg/s) are recovered by interpolation of the grid | notebook 02, `figures/fig06_*` |
| Figs. 7–8 | Case 1 gas-lift inlet pressure, 2 in and 3 in gas lines | reproduced; the plotted pressure is that of node 9 (gas-lines interconnection), not node 10 as the text says (errata A4) | notebook 02, `figures/fig07_*`, `fig08_*` |
| Tables 7–12 | Case 2 lines, profiles, nodes, specifications | reproduced with the values of `case_04.m` (same geometry remarks as Tables 2–6; Table 10 last vertex is 3005 m, not 4005 m — errata B3) | `cases/case2_multi_well.json`, notebook 04 |
| Figs. 9–10 | Case 2 isometric view and digraph | reproduced | notebook 04 |
| Figs. 11–12 | Case 2 downhole pressure contours (nodes 7/13 and 19/25) | reproduced; state vectors equal to `results_Subsea2_1.mat` to ~10⁻⁹ at the converged points (the MATLAB sweep did not converge at 8 of the 256 grid points, the Python one at 9) | notebook 04, `figures/fig11_*`, `fig12_*` |
| Fig. 13–14 | Case 1 maximum production vs well pressure (PI = ∞) | reproduced; SLSQP optima equal to the `fmincon` optima within 10⁻³ kg/s of oil except at 260 bar, where SLSQP finds a slightly better point (16.32 kg/s with 1.71 kg/s of gas vs 16.28 kg/s with 1.60 kg/s; objective +0.2 %) | notebook 03, `figures/fig13_*`, `fig14_*` |
| Fig. 15 | Case 1 maximum production vs productivity index | reproduced (oil within 10⁻⁴ kg/s, gas within 2·10⁻³ kg/s of the MATLAB optima) | notebook 03 |
| Table 13, Fig. 16 | Case 1 maximum profit scenarios | reproduced to all printed digits when evaluated at the published decision variables; re-optimised values within 3·10⁻⁴ relative | notebook 03, `results/eval/table13*.csv` |
| Tables 14–15 | Case 2 maximum production at 220 and 240 bar | reproduced (oil and gas within 6·10⁻³ kg/s of the MATLAB optima, pressures within 0.01 bar) | notebook 05 |
| Table 16 | design parameters | input; the polytropic coefficient actually used in Eq. (31) is 1.4, not the printed 1.25 (errata A6) | notebook 05 |
| Table 17 | optimal design decision variables | reproduced within 0.012 in (diameters) and 5·10⁻³ kg/s; the objective is flat near the optimum | notebook 05, `results/eval/table17.csv` |
| Tables 18–19 | optimal design equipment and costs | reproduced to all printed digits at the published decision variables (profit 210.82 / 273.26 × 10⁶ USD/yr); the printed numbers require the program's forms of Eqs. (30)–(33) (hp, D·a) and the node-n+4 discharge pressure (errata A2–A7) | notebook 05, `results/eval/table18_19*.csv` |

Evaluation summary (`make eval`, this machine): state vectors at 14 sweep points and at the
17 published optima agree with MATLAB to max |Δy| = 3·10⁻⁹ (scaled units); Tables 13, 18 and 19
evaluated at the published decision variables agree with the paper to a maximum relative
deviation of 3·10⁻³ (rounding of the printed values); the 17 re-optimised problems reach the
MATLAB objectives within 2·10⁻³ relative (one case, see Fig. 13, is better than published).

## Differences between the paper text and the original program

The program is the source of the published numbers (the saved `.mat` files reproduce all
printed digits of Tables 13–19). The package reproduces the program by default and offers the
paper's form as an option. Items that are errors of the paper are formalised in
[`errata.md`](errata.md); the list below summarises all of them with the corresponding switch.

| # | Topic | Program (default) | Paper text | Option |
|---|---|---|---|---|
| 1 | Well geometry (Tables 2, 5, 7, 11) | wellbore and flowline 4 in; annulus 1710 m; lower wellbore 338 m; gas-line roughness 0 | 5 in; 1500 m; 300 m; 100 µm (Table 7) | `case1_single_well(wellbore_diameter_in=5, annulus_length=1500, wellbore_length=1800, gas_line_roughness=100e-6)` |
| 2 | Annulus orientation | geometry written bottom-up, so the hydrostatic term acts as for upward flow (gas-side pressures ≈ 45 bar too high) | −90°, gas flowing down | `WellDesign(annulus_orientation="physical")`; `Network.validate()` reports the inconsistency of the faithful network |
| 3 | Annulus equivalent diameter | (D_o² + D_i²)/(D_o + D_i) = 8.29 in | hydraulic diameter D_o − D_i | `WellDesign(annulus_diameter_model="hydraulic")` |
| 4 | Compressor discharge pressure P_G | node n+4 (gas-lines interconnection, −1495 m) | node n+5 (platform, +20 m) | `CostModel(discharge_pressure_node="gas_inlet")` |
| 5 | Gas pipeline cost, Eq. (33) | diameter **multiplied** by the exponent (D·a) | D^a (McCoy, 2008), as printed | `CostModel(pipeline_cost_formulation="mccoy")` |
| 6 | Polytropic coefficient of Eq. (31) | γ = 1.4 | 1.25 (Table 16) | `CompressorParameters(polytropic_exponent=1.25)` |
| 7 | Power unit in Eqs. (30), (32) | hp (w_C × 1.341) | kW | `CostModel(compressor_cost_power_unit="kW")` |
| 8 | Equation counts | 3N + 4S: 66 (Case 1), 192 (Case 2) | 66 and 232 | — (`StateLayout.size`) |
| 9 | ODE tolerance | rtol 10⁻³ for the sweeps of Figs. 6–8, 10⁻² elsewhere; atol 10⁻⁶; max step 0.5 L | "3rd to 4th order Runge–Kutta" | `IntegrationOptions(rtol=...)`; the case files carry the right value per study item |
| 10 | Kinetic terms of Eqs. (11)–(13) | \|v\| A q / (2 ρ²), ~10⁻¹⁴ in scaled units | as printed (not dimensionally a pressure) | kept (negligible) |
| 11 | "Zero" flows at the sources | 10⁻⁶ kg/s for the absent phase | 0 | `cases.EPSILON_FLOW` |
| 12 | Decision-variable bounds | 5–20 / 1–4 kg/s (Case 1); 0.2–20 / 0.2–10 (Case 2 production); 0.5–20 / 0.1–4 and 1–4 in (design) | 0 lower bounds | `study.default_study` (case files) |
| 13 | Newton tolerance | ‖F‖ ≤ 10⁻⁹ | 10⁻¹⁰ | `NewtonOptions(tolerance=...)` |

Notebook 06 quantifies the effect of items 2–7 on the gas-side quantities of Tables 13 and 18
(compressor power 5–8 % lower with the physically oriented annulus and the platform node;
profit changes by less than 0.3 %; production-side results unchanged).

## Discrepancies that are not errors of the paper

* **Optimiser.** The paper used MATLAB `fmincon` (SQP, 2009); the package uses SciPy SLSQP with
  the same bounds, constraint and starting points. The objective surfaces are flat along the
  gas rate near the optimum, so the two stop at slightly different points with practically
  equal objectives (third decimal of the rates, second decimal of the diameters). At 260 bar
  (Fig. 13) SLSQP finds a better point than the published one.
* **Grid convergence.** The MATLAB sweeps have a few non-converged cells (NaN in the `.mat`
  files); the Python sweeps converge at the same cells except for isolated points, with the
  same serpentine warm-start strategy. `gaslift run` uses a continuation from lower rates when
  a cold start fails.
* **Digraph drawings** (Figs. 4, 5, 10) are schematic; the paper's figures were drawn by hand.

## Original MATLAB code (`legacy/matlab/`)

* `bbnet10/` is the simulator: `model_network.m` (residuals), `pipe_model.m` (Beggs & Brill
  and the integrand), `rigorous_model.m` / `nonrigorous_model.m` / `simplified_model.m`
  (integration variants), `initial_calculations.m` (linear estimate), `solver_sea.m` (Newton),
  `simulate.m`, `topology.m` / `topology_gas.m` (incidence matrices and well templates),
  `objective_function_v03.m` / `objective_nonlcon_v01.m` / `objective_simulation.m`
  (optimisation coupling and economics), `optimplot.m`, `outputfunction.m`, `stopfunction.m`.
* `bbnet10/cases/offshore_single_well/` and `offshore_multi_well/` hold the case definitions
  (`case_03.m`, `case_04.m`), the sweep drivers (`subsea_test_*.m`, `subsea2_test_*.m`) and
  the optimisation drivers (`optimize_subsea_v*.m`, `optimize_subsea2_v*.m`) that produced
  Figs. 6–16 and Tables 13–19. `RESUMO_CASOS_E_ARQUIVOS.txt` maps each driver to the result
  file it saved (e.g. Table 13 scenario 1 ↔ `result_opt_subsea_v41.mat`, Tables 17–19 ↔
  `result_opt_subsea2_v101/v102.mat`, Figs. 6–8 ↔ `results_Subsea_1_2.mat` /
  `results_Subsea_2.mat`, Figs. 11–12 ↔ `results_Subsea2_1.mat`).
* The code is kept read-only as reference; it needs MATLAB with the Optimization Toolbox.

## Reference data (`data/reference/`)

| File | Content | Source |
|---|---|---|
| `case1_grid_2in.json`, `case1_grid_3in.json` | rate sweeps of Case 1 with 2 in and 3 in gas lines: grids, downhole and node-9 pressures, full scaled state vectors, incidence/specification matrices, `ode45` rtol (10⁻³) | `results_Subsea_1_2.mat`, `results_Subsea_2.mat` |
| `case2_grid.json` | rate sweep of Case 2 (16 × 16 points, rtol 10⁻²) | `results_Subsea2_1.mat` |
| `matlab_optimization_results.json` | the 17 optima used in Figs. 13–16 and Tables 13–19: decision vectors, bounds, objective, reservoir data, prices, state vector and the `saidasOpt` economics (power, discharge pressure, costs) | `result_opt_subsea_v*.mat`, `result_opt_subsea2_v*.mat` |
| `published_tables.json` | Tables 13–19 and the numerical statements of Sections 4.1 and 6.1 as printed | the paper |

The `.mat` files are exported by `scripts/export_matlab_reference.py` from the raw-material
folder outside the repository. The tests (`tests/`) compare state vectors, economics and
optima with these files; `gaslift eval` writes the full comparison.
