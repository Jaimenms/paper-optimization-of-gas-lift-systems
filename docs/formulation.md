# Formulation: equation ↔ code map

Paper: de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010).
**Modeling, simulation and optimization of continuous gas lift systems for deepwater
offshore petroleum production.** *Journal of Petroleum Science and Engineering* 72 (3–4),
277–289. doi:10.1016/j.petrol.2010.03.028.

This note maps every section and equation of the paper to the modules and functions of
`gaslift` (`src/gaslift/`). SI units are used everywhere inside the package (Pa, m, kg/s);
bar and inches appear only at the API edges, in the tables and in the economic model, as in
the paper. The package reproduces the behaviour of the MATLAB program `bbnet10`
(`legacy/matlab/`) that generated the published numbers; where the paper text differs from
the program, the program is the default and the paper's variant is an option
(see [`reproduction.md`](reproduction.md) and [`errata.md`](errata.md)).

## Section 2 — Two-phase flow in a pipe (Eqs. 1–4) → `pipe_model.py`, `beggs_brill.py`, `fluids.py`

| Paper | Code |
|---|---|
| Table 1, fluid properties; polytropic gas ρ_G = ρ_G0 (P/P_0)^{1/γ}, incompressible liquid | `fluids.FluidProperties` (defaults are Table 1), `gas_density`, `gas_density_derivative` |
| Eq. (1), mixture momentum balance on the control volume of Fig. 1 | integrand assembled in `pipe_model.pressure_gradient` |
| Eq. (2), dP/dx = −(A ρ_s g sin θ + π D τ_m (1 + L_eq/L)) / (A + a_L q_L + a_G q_G) | `pipe_model.pressure_gradient` (returns dP/ds = L · dP/dx with s ∈ [0, 1]) |
| Eqs. (3)–(4), compressibility coefficients a_L, a_G with D_vsG evaluated numerically | same function: `comp = dρ_G/dP · v_sG/ρ_G`, `dhl_dvsg` by a forward difference of the hold-up |
| Beggs & Brill (1973) closure: flow pattern (Fig. 2), hold-up, friction | `beggs_brill.flow_pattern` (pattern map with logistic "soft if", as in the MATLAB code), `horizontal_holdup`, `inclination_factors`, `liquid_holdup`, `churchill_friction_factor` (Churchill, 1977), `two_phase_friction_multiplier`, `closure_relations` (τ_m, ρ_s, velocities, Froude number) |
| Pipe geometry: inclination θ, diameter D and roughness ε as functions of the axial position; equivalent length L_eq | `geometry.PipeGeometry` (3-D polyline, per-vertex diameter/roughness, `knots` with the `tinyl` transitions of the original code, `properties_at`), `pipe_model.PipeBatch` (padded arrays of all pipes) |
| Annulus treated as a pipe of equivalent diameter (footnote of Tables 2 and 7) | `geometry.annulus_equivalent_diameter` (`"matlab"`: (D_o² + D_i²)/(D_o + D_i), as the code; `"hydraulic"`: D_o − D_i) |

## Section 3.1–3.2 — Network variables and matrix notation (Eq. 5) → `network.py`

| Paper | Code |
|---|---|
| Digraph with N nodes and S lines; nodes are sources, clients, connections or fixed pressure | `network.Network` (`add_node`, `add_pipe`, `add_well`), `NodeRole`, `PipeCategory` |
| 3N + 4S variables: q_L, q_G, P_I, P_O (S each), P_W, w_L, w_G (N each); Eq. (5) y = [q_L, q_G, P_I, P_O, P_W, w_L, w_G] | `network.StateLayout` (slices `qL`, `qG`, `PI`, `PO`, `PW`, `wL`, `wG`; `size = 3N + 4S`); pressures are scaled by `ModelOptions.pressure_scale` (10⁷ Pa) |
| Incidence matrices M (±1), M1 (tails), M2 (heads) | `Network.incidence_matrix`, `start_matrix`, `end_matrix` (M1 selects the node at the beginning of each line and M2 the node at its end, as required by Eqs. (11)–(12) and as the MATLAB program does; the generation rules printed in Section 3.2 are swapped — errata B14) |
| Specification matrix E and vector e (N + S nodal specifications, at least one pressure) | `network.Specifications` (`set_pressure`, `set_liquid_flow`, `set_gas_flow`, `set_connection`, `set_liquid_source`, `set_gas_source`), `Specifications.matrix`, `validate` |
| Well sub-graph of Fig. 5 (6 nodes, 6 lines per well) | `network.WellDesign` → `Network.add_well` → `network.Well` (names of the nodes and pipes) |

## Section 3.3 — Non-linear model (Eqs. 6–14) → `model.py`

Residual vector `F(y)` assembled by `model.NetworkModel.residuals` in the order below;
`residuals_and_jacobian` adds the Jacobian (analytical blocks for the mass balances, the
nodal-momentum pressure terms and the specifications; forward finite differences, integrated
in one batch, for the derivatives of the momentum residuals with respect to q_L, q_G, P_I, P_O).

| Paper | Code |
|---|---|
| Eqs. (6)–(7), nodal mass balances M q_L + w_L = 0, M q_G + w_G = 0 (2N) | `f1 = M @ qL + wL`, `f2 = M @ qG + wG` |
| Eq. (8), pipe momentum balance P_O − P_I − ½(ΔP_1 + ΔP_2) = 0 (S) | `f3 = PO − PI − dp_mean / pressure_scale` in `NetworkModel._momentum` |
| Eqs. (9)–(10), ΔP_1 integrated from the inlet, ΔP_2 from the outlet | `pipe_model.integrate_pressure_drop` (`dp_forward`, `dp_backward`, `dp_mean`); both directions of all pipes in a single batch |
| Eqs. (11)–(12), nodal momentum balances M1ᵀ P_W − P_I + h_LI q_L + h_GI q_G = 0 and M2ᵀ P_W − P_O + h_LO q_L + h_GO q_G = 0 (2S) | `f42` (inlets) and `f41` (outlets) in `_momentum`; kinetic terms `kin_in`, `kin_out` |
| Eq. (13), h = 10⁻⁵ diag⁻¹(2 ρ ρ A) \|v\| | `kin = \|v\| A q / (2 ρ²) / pressure_scale` (kept as implemented in the MATLAB code; numerically negligible, see errata B12) |
| Eq. (14), specifications E y − e = 0 (N + S) | `f5 = E @ x − e` |
| Non-rigorous variant (gradient at mid-length, used as fall-back) | `ModelOptions.rigorous=False`, `pipe_model.midpoint_pressure_drop` |

## Section 3.4 — Linear model for the initial estimate (Eqs. 15–23) → `linear_model.py`

| Paper | Code |
|---|---|
| Eq. (15), P_I − P_O − F_L q_L − F_G q_G − d = 0 (laminar, phases independent, known hold-up H_L) | `linear_model.linear_initial_estimate` |
| Eqs. (16)–(17), F_L = 128 μ_L H_L L / (π ρ_L D⁴), F_G likewise with (1 − H_L) | `f_l`, `f_g` (divided by `pressure_scale`) |
| Eq. (18)?–(20), hydrostatic term d = (ρ_G g (1 − H_L) + ρ_L g H_L) Mᵀ z | `head` |
| Eqs. (21)–(23), linear system A y = c with the mass balances, Eq. (15), the nodal pressure identities and E y = e | matrix `a` and vector `c`, solved with `numpy.linalg.solve` |
| Hold-up assumed per pipe (0.5 for production lines and wellbores, 0 for gas lines and annuli, as in the MATLAB scripts) | `Pipe.initial_liquid_holdup` |

## Section 3.5 — Numerical methods → `integrator.py`, `solver.py`, `simulation.py`

| Paper | Code |
|---|---|
| ODE integration of Eq. (2) ("3rd to 4th order Runge–Kutta"; the program used MATLAB `ode45`) | `integrator.integrate_batch`: Dormand–Prince 5(4) with the step-size control of `ode45` (relative error against max(\|y\|, atol/rtol), growth/shrink factors 1.25/0.8, rejection factors 0.1 and 0.5, last step stretched up to 10 %), many independent problems with their own adaptive step |
| Tolerances (`odeset`: rtol 10⁻² or 10⁻³, atol 10⁻⁶, max step 0.5 L) | `pipe_model.IntegrationOptions` (defaults rtol 1e-2; the sweeps of Figs. 6–8 use 1e-3, see `study.default_study`) |
| Newton–Raphson with numerical Jacobian, ‖F‖ ≤ 10⁻¹⁰ (10⁻⁹ in the program) | `solver.newton_solve` (`NewtonOptions.tolerance = 1e-9`, bound clipping of the iterates), Jacobian in `NetworkModel.residuals_and_jacobian` (coarse perturbation 10⁻³ far from the solution, 10⁻⁵ near it) |
| Solution strategy: rigorous model from the linear estimate; non-rigorous model then rigorous refinement on failure | `simulation.Simulator.solve` / `_solve_from`, `simulation.simulate`; warm start from the previous solution for sweeps and optimisations |
| Results: nodal pressures, pipe flows, profiles along the pipes | `results.NetworkSolution` (`nodes_table`, `pipes_table`, `wells_table`, `profile`), `results.PipeProfiles` (`pipe_model.profile_states`) |

## Section 4 — Case studies (Tables 2–12, Figs. 3–12) → `cases.py`, `cases/*.json`, `sweeps.py`, `plotting.py`, `tables.py`

| Paper | Code |
|---|---|
| Case 1 (Section 4.1): Tables 2–6, 10 nodes / 9 lines, 66 equations | `cases.case1_single_well` = `cases/case1_single_well.json` (`cases.load_case("case1")`); `tables.pipes_table`, `profile_table`, `nodes_table`, `specifications_table` |
| Case 2 (Section 4.2): Tables 7–12, 28 nodes / 27 lines, 192 equations (paper: 232) | `cases.case2_multi_well` = `cases/case2_multi_well.json` |
| Node and line numbering of the paper | `CaseSetup.paper_node_numbers`, `paper_pipe_numbers` |
| Figs. 3 and 9, isometric projections | `plotting.plot_case_topology`, `plot_topology_3d` |
| Figs. 4, 5 and 10, digraphs | `plotting.plot_digraph` |
| Figs. 6–8 and 11–12, contour plots over injected gas × produced oil | `sweeps.sweep_rates` (serpentine traversal with warm start, as the original scripts), `plotting.plot_rate_contours`; grids defined in the `study.sweeps` block of the case files |
| "Gas lift inlet pressure" of Figs. 7–8 (node n+4, gas-lines interconnection, as in the program; the text says node n+5) | `RateSweep.gas_junction` (program) and `RateSweep.gas_inlet` (platform node) |

## Section 5 — Gas-lift optimisation (Eqs. 24–35) → `optimization.py`, `economics.py`

| Paper | Code |
|---|---|
| Eq. (24), maximise Ω_j subject to LB ≤ u ≤ UB and (P_well − P_downhole) − diag⁻¹(PI) Q_L = 0 | `optimization.GasLiftProblem` (`bounds`, `constraint_value` = P_res − P_downhole − Q_L/PI per well, `solve` with SciPy SLSQP, the SQP of the paper was MATLAB `fmincon`) |
| Decision variables Q_G, Q_L (and D_G for the design problem) | `GasLiftProblem.split`, `initial_vector`; `Objective.OPTIMAL_DESIGN` adds the gas-line diameters (`CaseSetup.set_gas_line_diameter`) |
| Reservoir model P_well, PI (PI = ∞ ⇒ P_downhole = P_well) | `optimization.Reservoir` |
| Eqs. (25)–(27), Ω_1 = R, Ω_2 = R − AGC, Ω_3 = R − AGC − ATCC | `optimization.Objective`, cost switches set by `GasLiftProblem._cost_for_objective` (`CostModel.include_fuel_gas_cost`, `include_capital_cost`); `EconomicEvaluation.profit` / `objective = −profit` |
| Eq. (28), ISBL = C_C + C_T + C_P; Eq. (29), ATCC = CCF · OSBL, CCF = 1/3, OSBL = 3 ISBL | `economics.evaluate_economics` (`isbl`, `osbl`, `annualized_capital_cost`; `CostModel.capital_charge_factor`, `osbl_factor`) |
| Eq. (30), C_C = 14.3·10⁻³ w_C^0.62 (Couper et al., 2005) | `economics.compressor_capital_cost` (power in hp as in the program, `CostModel.compressor_cost_power_unit="kW"` for the printed form) |
| Eq. (31), compression power of an n-stage machine with intercooling | `economics.compressor_power`, `CompressorParameters` (Table 16; γ = 1.4 as in the program, 1.25 printed) |
| Eq. (32), C_T = 7.19·10⁻³ w_C^0.81 | `economics.turbine_capital_cost` |
| Eq. (33), gas pipeline cost (McCoy, 2008) | `economics.pipeline_capital_cost` (`"matlab"`: D·a as in the program; `"mccoy"`: D^a as printed) |
| Eq. (34), AGC = 8.45·10⁻³ w_C c_GC | `economics.fuel_gas_volume` × `CostModel.gas_price` |
| Eq. (35), R = 198.34 Σ Q_L / ρ_L0 · c_L | `economics.oil_volume_rate` × `CostModel.oil_price` (198.34 = seconds per year / barrel / 10⁶) |
| Section 5.6, data exchanged between simulator and optimiser (Q_L, Q_G in; P_downhole, P_G out) | `GasLiftProblem.simulate` (sets the well rates, runs `Simulator`, evaluates `evaluate_economics`); P_G read at node n+4 (`CostModel.discharge_pressure_node="gas_junction"`, program) or n+5 (`"gas_inlet"`, text) |
| Section 5.7, SQP with bounds (the actual bounds of the runs are in the case files) | `study.default_study` → `cases/*.json`, `study.build_problem` |

## Section 6 — Optimisation results (Figs. 13–16, Tables 13–19) → `study.py`, `tables.py`, `figures.py`, `evaluation.py`

| Paper | Code |
|---|---|
| Section 6.1.1, Figs. 13–14 (PI = ∞, P_well 180–260 bar) | problems `fig13_*bar` of `cases/case1_single_well.json`; `tables.fig13_table`; `figures.export_figures` |
| Section 6.1.2, Fig. 15 (220 bar, PI 0.5–10) | problems `fig15_pi*`; `tables.fig15_table` |
| Section 6.1.3, Table 13 and Fig. 16 (maximum profit, oil price scenarios) | problems `table13_scenario*`; `tables.table13` |
| Section 6.2.1, Tables 14–15 (Case 2 maximum production) | problems `table14_maxprod_220bar`, `table15_maxprod_240bar`; `tables.table14_15` |
| Section 6.2.2, Tables 16–19 (optimal design) | problems `table17_19_design_*bar` (started from the maximum-production optimum with 2 in lines, as the paper); `tables.table17`, `table18`, `table19` |
| Published values of Tables 13–19 | `data/reference/published_tables.json`; MATLAB optima in `data/reference/matlab_optimization_results.json`; comparison by `evaluation.evaluate` (`gaslift eval`) |
