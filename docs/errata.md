# Errata

Formal erratum (draft) of the paper reimplemented in this repository, prepared from the
`gaslift` re-implementation. Every corrected number is reproduced by the tests and by
`gaslift eval`; the reproduction status per table and figure is in
[`reproduction.md`](reproduction.md) and the Portuguese working notes in
[`pt/paper_vs_code.md`](pt/paper_vs_code.md).

## Erratum to: de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010). Modeling, simulation and optimization of continuous gas lift systems for deepwater offshore petroleum production. *Journal of Petroleum Science and Engineering* 72 (3–4), 277–289. doi:10.1016/j.petrol.2010.03.028

## Scope and how locations are cited

While porting to Python the MATLAB program (`bbnet10`) that produced every number in the
paper, the published text and tables were compared, item by item, with the code and with
the saved result files (`.mat`). The saved results match all printed digits of Tables 13–19
and Figs. 6–8 and 11–16, so the program is unambiguously the source of the published
numbers. The items below are the places where the paper's text or tables differ from what
was actually computed, plus textual inconsistencies found on the way.

Locations are given as **section · page · line numbers of the typeset proof**
(Elsevier proof PETROL 1857, 13 pages; the page numbers of the published article are the
proof page plus 276, i.e. proof p. 4 = p. 280). Tables, figures and equations keep the same
numbers in the proof and in the published version.

Part A lists corrections that change reported numbers or the description of the model that
generated them. Part B lists textual and tabular inconsistencies with no effect on the
numbers. Part C gives reproducibility notes. Corrected values were recomputed with the
Python port (`gaslift`, repository notebooks 02–06 and `gaslift eval`), which reproduces the original
program to machine precision; where "corrected" values are given they are **evaluated at
the published decision variables**, not re-optimised.

---

## Part A — Corrections affecting the model description or the reported numbers

### A1. Well geometry used in the computations differs from Tables 2, 5, 7 and 11

**Location.** Section 4.1, Table 2 (p. 4, lines t2:7–t2:10) and Table 5 (p. 4, lines t5:9–t5:10); Section 4.2, Table 7 (p. 7, lines t7:7–t7:12) and Table 11 (p. 8, lines t11:12–t11:13).

**Published.** Wellbore lines (n, n+1) and well-to-production line (n−1): diameter 5 in; wellbore part 2: 1500 m; wellbore part 1: 300 m; gas-lift annulus: 1500 m; footnote: "annulus with internal diameter 5 in and external diameter 10 in". Node heights: annulus/wellbore interconnection −2900 m, downhole −3200 m. Table 7: roughness 100 µm for the annulus and the gas lines (n+2, n+3, n+4).

**Computed with.** Wellbore and flowline diameter **4 in**; wellbore part 2 (wellhead to injection point) **1710 m**; wellbore part 1 (injection point to reservoir) **338 m** (total 2048 m); annulus length **1710 m**; annulus inner diameter 4 in; node heights **−3110 m** (interconnection) and **−3448 m** (downhole); gas-line roughness **0 µm** in both cases (annulus 100 µm).

**Impact.** All simulation and optimisation results (Figs. 6–8, 11–16, Tables 13–19) correspond to the computed geometry. With the geometry printed in Table 2 the downhole pressure of Fig. 6 at 10 kg/s of oil and 2 kg/s of gas would be 178 bar instead of 214 bar. The tables should be corrected to the values above.

### A2. Orientation of the gas-lift annulus (Table 2, line 7; Table 7, line n+2)

**Location.** Table 2 (p. 4, line t2:10): "Gas lift annulus … 1500 … −90°"; Table 7 (p. 7, line t7:10).

**Issue.** The digraph orientation of the annulus (node n+3 at the wellhead to node n+1 at the injection point, i.e. gas flowing downwards, −90°) is correct, but the geometry fed to the integrator was written from the injection point up to the wellhead (+90°). The hydrostatic term of Eq. (2) was therefore evaluated as for **upward** flow along the annulus, and the pressure at the top of the annulus (node n+3) is **higher** than at the injection point by the gas column instead of lower.

**Impact.** Only the gas-side pressures are affected (the production side is specified by mass flow rates and does not feed back): the "gas lift inlet pressure" of Figs. 7–8, the discharge pressure P_G and compressor power of Table 13, the compressor and cost figures of Tables 18–19 and the fuel-gas cost in Fig. 16. At the operating point of Table 13 ("Max" scenario) the pressure at node 9 is 248.3 bar with the physically oriented annulus instead of 279.7 bar. See the recomputed values in A8.

### A3. Equivalent diameter of the annulus (footnote of Tables 2 and 7)

**Location.** Table 2 footnote a (p. 4, line t2:13); Table 7 footnote a (p. 7, line t7:13).

**Published.** "Hydraulic diameters: annulus with internal diameter 5 in and external diameter 10 in."

**Computed with.** D_eq = 4(A_outer + A_inner)/(P_outer + P_inner) = (D_o² + D_i²)/(D_o + D_i) = 8.29 in for D_o = 10 in and D_i = 4 in. The hydraulic diameter of an annulus is D_o − D_i = 6 in (5 in with the diameters printed in the footnote). The effect on the results is small (annulus friction is minor); the footnote should state the expression actually used, or the results should be recomputed with D_o − D_i.

### A4. Node at which the "compressor discharge" / "gas lift inlet pressure" is reported

**Location.** Section 4.1 (p. 7, lines 434–439): "the pressure at the compressor discharge or the pressure at the gas inlet (node 10) … presented in Fig. 7 … Fig. 8"; captions of Figs. 7 and 8 (pp. 5–6); Section 5.6 (p. 11, lines 605–606 and 613–616): "pressure at gas sources (injected gas pressure at the platform, P_G)"; Table 13 (p. 10, line t13:13) "Discharge pressure (P_G)"; Table 18 (p. 11, line t18:5).

**Issue.** The quantity plotted in Figs. 7–8 and used as P_G in Eq. (31) is the pressure at node **n+4** (node 9 in Case 1: interconnection of the gas lines at the base of the riser, elevation −1495 m), not at node **n+5** (node 10: gas inlet at the platform, +20 m). The two differ by the hydrostatic head of the 1515 m vertical gas line, 18 bar at the conditions of Table 13 (e.g. 261.9 bar at node 10 versus 279.7 bar at node 9 in the "Max" scenario). The text and captions should read "pressure at the gas-lines interconnection (node 9 / n+4)", or the results should be recomputed with the platform node.

### A5. Gas pipeline cost, Eq. (33)

**Location.** Section 5.4, Eq. (33) (p. 10, lines 554–563); Table 18 (p. 11, line t18:8) "Gas pipelines capital cost (C_P)".

**Issue.** Eq. (33) is printed correctly (diameter raised to the exponents 1.590, 0.940 and 0.791, as in McCoy, 2008), but the code multiplied the diameter by the exponent (D·a instead of D^a). The values of Table 18 (1.05 and 0.98 × 10⁶ USD) are reproduced only with the product form. With Eq. (33) as printed, C_P = **1.12** and **1.05** × 10⁶ USD for scenarios 1 and 2 (ISBL 6.85 and 6.53 × 10⁶ USD).

### A6. Polytropic coefficient used in the compressor power, Eq. (31)

**Location.** Section 5.2, Eq. (31) (p. 10, lines 538–547); Table 16 (p. 11, line t16:6): "Fuel gas polytrophic coefficient (γ) 1.25".

**Issue.** The compressor powers of Tables 13 and 18 are reproduced only with γ = **1.4**, not 1.25 as listed in Table 16 (γ = 1.25 is the coefficient of the gas density law of Table 1, used in the flow model). With γ = 1.25 the powers would be 944.3, 990.3, 1015.6 and 1042.9 kW in Table 13 (scenarios 1, 2, 3, Max) and 1696.4 and 1593.7 kW in Table 18. Table 16 should list γ = 1.4 for Eq. (31), or the results should be recomputed.

### A7. Units of the power in the cost correlations, Eqs. (30) and (32)

**Location.** Section 5.2 (p. 10, lines 535–537) and Section 5.3 (p. 10, lines 551–553): "w_C is the power of the machine in kW".

**Issue.** The code converts the power to brake horsepower (w_C × 1.341) before applying the Couper et al. (2005) correlations, and the values of Table 18 correspond to hp: C_C = 14.3·10⁻³ (1.341 w_C)^0.62 and C_T = 7.19·10⁻³ (1.341 w_C)^0.81. Evaluated with w_C in kW as printed, C_C = 1.49 and 1.43 × 10⁶ USD and C_T = 3.11 and 2.96 × 10⁶ USD (scenarios 1 and 2) instead of 1.79/1.72 and 3.95/3.76. The text should state hp.

### A8. Recomputed gas-side quantities (A2 + A3 + A4 together)

Evaluated at the published decision variables with the physically oriented annulus, hydraulic diameter D_o − D_i and discharge pressure at the platform node (hp-based cost correlations and the product form of Eq. (33) retained, so that only the gas-side items change):

| Table 13 scenario | P_G published (bar) | P_G corrected (bar) | w_C published (kW) | w_C corrected (kW) | AGC published → corrected (10⁶ USD/yr) |
|---|---|---|---|---|---|
| 1 | 269.18 | 221.85 | 1002.7 | 944.8 | 1.60 → 1.51 |
| 2 | 274.03 | 227.36 | 1051.8 | 993.3 | 1.68 → 1.58 |
| 3 | 276.74 | 230.43 | 1078.9 | 1020.2 | 1.72 → 1.63 |
| Max | 279.68 | 233.74 | 1108.1 | 1049.1 | 1.77 → 1.67 |

| Table 18/19 | Scenario 1 published | Scenario 1 corrected | Scenario 2 published | Scenario 2 corrected |
|---|---|---|---|---|
| Compressor discharge pressure (bar) | 235.44 | 183.38 | 258.02 | 202.15 |
| Compressor power (kW) | 1797.5 | 1659.8 | 1691.1 | 1567.0 |
| C_C (10⁶ USD) | 1.79 | 1.70 | 1.72 | 1.64 |
| C_T (10⁶ USD) | 3.95 | 3.70 | 3.76 | 3.53 |
| ISBL (10⁶ USD) | 6.79 | 6.45 | 6.46 | 6.16 |
| AGC (10⁶ USD/yr) | 2.86 | 2.65 | 2.69 | 2.50 |
| Profit (10⁶ USD/yr) | 210.82 | 211.37 | 273.26 | 273.76 |

The compression work depends weakly on the pressure ratio (exponent (γ−1)/(nγ) ≈ 0.1), so the
discharge-pressure corrections of 47–56 bar change the compressor power by 5–8 % and the
profit by less than 0.3 %. Downhole pressures, produced oil, injected gas, revenue and the
conclusions of Sections 6 and 7 are unchanged; the optimal decision variables would move
only marginally if re-optimised with the corrected model.

---

## Part B — Textual and tabular inconsistencies (no effect on the numbers)

### B1. Equation counts

**Location.** Section 4.1 (p. 7, lines 410–414); Section 4.2 (p. 8, lines 459–463).

**Published.** Case 1: "66 equations: 9 mass balance equations for liquid phase at lines, 9 mass balance equations for gas phase at lines, 9 momentum balance equations at pipes, 20 momentum balance equations at nodes and 19 specification equations." Case 2: "232 equations: 27 + 27 + 27 + 56 + 95".

**Correction.** The model of Section 3.3 has 2N + S + 2S + (N + S) = 3N + 4S equations. Case 1 (N = 10, S = 9): 66 = 10 + 10 + 9 + 18 + 19. Case 2 (N = 28, S = 27): **192** = 28 + 28 + 27 + 54 + 55 (the count of 232 and the decomposition 27/27/27/56/95 are not consistent with the formulation).

### B2. Length and cross-references of line 2 and line 8 in Tables 2 and 7

**Location.** Table 2 (p. 4, lines t2:5 and t2:11); Table 7 (p. 7, lines t7:5 and t7:11).

**Published.** Line 2 "Production line—part 2 … 4005 … See Table 5"; line 8 (n+3) "Gas line—part 1 … 4505 … See Table 6".

**Correction.** The length of line 2 is **2505 m** (4005 m is the distance from node 1, see Table 3). The inclination profiles are in Table **3** (line 2) and Table **4** (line 8) for Case 1, and Tables **8**, **9** (lines 8 and 14) and **10** (lines 20 and 26) for Case 2; Tables 5 and 6 are the nodal tables. The length of lines 20 and 26 is **3005 m** (see B3).

### B3. Table 10 (profile of lines 20 and 26)

**Location.** Table 10 (p. 8, lines t10:3–t10:6).

**Published.** Length 0, 1005, 2505, **4005** m; height −1495, −1395, −1395, **−1395** m.

**Correction.** The last vertex is the wellhead: length **3005** m and height **−1400** m (consistent with Table 16, L_G3 = L_G4 = 3005 + 1515 = 4520 m).

### B4. Table 17, third diameter row

**Location.** Table 17 (p. 11, line t17:14): "Pipe diameter (D_G2)—lines 20 and 21".

**Correction.** D_G**3**.

### B5. Well node numbers of Case 2

**Location.** Section 4.2 (p. 8, line 445): "nodes 5, 11, 17 and 13 correspond to wells".

**Correction.** Nodes 5, 11, 17 and **23**.

### B6. Node and pipe counts of Case 1

**Location.** Section 4.1 (p. 7, lines 394–398): "The topology of this network is composed by 9 nodes and 8 pipes … the fifth node corresponds to a well … contains 6 nodes and 5 pipes".

**Correction.** The complete network of Tables 2 and 5 has **10 nodes and 9 lines**; the condensed digraph of Fig. 4 has 5 nodes (4 network nodes plus the well) and 4 lines; the well sub-graph of Fig. 5 has 6 nodes and 6 lines (lines 4–9).

### B7. Equation cross-reference in Section 6.2.2

**Location.** Section 6.2.2 (p. 12, lines 697–698): "Economic parameters used in Eqs. (21), (23) and (24)".

**Correction.** Eqs. (21) and (23) belong to the linear initial-estimate model; the economic parameters enter Eqs. **(28)–(35)** (and Table 16).

### B8. "Different wellbore diameters" in the Abstract and Section 4

**Location.** Abstract (p. 1, lines 35–37); Section 4 (p. 7, lines 385–387): "as a function of the injected gas flow rate for different reservoir pressures and different wellbore diameters".

**Correction.** The diameter varied in Figs. 7–8 is that of the **gas supply lines** (2 in and 3 in); the wellbore diameter is fixed.

### B9. Table 4, x-coordinate row

**Location.** Table 4 (p. 4, line t4:5): "x-coordinate (m) 00004 0 0" (typesetting).

**Correction.** 0, 0, 0, 0, **400** (the last vertex is the wellhead at x = 400 m, y = −300 m).

### B10. Numerical methods statement

**Location.** Section 3.5 (p. 7, lines 377–381).

**Published.** "a 3rd to 4th order Runge–Kutta type method … convergence criterion was the norm of residue vector equal to 10⁻¹⁰".

**Correction.** The ODE solver used was MATLAB `ode45`, the Dormand–Prince **4(5)** pair, with relative tolerance 10⁻³ for the rate sweeps of Figs. 6–8 and 10⁻² for all other results, absolute tolerance 10⁻⁶ Pa and maximum step 0.5 of the pipe length; the Newton–Raphson stopping criterion was ‖F‖ ≤ **10⁻⁹**.

### B11. Bounds of the decision variables

**Location.** Section 5.7 (pp. 11–12, lines 619–623): "(iii) minimum gas injection equal to 0 kg/s and (iv) minimum oil production equal to 0 kg/s".

**Correction.** The runs used 5 ≤ Q_L ≤ 20 kg/s and 1 ≤ Q_G ≤ 4 kg/s (Case 1), 0.2 ≤ Q_L ≤ 20 and 0.2 ≤ Q_G ≤ 10 kg/s (Case 2, maximum production) and 0.5 ≤ Q_L ≤ 20, 0.1 ≤ Q_G ≤ 4 kg/s (Case 2, design). None of these bounds is active at the reported optima.

### B12. Kinetic terms of the nodal momentum balances, Eqs. (11)–(13)

**Location.** Section 3.3 (p. 5, lines 299–315).

**Note.** As printed and as implemented, h = 10⁻⁵ diag⁻¹(2 ρ × ρ × A)|v| multiplied by q is not dimensionally a pressure; its magnitude is of order 10⁻¹⁴ in the scaled pressure unit (10⁷ Pa), so the nodal balances reduce in practice to P_W = P_I and P_W = P_O. The term can be dropped without changing any result.

### B13. Zero flow specifications at the source nodes (Tables 6 and 12)

**Location.** Table 6 (p. 4, lines t6:6–t6:7); Table 12 (p. 8, lines t12:6–t12:7).

**Note.** The code specifies 10⁻⁶ kg/s, not 0, for the gas flow at the oil inlet and for the liquid flow at the gas inlet, to keep the Beggs and Brill correlation away from its single-phase limits. The effect is negligible.

### B14. Generation rules of the matrices M1 and M2

**Location.** Section 3.2 (p. 5, lines 247–250), rules (4) and (5); Eqs. (11)–(12).

**Published.** "M1 has elements equal to one at the positions where M is equal to one" (nodes at the *end* of a line) and "M2 … where M is equal to minus one" (nodes at the *beginning*).

**Correction.** In Eqs. (11)–(12) M1ᵀ P_W is compared with the inlet pressures P_I and M2ᵀ P_W with the outlet pressures P_O, so M1 must select the node at the **beginning** of each line (M = −1) and M2 the node at its **end** (M = +1); the program (`topology.m`) builds them this way. Rules (4) and (5) should be swapped.

---

## Part C — Reproducibility notes

* Results were generated with MATLAB program `bbnet10` and case scripts `case_03.m` (Case 1)
  and `case_04.m` (Case 2); the saved result files reproduce all printed digits of
  Tables 13–19 (e.g. Table 13 scenario 1 ↔ `result_opt_subsea_v41.mat`, Tables 17–19 ↔
  `result_opt_subsea2_v101/v102.mat`, Figs. 6–8 ↔ `results_Subsea_1_2.mat` and
  `results_Subsea_2.mat`, Figs. 11–12 ↔ `results_Subsea2_1.mat`).
* The Python re-implementation (`gaslift`) reproduces the MATLAB state vectors of all these
  files to about 10⁻¹⁰ (relative to the scaled unknowns) when run with the settings of
  B10, and the economic quantities of Tables 13, 18 and 19 to all printed digits. Its
  optimiser (SLSQP) recovers the published optima of Tables 13–15 to the third decimal and
  those of Table 17 to the second (the objective is flat near the optimum).
* Items already addressed in the authors' proof corrections (citation of Fig. 1 instead of
  Fig. 5 at line 147; the incomplete sentence at lines 441–443; keyword list; uncited
  references) are not repeated here.
