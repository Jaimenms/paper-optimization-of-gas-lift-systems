# Case files

The two case studies of de Souza et al. (2010) as editable JSON files, read by the
`gaslift` command line (`validate`, `run`, `figures`) and by `gaslift.cases.load_case`:

| File | Paper | Content |
|---|---|---|
| `case1_single_well.json` | Section 4.1 / 6.1, Tables 2-6, Figs. 3-8 and 13-16, Table 13 | one gas-lifted well, 4 km production line and 1.5 km riser (10 nodes, 9 lines) |
| `case2_multi_well.json` | Section 4.2 / 6.2, Tables 7-12 and 14-19, Figs. 9-12 | four wells sharing the production line and riser (28 nodes, 27 lines) |

The files are generated from the Python builders (`gaslift.cases.case1_single_well`,
`case2_multi_well`) and the study definitions (`gaslift.study.default_study`) by
`uv run python scripts/export_cases.py`; the tests check that file and builder agree.
The values follow the MATLAB scripts that produced the published numbers
(`legacy/matlab/bbnet10/cases/*/case_0{3,4}.m`), not the tables of the paper where the two
differ (see `docs/reproduction.md`).

## Fields

```jsonc
{
  "name": "case1",                       // short name (used for results/<name>/)
  "description": "...",
  "network": {
    "name": "...",
    "fluid": { "liquid_density": 850.0, "liquid_viscosity": 0.012, "gas_viscosity": 8e-5,
               "gas_reference_density": 80.0, "reference_pressure": 1e7,
               "gas_polytropic_coefficient": 1.25, "interfacial_tension": 0.071,
               "gas_normal_density": 1.2, "gas_molecular_weight": 25.0 },   // Table 1, SI units
    "nodes": [ { "name": "manifold_1", "x": 0.0, "y": 0.0, "z": -1400.0, "role": "connection" }, ... ],
    "pipes": [ { "name": "production_line_1", "start": "manifold_1", "end": "manifold_2",
                 "category": "production", "initial_liquid_holdup": 0.5,
                 "geometry": { "points": [[x, y, z], ...],      // vertices, inlet -> outlet, m
                               "diameters": [...], "roughness": [...],   // per vertex, m
                               "equivalent_length": 0.0, "transition_length": 1.0 } }, ... ],
    "wells": [ { "name": "well_1", "wellhead": "well_1:wellhead", ..., "design": { ... } } ]
  },
  "specifications": { "riser_top": { "pressure": 1500000.0 },           // Pa
                      "well_1:downhole": { "liquid_flow": 10.0, "gas_flow": 1e-06 },   // kg/s
                      "well_1:gas_inlet": { "liquid_flow": 1e-06, "gas_flow": 2.0 }, ... },
  "wells": ["well_1"],
  "paper_node_numbers": { "manifold_1": 1, ... },     // numbering of Tables 5 / 11
  "paper_pipe_numbers": { "production_line_1": 1, ... },
  "study": { ... }                                     // see below
}
```

* **Nodes** carry coordinates (m) and a role (`connection`, `outlet`, `wellhead`,
  `injection_point`, `downhole`, `annulus_top`, `gas_junction`, `gas_inlet`) used for
  reporting and plotting.
* **Pipes** are 3-D polylines with a diameter and roughness per vertex; the inclination of
  each segment follows from the vertices. `category` selects the plotting style
  (`production`, `flowline`, `wellbore`, `annulus`, `gas`). `initial_liquid_holdup` is the
  hold-up assumed by the linear initial-estimate model (Section 3.4).
* **Wells** record the six nodes and six pipes created from a `WellDesign` (Fig. 5) and keep
  the design (`wellbore_diameter`, `annulus_length`, `annulus_orientation`, ...). Editing a
  well geometry means editing its pipes/nodes (the `design` block is informative once the
  pipes exist).
* **Specifications** give, per node, the fixed pressure (Pa) and/or the external liquid and gas
  mass flow rates (kg/s). The model needs exactly N + S specified values and at least one
  pressure. Sources specify 1e-6 kg/s instead of 0 for the absent phase, as the original code.
* **Study** describes what `gaslift run` and `gaslift figures` do with the case:

```jsonc
"study": {
  "integration_rtol": 1e-2,              // default ode45 relative tolerance (Section 3.5)
  "operating_points": [                  // steady-state simulations (gaslift run)
    { "name": "matlab_grid_2in_10_3", "paper": "Fig. 6 / Fig. 7 sweep point",
      "q_oil": 10.0, "q_gas": 2.0, "rtol": 1e-3, "gas_line_diameter_in": 2.0,
      "reference": { "grid": "case1_grid_2in.json", "index": [10, 3] } }   // MATLAB state vector
  ],
  "sweeps": [                            // rate sweeps of Figs. 6-8 and 11-12 (gaslift figures)
    { "name": "fig06_fig07_2in", "paper": ["Fig. 6", "Fig. 7"],
      "q_gas": [0.5, ..., 4.0], "q_oil": [0.5, 1.0, ..., 20.0], "rtol": 1e-3 }
  ],
  "optimizations": [                     // optimisation problems (gaslift run)
    { "name": "table13_scenario1", "paper": ["Table 13", "Fig. 16"],
      "objective": "max_profit",         // max_production | max_profit | optimal_design
      "cost": { "oil_price": 12.7875 },  // CostModel fields (gaslift.economics)
      "reservoirs": [ { "pressure": 220.0, "productivity_index": 1.0 } ],   // bar, (kg/s)/bar; null = infinite
      "oil_bounds": [5.0, 20.0], "gas_bounds": [1.0, 4.0], "diameter_bounds": [1.0, 4.0],
      "x0": { "q_oil": 10.0, "q_gas": 2.0 },      // or {"from": "<problem name>", "diameters_in": 2.0}
      "reference": "case1_maxprofit_oil12.79",    // key of data/reference/matlab_optimization_results.json
      "published": { "table": "table13", "column": "1" } }
  ]
}
```

To study a different network, copy a file, edit nodes, pipes, specifications and the study
block, and run `uv run gaslift validate my_case.json` and `uv run gaslift run my_case.json`.
Networks can also be assembled in Python (`Network.add_node`, `add_pipe`, `add_well`) and
written with `gaslift.cases.save_case`.
