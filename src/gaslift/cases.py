"""Case studies of de Souza et al. (2010): single well (Case 1) and multiple wells (Case 2).

The builders reproduce the networks exactly as defined in the MATLAB scripts
``case_03.m`` and ``case_04.m`` that generated the published results. Where
those scripts differ from the values printed in the paper (wellbore diameter,
well segment lengths, gas-line roughness, annulus orientation), the script
values are used by default and the differences are documented in
``docs/reproduction.md`` (English) and ``docs/pt/paper_vs_code.md`` (Portuguese).
Use the keyword arguments to explore alternatives.

The same cases are stored as editable JSON files in ``cases/`` (see
``cases/README.md``); :func:`load_case` reads them and :func:`save_case` writes
them. ``load_case("case1")`` and ``load_case("case2")`` resolve to the files of
the paper cases.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .fluids import FluidProperties
from .network import Network, NodeRole, PipeCategory, Specifications, Well, WellDesign
from .simulation import SimulationOptions, Simulator
from .units import INCH

#: Tiny mass flow rate used instead of zero for the "absent" phase at sources
#: (keeps the Beggs & Brill correlation away from its singular single-phase limits).
EPSILON_FLOW = 1.0e-6

#: Riser outlet pressure before the choke valve, Pa (Table 6).
RISER_OUTLET_PRESSURE = 15.0e5

#: Directory with the JSON case files of the paper (``cases/`` at the repository root).
CASES_DIR = Path(__file__).resolve().parents[2] / "cases"

#: Short names of the paper cases and the JSON file of each one.
PAPER_CASE_FILES = {
    "case1": "case1_single_well.json",
    "case2": "case2_multi_well.json",
}


@dataclass
class CaseSetup:
    """A network together with its nodal specifications and well bookkeeping.

    ``study`` holds the JSON-serialisable description of the paper study attached to
    the case (operating points, rate sweeps and optimisation problems); it is
    interpreted by :mod:`gaslift.study`.
    """

    name: str
    network: Network
    specifications: Specifications
    wells: list[Well]
    paper_node_numbers: dict[str, int] = field(default_factory=dict)
    paper_pipe_numbers: dict[str, int] = field(default_factory=dict)
    description: str = ""
    study: dict = field(default_factory=dict)

    @property
    def n_wells(self) -> int:
        return len(self.wells)

    def set_well_rates(self, q_oil, q_gas) -> CaseSetup:
        """Specify produced oil (downhole) and injected gas (platform) rates, kg/s, per well."""
        q_oil = np.broadcast_to(np.asarray(q_oil, dtype=float), (self.n_wells,))
        q_gas = np.broadcast_to(np.asarray(q_gas, dtype=float), (self.n_wells,))
        for w, qo, qg in zip(self.wells, q_oil, q_gas, strict=True):
            self.specifications.set_liquid_flow(w.downhole, float(qo))
            self.specifications.set_gas_flow(w.gas_inlet, float(qg))
        return self

    def well_rates(self) -> tuple[np.ndarray, np.ndarray]:
        q_oil = np.array([self.specifications[w.downhole].liquid_flow for w in self.wells])
        q_gas = np.array([self.specifications[w.gas_inlet].gas_flow for w in self.wells])
        return q_oil, q_gas

    def set_gas_line_diameter(self, well: Well | str, diameter_m: float) -> CaseSetup:
        """Set the diameter (m) of both gas supply lines of a well."""
        w = well if isinstance(well, Well) else self.network.well(well)
        for pipe in w.gas_lines:
            self.network.set_pipe_diameter(pipe, diameter_m)
        return self

    def gas_line_length(self, well: Well | str) -> float:
        """Total length of the gas supply lines of a well, m."""
        w = well if isinstance(well, Well) else self.network.well(well)
        return float(sum(self.network.pipe(p).length for p in w.gas_lines))

    def simulator(self, options: SimulationOptions | None = None) -> Simulator:
        return Simulator(self.network, self.specifications, options)

    def copy(self) -> CaseSetup:
        """Independent copy (network geometry can then be modified freely)."""
        return CaseSetup.from_dict(self.to_dict())

    # ------------------------------------------------------------------ serialisation
    def to_dict(self) -> dict:
        """JSON-serialisable description of the case (network, specifications, study)."""
        return {
            "name": self.name,
            "description": self.description,
            "network": self.network.to_dict(),
            "specifications": self.specifications.to_dict(),
            "wells": [w.name for w in self.wells],
            "paper_node_numbers": dict(self.paper_node_numbers),
            "paper_pipe_numbers": dict(self.paper_pipe_numbers),
            "study": json.loads(json.dumps(self.study)),
        }

    @classmethod
    def from_dict(cls, data: dict) -> CaseSetup:
        net = Network.from_dict(data["network"])
        specs = Specifications(net)
        for name, s in data.get("specifications", {}).items():
            if s.get("pressure") is not None:
                specs.set_pressure(name, s["pressure"])
            if s.get("liquid_flow") is not None:
                specs.set_liquid_flow(name, s["liquid_flow"])
            if s.get("gas_flow") is not None:
                specs.set_gas_flow(name, s["gas_flow"])
        well_names = data.get("wells") or [w.name for w in net.wells]
        wells = [net.well(n) for n in well_names]
        return cls(
            data.get("name", net.name),
            net,
            specs,
            wells,
            {k: int(v) for k, v in data.get("paper_node_numbers", {}).items()},
            {k: int(v) for k, v in data.get("paper_pipe_numbers", {}).items()},
            data.get("description", ""),
            dict(data.get("study", {})),
        )


def save_case(case: CaseSetup, path: str | Path) -> Path:
    """Write a case to a JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(case.to_dict(), indent=1) + "\n")
    return path


def resolve_case_path(case: str | Path) -> Path:
    """Path of a case file: a JSON path, or a paper case name (``case1``/``case2``)."""
    name = str(case)
    if name in PAPER_CASE_FILES:
        return CASES_DIR / PAPER_CASE_FILES[name]
    return Path(case)


def load_case(case: str | Path) -> CaseSetup:
    """Load a case from a JSON file or by paper case name (``"case1"``, ``"case2"``)."""
    path = resolve_case_path(case)
    if not path.exists():
        raise FileNotFoundError(
            f"case file {path} not found (known paper cases: {', '.join(PAPER_CASE_FILES)})"
        )
    with open(path) as f:
        return CaseSetup.from_dict(json.load(f))


def paper_fluid() -> FluidProperties:
    """Fluid properties of Table 1."""
    return FluidProperties()


def _gathering_network(fluid: FluidProperties, name: str) -> Network:
    """Production line (two parts) and riser shared by both cases (Tables 2/3 and 7/8)."""
    net = Network(fluid, name)
    net.add_node("manifold_1", 0.0, 0.0, -1400.0)
    net.add_node("manifold_2", 0.0, 1500.0, -1400.0)
    net.add_node("riser_base", 0.0, 4000.0, -1500.0)
    net.add_node("riser_top", 0.0, 4000.0, 0.0, NodeRole.OUTLET)
    d6, rough = 6 * INCH, 100e-6
    net.add_pipe(
        "production_line_1",
        "manifold_1",
        "manifold_2",
        path=[(0, 0, -1400), (0, 1500, -1400)],
        diameter=d6,
        roughness=rough,
        category=PipeCategory.PRODUCTION,
    )
    net.add_pipe(
        "production_line_2",
        "manifold_2",
        "riser_base",
        path=[(0, 1500, -1400), (0, 3000, -1400), (0, 4000, -1500)],
        diameter=d6,
        roughness=rough,
        category=PipeCategory.PRODUCTION,
    )
    net.add_pipe(
        "riser",
        "riser_base",
        "riser_top",
        path=[(0, 4000, -1500), (0, 4000, 0)],
        diameter=d6,
        roughness=rough,
        category=PipeCategory.PRODUCTION,
    )
    return net


def _well_design(
    name: str,
    attach_to: str,
    wellhead: tuple[float, float, float],
    gas_line_1_path: Sequence[tuple[float, float, float]],
    *,
    gas_line_diameter: float,
    wellbore_diameter: float,
    annulus_outer_diameter: float,
    annulus_length: float,
    wellbore_length: float,
    gas_line_roughness: float,
    annulus_orientation: str,
    annulus_diameter_model: str,
) -> WellDesign:
    return WellDesign(
        name=name,
        attach_to=attach_to,
        wellhead=wellhead,
        wellbore_diameter=wellbore_diameter,
        wellbore_length=wellbore_length,
        annulus_length=annulus_length,
        annulus_outer_diameter=annulus_outer_diameter,
        gas_line_diameter=gas_line_diameter,
        gas_line_1_path=tuple(gas_line_1_path),
        gas_line_2_path=((0.0, 4000.0, 20.0), (0.0, 4000.0, -1495.0)),
        wellbore_roughness=100e-6,
        annulus_roughness=100e-6,
        gas_line_roughness=gas_line_roughness,
        annulus_orientation=annulus_orientation,
        annulus_diameter_model=annulus_diameter_model,
    )


def _apply_specifications(case: CaseSetup, q_oil, q_gas) -> None:
    specs = case.specifications
    specs.set_connection("manifold_1", "manifold_2", "riser_base")
    specs.set_pressure("riser_top", RISER_OUTLET_PRESSURE)
    for w in case.wells:
        specs.set_connection(w.wellhead, w.injection_point, w.annulus_top, w.gas_junction)
        specs.set_liquid_source(w.downhole, 0.0, EPSILON_FLOW)
        specs.set_gas_source(w.gas_inlet, 0.0, EPSILON_FLOW)
    case.set_well_rates(q_oil, q_gas)


def _number_elements(case: CaseSetup) -> None:
    case.paper_node_numbers = {n: i + 1 for i, n in enumerate(case.network.node_names)}
    case.paper_pipe_numbers = {p: i + 1 for i, p in enumerate(case.network.pipe_names)}


# Gas line 1 vertices (from the gas junction at the riser base towards the wellhead).
# The z coordinates carry the +5 m offset of the original scripts.
_GAS_LINE_1_LONG = (
    (0.0, 4000.0, -1495.0),
    (0.0, 3000.0, -1395.0),
    (0.0, 1500.0, -1395.0),
    (0.0, 0.0, -1395.0),
)
_GAS_LINE_1_SHORT = ((0.0, 4000.0, -1495.0), (0.0, 3000.0, -1395.0), (0.0, 1500.0, -1395.0))


def case1_single_well(
    q_oil: float = 10.0,
    q_gas: float = 2.0,
    *,
    gas_line_diameter_in: float = 2.0,
    wellbore_diameter_in: float = 4.0,
    annulus_outer_diameter_in: float = 10.0,
    annulus_length: float = 1710.0,
    wellbore_length: float = 2048.0,
    gas_line_roughness: float = 0.0,
    annulus_orientation: str = "matlab",
    annulus_diameter_model: str = "matlab",
    fluid: FluidProperties | None = None,
) -> CaseSetup:
    """Case 1 - offshore single well (paper Section 4.1, Tables 2-6, Figs. 3-5).

    Defaults reproduce ``case_03.m``. ``gas_line_diameter_in=3.0`` gives the
    variant of Fig. 8.
    """
    fluid = fluid or paper_fluid()
    net = _gathering_network(fluid, "Case 1 - offshore single well")
    design = _well_design(
        "well_1",
        "manifold_1",
        (400.0, -300.0, -1400.0),
        _GAS_LINE_1_LONG,
        gas_line_diameter=gas_line_diameter_in * INCH,
        wellbore_diameter=wellbore_diameter_in * INCH,
        annulus_outer_diameter=annulus_outer_diameter_in * INCH,
        annulus_length=annulus_length,
        wellbore_length=wellbore_length,
        gas_line_roughness=gas_line_roughness,
        annulus_orientation=annulus_orientation,
        annulus_diameter_model=annulus_diameter_model,
    )
    well = net.add_well(design)
    case = CaseSetup(
        "case1",
        net,
        Specifications(net),
        [well],
        description=(
            "Offshore single well with gas lift connected to a 4 km production line "
            "and a 1.5 km riser."
        ),
    )
    _apply_specifications(case, q_oil, q_gas)
    _number_elements(case)
    return case


def case2_multi_well(
    q_oil: Sequence[float] | float = (3.0, 3.0, 3.0, 3.0),
    q_gas: Sequence[float] | float = (1.0, 1.0, 1.0, 1.0),
    *,
    gas_line_diameter_in: Sequence[float] | float = 2.0,
    wellbore_diameter_in: float = 4.0,
    annulus_outer_diameter_in: float = 10.0,
    annulus_length: float = 1710.0,
    wellbore_length: float = 2048.0,
    gas_line_roughness: float = 0.0,
    annulus_orientation: str = "matlab",
    annulus_diameter_model: str = "matlab",
    fluid: FluidProperties | None = None,
) -> CaseSetup:
    """Case 2 - offshore multiple wells (paper Section 4.2, Tables 7-12, Figs. 9-10).

    Four wells: wells 1 and 2 feed manifold 1 (production line inlet), wells 3 and 4
    feed manifold 2. Defaults reproduce ``case_04.m``.
    """
    fluid = fluid or paper_fluid()
    net = _gathering_network(fluid, "Case 2 - offshore multiple wells")
    d_gas = np.broadcast_to(np.asarray(gas_line_diameter_in, dtype=float), (4,)) * INCH
    wellheads = [
        ("manifold_1", (400.0, -300.0, -1400.0), _GAS_LINE_1_LONG),
        ("manifold_1", (400.0, 300.0, -1400.0), _GAS_LINE_1_LONG),
        ("manifold_2", (400.0, 1200.0, -1400.0), _GAS_LINE_1_SHORT),
        ("manifold_2", (400.0, 1800.0, -1400.0), _GAS_LINE_1_SHORT),
    ]
    wells = []
    for i, (attach, head, gl1) in enumerate(wellheads):
        design = _well_design(
            f"well_{i + 1}",
            attach,
            head,
            gl1,
            gas_line_diameter=float(d_gas[i]),
            wellbore_diameter=wellbore_diameter_in * INCH,
            annulus_outer_diameter=annulus_outer_diameter_in * INCH,
            annulus_length=annulus_length,
            wellbore_length=wellbore_length,
            gas_line_roughness=gas_line_roughness,
            annulus_orientation=annulus_orientation,
            annulus_diameter_model=annulus_diameter_model,
        )
        wells.append(net.add_well(design))
    case = CaseSetup(
        "case2",
        net,
        Specifications(net),
        wells,
        description="Four gas-lifted wells sharing a production line and riser.",
    )
    _apply_specifications(case, q_oil, q_gas)
    _number_elements(case)
    return case
