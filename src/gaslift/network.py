"""Pipe-network topology (digraph), well templates and nodal specifications.

The network is a directed graph with ``N`` nodes and ``S`` pipes (lines). Nodes
carry the external mass flow rates (sources/clients) and the stagnant pressure;
pipes carry internal mass flow rates and inlet/outlet pressures. The state
vector of the two-phase model is ``y = [qL, qG, PI, PO, PW, wL, wG]`` (paper
Eq. 5), which :class:`StateLayout` indexes.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from enum import StrEnum

import numpy as np

from .fluids import FluidProperties
from .geometry import PipeGeometry, Point, annulus_equivalent_diameter


class PipeCategory(StrEnum):
    """Role of a pipe in a production system (used for reporting and plotting)."""

    PRODUCTION = "production"
    FLOWLINE = "flowline"
    WELLBORE = "wellbore"
    ANNULUS = "annulus"
    GAS = "gas"
    OTHER = "other"


class NodeRole(StrEnum):
    """Role of a node in a production system (used for reporting and plotting)."""

    CONNECTION = "connection"
    OUTLET = "outlet"
    WELLHEAD = "wellhead"
    INJECTION_POINT = "injection_point"
    DOWNHOLE = "downhole"
    ANNULUS_TOP = "annulus_top"
    GAS_JUNCTION = "gas_junction"
    GAS_INLET = "gas_inlet"


@dataclass(frozen=True)
class Node:
    name: str
    x: float
    y: float
    z: float
    role: NodeRole = NodeRole.CONNECTION

    @property
    def position(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


@dataclass(frozen=True)
class Pipe:
    name: str
    start: str
    """Name of the node at the beginning of the pipe (digraph tail)."""
    end: str
    """Name of the node at the end of the pipe (digraph head)."""
    geometry: PipeGeometry
    category: PipeCategory = PipeCategory.OTHER
    initial_liquid_holdup: float = 0.5
    """Liquid hold-up assumed by the linear model that generates the initial estimate."""

    @property
    def length(self) -> float:
        return self.geometry.length


# --------------------------------------------------------------------------- wells
@dataclass(frozen=True)
class WellDesign:
    """Template describing a gas-lifted well and its gas supply lines.

    A well expands into six nodes and six pipes (paper Fig. 5): the flowline from
    the wellhead to the gathering network, the upper and lower wellbore sections
    (above and below the gas injection point), the casing-tubing annulus that
    carries the lift gas down to the injection point, and two gas supply lines
    from the platform to the wellhead.

    All lengths are in metres and diameters in metres.
    """

    name: str
    attach_to: str
    """Name of the gathering-network node that receives the well flowline."""
    wellhead: Point
    wellbore_diameter: float
    wellbore_length: float
    """Developed length from the wellhead down to the reservoir (downhole node)."""
    annulus_length: float
    """Developed length from the wellhead down to the gas injection point."""
    annulus_outer_diameter: float
    gas_line_diameter: float
    gas_line_1_path: Sequence[Point]
    """Vertices of the subsea gas line, from the gas junction towards the wellhead
    (the wellhead itself is appended automatically)."""
    gas_line_2_path: Sequence[Point]
    """Vertices of the gas line from the platform (gas inlet) to the gas junction."""
    wellbore_roughness: float = 100e-6
    annulus_roughness: float = 100e-6
    gas_line_roughness: float = 0.0
    annulus_orientation: str = "matlab"
    """``"matlab"`` writes the annulus geometry bottom-up (as the original code did,
    which makes its hydrostatic term act as in upward flow); ``"physical"`` writes it
    top-down, consistent with the digraph orientation and gas flowing downwards."""
    annulus_diameter_model: str = "matlab"
    """See :func:`gaslift.geometry.annulus_equivalent_diameter`."""

    @property
    def annulus_equivalent_diameter(self) -> float:
        return annulus_equivalent_diameter(
            self.annulus_outer_diameter, self.wellbore_diameter, self.annulus_diameter_model
        )


@dataclass(frozen=True)
class Well:
    """Names of the nodes and pipes created for a well (see :class:`WellDesign`)."""

    name: str
    design: WellDesign
    wellhead: str
    injection_point: str
    downhole: str
    annulus_top: str
    gas_junction: str
    gas_inlet: str
    flowline: str
    wellbore_upper: str
    wellbore_lower: str
    annulus: str
    gas_line_1: str
    gas_line_2: str

    @property
    def nodes(self) -> tuple[str, ...]:
        return (
            self.wellhead,
            self.injection_point,
            self.downhole,
            self.annulus_top,
            self.gas_junction,
            self.gas_inlet,
        )

    @property
    def pipes(self) -> tuple[str, ...]:
        return (
            self.flowline,
            self.wellbore_upper,
            self.wellbore_lower,
            self.annulus,
            self.gas_line_1,
            self.gas_line_2,
        )

    @property
    def gas_lines(self) -> tuple[str, str]:
        return (self.gas_line_1, self.gas_line_2)


# --------------------------------------------------------------------------- state layout
@dataclass(frozen=True)
class StateLayout:
    """Index bookkeeping for the state vector ``[qL, qG, PI, PO, PW, wL, wG]``."""

    n_nodes: int
    n_pipes: int

    @property
    def size(self) -> int:
        return 3 * self.n_nodes + 4 * self.n_pipes

    @property
    def n_specifications(self) -> int:
        return self.n_nodes + self.n_pipes

    @property
    def qL(self) -> slice:
        return slice(0, self.n_pipes)

    @property
    def qG(self) -> slice:
        return slice(self.n_pipes, 2 * self.n_pipes)

    @property
    def PI(self) -> slice:
        return slice(2 * self.n_pipes, 3 * self.n_pipes)

    @property
    def PO(self) -> slice:
        return slice(3 * self.n_pipes, 4 * self.n_pipes)

    @property
    def PW(self) -> slice:
        return slice(4 * self.n_pipes, 4 * self.n_pipes + self.n_nodes)

    @property
    def wL(self) -> slice:
        s = 4 * self.n_pipes + self.n_nodes
        return slice(s, s + self.n_nodes)

    @property
    def wG(self) -> slice:
        s = 4 * self.n_pipes + 2 * self.n_nodes
        return slice(s, s + self.n_nodes)

    def index_PW(self, node: int) -> int:
        return self.PW.start + node

    def index_wL(self, node: int) -> int:
        return self.wL.start + node

    def index_wG(self, node: int) -> int:
        return self.wG.start + node


# --------------------------------------------------------------------------- network
class Network:
    """Directed pipe network with named nodes and pipes.

    Nodes and pipes keep their insertion order, which defines the numbering used
    by the incidence matrices and the state vector.
    """

    def __init__(self, fluid: FluidProperties | None = None, name: str = "") -> None:
        self.name = name
        self.fluid = fluid if fluid is not None else FluidProperties()
        self._nodes: list[Node] = []
        self._pipes: list[Pipe] = []
        self._node_index: dict[str, int] = {}
        self._pipe_index: dict[str, int] = {}
        self._wells: list[Well] = []

    # ------------------------------------------------------------------ properties
    @property
    def nodes(self) -> tuple[Node, ...]:
        return tuple(self._nodes)

    @property
    def pipes(self) -> tuple[Pipe, ...]:
        return tuple(self._pipes)

    @property
    def wells(self) -> tuple[Well, ...]:
        return tuple(self._wells)

    @property
    def n_nodes(self) -> int:
        return len(self._nodes)

    @property
    def n_pipes(self) -> int:
        return len(self._pipes)

    @property
    def layout(self) -> StateLayout:
        return StateLayout(self.n_nodes, self.n_pipes)

    @property
    def node_names(self) -> list[str]:
        return [n.name for n in self._nodes]

    @property
    def pipe_names(self) -> list[str]:
        return [p.name for p in self._pipes]

    def node(self, name: str) -> Node:
        return self._nodes[self.node_index(name)]

    def pipe(self, name: str) -> Pipe:
        return self._pipes[self.pipe_index(name)]

    def well(self, name: str) -> Well:
        for w in self._wells:
            if w.name == name:
                return w
        raise KeyError(f"unknown well {name!r}")

    def node_index(self, name: str) -> int:
        try:
            return self._node_index[name]
        except KeyError as exc:
            raise KeyError(f"unknown node {name!r}") from exc

    def pipe_index(self, name: str) -> int:
        try:
            return self._pipe_index[name]
        except KeyError as exc:
            raise KeyError(f"unknown pipe {name!r}") from exc

    @property
    def elevations(self) -> np.ndarray:
        return np.array([n.z for n in self._nodes], dtype=float)

    @property
    def pipe_lengths(self) -> np.ndarray:
        return np.array([p.length for p in self._pipes], dtype=float)

    # ------------------------------------------------------------------ construction
    def add_node(
        self, name: str, x: float, y: float, z: float, role: NodeRole = NodeRole.CONNECTION
    ) -> Node:
        if name in self._node_index:
            raise ValueError(f"node {name!r} already exists")
        node = Node(name, float(x), float(y), float(z), role)
        self._node_index[name] = len(self._nodes)
        self._nodes.append(node)
        return node

    def add_pipe(
        self,
        name: str,
        start: str,
        end: str,
        geometry: PipeGeometry | None = None,
        *,
        path: Iterable[Point] | None = None,
        diameter: float | None = None,
        roughness: float = 0.0,
        equivalent_length: float = 0.0,
        category: PipeCategory = PipeCategory.OTHER,
        initial_liquid_holdup: float = 0.5,
    ) -> Pipe:
        """Add a pipe from node ``start`` to node ``end``.

        Either a ready :class:`PipeGeometry` or ``path`` + ``diameter`` must be given.
        When ``path`` is omitted the pipe is a straight segment between the nodes.
        """
        if name in self._pipe_index:
            raise ValueError(f"pipe {name!r} already exists")
        n_start, n_end = self.node(start), self.node(end)
        if geometry is None:
            if diameter is None:
                raise ValueError("diameter is required when no geometry is given")
            pts = list(path) if path is not None else [n_start.position, n_end.position]
            geometry = PipeGeometry.from_points(pts, diameter, roughness, equivalent_length)
        pipe = Pipe(name, start, end, geometry, category, float(initial_liquid_holdup))
        self._pipe_index[name] = len(self._pipes)
        self._pipes.append(pipe)
        return pipe

    def add_well(self, design: WellDesign) -> Well:
        """Expand a :class:`WellDesign` into its six nodes and six pipes."""
        w = design.name
        attach = self.node(design.attach_to)
        wx, wy, wz = (float(v) for v in design.wellhead)
        gl1 = [tuple(float(v) for v in p) for p in design.gas_line_1_path]
        gl2 = [tuple(float(v) for v in p) for p in design.gas_line_2_path]
        if not gl1 or not gl2:
            raise ValueError("gas line paths must have at least one vertex")

        names = {
            "wellhead": f"{w}:wellhead",
            "injection_point": f"{w}:injection_point",
            "downhole": f"{w}:downhole",
            "annulus_top": f"{w}:annulus_top",
            "gas_junction": f"{w}:gas_junction",
            "gas_inlet": f"{w}:gas_inlet",
        }
        z_inj = wz - design.annulus_length
        z_down = wz - design.wellbore_length
        self.add_node(names["wellhead"], wx, wy, wz, NodeRole.WELLHEAD)
        self.add_node(names["injection_point"], wx, wy, z_inj, NodeRole.INJECTION_POINT)
        self.add_node(names["downhole"], wx, wy, z_down, NodeRole.DOWNHOLE)
        self.add_node(names["annulus_top"], wx, wy, wz, NodeRole.ANNULUS_TOP)
        self.add_node(names["gas_junction"], *gl1[0], NodeRole.GAS_JUNCTION)
        self.add_node(names["gas_inlet"], *gl2[0], NodeRole.GAS_INLET)

        wb_d, wb_r = design.wellbore_diameter, design.wellbore_roughness
        self.add_pipe(
            f"{w}:flowline",
            names["wellhead"],
            attach.name,
            path=[(wx, wy, wz), tuple(attach.position)],
            diameter=wb_d,
            roughness=wb_r,
            category=PipeCategory.FLOWLINE,
            initial_liquid_holdup=0.5,
        )
        self.add_pipe(
            f"{w}:wellbore_upper",
            names["injection_point"],
            names["wellhead"],
            path=[(wx, wy, z_inj), (wx, wy, wz)],
            diameter=wb_d,
            roughness=wb_r,
            category=PipeCategory.WELLBORE,
            initial_liquid_holdup=0.5,
        )
        self.add_pipe(
            f"{w}:wellbore_lower",
            names["downhole"],
            names["injection_point"],
            path=[(wx, wy, z_down), (wx, wy, z_inj)],
            diameter=wb_d,
            roughness=wb_r,
            category=PipeCategory.WELLBORE,
            initial_liquid_holdup=0.5,
        )
        if design.annulus_orientation == "matlab":
            annulus_path = [(wx, wy, z_inj), (wx, wy, wz)]
        elif design.annulus_orientation == "physical":
            annulus_path = [(wx, wy, wz), (wx, wy, z_inj)]
        else:
            raise ValueError(f"unknown annulus orientation {design.annulus_orientation!r}")
        self.add_pipe(
            f"{w}:annulus",
            names["annulus_top"],
            names["injection_point"],
            path=annulus_path,
            diameter=design.annulus_equivalent_diameter,
            roughness=design.annulus_roughness,
            category=PipeCategory.ANNULUS,
            initial_liquid_holdup=0.0,
        )
        self.add_pipe(
            f"{w}:gas_line_1",
            names["gas_junction"],
            names["annulus_top"],
            path=[*gl1, (wx, wy, wz)],
            diameter=design.gas_line_diameter,
            roughness=design.gas_line_roughness,
            category=PipeCategory.GAS,
            initial_liquid_holdup=0.0,
        )
        self.add_pipe(
            f"{w}:gas_line_2",
            names["gas_inlet"],
            names["gas_junction"],
            path=gl2,
            diameter=design.gas_line_diameter,
            roughness=design.gas_line_roughness,
            category=PipeCategory.GAS,
            initial_liquid_holdup=0.0,
        )
        well = Well(
            name=w,
            design=design,
            flowline=f"{w}:flowline",
            wellbore_upper=f"{w}:wellbore_upper",
            wellbore_lower=f"{w}:wellbore_lower",
            annulus=f"{w}:annulus",
            gas_line_1=f"{w}:gas_line_1",
            gas_line_2=f"{w}:gas_line_2",
            **names,
        )
        self._wells.append(well)
        return well

    def set_pipe_geometry(self, name: str, geometry: PipeGeometry) -> None:
        """Replace the geometry of an existing pipe (used for design variables)."""
        i = self.pipe_index(name)
        self._pipes[i] = replace(self._pipes[i], geometry=geometry)

    def set_pipe_diameter(self, name: str, diameter: float) -> None:
        """Give a pipe a uniform diameter (m)."""
        self.set_pipe_geometry(name, self.pipe(name).geometry.with_uniform_diameter(diameter))

    # ------------------------------------------------------------------ matrices
    def incidence_matrix(self) -> np.ndarray:
        """Incidence matrix ``M`` (N x S): -1 at the pipe tail node, +1 at its head node."""
        m = np.zeros((self.n_nodes, self.n_pipes))
        for j, p in enumerate(self._pipes):
            m[self.node_index(p.start), j] = -1.0
            m[self.node_index(p.end), j] = 1.0
        return m

    def start_matrix(self) -> np.ndarray:
        """``M1`` (N x S): 1 where a pipe starts at a node."""
        return (self.incidence_matrix() < 0).astype(float)

    def end_matrix(self) -> np.ndarray:
        """``M2`` (N x S): 1 where a pipe ends at a node."""
        return (self.incidence_matrix() > 0).astype(float)

    def downstream_path(self, pipe: str) -> list[str]:
        """Pipes traversed from ``pipe`` to the network outlet (tree networks only)."""
        m = self.incidence_matrix()
        j = self.pipe_index(pipe)
        path = [j]
        while True:
            head = np.flatnonzero(m[:, j] > 0)
            if head.size == 0:
                break
            nxt = np.flatnonzero(m[head[0], :] < 0)
            if nxt.size == 0:
                break
            j = int(nxt[0])
            path.append(j)
        return [self._pipes[k].name for k in path]

    # ------------------------------------------------------------------ validation
    def validate(self) -> list[str]:
        """Return human-readable warnings about the network definition."""
        issues: list[str] = []
        m = self.incidence_matrix()
        for i, node in enumerate(self._nodes):
            if np.count_nonzero(m[i, :] < 0) > 1:
                issues.append(
                    f"node {node.name!r} is a divergent node (two or more pipes leave it); "
                    "the nodal model of the paper requires convergent tree networks"
                )
        for p in self._pipes:
            tail, head = self.node(p.start), self.node(p.end)
            if not np.allclose(p.geometry.start, tail.position, atol=1e-6):
                issues.append(
                    f"pipe {p.name!r}: geometry starts at {tuple(p.geometry.start)} but its tail "
                    f"node {p.start!r} is at {tuple(tail.position)}"
                )
            if not np.allclose(p.geometry.end, head.position, atol=1e-6):
                issues.append(
                    f"pipe {p.name!r}: geometry ends at {tuple(p.geometry.end)} but its head "
                    f"node {p.end!r} is at {tuple(head.position)}"
                )
        if math.isclose(np.linalg.matrix_rank(m), 0):
            issues.append("network has no pipes")
        return issues

    # ------------------------------------------------------------------ serialisation
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "fluid": self.fluid.to_dict(),
            "nodes": [
                {"name": n.name, "x": n.x, "y": n.y, "z": n.z, "role": n.role.value}
                for n in self._nodes
            ],
            "pipes": [
                {
                    "name": p.name,
                    "start": p.start,
                    "end": p.end,
                    "category": p.category.value,
                    "initial_liquid_holdup": p.initial_liquid_holdup,
                    "geometry": p.geometry.to_dict(),
                }
                for p in self._pipes
            ],
            "wells": [
                {k: getattr(w, k) for k in w.__dataclass_fields__ if k != "design"}
                | {"design": {**w.design.__dict__}}
                for w in self._wells
            ],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Network:
        net = cls(FluidProperties.from_dict(data.get("fluid", {})), data.get("name", ""))
        for n in data["nodes"]:
            net.add_node(n["name"], n["x"], n["y"], n["z"], NodeRole(n.get("role", "connection")))
        for p in data["pipes"]:
            net.add_pipe(
                p["name"],
                p["start"],
                p["end"],
                PipeGeometry.from_dict(p["geometry"]),
                category=PipeCategory(p.get("category", "other")),
                initial_liquid_holdup=p.get("initial_liquid_holdup", 0.5),
            )
        for w in data.get("wells", []):
            design = WellDesign(**w["design"])
            net._wells.append(Well(design=design, **{k: v for k, v in w.items() if k != "design"}))
        return net

    def __repr__(self) -> str:
        return (
            f"Network(name={self.name!r}, nodes={self.n_nodes}, pipes={self.n_pipes}, "
            f"wells={len(self._wells)})"
        )


# --------------------------------------------------------------------------- specifications
@dataclass
class NodeSpecification:
    """Specified nodal variables; ``None`` means free (computed by the model)."""

    pressure: float | None = None
    """Stagnant pressure, Pa."""
    liquid_flow: float | None = None
    """External liquid mass flow rate, kg/s (positive entering the network)."""
    gas_flow: float | None = None
    """External gas mass flow rate, kg/s (positive entering the network)."""

    @property
    def count(self) -> int:
        return sum(v is not None for v in (self.pressure, self.liquid_flow, self.gas_flow))


class Specifications:
    """Nodal specifications of a network model (paper Section 3.2, matrix ``E`` and vector ``e``).

    A well-posed model needs exactly ``N + S`` specified nodal variables, including
    at least one pressure.
    """

    def __init__(self, network: Network) -> None:
        self._network = network
        self._specs: dict[str, NodeSpecification] = {
            n.name: NodeSpecification() for n in network.nodes
        }

    def copy(self) -> Specifications:
        new = Specifications(self._network)
        new._specs = {
            k: NodeSpecification(v.pressure, v.liquid_flow, v.gas_flow)
            for k, v in self._specs.items()
        }
        return new

    def __getitem__(self, node: str) -> NodeSpecification:
        return self._specs[node]

    def items(self):
        return self._specs.items()

    # ------------------------------------------------------------------ setters
    def set_pressure(self, node: str, pressure: float) -> Specifications:
        self._specs[node].pressure = float(pressure)
        return self

    def set_liquid_flow(self, node: str, flow: float) -> Specifications:
        self._specs[node].liquid_flow = float(flow)
        return self

    def set_gas_flow(self, node: str, flow: float) -> Specifications:
        self._specs[node].gas_flow = float(flow)
        return self

    def set_connection(self, *nodes: str) -> Specifications:
        """Mark nodes as pure connections (zero external liquid and gas flows)."""
        for node in nodes:
            self.set_liquid_flow(node, 0.0).set_gas_flow(node, 0.0)
        return self

    def set_liquid_source(
        self, node: str, liquid_flow: float, gas_flow: float = 0.0
    ) -> Specifications:
        return self.set_liquid_flow(node, liquid_flow).set_gas_flow(node, gas_flow)

    def set_gas_source(
        self, node: str, gas_flow: float, liquid_flow: float = 0.0
    ) -> Specifications:
        return self.set_liquid_flow(node, liquid_flow).set_gas_flow(node, gas_flow)

    def clear(self, node: str) -> Specifications:
        self._specs[node] = NodeSpecification()
        return self

    # ------------------------------------------------------------------ queries
    @property
    def count(self) -> int:
        return sum(s.count for s in self._specs.values())

    def validate(self) -> None:
        layout = self._network.layout
        if self.count != layout.n_specifications:
            raise ValueError(
                f"the model needs N + S = {layout.n_specifications} specifications, "
                f"got {self.count}"
            )
        if not any(s.pressure is not None for s in self._specs.values()):
            raise ValueError("at least one nodal pressure must be specified")

    def matrix(self, pressure_scale: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
        """Specification matrix ``E`` ((N+S) x (3N+4S)) and vector ``e``.

        Pressures are divided by ``pressure_scale`` to match the scaled state vector.
        """
        self.validate()
        layout = self._network.layout
        rows: list[tuple[int, float]] = []
        for i, node in enumerate(self._network.nodes):
            s = self._specs[node.name]
            if s.pressure is not None:
                rows.append((layout.index_PW(i), s.pressure / pressure_scale))
            if s.liquid_flow is not None:
                rows.append((layout.index_wL(i), s.liquid_flow))
            if s.gas_flow is not None:
                rows.append((layout.index_wG(i), s.gas_flow))
        e_mat = np.zeros((len(rows), layout.size))
        e_vec = np.zeros(len(rows))
        for r, (col, value) in enumerate(rows):
            e_mat[r, col] = 1.0
            e_vec[r] = value
        return e_mat, e_vec

    def to_dict(self) -> dict:
        return {
            name: {k: v for k, v in vars(s).items() if v is not None}
            for name, s in self._specs.items()
            if s.count
        }

    def __repr__(self) -> str:
        return f"Specifications({self.count} of {self._network.layout.n_specifications})"
