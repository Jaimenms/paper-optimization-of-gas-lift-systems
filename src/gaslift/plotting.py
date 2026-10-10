"""Matplotlib helpers reproducing the figures of the paper.

All functions accept an optional ``ax`` and return the axes they drew on, so
they compose into multi-panel figures inside the notebooks.
"""

from __future__ import annotations

from collections.abc import Sequence

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from .beggs_brill import flow_pattern
from .network import Network, NodeRole, PipeCategory
from .results import NetworkSolution
from .units import BAR

PRODUCTION_CATEGORIES = (PipeCategory.PRODUCTION, PipeCategory.FLOWLINE, PipeCategory.WELLBORE)
GAS_CATEGORIES = (PipeCategory.GAS, PipeCategory.ANNULUS)

_STYLE = {
    PipeCategory.PRODUCTION: dict(color="tab:blue", ls="-", lw=3, label="Oil/gas lines"),
    PipeCategory.FLOWLINE: dict(color="tab:blue", ls="-", lw=3, label="Oil/gas lines"),
    PipeCategory.WELLBORE: dict(color="tab:blue", ls="-.", lw=3, label="Wellbores"),
    PipeCategory.ANNULUS: dict(color="tab:red", ls=":", lw=3, label="Annuli"),
    PipeCategory.GAS: dict(color="tab:red", ls="--", lw=3, label="Gas lines"),
    PipeCategory.OTHER: dict(color="0.4", ls="-", lw=2, label="Other"),
}


# --------------------------------------------------------------------------- contour plots
def plot_rate_contours(
    q_gas: np.ndarray,
    q_oil: np.ndarray,
    values: np.ndarray,
    ax: Axes | None = None,
    levels: Sequence[float] | int | None = None,
    *,
    title: str = "",
    xlabel: str = "Injected gas mass flow rate (kg/s)",
    ylabel: str = "Produced oil mass flow rate (kg/s)",
    fmt: str = "%g",
    label_pad: float = 0.0,
) -> Axes:
    """Labelled black contour lines of ``values[i_oil, j_gas]`` (paper Figs. 6-8, 11-12)."""
    ax = ax or plt.gca()
    qg, qo = np.meshgrid(q_gas, q_oil)
    z = np.ma.masked_invalid(np.asarray(values, dtype=float))
    cs = ax.contour(qg, qo, z, levels=levels, colors="k", linewidths=1.5)
    ax.clabel(cs, inline=True, fmt=fmt, fontsize=9, inline_spacing=3 + label_pad)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return ax


# --------------------------------------------------------------------------- 3-D topology
def _box(ax, x0, x1, y0, y1, z0, z1, color, alpha):
    """Translucent box faces (sea water / sea bed)."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    faces = [
        [(x0, y0, z1), (x0, y1, z1), (x1, y1, z1), (x1, y0, z1)],
        [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
        [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
    ]
    coll = Poly3DCollection(faces, facecolor=color, edgecolor="0.3", linewidths=0.5, alpha=alpha)
    ax.add_collection3d(coll)


def _bounds(network: Network, seafloor: float, sea_level: float, pad: float = 0.05):
    pts = np.vstack([p.geometry.points for p in network.pipes])
    x0, x1 = pts[:, 0].min(), pts[:, 0].max()
    y0, y1 = pts[:, 1].min(), pts[:, 1].max()
    z0, z1 = min(pts[:, 2].min(), seafloor), max(pts[:, 2].max(), sea_level)
    dx, dy, dz = max(x1 - x0, 1.0), max(y1 - y0, 1.0), max(z1 - z0, 1.0)
    return (
        x0 - pad * dx,
        x1 + pad * dx,
        y0 - pad * dy,
        y1 + pad * dy,
        z0 - pad * dz,
        z1 + pad * dz,
    )


def plot_topology_3d(
    network: Network,
    ax: Axes | None = None,
    *,
    categories: Sequence[PipeCategory] | None = None,
    nodes: Sequence[str] | None = None,
    seafloor: float = -1500.0,
    sea_level: float = 0.0,
    view: tuple[float, float] = (28.0, -50.0),
    legend: bool = True,
    bounds: tuple[float, float, float, float, float, float] | None = None,
) -> Axes:
    """Isometric view of the pipes (paper Figs. 3 and 9) with sea and seabed volumes."""
    if ax is None:
        ax = plt.figure(figsize=(7, 6)).add_subplot(projection="3d")
    cats = tuple(categories) if categories is not None else tuple(PipeCategory)
    drawn: set[PipeCategory] = set()
    for p in network.pipes:
        if p.category not in cats:
            continue
        st = _STYLE[p.category]
        pts = p.geometry.points
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=st["color"], ls=st["ls"], lw=st["lw"])
        drawn.add(p.category)
    node_list = [network.node(n) for n in nodes] if nodes is not None else list(network.nodes)
    if node_list:
        ax.scatter(
            [n.x for n in node_list],
            [n.y for n in node_list],
            [n.z for n in node_list],
            color="k",
            s=40,
            depthshade=False,
        )
    x0, x1, y0, y1, z0, z1 = bounds if bounds is not None else _bounds(network, seafloor, sea_level)
    _box(ax, x0, x1, y0, y1, z0, seafloor, "tan", 0.25)
    _box(ax, x0, x1, y0, y1, seafloor, sea_level, "cyan", 0.15)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_zlim(z0, z1)
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    ax.set_box_aspect((max(dx, 0.3 * dy), dy, dz))
    ax.set_xticks([round(x0), round(x1)])
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.view_init(elev=view[0], azim=view[1])
    if legend:
        handles, seen = [], set()
        for c in PipeCategory:
            st = _STYLE[c]
            if c in drawn and st["label"] not in seen:
                handles.append(
                    Line2D([], [], color=st["color"], ls=st["ls"], lw=st["lw"], label=st["label"])
                )
                seen.add(st["label"])
        handles.append(Line2D([], [], color="k", marker="o", ls="", label="Nodes"))
        ax.legend(handles=handles, loc="upper left", fontsize=8)
    return ax


def plot_case_topology(network: Network, figsize=(15, 6.5)) -> Figure:
    """Two-panel isometric view: production system (left) and gas supply system (right)."""
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    prod_roles = {
        NodeRole.CONNECTION,
        NodeRole.OUTLET,
        NodeRole.WELLHEAD,
        NodeRole.INJECTION_POINT,
        NodeRole.DOWNHOLE,
    }
    gas_roles = {
        NodeRole.INJECTION_POINT,
        NodeRole.ANNULUS_TOP,
        NodeRole.GAS_JUNCTION,
        NodeRole.GAS_INLET,
    }
    prod_nodes = [n.name for n in network.nodes if n.role in prod_roles]
    gas_nodes = [n.name for n in network.nodes if n.role in gas_roles]
    bounds = _bounds(network, -1500.0, 0.0)
    plot_topology_3d(
        network, ax1, categories=PRODUCTION_CATEGORIES, nodes=prod_nodes, bounds=bounds
    )
    plot_topology_3d(network, ax2, categories=GAS_CATEGORIES, nodes=gas_nodes, bounds=bounds)
    ax1.set_title("Production system")
    ax2.set_title("Gas-lift supply system")
    fig.subplots_adjust(wspace=0.02, left=0.02, right=0.98)
    return fig


# --------------------------------------------------------------------------- digraph
def _schematic_positions(network: Network) -> dict[str, tuple[float, float]]:
    """Node positions for a convergent gathering network with well sub-graphs.

    Non-well nodes are placed on a horizontal line ordered by their distance (in
    pipes) to the outlet, with mostly vertical pipes (risers) drawn upwards. Each
    well is drawn below its attachment node as in Fig. 5 of the paper: wellhead,
    injection point and downhole in a column, the annulus top and the gas supply
    chain to the right.
    """
    well_nodes = {n for w in network.wells for n in w.nodes}
    trunk = [n for n in network.nodes if n.name not in well_nodes]
    m = network.incidence_matrix()
    # distance to outlet, following the (unique) outgoing pipe of each node
    rank: dict[str, int] = {}
    for n in trunk:
        i = network.node_index(n.name)
        steps, visited = 0, {i}
        while True:
            out = np.flatnonzero(m[i, :] < 0)
            if out.size == 0:
                break
            head = int(np.flatnonzero(m[:, out[0]] > 0)[0])
            if head in visited:
                break
            visited.add(head)
            i, steps = head, steps + 1
        rank[n.name] = steps
    max_rank = max(rank.values(), default=0)
    pos: dict[str, tuple[float, float]] = {}
    spacing = 3.0
    for n in trunk:
        pos[n.name] = ((max_rank - rank[n.name]) * spacing, 0.0)
    # risers: draw the head of a mostly vertical pipe above its tail
    for p in network.pipes:
        if p.start in pos and p.end in pos:
            v = p.geometry.end - p.geometry.start
            if abs(v[2]) > np.hypot(v[0], v[1]):
                pos[p.end] = (pos[p.start][0], pos[p.start][1] + 1.3)
    per_attach: dict[str, int] = {}
    for w in network.wells:
        attach = w.design.attach_to
        j = per_attach.get(attach, 0)
        per_attach[attach] = j + 1
        xa = pos[attach][0] if attach in pos else 0.0
        xw = xa - 0.9 + 1.4 * j
        pos[w.wellhead] = (xw, -0.7)
        pos[w.injection_point] = (xw, -1.6)
        pos[w.downhole] = (xw, -2.3)
        pos[w.annulus_top] = (xw + 0.55, -1.0)
        pos[w.gas_junction] = (xw + 0.95, -1.55)
        pos[w.gas_inlet] = (xw + 0.95, -2.3)
    return pos


def plot_digraph(
    network: Network,
    ax: Axes | None = None,
    *,
    node_numbers: dict[str, int] | None = None,
    pipe_numbers: dict[str, int] | None = None,
    highlight_wells: bool = True,
    layout: str = "schematic",
) -> Axes:
    """Digraph of the network (paper Figs. 4, 5 and 10).

    ``layout="schematic"`` draws wells as compact sub-graphs below the gathering
    line; ``layout="coordinates"`` uses the ``y`` coordinate and elevation of the
    nodes.
    """
    ax = ax or plt.gca()
    g = nx.DiGraph()
    for n in network.nodes:
        g.add_node(n.name)
    for p in network.pipes:
        g.add_edge(p.start, p.end, name=p.name, category=p.category)
    if layout == "schematic":
        pos = _schematic_positions(network)
    elif layout == "coordinates":
        pos = {n.name: (n.y, n.z) for n in network.nodes}
    else:
        raise ValueError(f"unknown layout {layout!r}")
    labels = {
        n: str(node_numbers.get(n, "")) if node_numbers else n.split(":")[-1] for n in g.nodes
    }
    colors = ["white"] * len(g.nodes)
    if highlight_wells:
        well_nodes = {nn for w in network.wells for nn in w.nodes}
        colors = ["lightyellow" if n in well_nodes else "white" for n in g.nodes]
    nx.draw_networkx_nodes(g, pos, ax=ax, node_size=450, node_color=colors, edgecolors="k")
    nx.draw_networkx_labels(g, pos, labels=labels, ax=ax, font_size=8)
    for cat in PipeCategory:
        edges = [(u, v) for u, v, d in g.edges(data=True) if d["category"] == cat]
        if not edges:
            continue
        st = _STYLE[cat]
        nx.draw_networkx_edges(
            g,
            pos,
            edgelist=edges,
            ax=ax,
            edge_color=st["color"],
            style=st["ls"],
            width=1.8,
            arrows=True,
            arrowsize=14,
            node_size=450,
        )
    if pipe_numbers:
        edge_labels = {(p.start, p.end): f"({pipe_numbers[p.name]})" for p in network.pipes}
        nx.draw_networkx_edge_labels(
            g, pos, edge_labels=edge_labels, ax=ax, font_size=7, rotate=False
        )
    handles = []
    for cat in (
        PipeCategory.PRODUCTION,
        PipeCategory.WELLBORE,
        PipeCategory.ANNULUS,
        PipeCategory.GAS,
    ):
        if any(p.category == cat for p in network.pipes):
            st = _STYLE[cat]
            handles.append(Line2D([], [], color=st["color"], ls=st["ls"], lw=2, label=st["label"]))
    if handles:
        ax.legend(handles=handles, loc="upper left", fontsize=8)
    ax.set_axis_off()
    ax.margins(0.08)
    return ax


# --------------------------------------------------------------------------- flow map
def plot_flow_pattern_map(ax: Axes | None = None, n: int = 400) -> Axes:
    """Horizontal flow pattern map of Beggs & Brill (paper Fig. 2)."""
    ax = ax or plt.gca()
    lam = np.logspace(-3, 0, n)
    fr = np.logspace(-1, 2, n)
    lam_g, fr_g = np.meshgrid(lam, fr)
    weights, _ = flow_pattern(lam_g, fr_g)
    idx = np.rint(weights.index).astype(int)
    cmap = ListedColormap(["white", "0.85", "0.6", "0.35"])
    mesh = ax.pcolormesh(lam_g, fr_g, idx, cmap=cmap, vmin=0.5, vmax=4.5, shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("No-slip liquid hold-up")
    ax.set_ylabel("Liquid Froude number")
    ax.set_title("Horizontal flow pattern map")
    cb = plt.colorbar(mesh, ax=ax, ticks=[1, 2, 3, 4])
    cb.ax.set_yticklabels(["Segregated", "Transition", "Intermittent", "Distributed"])
    return ax


# --------------------------------------------------------------------------- profiles
def plot_pressure_path(
    solution: NetworkSolution,
    pipes: Sequence[str],
    ax: Axes | None = None,
    *,
    label: str | None = None,
    reverse: bool = False,
) -> Axes:
    """Pressure along a chain of pipes (requires profiles)."""
    ax = ax or plt.gca()
    if solution.profiles is None:
        raise ValueError("solution has no profiles")
    offset = 0.0
    xs, ps = [], []
    for name in pipes:
        j = solution.network.pipe_index(name)
        s = solution.profiles.axial_position[j]
        p = solution.profiles.pressure[j] / BAR
        if reverse:
            s, p = s[::-1], p[::-1]
            s = s[0] - s
        xs.append(offset + s)
        ps.append(p)
        offset += s[-1]
    ax.plot(np.concatenate(xs), np.concatenate(ps), label=label)
    ax.set_xlabel("Distance along path (m)")
    ax.set_ylabel("Pressure (bar)")
    return ax


def plot_optimization_history(history, ax: Axes | None = None) -> Axes:
    """Objective value of every simulation evaluated during an optimisation."""
    ax = ax or plt.gca()
    obj = np.array([r.objective for r in history])
    conv = np.array([r.converged for r in history])
    k = np.arange(len(obj))
    ax.plot(k[conv], -obj[conv], "k.-", label="profit (converged)")
    if (~conv).any():
        ax.plot(k[~conv], np.zeros((~conv).sum()), "rx", label="failed simulation")
    ax.set_xlabel("Evaluation")
    ax.set_ylabel("Objective (10⁶ USD/yr)")
    ax.legend()
    return ax
