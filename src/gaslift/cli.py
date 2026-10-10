"""Command line interface: ``gaslift <verb> ...`` (see ``docs/conventions.md``).

Verbs: ``info``, ``validate``, ``run``, ``eval``, ``figures``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .cases import CASES_DIR, PAPER_CASE_FILES, load_case

PAPER = (
    "de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010). "
    "Modeling, simulation and optimization of continuous gas lift systems for deepwater "
    "offshore petroleum production. Journal of Petroleum Science and Engineering 72 (3-4), "
    "277-289. doi:10.1016/j.petrol.2010.03.028"
)


def _print(msg: str) -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------- verbs
def cmd_info(args: argparse.Namespace) -> int:
    _print(f"gaslift {__version__} - two-phase pipe-network simulation and gas-lift optimisation")
    _print("")
    _print("Paper:")
    _print(f"  {PAPER}")
    _print("")
    _print("Cases:")
    for name, file in PAPER_CASE_FILES.items():
        path = CASES_DIR / file
        if path.exists():
            case = load_case(path)
            lay = case.network.layout
            study = case.study
            root = CASES_DIR.parent
            shown = path.relative_to(root) if path.is_relative_to(root) else path
            _print(f"  {name:6s} {shown}  - {case.description}")
            _print(
                f"         {lay.n_nodes} nodes, {lay.n_pipes} pipes, {case.n_wells} well(s); "
                f"{len(study.get('operating_points', []))} operating points, "
                f"{len(study.get('sweeps', []))} sweeps, "
                f"{len(study.get('optimizations', []))} optimisation problems"
            )
        else:
            _print(f"  {name:6s} {path} (missing)")
    _print("")
    _print("Verbs: info, validate <case>, run [case1|case2|all] [--no-optimize], eval, figures")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    case = load_case(args.case)
    net = case.network
    lay = net.layout
    _print(f"{case.name}: {net.name}")
    _print(f"  nodes: {lay.n_nodes}, pipes: {lay.n_pipes}, wells: {case.n_wells}")
    _print(
        f"  equations/unknowns: 3N + 4S = {lay.size}; "
        f"specifications: {case.specifications.count} of N + S = {lay.n_specifications}"
    )
    for w in case.wells:
        s = case.specifications
        _print(
            f"  well {w.name}: oil {s[w.downhole].liquid_flow:g} kg/s at {w.downhole}, "
            f"gas {s[w.gas_inlet].gas_flow:g} kg/s at {w.gas_inlet}, "
            f"gas lines {net.pipe(w.gas_line_1).geometry.inlet_diameter / 0.0254:.3g} in"
        )
    ok = True
    try:
        case.specifications.validate()
    except ValueError as exc:
        _print(f"  ERROR: {exc}")
        ok = False
    issues = net.validate()
    annuli = sorted({w.name for w in case.wells if any(f"pipe '{w.annulus}'" in i for i in issues)})
    for issue in issues:
        if not any(f"pipe '{w.annulus}'" in issue for w in case.wells):
            _print(f"  WARNING: {issue}")
    if annuli:
        _print(
            f"  note: annulus geometry of {', '.join(annuli)} is written bottom-up "
            "(annulus_orientation='matlab', as in the original code; see docs/reproduction.md)"
        )
    study = case.study
    _print(
        f"  study: {len(study.get('operating_points', []))} operating points, "
        f"{len(study.get('sweeps', []))} sweeps, "
        f"{len(study.get('optimizations', []))} optimisation problems"
    )
    if args.simulate and ok:
        from .study import simulation_options

        rtol = float(study.get("integration_rtol", 1e-2))
        sol = case.simulator(simulation_options(rtol)).solve()
        _print(f"  simulation: {sol.status}")
        if sol.converged:
            for w in case.wells:
                _print(
                    f"    {w.name}: downhole {sol.pressure_bar(w.downhole):.2f} bar, "
                    f"gas inlet {sol.pressure_bar(w.gas_inlet):.2f} bar"
                )
        ok = ok and sol.converged
    _print("  OK" if ok else "  FAILED")
    return 0 if ok else 1


def _case_names(arg: str) -> list[str]:
    if arg in ("all", None):
        return list(PAPER_CASE_FILES)
    return [arg]


def cmd_run(args: argparse.Namespace) -> int:
    from .study import run_study

    out = Path(args.out)
    for name in _case_names(args.case):
        case = load_case(name)
        target = out / case.name
        summary = run_study(
            case, target, optimize=not args.no_optimize, only=args.only, progress=_print
        )
        _print(f"{case.name}: results written to {target} ({summary['elapsed_s']:.0f} s)")
    return 0


def cmd_eval(args: argparse.Namespace) -> int:
    from .evaluation import evaluate

    report = evaluate(args.results, args.out, progress=_print)
    _print(report.read_text())
    return 0


def cmd_figures(args: argparse.Namespace) -> int:
    from .figures import export_figures

    export_figures(args.out, results_dir=args.results, quick=args.quick, progress=_print)
    return 0


# --------------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gaslift",
        description=(
            "Steady-state two-phase pipe-network simulation and continuous gas-lift "
            "optimisation (de Souza et al., 2010)."
        ),
    )
    p.add_argument("--version", action="version", version=f"gaslift {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("info", help="paper reference, available cases, package version")
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("validate", help="check a case file (schema, topology, specifications)")
    s.add_argument("case", help="case JSON file, or 'case1' / 'case2'")
    s.add_argument("--simulate", action="store_true", help="also solve the steady state")
    s.set_defaults(func=cmd_validate)

    s = sub.add_parser(
        "run",
        help="simulate the paper operating points and solve the optimisation problems",
    )
    s.add_argument(
        "case",
        nargs="?",
        default="all",
        help="'case1', 'case2', 'all' (default) or a case JSON file",
    )
    s.add_argument("--out", default="results", help="output directory (default: results/)")
    s.add_argument(
        "--no-optimize",
        action="store_true",
        help="skip the optimisation problems (operating points and published optima only)",
    )
    s.add_argument(
        "--only",
        nargs="*",
        help="prefixes of the optimisation problems to run (e.g. fig13 table13)",
    )
    s.set_defaults(func=cmd_run)

    s = sub.add_parser(
        "eval", help="compare the results of `run` with the MATLAB results and the paper"
    )
    s.add_argument("--results", default="results", help="directory written by `run`")
    s.add_argument("--out", default="results/eval", help="output directory")
    s.set_defaults(func=cmd_eval)

    s = sub.add_parser("figures", help="export the figures and tables of the paper")
    s.add_argument("--out", default="figures", help="output directory (default: figures/)")
    s.add_argument(
        "--results",
        default="results",
        help="directory written by `run` (used for Figs. 13-16 and Tables 13-19)",
    )
    s.add_argument("--quick", action="store_true", help="coarser rate sweeps")
    s.set_defaults(func=cmd_figures)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
