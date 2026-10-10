"""Write the JSON case files of ``cases/`` from the Python case builders.

The files are the input of the ``gaslift`` command line (``validate``, ``run``,
``figures``). Re-run this script after changing ``gaslift.cases`` or
``gaslift.study.default_study``::

    uv run python scripts/export_cases.py
"""

from __future__ import annotations

from pathlib import Path

from gaslift.cases import PAPER_CASE_FILES, case1_single_well, case2_multi_well, save_case
from gaslift.study import default_study

ROOT = Path(__file__).resolve().parents[1]
BUILDERS = {"case1": case1_single_well, "case2": case2_multi_well}


def main() -> None:
    for name, builder in BUILDERS.items():
        case = builder()
        case.study = default_study(name)
        path = save_case(case, ROOT / "cases" / PAPER_CASE_FILES[name])
        print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
