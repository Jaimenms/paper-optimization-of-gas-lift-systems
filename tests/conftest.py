import json
from pathlib import Path

import pytest

REFERENCE_DIR = Path(__file__).resolve().parents[1] / "data" / "reference"


def load_reference(name: str) -> dict:
    with open(REFERENCE_DIR / name) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def case1_grid():
    return load_reference("case1_grid_2in.json")


@pytest.fixture(scope="session")
def case1_grid_3in():
    return load_reference("case1_grid_3in.json")


@pytest.fixture(scope="session")
def case2_grid():
    return load_reference("case2_grid.json")


@pytest.fixture(scope="session")
def matlab_optima():
    return load_reference("matlab_optimization_results.json")
