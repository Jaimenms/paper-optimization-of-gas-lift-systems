# gaslift - series conventions: see docs/conventions.md
.PHONY: setup lint format test test-all run eval figures notebooks clean

setup:            ## create .venv with uv and install the package with all extras
	uv sync --all-extras

lint:             ## static checks
	uv run ruff check src tests scripts
	uv run ruff format --check src tests scripts

format:           ## format the code
	uv run ruff format src tests scripts
	uv run ruff check --fix src tests scripts

test:             ## fast tests (seconds)
	uv run pytest

test-all:         ## all tests, including the slow ones (optimisations, minutes)
	uv run pytest -m ""

run:              ## run both paper cases (simulations + optimisations, ~10 min) and write results/
	uv run gaslift run --out results

eval:             ## compare results/ with the MATLAB results and the published tables -> results/eval/
	uv run gaslift eval --results results --out results/eval

figures:          ## export the figures and tables of the paper to figures/
	uv run gaslift figures --out figures --results results

notebooks:        ## regenerate and execute the notebooks (~15 min)
	uv run python scripts/build_notebooks.py --execute

clean:            ## remove caches and build artefacts
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
	find . -name ".ipynb_checkpoints" -type d -prune -exec rm -rf {} +
