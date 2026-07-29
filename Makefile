.PHONY: install data analysis figures test clean all

PYTHON ?= python3

all: analysis

install:
	$(PYTHON) -m pip install -e ".[dev]"

## Download 1HFR and 1KMV from RCSB into data/pdb/
data:
	$(PYTHON) -c "from pathlib import Path; \
	from plinter.structures import TARGETS, download_structure; \
	[download_structure(t.pdb_id, Path('data/pdb')/t.filename) for t in TARGETS]"

## Run the contact analysis and write results/
analysis: data
	$(PYTHON) -m plinter.cli

## Render PyMOL binding-site figures (requires a PyMOL installation)
figures: data
	pymol -cq pymol/render_binding_sites.py

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*.csv results/*.json results/*.png results/figures
	find . -name __pycache__ -type d -exec rm -rf {} +
