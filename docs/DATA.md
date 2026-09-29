# Data sources and licences

The MIT licence covers **the code in this repository only**. The data it retrieves, and the
third-party tools it invokes, carry their own terms.

| Source | Used for | Licence |
| --- | --- | --- |
| [RCSB PDB](https://www.rcsb.org/pages/usage-policy) | Structures 1HFR and 1KMV; chemical component definitions | CC0 / public domain |
| [Biopython](https://github.com/biopython/biopython/blob/master/LICENSE.rst) | Structure parsing, KD-tree neighbour search | Biopython Licence (BSD-like) |
| [RDKit](https://github.com/rdkit/rdkit/blob/master/license.txt) | Ligand descriptors | BSD-3-Clause |
| [PyMOL](https://github.com/schrodinger/pymol-open-source) | Molecular graphics (optional) | Open-source PyMOL licence |
| [PLIP](https://github.com/pharmai/plip) | Independent interaction profiling (optional) | GPL-2.0-only |
| [Open Babel](https://github.com/openbabel/openbabel) | Protonation inside PLIP (optional) | GPL-2.0-only |

**What this repository ships.** `data/pdb/` contains the two PDB entries and their chemical
component records, unmodified, as a convenience so the analysis runs offline. Those are
CC0 / public domain, so redistribution carries no conditions. Everything under `results/`
is generated from them by the code here. PLIP and Open Babel are imported only when installed
and are not redistributed here.
