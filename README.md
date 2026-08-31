# Protein–Ligand Interaction Analysis of Human DHFR Antifolate Complexes
Comparative contact analysis of a classical and a lipophilic antifolate bound to human dihydrofolate reductase.

[![Pipeline](https://github.com/aposfys/protein-ligand-interaction-pymol/actions/workflows/pipeline.yml/badge.svg)](https://github.com/aposfys/protein-ligand-interaction-pymol/actions/workflows/pipeline.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

The pipeline maps every heavy-atom contact each inhibitor makes with the binding site, classifies the interactions, and quantifies which pocket residues are shared versus inhibitor-specific. Structures are [1HFR](https://www.rcsb.org/structure/1HFR) (2.10 Å, hDHFR + NADPH + **MOT**) and [1KMV](https://www.rcsb.org/structure/1KMV) (1.05 Å, hDHFR + NADPH + **LII**).

<p align="center">
  <img src="results/figures/overlay_MOT_vs_LII.png" width="720" alt="MOT and LII superposed in the hDHFR active site">
</p>

### Results

| Metric | MOT (1HFR) | LII (1KMV) |
| --- | ---: | ---: |
| Contacts ≤ 5 Å | 330 | 274 |
| Residues contacted | 21 | 18 |
| Hydrogen bonds | 10 | 4 |
| Salt bridges | 7 | 2 |
| Hydrophobic contacts | 66 | 74 |
| Shortest contact | 2.25 Å (Arg70) | 2.77 Å (Glu30) |

Both inhibitors are anchored by the same conserved pharmacophore — a bidentate salt bridge from the 2,4-diaminopyrimidine head to **Glu30**, plus hydrogen bonds to the **Ile7** and **Val115** backbone carbonyls — and share 17 pocket residues. They then diverge: MOT's L-glutamate tail contributes two ionised carboxylates and reaches **Arg70** and **Asn64**; LII has no anion to pair with the guanidinium and compensates with hydrophobic packing against **Phe34** and **Pro61**. That trade-off is the design lever.

Every residue reported is checked against the depositors' SITE records: **100% recall in both structures** (14/14 and 10/10), with 7 and 8 additional residues from the 5.0 Å cutoff.

Cα RMSD between the two structures is 0.51 Å, so the differences are ligand-driven rather than conformational. The entries differ in resolution, so sub-0.1 Å distance differences should not be over-interpreted.

### Quick start

```
pip install -e ".[dev]"   # or: pip install -r requirements.txt
make analysis             # downloads PDB files, writes results/
make figures              # PyMOL renders (requires a PyMOL install)
make test
```

The base install needs only Biopython; `pip install -e ".[chemistry]"` adds RDKit for the ligand descriptors.

### More

- [Method, design decisions and output files](docs/METHOD.md)
- [Data sources and licences](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
