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

Every residue reported is checked against two independent references. Against the depositors' SITE records: **100% recall in both structures** (14/14 and 10/10), with 7 and 8 additional residues from the 5.0 Å cutoff. Against [PLIP](https://doi.org/10.1093/nar/gkab294), which types interactions from perceived chemistry rather than from distance rules: **Arg70 and Asn64 appear for MOT and for neither does LII, and Glu30 anchors both** — the discriminating result, reproduced without this pipeline's criteria.

PLIP also finds two classes this pipeline cannot: **π-stacking** with Phe31 in both complexes, which has no geometric rule here, and **water-mediated bridges**, which are invisible in principle because waters are stripped before contact detection. Both blind spots are asserted by tests. PLIP's residue counts (9 and 5) are smaller than this pipeline's (21 and 18) because it reports only typed interactions rather than every contact under a cutoff; the two are not comparable.

Cα RMSD between the two structures is 0.51 Å, so the differences are ligand-driven rather than conformational. The entries differ in resolution, so sub-0.1 Å distance differences should not be over-interpreted.

### Quick start

```
pip install -e ".[dev]"   # or: pip install -r requirements.txt
make analysis             # downloads PDB files, writes results/
make figures              # PyMOL renders (requires a PyMOL install)
make test
```

The base install needs only Biopython; `pip install -e ".[chemistry]"` adds RDKit for the ligand descriptors.

### Prior work

Human DHFR and its antifolates are among the most thoroughly characterised protein–ligand
systems in structural biology, and both entries analysed here have been in the PDB for over
two decades. The 2,4-diaminopyrimidine–Glu30 pharmacophore, the role of Arg70 in binding the
glutamate tail of classical antifolates, and the lipophilic-analogue design rationale are all
long-established. **Nothing in the results above is a new observation about DHFR.**

PLIP (Adasme et al., *Nucleic Acids Research* 2021) is the standard tool for the profiling
this pipeline does, and it is used here as an independent check rather than being
reimplemented.

What this repository is: a small, fully reproducible contact-analysis pipeline whose every
reported residue is checked against two independent references, whose blind spots
(π-stacking, water-mediated bridges) are asserted by tests rather than left implicit, and
whose binding-site finding is reproduced by a different route in a different repository
([`dhfr-campaign`](https://github.com/aposfys/dhfr-campaign)). It is a teaching and
verification artefact, not a research contribution.

### More

- [Method, design decisions and output files](docs/METHOD.md)
- [Data sources and licences](docs/DATA.md)

---

Apostolos Fysekidis · [MIT Licence](LICENSE)
