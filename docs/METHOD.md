# Method

1. **Retrieve** both entries from RCSB (cached in `data/pdb/`).
2. **Define the environment.** Waters, the NADPH cofactor and crystallisation additives are
   removed. Every remaining heavy atom is indexed in a KD-tree.
3. **Find contacts.** For each ligand heavy atom, `Bio.PDB.NeighborSearch` returns all
   environment atoms within 5.0 Å.
4. **Classify** each contact geometrically, as is standard for X-ray structures deposited
   without hydrogens:

   | Interaction | Criterion |
   | --- | --- |
   | Salt bridge | N/O ↔ formally charged side-chain N/O, ≤ 4.0 Å |
   | Hydrogen bond | N/O ↔ N/O, ≤ 3.5 Å |
   | Hydrophobic | C ↔ C, ≤ 4.5 Å |
   | van der Waals | any other pair ≤ 5.0 Å |

5. **Compare** the two fingerprints: shared residues, shared hydrogen-bond partners, and
   residues unique to each ligand.

## The chemistry explains the binding modes

RDKit descriptors over the PDB Chemical Component Dictionary SMILES turn the structural
observation into a mechanistic one:

| | MOT (classical) | LII (lipophilic) |
| --- | ---: | ---: |
| Molecular weight | 442.4 | 337.4 |
| cLogP | +1.07 | **+2.69** |
| Topological polar surface area | 197.9 Å² | **109.2 Å²** |
| H-bond donors / acceptors | 5 / 9 | 2 / 7 |
| Rotatable bonds | 9 | 4 |
| **Carboxylic acids** | **2** | **0** |
| Primary aromatic amines | 2 | 2 |

The decisive rows are the last two. Both ligands carry the 2,4-diaminopyrimidine head group
— the two primary aromatic amines that hydrogen-bond to Glu30 — which is why both are
anchored identically. But MOT's L-glutamate tail contributes two carboxylates, ionised at
physiological pH, and LII has none. That is the physical reason LII cannot reach Arg70.

## Validation against the depositors' annotation

| | Annotated residues | Recovered | Additional |
| --- | ---: | ---: | ---: |
| 1HFR / MOT | 14 | **14 (100%)** | 7 |
| 1KMV / LII | 10 | **10 (100%)** | 8 |

Linking a SITE record to its ligand requires the `REMARK 800 SITE_DESCRIPTION` line:
selecting instead the site whose residue list *mentions* the ligand returns the NADPH pocket
in both of these entries, because adjacent ligands appear in each other's site lists. A test
pins that distinction.

## Independent profiling with PLIP

The SITE records say *which residues* line the pocket; they say nothing about what kind of
interaction each one makes. Those assignments come from this pipeline's own distance rules,
so checking them needs a second implementation that types interactions from chemistry rather
than from geometry alone.

[PLIP](https://doi.org/10.1093/nar/gkab294) (Adasme et al., *Nucleic Acids Research* 2021)
is that reference. It protonates the entry internally with OpenBabel, so unlike ProLIF it
runs on a deposited X-ray file as-is — which matters here because neither structure carries
explicit hydrogens, and `prolif_residues` refuses both by design.

| | Residues with a typed interaction | H-bond | Hydrophobic | Salt bridge | π-stacking | Water bridge |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1HFR / MOT | 9 | 4 | 3 | 2 | 1 | 2 |
| 1KMV / LII | 5 | 2 | 1 | 1 | 1 | 1 |

**PLIP reproduces the discriminating result by an independent route.** Arg70 and Asn64 appear
for MOT and neither appears for LII, and Glu30 anchors both — the same conclusion the contact
analysis reaches, reached without its distance rules. A test asserts each of these rather than
leaving them to be read off a table.

**PLIP also finds two interaction classes this pipeline cannot see.**

- **π-stacking** with Phe31 in both complexes. There is no geometric rule for ring-centroid
  geometry in `contacts.py`, so these contacts are currently absorbed into the hydrophobic
  and van der Waals counts.
- **Water-mediated bridges**, two for MOT and one for LII. Waters are stripped before the
  KD-tree is built, so a bridging water is invisible in principle rather than merely missed.

Both are asserted by tests, so the blind spots cannot quietly close without the documentation
changing with them.

**The two residue counts are not comparable and are not presented as if they were.** This
pipeline reports every heavy-atom contact within 5.0 Å — 21 and 18 residues. PLIP reports
only chemically typed interactions — 9 and 5. A smaller number here is not disagreement; the
two methods answer different questions, and a test pins the direction of the inequality.

## Design decisions

Each of these is a place where a contact analysis can quietly go wrong:

- **Crystallisation additives are excluded from the binding site.** 1KMV contains DMSO
  (`DMS 203`) within 3 Å of LII. DMSO is a cryoprotectant, not part of the protein, and
  counting it as a binding-site partner invents an interaction that does not exist in
  solution. `CRYSTALLISATION_ADDITIVES` in [`contacts.py`](../src/plinter/contacts.py) lists
  what is filtered.
- **KD-tree neighbour search, not a nested loop.** Testing every ligand × environment atom
  pair is O(n×m) and, through PyMOL's `cmd.get_distance`, slow enough to be measured in
  minutes. `Bio.PDB.NeighborSearch` gives identical distances in under a second and removes
  the dependency on a running PyMOL session for the numerical work.
- **Interaction typing beyond "hydrogen bond: yes/no".** Salt bridges are separated from
  neutral hydrogen bonds and hydrophobic contacts are identified explicitly. This is what
  makes the two binding modes distinguishable; a boolean H-bond flag cannot express it.
- **Descriptors come from curated SMILES, not from coordinates.** Bond orders and
  protonation inferred from a 2.10 Å electron-density model would be guesses.
- **The pocket is checked against an external reference.** A geometric analysis that agrees
  with nobody is hard to trust.
- **Everything regenerates from source.** 19 tests assert the invariants — cutoffs
  respected, cofactor and solvent excluded, the conserved Glu30 anchor recovered, and full
  recall of the annotated binding site.

## CLI

```
python -m plinter.cli --cutoff 4.5 --results-dir results/strict
python -m plinter.cli --keep-cofactor    # include NADPH in the environment
python -m plinter.cli --no-chemistry     # skip RDKit descriptors
python -m plinter.cli --no-validation    # skip the SITE-record comparison
```

## Repository layout

```
src/plinter/
  contacts.py     KD-tree contact detection and interaction classification
  compare.py      Binding-site fingerprints and their intersection
  structures.py   RCSB retrieval and PDB parsing
  chemistry.py    RDKit descriptors from the PDB Chemical Component Dictionary
  validation.py   SITE-record and ProLIF cross-checks
  report.py       CSV / JSON / Markdown output
  plots.py        Matplotlib figures
  cli.py          Command-line entry point
pymol/
  render_binding_sites.py   Headless PyMOL figure rendering
tests/            pytest suite
results/          Generated tables and figures
```

## Output files

| File | Contents |
| --- | --- |
| `results/contacts_1hfr_MOT.csv` | Every MOT contact: atoms, residue, distance, interaction type |
| `results/contacts_1kmv_LII.csv` | The same for LII |
| `results/residue_comparison.csv` | One row per pocket residue, both ligands' closest approach |
| `results/summary.json` | Per-ligand counts and the shared/unique residue partition |
| `results/ligand_properties.csv` | RDKit descriptors and functional-group counts |
| `results/validation.csv` | Agreement with the depositors' SITE records |
| `results/extras.json` | Chemistry and validation, machine-readable |
| `results/figures/` | PyMOL renders of each site and the superposition |

## References

1. Ferreira de Freitas, R. & Schapira, M. (2017). A systematic analysis of atomic
   protein–ligand interactions in the PDB. *MedChemComm* **8**, 1970–1981.
2. Cody, V., Galitsky, N., Luft, J. R., Pangborn, W., Blakley, R. L. & Gangjee, A. (1998).
   Comparison of ternary crystal complexes of F31 variants of human dihydrofolate reductase
   with NADPH and a classical antitumor furopyrimidine. *Anti-Cancer Drug Design* **13**,
   307–315. PDB **1HFR**.
3. Klon, A. E., Héroux, A., Ross, L. J., Pathak, V., Johnson, C. A., Piper, J. R. &
   Borhani, D. W. (2002). Atomic structures of human dihydrofolate reductase complexed with
   NADPH and two lipophilic antifolates at 1.09 Å and 1.05 Å resolution. *Journal of
   Molecular Biology* **320**, 677–693. PDB **1KMV**.
4. Cock, P. J. A. *et al.* (2009). Biopython: freely available Python tools for
   computational molecular biology and bioinformatics. *Bioinformatics* **25**, 1422–1423.
5. Schrödinger, LLC. The PyMOL Molecular Graphics System.
6. Landrum, G. RDKit: Open-source cheminformatics. https://www.rdkit.org
7. Bouysset, C. & Fiorucci, S. (2021). ProLIF: a library to encode molecular interactions as
   fingerprints. *Journal of Cheminformatics* **13**, 72.
