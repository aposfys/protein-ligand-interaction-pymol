"""Independent checks on the contact detection.

A geometric contact analysis is only trustworthy if an independent
implementation, and the people who solved the structure, agree with it. Two
references are used:

* **ProLIF**, a community interaction-fingerprint library built on RDKit and
  MDAnalysis. It applies its own geometric criteria, including angular terms
  this pipeline does not model, so agreement is meaningful and disagreement is
  informative.
* **The depositors' SITE records**, the binding-site residues annotated in the
  PDB entry itself.

Neither is ground truth. The point is to make the differences explicit rather
than to claim a single correct answer.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

# SITE records list residues in fixed-width groups of four across each line.
SITE_RESIDUE = re.compile(r"([A-Z0-9]{3})\s+([A-Z])\s*(-?\d+)")

# Species that are not part of the protein and are excluded from comparison, so
# the two methods are judged on the same universe of residues.
NON_PROTEIN = frozenset(
    {"HOH", "DOD", "NDP", "DMS", "SO4", "GOL", "EDO", "MPD", "PEG", "ACT", "PO4"}
)


@dataclass(frozen=True)
class Agreement:
    """Overlap between this pipeline's pocket and a reference set."""

    reference: str
    ours: frozenset[str]
    theirs: frozenset[str]

    @property
    def shared(self) -> frozenset[str]:
        return self.ours & self.theirs

    @property
    def only_ours(self) -> frozenset[str]:
        return self.ours - self.theirs

    @property
    def only_theirs(self) -> frozenset[str]:
        return self.theirs - self.ours

    @property
    def jaccard(self) -> float:
        union = self.ours | self.theirs
        return len(self.shared) / len(union) if union else 0.0

    @property
    def recall_of_reference(self) -> float:
        """Fraction of the reference's residues this pipeline also found."""
        return len(self.shared) / len(self.theirs) if self.theirs else 0.0

    def as_dict(self) -> dict[str, object]:
        return {
            "reference": self.reference,
            "residues_ours": len(self.ours),
            "residues_reference": len(self.theirs),
            "shared": sorted(self.shared),
            "only_ours": sorted(self.only_ours),
            "only_reference": sorted(self.only_theirs),
            "jaccard": round(self.jaccard, 3),
            "recall_of_reference": round(self.recall_of_reference, 3),
        }


def site_record_residues(pdb: Path, ligand: str) -> frozenset[str]:
    """Binding-site residues the depositors annotated for a ligand.

    The SITE records themselves carry no ligand label; the association lives in
    the matching ``REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE <lig>``
    line. Selecting instead the site whose residue list *mentions* the ligand
    returns the wrong pocket, because a neighbouring ligand's site lists it too:
    in both of these entries that picks the NADPH site rather than the
    antifolate site.
    """
    text = pdb.read_text(encoding="utf-8")
    lines = text.splitlines()

    # REMARK 800 maps each site identifier to the ligand it surrounds.
    site_for_ligand: str | None = None
    current_id: str | None = None
    for line in lines:
        if not line.startswith("REMARK 800"):
            continue
        body = line[10:].strip()
        if body.startswith("SITE_IDENTIFIER:"):
            current_id = body.split(":", 1)[1].strip()
        elif body.startswith("SITE_DESCRIPTION:") and current_id:
            description = body.split(":", 1)[1].strip().upper()
            if re.search(rf"\bRESIDUE\s+{re.escape(ligand.upper())}\b", description):
                site_for_ligand = current_id
                break

    if site_for_ligand is None:
        return frozenset()

    residues: set[str] = set()
    for line in lines:
        if not line.startswith("SITE") or line[11:14].strip() != site_for_ligand:
            continue
        for resname, _chain, resseq in SITE_RESIDUE.findall(line[18:]):
            if resname not in NON_PROTEIN and resname != ligand.upper():
                residues.add(f"{resname}{int(resseq)}")

    return frozenset(residues)


def has_hydrogens(pdb: Path) -> bool:
    """True if the file contains explicit hydrogen atoms."""
    for line in pdb.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[76:78].strip().upper() == "H":
            return True
    return False


def prolif_residues(
    pdb: Path, ligand: str, cutoff: float = 5.0
) -> tuple[frozenset[str], dict[str, int]]:
    """Binding-site residues and interaction counts according to ProLIF.

    ProLIF assigns hydrogen bonds and hydrophobic contacts from explicit
    hydrogens and perceived bond orders. X-ray structures are deposited without
    hydrogens, and on such a file ProLIF silently degrades to reporting a
    handful of van der Waals contacts -- which would look like disagreement with
    this pipeline when it is really an unprepared input.

    Protonate first, for example with ``pdbfixer`` or ``reduce``, and pass the
    prepared structure here.

    Raises:
        ValueError: if the structure has no hydrogens, or lacks the ligand.
    """
    import MDAnalysis as mda
    import prolif

    if not has_hydrogens(pdb):
        raise ValueError(
            f"{pdb.name} has no explicit hydrogens. ProLIF needs a protonated "
            "structure; add hydrogens with pdbfixer or reduce first. The "
            "SITE-record comparison works on the deposited file as-is."
        )

    universe = mda.Universe(str(pdb))

    ligand_selection = universe.select_atoms(f"resname {ligand.upper()}")
    if not ligand_selection:
        raise ValueError(f"ligand {ligand!r} not found in {pdb.name}")

    excluded = " ".join(sorted(NON_PROTEIN))
    pocket_selection = universe.select_atoms(
        f"protein and not resname {excluded} and "
        f"byres (around {cutoff} resname {ligand.upper()})"
    )

    ligand_mol = prolif.Molecule.from_mda(ligand_selection, force=True)
    pocket_mol = prolif.Molecule.from_mda(pocket_selection, force=True)

    fingerprint = prolif.Fingerprint()
    fingerprint.run_from_iterable([ligand_mol], pocket_mol, progress=False)
    frame = fingerprint.to_dataframe()

    residues: set[str] = set()
    counts: dict[str, int] = {}
    for column in frame.columns:
        _ligand_residue, protein_residue, interaction = column
        if not bool(frame[column].iloc[0]):
            continue
        label = str(protein_residue)
        # ProLIF labels residues as e.g. "GLU30.A"; drop the chain suffix.
        residues.add(label.split(".")[0])
        counts[interaction] = counts.get(interaction, 0) + 1

    return frozenset(residues), counts


def compare_to_reference(
    ours: Sequence[str], theirs: frozenset[str], reference: str
) -> Agreement:
    return Agreement(reference=reference, ours=frozenset(ours), theirs=theirs)


def format_agreement(agreement: Agreement) -> str:
    """Render an agreement as a short Markdown block."""
    lines = [
        f"**{agreement.reference}** — Jaccard {agreement.jaccard:.2f}, "
        f"recovers {agreement.recall_of_reference:.0%} of the reference residues",
        "",
        f"- shared ({len(agreement.shared)}): "
        + (", ".join(sorted(agreement.shared)) or "none"),
        f"- only this pipeline ({len(agreement.only_ours)}): "
        + (", ".join(sorted(agreement.only_ours)) or "none"),
        f"- only the reference ({len(agreement.only_theirs)}): "
        + (", ".join(sorted(agreement.only_theirs)) or "none"),
    ]
    return "\n".join(lines)
