"""Detection and classification of protein-ligand atomic contacts.

Contacts are found with a KD-tree neighbour search over heavy atoms, then
classified with the geometric heuristics that are standard for X-ray
structures deposited without hydrogens (Ferreira de Freitas & Schapira, 2017).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable, Iterator, Sequence

from Bio.PDB import NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure

# Heavy-atom distance below which two atoms are recorded as being in contact.
CONTACT_CUTOFF = 5.0
# Donor/acceptor distance below which a polar contact is called a hydrogen bond.
HBOND_CUTOFF = 3.5
# Carbon-carbon distance below which a contact is called hydrophobic.
HYDROPHOBIC_CUTOFF = 4.5
# Charged-group distance below which a polar contact is also called a salt bridge.
SALT_BRIDGE_CUTOFF = 4.0

# Elements able to donate or accept a hydrogen bond. Hydrogens are absent from
# most X-ray depositions, so donor/acceptor identity is inferred from element.
POLAR_ELEMENTS = frozenset({"N", "O"})

# Formally charged side-chain atoms, used to flag salt bridges.
CATIONIC_ATOMS = {
    ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"),
    ("LYS", "NZ"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}
ANIONIC_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"),
    ("GLU", "OE1"), ("GLU", "OE2"),
}

# Non-biological species routinely present in crystals: cryoprotectants,
# precipitants and buffer components. They are not part of the binding site and
# are excluded from the protein environment by default.
CRYSTALLISATION_ADDITIVES = frozenset(
    {
        "HOH", "DOD",                      # water
        "DMS", "GOL", "EDO", "MPD", "PEG", # cryoprotectants
        "SO4", "PO4", "ACT", "CIT", "TRS", # buffer / precipitant ions
        "MES", "EPE", "FMT", "IMD", "NO3",
        "CL", "NA", "K", "MG", "CA", "ZN",
    }
)


@dataclass(frozen=True)
class Contact:
    """A single ligand-atom / environment-atom pair within the cutoff."""

    ligand_atom: str
    ligand_element: str
    residue: str          # e.g. "GLU30"
    residue_name: str
    residue_seq: int
    chain: str
    protein_atom: str
    protein_element: str
    distance: float
    interaction: str      # hydrogen bond | salt bridge | hydrophobic | van der Waals
    hydrogen_bond: bool

    def as_row(self) -> dict:
        return asdict(self)


def _element(atom: Atom) -> str:
    """Return the element symbol, falling back to the PDB atom-name convention."""
    element = (atom.element or "").strip().upper()
    if element:
        return element
    # Some older PDB files leave columns 77-78 blank; the first alphabetic
    # character of the atom name is the element for all standard residues.
    return next((c for c in atom.get_name() if c.isalpha()), "").upper()


def _is_hydrogen(atom: Atom) -> bool:
    return _element(atom) in {"H", "D"}


def _classify(
    ligand_atom: Atom,
    partner_atom: Atom,
    partner_residue: str,
    distance: float,
) -> tuple[str, bool]:
    """Assign an interaction type to a contact and flag hydrogen bonds."""
    lig_element = _element(ligand_atom)
    par_element = _element(partner_atom)

    both_polar = lig_element in POLAR_ELEMENTS and par_element in POLAR_ELEMENTS
    if both_polar and distance <= HBOND_CUTOFF:
        partner_key = (partner_residue, partner_atom.get_name())
        charged_partner = partner_key in CATIONIC_ATOMS or partner_key in ANIONIC_ATOMS
        if charged_partner and distance <= SALT_BRIDGE_CUTOFF:
            return "salt bridge", True
        return "hydrogen bond", True

    if lig_element == "C" and par_element == "C" and distance <= HYDROPHOBIC_CUTOFF:
        return "hydrophobic", False

    return "van der Waals", False


def iter_ligand_atoms(structure: Structure, ligand_resname: str) -> Iterator[Atom]:
    """Yield every heavy atom belonging to the named ligand."""
    for residue in structure.get_residues():
        if residue.get_resname().strip() == ligand_resname:
            for atom in residue:
                if not _is_hydrogen(atom):
                    yield atom


def _environment_atoms(
    structure: Structure,
    exclude_resnames: Iterable[str],
) -> list[Atom]:
    """Collect the heavy atoms that make up the binding-site environment."""
    excluded = {name.strip().upper() for name in exclude_resnames}
    return [
        atom
        for residue in structure.get_residues()
        if residue.get_resname().strip().upper() not in excluded
        for atom in residue
        if not _is_hydrogen(atom)
    ]


def find_contacts(
    structure: Structure,
    ligand_resname: str,
    exclude_resnames: Sequence[str] = (),
    cutoff: float = CONTACT_CUTOFF,
) -> list[Contact]:
    """Find every heavy-atom contact between a ligand and its environment.

    Args:
        structure: A parsed Biopython structure.
        ligand_resname: PDB chemical component ID of the ligand, e.g. ``"MOT"``.
        exclude_resnames: Residue names to remove from the environment, on top of
            the ligand itself and the default crystallisation additives.
        cutoff: Maximum heavy-atom separation, in angstroms.

    Returns:
        Contacts sorted by increasing distance.
    """
    ligand_atoms = list(iter_ligand_atoms(structure, ligand_resname))
    if not ligand_atoms:
        raise ValueError(
            f"Ligand {ligand_resname!r} not found in {structure.get_id()!r}"
        )

    excluded = set(CRYSTALLISATION_ADDITIVES) | {ligand_resname.upper()}
    excluded.update(name.strip().upper() for name in exclude_resnames)

    environment = _environment_atoms(structure, excluded)
    search = NeighborSearch(environment)

    contacts: list[Contact] = []
    for ligand_atom in ligand_atoms:
        for partner in search.search(ligand_atom.coord, cutoff):
            distance = float(ligand_atom - partner)
            residue = partner.get_parent()
            resname = residue.get_resname().strip()
            resseq = residue.get_id()[1]
            interaction, is_hbond = _classify(ligand_atom, partner, resname, distance)
            contacts.append(
                Contact(
                    ligand_atom=ligand_atom.get_name(),
                    ligand_element=_element(ligand_atom),
                    residue=f"{resname}{resseq}",
                    residue_name=resname,
                    residue_seq=resseq,
                    chain=residue.get_parent().get_id(),
                    protein_atom=partner.get_name(),
                    protein_element=_element(partner),
                    distance=round(distance, 2),
                    interaction=interaction,
                    hydrogen_bond=is_hbond,
                )
            )

    contacts.sort(key=lambda c: c.distance)
    return contacts


def contacting_residues(contacts: Sequence[Contact]) -> dict[str, float]:
    """Map each contacted residue to its shortest contact distance."""
    closest: dict[str, float] = {}
    for contact in contacts:
        key = contact.residue
        if contact.distance < closest.get(key, float("inf")):
            closest[key] = contact.distance
    return dict(sorted(closest.items(), key=lambda item: item[1]))


def hydrogen_bonds(contacts: Sequence[Contact]) -> list[Contact]:
    """Return only the contacts that satisfy the hydrogen-bond criterion."""
    return [c for c in contacts if c.hydrogen_bond]
