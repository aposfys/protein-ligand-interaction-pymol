"""Comparison of the binding-site fingerprints of two ligands."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from .contacts import Contact, contacting_residues, hydrogen_bonds


@dataclass(frozen=True)
class Fingerprint:
    """The set of residues one ligand touches, split by interaction type."""

    ligand: str
    pdb_id: str
    contacts: tuple[Contact, ...]

    @property
    def residues(self) -> set[str]:
        return {c.residue for c in self.contacts}

    @property
    def hbond_residues(self) -> set[str]:
        return {c.residue for c in self.contacts if c.hydrogen_bond}

    @property
    def hydrophobic_residues(self) -> set[str]:
        return {c.residue for c in self.contacts if c.interaction == "hydrophobic"}

    @property
    def closest_approach(self) -> dict[str, float]:
        return contacting_residues(self.contacts)

    def summary(self) -> dict[str, object]:
        return {
            "pdb_id": self.pdb_id,
            "ligand": self.ligand,
            "ligand_atoms_in_contact": len({c.ligand_atom for c in self.contacts}),
            "contacts": len(self.contacts),
            "residues_contacted": len(self.residues),
            "hydrogen_bonds": len(hydrogen_bonds(self.contacts)),
            "salt_bridges": sum(1 for c in self.contacts if c.interaction == "salt bridge"),
            "hydrophobic_contacts": sum(
                1 for c in self.contacts if c.interaction == "hydrophobic"
            ),
            "shortest_contact_A": min((c.distance for c in self.contacts), default=None),
        }


@dataclass(frozen=True)
class Comparison:
    """Shared and ligand-specific binding-site residues."""

    first: Fingerprint
    second: Fingerprint

    @property
    def shared_residues(self) -> set[str]:
        return self.first.residues & self.second.residues

    @property
    def shared_hbond_residues(self) -> set[str]:
        return self.first.hbond_residues & self.second.hbond_residues

    def unique_to(self, fingerprint: Fingerprint) -> set[str]:
        other = self.second if fingerprint is self.first else self.first
        return fingerprint.residues - other.residues

    def residue_table(self) -> list[dict[str, object]]:
        """One row per residue contacted by either ligand, with both distances."""
        first_closest = self.first.closest_approach
        second_closest = self.second.closest_approach

        rows: list[dict[str, object]] = []
        for residue in sorted(
            self.first.residues | self.second.residues, key=_residue_sort_key
        ):
            in_first = residue in first_closest
            in_second = residue in second_closest
            rows.append(
                {
                    "residue": residue,
                    f"{self.first.ligand}_min_distance_A": first_closest.get(residue),
                    f"{self.second.ligand}_min_distance_A": second_closest.get(residue),
                    f"{self.first.ligand}_hbond": residue in self.first.hbond_residues,
                    f"{self.second.ligand}_hbond": residue in self.second.hbond_residues,
                    "shared": in_first and in_second,
                }
            )
        return rows


def _residue_sort_key(residue: str) -> tuple[int, str]:
    """Sort residue labels such as ``GLU30`` by sequence number."""
    digits = "".join(c for c in residue if c.isdigit())
    return (int(digits) if digits else 0, residue)


def build_fingerprint(pdb_id: str, ligand: str, contacts: Sequence[Contact]) -> Fingerprint:
    return Fingerprint(ligand=ligand, pdb_id=pdb_id, contacts=tuple(contacts))
