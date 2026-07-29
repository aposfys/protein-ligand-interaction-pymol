"""Rendering of contact tables and the comparative summary."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Sequence

from .compare import Comparison, Fingerprint
from .contacts import Contact

CONTACT_FIELDS = [
    "ligand_atom",
    "ligand_element",
    "residue",
    "residue_name",
    "residue_seq",
    "chain",
    "protein_atom",
    "protein_element",
    "distance",
    "interaction",
    "hydrogen_bond",
]


def write_contacts_csv(contacts: Sequence[Contact], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONTACT_FIELDS)
        writer.writeheader()
        writer.writerows(contact.as_row() for contact in contacts)
    return path


def write_rows_csv(rows: Sequence[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_summary_json(
    fingerprints: Sequence[Fingerprint], comparison: Comparison, path: Path
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "per_ligand": [fp.summary() for fp in fingerprints],
        "comparison": {
            "shared_residues": sorted(comparison.shared_residues),
            "shared_hydrogen_bond_residues": sorted(comparison.shared_hbond_residues),
            f"unique_to_{comparison.first.ligand}": sorted(
                comparison.unique_to(comparison.first)
            ),
            f"unique_to_{comparison.second.ligand}": sorted(
                comparison.unique_to(comparison.second)
            ),
        },
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def format_top_contacts(contacts: Sequence[Contact], limit: int = 10) -> str:
    """Render the shortest contacts as a Markdown table."""
    header = (
        "| Ligand atom | Residue / atom | Distance (Å) | Interaction |\n"
        "| --- | --- | --- | --- |\n"
    )
    rows = "".join(
        f"| {c.ligand_atom} | {c.residue}/{c.protein_atom} | {c.distance:.2f} "
        f"| {c.interaction} |\n"
        for c in contacts[:limit]
    )
    return header + rows


def format_comparison(comparison: Comparison) -> str:
    """Render the shared/unique residue breakdown as Markdown."""
    first, second = comparison.first, comparison.second
    lines = [
        f"**Shared binding-site residues ({len(comparison.shared_residues)}):** "
        + ", ".join(sorted(comparison.shared_residues)),
        "",
        "**Shared hydrogen-bond partners:** "
        + (", ".join(sorted(comparison.shared_hbond_residues)) or "none"),
        "",
        f"**Unique to {first.ligand} ({first.pdb_id}):** "
        + (", ".join(sorted(comparison.unique_to(first))) or "none"),
        "",
        f"**Unique to {second.ligand} ({second.pdb_id}):** "
        + (", ".join(sorted(comparison.unique_to(second))) or "none"),
    ]
    return "\n".join(lines)
