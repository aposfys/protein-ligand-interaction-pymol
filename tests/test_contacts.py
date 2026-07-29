"""Tests for contact detection and classification."""

from __future__ import annotations

from pathlib import Path

import pytest

from plinter.compare import Comparison, build_fingerprint
from plinter.contacts import (
    CRYSTALLISATION_ADDITIVES,
    contacting_residues,
    find_contacts,
    hydrogen_bonds,
)
from plinter.structures import COFACTOR, TARGETS, load_structure

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "pdb"


@pytest.fixture(scope="module")
def fingerprints():
    prints = []
    for target in TARGETS:
        path = DATA_DIR / target.filename
        if not path.exists():
            pytest.skip(f"{path} missing; run `make data` first")
        contacts = find_contacts(
            load_structure(path),
            ligand_resname=target.ligand,
            exclude_resnames=(COFACTOR,),
        )
        prints.append(build_fingerprint(target.pdb_id, target.ligand, contacts))
    return prints


def test_every_ligand_finds_contacts(fingerprints):
    for fingerprint in fingerprints:
        assert fingerprint.contacts, f"no contacts for {fingerprint.ligand}"


def test_contacts_are_sorted_by_distance(fingerprints):
    for fingerprint in fingerprints:
        distances = [c.distance for c in fingerprint.contacts]
        assert distances == sorted(distances)


def test_contacts_respect_the_cutoff(fingerprints):
    for fingerprint in fingerprints:
        assert all(c.distance <= 5.0 for c in fingerprint.contacts)


def test_cofactor_and_solvent_are_excluded(fingerprints):
    """NADPH, waters and cryoprotectants must not appear as binding-site residues."""
    for fingerprint in fingerprints:
        names = {c.residue_name for c in fingerprint.contacts}
        assert COFACTOR not in names
        assert not (names & CRYSTALLISATION_ADDITIVES)


def test_ligand_never_contacts_itself(fingerprints):
    for fingerprint in fingerprints:
        assert fingerprint.ligand not in {c.residue_name for c in fingerprint.contacts}


def test_hydrogen_bonds_are_polar_and_short(fingerprints):
    for fingerprint in fingerprints:
        for contact in hydrogen_bonds(fingerprint.contacts):
            assert contact.distance <= 3.5
            assert contact.ligand_element in {"N", "O"}
            assert contact.protein_element in {"N", "O"}


def test_glu30_anchors_both_inhibitors(fingerprints):
    """Glu30 is the conserved hDHFR anchor for the 2,4-diaminopyrimidine ring."""
    for fingerprint in fingerprints:
        assert "GLU30" in fingerprint.hbond_residues


def test_closest_approach_matches_raw_contacts(fingerprints):
    for fingerprint in fingerprints:
        closest = contacting_residues(fingerprint.contacts)
        for contact in fingerprint.contacts:
            assert closest[contact.residue] <= contact.distance


def test_comparison_partitions_residues(fingerprints):
    comparison = Comparison(first=fingerprints[0], second=fingerprints[1])
    first, second = comparison.first, comparison.second

    assert comparison.shared_residues == first.residues & second.residues
    assert comparison.unique_to(first).isdisjoint(second.residues)
    assert comparison.unique_to(second).isdisjoint(first.residues)

    partitioned = (
        comparison.shared_residues
        | comparison.unique_to(first)
        | comparison.unique_to(second)
    )
    assert partitioned == first.residues | second.residues


def test_residue_table_covers_every_residue(fingerprints):
    comparison = Comparison(first=fingerprints[0], second=fingerprints[1])
    rows = comparison.residue_table()
    assert {row["residue"] for row in rows} == (
        comparison.first.residues | comparison.second.residues
    )


def test_missing_ligand_raises():
    path = DATA_DIR / TARGETS[0].filename
    if not path.exists():
        pytest.skip(f"{path} missing; run `make data` first")
    with pytest.raises(ValueError, match="not found"):
        find_contacts(load_structure(path), ligand_resname="XYZ")
