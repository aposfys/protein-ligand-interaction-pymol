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
        comparison.shared_residues | comparison.unique_to(first) | comparison.unique_to(second)
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


# --- ligand chemistry --------------------------------------------------------


@pytest.fixture(scope="module")
def ligand_profiles():
    from plinter import chemistry

    cache = Path(__file__).resolve().parents[1] / "data"
    try:
        return {target.ligand: chemistry.profile(target.ligand, cache) for target in TARGETS}
    except Exception as error:  # offline and uncached
        pytest.skip(f"chemical component data unavailable: {error}")


def test_molecular_weights_match_the_deposited_formulae(ligand_profiles):
    assert ligand_profiles["MOT"].molecular_weight == pytest.approx(442.4, abs=0.5)
    assert ligand_profiles["LII"].molecular_weight == pytest.approx(337.4, abs=0.5)


def test_the_lipophilic_antifolate_is_more_lipophilic(ligand_profiles):
    """LII is described as lipophilic; cLogP and TPSA should both agree."""
    assert ligand_profiles["LII"].logp > ligand_profiles["MOT"].logp
    assert ligand_profiles["LII"].tpsa < ligand_profiles["MOT"].tpsa


def test_only_the_classical_antifolate_carries_carboxylates(ligand_profiles):
    """The glutamate tail is what reaches the Arg70 subsite; LII has none."""
    from plinter import chemistry

    assert chemistry.ionisable_groups(ligand_profiles["MOT"])["carboxylic_acid"] == 2
    assert chemistry.ionisable_groups(ligand_profiles["LII"])["carboxylic_acid"] == 0


def test_both_share_the_diaminopyrimidine_amines(ligand_profiles):
    """The 2,4-diamino head group is the conserved recognition motif."""
    from plinter import chemistry

    for profile in ligand_profiles.values():
        assert chemistry.ionisable_groups(profile)["primary_aromatic_amine"] == 2


# --- validation --------------------------------------------------------------


def test_site_records_are_linked_through_remark_800():
    """Matching on residue names alone selects the NADPH site, not the ligand's."""
    from plinter import validation

    residues = validation.site_record_residues(DATA_DIR / "1kmv.pdb", "LII")
    assert "GLU30" in residues  # in the antifolate site
    assert "ASP21" not in residues  # in the NADPH site


def test_missing_ligand_site_returns_empty():
    from plinter import validation

    assert validation.site_record_residues(DATA_DIR / "1hfr.pdb", "XYZ") == frozenset()


def test_pipeline_recovers_every_annotated_residue(fingerprints):
    """No false negatives against the depositors' own annotation."""
    from plinter import validation

    for target, fingerprint in zip(TARGETS, fingerprints, strict=True):
        reference = validation.site_record_residues(DATA_DIR / target.filename, target.ligand)
        assert reference, f"no SITE record for {target.ligand}"
        agreement = validation.compare_to_reference(
            sorted(fingerprint.residues), reference, target.pdb_id
        )
        assert agreement.recall_of_reference == 1.0
        assert not agreement.only_theirs


def test_prolif_refuses_unprotonated_structures():
    """Silently degrading on a hydrogen-free crystal structure would look like
    disagreement with this pipeline when it is really an unprepared input."""
    from plinter import validation

    assert not validation.has_hydrogens(DATA_DIR / "1kmv.pdb")
    pytest.importorskip("prolif")
    with pytest.raises(ValueError, match="no explicit hydrogens"):
        validation.prolif_residues(DATA_DIR / "1kmv.pdb", "LII")


@pytest.fixture(scope="module")
def plip_profiles():
    """PLIP's typed interaction profile for both complexes."""
    pytest.importorskip("plip")
    from plinter import validation

    profiles = {}
    for target in TARGETS:
        path = DATA_DIR / target.filename
        if not path.exists():
            pytest.skip(f"{path} missing; run `make data` first")
        profiles[target.ligand] = validation.plip_interactions(path, target.ligand)
    return profiles


def test_plip_runs_on_the_deposited_structures(plip_profiles):
    """PLIP protonates internally, so it accepts the very files ProLIF refuses.

    That is why it is the reference that actually executes here: every entry in
    this project is a crystal structure deposited without hydrogens.
    """
    from plinter import validation

    for target in TARGETS:
        assert not validation.has_hydrogens(DATA_DIR / target.filename)
        assert plip_profiles[target.ligand].residues


def test_plip_independently_reproduces_the_arg70_discrimination(plip_profiles):
    """The headline finding, confirmed by an independent implementation.

    MOT's ionised glutamate tail reaches Arg70; LII has no anion to pair with
    the guanidinium. PLIP types interactions from perceived chemistry rather
    than from this pipeline's distance rules, so agreement here is not
    circular.
    """
    assert "ARG70" in plip_profiles["MOT"].residues
    assert "ARG70" not in plip_profiles["LII"].residues


def test_plip_confirms_the_shared_glu30_anchor(plip_profiles):
    """Both antifolates carry the 2,4-diaminopyrimidine head that pairs with
    Glu30, so the anchor must appear in both profiles."""
    for ligand in ("MOT", "LII"):
        assert "GLU30" in plip_profiles[ligand].residues
        assert plip_profiles[ligand].counts.get("salt_bridge", 0) >= 1


def test_plip_finds_interaction_classes_this_pipeline_cannot(plip_profiles):
    """Asserted rather than mentioned, so the blind spot cannot quietly close.

    Waters are stripped before contact detection, so a water bridge is
    invisible in principle; pi-stacking has no geometric rule at all.
    """
    from plinter import validation

    for ligand in ("MOT", "LII"):
        counts = plip_profiles[ligand].counts
        unmodelled = [n for n in validation.UNMODELLED_BY_US if n in counts]
        assert "pi_stacking" in unmodelled
        assert "water_bridge" in unmodelled


def test_plip_places_each_pi_stack_on_a_different_phenylalanine(plip_profiles):
    """MOT's benzoyl ring stacks edge-on against Phe34, LII's ring face-on
    against Phe31. Both residues appear in both profiles, so the partner has to
    be read from the stacking interaction itself."""
    assert plip_profiles["MOT"].pi_stacking_residues == {"PHE34"}
    assert plip_profiles["LII"].pi_stacking_residues == {"PHE31"}


def test_plip_reports_fewer_residues_than_a_distance_cutoff(plip_profiles, fingerprints):
    """The two methods answer different questions and must not be conflated.

    This pipeline reports every heavy-atom contact within 5 A; PLIP reports
    only chemically typed interactions. PLIP's pocket is therefore a subset in
    size, and quoting the two counts side by side as if they measured the same
    thing would be wrong.
    """
    for target, fingerprint in zip(TARGETS, fingerprints, strict=True):
        assert len(plip_profiles[target.ligand].residues) < len(fingerprint.residues)
