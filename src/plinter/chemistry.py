"""Physicochemical profiling of the two inhibitors.

The contact analysis says *where* each ligand binds. This module says *what
each ligand is*, and the two together explain the binding modes: a classical
antifolate carries an ionised glutamate tail that reaches a charged subsite,
while a lipophilic analogue trades those contacts for permeability.

Descriptors come from RDKit over the canonical SMILES in the PDB Chemical
Component Dictionary, rather than from the deposited coordinates, so bond
orders and protonation are the curated ones rather than inferred from geometry.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

CHEMCOMP_API = "https://data.rcsb.org/rest/v1/core/chemcomp/{ligand}"

# Preference order for the SMILES flavour: the stereochemistry-bearing CACTVS
# string first, then any other canonical form.
SMILES_PREFERENCE = ("SMILES_CANONICAL", "SMILES")

# Lipinski's rule of five, for reference rather than as a pass/fail gate.
LIPINSKI = {
    "molecular_weight": 500,
    "logp": 5,
    "hydrogen_bond_donors": 5,
    "hydrogen_bond_acceptors": 10,
}


@dataclass(frozen=True)
class LigandProfile:
    """RDKit descriptors for one PDB chemical component."""

    ligand: str
    name: str
    formula: str
    smiles: str
    molecular_weight: float
    logp: float  # Crippen cLogP: lipophilicity
    tpsa: float  # topological polar surface area
    hydrogen_bond_donors: int
    hydrogen_bond_acceptors: int
    rotatable_bonds: int
    aromatic_rings: int
    formal_charge: int
    fraction_sp3: float
    heavy_atoms: int

    def as_dict(self) -> dict:
        return asdict(self)

    @property
    def lipinski_violations(self) -> list[str]:
        return [name for name, limit in LIPINSKI.items() if getattr(self, name) > limit]


def fetch_chemcomp(ligand: str, cache: Path) -> dict:
    """Cache the RCSB chemical component record for a ligand."""
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists() and cache.stat().st_size > 0:
        return json.loads(cache.read_text(encoding="utf-8"))

    with urllib.request.urlopen(
        CHEMCOMP_API.format(ligand=ligand.upper()), timeout=60
    ) as response:
        payload = response.read()
    cache.write_bytes(payload)
    return json.loads(payload)


def _smiles_from(record: dict) -> str:
    descriptors = record.get("pdbx_chem_comp_descriptor", [])
    by_type: dict[str, str] = {}
    for entry in descriptors:
        kind = entry.get("type", "")
        if kind in SMILES_PREFERENCE and kind not in by_type:
            by_type[kind] = entry["descriptor"]

    for kind in SMILES_PREFERENCE:
        if kind in by_type:
            return by_type[kind]
    raise ValueError(f"no SMILES in chemical component {record.get('rcsb_id')}")


def profile(ligand: str, cache_dir: Path) -> LigandProfile:
    """Compute descriptors for one ligand by its PDB chemical component ID."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")

    record = fetch_chemcomp(ligand, cache_dir / f"chemcomp_{ligand.upper()}.json")
    chem_comp = record.get("chem_comp", {})
    smiles = _smiles_from(record)

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse the SMILES for {ligand}: {smiles}")

    return LigandProfile(
        ligand=ligand.upper(),
        name=(chem_comp.get("name") or "").title(),
        formula=chem_comp.get("formula", ""),
        smiles=smiles,
        molecular_weight=round(Descriptors.MolWt(molecule), 2),
        logp=round(Crippen.MolLogP(molecule), 2),
        tpsa=round(rdMolDescriptors.CalcTPSA(molecule), 2),
        hydrogen_bond_donors=Lipinski.NumHDonors(molecule),
        hydrogen_bond_acceptors=Lipinski.NumHAcceptors(molecule),
        rotatable_bonds=Lipinski.NumRotatableBonds(molecule),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(molecule),
        formal_charge=Chem.GetFormalCharge(molecule),
        fraction_sp3=round(rdMolDescriptors.CalcFractionCSP3(molecule), 3),
        heavy_atoms=molecule.GetNumHeavyAtoms(),
    )


def ionisable_groups(profile_: LigandProfile) -> dict[str, int]:
    """Count acidic and basic groups, which drive the charged contacts.

    A carboxylate is deprotonated at physiological pH and is what reaches a
    charged arginine subsite; a neutral analogue has nothing to put there.
    """
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(profile_.smiles)
    patterns = {
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "primary_aromatic_amine": "[NX3;H2][c]",
        "amide": "[NX3][CX3](=[OX1])",
    }
    return {
        name: len(molecule.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        for name, smarts in patterns.items()
    }


def compare(first: LigandProfile, second: LigandProfile) -> dict[str, object]:
    """Contrast two ligands on the properties that separate their binding modes."""
    return {
        "ligands": [first.ligand, second.ligand],
        "lipophilicity_difference_logp": round(second.logp - first.logp, 2),
        "polar_surface_difference_tpsa": round(second.tpsa - first.tpsa, 2),
        "formal_charge": {
            first.ligand: first.formal_charge,
            second.ligand: second.formal_charge,
        },
        "rotatable_bonds": {
            first.ligand: first.rotatable_bonds,
            second.ligand: second.rotatable_bonds,
        },
        "lipinski_violations": {
            first.ligand: first.lipinski_violations,
            second.ligand: second.lipinski_violations,
        },
    }
