"""Retrieval and parsing of the PDB entries analysed in this project."""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


@dataclass(frozen=True)
class Target:
    """One protein-ligand complex to analyse."""

    pdb_id: str
    ligand: str
    inhibitor_name: str

    @property
    def filename(self) -> str:
        return f"{self.pdb_id.lower()}.pdb"


# Two ternary complexes of human dihydrofolate reductase (hDHFR) with NADPH and
# a lipophilic antifolate. The shared target and cofactor make the two
# inhibitors directly comparable.
TARGETS: tuple[Target, ...] = (
    Target("1HFR", "MOT", "classical furo[2,3-d]pyrimidine antifolate, L-glutamate tail"),
    Target("1KMV", "LII", "SRI-9662, lipophilic pyrido[2,3-d]pyrimidine antifolate"),
)

# NADPH is a cofactor rather than part of the binding site under study, and is
# excluded from the protein environment following the original assignment.
COFACTOR = "NDP"


def download_structure(pdb_id: str, destination: Path) -> Path:
    """Download a PDB entry from RCSB unless it is already present."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return destination

    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id.upper())
    with urllib.request.urlopen(url, timeout=60) as response:
        destination.write_bytes(response.read())
    return destination


def load_structure(path: Path) -> Structure:
    """Parse a PDB file, keeping only the first model."""
    parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem.upper(), str(path))
