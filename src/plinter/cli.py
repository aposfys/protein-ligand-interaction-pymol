"""Command-line entry point for the hDHFR antifolate contact analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import chemistry, validation
from .compare import Comparison, build_fingerprint
from .contacts import CONTACT_CUTOFF, find_contacts, hydrogen_bonds
from .report import (
    format_comparison,
    format_top_contacts,
    write_contacts_csv,
    write_rows_csv,
    write_summary_json,
)
from .structures import COFACTOR, TARGETS, download_structure, load_structure


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="plinter",
        description=(
            "Compare the binding-site interactions of two antifolate inhibitors "
            "in ternary complexes of human dihydrofolate reductase."
        ),
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/pdb"),
        help="Directory holding (or receiving) the PDB files. Default: data/pdb",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV/JSON/figure output. Default: results",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=CONTACT_CUTOFF,
        help=f"Heavy-atom contact cutoff in angstroms. Default: {CONTACT_CUTOFF}",
    )
    parser.add_argument(
        "--keep-cofactor",
        action="store_true",
        help="Include the NADPH cofactor in the binding-site environment.",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip figure generation (avoids the matplotlib dependency).",
    )
    parser.add_argument(
        "--no-chemistry",
        action="store_true",
        help="Skip the RDKit ligand profiling stage.",
    )
    parser.add_argument(
        "--no-validation",
        action="store_true",
        help="Skip the comparison against the depositors' SITE records.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    excluded = () if args.keep_cofactor else (COFACTOR,)

    fingerprints = []
    for target in TARGETS:
        path = download_structure(target.pdb_id, args.data_dir / target.filename)
        structure = load_structure(path)
        contacts = find_contacts(
            structure,
            ligand_resname=target.ligand,
            exclude_resnames=excluded,
            cutoff=args.cutoff,
        )
        fingerprint = build_fingerprint(target.pdb_id, target.ligand, contacts)
        fingerprints.append(fingerprint)

        write_contacts_csv(
            contacts,
            args.results_dir / f"contacts_{target.pdb_id.lower()}_{target.ligand}.csv",
        )

        print(f"\n## {target.pdb_id} — {target.ligand} ({target.inhibitor_name})")
        print(
            f"{len(contacts)} contacts within {args.cutoff} Å across "
            f"{len(fingerprint.residues)} residues; "
            f"{len(hydrogen_bonds(contacts))} hydrogen bonds.\n"
        )
        print(format_top_contacts(contacts))

    comparison = Comparison(first=fingerprints[0], second=fingerprints[1])
    write_rows_csv(comparison.residue_table(), args.results_dir / "residue_comparison.csv")
    write_summary_json(fingerprints, comparison, args.results_dir / "summary.json")

    print("\n## Comparative analysis\n")
    print(format_comparison(comparison))

    extras: dict[str, object] = {}

    if not args.no_validation:
        print("\n## Validation against the depositors' annotation\n")
        agreements = []
        for target, fingerprint in zip(TARGETS, fingerprints, strict=True):
            reference = validation.site_record_residues(
                args.data_dir / target.filename, target.ligand
            )
            if not reference:
                print(f"{target.pdb_id}: no SITE record for {target.ligand}")
                continue
            agreement = validation.compare_to_reference(
                sorted(fingerprint.residues),
                reference,
                f"{target.pdb_id} SITE records",
            )
            agreements.append(agreement.as_dict())
            print(validation.format_agreement(agreement))
            print()
        extras["validation"] = agreements
        write_rows_csv(agreements, args.results_dir / "validation.csv")

    if not args.no_chemistry:
        print("\n## Ligand physicochemical profile\n")
        profiles = [chemistry.profile(target.ligand, args.data_dir) for target in TARGETS]
        rows = []
        for entry in profiles:
            groups = chemistry.ionisable_groups(entry)
            print(
                f"{entry.ligand}: MW {entry.molecular_weight}, "
                f"cLogP {entry.logp:+.2f}, TPSA {entry.tpsa} A^2, "
                f"HBD {entry.hydrogen_bond_donors}, HBA {entry.hydrogen_bond_acceptors}, "
                f"rotatable {entry.rotatable_bonds}, "
                f"carboxylates {groups['carboxylic_acid']}"
            )
            rows.append({**entry.as_dict(), **groups})
        write_rows_csv(rows, args.results_dir / "ligand_properties.csv")
        extras["chemistry"] = {
            "profiles": [entry.as_dict() for entry in profiles],
            "comparison": chemistry.compare(*profiles),
        }

    if extras:
        import json

        (args.results_dir / "extras.json").write_text(
            json.dumps(extras, indent=2) + "\n", encoding="utf-8"
        )

    if not args.no_figures:
        try:
            from .plots import plot_closest_approach, plot_interaction_profile
        except ImportError:
            print("\nmatplotlib is not installed; skipping figures.", file=sys.stderr)
        else:
            plot_closest_approach(comparison, args.results_dir / "closest_approach.png")
            plot_interaction_profile(comparison, args.results_dir / "interaction_profile.png")

    print(f"\nResults written to {args.results_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
