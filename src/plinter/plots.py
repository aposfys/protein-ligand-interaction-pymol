"""Figures summarising the two binding-site fingerprints."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .compare import Comparison  # noqa: E402

INTERACTION_COLOURS = {
    "salt bridge": "#b4413c",
    "hydrogen bond": "#1f77b4",
    "hydrophobic": "#e8a33d",
    "van der Waals": "#b8c4cc",
}


FIRST_COLOUR = "#1f77b4"
SECOND_COLOUR = "#e8a33d"


def plot_closest_approach(comparison: Comparison, path: Path) -> Path:
    """Dumbbell chart of each ligand's closest approach to every pocket residue.

    A dumbbell keeps the eye on the *gap* between the two ligands rather than on
    bar heights, and reads correctly on an axis that does not start at zero.
    """
    first, second = comparison.first, comparison.second
    first_closest = first.closest_approach
    second_closest = second.closest_approach

    residues = sorted(
        first.residues | second.residues,
        key=lambda r: min(
            first_closest.get(r, float("inf")), second_closest.get(r, float("inf"))
        ),
        reverse=True,
    )

    fig, ax = plt.subplots(figsize=(8, max(5, len(residues) * 0.28)))

    for y, residue in enumerate(residues):
        a = first_closest.get(residue)
        b = second_closest.get(residue)
        if a is not None and b is not None:
            ax.plot([a, b], [y, y], color="#c8ced4", linewidth=1.6, zorder=1)
        if a is not None:
            ax.scatter(a, y, s=46, color=FIRST_COLOUR, zorder=2)
        if b is not None:
            ax.scatter(b, y, s=46, color=SECOND_COLOUR, zorder=2)

    ax.axvline(3.5, color="#b4413c", linestyle="--", linewidth=1, zorder=0)
    ax.text(
        3.45, len(residues) - 0.4, "hydrogen-bond cutoff",
        ha="right", va="top", fontsize=8, color="#b4413c", rotation=90,
    )

    ax.set_yticks(range(len(residues)))
    ax.set_yticklabels(residues, fontsize=9)
    ax.set_ylim(-0.8, len(residues) - 0.2)
    ax.set_xlabel("Closest heavy-atom approach (Å)")
    ax.set_title(
        "Binding-site residues of human DHFR\ncontacted by each antifolate", fontsize=11
    )
    ax.scatter([], [], s=46, color=FIRST_COLOUR, label=f"{first.ligand} ({first.pdb_id})")
    ax.scatter(
        [], [], s=46, color=SECOND_COLOUR, label=f"{second.ligand} ({second.pdb_id})"
    )
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#eef1f3", zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_interaction_profile(comparison: Comparison, path: Path) -> Path:
    """Stacked bars counting each interaction type per ligand."""
    fingerprints = [comparison.first, comparison.second]
    types = list(INTERACTION_COLOURS)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    bottoms = [0, 0]
    for interaction in types:
        counts = [
            sum(1 for c in fp.contacts if c.interaction == interaction)
            for fp in fingerprints
        ]
        ax.bar(
            [f"{fp.ligand}\n{fp.pdb_id}" for fp in fingerprints],
            counts,
            bottom=bottoms,
            label=interaction,
            color=INTERACTION_COLOURS[interaction],
        )
        bottoms = [b + c for b, c in zip(bottoms, counts)]

    ax.set_ylabel("Heavy-atom contacts within 5 Å")
    ax.set_title("Interaction-type composition")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path
