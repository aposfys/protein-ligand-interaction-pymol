"""Render publication figures of both hDHFR binding sites in PyMOL.

Run headless from the repository root:

    pymol -cq pymol/render_binding_sites.py

Writes one PNG per complex plus a superposed overlay into ``results/figures/``.
"""

import os

from pymol import cmd

# Carbon colours match the matplotlib figures produced by plinter.plots.
COMPLEXES = [
    {"pdb_id": "1hfr", "ligand": "MOT", "carbon": "0x1F77B4"},
    {"pdb_id": "1kmv", "ligand": "LII", "carbon": "0xE8A33D"},
]
COFACTOR = "NDP"
POCKET_RADIUS = 5.0
FIGURE_DIR = os.path.join("results", "figures")
DATA_DIR = os.path.join("data", "pdb")


def setup_scene():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_transparency", 0.55)
    cmd.set("label_size", 16)
    cmd.set("label_color", "black")
    cmd.set("dash_color", "grey40")
    cmd.set("dash_gap", 0.3)
    cmd.set("dash_width", 2.0)


def load_complex(pdb_id):
    """Load a local copy if present, otherwise fetch from RCSB."""
    local = os.path.join(DATA_DIR, f"{pdb_id}.pdb")
    if os.path.exists(local):
        cmd.load(local, pdb_id)
    else:
        cmd.fetch(pdb_id, name=pdb_id, async_=0)
    cmd.remove(f"{pdb_id} and solvent")
    cmd.remove(f"{pdb_id} and hydro")


def style_binding_site(pdb_id, ligand, carbon_colour):
    """Show the pocket residues as sticks and dash every polar contact."""
    pocket = f"{pdb_id}_pocket"
    lig = f"{pdb_id}_lig"

    cmd.select(lig, f"{pdb_id} and resn {ligand}")
    cmd.select(
        pocket,
        f"byres ({lig} around {POCKET_RADIUS}) and {pdb_id} "
        f"and polymer and not resn {ligand}+{COFACTOR}",
    )

    cmd.hide("everything", pdb_id)
    cmd.show("cartoon", f"{pdb_id} and polymer")
    cmd.color("grey80", f"{pdb_id} and polymer")

    cmd.show("sticks", f"{pocket} and not (name C+N+O)")
    cmd.color("palecyan", pocket)
    cmd.util.cnc(pocket)

    cmd.show("sticks", lig)
    cmd.color(carbon_colour, f"{lig} and elem C")
    cmd.util.cnc(lig)

    cmd.distance(f"{pdb_id}_hbonds", lig, pocket, 3.5, mode=2)
    cmd.label(f"{pocket} and name CA", '"%s%s" % (resn, resi)')

    cmd.orient(lig)
    cmd.zoom(lig, 4.0)


def render(name):
    os.makedirs(FIGURE_DIR, exist_ok=True)
    cmd.png(os.path.join(FIGURE_DIR, f"{name}.png"), width=1600, height=1200, dpi=300, ray=1)


def render_individual_sites():
    for entry in COMPLEXES:
        setup_scene()
        load_complex(entry["pdb_id"])
        style_binding_site(entry["pdb_id"], entry["ligand"], entry["carbon"])
        render(f"{entry['pdb_id']}_{entry['ligand']}_binding_site")


def render_overlay():
    """Superpose both complexes to compare the two inhibitors in one pocket."""
    setup_scene()
    for entry in COMPLEXES:
        load_complex(entry["pdb_id"])

    reference, mobile = COMPLEXES[0]["pdb_id"], COMPLEXES[1]["pdb_id"]
    rms = cmd.align(f"{mobile} and polymer", f"{reference} and polymer")[0]
    print(f"Cα alignment RMSD {reference} vs {mobile}: {rms:.2f} A")

    cmd.hide("everything")
    cmd.show("cartoon", f"{reference} and polymer")
    cmd.color("grey85", f"{reference} and polymer")

    for entry in COMPLEXES:
        lig = f"{entry['pdb_id']}_lig"
        cmd.select(lig, f"{entry['pdb_id']} and resn {entry['ligand']}")
        cmd.show("sticks", lig)
        cmd.color(entry["carbon"], f"{lig} and elem C")
        cmd.util.cnc(lig)

    shared = f"byres ({reference}_lig around {POCKET_RADIUS}) and {reference} and polymer"
    cmd.select("shared_pocket", shared)
    cmd.show("sticks", "shared_pocket and not (name C+N+O)")
    cmd.color("palecyan", "shared_pocket")
    cmd.util.cnc("shared_pocket")

    cmd.orient(f"{reference}_lig or {mobile}_lig")
    cmd.zoom(f"{reference}_lig or {mobile}_lig", 4.0)
    render("overlay_MOT_vs_LII")


if __name__ == "pymol":
    render_individual_sites()
    render_overlay()
    print(f"Figures written to {FIGURE_DIR}/")
