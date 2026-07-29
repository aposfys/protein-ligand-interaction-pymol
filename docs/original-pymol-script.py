from pymol import cmd

def measure_interactions(structure):
    # Select ligand and nearby residues, excluding ligand, water, and NDP
    ligand = "MOT" if structure == "1HFR" else "LII"
    cmd.select(f"{ligand}_lig", f"resn {ligand}")
    cmd.select("nearby", f"byres ({ligand}_lig around 5.0) and not (resn {ligand} or resn HOH or resn NDP)")
    
    # List to store measurements
    measurements = []
    
    # Measure distances
    model_lig = cmd.get_model(f"{ligand}_lig")
    model_nearby = cmd.get_model("nearby")
    for atom1 in model_lig.atom:
        for atom2 in model_nearby.atom:
            distance = cmd.get_distance(f"id {atom1.id}", f"id {atom2.id}")
            
            if distance <= 5.0:
                # Check for possible hydrogen bonds
                possible_hbond = "Yes" if (atom1.name[0] in ["O", "N"] and 
                                           atom2.name[0] in ["O", "N"] and 
                                           distance <= 3.5) else "No"
                measurements.append({
                    f"{ligand} Atoms": f"{atom1.name}",
                    "Protein Atoms": f"{atom2.resn}{atom2.resi}/{atom2.name}",
                    "Distance (Å)": round(distance, 2),
                    "Possible Hydrogen Bonds": possible_hbond
                })
    
    # Sort measurements by distance
    measurements.sort(key=lambda x: x["Distance (Å)"])
    
    # Print results
    print(f"\nResults for {structure}:")
    print(f"{ligand} Atoms | Protein Atoms | Distance (Å) | Possible Hydrogen Bonds")
    print("-" * 75)
    for m in measurements:
        print(f"{m[f'{ligand} Atoms']:<10} | {m['Protein Atoms']:<17} | {m['Distance (Å)']:<13} | {m['Possible Hydrogen Bonds']}")

# Execute the function for both structures
measure_interactions("1HFR")
measure_interactions("1KMV")
