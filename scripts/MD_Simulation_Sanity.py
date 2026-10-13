import os
import numpy as np
import mdtraj as md

# =====================================================================
# 1. MASTER CONFIGURATION (All Three Studies)
# =====================================================================
STUDIES = {
    "{Comparison_1}": {
        "base_dir": "{add_base_directory_path_here}",
        "paths": {
            "{Control_1}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "{Control_2}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "S3-DOCK": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"]
        }
    },
    "{Comparison_2}": {
        "base_dir": "{add_base_directory_path_here}",
        "paths": {
            "{Control_1}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "{Control_2}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "S3-DOCK": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"]
        }
    },
    "{Comparison_3}": {
        "base_dir": "{add_base_directory_path_here}",
        "paths": {
            "{Control_1}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "{Control_2}": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"],
            "S3-DOCK": ["{add_directory_path_here_1}", "{add_directory_path_here_2}", "{add_directory_path_here_3}"]
        }
    }
}

# =====================================================================
# 2. PBC SANITY CHECK ENGINE
# =====================================================================
def run_pbc_sanity_check(folder):
    """Finds the traj/top files and calculates unbiased PBC min distances."""
    md_base = os.path.join(folder, "md_simulations")
    if not os.path.exists(md_base): 
        return "ERROR: No md_simulations folder found."
        
    subdirs = [os.path.join(md_base, d) for d in os.listdir(md_base) if os.path.isdir(os.path.join(md_base, d))]
    if not subdirs: 
        return "ERROR: No inner MD folder found."
        
    actual_md_dir = subdirs[0]
    
    # Try .pdb first, fallback to .cif
    top_file = os.path.join(actual_md_dir, "topology_template.pdb")
    if not os.path.exists(top_file):
        top_file = os.path.join(actual_md_dir, "topology_template.cif")
        
    traj_file = os.path.join(actual_md_dir, "trajectory.nc")
    
    if not os.path.exists(top_file) or not os.path.exists(traj_file):
        return "ERROR: Missing topology or trajectory file."

    # Load trajectory and isolate protein
    try:
        traj = md.load(traj_file, top=top_file)
    except Exception as e:
        return f"ERROR loading trajectory: {e}"

    dry_idx = traj.topology.select("protein")
    traj_dry = traj.atom_slice(dry_idx)
    
    n_chains = traj_dry.topology.n_chains
    if n_chains < 2: return "ERROR: Less than 2 chains found."
    
    # Select Receptor and Ligand
    rec_idx = traj_dry.topology.select(f"not chainid {n_chains - 1}")
    lig_idx = traj_dry.topology.select(f"chainid {n_chains - 1}")
    
    box = traj_dry.unitcell_lengths
    if box is None:
        return "ERROR: No unit cell info in this trajectory — cannot apply minimum image convention."
        
    min_dists = []
    # Calculate PBC-corrected minimum distance for every frame
    for i in range(traj_dry.n_frames):
        rec_xyz = traj_dry.xyz[i, rec_idx]
        lig_xyz = traj_dry.xyz[i, lig_idx]
        
        diff = rec_xyz[:, None, :] - lig_xyz[None, :, :]
        diff -= box[i] * np.round(diff / box[i])   # Minimum image convention
        
        d = np.linalg.norm(diff, axis=-1)
        min_dists.append(d.min() * 10)  # nm -> Å

    min_dists = np.array(min_dists)
    
    min_val = min_dists.min()
    max_val = min_dists.max()
    fly_away_frac = (min_dists > 10).mean() * 100 # percentage > 10 Å
    
    status = "✅ PASS" if fly_away_frac == 0.0 else "❌ WARN (Peptide drifting)"
    return f"{status} | Range: {min_val:.1f}–{max_val:.1f} Å | Frames >10Å: {fly_away_frac:.1f}%"

# =====================================================================
# 3. EXECUTION OVER ALL STUDIES
# =====================================================================
print("=========================================================")
print("          MASTER PBC-CORRECTED SANITY CHECK              ")
print("=========================================================\n")

for study_name, config in STUDIES.items():
    print(f"► STUDY: {study_name.replace('_', ' ').upper()}")
    base_dir = config["base_dir"]
    
    for cond, subpaths in config["paths"].items():
        print(f"  Condition: {cond}")
        for sub in subpaths:
            full_path = os.path.join(base_dir, sub)
            rep_name = os.path.basename(sub)
            result = run_pbc_sanity_check(full_path)
            print(f"    [{rep_name}] -> {result}")
    print("-" * 57)

print("Sanity check complete")