import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.ndimage import gaussian_filter1d
import mdtraj as md

# =====================================================================
# 1. CONFIGURATION
# =====================================================================
TARGET_NAME = "MDM2_vs_Control_Comparison_Boltz" 
BASE_DIR = "/home/aicd3-2/Desktop/nviswanathan/S3-Dock/results/MAIN/RESULTS_Manuscript_S3-DOCK/BOLTZ/MDM2_COMPARISON/"

REPLICATE_SUBPATHS = {
    "Negative_Control": [
        "MDM2_Negative/MDM2_Negative_rep_1",
        "MDM2_Negative/MDM2_Negative_rep_2",
        "MDM2_Negative/MDM2_Negative_rep_3",
    ],
    "Positive_Control": [
        "MDM2_positive/MDM2_positive_rep_1",
        "MDM2_positive/MDM2_positive_rep_2",
        "MDM2_positive/MDM2_positive_rep_3",
    ],
    "S3-DOCK": [
        "MDM2_OUR/MDM2_OURS_rep_1",
        "MDM2_OUR/MDM2_OURS_rep_2",
        "MDM2_OUR/MDM2_OURS_rep_3",
    ],
}

REPLICATES = {cond: [os.path.join(BASE_DIR, sub) for sub in subs] for cond, subs in REPLICATE_SUBPATHS.items()}
OUTPUT_FIG = f"{TARGET_NAME}_4Panel_Master.png"
COLORS = {"Negative_Control": "#AB0B0B", "Positive_Control": "#24A355", "S3-DOCK": "#6A0E86"}


# =====================================================================
# 2. DYNAMIC EXTRACTION ENGINE (PBC-SAFE & INDEX-MAPPED)
# =====================================================================
def extract_md_data(folder):
    """Extracts RMSD, 1D Target Contacts, and Time-Series H-Bonds safely."""
    md_base = os.path.join(folder, "md_simulations")
    if not os.path.exists(md_base): return None, None, None
        
    subdirs = [os.path.join(md_base, d) for d in os.listdir(md_base) if os.path.isdir(os.path.join(md_base, d))]
    if not subdirs: return None, None, None
        
    actual_md_dir = subdirs[0]
    top_file = os.path.join(actual_md_dir, "topology_template.pdb")
    traj_file = os.path.join(actual_md_dir, "trajectory.nc")
    
    if not os.path.exists(top_file) or not os.path.exists(traj_file):
        return None, None, None
        
    print(f"  -> Loading & unwrapping trajectory from {os.path.basename(actual_md_dir)}...")
    raw_traj = md.load(traj_file, top=top_file)
    
    # Filter out solvent/ions for fast, stable processing
    protein_atoms = raw_traj.topology.select("protein")
    traj = raw_traj.atom_slice(protein_atoms)
    
    n_chains = traj.topology.n_chains
    if n_chains < 2: return None, None, None
    rec_chain_idx, lig_chain_idx = 0, n_chains - 1

    rec_atoms = traj.topology.select(f"chainid {rec_chain_idx}")
    lig_atoms = traj.topology.select(f"chainid {lig_chain_idx}")
    
    # PBC Fix: Map to proper MDTraj Atom objects to avoid the ValueError
    traj.make_molecules_whole(inplace=True)
    rec_indices = set(rec_atoms.tolist())
    anchor_mols = [mol for mol in traj.topology.find_molecules() if any(a.index in rec_indices for a in mol)]
    traj.image_molecules(anchor_molecules=anchor_mols, inplace=True)
    
    rec_bb_idx = traj.topology.select(f"chainid {rec_chain_idx} and backbone")
    traj.superpose(traj, 0, atom_indices=rec_bb_idx)
    
    # 1. RMSD Extraction
    pep_bb_idx = traj.topology.select(f"chainid {lig_chain_idx} and backbone")
    rmsd = md.rmsd(traj, traj, 0, atom_indices=pep_bb_idx) * 10 
    
    traj_eq = traj[len(traj)//2:] 
    
    # 2. 1D Target Contact Profile
    rec_res = [r.index for r in traj_eq.topology.chain(rec_chain_idx).residues]
    pep_res = [r.index for r in traj_eq.topology.chain(lig_chain_idx).residues]
    
    if len(rec_res) == 0 or len(pep_res) == 0:
        return rmsd, None, None 
        
    pairs = np.array([[r, p] for r in rec_res for p in pep_res])
    distances, _ = md.compute_contacts(traj_eq, contacts=pairs, scheme='closest-heavy', periodic=True)
    contact_freq = (distances < 0.45).mean(axis=0)
    
    contact_matrix = contact_freq.reshape(len(rec_res), len(pep_res))
    target_contact_profile = contact_matrix.sum(axis=1)
        
    # 3. Dynamic Time-Series H-Bonds
    hbonds = md.baker_hubbard(traj_eq, periodic=True)
    hbond_counts_per_frame = np.zeros(traj_eq.n_frames)
    
    if len(hbonds) > 0:
        da_distances = md.compute_distances(traj_eq, hbonds[:, [0,2]], periodic=True)
        dha_angles = md.compute_angles(traj_eq, hbonds, periodic=True)
        is_hbond = (da_distances < 0.35) & (dha_angles > np.deg2rad(120))
        
        rec_atom_set = set(rec_atoms)
        lig_atom_set = set(lig_atoms)
        
        inter_mask = []
        for d, h, a in hbonds:
            inter_mask.append((d in rec_atom_set and a in lig_atom_set) or (d in lig_atom_set and a in rec_atom_set))
        
        inter_mask = np.array(inter_mask)
        if np.any(inter_mask):
            inter_is_hbond = is_hbond[:, inter_mask]
            hbond_counts_per_frame = inter_is_hbond.sum(axis=1)

    return rmsd, target_contact_profile, hbond_counts_per_frame

def extract_mmgbsa(folder, cond, verbose=True):
    """Pulls energy data, explicitly skipping the summary header to target frame data."""
    path = os.path.join(folder, "mmgbsa_results", "**", "*.csv")
    matches = glob.glob(path, recursive=True)
    if not matches:
        if verbose: print(f"  [WARN] No MMGBSA csv found in {folder}")
        return []
        
    try:
        df = pd.read_csv(matches[0], skiprows=3)
        target_col = "Frame_Delta_G_kcal_mol"
        
        if target_col in df.columns:
            vals = pd.to_numeric(df[target_col], errors='coerce').dropna().tolist()
            return vals[len(vals)//2:]
        else:
            return []
            
    except Exception as e:
        if verbose: print(f"  [ERROR] [{os.path.basename(folder)}] Failed to parse CSV: {e}")
        return []

# =====================================================================
# 3. RUN EXTRACTION & PLOT
# =====================================================================
print("--- Extracting Physics Data from Raw Folders ---")
fig = plt.figure(figsize=(18, 12), dpi=300)
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

all_mmgbsa = {cond: [] for cond in REPLICATES}
rep_mmgbsa_means = {cond: [] for cond in REPLICATES} # Added for Mean of Means Table Math
all_rmsd = {cond: [] for cond in REPLICATES}
all_contacts = {cond: [] for cond in REPLICATES}
all_hbonds = {cond: [] for cond in REPLICATES}

for cond, folders in REPLICATES.items():
    print(f"\nScanning {cond}...")
    for f in folders:
        energies = extract_mmgbsa(f, cond, verbose=True)
        if energies: 
            all_mmgbsa[cond].extend(energies)
            rep_mmgbsa_means[cond].append(np.mean(energies)) # Store the per-replicate mean
        
        rmsd, t_profile, hb_time = extract_md_data(f)
        if rmsd is not None: all_rmsd[cond].append(rmsd)
        if t_profile is not None: all_contacts[cond].append(t_profile)
        if hb_time is not None: all_hbonds[cond].append(hb_time)

# ---- Plot A: MMGBSA Violin Plot (Last 50ns) ----
axA = fig.add_subplot(gs[0, 0])
v_data = [all_mmgbsa[c] for c in REPLICATES if all_mmgbsa[c]]
v_labels = [c for c in REPLICATES if all_mmgbsa[c]]
if v_data:
    parts = axA.violinplot(v_data, showmeans=True, showextrema=True, points=1000, bw_method='scott')
    
    for i, pc in enumerate(parts["bodies"]): 
        color = COLORS.get(v_labels[i], "gray")
        pc.set_facecolor(color)
        pc.set_alpha(0.65)
        pc.set_edgecolor("black")
        pc.set_linewidth(1.5)
        
    parts['cmeans'].set_color('black')
    parts['cmeans'].set_linewidth(2)
    parts['cbars'].set_color('black')
    parts['cmaxes'].set_color('black')
    parts['cmins'].set_color('black')
    
    for i, data in enumerate(v_data):
        jitter = np.random.normal(i + 1, 0.05, size=len(data))
        axA.scatter(jitter, data, color=COLORS.get(v_labels[i], "gray"), alpha=0.15, s=3, zorder=1)

    axA.set_xticks(range(1, len(v_labels) + 1))
    axA.set_xticklabels(v_labels, rotation=10, fontsize=11, fontweight="bold")
    axA.set_ylabel("MM/GBSA ΔG (kcal/mol)", fontweight="bold")
    axA.set_title("A. Binding Free Energy (Last 50ns Equilibration)", fontweight="bold")
    axA.grid(axis='y', alpha=0.3)

# ---- Plot B: RMSD Mean +/- SD (Full 100ns) ----
axB = fig.add_subplot(gs[0, 1])
for cond in REPLICATES:
    runs = all_rmsd[cond]
    if not runs: continue
    min_len = min(len(r) for r in runs)
    matrix = np.vstack([r[:min_len] for r in runs])
    time = np.linspace(0, 100, min_len) 
    
    mean, sd = matrix.mean(axis=0), matrix.std(axis=0)
    smoothed_mean = gaussian_filter1d(mean, sigma=5)
    smoothed_sd = gaussian_filter1d(sd, sigma=5)
    
    axB.plot(time, mean, color=COLORS[cond], linewidth=1, alpha=0.3)
    axB.plot(time, smoothed_mean, label=cond, color=COLORS[cond], linewidth=2.5)
    axB.fill_between(time, smoothed_mean - smoothed_sd, smoothed_mean + smoothed_sd, color=COLORS[cond], alpha=0.15)

axB.set_xlabel("Time (ns)")
axB.set_ylabel("Peptide Backbone RMSD (Å)")
axB.set_title("B. Structural Stability (Smoothed, n=3)", fontweight="bold")
axB.legend(loc="upper left")

# ---- Plot C: Dynamic Time-Series H-Bond Trace (Last 50ns) ----
axC = fig.add_subplot(gs[1, 0])
for cond in REPLICATES:
    runs = all_hbonds[cond]
    if not runs: continue
    min_len = min(len(r) for r in runs)
    matrix = np.vstack([r[:min_len] for r in runs])
    time = np.linspace(50, 100, min_len)
    
    mean = matrix.mean(axis=0)
    smoothed_mean = gaussian_filter1d(mean, sigma=5)
    
    axC.plot(time, mean, color=COLORS[cond], linewidth=1, alpha=0.2)
    axC.plot(time, smoothed_mean, label=f"{cond} (Avg: {mean.mean():.1f})", color=COLORS[cond], linewidth=2.5)

axC.set_xlabel("Time (ns)")
axC.set_ylabel("Active Intermolecular H-Bonds")
axC.set_title("C. Dynamic H-Bond Network Integrity (Last 50ns)", fontweight="bold")
axC.grid(True, alpha=0.3)
axC.legend(loc="upper left")

# ---- Plot D: Target Hotspot Engagement 1D Line Profile ----
axD = fig.add_subplot(gs[1, 1])
global_min_idx = float('inf')
global_max_idx = float('-inf')

for cond in REPLICATES:
    if not all_contacts[cond]: continue
    
    mean_profile = np.mean(all_contacts[cond], axis=0)
    target_positions = np.arange(len(mean_profile))
    
    axD.plot(target_positions, mean_profile, linewidth=2.0, alpha=0.85, color=COLORS[cond], label=cond)
    axD.fill_between(target_positions, 0, mean_profile, color=COLORS[cond], alpha=0.3)
    
    active_idx = np.where(mean_profile > 0.05)[0]
    if len(active_idx) > 0:
        global_min_idx = min(global_min_idx, active_idx[0])
        global_max_idx = max(global_max_idx, active_idx[-1])

if global_min_idx < float('inf') and global_max_idx > float('-inf'):
    axD.set_xlim(max(0, global_min_idx - 5), global_max_idx + 5)

axD.set_xlabel("Target Residue Index")
axD.set_ylabel("Total Peptide Contacts (Summed)")
axD.set_title("D. Target Hotspot Engagement Profile (Last 50ns)", fontweight="bold")
axD.grid(True, alpha=0.3)
axD.legend(loc="upper left")

plt.suptitle(f"{TARGET_NAME}: Master Analysis Figure (n=3)", fontsize=15, fontweight="bold")
plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches="tight")
print(f"\n✅ SUCCESS: Extracted data & saved unified master figure -> {OUTPUT_FIG}")

# =====================================================================
# 4. GENERATE SUMMARY TABLE (FOR MANUSCRIPT)
# =====================================================================
print("\n--- Generating Paper-Ready Summary Table ---")
table_data = []

for cond in REPLICATES.keys():
    row = {"Condition": cond}
    
    # 1. Aggregate MM/GBSA (Strict Mean of Means for proper variance)
    if rep_mmgbsa_means[cond]:
        mean_dg = np.mean(rep_mmgbsa_means[cond])
        std_dg = np.std(rep_mmgbsa_means[cond])
        row["MM/GBSA ΔG (kcal/mol)"] = f"{mean_dg:.2f} ± {std_dg:.2f}"
    else:
        row["MM/GBSA ΔG (kcal/mol)"] = "N/A"
        
    # 2. Aggregate RMSD (Mean ± SD of the full trajectory across all replicates)
    if all_rmsd[cond]:
        flat_rmsd = np.concatenate(all_rmsd[cond])
        mean_rmsd = np.mean(flat_rmsd)
        std_rmsd = np.std(flat_rmsd)
        row["Peptide RMSD (Å)"] = f"{mean_rmsd:.2f} ± {std_rmsd:.2f}"
    else:
        row["Peptide RMSD (Å)"] = "N/A"
        
    # 3. Aggregate H-Bonds (Mean ± SD of the last 50ns across all replicates)
    if all_hbonds[cond]:
        flat_hbonds = np.concatenate(all_hbonds[cond])
        mean_hb = np.mean(flat_hbonds)
        std_hb = np.std(flat_hbonds)
        row["Active H-Bonds"] = f"{mean_hb:.1f} ± {std_hb:.1f}"
    else:
        row["Active H-Bonds"] = "N/A"

    # 4. Aggregate Plot D Metrics (Target Hotspot Engagement)
    if all_contacts[cond]:
        # Average the 1D profile across replicates first
        mean_profile = np.mean(all_contacts[cond], axis=0) 
        
        total_mass = np.sum(mean_profile)
        engaged_res = np.sum(mean_profile > 0.1)
        peak_strength = np.max(mean_profile) if len(mean_profile) > 0 else 0.0
        
        row["Total Contact Mass"] = f"{total_mass:.1f}"
        row["Engaged Residues (>0.1)"] = f"{int(engaged_res)}"
        row["Peak Anchor Strength"] = f"{peak_strength:.2f}"
    else:
        row["Total Contact Mass"] = "N/A"
        row["Engaged Residues (>0.1)"] = "N/A"
        row["Peak Anchor Strength"] = "N/A"
        
    table_data.append(row)

# Convert to DataFrame, save as CSV, and print to console as Markdown
df_summary = pd.DataFrame(table_data)
csv_output_path = os.path.join(BASE_DIR, f"{TARGET_NAME}_Summary_Table.csv")
df_summary.to_csv(csv_output_path, index=False)

print(f"✅ SUCCESS: Saved summary table CSV -> {csv_output_path}\n")
print("📋 Copy this directly into your manuscript's Results section:\n")
print(df_summary.to_markdown(index=False))
print("\n" + "="*70)