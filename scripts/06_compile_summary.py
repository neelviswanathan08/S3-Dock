import sys
import os
import yaml
import csv
import warnings
import math
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg') # CRITICAL for concurrent background plotting
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial.distance import cdist
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ["PYTHONUNBUFFERED"] = "1"
warnings.filterwarnings('ignore')

print("====================================================", flush=True)
print("[PHASE 6] AUTOMATED GLOBAL DISCOVERY COMPILER (CONCURRENT)", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
md_dir = os.path.join(run_dir, "md_simulations")
mmpbsa_out_dir = os.path.join(run_dir, "mmgbsa_results")
final_dir = os.path.join(run_dir, "top_designs")
summary_dir = os.path.join(run_dir, "final_summary")

os.makedirs(summary_dir, exist_ok=True)

if not os.path.exists(md_dir) or not os.listdir(md_dir):
    print("❌ ERROR: Discovery components missing. Cannot compile master report.", flush=True)
    sys.exit(1)

CONTACT_CUTOFF = 4.5

# ==========================================
# 1. BIOPHYSICAL SCALES & NATIVE ENGINES (NO BIOPYTHON)
# ==========================================
HYDROPHOBIC_RES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'PRO', 'GLY', 'CYS', 'CYX'}

HYDROPHOBICITY_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

EISENBERG_SCALE = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
    'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
    'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
    'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
}

def calculate_residue_counts(sequence):
    counts = {"pos": 0, "neg": 0, "polar": 0, "arom": 0, "hydro": 0}
    for aa in sequence:
        if aa in "RKH": counts["pos"] += 1
        if aa in "DE": counts["neg"] += 1
        if aa in "STNQCY": counts["polar"] += 1
        if aa in "FYW": counts["arom"] += 1
        if aa in "AILMFVWY": counts["hydro"] += 1
    return counts

def calculate_hydrophobic_moment(sequence, angle=100):
    rad_per_deg = math.pi / 180.0
    sum_sin, sum_cos = 0, 0
    for i, aa in enumerate(sequence):
        hp = EISENBERG_SCALE.get(aa, 0.0)
        current_angle = i * angle * rad_per_deg
        sum_sin += hp * math.sin(current_angle)
        sum_cos += hp * math.cos(current_angle)
    n = len(sequence) if len(sequence) > 0 else 1
    return math.sqrt(sum_sin ** 2 + sum_cos ** 2) / n

def calculate_native_charge(sequence, pH=7.4):
    """Native Python replacement for Biopython's ProteinAnalysis.charge_at_pH()"""
    if not sequence: return 0.0
    pKa_pos = {'K': 10.53, 'R': 12.48, 'H': 6.00, 'N_term': 9.69}
    pKa_neg = {'D': 3.65, 'E': 4.25, 'C': 8.18, 'Y': 10.07, 'C_term': 2.34}
    
    charge = 0.0
    # Positive groups (including N-terminus)
    for aa, pka in pKa_pos.items():
        count = sequence.count(aa) if aa != 'N_term' else 1
        if count > 0:
            charge += count * (10**(pka - pH)) / (1 + 10**(pka - pH))
            
    # Negative groups (including C-terminus)
    for aa, pka in pKa_neg.items():
        count = sequence.count(aa) if aa != 'C_term' else 1
        if count > 0:
            charge -= count / (1 + 10**(pka - pH))
            
    return charge

# -------------------------------------------------------------------
# INGEST PHASE 2 MASTER METRICS
# -------------------------------------------------------------------
phase2_metrics = {}
master_csv_phase2 = os.path.join(final_dir, f"caps_2_{run_name}_master_metrics.csv")
if os.path.exists(master_csv_phase2):
    print(f"📦 Discovered Phase 2 Master Metrics. Ingesting Rosetta & Sequence data...", flush=True)
    with open(master_csv_phase2, 'r') as p2_file:
        reader = csv.reader(p2_file)
        header = next(reader, None) 
        if header:
            for row in reader:
                if row:
                    model_id = row[0]
                    phase2_metrics[model_id] = row
else:
    print(f"⚠️ WARNING: Phase 2 Metrics File not discovered at {master_csv_phase2}.", flush=True)

# -------------------------------------------------------------------
# CONCURRENT WORKER FUNCTION
# -------------------------------------------------------------------
def compile_candidate(folder):
    model_md_path = os.path.join(md_dir, folder)
    pdb_template_path = os.path.join(model_md_path, "topology_template.pdb")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    mmpbsa_csv = os.path.join(mmpbsa_out_dir, f"{folder}_mmpbsa.csv")
    
    if not os.path.exists(pdb_template_path) or not os.path.exists(nc_path):
        return None
        
    print(f"\n📊 Compiling Advanced Biophysics Reports for Design Vector: {folder}", flush=True)
    
    try:
        raw_traj = md.load(nc_path, top=pdb_template_path)
        total_frames = raw_traj.n_frames
        
        start_frame = config.get('energy_start_frame', total_frames // 2)
        end_frame = config.get('energy_end_frame', total_frames)
        frame_interval = config.get('energy_interval', 1)
        
        if start_frame >= total_frames: start_frame = 0
        if end_frame == -1 or end_frame > total_frames or end_frame <= start_frame: end_frame = total_frames
        
        traj = raw_traj[start_frame:end_frame:frame_interval]
        traj.image_molecules(inplace=True)
        
        dry_idx = traj.topology.select('protein')
        traj_dry = traj.atom_slice(dry_idx)
        
        ligand_chain_index = traj_dry.topology.n_chains - 1
        rec_idx = traj_dry.topology.select(f'not chainid {ligand_chain_index}')
        lig_idx = traj_dry.topology.select(f'chainid {ligand_chain_index}')
        
        traj_rec = traj_dry.atom_slice(rec_idx)
        traj_lig = traj_dry.atom_slice(lig_idx)
        
        # 1. RMSD
        bb_selection = traj_dry.topology.select(f"chainid {ligand_chain_index} and backbone")
        rmsd_values = md.rmsd(traj_dry, traj_dry[0], atom_indices=bb_selection) * 10.0 
        mean_rmsd = np.mean(rmsd_values)
        
        # 2. H-Bonds
        hbond_data = md.wernet_nilsson(traj_dry)
        hbond_counts = [len(frame) for frame in hbond_data]
        mean_hbonds = np.mean(hbond_counts)
        
        # 3. BSA
        sasa_comp = np.sum(md.shrake_rupley(traj_dry), axis=1)
        sasa_rec = np.sum(md.shrake_rupley(traj_rec), axis=1)
        sasa_lig = np.sum(md.shrake_rupley(traj_lig), axis=1)
        bsa_values = (sasa_rec + sasa_lig - sasa_comp) * 100.0 
        mean_bsa = np.mean(bsa_values)
        
        # 4. Phase 5 Thermodynamics (🚨 SD FIX APPLIED HERE)
        dg_final, dg_std = 0.0, 0.0
        if os.path.exists(mmpbsa_csv):
            with open(mmpbsa_csv, 'r') as f:
                reader = csv.reader(f)
                next(reader) 
                for row in reader:
                    if len(row) >= 6:
                        dg_final = float(row[4])
                        dg_std = float(row[5]) # Directly extract the SD instead of recalculating it
        
        # 5. Data Alignment & Biophysical Native Calculations
        p2_data = phase2_metrics.get(folder, ["N/A"] * 12)
        if len(p2_data) < 12:
            p2_data = [folder] + list(p2_data[1:]) + ["N/A"] * (12 - len(p2_data))
            
        seq = p2_data[2]
        if seq == "N/A" or not seq:
            seq = ""
            
        # Native math executions
        counts = calculate_residue_counts(seq)
        kd_scores = [HYDROPHOBICITY_SCALE.get(aa, 0.0) for aa in seq]
        avg_kd = sum(kd_scores) / len(seq) if seq else 0.0
        net_charge = calculate_native_charge(seq, 7.4)
        hmoment = calculate_hydrophobic_moment(seq)
            
        compiled_row = [
            folder,                   # Design_Target
            p2_data[1],               # Seed_ID
            p2_data[2],               # Binder_Sequence
            p2_data[3],               # Seq_Length
            p2_data[4],               # MW_kDa
            p2_data[7],               # Boltz_ipTM
            p2_data[8],               # Complex_TM_Score
            p2_data[10],              # Rosetta_dG_separated
            p2_data[11],              # Phase_2_Status
            round(dg_final, 2) if dg_final != 0.0 else "N/A", 
            round(dg_std, 2) if dg_std != 0.0 else "N/A",     
            round(mean_rmsd, 2),      
            round(mean_hbonds, 1),    
            round(mean_bsa, 1),
            round(net_charge, 2),     # NEW: Native Charge
            counts["pos"],            # NEW: Positive Res
            counts["neg"],            # NEW: Negative Res
            counts["polar"],          # NEW: Polar Res
            counts["arom"],           # NEW: Aromatic Res
            counts["hydro"],          # NEW: Hydrophobic Res
            round(avg_kd, 2),         # NEW: Avg Hydrophobicity KD
            round(hmoment, 3)         # NEW: Hydrophobic Moment uH
        ]
        
        # -------------------------------------------------------------------
        # GENERATE HIGH-RES FIGURES
        # -------------------------------------------------------------------
        all_chains = list(traj_dry.topology.chains)
        ligand_chain = all_chains[ligand_chain_index]
        k_residues = list(ligand_chain.residues)
        k_res_labels = [f"{res.name}{res.resSeq}" for res in k_residues]
        
        # A. Binding Hotspots Histogram
        contact_counts = []
        fibril_atom_indices = [a.index for a in traj_rec.topology.atoms]
        for res_k in k_residues:
            atoms_k = [a.index for a in res_k.atoms]
            res_contacts_per_frame = []
            for f_idx in range(traj_dry.n_frames):
                c_k = traj_dry.xyz[f_idx, atoms_k]
                c_f = traj_dry.xyz[f_idx, fibril_atom_indices]
                if len(c_k) > 0 and len(c_f) > 0:
                    dists = cdist(c_k, c_f) * 10.0 
                    num_contacts = np.sum(dists < CONTACT_CUTOFF)
                    res_contacts_per_frame.append(num_contacts)
                else:
                    res_contacts_per_frame.append(0)
            contact_counts.append(np.mean(res_contacts_per_frame))

        plt.figure(figsize=(12, 6))
        colors = ['tab:orange' if res.name in HYDROPHOBIC_RES else 'tab:blue' for res in k_residues]
        bars = plt.bar(k_res_labels, contact_counts, color=colors, edgecolor='black', zorder=3)
        for bar in bars:
            height = bar.get_height()
            if height > 1.0:
                plt.text(bar.get_x() + bar.get_width()/2., height + 1.0, f'{height:.1f}', ha='center', va='bottom', fontsize=9, color='black', fontweight='bold', rotation=90)
        
        plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
        plt.title(f"Biochemical Binding Hotspots: ({folder})", fontsize=14, pad=15)
        plt.xlabel("Peptide Residue", fontsize=12, fontweight='bold', labelpad=10)
        plt.ylabel(f"Average Atomic Contacts (< {CONTACT_CUTOFF} Å)", fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.ylim(0, max(contact_counts) * 1.15 if len(contact_counts) > 0 else 10)
        
        hydrophobic_patch = mpatches.Patch(color='tab:orange', label='Hydrophobic')
        hydrophilic_patch = mpatches.Patch(color='tab:blue', label='Polar / Charged')
        plt.legend(handles=[hydrophobic_patch, hydrophilic_patch], loc='upper right')
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, f"{folder}_binding_hotspots.png"), dpi=300)
        plt.close()

        # B. Interaction Network Map
        rec_chains = [c for c in all_chains if c.index != ligand_chain.index]
        fibril_residues, fib_res_labels = [], []
        for chain_obj in rec_chains:
            c_id = getattr(chain_obj, 'id', f"C{chain_obj.index}")
            for r in chain_obj.residues:
                fibril_residues.append(r)
                fib_res_labels.append(f"{r.name}{r.resSeq}_{c_id}")
                         
        global_dist_matrix = np.zeros((len(k_residues), len(fibril_residues)))
        for y, res_k in enumerate(k_residues):
            atoms_k = [a.index for a in res_k.atoms]
            for x, res_f in enumerate(fibril_residues):
                atoms_f = [a.index for a in res_f.atoms]
                frame_dists = []
                for f_idx in range(traj_dry.n_frames):
                    c_k = traj_dry.xyz[f_idx, atoms_k]
                    c_f = traj_dry.xyz[f_idx, atoms_f]
                    if len(c_k) > 0 and len(c_f) > 0:
                        frame_dists.append(np.min(cdist(c_k, c_f)))
                    else:
                        frame_dists.append(np.nan)
                global_dist_matrix[y, x] = np.nanmean(frame_dists) * 10.0 

        active_y, active_x = np.where(global_dist_matrix <= CONTACT_CUTOFF)
        unique_active_x = np.unique(active_x) 
        fig, ax = plt.subplots(figsize=(12, 10))

        for i, y_idx in enumerate(range(len(k_residues))):
            y_pos = len(k_residues) - i
            node_color = 'tab:orange' if k_residues[y_idx].name in HYDROPHOBIC_RES else 'tab:blue'
            ax.scatter(0, y_pos, color=node_color, s=150, zorder=3, edgecolor='black')
            ax.text(-0.05, y_pos, k_res_labels[y_idx], ha='right', va='center', fontsize=11, fontweight='bold')

        if len(unique_active_x) > 0:
            for j, x_idx in enumerate(unique_active_x):
                y_pos = len(k_residues) - (j * (len(k_residues) / max(1, len(unique_active_x) - 1))) if len(unique_active_x) > 1 else len(k_residues)/2
                node_color = 'tab:orange' if fibril_residues[x_idx].name in HYDROPHOBIC_RES else 'tab:blue'
                ax.scatter(1, y_pos, color=node_color, s=150, zorder=3, edgecolor='black')
                ax.text(1.05, y_pos, fib_res_labels[x_idx], ha='left', va='center', fontsize=10, fontweight='bold')
                
                connected_peptides = np.where(global_dist_matrix[:, x_idx] <= CONTACT_CUTOFF)[0]
                for pep_idx in connected_peptides:
                    pep_y_pos = len(k_residues) - pep_idx
                    dist_val = global_dist_matrix[pep_idx, x_idx]
                    line_width = max(0.5, 5.0 - (dist_val - 2.0))
                    ax.plot([0, 1], [pep_y_pos, y_pos], color='gray', alpha=0.4, linewidth=line_width, zorder=1)
        else:
            ax.text(0.5, len(k_residues)/2, "NO CONTACTS FOUND BELOW CUTOFF", ha='center', va='center', fontsize=12, color='red', fontweight='bold')

        ax.set_xlim(-0.5, 1.5)
        ax.axis('off')
        plt.legend(handles=[hydrophobic_patch, hydrophilic_patch], loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2)
        plt.title(f"Biochemical Interaction Network: {folder}\n(Observed Contacts < {CONTACT_CUTOFF}Å)", pad=30, fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, f"{folder}_interaction_network.png"), dpi=300)
        plt.close()

        # C. Time-Series Stability Dashboard
        time_axis = np.arange(len(rmsd_values)) * frame_interval
        fig, ax1 = plt.subplots(figsize=(10, 5))
        color = 'tab:gray'
        ax1.set_xlabel('Sampled Simulation Frame Index', fontsize=10)
        ax1.set_ylabel('Backbone RMSD (Å)', color=color, fontsize=10)
        ax1.plot(time_axis, rmsd_values, color=color, linewidth=2, label='Peptide Backbone RMSD')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, alpha=0.3)
        ax2 = ax1.twinx()
        color = 'tab:green'
        ax2.set_ylabel('Active Inter-Chain H-Bonds', color=color, fontsize=10)
        ax2.plot(time_axis, hbond_counts, color=color, linewidth=1.5, alpha=0.7, linestyle=':', label='Hydrogen Bonds')
        ax2.tick_params(axis='y', labelcolor=color)
        plt.title(f"Dynamic Structural Trace Dashboard: {folder}", fontsize=12, pad=10)
        fig.tight_layout()
        plt.savefig(os.path.join(summary_dir, f"{folder}_stability_profile.png"), dpi=300)
        plt.close()
        
        print(f"   ✅ [SUCCESS] Graphic processing and biophysics mapping complete for {folder}.", flush=True)
        return compiled_row

    except Exception as err:
        print(f"   ❌ [METRIC RUNTIME WARNING] Visual compilation failed for {folder}: {err}", flush=True)
        return None

# -------------------------------------------------------------------
# EXECUTE CONCURRENT GENERATION
# -------------------------------------------------------------------
master_data = []
folders_to_process = [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]

if folders_to_process:
    max_workers = min(os.cpu_count() or 4, len(folders_to_process))
    print(f"\n[INFO] Launching Phase 6 multiprocessing engine with {max_workers} active workers...", flush=True)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(compile_candidate, folder): folder for folder in folders_to_process}
        for future in as_completed(futures):
            result = future.result()
            if result:
                master_data.append(result)

# -------------------------------------------------------------------
# EXPORT MASTER REPORT
# -------------------------------------------------------------------
output_report_csv = os.path.join(summary_dir, "pipeline_summary_report.csv")
with open(output_report_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Design_Target", "Seed_ID", "Binder_Sequence", "Seq_Length", "MW_kDa", 
        "Boltz_ipTM", "Complex_TM_Score", "Rosetta_dG_kcal_mol", "Phase_2_Status",
        "MMGBSA_Delta_G_kcal_mol", "MMGBSA_Std_Dev", "Mean_Backbone_RMSD_A", "Mean_Hydrogen_Bonds", "Buried_Surface_Area_A2",
        "Net_Charge_pH7.4", "Pos_Res", "Neg_Res", "Polar_Res", "Aromatic_Res", "Hydrophobic_Res", "Avg_Hydrophobicity_KD", "Hydrophobic_Moment_uH"
    ])
    writer.writerows(master_data)

print("\n----------------------------------------------------", flush=True)
print(f"[SUCCESS] Merged Master CSV Report logged at: {output_report_csv}", flush=True)
print("[PHASE 6 COMPLETE] ALL ASSETS INTEGRATED SUCCESSFULLY! THE DISCOVERY PIPELINE IS COMPLETE.", flush=True)