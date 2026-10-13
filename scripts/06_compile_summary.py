import sys
import os
import yaml
import csv
import warnings
import math
import argparse
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg') # CRITICAL for concurrent background plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter1d
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ["PYTHONUNBUFFERED"] = "1"
warnings.filterwarnings('ignore')

print("====================================================", flush=True)
print("[PHASE 6] AUTOMATED GLOBAL DISCOVERY COMPILER (CONCURRENT)", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

# ---------------------------------------------------------------------------
# 0. CLI & DYNAMIC ROUTER
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--config", type=str, default="config.yaml", 
                    help="Name of the specific config YAML file to use.")
cli_args, unknown = parser.parse_known_args()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

if os.path.isabs(cli_args.config):
    CONFIG_PATH = cli_args.config
else:
    CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', cli_args.config))

if not os.path.exists(CONFIG_PATH):
    print(f"❌ [FATAL ERROR] Config file not found at: {CONFIG_PATH}", flush=True)
    sys.exit(1)

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

# ==========================================
# 1. BIOPHYSICAL SCALES
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

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
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
    if not sequence: return 0.0
    pKa_pos = {'K': 10.53, 'R': 12.48, 'H': 6.00, 'N_term': 9.69}
    pKa_neg = {'D': 3.65, 'E': 4.25, 'C': 8.18, 'Y': 10.07, 'C_term': 2.34}
    
    charge = 0.0
    for aa, pka in pKa_pos.items():
        count = sequence.count(aa) if aa != 'N_term' else 1
        if count > 0:
            charge += count * (10**(pka - pH)) / (1 + 10**(pka - pH))
            
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
    print(f"⚠️ [INFO] Phase 2 Metrics File not discovered at {master_csv_phase2}.", flush=True)

# -------------------------------------------------------------------
# CONCURRENT WORKER FUNCTION (WITH 4-PANEL INDIVIDUAL REPORT CARD)
# -------------------------------------------------------------------
def compile_candidate(folder):
    model_md_path = os.path.join(md_dir, folder)
    pdb_template_path = os.path.join(model_md_path, "topology_template.pdb")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    mmpbsa_csv = os.path.join(mmpbsa_out_dir, f"{folder}_mmpbsa.csv")
    
    if not os.path.exists(pdb_template_path) or not os.path.exists(nc_path):
        return None
        
    print(f"\n📊 Compiling 4-Panel Report Card for: {folder}", flush=True)
    
    try:
        raw_traj = md.load(nc_path, top=pdb_template_path)
        total_frames = raw_traj.n_frames
        
        # Windowing configuration
        start_frame = config.get('analysis_start_frame', config.get('energy_start_frame', 0))
        end_frame = config.get('analysis_end_frame', config.get('energy_end_frame', -1))
        frame_interval = config.get('analysis_stride', config.get('energy_interval', 1))
        
        if start_frame >= total_frames: 
            start_frame = 0
        if end_frame == -1 or end_frame > total_frames or end_frame <= start_frame: 
            end_frame = total_frames
            
        # 1. Filter out solvent/ions first for fast, stable processing
        protein_atoms = raw_traj.topology.select("protein")
        traj = raw_traj.atom_slice(protein_atoms)
        
        # 2. Identify chains safely
        n_chains = traj.topology.n_chains
        if n_chains < 2: return None
        rec_chain_idx, lig_chain_idx = 0, n_chains - 1

        rec_atoms = traj.topology.select(f"chainid {rec_chain_idx}")
        lig_atoms = traj.topology.select(f"chainid {lig_chain_idx}")

        # 3. PBC Fix using the robust Atom mapping
        traj.make_molecules_whole(inplace=True)
        rec_indices = set(rec_atoms.tolist())
        anchor_mols = [mol for mol in traj.topology.find_molecules() if any(a.index in rec_indices for a in mol)]
        traj.image_molecules(anchor_molecules=anchor_mols, inplace=True)
        
        # 4. Superpose receptor backbone
        rec_bb = traj.topology.select(f"chainid {rec_chain_idx} and backbone")
        traj.superpose(traj, 0, atom_indices=rec_bb)
        
        # Slice to the configured window
        traj_eq = traj[start_frame:end_frame:frame_interval]
        time_axis = np.linspace(start_frame, end_frame, traj_eq.n_frames)
        
        # --- A. Calculate RMSD ---
        lig_bb = traj.topology.select(f"chainid {lig_chain_idx} and backbone")
        rmsd_values = md.rmsd(traj_eq, traj_eq, 0, atom_indices=lig_bb) * 10.0
        mean_rmsd = np.mean(rmsd_values)
        
        # --- B. Calculate Interface Contacts (1D Profile) ---
        rec_res = [r.index for r in traj_eq.topology.chain(rec_chain_idx).residues]
        pep_res = [r.index for r in traj_eq.topology.chain(lig_chain_idx).residues]
        
        pairs = np.array([[r, p] for r in rec_res for p in pep_res])
        distances, _ = md.compute_contacts(traj_eq, contacts=pairs, scheme='closest-heavy', periodic=True)
        contact_freq = (distances < 0.45).mean(axis=0)
        
        contact_matrix = contact_freq.reshape(len(rec_res), len(pep_res))
        target_contact_profile = contact_matrix.sum(axis=1)
        
        # --- C. Calculate Dynamic H-Bonds ---
        hbonds = md.baker_hubbard(traj_eq, periodic=True)
        hbond_counts_per_frame = np.zeros(traj_eq.n_frames)
        mean_hbonds = 0.0
        
        if len(hbonds) > 0:
            da_distances = md.compute_distances(traj_eq, hbonds[:, [0,2]], periodic=True)
            dha_angles = md.compute_angles(traj_eq, hbonds, periodic=True)
            is_hbond = (da_distances < 0.35) & (dha_angles > np.deg2rad(120))
            
            rec_atom_set, lig_atom_set = set(rec_atoms), set(lig_atoms)
            inter_mask = []
            for d, h, a in hbonds:
                inter_mask.append((d in rec_atom_set and a in lig_atom_set) or (d in lig_atom_set and a in rec_atom_set))
            
            inter_mask = np.array(inter_mask)
            if np.any(inter_mask):
                inter_is_hbond = is_hbond[:, inter_mask]
                hbond_counts_per_frame = inter_is_hbond.sum(axis=1)
                mean_hbonds = np.mean(hbond_counts_per_frame)
                
        # --- D. Extract MMGBSA Energy ---
        dg_final, dg_std = 0.0, 0.0
        window_energies = []
        if os.path.exists(mmpbsa_csv):
            with open(mmpbsa_csv, 'r') as f:
                lines = f.readlines()
            start_idx = -1
            for i, line in enumerate(lines):
                if "Absolute_Frame_Index" in line or "Frame_Delta_G_kcal_mol" in line:
                    start_idx = i + 1
                    break
            if start_idx != -1:
                for line in lines[start_idx:]:
                    if not line.strip(): continue
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            energy = float(parts[1])
                            window_energies.append(energy)
                        except ValueError:
                            continue
                if window_energies:
                    dg_final = np.mean(window_energies)
                    dg_std = np.std(window_energies)

        # =====================================================================
        # GENERATE 4-PANEL INDIVIDUAL REPORT CARD
        # =====================================================================
        fig = plt.figure(figsize=(16, 10), dpi=300)
        gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)
        plot_color = '#1565C0' # Default S3-Dock Blue
        
        # Panel A: Violin Plot (Energy)
        axA = fig.add_subplot(gs[0, 0])
        if window_energies:
            parts = axA.violinplot([window_energies], showmeans=True, showextrema=True)
            for pc in parts['bodies']:
                pc.set_facecolor(plot_color)
                pc.set_alpha(0.65)
                pc.set_edgecolor('black')
            for partname in ['cmeans', 'cbars', 'cmaxes', 'cmins']:
                parts[partname].set_color('black')
            
            jitter = np.random.normal(1, 0.04, size=len(window_energies))
            axA.scatter(jitter, window_energies, color=plot_color, alpha=0.15, s=3, zorder=1)
            axA.set_xticks([1])
            axA.set_xticklabels([folder], fontweight="bold")
        axA.set_ylabel("MM/GBSA ΔG (kcal/mol)", fontweight="bold")
        axA.set_title("A. Binding Free Energy Distribution", fontweight="bold")
        axA.grid(axis='y', alpha=0.3)

        # Panel B: Backbone RMSD Profile
        axB = fig.add_subplot(gs[0, 1])
        smoothed_rmsd = gaussian_filter1d(rmsd_values, sigma=4)
        axB.plot(time_axis, rmsd_values, color=plot_color, alpha=0.3, linewidth=1)
        axB.plot(time_axis, smoothed_rmsd, color=plot_color, linewidth=2.5, label=f"Avg RMSD: {mean_rmsd:.2f} Å")
        axB.set_xlabel("Time (Frames)", fontweight="bold")
        axB.set_ylabel("Backbone RMSD (Å)", fontweight="bold")
        axB.set_title("B. Structural Stability Profile", fontweight="bold")
        axB.legend(loc="upper left")
        axB.grid(True, alpha=0.3)

        # Panel C: Active H-Bonds
        axC = fig.add_subplot(gs[1, 0])
        smoothed_hb = gaussian_filter1d(hbond_counts_per_frame, sigma=4)
        axC.plot(time_axis, hbond_counts_per_frame, color=plot_color, alpha=0.3, linewidth=1)
        axC.plot(time_axis, smoothed_hb, color=plot_color, linewidth=2.5, label=f"Avg H-Bonds: {mean_hbonds:.1f}")
        axC.set_xlabel("Time (Frames)", fontweight="bold")
        axC.set_ylabel("Active Intermolecular H-Bonds", fontweight="bold")
        axC.set_title("C. Dynamic H-Bond Network Integrity", fontweight="bold")
        axC.legend(loc="upper left")
        axC.grid(True, alpha=0.3)

        # Panel D: 1D Target Contact Profile
        axD = fig.add_subplot(gs[1, 1])
        target_positions = np.arange(len(target_contact_profile))
        axD.plot(target_positions, target_contact_profile, color=plot_color, linewidth=2.0)
        axD.fill_between(target_positions, 0, target_contact_profile, color=plot_color, alpha=0.3)
        
        active_idx = np.where(target_contact_profile > 0.05)[0]
        if len(active_idx) > 0:
            axD.set_xlim(max(0, active_idx[0] - 5), active_idx[-1] + 5)
        
        axD.set_xlabel("Target Residue Index", fontweight="bold")
        axD.set_ylabel("Total Peptide Contacts (Summed)", fontweight="bold")
        axD.set_title("D. Target Hotspot Engagement Profile", fontweight="bold")
        axD.grid(True, alpha=0.3)

        plt.suptitle(f"Individual Simulation Report: {folder}", fontsize=16, fontweight="bold")
        plt.savefig(os.path.join(summary_dir, f"{folder}_4Panel_Report.png"), dpi=300, bbox_inches="tight")
        plt.close()
        
        # --- DATA COMPILATION FOR CSV ROW ---
        base_seed = folder.replace("_DEEP", "").split("_M")[0]
        p2_data = None
        for row in phase2_metrics.values():
            if len(row) > 1 and row[1] == base_seed:
                p2_data = row
                break
                
        if p2_data:
            seq = p2_data[2]
            if seq == "N/A" or not seq: seq = ""
            compiled_row = [
                folder, p2_data[1], p2_data[2], p2_data[3], p2_data[4],
                p2_data[7], p2_data[8], p2_data[10], p2_data[11],
                round(dg_final, 2) if dg_final != 0.0 else "N/A", 
                round(dg_std, 2) if dg_std != 0.0 else "N/A",     
                round(mean_rmsd, 2), round(mean_hbonds, 1), 
                "N/A", # BSA bypassed for speed in this version
                round(calculate_native_charge(seq, 7.4), 2),     
                calculate_residue_counts(seq)["pos"], calculate_residue_counts(seq)["neg"],             
                calculate_residue_counts(seq)["polar"], calculate_residue_counts(seq)["arom"],            
                calculate_residue_counts(seq)["hydro"],          
                round(sum([HYDROPHOBICITY_SCALE.get(aa, 0.0) for aa in seq]) / len(seq), 2) if seq else 0.0,         
                round(calculate_hydrophobic_moment(seq), 3)         
            ]
            return compiled_row
        return None

    except Exception as err:
        print(f"   ❌ [ERROR] Compilation failed for {folder}: {err}", flush=True)
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
print("[PHASE 6 COMPLETE] ALL INDIVIDUAL REPORT CARDS GENERATED SUCCESSFULLY.", flush=True)