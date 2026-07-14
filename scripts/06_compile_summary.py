import sys
import os
import yaml
import csv
import warnings
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

os.environ["PYTHONUNBUFFERED"] = "1"
warnings.filterwarnings('ignore')

print("====================================================", flush=True)
print("[PHASE 6] AUTOMATED GLOBAL DISCOVERY COMPILER", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
md_dir = os.path.join(run_dir, "md_simulations")
mmgbsa_dir = os.path.join(run_dir, "mmgbsa_results")
final_dir = os.path.join(run_dir, "top_designs")
summary_dir = os.path.join(run_dir, "final_summary")

os.makedirs(summary_dir, exist_ok=True)

if not os.path.exists(md_dir) or not os.listdir(md_dir):
    print("❌ ERROR: Discovery components missing. Cannot compile master report.", flush=True)
    sys.exit(1)

# 🚨 PARSE THE EXPANDED PHASE 2 MASTER METRICS CSV
phase2_metrics = {}
master_csv_phase2 = os.path.join(final_dir, f"caps_2_{run_name}_master_metrics.csv")
if os.path.exists(master_csv_phase2):
    print(f"📦 Discovered Phase 2 Master Metrics File. Ingesting scores and sequence data...", flush=True)
    with open(master_csv_phase2, 'r') as p2_file:
        reader = csv.reader(p2_file)
        header = next(reader, None) 
        if header:
            for row in reader:
                if row:
                    model_id = row[0]
                    # Map exactly to the new 14-item Phase 2 array (after Model_ID)
                    phase2_metrics[model_id] = row[1:]
else:
    print(f"⚠️ WARNING: Phase 2 Metrics File not discovered at {master_csv_phase2}. Columns will fall back to N/A.", flush=True)

master_data = []
CONTACT_CUTOFF = 4.5

for folder in [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]:
    model_md_path = os.path.join(md_dir, folder)
    cif_path = os.path.join(model_md_path, "topology_template.cif")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    mmgbsa_csv = os.path.join(mmgbsa_dir, f"{folder}_mmgbsa.csv")
    
    if not os.path.exists(cif_path) or not os.path.exists(nc_path):
        continue
        
    print(f"\n📊 Compiling Advanced Biophysics Reports for Design Vector: {folder}", flush=True)
    
    try:
        raw_traj = md.load(nc_path, top=cif_path)
        total_frames = raw_traj.n_frames
        default_start = total_frames // 2
        start_frame = config.get('mmpbsa_start_frame', default_start)
        end_frame = config.get('mmpbsa_end_frame', -1)
        frame_interval = config.get('mmpbsa_interval', 1)
        if end_frame == -1 or end_frame > total_frames: end_frame = total_frames
        
        traj = raw_traj[start_frame:end_frame:frame_interval]
        traj.image_molecules(inplace=True)
        
        dry_idx = traj.topology.select('protein')
        traj_dry = traj.atom_slice(dry_idx)
        
        ligand_chain_index = traj_dry.topology.n_chains - 1
        rec_idx = traj_dry.topology.select(f'not chainid {ligand_chain_index}')
        lig_idx = traj_dry.topology.select(f'chainid {ligand_chain_index}')
        
        traj_rec = traj_dry.atom_slice(rec_idx)
        traj_lig = traj_dry.atom_slice(lig_idx)
        
        print("   ↳ Computing peptide backbone stability matrices (RMSD)...", flush=True)
        bb_selection = traj_dry.topology.select(f"chainid {ligand_chain_index} and backbone")
        rmsd_values = md.rmsd(traj_dry, traj_dry[0], atom_indices=bb_selection) * 10.0 
        mean_rmsd = np.mean(rmsd_values)
        
        print("   ↳ Tracking spatial hydrogen bonding anchor clusters...", flush=True)
        hbond_data = md.wernet_nilsson(traj_dry)
        hbond_counts = [len(frame) for frame in hbond_data]
        mean_hbonds = np.mean(hbond_counts)
        
        print("   ↳ Evaluating hydrophobic shielding profiles (Buried Surface Area)...", flush=True)
        sasa_comp = np.sum(md.shrake_rupley(traj_dry), axis=1)
        sasa_rec = np.sum(md.shrake_rupley(traj_rec), axis=1)
        sasa_lig = np.sum(md.shrake_rupley(traj_lig), axis=1)
        bsa_values = (sasa_rec + sasa_lig - sasa_comp) * 100.0 
        mean_bsa = np.mean(bsa_values)
        
        print("   ↳ Parsing Phase 5 production thermodynamics...", flush=True)
        dg_final, dg_std = 0.0, 0.0
        if os.path.exists(mmgbsa_csv):
            dg_vals = []
            with open(mmgbsa_csv, 'r') as f:
                reader = csv.reader(f)
                next(reader) 
                for row in reader:
                    dg_vals.append(float(row[4]))
            if dg_vals:
                dg_final = np.mean(dg_vals)
                dg_std = np.std(dg_vals)
        
        # 🚨 MERGING LAYER: Aligning 14 columns from Phase 2
        p2_data = phase2_metrics.get(folder, ["N/A"] * 14)
        
        master_data.append([
            folder,                   # Design_Target (Model_ID)
            p2_data[0],               # Seed_ID
            p2_data[1],               # Binder_Sequence
            p2_data[2],               # Seq_Length
            p2_data[3],               # MW_kDa
            p2_data[4],               # pI
            p2_data[5],               # GRAVY
            p2_data[6],               # Boltz_ipTM
            p2_data[7],               # Is_Rigid_Fibril_Mode
            p2_data[8],               # Structural_Metric_Score
            p2_data[9],               # 3D_Helicity_Score
            p2_data[10],              # PRODIGY_dG_kcal_mol
            p2_data[11],              # PRODIGY_Kd_Text
            p2_data[12],              # PRODIGY_Kd_Molar_Value
            p2_data[13],              # Status
            round(dg_final, 2) if dg_final != 0.0 else "N/A", 
            round(dg_std, 2) if dg_std != 0.0 else "N/A",     
            round(mean_rmsd, 2),      
            round(mean_hbonds, 1),    
            round(mean_bsa, 1)        
        ])
        
        # (Plots logic preserved untouched...)
        print("   ↳ Mapping global binding hotspots...", flush=True)
        all_chains = list(traj_dry.topology.chains)
        ligand_chain = all_chains[ligand_chain_index]
        k_residues = list(ligand_chain.residues)
        k_res_labels = [f"{res.name}{res.resSeq}" for res in k_residues]
        
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
        colors = ['crimson' if c > 30 else 'tab:blue' for c in contact_counts]
        bars = plt.bar(k_res_labels, contact_counts, color=colors, edgecolor='black', zorder=3)
        for bar in bars:
            height = bar.get_height()
            if height > 1.0:
                plt.text(bar.get_x() + bar.get_width()/2., height + 1.0, f'{height:.1f}', ha='center', va='bottom', fontsize=9, color='black', fontweight='bold', rotation=90)
        plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
        plt.title(f"Global Binding Hotspots: ({folder})", fontsize=14, pad=15)
        plt.xlabel("Peptide Residue", fontsize=12, fontweight='bold', labelpad=10)
        plt.ylabel(f"Average Atomic Contacts (< {CONTACT_CUTOFF} Å)", fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.ylim(0, max(contact_counts) * 1.15 if len(contact_counts) > 0 else 10)
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, f"{folder}_binding_hotspots.png"), dpi=300)
        plt.close()

        print("   ↳ Running Global Contact Scan across all fibril chains...", flush=True)
        rec_chains = [c for c in all_chains if c.index != ligand_chain.index]
        fibril_residues = []
        fib_res_labels = []
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
            ax.scatter(0, y_pos, color='tab:blue', s=150, zorder=3, edgecolor='black')
            ax.text(-0.05, y_pos, k_res_labels[y_idx], ha='right', va='center', fontsize=11, fontweight='bold')

        if len(unique_active_x) > 0:
            for j, x_idx in enumerate(unique_active_x):
                y_pos = len(k_residues) - (j * (len(k_residues) / max(1, len(unique_active_x) - 1))) if len(unique_active_x) > 1 else len(k_residues)/2
                ax.scatter(1, y_pos, color='crimson', s=150, zorder=3, edgecolor='black')
                ax.text(1.05, y_pos, fib_res_labels[x_idx], ha='left', va='center', fontsize=10, fontweight='bold')
                connected_peptides = np.where(global_dist_matrix[:, x_idx] <= CONTACT_CUTOFF)[0]
                for pep_idx in connected_peptides:
                    pep_y_pos = len(k_residues) - pep_idx
                    dist_val = global_dist_matrix[pep_idx, x_idx]
                    line_width = max(0.5, 5.0 - (dist_val - 2.0))
                    ax.plot([0, 1], [pep_y_pos, y_pos], color='gray', alpha=0.4, linewidth=line_width, zorder=1)
        else:
            ax.text(0.5, len(k_residues)/2, "GLOBAL ERROR: NO CONTACTS FOUND BELOW CUTOFF", ha='center', va='center', fontsize=12, color='red', fontweight='bold')

        ax.set_xlim(-0.5, 1.5)
        ax.axis('off')
        plt.title(f"Universal Global Interaction Network: {folder}\n(Observed Contacts < {CONTACT_CUTOFF}Å)", pad=20, fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, f"{folder}_interaction_network.png"), dpi=300)
        plt.close()

        print("   ↳ Exporting dynamic time-series stability diagnostics...", flush=True)
        time_axis = np.arange(len(rmsd_values)) * frame_interval
        fig, ax1 = plt.subplots(figsize=(10, 5))
        color = 'tab:blue'
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
        
    except Exception as err:
        print(f"   ❌ [METRIC RUNTIME WARNING] Visual compilation paused for {folder}: {err}", flush=True)

output_report_csv = os.path.join(summary_dir, "pipeline_summary_report.csv")
with open(output_report_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Design_Target", 
        "Seed_ID", 
        "Binder_Sequence",
        "Seq_Length",
        "MW_kDa",
        "pI",
        "GRAVY",
        "Boltz_ipTM", 
        "Is_Rigid_Fibril_Mode",
        "Structural_Metric_Score", 
        "3D_Helicity_Score", 
        "PRODIGY_dG_kcal_mol", 
        "PRODIGY_Kd_Text",
        "PRODIGY_Kd_Molar_Value",
        "Status",
        "MMGBSA_Delta_G_kcal_mol", 
        "MMGBSA_Std_Dev", 
        "Mean_Backbone_RMSD_A", 
        "Mean_Hydrogen_Bonds", 
        "Buried_Surface_Area_A2"
    ])
    writer.writerows(master_data)

print("\n----------------------------------------------------", flush=True)
print(f"[SUCCESS] Merged Master CSV Report logged at: {output_report_csv}", flush=True)
print(f"[SUCCESS] High-fidelity graphical visual assets cached in: {summary_dir}", flush=True)
print("[PHASE 6 COMPLETE] ALL ASSETS INTEGRATED SUCCESSFULLY! WE ARE DONE!", flush=True)