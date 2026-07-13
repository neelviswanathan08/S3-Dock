import sys
import os
import yaml
import re
import csv
import time
import shutil
import subprocess
import json
import warnings
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Superimposer import Superimposer
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

print("====================================================")
print("🚀 S3-DOCK: FACTORY RUN RUNNING...")
print("====================================================")

# 1. Setup absolute paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
final_dir = os.path.join(run_dir, "top_designs")
fasta_bridge_path = os.path.join(run_dir, "library.fasta")
master_csv = os.path.join(final_dir, f"caps_2_{run_name}_master_metrics.csv")

target_cif = os.path.abspath(os.path.join(SCRIPT_DIR, "..", config['target_cif']))
binder_id = config['binder_chain_id']
target_chains = list(config['target_chains_and_sequences'].keys())

def get_chain_centroids(model, chains):
    centroids = {}
    coords_cache = {}
    for cid in chains:
        if cid in model:
            atoms = [res['CA'] for res in model[cid] if 'CA' in res]
            if atoms:
                coords_cache[cid] = atoms
                centroids[cid] = np.mean([a.coord for a in atoms], axis=0)
    return centroids, coords_cache

def check_smart_holistic_similarity(pred_cif_path, ref_cif_path):
    try:
        parser = MMCIFParser(QUIET=True)
        ref_model = parser.get_structure("ref", ref_cif_path)[0]
        pred_model = parser.get_structure("pred", pred_cif_path)[0]
        
        ref_centroids, ref_atoms_cache = get_chain_centroids(ref_model, target_chains)
        pred_centroids, pred_atoms_cache = get_chain_centroids(pred_model, target_chains)
        
        if len(ref_centroids) != len(target_chains) or len(pred_centroids) != len(target_chains):
            return 999.0, 999
            
        ref_cids = list(ref_centroids.keys())
        pred_cids = list(pred_centroids.keys())
        
        best_rms = 999.0
        best_bad_atoms = 999
        
        for trial_p_cid in pred_cids:
            sup_trial = Superimposer()
            sup_trial.set_atoms(ref_atoms_cache[ref_cids[len(ref_cids)//2]], pred_atoms_cache[trial_p_cid]) 
            rot_trial, tran_trial = sup_trial.rotran
            
            rotated_p_centroids = {}
            for p_cid in pred_cids:
                rotated_p_centroids[p_cid] = np.dot(pred_centroids[p_cid], rot_trial) + tran_trial
                
            cost_matrix = np.zeros((len(target_chains), len(target_chains)))
            for i, r_cid in enumerate(ref_cids):
                for j, p_cid in enumerate(pred_cids):
                    cost_matrix[i, j] = np.linalg.norm(ref_centroids[r_cid] - rotated_p_centroids[p_cid])
                    
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            all_ref_atoms = []
            all_pred_atoms = []
            
            for r_idx, p_idx in zip(row_ind, col_ind):
                r_cid = ref_cids[r_idx]
                p_cid = pred_cids[p_idx]
                for r_atom, p_atom in zip(ref_atoms_cache[r_cid], pred_atoms_cache[p_cid]):
                    all_ref_atoms.append(r_atom)
                    all_pred_atoms.append(p_atom)
                    
            sup_final = Superimposer()
            sup_final.set_atoms(all_ref_atoms, all_pred_atoms)
            
            rot_f, tran_f = sup_final.rotran
            bad_atoms = 0
            
            for r_atom, p_atom in zip(all_ref_atoms, all_pred_atoms):
                t_coord = np.dot(p_atom.coord, rot_f) + tran_f
                if np.linalg.norm(r_atom.coord - t_coord) > config['dist_threshold']:
                    bad_atoms += 1
                    
            if sup_final.rms < best_rms:
                best_rms = sup_final.rms
                best_bad_atoms = bad_atoms
                
        return best_rms, best_bad_atoms
    except Exception:
        return 999.0, 999

def calculate_3d_helicity_score(pred_cif_path, chain_id=binder_id):
    try:
        parser = MMCIFParser(QUIET=True)
        model = parser.get_structure("pred", pred_cif_path)[0]
        if chain_id not in model: return 0.0
        ca_atoms = [r['CA'].coord for r in model[chain_id].get_residues() if 'CA' in r]
        if len(ca_atoms) < 5: return 0.0
        helical_bonds = 0
        total_checks = len(ca_atoms) - 4
        for i in range(total_checks):
            dist = np.linalg.norm(ca_atoms[i] - ca_atoms[i+4])
            if 5.3 <= dist <= 7.2: helical_bonds += 1
        return round(helical_bonds / total_checks, 3)
    except: return 0.0

def run_prodigy_scoring(cif_path, model_id):
    try:
        receptor_str = ",".join(target_chains)
        prodigy_bin = os.path.join(os.path.dirname(sys.executable), "prodigy")
        prodigy_cmd = f"{prodigy_bin} {cif_path} --selection {receptor_str} {binder_id}"
        result = subprocess.run(prodigy_cmd, shell=True, capture_output=True, text=True)
        delta_g, kd = 999.0, "N/A"
        for line in result.stdout.split('\n'):
            if "predicted binding affinity" in line.lower():
                try: delta_g = float(line.split()[-1])
                except: pass
            if "predicted dissociation constant" in line.lower(): kd = line.split()[-1]
        return delta_g, kd
    except: return 999.0, "N/A"

def parse_fasta(path):
    seqs = []
    current_id, current_seq = None, []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id: seqs.append((current_id, "".join(current_seq)))
                current_id = re.sub(r'[^\w\-_]', '_', line[1:].split()[0])
                current_seq = []
            else: current_seq.append(line.strip())
        if current_id: seqs.append((current_id, "".join(current_seq)))
    return seqs

def format_time(seconds):
    hours, remainder = divmod(int(seconds), 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours}h {minutes}m {secs}s"

if not os.path.exists(fasta_bridge_path):
    print(f"❌ ERROR: Could not find bridge library at {fasta_bridge_path}. Stage 1 failed.")
    exit(1)

jobs = parse_fasta(fasta_bridge_path)

if not os.path.exists(master_csv):
    os.makedirs(final_dir, exist_ok=True)
    with open(master_csv, "w", newline="") as f:
        csv.writer(f).writerow(["Model_ID", "Seed_ID", "Boltz_ipTM", "Global_Backbone_RMSD", "Dislocated_Atoms", "3D_Helicity_Score", "PRODIGY_dG_kcal_mol", "PRODIGY_Kd_M"])

print(" Loading machine learning libraries (this takes a few seconds)...", flush=True)
print(" Libraries loaded. Starting Discovery Queue...\n", flush=True)

start_time = time.time()

for idx, (seq_id, sequence) in enumerate(jobs, 1):
    
    # Calculate ETA for the UI
    if idx == 1:
        eta_str = "Calculating... (Waiting for first fold to calibrate speed)"
    else:
        avg_time = (time.time() - start_time) / (idx - 1)
        eta = avg_time * (len(jobs) - idx + 1)
        eta_str = f"~{format_time(eta)} remaining for entire run"

    print(f" [{idx}/{len(jobs)}] DISCOVERY QUEUE: {seq_id}", flush=True)
    print(f"    Sequence : {sequence}", flush=True)
    print(f"    Length   : {len(sequence)} AA", flush=True)
    print(f"    Est. Time: {eta_str}", flush=True)
    print(f"     Status  : Booting Boltz-2 Deep Learning Engine (Folding in progress... this takes time!)", flush=True)
    
    design_dir = os.path.join(run_dir, seq_id)
    os.makedirs(design_dir, exist_ok=True)
    yaml_filename = os.path.join(design_dir, f"{seq_id}.yaml")
    
    yaml_lines = ["version: 1", "sequences:"]
    for chain, seq in config['target_chains_and_sequences'].items():
        yaml_lines.append(f"  - protein:\n      id: {chain}\n      sequence: '{seq}'")
    yaml_lines.append(f"  - protein:\n      id: {binder_id}\n      sequence: '{sequence}'")
    yaml_lines.append("\ntemplates:")
    for chain in target_chains: yaml_lines.append(f"  - cif: '{target_cif}'\n    chain_id: '{chain}'\n    template_id: '{chain}'")
    
    if 'pocket_contacts' in config and config['pocket_contacts']:
        yaml_lines.append("\nconstraints:\n  - pocket:\n" + f"      binder: {binder_id}\n      contacts:")
        for contact in config['pocket_contacts']: yaml_lines.append(f"        - [{contact[0]}, {contact[1]}]")
        yaml_lines.append(f"      max_distance: {config.get('max_distance_threshold', 5.0)}")

    with open(yaml_filename, "w") as f:
        f.write("\n".join(yaml_lines))
        
    boltz_bin = os.path.join(os.path.dirname(sys.executable), "boltz")
    samples = config.get('samples_per_seed', 1)
    
    boltz_cmd = f"{boltz_bin} predict {yaml_filename} --use_msa_server --use_potentials --out_dir {design_dir} --recycling_steps 10 --diffusion_samples {samples} --override"
    
    # Run Boltz silently but capture output to catch errors
    process = subprocess.run(boltz_cmd, shell=True, capture_output=True, text=True)
    
    if process.returncode != 0:
        print(f"\n [CRITICAL BOLTZ CRASH] Boltz failed to execute for {seq_id}!\n")
        print(f"--- ERROR LOG START ---\n{process.stderr}\n--- ERROR LOG END ---", flush=True)
        continue
    
    generated_models = []
    # 🚨 SEED-LEVEL SMART RESUME: Check if Boltz already ran for this specific seed
    for root, _, files in os.walk(design_dir):
        for file in files:
            if file.endswith(".cif") and file != os.path.basename(target_cif):
                generated_models.append(os.path.join(root, file))
    
    valid_candidates = []
    print(f"    Folding Complete! Evaluating {len(generated_models)} generated models for {seq_id}:", flush=True)
    
    for best_cif_path in generated_models:
        filename = os.path.basename(best_cif_path)
        model_num_match = re.search(r'model_(\d+)', filename)
        model_num = model_num_match.group(1) if model_num_match else "0"
        model_id = f"{seq_id}_M{model_num}"
        
        iptm_score = 0.0
        for root, _, files in os.walk(design_dir):
            for file in files:
                if file.endswith(".json") and "confidence" in file.lower():
                    try:
                        with open(os.path.join(root, file), 'r') as jf: data = json.load(jf)
                        if "models" in data and len(data["models"]) > int(model_num): iptm_score = float(data["models"][int(model_num)].get("iptm", 0.0))
                        elif f"model_{model_num}" in data: iptm_score = float(data[f"model_{model_num}"].get("iptm", 0.0))
                    except: pass

        global_rmsd, bad_atoms = check_smart_holistic_similarity(best_cif_path, target_cif)
        if global_rmsd > config['rmsd_cutoff'] or bad_atoms > config['max_bad_atoms']:
            continue
            
        helicity_score = calculate_3d_helicity_score(best_cif_path)
        dG, Kd = run_prodigy_scoring(best_cif_path, model_id)
        
        valid_candidates.append({
            'path': best_cif_path, 'id': model_id, 'iptm': iptm_score,
            'rmsd': global_rmsd, 'bad_atoms': bad_atoms, 'helicity': helicity_score,
            'dG': dG, 'Kd': Kd
        })

    if valid_candidates:
        valid_candidates.sort(key=lambda x: x['dG'])
        top_design = valid_candidates[0]
        
        print(f"    TOP MODEL: {top_design['id']} | dG: {top_design['dG']} kcal/mol | ipTM: {top_design['iptm']:.2f} | Kd: {top_design['Kd']}")

        with open(master_csv, "a", newline="") as f:
            csv.writer(f).writerow([top_design['id'], seq_id, top_design['iptm'], round(top_design['rmsd'], 2), top_design['bad_atoms'], top_design['helicity'], top_design['dG'], top_design['Kd']])
            
        shutil.copy(top_design['path'], os.path.join(final_dir, f"{top_design['id']}_best.cif"))
    else:
        print(f"    No trajectories passed structural validation limits.")

    print("----------------------------------------------------", flush=True)

if not os.path.exists(final_dir) or not os.listdir(final_dir):
    print(" PHASE 2 COMPLETE: No sequences survived the structural gatekeepers.", flush=True)