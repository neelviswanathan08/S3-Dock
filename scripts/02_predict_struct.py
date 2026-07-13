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

print("====================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 2 - GLOBAL STRUCTURAL SAMPLING MATRIX", flush=True)
print("====================================================", flush=True)

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
rng_seed = config.get('rng_seed', 42)

samples = config.get('samples_per_seed', 1)
is_benchmark = config.get('benchmark_mode', False)

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
        
        # 🚨 UNIVERSAL UPGRADE: Dynamically detect whatever receptor chains exist in the files
        all_pred_chains = [c.id for c in pred_model.get_chains() if c.id != config['binder_chain_id']]
        all_ref_chains = [c.id for c in ref_model.get_chains() if c.id != config['binder_chain_id']]
        
        ref_centroids, ref_atoms_cache = get_chain_centroids(ref_model, all_ref_chains)
        pred_centroids, pred_atoms_cache = get_chain_centroids(pred_model, all_pred_chains)
        
        # If either file is corrupt or unreadable, fail safely
        if not ref_centroids or not pred_centroids: 
            return 999.0, 999
            
        # Extract CA atoms from the first available receptor chain for a raw alignment
        ref_cids = list(ref_centroids.keys())
        pred_cids = list(pred_centroids.keys())
        
        sup_final = Superimposer()
        # Pair up only as many atoms as the shorter chain contains to prevent zip crashes
        min_atoms = min(len(ref_atoms_cache[ref_cids[0]]), len(pred_atoms_cache[pred_cids[0]]))
        
        sup_final.set_atoms(ref_atoms_cache[ref_cids[0]][:min_atoms], pred_atoms_cache[pred_cids[0]][:min_atoms])
        
        bad_atoms = 0
        rot_f, tran_f = sup_final.rotran
        for r_atom, p_atom in zip(ref_atoms_cache[ref_cids[0]][:min_atoms], pred_atoms_cache[pred_cids[0]][:min_atoms]):
            t_coord = np.dot(p_atom.coord, rot_f) + tran_f
            if np.linalg.norm(r_atom.coord - t_coord) > config['dist_threshold']: 
                bad_atoms += 1
                
        return sup_final.rms, bad_atoms
    except Exception as e: 
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

jobs = parse_fasta(fasta_bridge_path)
global_candidates = []

for idx, (seq_id, sequence) in enumerate(jobs, 1):
    print(f"[QUEUE] [{idx}/{len(jobs)}] DISCOVERY QUEUE: {seq_id}", flush=True)
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

    with open(yaml_filename, "w") as f: f.write("\n".join(yaml_lines))
        
    generated_models = []
    # 🚨 SEED-LEVEL SMART RESUME: Check if Boltz already ran for this specific seed
    for root, _, files in os.walk(design_dir):
        for file in files:
            if file.endswith(".cif") and file != os.path.basename(target_cif): 
                generated_models.append(os.path.join(root, file))
                
    if generated_models:
        print(f"   [SMART RESUME] Found existing Boltz outputs for {seq_id}. Skipping Boltz-2 Prediction.", flush=True)
    else:
        current_track_seed = rng_seed + idx
        boltz_bin = os.path.join(os.path.dirname(sys.executable), "boltz")
        boltz_cmd = f"{boltz_bin} predict {yaml_filename} --use_msa_server --use_potentials --out_dir {design_dir} --recycling_steps 10 --diffusion_samples {samples} --seed {current_track_seed} --override"
        
        process = subprocess.run(boltz_cmd, shell=True, capture_output=True, text=True)
        if process.returncode != 0: continue
        
        for root, _, files in os.walk(design_dir):
            for file in files:
                if file.endswith(".cif") and file != os.path.basename(target_cif): 
                    generated_models.append(os.path.join(root, file))
            
    generated_models.sort(key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
    
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
        print(f"DEBUG {model_id} -> RMSD: {global_rmsd:.2f} (Cutoff: {config['rmsd_cutoff']}) | Bad Atoms: {bad_atoms} (Max: {config['max_bad_atoms']})")
        if global_rmsd > config['rmsd_cutoff'] or bad_atoms > config['max_bad_atoms']:
            if is_benchmark:
                print(f"   -> [REJECTED] {model_id} failed structural filters.", flush=True)
            continue
            
        helicity_score = calculate_3d_helicity_score(best_cif_path)
        dG, Kd = run_prodigy_scoring(best_cif_path, model_id)
        
        if dG != 999.0:
            global_candidates.append({
                'path': best_cif_path, 'id': model_id, 'seq_id': seq_id, 'iptm': iptm_score,
                'rmsd': global_rmsd, 'bad_atoms': bad_atoms, 'helicity': helicity_score,
                'dG': dG, 'Kd': Kd
            })
            if is_benchmark:
                print(f"   -> [EVALUATED] {model_id} | dG: {dG} kcal/mol | Kd: {Kd} | RMSD: {global_rmsd:.2f} A", flush=True)

    print("----------------------------------------------------", flush=True)

# =====================================================================
# GLOBAL TOP-TIER CHAMPION SELECTION
# =====================================================================
if os.path.exists(final_dir): shutil.rmtree(final_dir)
os.makedirs(final_dir, exist_ok=True)

if global_candidates:
    global_candidates.sort(key=lambda x: float(x['dG']))
    top_design = global_candidates[0]
    
    print(f"\n[GLOBAL CHAMPION] Absolute Best Target Candidate: {top_design['id']}", flush=True)
    print(f"                  dG: {top_design['dG']} kcal/mol | Kd: {top_design['Kd']}\n", flush=True)

    with open(master_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Model_ID", "Seed_ID", "Boltz_ipTM", "Global_Backbone_RMSD", "Dislocated_Atoms", "3D_Helicity_Score", "PRODIGY_dG_kcal_mol", "PRODIGY_Kd_M"])
        writer.writerow([top_design['id'], top_design['seq_id'], top_design['iptm'], round(top_design['rmsd'], 2), top_design['bad_atoms'], top_design['helicity'], top_design['dG'], top_design['Kd']])
        
    shutil.copy(top_design['path'], os.path.join(final_dir, f"{top_design['id']}_best.cif"))
else:
    print(f"\n[WARNING] No configurations passed project validation parameters.", flush=True)