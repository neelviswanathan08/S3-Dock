import sys
import os
import yaml
import re
import csv
import shutil
import subprocess
import json
import warnings
import tempfile
import numpy as np
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

print("====================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 2 - STRUCTURAL SAMPLING SCREENING MATRIX", flush=True)
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
is_fibril_mode = config.get('is_fibril', True) # True = Rigid RMSD/Atom checks, False = US-Align/TM-Score

if config.get('benchmark_mode', False):
    samples = 1000
    print("[BENCHMARK ROUTINE] Matrix forced to 1000 structural diffusion variations.", flush=True)
else:
    samples = config.get('samples_per_seed', 1)

# ---------------------------------------------------------------------------
# STRUCTURE SELECTORS & AUXILIARY HELPERS
# ---------------------------------------------------------------------------
class TargetSelect(Select):
    """Filter to isolate the target receptor and drop the designer binder"""
    def accept_chain(self, chain):
        return chain.id != binder_id

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

# ---------------------------------------------------------------------------
# METHOD 1: RIGID COORDINATE RMSD & DISLOCATED ATOMS CHECK (For Symmetrical Fibrils)
# ---------------------------------------------------------------------------
def check_smart_holistic_similarity(pred_cif_path, ref_cif_path):
    try:
        from Bio.PDB.Superimposer import Superimposer
        parser = MMCIFParser(QUIET=True)
        ref_model = parser.get_structure("ref", ref_cif_path)[0]
        pred_model = parser.get_structure("pred", pred_cif_path)[0]
        ref_centroids, ref_atoms_cache = get_chain_centroids(ref_model, target_chains)
        pred_centroids, pred_atoms_cache = get_chain_centroids(pred_model, target_chains)
        
        if len(ref_centroids) != len(target_chains) or len(pred_centroids) != len(target_chains): 
            return 999.0, 999
            
        ref_cids, pred_cids = list(ref_centroids.keys()), list(pred_centroids.keys())
        best_rms, best_bad_atoms = 999.0, 999
        
        for trial_p_cid in pred_cids:
            sup_trial = Superimposer()
            sup_trial.set_atoms(ref_atoms_cache[ref_cids[len(ref_cids)//2]], pred_atoms_cache[trial_p_cid]) 
            rot_trial, tran_trial = sup_trial.rotran
            rotated_p_centroids = {p_cid: np.dot(pred_centroids[p_cid], rot_trial) + tran_trial for p_cid in pred_cids}
            cost_matrix = np.zeros((len(target_chains), len(target_chains)))
            for i, r_cid in enumerate(ref_cids):
                for j, p_cid in enumerate(pred_cids):
                    cost_matrix[i, j] = np.linalg.norm(ref_centroids[r_cid] - rotated_p_centroids[p_cid])
            from scipy.optimize import linear_sum_assignment
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            all_ref_atoms, all_pred_atoms = [], []
            for r_idx, p_idx in zip(row_ind, col_ind):
                r_cid, p_cid = ref_cids[r_idx], pred_cids[p_idx]
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
                best_rms, best_bad_atoms = sup_final.rms, bad_atoms
        return best_rms, best_bad_atoms
    except: 
        return 999.0, 999

# ---------------------------------------------------------------------------
# METHOD 2: TAIL-TOLERANT US-ALIGN TM-SCORE EVALUATION (For Globular Monomers)
# ---------------------------------------------------------------------------
def calculate_oligomeric_tm_score(pred_cif_path, ref_cif_path):
    try:
        usalign_bin = None
        for binary_name in ["USalign", "usalign", "US-align", "us-align"]:
            found_path = shutil.which(binary_name)
            if found_path:
                usalign_bin = found_path
                break
        if not usalign_bin:
            for fallback in ["../USalign", "./USalign", "../usalign", "./usalign"]:
                test_path = os.path.abspath(os.path.join(SCRIPT_DIR, fallback))
                if os.path.exists(test_path):
                    usalign_bin = test_path
                    break
        if not usalign_bin:
            print("   -> [CRITICAL ERROR] US-align executable not found in PATH or directory!", flush=True)
            return 0.0

        ram_disk_dir = "/dev/shm"
        if os.path.exists(ram_disk_dir) and os.access(ram_disk_dir, os.W_OK):
            temp_dir = tempfile.mkdtemp(dir=ram_disk_dir)
        else:
            temp_dir = tempfile.mkdtemp()

        parser = MMCIFParser(QUIET=True)
        pred_struct = parser.get_structure("pred", pred_cif_path)
        
        temp_cif = os.path.join(temp_dir, "temp_target_only.cif")
        io = MMCIFIO()
        io.set_structure(pred_struct)
        io.save(temp_cif, TargetSelect())
        
        cmd = [usalign_bin, temp_cif, ref_cif_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        tm_score = 0.0
        lines = result.stdout.split('\n')
        for line in lines:
            if "TM-score=" in line:
                try:
                    score = float(line.split()[1])
                    if "Structure_2" in line:
                        tm_score = score
                        break
                    elif score > tm_score:
                        tm_score = score
                except:
                    pass

        shutil.rmtree(temp_dir, ignore_errors=True)
        return tm_score
    except Exception as e:
        print(f"   -> [ERROR] US-Align engine crash or exception: {e}", flush=True)
        return 0.0

# ---------------------------------------------------------------------------
# AUXILIARY ANALYSIS ENGINES
# ---------------------------------------------------------------------------
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
        delta_g, kd_str = 999.0, "N/A"
        kd_molar = 999.0
        
        for line in result.stdout.split('\n'):
            if "predicted binding affinity" in line.lower():
                try: 
                    delta_g = float(line.split()[-1])
                except: 
                    pass
            if "predicted dissociation constant" in line.lower(): 
                kd_str = line.split()[-1]
                
        # Natively calculate molar Kd from Delta G to avoid parsing scientific notation variations
        if delta_g != 999.0:
            # Kd = exp(dG / (R * T))
            # R = 0.0019872 kcal/(mol*K), T = 298.15 K
            kd_molar = np.exp(delta_g / (0.0019872 * 298.15))
            
        return delta_g, kd_str, kd_molar
    except: 
        return 999.0, "N/A", 999.0

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

if not os.path.exists(fasta_bridge_path):
    print(f"[ERROR] Could not find library file at {fasta_bridge_path}.", flush=True)
    exit(1)

# Initialize global candidate pool to track structural and affinity standouts
global_candidates = []

jobs = parse_fasta(fasta_bridge_path)
if not os.path.exists(master_csv):
    os.makedirs(final_dir, exist_ok=True)
    with open(master_csv, "w", newline="") as f:
        csv.writer(f).writerow([
            "Model_ID", "Seed_ID", "Boltz_ipTM", "Is_Rigid_Fibril_Mode", 
            "Structural_Metric_Score", "3D_Helicity_Score", 
            "PRODIGY_dG_kcal_mol", "PRODIGY_Kd_Text", "PRODIGY_Kd_Molar_Value", "Status"
        ])

for idx, (seq_id, sequence) in enumerate(jobs, 1):
    print(f"[INFO] Processing Target Workspace Component: {seq_id}", flush=True)
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
        
    boltz_bin = os.path.join(os.path.dirname(sys.executable), "boltz")
    
    generated_models = []
    for root, _, files in os.walk(design_dir):
        for file in files:
            if file.endswith(".cif") and file != os.path.basename(target_cif) and "_target_only" not in file: 
                generated_models.append(os.path.join(root, file))
                
    if not generated_models:
        print(f"[INFO] Booting Boltz-2 Engine to generate {samples} structural samples...", flush=True)
        boltz_cmd = f"{boltz_bin} predict {yaml_filename} --use_msa_server --use_potentials --out_dir {design_dir} --recycling_steps 10 --diffusion_samples {samples} --override"
        process = subprocess.run(boltz_cmd, shell=True, capture_output=True, text=True)
        if process.returncode != 0: continue
        
        for root, _, files in os.walk(design_dir):
            for file in files:
                if file.endswith(".cif") and file != os.path.basename(target_cif) and "_target_only" not in file: 
                    generated_models.append(os.path.join(root, file))
    else:
        print(f"[INFO] Outputs found for {seq_id}. Skipping Boltz-2 prediction.", flush=True)

    print(f"[INFO] Screening generated models via dynamic validation matrix...", flush=True)
    
    tm_cutoff = config.get('tm_threshold', 0.75)
    
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

        # SCREEN ALL GENERATED STRUCTURES VIA THERMODYNAMICS FIRST (PRODIGY)
        helicity_score = calculate_3d_helicity_score(best_cif_path)
        dG, Kd_text, Kd_molar = run_prodigy_scoring(best_cif_path, model_id)
        
        # DUAL-MODE VALIDATION FILTER (True = Rigid RMSD/Atom checks, False = US-Align/TM-Score)
        passed_structure_filter = True
        
        if is_fibril_mode:
            # Fibril Mode: Brutal Coordinate Hungarian RMSD & Bad Atom Filtering
            global_rmsd, bad_atoms = check_smart_holistic_similarity(best_cif_path, target_cif)
            score_metric = global_rmsd
            max_bad_atoms = config.get('max_bad_atoms', 350)
            rmsd_cutoff = config.get('rmsd_cutoff', 3.5)
            
            if global_rmsd > rmsd_cutoff or bad_atoms > max_bad_atoms:
                passed_structure_filter = False
                status_msg = f"REJECTED: RMSD {global_rmsd:.2f}A (cutoff {rmsd_cutoff}A), Bad Atoms {bad_atoms} (cutoff {max_bad_atoms})"
                print(f"   -> [REJECTED] {model_id} failed backbone configuration metrics ({status_msg})", flush=True)
            else:
                status_msg = f"PASSED (RMSD: {global_rmsd:.2f}A)"
        else:
            # Globular Mode: Tail-tolerant US-Align TM-Score
            score_metric = calculate_oligomeric_tm_score(best_cif_path, target_cif)
            if score_metric < tm_cutoff:
                passed_structure_filter = False
                status_msg = f"REJECTED: TM-Score {score_metric:.3f} < {tm_cutoff}"
                print(f"   -> [REJECTED] {model_id} failed complex stack validation ({status_msg})", flush=True)
            else:
                status_msg = f"PASSED (TM-Score: {score_metric:.3f})"

        # Export ALL structural metrics and PRODIGY energy affinities immediately to Master Metrics CSV
        with open(master_csv, "a", newline="") as f:
            csv.writer(f).writerow([
                model_id, seq_id, iptm_score, is_fibril_mode, 
                round(score_metric, 3), helicity_score, dG, Kd_text, f"{Kd_molar:.3e}", status_msg
            ])

        # If passed the structural filter gate, add to the pool for final global champion sorting
        if passed_structure_filter and dG != 999.0:
            global_candidates.append({
                'path': best_cif_path, 'id': model_id, 'seq_id': seq_id, 'iptm': iptm_score, 
                'metric': score_metric, 'helicity': helicity_score, 'dG': dG, 'Kd_text': Kd_text, 'Kd_molar': Kd_molar
            })
            metric_label = "RMSD" if is_fibril_mode else "TM-Score"
            print(f"   -> [EVALUATED] {model_id} | dG: {dG} kcal/mol | Kd: {Kd_text} ({Kd_molar:.3e} M) | {metric_label}: {score_metric:.3f}", flush=True)

    print("----------------------------------------------------", flush=True)

# =====================================================================
# GLOBAL CHAMPION SELECTION (EXECUTES ONLY ONCE AT THE END OF PIPELINE)
# =====================================================================
if os.path.exists(final_dir): 
    shutil.rmtree(final_dir)
os.makedirs(final_dir, exist_ok=True)

if global_candidates:
    # Sort ALL structurally valid candidates across ALL seeds strictly by dG (strongest affinity)
    global_candidates.sort(key=lambda x: float(x['dG']))
    top_design = global_candidates[0]
    metric_label = "RMSD" if is_fibril_mode else "TM-Score"
    
    print(f"\n[GLOBAL CHAMPION] Selected Absolute Best Candidate from total library run: {top_design['id']}", flush=True)
    print(f"                  dG: {top_design['dG']} kcal/mol | Kd: {top_design['Kd_text']} ({top_design['Kd_molar']:.3e} M) | {metric_label}: {top_design['metric']:.3f}\n", flush=True)
    
    # Push the absolute champion forward to the final directory for HADDOCK refinement
    shutil.copy(top_design['path'], os.path.join(final_dir, f"{top_design['id']}_best.cif"))
else:
    print(f"\n[WARNING] No structural variations across any seeds passed validation gates.", flush=True)
