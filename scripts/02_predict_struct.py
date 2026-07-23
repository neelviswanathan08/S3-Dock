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
import time
import glob
from scipy.spatial.distance import cdist
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import BiopythonWarning

# Initialize SOTA Rosetta Engine (Runs once silently)
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
pyrosetta.init("-mute all")

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

samples = config.get('samples_per_seed', 1)

if config.get('benchmark_mode', False):
    print(f"[BENCHMARK ROUTINE] Processing 1 prediction per seed across fasta library...", flush=True)

# ---------------------------------------------------------------------------
# STRUCTURE SELECTORS & AUXILIARY HELPERS
# ---------------------------------------------------------------------------
class TargetSelect(Select):
    def accept_chain(self, chain):
        return chain.id != binder_id

# ---------------------------------------------------------------------------
# UNIVERSAL US-ALIGN ENGINE (MONOMER & MULTIMER)
# ---------------------------------------------------------------------------
def calculate_universal_tm_score(pred_cif_path, ref_cif_path):
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
            print("   -> [CRITICAL ERROR] US-align executable not found in PATH!", flush=True)
            return 0.0

        ram_disk_dir = "/dev/shm"
        temp_dir = tempfile.mkdtemp(dir=ram_disk_dir) if os.path.exists(ram_disk_dir) else tempfile.mkdtemp()
        
        parser = MMCIFParser(QUIET=True)
        pred_struct = parser.get_structure("pred", pred_cif_path)
        temp_cif = os.path.join(temp_dir, "temp_target_only.cif")
        
        io = MMCIFIO()
        io.set_structure(pred_struct)
        io.save(temp_cif, TargetSelect())
        
        cmd = [usalign_bin, temp_cif, ref_cif_path, "-mm", "1"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        tm_score = 0.0
        lines = result.stdout.split('\n')
        for line in lines:
            if "TM-score=" in line and "normalized by length of Structure_2" in line:
                try:
                    tm_score = float(line.split()[1])
                except:
                    pass
        
        if tm_score == 0.0:
            for line in lines:
                if "TM-score=" in line:
                    try:
                        tm_score = max(tm_score, float(line.split()[1]))
                    except:
                        pass

        shutil.rmtree(temp_dir, ignore_errors=True)
        return tm_score

    except Exception as e:
        print(f"   -> [ERROR] US-Align engine exception: {e}", flush=True)
        return 0.0

# ---------------------------------------------------------------------------
# BLIND PREDICTIVE POCKET FILTER
# ---------------------------------------------------------------------------
def check_interface_contacts(cif_path, binder_id, pocket_contacts, dist_cutoff=5.0, min_contacts=3):
    try:
        parser = MMCIFParser(QUIET=True)
        model = parser.get_structure("pred", cif_path)[0]
        
        if binder_id not in model: return False
        
        binder_atoms = [atom.coord for res in model[binder_id] for atom in res]
        if not binder_atoms: return False
        binder_coords = np.array(binder_atoms)
        
        active_contacts = 0
        for chain_id, res_id in pocket_contacts:
            if chain_id in model and res_id in model[chain_id]:
                target_atoms = [atom.coord for atom in model[chain_id][res_id]]
                if not target_atoms: continue
                target_coords = np.array(target_atoms)
                dists = cdist(binder_coords, target_coords)
                if np.min(dists) <= dist_cutoff:
                    active_contacts += 1
                    
        return active_contacts >= min_contacts
    except Exception as e:
        print(f"   -> [WARNING] Pocket contact filter failed to parse: {e}", flush=True)
        return False

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

def run_rosetta_scoring(cif_path):
    """
    SOTA Rosetta Interface Analyzer.
    Replaces PRODIGY. Calculates thermodynamic penalty of separation.
    """
    try:
        pose = pyrosetta.pose_from_file(cif_path)
        
        # Build interface string (e.g., "AB_C")
        receptor_str = "".join(target_chains)
        interface_str = f"{receptor_str}_{binder_id}"
        
        iam = InterfaceAnalyzerMover(interface_str)
        iam.set_compute_packstat(False) # Optimization for speed
        iam.set_compute_interface_energy(True)
        iam.apply(pose)
        
        dG_separated = iam.get_interface_dG()
        return round(dG_separated, 3)
    except Exception as e:
        print(f"   -> [WARNING] Rosetta scoring failed: {e}", flush=True)
        return 999.0

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

global_candidates = []
jobs = parse_fasta(fasta_bridge_path)

if not os.path.exists(master_csv):
    os.makedirs(final_dir, exist_ok=True)
    with open(master_csv, "w", newline="") as f:
        csv.writer(f).writerow([
            "Model_ID", "Seed_ID", "Binder_Sequence", "Seq_Length", "MW_kDa", "pI", "GRAVY",
            "Boltz_ipTM", "Complex_TM_Score", "3D_Helicity_Score", 
            "Rosetta_dG_separated", "Status"
        ])

enforce_boltz = config.get('enforce_boltz_constraints', False)
require_pocket_filter = config.get('require_pocket_filter', True)
min_pocket_contacts = config.get('min_pocket_contacts', 3)
iptm_cutoff = config.get('iptm_threshold', 0.4)

use_msa = config.get('use_msa_server', True)
cached_msa_path = os.path.abspath(os.path.join(run_dir, "target_cached_msa.a3m"))

for idx, (seq_id, sequence) in enumerate(jobs, 1):
    print(f"[INFO] Processing Target Workspace Component: {seq_id}", flush=True)
    design_dir = os.path.join(run_dir, seq_id)
    os.makedirs(design_dir, exist_ok=True)
    yaml_filename = os.path.join(design_dir, f"{seq_id}.yaml")
    
    # --- CALCULATE SEQUENCE PROPERTIES ---
    try:
        pa = ProteinAnalysis(sequence.replace("X", ""))
        seq_len = len(sequence)
        mw = round(pa.molecular_weight() / 1000.0, 2)
        pi = round(pa.isoelectric_point(), 2)
        gravy = round(pa.gravy(), 3)
    except:
        seq_len, mw, pi, gravy = len(sequence), 0.0, 0.0, 0.0
    
    yaml_lines = ["version: 1", "sequences:"]
    for chain, seq in config['target_chains_and_sequences'].items():
        yaml_lines.append(f"  - protein:\n      id: {chain}\n      sequence: '{seq}'")
        if use_msa and os.path.exists(cached_msa_path):
            yaml_lines.append(f"      msa: '{cached_msa_path}'")
            
    yaml_lines.append(f"  - protein:\n      id: {binder_id}\n      sequence: '{sequence}'")
    yaml_lines.append("\ntemplates:")
    for chain in target_chains: yaml_lines.append(f"  - cif: '{target_cif}'\n    chain_id: '{chain}'\n    template_id: '{chain}'")
    
    if 'pocket_contacts' in config and config['pocket_contacts'] and enforce_boltz:
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
        boltz_cmd = f"{boltz_bin} predict {yaml_filename} --use_potentials --out_dir {design_dir} --recycling_steps 10 --diffusion_samples {samples} --override"
        
        if not use_msa or os.path.exists(cached_msa_path):
            boltz_cmd += " --use_msa_server False"
        else:
            boltz_cmd += " --use_msa_server"
        
        process = subprocess.run(boltz_cmd, shell=True, capture_output=True, text=True)
        
        if process.returncode != 0: 
            print(f"   -> [CRITICAL ERROR] Boltz engine crashed for {seq_id}!", flush=True)
            print(f"   -> [ERROR DETAILS] {process.stderr.strip()}", flush=True)
            continue
        
        for root, _, files in os.walk(design_dir):
            for file in files:
                if file.endswith(".cif") and file != os.path.basename(target_cif) and "_target_only" not in file: 
                    generated_models.append(os.path.join(root, file))
                    
        if not generated_models:
            print(f"   -> [FATAL ERROR] Boltz ran, but NO .cif files were generated for {seq_id}!", flush=True)
            print(f"   -> [BOLTZ LOGS]: {process.stderr.strip()} | {process.stdout.strip()}", flush=True)
            continue

        # HARVEST THE MSA CACHE
        if use_msa and not os.path.exists(cached_msa_path):
            found_a3m = glob.glob(f"{design_dir}/**/msa/*.a3m", recursive=True)
            if found_a3m:
                shutil.copy(found_a3m[0], cached_msa_path)
                print(f"   -> [SMART CACHE] Successfully harvested MSA for target! Saved to {cached_msa_path}", flush=True)
                time.sleep(15)
    else:
        print(f"[INFO] Outputs found for {seq_id}. Skipping Boltz-2 prediction.", flush=True)

    print(f"[INFO] Screening generated models via dynamic validation matrix...", flush=True)
    
    for best_cif_path in generated_models:
        filename = os.path.basename(best_cif_path)
        model_num_match = re.search(r'model_(\d+)', filename)
        model_num = model_num_match.group(1) if model_num_match else "0"
        model_id = f"{seq_id}_M{model_num}"
        
        iptm_score = 0.0
        for root, _, files in os.walk(design_dir):
            for file in files:
                if file.endswith(".json"):
                    if samples > 1 and f"model_{model_num}" not in file.lower():
                        continue
                    try:
                        with open(os.path.join(root, file), 'r') as jf:
                            content = jf.read()
                            match = re.search(r'"iptm"\s*:\s*([0-9.]+)', content, re.IGNORECASE)
                            if match:
                                iptm_score = float(match.group(1))
                            else:
                                match_comp = re.search(r'"complex_iptm"\s*:\s*([0-9.]+)', content, re.IGNORECASE)
                                if match_comp:
                                    iptm_score = float(match_comp.group(1))
                    except:
                        pass
        
        if iptm_score == 0.0:
            print(f"   -> [WARNING] ipTM score not found in JSON outputs for {model_id}. Defaulting to 0.0", flush=True)

        helicity_score = calculate_3d_helicity_score(best_cif_path)
        
        passed_structure_filter = True
        status_msg = ""
        
        # 1. UNIVERSAL AI CONFIDENCE GATE
        if iptm_score > 0.0 and iptm_score < iptm_cutoff:
            passed_structure_filter = False
            status_msg = f"REJECTED: Hallucination (ipTM {iptm_score:.2f} < {iptm_cutoff})"
            print(f"   -> [REJECTED] {model_id} interaction does not meet Boltz confidence thresholds ({status_msg})", flush=True)

        # 2. BLIND PREDICTIVE POCKET FILTER
        if passed_structure_filter and require_pocket_filter and 'pocket_contacts' in config:
            passed_pocket = check_interface_contacts(
                best_cif_path, binder_id, config['pocket_contacts'], 
                config.get('max_distance_threshold', 5.0), min_pocket_contacts
            )
            if not passed_pocket:
                passed_structure_filter = False
                status_msg = "REJECTED: Off-Target Binding (Failed Pocket Filter)"
                print(f"   -> [REJECTED] {model_id} spontaneously bound to the wrong region.", flush=True)

        # 3. UNIVERSAL QUATERNARY STRUCTURE GATE
        score_metric = 0.0
        if passed_structure_filter:
            complex_tm_score = calculate_universal_tm_score(best_cif_path, target_cif)
            score_metric = complex_tm_score
            tm_threshold = config.get('tm_threshold', 0.65)
            
            if complex_tm_score < tm_threshold:
                passed_structure_filter = False
                status_msg = f"REJECTED: Target Fold Distorted (TM-Score: {complex_tm_score:.3f} < {tm_threshold})"
                print(f"   -> [REJECTED] {model_id} distorted the native quaternary structure ({status_msg})", flush=True)
            else:
                status_msg = f"PASSED (TM-Score: {complex_tm_score:.3f})"

        # 4. ROSETTA SCORING (Replaced PRODIGY)
        dG = 999.0
        if passed_structure_filter:
            dG = run_rosetta_scoring(best_cif_path)

        # Export metrics immediately
        with open(master_csv, "a", newline="") as f:
            csv.writer(f).writerow([
                model_id, seq_id, sequence, seq_len, mw, pi, gravy,
                iptm_score, round(score_metric, 3), helicity_score, dG, status_msg
            ])

        if passed_structure_filter and dG != 999.0:
            global_candidates.append({
                'path': best_cif_path, 'id': model_id, 'seq_id': seq_id, 'iptm': iptm_score, 
                'metric': score_metric, 'helicity': helicity_score, 'dG': dG
            })
            print(f"   -> [EVALUATED] {model_id} | Rosetta dG: {dG} kcal/mol | ipTM: {iptm_score:.3f} | TM-Score: {score_metric:.3f}", flush=True)

    print("----------------------------------------------------", flush=True)

# =====================================================================
# GLOBAL CHAMPION SELECTION & DEEP DIVE PROTOCOL
# =====================================================================
if not os.path.exists(final_dir):
    os.makedirs(final_dir, exist_ok=True)

if global_candidates:
    # 1. DEDUPLICATE: Get only the absolute best score for each unique sequence
    unique_seq_map = {}
    for cand in global_candidates:
        sid = cand['seq_id']
        if sid not in unique_seq_map or cand['dG'] < unique_seq_map[sid]['dG']:
            unique_seq_map[sid] = cand
            
    # 2. Sort unique sequences by their best dG (lowest is best)
    ranked_unique_seqs = sorted(unique_seq_map.values(), key=lambda x: float(x['dG']))
    top_k = min(config.get("top_k", 3), len(ranked_unique_seqs))
    top_seeds = ranked_unique_seqs[:top_k]
    
    is_benchmark = config.get('benchmark_mode', False)
    
    if is_benchmark:
        print(f"\n[BENCHMARK MODE] Selected Top {top_k} Unique Candidates (Bypassing Deep Dive):", flush=True)
        for i, candidate in enumerate(top_seeds):
            print(f"   #{i+1}: {candidate['id']} | Rosetta dG: {candidate['dG']} kcal/mol", flush=True)
            shutil.copy(candidate['path'], os.path.join(final_dir, f"{candidate['id']}_best.cif"))
        print("\n", flush=True)
        
    else:
        deep_dive_samples = config.get('deep_dive_samples', 10)
        print(f"\n[DEEP DIVE INITIATED] Re-sampling the Top {top_k} sequences with {deep_dive_samples} structures each...", flush=True)
        
        # Create a dictionary to easily grab the raw sequence string for the YAML
        job_dict = dict(jobs)
        boltz_bin = os.path.join(os.path.dirname(sys.executable), "boltz")
        
        for i, candidate in enumerate(top_seeds):
            seq_id = candidate['seq_id']
            sequence = job_dict[seq_id]
            print(f"\n--- [DEEP DIVE] Candidate #{i+1}: {seq_id} (Initial dG: {candidate['dG']}) ---", flush=True)
            
            deep_dir = os.path.join(run_dir, f"{seq_id}_DEEP")
            os.makedirs(deep_dir, exist_ok=True)
            yaml_filename = os.path.join(deep_dir, f"{seq_id}_deep.yaml")
            
            # Rebuild YAML for Deep Dive
            yaml_lines = ["version: 1", "sequences:"]
            for chain, seq in config['target_chains_and_sequences'].items():
                yaml_lines.append(f"  - protein:\n      id: {chain}\n      sequence: '{seq}'")
                if use_msa and os.path.exists(cached_msa_path):
                    yaml_lines.append(f"      msa: '{cached_msa_path}'")
            yaml_lines.append(f"  - protein:\n      id: {binder_id}\n      sequence: '{sequence}'")
            yaml_lines.append("\ntemplates:")
            for chain in target_chains: 
                yaml_lines.append(f"  - cif: '{target_cif}'\n    chain_id: '{chain}'\n    template_id: '{chain}'")
            if 'pocket_contacts' in config and config['pocket_contacts'] and enforce_boltz:
                yaml_lines.append("\nconstraints:\n  - pocket:\n" + f"      binder: {binder_id}\n      contacts:")
                for contact in config['pocket_contacts']: yaml_lines.append(f"        - [{contact[0]}, {contact[1]}]")
                yaml_lines.append(f"      max_distance: {config.get('max_distance_threshold', 5.0)}")
            
            with open(yaml_filename, "w") as f: f.write("\n".join(yaml_lines))
            
            # Execute Boltz-2 Deep Dive
            print(f"   -> Booting Boltz-2 for {deep_dive_samples} exhaustive samples...", flush=True)
            boltz_cmd = f"{boltz_bin} predict {yaml_filename} --use_potentials --out_dir {deep_dir} --recycling_steps 10 --diffusion_samples {deep_dive_samples} --override"
            if not use_msa or os.path.exists(cached_msa_path):
                boltz_cmd += " --use_msa_server False"
            else:
                boltz_cmd += " --use_msa_server"
            
            subprocess.run(boltz_cmd, shell=True, capture_output=True, text=True)
            
            # Evaluate the Deep Dive outputs
            deep_candidates = []
            for root, _, files in os.walk(deep_dir):
                for file in files:
                    if file.endswith(".cif") and file != os.path.basename(target_cif) and "_target_only" not in file:
                        cif_path = os.path.join(root, file)
                        passed = True
                        
                        # Apply Geometric Filters
                        if passed and require_pocket_filter and 'pocket_contacts' in config:
                            if not check_interface_contacts(cif_path, binder_id, config['pocket_contacts'], config.get('max_distance_threshold', 5.0), min_pocket_contacts):
                                passed = False
                        
                        if passed:
                            tm = calculate_universal_tm_score(cif_path, target_cif)
                            if tm < config.get('tm_threshold', 0.65):
                                passed = False
                        
                        # Apply Thermodynamic Filter
                        if passed:
                            dG = run_rosetta_scoring(cif_path)
                            if dG != 999.0:
                                deep_candidates.append({'path': cif_path, 'dG': dG})
            
            # Select the absolute best from the Deep Dive
            if deep_candidates:
                deep_candidates.sort(key=lambda x: float(x['dG']))
                best_deep = deep_candidates[0]
                survival_rate = (len(deep_candidates) / deep_dive_samples) * 100
                
                print(f"   -> [METRICS] {len(deep_candidates)}/{deep_dive_samples} models passed rigid geometry constraints ({survival_rate:.1f}% Survival).", flush=True)
                print(f"   -> [CHAMPION SECURED] Deep Dive Best dG: {best_deep['dG']} kcal/mol", flush=True)
                shutil.copy(best_deep['path'], os.path.join(final_dir, f"{seq_id}_DEEP_best.cif"))
            else:
                print(f"   -> [FATAL FLUKE] All {deep_dive_samples} deep dive models failed. The initial score was a geometric anomaly. Candidate discarded.", flush=True)
                
    print("\n[PHASE 2/3 COMPLETE] Ready for Phase 4 Molecular Dynamics.", flush=True)
else:
    print(f"\n[WARNING] No structural variations across any seeds passed validation gates.", flush=True)