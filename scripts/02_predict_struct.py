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
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

print("====================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 2 - STRUCTURAL SAMPLING MATRIX (US-ALIGN ENGINE)", flush=True)
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

if config.get('benchmark_mode', False):
    samples = 1000
    print("[BENCHMARK ROUTINE] Matrix forced to 1000 structural diffusion variations.", flush=True)
else:
    samples = config.get('samples_per_seed', 1)

# ---------------------------------------------------------------------------
# REPRODUCIBLE & ULTRA-FAST US-ALIGN WRAPPER (RAM-DISK OPTIMIZED)
# ---------------------------------------------------------------------------
class TargetSelect(Select):
    """Biopython filter to strip the binder out before sending to US-Align"""
    def accept_chain(self, chain):
        return chain.id != binder_id

def calculate_oligomeric_tm_score(pred_cif_path, ref_cif_path):
    try:
        # 1. Fully Reproducible Auto-Detection of the US-Align Binary
        # Scans system path (Conda/Brew) first, then falls back to local compiles
        usalign_bin = None
        for binary_name in ["USalign", "usalign", "US-align", "us-align"]:
            found_path = shutil.which(binary_name)
            if found_path:
                usalign_bin = found_path
                break
                
        if not usalign_bin:
            # Fallback local checks
            for fallback in ["../USalign", "./USalign", "../usalign", "./usalign"]:
                test_path = os.path.abspath(os.path.join(SCRIPT_DIR, fallback))
                if os.path.exists(test_path):
                    usalign_bin = test_path
                    break
                    
        if not usalign_bin:
            print("   -> [CRITICAL ERROR] US-align executable not found in PATH or project directory!", flush=True)
            print("      Please run: 'conda install -c bioconda usalign' or compile locally.", flush=True)
            return 0.0

        # 2. Lightning-Speed Shared Memory (RAM-Disk) Temp Space Selection
        # On Linux, /dev/shm is a virtual memory RAM disk. Writing here bypasses disk storage.
        ram_disk_dir = "/dev/shm"
        if os.path.exists(ram_disk_dir) and os.access(ram_disk_dir, os.W_OK):
            temp_dir = tempfile.mkdtemp(dir=ram_disk_dir)
        else:
            temp_dir = tempfile.mkdtemp() # Standard fallback if RAM disk is unavailable/restricted

        # 3. Strip the binder and write the target-only structure straight to RAM
        parser = MMCIFParser(QUIET=True)
        pred_struct = parser.get_structure("pred", pred_cif_path)
        
        temp_cif = os.path.join(temp_dir, "temp_target_only.cif")
        io = MMCIFIO()
        io.set_structure(pred_struct)
        io.save(temp_cif, TargetSelect())
        
        # 4. Execute the binary via optimized subprocess call
        cmd = [usalign_bin, temp_cif, ref_cif_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        tm_score = 0.0
        for line in result.stdout.split('\n'):
            if "TM-score=" in line:
                try:
                    score = float(line.split()[1])
                    # Ensure we normalize based on the reference structure
                    if "Structure_2" in line:
                        tm_score = score
                        break
                    elif score > tm_score:
                        tm_score = score
                except:
                    pass
                    
        # 5. Clean up RAM space
        shutil.rmtree(temp_dir, ignore_errors=True)
            
        return tm_score
    except Exception as e:
        print(f"   -> [ERROR] US-Align execution failed: {e}", flush=True)
        return 0.0

# ---------------------------------------------------------------------------
# AUXILIARY SCORING ENGINES
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

if not os.path.exists(fasta_bridge_path):
    print(f"[ERROR] Could not find library file at {fasta_bridge_path}.", flush=True)
    exit(1)

jobs = parse_fasta(fasta_bridge_path)
if not os.path.exists(master_csv):
    os.makedirs(final_dir, exist_ok=True)
    with open(master_csv, "w", newline="") as f:
        csv.writer(f).writerow(["Model_ID", "Seed_ID", "Boltz_ipTM", "Oligomeric_TM_Score", "3D_Helicity_Score", "PRODIGY_dG_kcal_mol", "PRODIGY_Kd_M"])

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

    valid_candidates = []
    print(f"[INFO] Screening all {len(generated_models)} generated models via US-Align & PRODIGY Matrix...", flush=True)
    
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

        # 🚨 TRUE RAM-DISK ACCELERATED US-ALIGN SCORE CALCULATION 
        tm_score = calculate_oligomeric_tm_score(best_cif_path, target_cif)
        
        if tm_score < tm_cutoff: 
            print(f"   -> [REJECTED] {model_id} failed TM-Score threshold (Score: {tm_score:.3f} < {tm_cutoff})", flush=True)
            continue
            
        helicity_score = calculate_3d_helicity_score(best_cif_path)
        dG, Kd = run_prodigy_scoring(best_cif_path, model_id)
        
        if dG != 999.0:
            valid_candidates.append({'path': best_cif_path, 'id': model_id, 'iptm': iptm_score, 'tm_score': tm_score, 'helicity': helicity_score, 'dG': dG, 'Kd': Kd})
            print(f"   -> [EVALUATED] {model_id} | dG: {dG} kcal/mol | TM-Score: {tm_score:.3f}", flush=True)

    if valid_candidates:
        valid_candidates.sort(key=lambda x: float(x['dG']))
        top_design = valid_candidates[0]
        print(f"[CHAMPION] Selected Best Conformational State: {top_design['id']} | dG: {top_design['dG']} kcal/mol | TM-Score: {top_design['tm_score']:.3f}", flush=True)
        with open(master_csv, "a", newline="") as f:
            csv.writer(f).writerow([top_design['id'], seq_id, top_design['iptm'], round(top_design['tm_score'], 3), top_design['helicity'], top_design['dG'], top_design['Kd']])
        shutil.copy(top_design['path'], os.path.join(final_dir, f"{top_design['id']}_best.cif"))
    else:
        print(f"[WARNING] No structural variations passed validation gates.", flush=True)
    print("----------------------------------------------------", flush=True)
