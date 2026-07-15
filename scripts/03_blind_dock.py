import os
import yaml
import subprocess
import glob
import gzip
from Bio.PDB import MMCIFParser, PDBIO

print("====================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 3 - UNBIASED BLIND GLOBAL DOCKING", flush=True)
print("====================================================", flush=True)

# 1. Setup absolute paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

if not config.get('run_haddock', False):
    print("[SKIP] HADDOCK Global Docking disabled in config.yaml. Skipping Phase 3.", flush=True)
    exit(0)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
final_dir = os.path.join(run_dir, "top_designs")

binder_id = config['binder_chain_id']
target_chains = list(config['target_chains_and_sequences'].keys())

# --- ADVANCED SANITIZATION LAYER: Prevents CNS Overlap Crashes ---
def extract_and_sanitize_target(cif_path, output_pdb_path, chains_to_keep):
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("struct", cif_path)
        
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        from Bio.PDB.Chain import Chain
        
        new_struct = Structure("clean_target")
        new_model = Model(0)
        flattened_chain = Chain("A") # Merge all parts into a clean Chain A
        
        residue_counter = 1
        for model in structure:
            for chain in model:
                if chain.id in chains_to_keep:
                    for residue in chain:
                        if residue.id[0] == " ": # Skip heteroatoms/water
                            res_copy = residue.copy()
                            res_copy.id = (" ", residue_counter, " ")
                            flattened_chain.add(res_copy)
                            residue_counter += 1
                            
        new_model.add(flattened_chain)
        new_struct.add(new_model)
        
        io = PDBIO()
        io.set_structure(new_struct)
        io.save(output_pdb_path)
        return True
    except Exception as e:
        print(f"   [ERROR] Target Sanitization Error: {e}", flush=True)
        return False

def extract_and_sanitize_binder(cif_path, output_pdb_path, binder_chain_id):
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("struct", cif_path)
        
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        from Bio.PDB.Chain import Chain
        
        new_struct = Structure("clean_binder")
        new_model = Model(0)
        flattened_chain = Chain("B") # Keep binder isolated on Chain B
        
        residue_counter = 1
        for model in structure:
            if binder_chain_id in model:
                for residue in model[binder_chain_id]:
                    if residue.id[0] == " ":
                        res_copy = residue.copy()
                        res_copy.id = (" ", residue_counter, " ")
                        flattened_chain.add(res_copy)
                        residue_counter += 1
                        
        new_model.add(flattened_chain)
        new_struct.add(new_model)
        
        io = PDBIO()
        io.set_structure(new_struct)
        io.save(output_pdb_path)
        return True
    except Exception as e:
        print(f"   [ERROR] Binder Sanitization Error: {e}", flush=True)
        return False

if not os.path.exists(final_dir):
    print(f"[ERROR] No 'top_designs' directory found at {final_dir}.", flush=True)
    exit(1)

archived_models = [os.path.join(final_dir, f) for f in os.listdir(final_dir) if f.endswith("_best.cif")]

if not archived_models:
    print("[WARNING] No valid structural candidates found in top_designs folder to dock.", flush=True)
    exit(0)

print(f"[INFO] Found {len(archived_models)} winning candidates to push through global blind dock.", flush=True)

local_haddock_bin = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'envs', 'haddock_env', 'bin', 'haddock3'))
local_prodigy_bin = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'envs', 'boltz_env', 'bin', 'prodigy'))

# GRAB SEEDS FROM CONFIG 
haddock_seed = config.get('rnd_haddock_seed', config.get('rng_seed', 42))

for cif_path in archived_models:
    model_base_name = os.path.basename(cif_path).replace("_best.cif", "")
    print(f"\n[QUEUE] Blind Docking Matrix Setup for: {model_base_name}...", flush=True)
    
    haddock_work_dir = os.path.join(run_dir, "haddock_runs", model_base_name)
    os.makedirs(haddock_work_dir, exist_ok=True)
    
    target_pdb_path = os.path.join(haddock_work_dir, "target_complex.pdb")
    peptide_pdb_path = os.path.join(haddock_work_dir, "candidate_peptide.pdb")
    
    print("   -> Running Advanced Sequential Renumbering on PDB inputs...", flush=True)
    extract_and_sanitize_target(cif_path, target_pdb_path, target_chains)
    extract_and_sanitize_binder(cif_path, peptide_pdb_path, binder_id)
    
    haddock_cfg_path = os.path.join(haddock_work_dir, "blind_dock.toml")
    haddock_output_dir = os.path.join(haddock_work_dir, "haddock3_output")
    
    # 🚨 FIX: MOVED iniseed INTO THE PHYSICS MODULES 🚨
    haddock_toml_content = f"""run_dir = "{haddock_output_dir}"
mode = "local"
ncores = 2

molecules = [
    "{target_pdb_path}",
    "{peptide_pdb_path}"
]

[topoaa]
autohis = true

[rigidbody]
sampling = {config.get('haddock_sampling', 1000)}
cmrest = true
iniseed = {haddock_seed}

[seletop]
select = {config.get('haddock_seletop', 200)}

[flexref]
tolerance = 20
iniseed = {haddock_seed}

[emref]

[clustfcc]

[seletopclusts]
top_clusters = {config.get('haddock_top_clusters', 10)}

[caprieval]
"""
    with open(haddock_cfg_path, "w") as hf:
        hf.write(haddock_toml_content)
        
    print(f"   -> Docking... Running full 7-stage HADDOCK3 pipeline Matrix with seed [{haddock_seed}]...", flush=True)
    subprocess.run(f"{local_haddock_bin} {haddock_cfg_path}", shell=True)
    print(f"   [SUCCESS] HADDOCK3 Complete!", flush=True)

    print(f"   -> Evaluating Top Docked Pose with PRODIGY...", flush=True)
    haddock_models = glob.glob(os.path.join(haddock_output_dir, "**", "*_1.pdb*"), recursive=True)
    if not haddock_models:
        haddock_models = glob.glob(os.path.join(haddock_output_dir, "**", "*.pdb*"), recursive=True)
        
    if haddock_models:
        top_docked_pdb = sorted(haddock_models)[0]
        
        if top_docked_pdb.endswith(".gz"):
            unzipped_pdb = top_docked_pdb.replace(".gz", "")
            with gzip.open(top_docked_pdb, 'rb') as f_in:
                with open(unzipped_pdb, 'wb') as f_out:
                    f_out.write(f_in.read())
            top_docked_pdb = unzipped_pdb
        
        prodigy_cmd = f"{local_prodigy_bin} {top_docked_pdb} --selection A B"
        result = subprocess.run(prodigy_cmd, shell=True, capture_output=True, text=True)
        
        delta_g = 999.0
        kd = "N/A"
        for line in result.stdout.split('\n'):
            if "predicted binding affinity" in line.lower():
                try:
                    delta_g = float(line.split()[-1])
                except:
                    pass
            if "predicted dissociation constant" in line.lower():
                kd = line.split()[-1]
                
        print(f"   [SUCCESS] DOCKED AFFINITY (PRODIGY): {delta_g} kcal/mol | Kd: {kd}", flush=True)
    else:
        print("   [ERROR] HADDOCK3 failed to produce structural outputs. Review the stage logs above.", flush=True)