import sys
import os
import yaml
import csv
import numpy as np
import mdtraj as md
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from concurrent.futures import ProcessPoolExecutor, as_completed

# ⚡ Force immediate terminal output
os.environ["PYTHONUNBUFFERED"] = "1"

print("====================================================", flush=True)
print("[PHASE 5] NATIVE OPENMM MM-GBSA ENGINE (CONCURRENT)", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
md_dir = os.path.join(run_dir, "md_simulations")
mmgbsa_out_dir = os.path.join(run_dir, "mmgbsa_results") 

os.makedirs(mmgbsa_out_dir, exist_ok=True)

if not os.path.exists(md_dir) or not os.listdir(md_dir):
    print("❌ ERROR: No active production coordinates discovered. Exiting.", flush=True)
    sys.exit(1)

run_concurrently = config.get("concurrent_md_execution", False)

# -------------------------------------------------------------------
# INDIVIDUAL MM-GBSA WORKER FUNCTION
# -------------------------------------------------------------------
def process_mmgbsa(folder):
    model_md_path = os.path.join(md_dir, folder)
    cif_path = os.path.join(model_md_path, "topology_template.cif")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    output_csv = os.path.join(mmgbsa_out_dir, f"{folder}_mmpbsa.csv")
    
    if os.path.exists(output_csv):
        return f" ⏩ [SMART RESUME] Existing MM-GBSA data found for {folder}. Skipping."
        
    if not os.path.exists(cif_path) or not os.path.exists(nc_path):
        return f" ⚠️ [WARNING] Missing topology or trajectory for {folder}. Skipping."
    
    try:
        # Load trajectory
        raw_traj = md.load(nc_path, top=cif_path)
        total_frames = raw_traj.n_frames
        
        default_start = total_frames // 2
        start_frame = config.get('mmpbsa_start_frame', default_start)
        end_frame = config.get('mmpbsa_end_frame', -1)
        frame_interval = config.get('mmpbsa_interval', 1)
        
        if end_frame == -1 or end_frame > total_frames:
            end_frame = total_frames
            
        traj = raw_traj[start_frame:end_frame:frame_interval]
        
        # Safe Fibril Unwrap
        traj.make_molecules_whole(inplace=True)
        protein_alignment_idx = traj.topology.select("protein")
        traj.superpose(traj, 0, atom_indices=protein_alignment_idx)
        
        dry_idx = traj.topology.select('protein')
        traj_dry = traj.atom_slice(dry_idx)
        
        ligand_chain_index = traj_dry.topology.n_chains - 1
        rec_idx = traj_dry.topology.select(f'not chainid {ligand_chain_index}')
        lig_idx = traj_dry.topology.select(f'chainid {ligand_chain_index}')
        
        traj_rec = traj_dry.atom_slice(rec_idx)
        traj_lig = traj_dry.atom_slice(lig_idx)
        
        top_comp = traj_dry.topology.to_openmm()
        top_rec = traj_rec.topology.to_openmm()
        top_lig = traj_lig.topology.to_openmm()
        
        # MUST initialize Forcefield and Platform INSIDE the concurrent function
        ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
        platform = mm.Platform.getPlatformByName('CPU')
        
        sys_comp = ff.createSystem(top_comp, nonbondedMethod=app.NoCutoff)
        sys_rec = ff.createSystem(top_rec, nonbondedMethod=app.NoCutoff)
        sys_lig = ff.createSystem(top_lig, nonbondedMethod=app.NoCutoff)
        
        ctx_comp = mm.Context(sys_comp, mm.VerletIntegrator(1.0), platform)
        ctx_rec = mm.Context(sys_rec, mm.VerletIntegrator(1.0), platform)
        ctx_lig = mm.Context(sys_lig, mm.VerletIntegrator(1.0), platform)
        
        energies = []
        
        for i in range(traj.n_frames):
            ctx_comp.setPositions(traj_dry.xyz[i])
            ctx_rec.setPositions(traj_rec.xyz[i])
            ctx_lig.setPositions(traj_lig.xyz[i])
            
            e_comp = ctx_comp.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            e_rec = ctx_rec.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            e_lig = ctx_lig.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            
            dg = e_comp - (e_rec + e_lig)
            energies.append(dg)
                    
        dg_final = np.mean(energies)
        dg_std = np.std(energies)
        
        # Write Phase 6 CSV
        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Model", "Method", "Gas", "Solv", "Delta_G_kcal_mol", "Std_Dev"])
            writer.writerow([folder, "MM-GBSA", "0.0", "0.0", round(dg_final, 2), round(dg_std, 2)])
        
        return f" ✅ FINAL PRODUCTION MM-GBSA ΔG for {folder}: {dg_final:.2f} ± {dg_std:.2f} kcal/mol"
        
    except Exception as e:
        if os.path.exists(output_csv):
            os.remove(output_csv)
        return f" ❌ [ERROR] Calculation failed for {folder}: {e}"

# -------------------------------------------------------------------
# THE EXECUTION ENGINE
# -------------------------------------------------------------------
folders = [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]

if run_concurrently:
    print(f"\n[INFO] INITIATING CONCURRENT MM-GBSA CALCULATIONS...", flush=True)
    # CPU bound: Safely runs up to 3 models simultaneously 
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(process_mmgbsa, folder): folder for folder in folders}
        for future in as_completed(futures):
            print(future.result(), flush=True)
else:
    print(f"\n[INFO] Initiating Sequential MM-GBSA Calculations...", flush=True)
    for folder in folders:
        print(f" -> Processing {folder}...", flush=True)
        print(process_mmgbsa(folder), flush=True)

print("\n----------------------------------------------------", flush=True)
print("[PHASE 5 COMPLETE] ALL MM-GBSA ENERGY COEFFICIENTS SYNCED!", flush=True)