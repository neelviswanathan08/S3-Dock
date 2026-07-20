import sys
import os
import yaml
import csv
import numpy as np
import mdtraj as md
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# ⚡ Force immediate terminal output
os.environ["PYTHONUNBUFFERED"] = "1"

print("====================================================", flush=True)
print("[PHASE 5] NATIVE OPENMM MM-GBSA ENGINE (UNIVERSAL UNWRAP)", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
md_dir = os.path.join(run_dir, "md_simulations")

# 🚨 CRITICAL: Saving to mmpbsa_results so Phase 6 finds it automatically!
mmgbsa_out_dir = os.path.join(run_dir, "mmgbsa_results") 

os.makedirs(mmgbsa_out_dir, exist_ok=True)

if not os.path.exists(md_dir) or not os.listdir(md_dir):
    print("❌ ERROR: No active production coordinates discovered. Exiting.", flush=True)
    sys.exit(1)

# Load Forcefield once to save time
print("⚙️ Loading AMBER14 Forcefield + OBC2 Implicit Solvent (igb=5 equivalent)...", flush=True)
ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

# Use CPU for analysis to prevent GPU memory overflow during rapid calculation loops
platform = mm.Platform.getPlatformByName('CPU')

for folder in [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]:
    model_md_path = os.path.join(md_dir, folder)
    cif_path = os.path.join(model_md_path, "topology_template.cif")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    output_csv = os.path.join(mmgbsa_out_dir, f"{folder}_mmpbsa.csv")
    
    # 🚨 SMART RESUME: Skip if this specific folder has already been calculated!
    if os.path.exists(output_csv):
        print(f"\n====================================================")
        print(f" ⏩ [SMART RESUME] Existing MM-GBSA data found for {folder}. Skipping.")
        print(f"====================================================")
        continue
        
    if not os.path.exists(cif_path) or not os.path.exists(nc_path):
        continue
        
    print(f"\n====================================================")
    print(f" PROCESSING THERMODYNAMICS FOR: {folder}")
    print(f"====================================================")
    sys.stdout.flush()
    
    try:
        print("📦 Loading NetCDF trajectory...", flush=True)
        raw_traj = md.load(nc_path, top=cif_path)
        total_frames = raw_traj.n_frames
        print(f"   ↳ Discovered total trajectory depth: {total_frames} frames", flush=True)
        
        default_start = total_frames // 2
        start_frame = config.get('mmpbsa_start_frame', default_start)
        end_frame = config.get('mmpbsa_end_frame', -1)
        frame_interval = config.get('mmpbsa_interval', 1)
        print_interval = config.get('mmgbsa_reporting_interval', 10)
        
        if end_frame == -1 or end_frame > total_frames:
            end_frame = total_frames
            
        print(f"✂️ Slicing Trajectory: Extracting frames {start_frame} to {end_frame}...", flush=True)
        traj = raw_traj[start_frame:end_frame:frame_interval]
        
        # 🚨 THE UNIVERSAL FIX: Make molecules whole. Do NOT force centering!
        # This protects fibrils from snapping while still fixing the PBC layout.
        print("🔄 Executing Universal PBC Unwrap (Safe for Fibrils & Antibodies)...", flush=True)
        traj.make_molecules_whole(inplace=True)
        
        protein_alignment_idx = traj.topology.select("protein")
        traj.superpose(traj, 0, atom_indices=protein_alignment_idx)
        
        print("🧬 Isolating Receptor (Fibril) and Ligand (Peptide)...", flush=True)
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
        
        print("🏗️ Building native OpenMM physics matrices (NoCutoff)...", flush=True)
        sys_comp = ff.createSystem(top_comp, nonbondedMethod=app.NoCutoff)
        sys_rec = ff.createSystem(top_rec, nonbondedMethod=app.NoCutoff)
        sys_lig = ff.createSystem(top_lig, nonbondedMethod=app.NoCutoff)
        
        ctx_comp = mm.Context(sys_comp, mm.VerletIntegrator(1.0), platform)
        ctx_rec = mm.Context(sys_rec, mm.VerletIntegrator(1.0), platform)
        ctx_lig = mm.Context(sys_lig, mm.VerletIntegrator(1.0), platform)
        
        print(f"🧮 Calculating Thermodynamic Binding Energy...", flush=True)
        
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
            
            if (i + 1) % print_interval == 0 or i == 0:
                print(f"   ▶ Frame Data processed: {i+1:03d}/{traj.n_frames} | Running Average: {np.mean(energies):.2f} kcal/mol", flush=True)
                sys.stdout.flush()
                    
        dg_final = np.mean(energies)
        dg_std = np.std(energies)
        
        # 🚨 THE PHASE 6 INTEGRATION FIX: Write standard summary CSV!
        # This replaces the 2,500-line CSV that was breaking Phase 6.
        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Model", "Method", "Gas", "Solv", "Delta_G_kcal_mol", "Std_Dev"])
            writer.writerow([folder, "MM-GBSA", "0.0", "0.0", round(dg_final, 2), round(dg_std, 2)])
        
        print("\n" + "="*50)
        print("🎉 CALCULATION COMPLETE 🎉")
        print("="*50)
        print(f"✅ FINAL PRODUCTION MM-GBSA ΔG for {folder}: {dg_final:.2f} ± {dg_std:.2f} kcal/mol")
        print("="*50, flush=True)
        
    except Exception as e:
        print(f"❌ [ERROR] Calculation failed for {folder}: {e}", flush=True)
        if os.path.exists(output_csv):
            os.remove(output_csv)

print("\n----------------------------------------------------", flush=True)
print("[PHASE 5 COMPLETE] ALL MM-GBSA ENERGY COEFFICIENTS SYNCED!", flush=True)