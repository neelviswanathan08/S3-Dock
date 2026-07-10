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
print("[PHASE 5] NATIVE OPENMM MM-GBSA ENGINE", flush=True)
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

# Load Forcefield once to save time
print("[SYSTEM] Loading AMBER14 Forcefield + OBC2 Implicit Solvent (igb=5)...", flush=True)
ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

# Use CPU for analysis to prevent GPU memory overflow during rapid calculation loops
platform = mm.Platform.getPlatformByName('CPU')

for folder in [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]:
    model_md_path = os.path.join(md_dir, folder)
    cif_path = os.path.join(model_md_path, "topology_template.cif")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    
    if not os.path.exists(cif_path) or not os.path.exists(nc_path):
        continue
        
    print(f"\n====================================================")
    print(f"🧮 PROCESSING THERMODYNAMICS FOR: {folder}")
    print(f"====================================================")
    sys.stdout.flush()
    
    try:
        print("📦 Loading NetCDF trajectory...", flush=True)
        traj = md.load(nc_path, top=cif_path)
        
        print("💧 Stripping explicit water and isolating chains...", flush=True)
        # Select ONLY protein (strips all water, Na+, Cl-, etc.)
        dry_idx = traj.topology.select('protein')
        traj_dry = traj.atom_slice(dry_idx)
        
        # In Boltz, the designed peptide is always the LAST chain.
        # We dynamically find the index of the last chain to isolate it.
        ligand_chain_index = traj_dry.topology.n_chains - 1
        
        rec_idx = traj_dry.topology.select(f'not chainid {ligand_chain_index}')
        lig_idx = traj_dry.topology.select(f'chainid {ligand_chain_index}')
        
        traj_rec = traj_dry.atom_slice(rec_idx)
        traj_lig = traj_dry.atom_slice(lig_idx)
        
        print("⚙️ Converting to OpenMM Native Topologies...", flush=True)
        top_comp = traj_dry.topology.to_openmm()
        top_rec = traj_rec.topology.to_openmm()
        top_lig = traj_lig.topology.to_openmm()
        
        print("🏗️ Building physics systems (NoCutoff)...", flush=True)
        sys_comp = ff.createSystem(top_comp, nonbondedMethod=app.NoCutoff)
        sys_rec = ff.createSystem(top_rec, nonbondedMethod=app.NoCutoff)
        sys_lig = ff.createSystem(top_lig, nonbondedMethod=app.NoCutoff)
        
        ctx_comp = mm.Context(sys_comp, mm.VerletIntegrator(1.0), platform)
        ctx_rec = mm.Context(sys_rec, mm.VerletIntegrator(1.0), platform)
        ctx_lig = mm.Context(sys_lig, mm.VerletIntegrator(1.0), platform)
        
        print(f"▶️ Calculating Binding Free Energy across {traj.n_frames} frames...", flush=True)
        
        output_csv = os.path.join(mmgbsa_out_dir, f"{folder}_mmgbsa.csv")
        energies = []
        
        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Frame", "Complex_E", "Receptor_E", "Ligand_E", "Estimated_Binding_Free_Energy_kcal_mol"])
            
            for i in range(traj.n_frames):
                # Apply coordinates directly from the sliced trajectories
                ctx_comp.setPositions(traj_dry.xyz[i])
                ctx_rec.setPositions(traj_rec.xyz[i])
                ctx_lig.setPositions(traj_lig.xyz[i])
                
                # Extract Potential Energy (kcal/mol)
                e_comp = ctx_comp.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                e_rec = ctx_rec.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                e_lig = ctx_lig.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
                
                # Calculate Binding Free Energy: ΔG = G_complex - (G_receptor + G_ligand)
                dg = e_comp - (e_rec + e_lig)
                energies.append(dg)
                
                writer.writerow([i, round(e_comp, 2), round(e_rec, 2), round(e_lig, 2), round(dg, 2)])
                
                # Live Terminal Feedback exactly as you requested!
                if (i + 1) % 10 == 0 or i == 0:
                    print(f"   ↳ Frame {i+1:03d}/{traj.n_frames} | Running Average: {np.mean(energies):.2f} kcal/mol", flush=True)
                    sys.stdout.flush()
                    
        dg_final = np.mean(energies)
        dg_std = np.std(energies)
        print(f"✅ FINAL MM-GBSA ΔG for {folder}: {dg_final:.2f} ± {dg_std:.2f} kcal/mol", flush=True)
        
    except Exception as e:
        print(f"❌ [ERROR] Calculation failed for {folder}: {e}", flush=True)

print("\n----------------------------------------------------", flush=True)
print("[PHASE 5 COMPLETE] ALL MM-GBSA ENERGY COEFFICIENTS SYNCED!", flush=True)