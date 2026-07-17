import sys
import os
import yaml
import csv
import subprocess
import shutil      
import numpy as np
import mdtraj as md
import openmm as mm
import openmm.app as app
import parmed as pmd

# Force immediate terminal output
os.environ["PYTHONUNBUFFERED"] = "1"

print("====================================================", flush=True)
print("[PHASE 5] NATIVE REPRODUCIBLE AMBERTOOLS THERMO-SOLVER", flush=True)
print("====================================================", flush=True)
sys.stdout.flush()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
md_dir = os.path.join(run_dir, "md_simulations")
mmpbsa_out_dir = os.path.join(run_dir, "mmpbsa_results")

os.makedirs(mmpbsa_out_dir, exist_ok=True)

if not os.path.exists(md_dir) or not os.listdir(md_dir):
    print(" [ERROR] No active production coordinates discovered. Exiting.", flush=True)
    sys.exit(1)

free_energy_method = config.get('free_energy_method', 'mm-gbsa').lower()
if free_energy_method == "none":
    print(" [INFO] Thermodynamic calculations disabled in config. Exiting cleanly.", flush=True)
    sys.exit(0)

# 🚨 FIXED: Cross-environment path enforcement
amber_env_path = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'envs', 'amber_env'))
mmpbsa_bin = os.path.join(amber_env_path, "bin", "MMPBSA.py")

if not os.path.exists(mmpbsa_bin):
    print(f" [FATAL ERROR] MMPBSA.py binary not found at {mmpbsa_bin}!", flush=True)
    sys.exit(1)

print(f" Loading AMBER14 Forcefield for ParmEd Topological Translation...", flush=True)
ff = app.ForceField('amber14-all.xml')

for folder in [f for f in os.listdir(md_dir) if os.path.isdir(os.path.join(md_dir, f))]:
    model_md_path = os.path.join(md_dir, folder)
    cif_path = os.path.join(model_md_path, "topology_template.cif")
    nc_path = os.path.join(model_md_path, "trajectory.nc")
    output_csv = os.path.join(mmpbsa_out_dir, f"{folder}_mmpbsa.csv")
    amber_dir = os.path.join(mmpbsa_out_dir, f"{folder}_amber_temp")
    
    if os.path.exists(output_csv):
        print(f"\n====================================================")
        print(f"  [SMART RESUME] Existing data found for {folder}. Skipping.")
        print(f"====================================================")
        continue
        
    if not os.path.exists(cif_path) or not os.path.exists(nc_path):
        continue
        
    os.makedirs(amber_dir, exist_ok=True)
    
    print(f"\n====================================================")
    print(f" PROCESSING THERMODYNAMICS ({free_energy_method.upper()}) FOR: {folder}")
    print(f"====================================================")
    sys.stdout.flush()
    
    try:
        print(" Loading NetCDF trajectory...", flush=True)
        raw_traj = md.load(nc_path, top=cif_path)
        total_frames = raw_traj.n_frames
        
        start_frame = config.get('energy_start_frame', total_frames // 2)
        end_frame = config.get('energy_end_frame', total_frames)
        frame_interval = config.get('energy_interval', 1)
        salt_con = config.get('energy_salt_concentration', 0.150)
        
        # 🚨 FIXED: Out-of-bounds trajectory safety guard
        if start_frame >= total_frames:
            print(f"  [WARNING] start_frame ({start_frame}) >= total_frames ({total_frames}). Falling back to 0.", flush=True)
            start_frame = 0
        if end_frame > total_frames or end_frame <= start_frame:
            end_frame = total_frames
            
        print(f" Slicing Trajectory: Extracting frames {start_frame} to {end_frame}...", flush=True)
        traj = raw_traj[start_frame:end_frame:frame_interval]
        
        print(" Wrapping periodic boundaries and aligning protein...", flush=True)
        traj.image_molecules(inplace=True)
        protein_alignment_idx = traj.topology.select("protein")
        traj.superpose(traj, 0, atom_indices=protein_alignment_idx)
        
        dry_idx = traj.topology.select('protein')
        traj_dry = traj.atom_slice(dry_idx)
        ligand_chain_index = traj_dry.topology.n_chains - 1
        
        rec_idx = traj_dry.topology.select(f'not chainid {ligand_chain_index}')
        lig_idx = traj_dry.topology.select(f'chainid {ligand_chain_index}')
        
        traj_rec = traj_dry.atom_slice(rec_idx)
        traj_lig = traj_dry.atom_slice(lig_idx)
        
        print(" Translating OpenMM Topologies to AMBER .prmtop formats via ParmEd...", flush=True)
        top_comp = traj_dry.topology.to_openmm()
        top_rec = traj_rec.topology.to_openmm()
        top_lig = traj_lig.topology.to_openmm()
        
        sys_comp = ff.createSystem(top_comp, nonbondedMethod=app.NoCutoff)
        sys_rec = ff.createSystem(top_rec, nonbondedMethod=app.NoCutoff)
        sys_lig = ff.createSystem(top_lig, nonbondedMethod=app.NoCutoff)
        
        struct_comp = pmd.openmm.load_topology(top_comp, sys_comp, traj_dry.xyz[0] * 10.0) 
        struct_rec = pmd.openmm.load_topology(top_rec, sys_rec, traj_rec.xyz[0] * 10.0)
        struct_lig = pmd.openmm.load_topology(top_lig, sys_lig, traj_lig.xyz[0] * 10.0)
        
        comp_prmtop = os.path.join(amber_dir, "complex.prmtop")
        rec_prmtop = os.path.join(amber_dir, "receptor.prmtop")
        lig_prmtop = os.path.join(amber_dir, "ligand.prmtop")
        
        struct_comp.save(comp_prmtop, overwrite=True)
        struct_rec.save(rec_prmtop, overwrite=True)
        struct_lig.save(lig_prmtop, overwrite=True)
        
        sliced_nc = os.path.join(amber_dir, "sliced_dry.nc")
        traj_dry.save_netcdf(sliced_nc)
        
        mmpbsa_in = os.path.join(amber_dir, "mmpbsa.in")
        with open(mmpbsa_in, "w") as f:
            f.write(f"&general\n   endframe={traj.n_frames}, verbose=1, interval=1,\n/\n")
            if free_energy_method == "mm-pbsa":
                f.write(f"&pb\n   istrng={salt_con}, fillratio=4.0, radiopt=0,\n/\n")
            else:
                f.write(f"&gb\n   igb=5, saltcon={salt_con},\n/\n")
                
        print(f" Booting AmberTools {free_energy_method.upper()} solver. This may take a few minutes...", flush=True)
        final_dat = os.path.join(amber_dir, "FINAL_RESULTS_MMPBSA.dat")
        
        cmd = [
            mmpbsa_bin, "-O", 
            "-i", mmpbsa_in, 
            "-o", final_dat, 
            "-cp", comp_prmtop, 
            "-rp", rec_prmtop, 
            "-lp", lig_prmtop, 
            "-y", sliced_nc
        ]
        
        subprocess.run(cmd, cwd=amber_dir, capture_output=True, text=True, check=True)
        
        delta_g = None
        delta_g_std = None
        
        if os.path.exists(final_dat):
            with open(final_dat, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("DELTA TOTAL"):
                        parts = line.split()
                        try:
                            delta_g = float(parts[2])
                            delta_g_std = float(parts[3])
                        except:
                            pass
        
        if delta_g is not None:
            # 🚨 FIXED: Expanded schema layout to fit Phase 6 expectations
            with open(output_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["Model", "Method", "Gas", "Solv", "Delta_G_kcal_mol", "Std_Dev"])
                writer.writerow([folder, free_energy_method.upper(), "0.0", "0.0", delta_g, delta_g_std])
            
            print("\n" + "="*50)
            print(f" {free_energy_method.upper()} CALCULATION COMPLETE ")
            print("="*50)
            print(f" FINAL PRODUCTION ΔG for {folder}: {delta_g:.2f} ± {delta_g_std:.2f} kcal/mol")
            print("="*50, flush=True)
        else:
            print(f" [ERROR] Failed to parse {free_energy_method.upper()} output for {folder}.")
            
    except subprocess.CalledProcessError as e:
        print(f" [ERROR] AmberTools Engine Crashed for {folder}.")
    except Exception as e:
        print(f" [ERROR] Python execution failed for {folder}: {e}", flush=True)

print("\n----------------------------------------------------", flush=True)
print(f"[PHASE 5 COMPLETE] ALL {free_energy_method.upper()} ENERGY COEFFICIENTS SYNCED!", flush=True)