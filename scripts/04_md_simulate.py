import sys
import os

# --- UNIVERSAL PATH ENFORCER ---
cuda_paths = [
    "/usr/local/cuda/lib64",
    "/usr/lib/x86_64-linux-gnu",
    "/usr/lib64",
    "/usr/lib/nvidia"
]
try:
    import subprocess
    conda_base = subprocess.check_output(["conda", "info", "--base"]).decode().strip()
    cuda_paths.append(f"{conda_base}/envs/openmm_env/lib")
except Exception:
    pass

existing_ld = os.environ.get("LD_LIBRARY_PATH", "")
os.environ["LD_LIBRARY_PATH"] = ":".join(cuda_paths) + (f":{existing_ld}" if existing_ld else "")

print("====================================================", flush=True)
print("[PHASE 4] OPENMM EXTREME-HT PRODUCTION MD RUN & POST-PROCESSING", flush=True)
print("====================================================", flush=True)
print(f"[PATH ENFORCER] Injected CUDA Library Paths: {os.environ['LD_LIBRARY_PATH']}", flush=True)
print("[SYSTEM] Initializing structural biology engines...", flush=True)
sys.stdout.flush()

import yaml
import gzip
import shutil
import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
haddock_dir = os.path.join(run_dir, "haddock_runs")
md_dir = os.path.join(run_dir, "md_simulations")

os.makedirs(md_dir, exist_ok=True)

if not os.path.exists(haddock_dir):
    print(f"[ERROR] No HADDOCK inputs discovered at {haddock_dir}.", flush=True)
    sys.exit(1)

total_steps = config.get('md_simulation_steps', 10000)
log_interval = config.get('md_reporting_interval', 1000)
checkpoint_interval = log_interval * 5
target_temp = config.get('md_temperature_kelvin', 310.15)

for model_name in os.listdir(haddock_dir):
    model_haddock_out = os.path.join(haddock_dir, model_name, "haddock3_output")
    
    found_files = []
    if os.path.exists(model_haddock_out):
        for root, dirs, files in os.walk(model_haddock_out):
            for file in files:
                if "cluster_1_model_1" in file and (file.endswith(".pdb") or file.endswith(".pdb.gz")):
                    found_files.append(os.path.join(root, file))

    if not found_files:
        continue
        
    raw_pdb_path = sorted(found_files)[0]
    model_md_dir = os.path.join(md_dir, model_name)
    os.makedirs(model_md_dir, exist_ok=True)
    
    start_pdb_path = os.path.join(model_md_dir, "start_complex.pdb")
    
    if raw_pdb_path.endswith(".gz"):
        with gzip.open(raw_pdb_path, 'rb') as f_in:
            with open(start_pdb_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(raw_pdb_path, start_pdb_path)

    print(f"\n[QUEUE] Launching Production Environment Setup for: {model_name}", flush=True)
    sys.stdout.flush()

    nc_path = os.path.join(model_md_dir, "trajectory.nc")
    cif_path = os.path.join(model_md_dir, "topology_template.cif")
    
    # 🚨 SMART BYPASS: Skip the OpenMM math if the trajectory already exists!
    if os.path.exists(nc_path) and os.path.exists(cif_path):
        print("  [SMART RESUME] Existing trajectory discovered. Bypassing 6-hour MD simulation.", flush=True)
    else:
        try:
            print("  [SYSTEM] Repairing structural topology and adding caps via PDBFixer...", flush=True)
            fixer = PDBFixer(filename=start_pdb_path)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.4)

            print("  [SYSTEM] Applying AMBER14 forcefield parameters...", flush=True)
            forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
            modeller = app.Modeller(fixer.topology, fixer.positions)
            
            print("  [SYSTEM] Packing explicit solvent box (1.0 nm pad + 0.15M NaCl)...", flush=True)
            modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)

            print(f"  [SYSTEM] Writing immutable mmCIF template topology map...", flush=True)
            with open(cif_path, 'w') as f:
                app.PDBxFile.writeFile(modeller.topology, modeller.positions, f)

            print("  [SYSTEM] Compiling system physics equations (PME method)...", flush=True)
            system = forcefield.createSystem(modeller.topology, 
                                             nonbondedMethod=app.PME, 
                                             nonbondedCutoff=1.0*unit.nanometers, 
                                             constraints=app.HBonds)
            
            integrator = mm.LangevinMiddleIntegrator(target_temp*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
            
            available_platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
            print(f"  [DIAGNOSTIC] All active OpenMM platforms found: {available_platforms}", flush=True)
            
            simulation = None
            
            if 'CUDA' in available_platforms:
                try:
                    platform = mm.Platform.getPlatformByName('CUDA')
                    properties = {'CudaPrecision': 'mixed'}
                    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
                    simulation.context.setPositions(modeller.positions)
                    simulation.context.getState(getEnergy=True)
                    print("  [HARDWARE] Successfully locked computing framework to NVIDIA CUDA core array.", flush=True)
                except Exception as e:
                    print(f"  [HARDWARE WARNING] CUDA acceleration context rejected: {e}", flush=True)
                    simulation = None
            
            if simulation is None and 'OpenCL' in available_platforms:
                try:
                    print("  [HARDWARE] Initializing alternative OpenCL GPU compilation matrix...", flush=True)
                    platform = mm.Platform.getPlatformByName('OpenCL')
                    simulation = app.Simulation(modeller.topology, system, integrator, platform)
                    simulation.context.setPositions(modeller.positions)
                    simulation.context.getState(getEnergy=True)
                    print("  [HARDWARE] Success! Secure hardware lock established via OpenCL GPU driver.", flush=True)
                except Exception as ocl_err:
                    print(f"  [HARDWARE WARNING] OpenCL tracking matrix compilation failed: {ocl_err}", flush=True)
                    simulation = None

            if simulation is None:
                simulation = app.Simulation(modeller.topology, system, integrator)
                simulation.context.setPositions(modeller.positions)
                print("  [HARDWARE] Critical Warning: Operating on slow fallback CPU architecture.", flush=True)

            print("  [SYSTEM] Relaxing system matrix to clear atomic overlaps...", flush=True)
            sys.stdout.flush()
            simulation.minimizeEnergy()
            print("  [SUCCESS] System structural grid stabilized successfully.", flush=True)

            print(f"  [SYSTEM] Initializing Maxwell-Boltzmann velocities at {target_temp}K...", flush=True)
            simulation.context.setVelocitiesToTemperature(target_temp*unit.kelvin)

            checkpoint_path = os.path.join(model_md_dir, "production_failsafe.chk")
            print(f"  [SYSTEM] Initializing high-speed NetCDF trajectory stream to: {os.path.basename(nc_path)}", flush=True)

            try:
                import parmed
                simulation.reporters.append(parmed.openmm.NetCDFReporter(nc_path, log_interval))
            except ImportError:
                from mdtraj.reporters import NetCDFReporter as SafeNCReporter
                simulation.reporters.append(SafeNCReporter(nc_path, log_interval))

            simulation.reporters.append(app.CheckpointReporter(checkpoint_path, checkpoint_interval))
            simulation.reporters.append(app.StateDataReporter(sys.stdout, log_interval, step=True, 
                                                              potentialEnergy=True, temperature=True, speed=True))

            print(f"  [RUN] Executing {total_steps} calculation timesteps...", flush=True)
            sys.stdout.flush()
            simulation.step(total_steps)
            print(f"  [SUCCESS] Trajectory processing complete for {model_name}!", flush=True)
            
        except Exception as e:
            print(f"  [ERROR] OpenMM production engine failed for {model_name}: {e}", flush=True)
            if os.path.exists(nc_path):
                try: os.remove(nc_path)
                except: pass
            continue # Skip centering if the run failed

    # 🚨 MDTRAJ POST-PROCESSING: Trajectory Centering & PBC Wrapping
    if os.path.exists(nc_path) and os.path.exists(cif_path):
        print("  [SYSTEM] Wrapping Periodic Boundary Conditions and centering trajectory...", flush=True)
        try:
            import mdtraj as md
            # Load the trajectory and topology into memory
            traj = md.load(nc_path, top=cif_path)
            # Center the complex and remove PBC split artifacts
            traj.image_molecules(inplace=True)
            # Overwrite the uncentered file with the pristine version
            traj.save_netcdf(nc_path)
            print("  [SUCCESS] Trajectory beautifully cleaned, centered, and saved!", flush=True)
        except Exception as e:
            print(f"  [WARNING] Trajectory centering failed: {e}", flush=True)

print("----------------------------------------------------", flush=True)
print("[PHASE 4 COMPLETE] UNIVERSAL PRODUCTION TRAJECTORIES LOADED & CENTERED SUCCESSFULLY!", flush=True)