import sys
import os
import yaml
import gzip
import shutil
import subprocess
import argparse

# --- UNIVERSAL PATH ENFORCER ---
cuda_paths = [
    "/usr/local/cuda/lib64",
    "/usr/lib/x86_64-linux-gnu",
    "/usr/lib64",
    "/usr/lib/nvidia"
]
try:
    conda_base = subprocess.check_output(["conda", "info", "--base"]).decode().strip()
    cuda_paths.append(f"{conda_base}/envs/openmm_env/lib")
except Exception:
    pass

existing_ld = os.environ.get("LD_LIBRARY_PATH", "")
os.environ["LD_LIBRARY_PATH"] = ":".join(cuda_paths) + (f":{existing_ld}" if existing_ld else "")

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

print("====================================================", flush=True)
print("[PHASE 4] OPENMM EXTREME-HT PRODUCTION MD RUN", flush=True)
print("====================================================", flush=True)
print(f"[PATH ENFORCER] Injected CUDA Library Paths: {os.environ['LD_LIBRARY_PATH']}", flush=True)

# ---------------------------------------------------------------------------
# 0. CLI & DYNAMIC ROUTER
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--config", type=str, default="config.yaml", 
                    help="Name of the specific config YAML file to use.")
cli_args, unknown = parser.parse_known_args()

# 1. Setup absolute paths using dynamic config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

if os.path.isabs(cli_args.config):
    CONFIG_PATH = cli_args.config
else:
    CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', cli_args.config))

if not os.path.exists(CONFIG_PATH):
    print(f"❌ [FATAL ERROR] Config file not found at: {CONFIG_PATH}", flush=True)
    sys.exit(1)

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
haddock_dir = os.path.join(run_dir, "haddock_runs")
final_dir = os.path.join(run_dir, "top_designs")
md_dir = os.path.join(run_dir, "md_simulations")
os.makedirs(md_dir, exist_ok=True)

run_haddock = config.get('run_haddock', False)
extend_sim = config.get('extend_simulation', False)
input_structures = {}

# -------------------------------------------------------------------
# DYNAMIC ROUTING LOGIC (UPDATED FOR CRYSTAL OVERRIDE)
# -------------------------------------------------------------------
is_pdb_override = config.get('benchmark_pdb_exist', False)
is_cif_override = config.get('benchmark_cif_exist', False)

if is_pdb_override and config.get('benchmark_pdb_path'):
    print(f"[INFO] CRYSTAL OVERRIDE ACTIVE. Sourcing empirical PDB...", flush=True)
    raw_path = config.get('benchmark_pdb_path')
    if not os.path.isabs(raw_path):
        raw_path = os.path.abspath(os.path.join(SCRIPT_DIR, '..', raw_path))
    if not os.path.exists(raw_path):
        print(f"[ERROR] Crystal PDB not found at {raw_path}", flush=True)
        sys.exit(1)
    model_name = os.path.basename(raw_path).replace(".pdb", "")
    input_structures[model_name] = raw_path

elif is_cif_override and config.get('benchmark_cif_path'):
    print(f"[INFO] CRYSTAL OVERRIDE ACTIVE. Sourcing empirical CIF...", flush=True)
    raw_path = config.get('benchmark_cif_path')
    if not os.path.isabs(raw_path):
        raw_path = os.path.abspath(os.path.join(SCRIPT_DIR, '..', raw_path))
    if not os.path.exists(raw_path):
        print(f"[ERROR] Crystal CIF not found at {raw_path}", flush=True)
        sys.exit(1)
    model_name = os.path.basename(raw_path).replace(".cif", "")
    input_structures[model_name] = raw_path

elif run_haddock:
    print("[INFO] HADDOCK mode active. Sourcing initial coordinates from docking clusters...", flush=True)
    if not os.path.exists(haddock_dir):
        print(f"[ERROR] No HADDOCK inputs discovered at {haddock_dir}.", flush=True)
        sys.exit(1)
    for model_name in os.listdir(haddock_dir):
        model_haddock_out = os.path.join(haddock_dir, model_name, "haddock3_output")
        found_files = []
        if os.path.exists(model_haddock_out):
            for root, dirs, files in os.walk(model_haddock_out):
                for file in files:
                    if "cluster_1_model_1" in file and (file.endswith(".pdb") or file.endswith(".pdb.gz")):
                        found_files.append(os.path.join(root, file))
        if found_files:
            input_structures[model_name] = sorted(found_files)[0]
else:
    print("[INFO] HADDOCK bypass active. Sourcing initial coordinates directly from Boltz-2 PDBs/CIFs...", flush=True)
    if not os.path.exists(final_dir):
        print(f"[ERROR] No top_designs folder discovered at {final_dir}.", flush=True)
        sys.exit(1)
    for file in os.listdir(final_dir):
        if file.endswith(".pdb") or file.endswith(".cif"):
            model_name = file.replace("_best.pdb", "").replace(".pdb", "").replace("_best.cif", "").replace(".cif", "")
            input_structures[model_name] = os.path.join(final_dir, file)

if not input_structures:
    print("[ERROR] No valid starting structures found to simulate. Exiting Phase 4.", flush=True)
    sys.exit(1)

total_steps = config.get('md_simulation_steps', 10000)
log_interval = config.get('md_reporting_interval', 1000)
checkpoint_interval = log_interval * 5
target_temp = config.get('md_temperature_kelvin', 310.15)
salt = config.get('md_NaCl_concentration', 0.150)

# -------------------------------------------------------------------
# INDIVIDUAL SIMULATION FUNCTION
# -------------------------------------------------------------------
def run_simulation(model_name, raw_pdb_path):
    model_md_dir = os.path.join(md_dir, model_name)
    os.makedirs(model_md_dir, exist_ok=True)
    
    # Safely handle extensions (PDBFixer needs the right extension)
    ext = ".cif" if raw_pdb_path.endswith(".cif") else ".pdb"
    start_structure_path = os.path.join(model_md_dir, f"start_complex{ext}")
    
    nc_path = os.path.join(model_md_dir, "trajectory.nc")
    part1_nc = os.path.join(model_md_dir, "trajectory_part1.nc")
    part2_nc = os.path.join(model_md_dir, "trajectory_part2.nc")
    
    cif_path = os.path.join(model_md_dir, "topology_template.cif")
    pdb_template_path = os.path.join(model_md_dir, "topology_template.pdb")
    log_path = os.path.join(model_md_dir, "openmm_production.log")
    checkpoint_path = os.path.join(model_md_dir, "production_failsafe.chk")
    
    is_resuming = False
    
    # Check if we should extend the simulation
    if os.path.exists(nc_path) or os.path.exists(part1_nc):
        if extend_sim and os.path.exists(checkpoint_path):
            is_resuming = True
            print(f"   -> [EXTEND MODE] Found checkpoint for {model_name}. Attempting to resume.", flush=True)
        else:
            return f"[SMART RESUME] {model_name} trajectory already exists. Bypassing."

    if not is_resuming:
        if raw_pdb_path.endswith(".gz"):
            with gzip.open(raw_pdb_path, 'rb') as f_in:
                with open(start_structure_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy(raw_pdb_path, start_structure_path)

    try:
        # Build the system exactly as before to ensure topology matches
        fixer = PDBFixer(filename=start_structure_path)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = app.Modeller(fixer.topology, fixer.positions)
        
        # Only add solvent and write templates if it's a fresh run to avoid slight atom re-orderings
        if not is_resuming:
            modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=salt*unit.molar)
            with open(pdb_template_path, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            with open(cif_path, 'w') as f:
                app.PDBxFile.writeFile(modeller.topology, modeller.positions, f)
        else:
            # If resuming, load the exact topology we saved before
            pdb = app.PDBFile(pdb_template_path)
            modeller = app.Modeller(pdb.topology, pdb.positions)

        system = forcefield.createSystem(modeller.topology, 
                                         nonbondedMethod=app.PME, 
                                         nonbondedCutoff=1.0*unit.nanometers, 
                                         constraints=app.HBonds)
        
        integrator = mm.LangevinMiddleIntegrator(target_temp*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
        
        available_platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
        simulation = None
        
        if 'CUDA' in available_platforms:
            try:
                platform = mm.Platform.getPlatformByName('CUDA')
                properties = {'CudaPrecision': 'mixed'}
                simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
            except:
                simulation = None
        
        if simulation is None and 'OpenCL' in available_platforms:
            try:
                platform = mm.Platform.getPlatformByName('OpenCL')
                simulation = app.Simulation(modeller.topology, system, integrator, platform)
            except:
                simulation = None

        if simulation is None:
            simulation = app.Simulation(modeller.topology, system, integrator)

        # -------------------------------------------------------------
        # THE EXTENSION & SAFE FILE LOGIC
        # -------------------------------------------------------------
        if is_resuming:
            simulation.loadCheckpoint(checkpoint_path)
            current_step = simulation.context.getState().getStepCount()
            steps_to_run = total_steps - current_step
            if steps_to_run <= 0:
                return f"   -> [COMPLETE] {model_name} has already reached {total_steps} steps. Nothing to extend."
            print(f"   -> Resuming from step {current_step}. Running {steps_to_run} additional steps...", flush=True)
            
            # Safely move the 50ns file so we don't overwrite it
            if os.path.exists(nc_path):
                shutil.move(nc_path, part1_nc)
            
            # Output the new 50ns to part 2
            active_nc_output = part2_nc
        else:
            simulation.context.setPositions(modeller.positions)
            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(target_temp*unit.kelvin)
            steps_to_run = total_steps
            active_nc_output = nc_path

        # -------------------------------------------------------------
        # THE REPORTERS (NO APPEND ON NC TO PREVENT CRASHES)
        # -------------------------------------------------------------
        try:
            from mdtraj.reporters import NetCDFReporter
            simulation.reporters.append(NetCDFReporter(active_nc_output, log_interval))
        except ImportError:
            import parmed
            simulation.reporters.append(parmed.openmm.NetCDFReporter(active_nc_output, log_interval))

        simulation.reporters.append(app.CheckpointReporter(checkpoint_path, checkpoint_interval))
        
        # It IS safe to append text logs
        mode = 'a' if is_resuming else 'w'
        with open(log_path, mode) as log_file:
            simulation.reporters.append(app.StateDataReporter(log_file, log_interval, step=True, 
                                                              potentialEnergy=True, temperature=True, speed=True,
                                                              append=is_resuming))
            simulation.step(steps_to_run)
            
        # -------------------------------------------------------------
        # THE POST-SIMULATION STITCHER (Combines Part 1 and 2)
        # -------------------------------------------------------------
        if is_resuming:
            print(f"   -> [MERGE] Stitching Part 1 and Part 2 trajectories safely into 100ns master...", flush=True)
            import mdtraj as md
            
            # MDTraj automatically concatenates them seamlessly when passed as a list
            if os.path.exists(part1_nc) and os.path.exists(part2_nc):
                combined_traj = md.load([part1_nc, part2_nc], top=pdb_template_path)
                combined_traj.save_netcdf(nc_path)
                print(f"   -> [MERGE COMPLETE] Master trajectory {nc_path} is fully updated.", flush=True)
            
        return f"[SUCCESS] Extension complete for {model_name}!"

    except Exception as e:
        return f"[ERROR] OpenMM production engine failed for {model_name}: {e}"

# -------------------------------------------------------------------
# THE EXECUTION ENGINE (STRICTLY SEQUENTIAL)
# -------------------------------------------------------------------
print(f"\n[INFO] Initiating STRICTLY SEQUENTIAL MD Executions...", flush=True)
print(f"[NOTE] Live tracking data is being routed to: results/{run_name}/md_simulations/<model>/openmm_production.log", flush=True)

for model_name, raw_pdb_path in input_structures.items():
    print(f"   -> Starting {model_name}...", flush=True)
    result = run_simulation(model_name, raw_pdb_path)
    print(result, flush=True)

print("----------------------------------------------------", flush=True)
print("[PHASE 4 COMPLETE] UNIVERSAL PRODUCTION TRAJECTORIES LOADED SUCCESSFULLY!", flush=True)