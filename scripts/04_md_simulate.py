import sys
import os
import yaml
import gzip
import shutil
import subprocess

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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
haddock_dir = os.path.join(run_dir, "haddock_runs")
final_dir = os.path.join(run_dir, "top_designs")
md_dir = os.path.join(run_dir, "md_simulations")
os.makedirs(md_dir, exist_ok=True)

run_haddock = config.get('run_haddock', False)
input_structures = {}

# -------------------------------------------------------------------
# DYNAMIC ROUTING LOGIC
# -------------------------------------------------------------------
if run_haddock:
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
    print("[INFO] HADDOCK bypass active. Sourcing initial coordinates directly from Boltz-2 PDBs...", flush=True)
    if not os.path.exists(final_dir):
        print(f"[ERROR] No top_designs folder discovered at {final_dir}.", flush=True)
        sys.exit(1)
    for file in os.listdir(final_dir):
        if file.endswith(".pdb"):
            model_name = file.replace("_best.pdb", "").replace(".pdb", "")
            input_structures[model_name] = os.path.join(final_dir, file)

if not input_structures:
    print("[ERROR] No valid starting .pdb structures found to simulate. Exiting Phase 4.", flush=True)
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
    
    start_pdb_path = os.path.join(model_md_dir, "start_complex.pdb")
    nc_path = os.path.join(model_md_dir, "trajectory.nc")
    cif_path = os.path.join(model_md_dir, "topology_template.cif")
    pdb_template_path = os.path.join(model_md_dir, "topology_template.pdb")
    log_path = os.path.join(model_md_dir, "openmm_production.log")
    
    if os.path.exists(nc_path):
        return f"[SMART RESUME] {model_name} trajectory already exists. Bypassing."

    if raw_pdb_path.endswith(".gz"):
        with gzip.open(raw_pdb_path, 'rb') as f_in:
            with open(start_pdb_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(raw_pdb_path, start_pdb_path)

    try:
        fixer = PDBFixer(filename=start_pdb_path)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=salt*unit.molar)

        with open(pdb_template_path, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        with open(cif_path, 'w') as f:
            app.PDBxFile.writeFile(modeller.topology, modeller.positions, f)

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

        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(target_temp*unit.kelvin)

        checkpoint_path = os.path.join(model_md_dir, "production_failsafe.chk")
        
        try:
            import parmed
            simulation.reporters.append(parmed.openmm.NetCDFReporter(nc_path, log_interval))
        except ImportError:
            from mdtraj.reporters import NetCDFReporter as SafeNCReporter
            simulation.reporters.append(SafeNCReporter(nc_path, log_interval))

        simulation.reporters.append(app.CheckpointReporter(checkpoint_path, checkpoint_interval))
        
        # ROUTE LIVE LOGS TO A TEXT FILE TO PREVENT TERMINAL SPAM
        with open(log_path, 'w') as log_file:
            simulation.reporters.append(app.StateDataReporter(log_file, log_interval, step=True, 
                                                              potentialEnergy=True, temperature=True, speed=True))
            simulation.step(total_steps)
            
        return f"[SUCCESS] Trajectory processing complete for {model_name}!"

    except Exception as e:
        if os.path.exists(nc_path):
            try: os.remove(nc_path)
            except: pass
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