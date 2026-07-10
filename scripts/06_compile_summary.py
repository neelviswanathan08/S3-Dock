import sys
import os
import yaml
import shutil
import pandas as pd

print("====================================================", flush=True)
print("[PHASE 6] UNIVERSAL METRICS & ASSET CONSOLIDATION", flush=True)
print("====================================================", flush=True)

# 1. Setup absolute paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml')), 'r') as file: 
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))

top_designs_dir = os.path.join(run_dir, "top_designs")
haddock_dir = os.path.join(run_dir, "haddock_runs")
md_dir = os.path.join(run_dir, "md_simulations")
mmgbsa_dir = os.path.join(run_dir, "mmgbsa_results")
summary_out_dir = os.path.join(run_dir, "final_summary")

os.makedirs(summary_out_dir, exist_ok=True)

# 2. Locate the Phase 2 Master Metrics CSV (Contains your Initial PRODIGY values)
master_csv_path = os.path.join(top_designs_dir, f"caps_2_{run_name}_master_metrics.csv")

if not os.path.exists(master_csv_path):
    print(f"❌ ERROR: Could not find Phase 2 metrics at {master_csv_path}", flush=True)
    sys.exit(1)

# 3. Read the Master CSV via Pandas
df = pd.read_csv(master_csv_path)

# Add our new comprehensive statistical columns for Phase 5 Data
df["MMGBSA_Mean_kcal_mol"] = "N/A"
df["MMGBSA_StdDev_kcal_mol"] = "N/A"
df["MMGBSA_Min_kcal_mol"] = "N/A"
df["MMGBSA_Max_kcal_mol"] = "N/A"

print("\n[COMPILE] Merging Phase 1/2 PRODIGY values with comprehensive MM-GBSA statistics...", flush=True)

# 4. Process each candidate to collect physical files and MM-GBSA stats
for index, row in df.iterrows():
    model_id = str(row["Model_ID"])
    print(f"    -> Formatting data & assets for: {model_id}", flush=True)
    
    dest_boltz_cif = f"{model_id}_01_boltz_design.cif"
    dest_haddock_pdb = f"{model_id}_02_haddock_docked.pdb"
    dest_nc = f"{model_id}_03_md_trajectory.nc"
    dest_cif = f"{model_id}_04_md_topology.cif"

    # --- Find and Copy Boltz CIF ---
    src_boltz = os.path.join(top_designs_dir, f"{model_id}_best.cif")
    if os.path.exists(src_boltz):
        shutil.copy2(src_boltz, os.path.join(summary_out_dir, dest_boltz_cif))

    # --- Find and Copy HADDOCK PDB ---
    haddock_out = os.path.join(haddock_dir, model_id, "haddock3_output")
    if os.path.exists(haddock_out):
        for r, d, f in os.walk(haddock_out):
            for file in f:
                if "cluster_1_model_1" in file and file.endswith(".pdb"):
                    shutil.copy2(os.path.join(r, file), os.path.join(summary_out_dir, dest_haddock_pdb))
                    break

    # --- Find and Copy MD Files ---
    src_nc = os.path.join(md_dir, model_id, "trajectory.nc")
    if os.path.exists(src_nc): 
        shutil.copy2(src_nc, os.path.join(summary_out_dir, dest_nc))
    
    src_topology = os.path.join(md_dir, model_id, "topology_template.cif")
    if os.path.exists(src_topology): 
        shutil.copy2(src_topology, os.path.join(summary_out_dir, dest_cif))

    # --- Calculate Comprehensive MM-GBSA Statistics ---
    mmgbsa_csv = os.path.join(mmgbsa_dir, f"{model_id}_mmgbsa.csv")
    if os.path.exists(mmgbsa_csv):
        try:
            mmgbsa_df = pd.read_csv(mmgbsa_csv)
            energies = mmgbsa_df["Estimated_Binding_Free_Energy_kcal_mol"]
            
            # Map the totality of the energy landscape
            df.at[index, "MMGBSA_Mean_kcal_mol"] = round(energies.mean(), 2)
            df.at[index, "MMGBSA_StdDev_kcal_mol"] = round(energies.std(), 2)
            df.at[index, "MMGBSA_Min_kcal_mol"] = round(energies.min(), 2)
            df.at[index, "MMGBSA_Max_kcal_mol"] = round(energies.max(), 2)
            
        except Exception as e:
            print(f"       [WARNING] Could not parse MM-GBSA for {model_id}: {e}", flush=True)

# 5. Export cleanly formatted Final CSV
final_csv_path = os.path.join(summary_out_dir, "pipeline_summary_report.csv")

# Force Semicolons and a UTF-8 BOM so Excel natively reads the columns
df.to_csv(final_csv_path, index=False, sep=";", encoding="utf-8-sig")

print("\n----------------------------------------------------", flush=True)
print(f"[SUCCESS] Isolated summary folder generated at: {summary_out_dir}", flush=True)
print(f"[SUCCESS] Data beautifully formatted into: pipeline_summary_report.csv", flush=True)
print("[PHASE 6 COMPLETE] WE ARE DONE!", flush=True)