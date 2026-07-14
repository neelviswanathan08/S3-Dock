# 🧬 S3-DOCK: Sequence-Structure-Simulate Docking Pipeline

**S3-Dock** is an automated, high-throughput computational biophysics pipeline designed for the end-to-end discovery, folding, docking, simulation, and validation of designer peptides against target proteins.

The platform orchestrates evolutionary sequence design, deep-learning structural prediction (**Boltz-2**), thermodynamic filtering (**PRODIGY**), blind docking (**HADDOCK3**), molecular dynamics (**OpenMM**), and binding energy calculations (**MM-GBSA**).

---

## ✨ Key Features

* **Dual-Mode Structural Validation:** Switches between rigid Coordinate RMSD/Atom checks (for symmetrical fibrils) and US-Align TM-Score comparison (for globular monomers/hetero-complexes), matched to the geometry of the target.
* **Smart Resume Architecture:** Every phase checks for existing outputs before running. If a prediction, docking run, or MD simulation already exists for a given candidate, S3-Dock skips the heavy computation and moves to post-processing.
* **Automated MDTraj Centering:** Repairs periodic boundary condition (PBC) splitting artifacts in trajectories before analysis, so MM-GBSA energies and visualization output are computed on whole, correctly wrapped molecules.
* **Thermodynamic Gatekeeping:** Integrates PRODIGY to screen structural predictions by binding affinity (ΔG and molar Kd) before candidates advance to the more expensive docking and MD stages.

---

## 🏗️ Pipeline Architecture

S3-Dock executes over six specialized phases. The pipeline is fully modular — run it end-to-end or trigger individual stages.

* **Phase 1: Sequence Generation** (`01_generate_seq.py`)
  Simulated-annealing sequence optimizer that designs peptides matching a target secondary-structure H/P pattern (helix/beta-sheet), optimized toward target net charge, Boman index, and hydrophobic moment.
* **Phase 2: Structural Prediction & Filtering** (`02_predict_struct.py`)
  Runs ML-driven 3D structure prediction via **Boltz-2**. Validates predicted complexes against the reference target using RMSD (fibril mode) or US-Align TM-score (globular mode), then scores surviving candidates via **PRODIGY** and selects the top binder by affinity.
* **Phase 3: Blind Docking** (`03_blind_dock.py`)
  Performs blind rigid-body ensemble docking of the selected sequence against the receptor using **HADDOCK3**.
* **Phase 4: Molecular Dynamics** (`04_md_simulate.py`)
  Runs physical relaxation/production MD via **OpenMM**, followed by PBC wrapping and trajectory centering via **MDTraj**.
* **Phase 5: Thermodynamic Energy Calculation** (`05_calculate_energy.py`)
  Computes single-trajectory MM-GBSA binding energies across the production window of the MD trajectory.
* **Phase 6: Universal Summary Compilation** (`06_compile_summary.py`)
  Aggregates per-candidate results into interaction networks, stability traces, contact hotspots, and a master CSV report.

---

## 📂 Repository Layout

```text
S3-Dock/
├── config.yaml               # Master pipeline runtime configuration
├── run_pipeline.sh           # Wrapper to execute the full pipeline
├── conda_env.sh              # Helper script to build conda environments
├── envs/                     # Environment YAMLs keeping heavy deps isolated
│   ├── boltz_env.yml         # (Boltz-2 / deep learning dependencies)
│   ├── haddock_env.yml       # (HADDOCK3 / Biopython dependencies)
│   └── openmm_env.yml        # (OpenMM / MDTraj / AmberTools dependencies)
├── inputs/                   # Target structures (CIF format)
│   └── target_receptor.cif
├── results/                  # Pipeline outputs (logs, trajectories, scores)
│   └── [run_folder_name]/
└── scripts/                  # Core pipeline modules
    ├── 01_generate_seq.py
    ├── 02_predict_struct.py
    ├── 03_blind_dock.py
    ├── 04_md_simulate.py
    ├── 05_calculate_energy.py
    └── 06_compile_summary.py
```

---

## 💻 Requirements & Installation

### Prerequisites

* Linux-based OS (Ubuntu, CentOS, RHEL, etc.)
* Conda (Miniconda / Anaconda)
* NVIDIA GPU recommended for Boltz-2 (Phase 2) and MD (Phase 4)

### Environment Setup

S3-Dock isolates its environments to avoid dependency conflicts between PyTorch, HADDOCK, and OpenMM.

**Quick setup:**

```bash
chmod +x conda_env.sh
./conda_env.sh
```

**Manual setup:**

```bash
conda env create -f envs/openmm_env.yml
conda env create -f envs/boltz_env.yml
conda env create -f envs/haddock_env.yml
```

---

## ⚙️ Configuration

The pipeline is controlled via `config.yaml`. Key parameters:

* **Workspace & Targets:** `run_folder_name`, and `target_cif` pointing to your receptor file in `inputs/`.
* **Target Mode (`is_fibril`):** `true` for RMSD-based validation (symmetrical fibrils/amyloids), `false` for TM-score validation (globular monomers/hetero-complexes).
* **Peptide Constraints:** `peptide_length`, `secondary_structure`, `boman_index_min`, `target_net_charge`, `hydrophobic_moment_target`.
* **Optimizer Weights:** `charge_weight`, `boman_weight`, `moment_weight` — control how strongly each property is enforced during sequence search.
* **Throughput Control:** `benchmark_mode: true` with `benchmark_samples` to clone a known control sequence for testing; `benchmark_mode: false` with `num_seeds` to run the evolutionary optimizer normally.
* **Molecular Dynamics:** `md_simulation_steps` and related interval settings control OpenMM production depth and MM-GBSA sampling frequency.

---

## 🚀 Usage

Place your target structure (e.g., `target_receptor.cif`) into `inputs/` and configure `config.yaml`.

### 1. Run the full pipeline

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

### 2. Run individual stages (modular mode)

Each phase requires its matching conda environment:

```bash
# Phase 1 & 2 (Boltz environment)
./envs/boltz_env/bin/python scripts/01_generate_seq.py
./envs/boltz_env/bin/python scripts/02_predict_struct.py

# Phase 3 (HADDOCK environment)
./envs/haddock_env/bin/python scripts/03_blind_dock.py

# Phases 4, 5 & 6 (OpenMM environment)
./envs/openmm_env/bin/python scripts/04_md_simulate.py
./envs/openmm_env/bin/python scripts/05_calculate_energy.py
./envs/openmm_env/bin/python scripts/06_compile_summary.py
```

---

## 📊 Outputs & Analytics

All generated files are written to `results/[run_folder_name]/`.

After Phase 6 completes, `results/[run_folder_name]/final_summary/` contains:

* **`pipeline_summary_report.csv`** — merged sequence properties, Boltz confidence metrics, PRODIGY affinity/Kd, MM-GBSA energies, structural validation scores, and per-candidate stats.
* **`[Seed]_stability_profile.png`** — backbone RMSD and inter-chain hydrogen bond stability over the trajectory.
* **`[Seed]_binding_hotspots.png`** — peptide residues acting as binding anchors.
* **`[Seed]_interaction_network.png`** — bipartite interaction network of atomic contacts (<4.5 Å) between peptide and receptor.

---

## 📜 License & Citation

This project is licensed under the MIT License. If you use S3-Dock in your research, please consider citing this repository.
