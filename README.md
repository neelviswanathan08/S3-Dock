#  S3-DOCK: Sequence-Structure-Simulate Docking Pipeline

S3-Dock is an automated, computational structural pipeline designed for folding, docking, simulating, and validating optimized peptide candidates for binding with target proteins.

The platform connects evolutionary sequence design, deep-learning structural prediction (Boltz-2), blind docking of designed sequence (HADDOCK), solvent molecular dynamics (Open-MM) and thermodynamic energy (MM-GBSA/MM-PBSA) of protein complex into a unified, fault-tolerant workflow.

##  Table of Contents

- [Pipeline Architecture](#️-pipeline-architecture)
- [Repository Layout](#-repository-layout)
- [Requirements & Installation](#-requirements--installation)
- [Configuration](#️-configuration)
- [Usage](#-usage)
- [Outputs & Analytics](#-outputs--analytics)
- [License & Citation](#-license--citation)

##  Pipeline Architecture

S3-Dock executes over six specialized phases. The pipeline is fully modular, allowing you to run end-to-end or trigger individual stages as needed.

**Phase 1: Sequence Generation** (`01_generate_seq.py`)
Utilizes an evolutionary Monte Carlo simulated algorithm. It designs sequences matching specific secondary structure propensities (helices/sheets) under hydrophobic/polar (H/P) patterns, optimized for target net charge, Boman solubility indices, and hydrophobic moments. Desired user features are implemented in the first stage of the pipeline.

**Phase 2: Structural Prediction & Filtering** (`02_predict_struct.py`)
Triggers machine-learning-driven 3D structural folding (via Boltz-2). Allows for user to request binding to particular region of protein target. Implements a local pocket alignment filter to automatically discard misfolded candidates without manual design tracking. Following prediction, top candidate is chosen based on strongest binding affinity with the target protein (via PRODIGY). 

**Phase 3: Blind Docking** (`03_blind_dock.py`)
Prepares target protein and designed sequence's topologies and performs blind rigid-body ensemble docking against the receptor pocket (via HADDOCK). 

**Phase 4: Molecular Dynamics** (`04_md_simulate.py`)
Deploys massive explicit-solvent physical relaxation steps (via OpenMM) under fixed temperature/pressure parameters to observe the true time-series stability of the target-peptide complex.

**Phase 5: Thermodynamic Energy Calculation** (`05_calculate_energy.py`)
Calculates binding and interaction energies (MM-GBSA/MM-PBSA) across the stable production window of the MD trajectory.

**Phase 6: Universal Summary Compilation** (`06_compile_summary.py`)
An automated post-processing layer that aggregates scores and generates publication-ready graphical assets: bipartite interaction networks, dynamic stability traces, global contact hotspots, and the master CSV report.

##  Repository Layout

```
S3-Dock/
├── config.yaml               # Master pipeline runtime configuration profiles
├── run_pipeline.sh           # Universal wrapper to execute the full pipeline
├── conda_env.sh              # Helper script to construct conda environments
├── envs/                     # Environment YAMLs keeping heavy deps isolated
│   ├── boltz_env.yml         # (Boltz/Deep Learning dependencies)
│   ├── haddock_env.yml       # (HADDOCK3 dependencies)
│   └── openmm_env.yml        # (OpenMM/MDTraj dependencies)
├── inputs/                   # Directory for example input structures (CIF format)
│   ├── target_receptor.cif   
│   ├── 2BEG.cif              
│   └── 2BEG_fix.cif          
├── results/                  # Dynamic pipeline outputs (logs, trajectories, scores)
│   └── [run_folder_name]/    
└── scripts/                  # Core pipeline modules
    ├── 01_generate_seq.py    # Phase 1: Sequence / Ligand generation
    ├── 02_predict_struct.py  # Phase 2: 3D structure modeling & filtering
    ├── 03_blind_dock.py      # Phase 3: High-throughput blind docking
    ├── 04_md_simulate.py     # Phase 4: OpenMM molecular dynamics
    ├── 05_calculate_energy.py# Phase 5: Interaction energy computation
    └── 06_compile_summary.py # Phase 6: Aggregate analytics & report generation
```

##  Requirements & Installation

### Prerequisites

- Linux-based OS
- Conda (Miniconda / Anaconda)
- Python 3.8+
- GPU highly recommended for structural prediction (Phase 2) and Molecular Dynamics (Phase 4).

### Environment Setup

S3-Dock splits its environments to keep heavy dependencies (OpenMM, HADDOCK, Boltz) safely isolated and prevent package conflicts.

**Quick Setup:**
Run the provided helper script to automatically build all required environments:

```bash
chmod +x conda_env.sh
bash conda_env.sh
```

**Manual Setup:**
Alternatively, you can create environments individually:

```bash
conda env create -f envs/openmm_env.yml
conda env create -f envs/boltz_env.yml
conda env create -f envs/haddock_env.yml
```

##  Configuration

The entire pipeline is controlled via the `config.yaml` file. Before running, edit this file to suit your experimental needs. Key parameters include:

- **Workspace & Targets:** Set `run_folder_name` and point `target_cif` to your receptor file in the `inputs/` directory.
- **Peptide Constraints:** Define `peptide_length`, `secondary_structure`, target `boman_index_min`, and `net_charge`.
- **Weights of Features:** Change values for `charge_weight`, `boman_weight`, and `moment_weight` to increase or decrease impact of parameter.
- **Universal Structural Pocket Constraints:** Define the `pocket_contacts` using `[Integer, "Chain"]` format to allow the Phase 2 filter to isolate the binding cleft accurately.
- **Molecular Dynamics:** Configure `md_simulation_steps` and thermodynamics intervals to control your OpenMM production depth.

## Usage

Place your target structures (e.g., `target_receptor.cif`) into the `inputs/` directory.

### 1. Run the Full Pipeline

To execute all phases consecutively from end-to-end based on your `config.yaml`:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

### 2. Run Individual Stages (Modular Mode)

S3-Dock is fully modular. You can safely run, re-run, or test any single stage independently:

```bash
python3 scripts/01_generate_seq.py --config config.yaml
python3 scripts/02_predict_struct.py --config config.yaml
python3 scripts/03_blind_dock.py --config config.yaml
python3 scripts/04_md_simulate.py --config config.yaml
python3 scripts/05_calculate_energy.py --config config.yaml
python3 scripts/06_compile_summary.py --config config.yaml
```

> **Note:** Ensure the proper conda environment is activated for the specific script if running manually outside of the `run_pipeline.sh` wrapper.

##  Outputs & Analytics

All generated files are written directly to `results/[run_folder_name]/`.

Upon successful completion of Phase 6, navigate to `results/[run_folder_name]/final_summary/` to access your comprehensive analytics suite:

-  **`pipeline_summary_report.csv`** — The master dataset containing merged sequence logic, Boltz metrics, MM-GBSA energies, global RMSDs, and hydrogen bond counts for every tested seed.
-  **`[Seed]_stability_profile.png`** — Time-series dynamic traces showing Backbone RMSD alongside inter-chain hydrogen bond stability.
-  **`[Seed]_binding_hotspots.png`** — Bar charts identifying exactly which peptide residues serve as binding anchors.
-  **`[Seed]_interaction_network.png`** — A universal bipartite interaction network wiring diagram mapping atomic contacts (<4.5 Å) between the peptide and the target receptor.

##  License & Citation

This project is licensed under the MIT License. If you use S3-Dock in your research, please consider citing this repository.