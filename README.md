# 🧬 S3-Dock

> **A Feature-Based-to-Structure Pipeline for Computational Peptide Design**

![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)
![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)
![CUDA](https://img.shields.io/badge/CUDA-12.x-success.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

---

# Overview

**S3-Dock** is a modular computational workflow for **feature-guided peptide design and structural evaluation**.

The pipeline was originally developed to investigate peptide inhibitors targeting the **α-synuclein rod polymorph fibril**, with the goal of integrating sequence optimization, AI-based structural prediction, structural filtering, molecular docking, molecular dynamics simulations, and trajectory-based analysis into a single reproducible workflow.

Unlike many workflows that begin by generating structures directly, S3-Dock follows a **feature-based-to-structure** design philosophy. Candidate peptides are first generated to satisfy user-defined physicochemical objectives before undergoing progressively more computationally intensive structural evaluation.

S3-Dock should be viewed as a research workflow for computational hypothesis generation. It does **not** introduce new prediction models, nor is it intended to replace existing protein design software. Instead, it provides a configurable framework for combining several established computational tools into a single automated pipeline.

---

# Scope & Intended Use

> **Please read this section before using S3-Dock.**

S3-Dock was **originally developed to support computational studies of peptide inhibitors targeting the α-synuclein rod polymorph fibril**.

The default workflow, configuration examples, filtering criteria, and scoring strategy were developed with this application in mind.

Although the pipeline is modular and individual stages may be adapted to other peptide–protein interaction problems, **its default parameters should not be assumed to generalize to unrelated biological systems without additional validation.**

S3-Dock is therefore best viewed as an **alternative workflow** for situations where users wish to begin from predefined biophysical characteristics—such as charge, hydrophobic moment, or amphipathicity—and subsequently evaluate those candidates through progressively more detailed structural analyses.

Users intending to study different targets are encouraged to modify the configuration files, filtering thresholds, and validation strategy to match the biological system under investigation.

---

# Design Philosophy

Many modern peptide design workflows begin with structural generation or end-to-end AI models.

S3-Dock instead adopts a **feature-first** strategy.

Rather than attempting to directly generate binders, the workflow first constructs peptide sequences that satisfy predefined physicochemical objectives before evaluating those sequences through multiple independent structural and thermodynamic analyses.

This staged approach allows inexpensive computational filters to remove poor candidates before more computationally demanding simulations are performed.

The pipeline emphasizes:

- Modular workflow design
- Configuration-driven execution
- Reproducibility
- Checkpoint-based recovery
- Compatibility with established structural biology software

---

# Workflow Overview

```
Input Configuration
        │
        ▼
Feature-Based Sequence Optimization
        │
        ▼
Sequence Library
        │
        ▼
Boltz-2 Structure Prediction
        │
        ▼
Structural Filtering
        │
        ▼
(Optional)
HADDOCK3 Docking
        │
        ▼
OpenMM Molecular Dynamics
        │
        ▼
MM-GBSA Analysis
        │
        ▼
Result Aggregation
```

---

# Pipeline Architecture

S3-Dock is organized into six primary computational stages.

---

# Phase 1 — Feature-Based Sequence Generation

Novel peptide sequences are generated using a **Simulated Annealing Monte Carlo** optimization strategy.

Instead of directly optimizing structural predictions, candidate sequences are evaluated against predefined physicochemical objectives using the `peptides` Python package.

Typical optimization targets include:

- Net Charge
- Boman Index
- Hydrophobic Moment
- Secondary Structure Preference
  - α-Helix
  - β-Sheet

Known benchmark or experimentally characterized sequences may optionally bypass this optimization stage and proceed directly to structural prediction.

---

# Phase 2 — Structure Prediction & Structural Filtering

Candidate peptides are folded against the target receptor using **Boltz-2**.

Each predicted complex then undergoes several independent filtering stages before advancing further into the workflow.

## Structural Confidence

Low-confidence predictions may be removed using metrics such as:

- ipTM
- pTM
- model confidence estimates

---

## Interface Validation

Predicted complexes are evaluated to determine whether peptide binding occurs at the intended interface rather than an unrelated surface.

---

## Structural Integrity

US-align is used to evaluate whether peptide binding substantially perturbs the native receptor structure.

Typical metrics include:

- TM-score
- RMSD

---

## Interface Energetics

PyRosetta's **InterfaceAnalyzerMover** may be used to estimate interface energetics as an additional filtering criterion.

Depending on the workflow configuration, promising candidates may optionally undergo additional structural sampling.

---

# Phase 3 — Optional Docking Validation

This stage provides an independent assessment of predicted complexes.

Docking is performed using **HADDOCK3**.

Docked complexes may subsequently be evaluated using:

- PRODIGY
- Predicted binding affinity
- Estimated dissociation constant

Because Boltz-2 already predicts multimeric complexes, this stage is optional and may be disabled through the configuration file.

---

# Phase 4 — Molecular Dynamics Simulation

Selected peptide–protein complexes are simulated using **OpenMM** with CUDA acceleration.

Typical simulation workflow:

- Explicit TIP3P solvent
- AMBER14 force field
- Energy minimization
- Equilibration
- Production molecular dynamics

Checkpoint files are generated during execution, allowing interrupted simulations to resume without restarting.

---

# Phase 5 — MM-GBSA Analysis

Following production simulations, explicit solvent is removed and trajectories are re-evaluated using an implicit solvent model.

Typical analyses include:

- MM-GBSA energy estimation
- Frame averaging
- Trajectory slicing
- Running ΔG estimation

These values provide approximate energetic descriptors and should be interpreted accordingly.

---

# Phase 6 — Result Aggregation

Outputs generated throughout the workflow are combined into a single summary directory.

Typical outputs include:

- Sequence descriptors
- Structural metrics
- Docking statistics
- Molecular dynamics analyses
- MM-GBSA estimates

Several analysis figures are also generated automatically, including:

- Backbone RMSD
- Hydrogen Bond Profiles
- Residue Contact Maps
- Binding Hotspot Histograms
- Interface Interaction Networks

Trajectory files compatible with ChimeraX are exported for visualization.

---

# Hardware Requirements

S3-Dock was developed and tested on Linux systems equipped with NVIDIA GPUs.

Recommended environment:

| Component | Recommendation |
|------------|---------------|
| Operating System | Ubuntu 22.04+ |
| Python | 3.10+ |
| CUDA Toolkit | 12.x |
| GPU | CUDA-capable NVIDIA GPU |
| Package Manager | Miniconda / Anaconda |

Because multiple third-party packages have incompatible dependencies, S3-Dock uses isolated Conda environments.

---

# Installation

Clone the repository:

```bash
git clone https://github.com/neelviswanathan08/S3-Dock.git

cd S3-Dock
```

---

# Conda Environments

The workflow assumes separate environments for different software packages.

Example environments:

```
envs/

├── boltz_env
├── openmm_env
└── haddock_env
```

Example packages include:

### boltz_env

- Boltz-2
- PyTorch
- PyRosetta
- BioPython

### openmm_env

- OpenMM
- MDTraj
- PDBFixer

### haddock_env

- HADDOCK3

---

# Configuration

Pipeline behavior is controlled through a YAML configuration file.

Example:

```yaml
run_folder_name: Ultra_aSyn_candidate

peptide_length: 22

secondary_structure: helix

target_net_charge: 1.993

hydrophobic_moment_target: 0.678

target_cif: inputs/aSyn_fibril.cif

binder_chain_id: K

md_simulation_steps: 50000000

md_reporting_interval: 10000

mmpbsa_start_frame: 2500

mmpbsa_end_frame: 5000

mmpbsa_interval: 1
```

Because all major settings are configuration driven, multiple independent targets may be processed without modifying the source code.

---

# Running the Pipeline

Single run:

```bash
./run_pipeline.sh config_targetA.yaml
```

Multiple concurrent runs:

```bash
./run_pipeline.sh config_targetA.yaml &

./run_pipeline.sh config_targetB.yaml &
```

Completed workflow stages are automatically detected, allowing interrupted runs to resume without repeating completed computations.

---

# Output Directory

Results are organized by run name.

```
results/

└── <run_folder_name>/
```

Typical outputs include:

```
pipeline_summary_report.csv

binding_hotspots.png

interaction_network.png

stability_profile.png

trajectory.dcd

topology.pdb
```

Intermediate files generated throughout the workflow are retained to facilitate reproducibility and downstream analysis.

---

# Repository Philosophy

S3-Dock intentionally integrates existing computational methods rather than replacing them.

The pipeline emphasizes reproducibility and automation by combining established structural biology software into a configurable workflow.

Users are encouraged to adapt individual stages to their own systems, modify filtering criteria where appropriate, and independently validate any computational predictions.

---

# Current Limitations

S3-Dock is an active research workflow.

Current limitations include:

- The default workflow was developed for peptide inhibitor design against the α-synuclein rod polymorph fibril.
- Performance on unrelated targets has not been comprehensively evaluated.
- MM-GBSA provides approximate energetic estimates rather than rigorous binding free energies.
- Structural predictions remain dependent on the assumptions and limitations of the underlying prediction software.
- Computational predictions should be interpreted as hypotheses requiring experimental validation.

---

# Reproducibility

The workflow emphasizes reproducible computational experiments through:

- Configuration-driven execution
- Isolated software environments
- Checkpoint-based recovery
- Preservation of intermediate outputs
- Modular pipeline stages

---

# Citation

If S3-Dock contributes to published work, please cite this repository together with the original publications describing the external software used throughout the workflow, including Boltz-2, HADDOCK3, OpenMM, PyRosetta, MDTraj, and other dependencies.

---

# License

This project is distributed under the MIT License.

---

# Acknowledgements

S3-Dock builds upon numerous open-source software projects developed by the structural biology and computational chemistry communities, including:

- Boltz-2
- HADDOCK3
- OpenMM
- PyRosetta
- MDTraj
- BioPython
- peptides
- PDBFixer
- ChimeraX

The authors gratefully acknowledge the developers and maintainers of these projects.

---

# Disclaimer

S3-Dock is intended for **research purposes only**.

The workflow generates computational predictions intended to support hypothesis generation and exploratory computational studies. Outputs should not be interpreted as evidence of biological activity or therapeutic efficacy without appropriate experimental validation.
