# 🧬 S3-Dock
### **Sequence-to-Structure-to-Simulation Pipeline for Feature-Based Peptide Design**

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)]()
[![CUDA](https://img.shields.io/badge/CUDA-12.x-success.svg)]()

---

## Overview

**S3-Dock** is an automated computational pipeline that bridges **1D peptide feature optimization**, **3D structure prediction**, and **4D molecular dynamics simulations** into a unified, reproducible workflow.

Rather than designing peptides solely from sequence or solely from structural information, S3-Dock follows a **feature-based-to-structure philosophy**, where candidate peptides are first engineered to satisfy desired physicochemical properties before undergoing increasingly rigorous structural and thermodynamic validation.

The pipeline integrates modern AI-based structure prediction with classical molecular modeling and molecular dynamics to provide a systematic framework for evaluating peptide binders.

---

# ⚠️ Scope & Original Use Case

S3-Dock was **originally developed to design peptide inhibitors targeting the α-synuclein rod polymorph fibril**.

The software was designed around a **feature-based-to-structure pipeline**, where peptide sequences are generated from target biophysical features (e.g., hydrophobicity, charge, and amphipathicity) and subsequently evaluated through structural prediction, docking, molecular dynamics, and free-energy analysis.

Although the pipeline is highly modular and can be adapted to other peptide–protein interaction problems, **its default assumptions, scoring metrics, and validation strategy were developed and tuned for amyloid fibril inhibitor design**.

Accordingly, users intending to apply S3-Dock to other biological systems should carefully validate target-specific parameters, structural filters, and scoring thresholds before interpreting results.

S3-Dock should therefore be viewed as a **generalizable computational framework**, rather than a universally validated peptide discovery platform.

---

# Design Philosophy

S3-Dock follows a sequential validation strategy:

```
Feature Optimization
        ↓
Sequence Generation
        ↓
AI Structure Prediction
        ↓
Structural Filtering
        ↓
Docking Validation
        ↓
Molecular Dynamics
        ↓
Free Energy Analysis
        ↓
Candidate Ranking
```

Rather than relying on a single predictive model, each stage attempts to remove poor candidates before more computationally expensive simulations are performed.

---

# Pipeline Architecture

The pipeline is divided into six major computational stages.

---

# Phase 1 — Sequence Generation Engine

Novel peptide sequences are generated using a **Simulated Annealing Monte Carlo optimizer**.

Each candidate is evaluated against user-defined physicochemical objectives using the `peptides` Python package.

### Optimization targets include

- Net Charge
- Boman Index
- Hydrophobic Moment
- Secondary Structure Preference
  - α-Helix
  - β-Sheet

Known benchmark peptides (such as experimentally characterized control sequences) may optionally bypass optimization and proceed directly to structural prediction.

---

# Phase 2 — Structural Sampling & Screening

Candidate peptides are folded against the target receptor using **Boltz-2**.

Each predicted complex undergoes several independent filtering stages.

## Structural Confidence

Removes low-confidence AI predictions using:

- ipTM
- pTM
- confidence metrics

## Pocket Validation

Verifies that the peptide binds within the intended target interface rather than an arbitrary surface.

## Quaternary Structure Validation

Uses structural alignment to ensure that peptide binding does not significantly distort the receptor.

Typical metrics include:

- TM-score
- RMSD

## Thermodynamic Filtering

Interface energetics are estimated using **PyRosetta InterfaceAnalyzerMover** to approximate binding quality.

Top-performing candidates may optionally enter an exhaustive re-sampling workflow for additional structural exploration.

---

# Phase 3 — Blind Docking (Optional)

This stage provides an independent validation of AI-predicted complexes.

Docking is performed using **HADDOCK3**, followed by scoring using:

- PRODIGY
- Binding affinity estimates
- Dissociation constant estimation

This stage can be disabled when Boltz-2 multimer predictions are considered sufficiently reliable for downstream simulations.

---

# Phase 4 — Molecular Dynamics

Production molecular dynamics simulations are performed using **OpenMM** with CUDA acceleration.

Typical simulation workflow includes:

- Explicit TIP3P solvent
- AMBER14 force field
- Energy minimization
- Equilibration
- Production trajectory

The workflow includes:

- automatic checkpointing
- resume functionality
- GPU memory protection
- fault recovery

allowing interrupted simulations to continue without restarting.

---

# Phase 5 — MM-GBSA Thermodynamics

Following production simulations, explicit solvent is removed and trajectories are re-evaluated using an implicit solvent model.

Energy analysis includes:

- MM-GBSA
- Running ΔG estimation
- Frame averaging
- Trajectory slicing

This stage provides an approximate estimate of peptide binding energetics over the course of the simulation.

---

# Phase 6 — Automated Discovery Compiler

The final stage aggregates results from every previous phase into a unified report.

Outputs include:

- sequence descriptors
- structural metrics
- docking statistics
- molecular dynamics analyses
- thermodynamic estimates

Publication-ready visualizations are automatically generated, including:

- Backbone RMSD
- Hydrogen Bond Profiles
- Residue Contact Maps
- Binding Hotspots
- Interface Interaction Networks

Cleaned topology and trajectory files are exported for direct visualization in ChimeraX.

---

# Hardware Requirements

S3-Dock was developed on Linux workstations equipped with NVIDIA GPUs.

Recommended specifications include:

| Component | Recommendation |
|------------|---------------|
| Operating System | Ubuntu 22.04+ |
| GPU | NVIDIA CUDA-capable GPU |
| CUDA Toolkit | 12.x |
| Python | 3.10+ |
| Conda | Miniconda / Anaconda |

Because several software packages require conflicting dependencies, S3-Dock utilizes **multiple isolated Conda environments**.

---

# Repository Structure

```text
S3-Dock/
│
├── configs/
├── inputs/
├── scripts/
├── workflows/
├── environments/
├── results/
├── docs/
├── utilities/
├── run_pipeline.sh
└── README.md
```

---

# Installation

Clone the repository:

```bash
git clone https://github.com/neelviswanathan08/S3-Dock.git
cd S3-Dock
```

---

## Environment Setup

Independent Conda environments are recommended.

Example environments include:

```text
boltz_env
```

Contains:

- Boltz-2
- PyTorch
- PyRosetta
- BioPython

```text
openmm_env
```

Contains:

- OpenMM
- MDTraj
- PDBFixer

```text
haddock_env
```

Contains:

- HADDOCK3

---

# Configuration

The pipeline is controlled through a single YAML configuration file.

Example:

```yaml
run_folder_name: "Ultra_aSyn_candidate"

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

---

# Running S3-Dock

Single run:

```bash
./run_pipeline.sh config_targetA.yaml
```

Multiple concurrent targets:

```bash
./run_pipeline.sh config_targetA.yaml &

./run_pipeline.sh config_targetB.yaml &
```

Completed stages are automatically detected, allowing interrupted runs to resume without repeating finished computations.

---

# Output Directory

All results are written to

```text
results/<run_folder_name>/
```

The final summary directory contains:

```text
pipeline_summary_report.csv

binding_hotspots.png

interaction_network.png

stability_profile.png

trajectory.dcd

topology.pdb
```

Additional intermediate outputs from each computational phase are preserved to facilitate reproducibility and downstream analysis.

---

# Workflow Summary

```text
Input Configuration
        │
        ▼
Sequence Optimization
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
MM-GBSA Free Energy Analysis
        │
        ▼
Automated Candidate Ranking
        │
        ▼
Publication-Ready Outputs
```

---

# Current Limitations

S3-Dock is an active research pipeline.

Current limitations include:

- The default workflow was optimized for peptide inhibitor design rather than general small-molecule docking.
- Validation has primarily focused on α-synuclein fibril interfaces.
- MM-GBSA provides approximate free-energy estimates and should not be interpreted as rigorous binding free energies.
- Predicted structures remain dependent on the underlying AI structure prediction models.
- Biological activity ultimately requires experimental validation.

---

# Reproducibility

S3-Dock emphasizes reproducible computational experiments by:

- using configuration-driven workflows
- avoiding hard-coded paths
- preserving intermediate outputs
- supporting resumable execution
- maintaining isolated software environments

---

# Citation

If S3-Dock contributes to published research, please cite this repository along with the underlying software packages used throughout the pipeline (e.g., Boltz-2, HADDOCK3, OpenMM, PyRosetta, MDTraj, and associated dependencies).

A formal manuscript describing S3-Dock is currently in preparation.

---

# License

This project is distributed under the MIT License.

---

# Acknowledgements

S3-Dock builds upon numerous open-source scientific software projects, including:

- Boltz-2
- OpenMM
- HADDOCK3
- PyRosetta
- MDTraj
- BioPython
- peptides
- PDBFixer
- ChimeraX

The authors gratefully acknowledge the developers and maintainers of these projects, whose work made this pipeline possible.

---

## Disclaimer

S3-Dock is intended as a **research software package** for computational peptide design and structural evaluation.

Predictions generated by this software are computational hypotheses and **should not be interpreted as experimental validation or evidence of biological efficacy**. All candidate molecules should undergo appropriate experimental characterization before any biological conclusions are drawn.
