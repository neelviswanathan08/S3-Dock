# S3-Dock: Reproducible Workflow for Peptide Discovery and Structural Evaluation

S3-Dock is a computational workflow for the in silico generation, structural screening, docking, and molecular-dynamics refinement of peptide candidates against a target protein or receptor assembly. The pipeline integrates sequence generation, structure prediction, pocket-aware filtering, blind docking, molecular dynamics simulation, and energetic analysis within a reproducible computational framework designed for structured peptide-discovery studies.

This repository supports end-to-end computational analysis in a controlled and traceable manner. Candidate peptide libraries are generated, filtered using biophysical and structural criteria, evaluated through docking and simulation, and compiled into summary outputs suitable for downstream comparative analysis and manuscript reporting.

---

## Overview

The workflow is organized as a modular, end-to-end pipeline comprising the following stages:

1. Generation of peptide candidate libraries under defined physicochemical constraints;
2. Structural screening and validation of candidate designs;
3. Pocket-aware filtering based on contact geometry and interface plausibility;
4. Blind docking of selected candidates to the target receptor;
5. Molecular-dynamics refinement and trajectory analysis;
6. Energetic evaluation of simulated complexes;
7. Compilation of final results into summary tables and figures for comparative analysis.

The repository is intended for reproducible, research-focused evaluation of peptide–protein interaction hypotheses under the software stack and operating conditions specified in the environment files. It is not intended as a universal computational platform for all biomolecular systems without configuration, validation, and domain-specific adaptation.

---

## Scientific workflow

### Phase 1: Sequence Generation
File: `scripts/01_generate_seq.py`

This phase generates a peptide library under user-defined constraints, including peptide length, secondary-structure pattern, net charge, Boman index, and hydrophobic moment. The optimizer employs a stochastic search strategy to evolve sequences toward target physicochemical properties while maintaining a sequence pattern consistent with the desired fold.

### Phase 2: Structural Prediction and Screening
File: `scripts/02_predict_struct.py`

This phase evaluates each candidate sequence for structural plausibility using structure prediction and validation routines. It assesses peptide–target geometry, performs pocket-aware contact checks, and filters out designs that fail interface or structural quality criteria.

### Phase 3: Blind Docking
File: `scripts/03_blind_dock.py`

Surviving designs are evaluated by blind global docking. This stage explores plausible binding poses between the peptide and target without a predefined binding orientation, generating candidate docked complexes for downstream refinement.

### Phase 4: Molecular Dynamics Simulation
File: `scripts/04_md_simulate.py`

Candidate complexes are refined using explicit-solvent molecular dynamics in a reproducible simulation setup. The workflow includes geometry preparation, trajectory generation, and quality-control checks for dynamic stability.

### Phase 5: Energy Calculation
File: `scripts/05_calculate_energy.py`

Binding and interaction energetics are evaluated on the simulated trajectory, yielding quantitative measures relevant to stability and interaction strength.

### Phase 6: Summary Compilation
File: `scripts/06_compile_summary.py`

The final stage compiles per-candidate data into summary plots and tabular outputs. This facilitates comparative analysis across designs and provides a structured basis for manuscript-ready reporting.

---

## Repository structure

```text
S3-Dock/
├── README.md
├── config.yaml
├── conda_env.sh
├── run_pipeline.sh
├── inputs/
│   ├── MDM2.cif
│   ├── PSD-95.cif
│   └── ...
├── results/
│   └── <run_name>/
│       ├── library.fasta
│       ├── top_designs/
│       ├── md_simulations/
│       ├── mmgbsa_results/
│       └── final_summary/
├── scripts/
│   ├── 01_generate_seq.py
│   ├── 02_predict_struct.py
│   ├── 03_blind_dock.py
│   ├── 04_md_simulate.py
│   ├── 05_calculate_energy.py
│   ├── 06_compile_summary.py
│   ├── MD_Simulation_Sanity.py
│   ├── Sequence_Property_check.py
│   ├── Spearman_Correlation.py
│   └── ...
├── envs/
│   ├── blind_dock_env.yml
│   ├── boltz_env.yml
│   ├── openmm_env.yml
│   └── ...
├── Control_yaml/
│   ├── CRIPT_Boltz/
│   ├── MDM2_Comparison/
│   └── MDM2_p53_Boltz/
├── test_benchmarks/
│   ├── CRIPT_negative/
│   ├── CRIPT_positive/
│   ├── MDM2_Bindcraft/
│   ├── MDM2_Colab/
│   ├── MDM2_Negative/
│   └── MDM2_Positive/
├── Compilation_files/
│   ├── CRIPT_VS_CONTROL_BOLTZ.py
│   ├── CRIPT_VS_CONTROL_CRYSTAL.py
│   ├── MDM2_PIPELINE_COMP.py
│   ├── MDM2_VS_CONTROL_BOLTZ.py
│   └── MDM2_VS_CONTROL_CRYSTAL.py
└── ...
```

---

## Requirements

### Hardware

- Linux-based operating system;
- CUDA-capable NVIDIA GPU strongly recommended for structure-prediction and molecular-dynamics stages;
- adequate storage for intermediate trajectory files, docking outputs, and analysis results.

### Software

- Conda or Miniconda;
- Python 3.x compatible with the project environment files;
- GPU-enabled ML stack for the Boltz/structure-prediction stage;
- MD and docking dependencies defined in the repository environment specifications.

---

## Environment Setup

The repository requires dedicated conda environments for the primary computational stages. The corresponding environment files are provided in the `envs/` directory.

### Recommended Setup

```bash
chmod +x conda_env.sh
./conda_env.sh
```

This script creates the local environment prefixes required by the workflow, including the Boltz, blind docking, and OpenMM stages.

### Manual Environment Creation

If preferred, the environments may be created explicitly:

```bash
conda env create --prefix envs/boltz_env -f envs/boltz_env.yml
conda env create --prefix envs/blind_dock_env -f envs/blind_dock_env.yml
conda env create --prefix envs/openmm_env -f envs/openmm_env.yml
```

> The environment names and prefixes should match those used by the pipeline scripts to maintain reproducibility.

---

## Configuration

The central configuration file is `config.yaml`. It defines the run name, peptide design constraints, target receptor, docking settings, molecular dynamics parameters, and summary behavior.

Key configuration categories include:

- `run_folder_name`: output directory under `results/`;
- `num_seeds`: number of candidate sequences to generate;
- `peptide_length`: target peptide length;
- `secondary_structure`: helix or beta-sheet design mode;
- `peptide_property`: amphipathic, hydrophobic, or hydrophilic;
- `target_cif`: receptor structure used for docking and analysis;
- `target_chains_and_sequences`: receptor sequence metadata;
- `benchmark_mode`: optional benchmark or control-sequence mode;
- `md_simulation_steps`: molecular-dynamics production length;
- `analysis_start_frame`, `analysis_end_frame`, `analysis_stride`: analysis window configuration;
- `rng_seed`: reproducibility seed for stochastic generation.

A typical run is configured by modifying the target receptor path, peptide constraints, and output folder name before executing the pipeline.

---

## Full Reproducible Execution

The recommended approach is to run the complete workflow through the project runner:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh config.yaml
```

The runner executes the full computational sequence and automatically skips tasks already present in the output directory.

### End-to-End Workflow

The pipeline proceeds in the following order:

1. Sequence library generation;
2. Structure prediction and screening;
3. Blind docking;
4. Molecular dynamics simulation;
5. Energetic analysis;
6. Final summary generation.

---

## Modular Execution

For debugging or targeted validation, individual stages may be executed separately:

```bash
./envs/boltz_env/bin/python scripts/01_generate_seq.py --config config.yaml
./envs/boltz_env/bin/python scripts/02_predict_struct.py --config config.yaml
./envs/blind_dock_env/bin/python scripts/03_blind_dock.py --config config.yaml
./envs/openmm_env/bin/python scripts/04_md_simulate.py --config config.yaml
./envs/openmm_env/bin/python scripts/05_calculate_energy.py --config config.yaml
./envs/openmm_env/bin/python scripts/06_compile_summary.py --config config.yaml
```

This modular mode is useful for evaluating a single stage, comparing parameter sets, or diagnosing failures in a reproducible workflow.

---

## Reproducibility Practices

For manuscript-grade reproducibility, the following information should be recorded for each run:

- Exact `config.yaml` used;
- `run_folder_name`;
- RNG seed value in the configuration;
- Target input CIF file;
- Environment creation commands and package versions;
- Date and time of execution;
- Final result directory under `results/`.

The workflow is designed to preserve reproducibility by using explicit seeds and by storing results within a run-specific directory.

---

## Software Stack and Version Constraints

The repository is configured around a specific, versioned computational stack defined in the environment files. The workflow is intended to be reproducible only under the exact package combinations described below; results obtained with alternative versions or mixed environments should be considered non-equivalent unless separately validated.

### Boltz / Structure Prediction Environment

Environment file: `envs/boltz_env.yml`

- Python: `3.12.12`
- PyTorch: `2.9.1`
- Boltz: `2.2.1`
- Biopython: `1.84`
- RDKit: `2025.9.1`
- scikit-learn: `1.6.1`
- SciPy: `1.13.1`
- NumPy: `<2.0`
- `peptides` package
- `PyRosetta` (`pyrosetta`)

This environment is the primary engine for sequence optimization and structural screening. It is intended for Linux systems with NVIDIA CUDA support and is most reliable when GPU acceleration is available for the Boltz model and the underlying ML stack.

### Blind Docking Environment

Environment file: `envs/blind_dock_env.yml`

- Python: `3.10.20`
- HADDOCK3: `2026.5.0`
- Biopython: `1.87`
- NumPy: `2.2.6`
- SciPy: `1.15.3`
- pandas: `2.3.3`
- PRODIGY-prot: `2.3.0`
- PRODIGY-lig: `1.1.4`
- freesasa: `2.2.1`

This environment implements the docking and thermodynamic screening layer used by the repository.

### OpenMM / Molecular Dynamics Environment

Environment file: `envs/openmm_env.yml`

- Python: `3.12`
- OpenMM: `8.5.2`
- MDTraj: `1.11.1`
- PDBFixer: `1.12`
- NumPy: environment-managed variant
- SciPy: environment-managed variant
- matplotlib: `3.11.0`
- seaborn: `0.13.2`
- netCDF4: `1.7.1.post1`

This environment supports explicit-solvent molecular dynamics refinement and downstream energetic analysis. It is among the most hardware-sensitive components of the workflow and is intended for execution on a Linux workstation or cluster node with a compatible NVIDIA CUDA stack.

---

## Limitations and Operating Constraints

This workflow is designed for a defined operating envelope and should be interpreted within that scope.

- Linux is the supported operating environment for the current setup;
- CUDA-capable GPU hardware is strongly recommended for the Boltz and molecular-dynamics stages; CPU-only execution may be possible for limited tests but is not the intended production mode;
- the receptor structure should be supplied as a reasonably complete and well-resolved CIF file;
- the pipeline is intended primarily for peptide–protein interaction studies and is not a general-purpose solution for all biomolecular systems without additional configuration and validation;
- the docking and MD modules are expected to perform most reliably when the target has a well-defined binding pocket or interface region;
- the implementation depends on the exact environment files and package versions included in the repository; mixing or substituting versions outside these files is not guaranteed to reproduce the same behavior;
- some benchmarking and bypass controls are tuned for specific internal validation workflows and may require reconfiguration for unrelated targets or study designs.

### Conditions Under Which the Workflow Is Expected to Perform Most Reliably

The workflow is expected to be most reliable when the following conditions are met:

1. the receptor is a well-resolved protein structure with a defined binding interface or pocket;
2. the peptide length and design constraints are compatible with a short-to-mid-length peptide discovery setting, typically within the range of approximately 8–20 residues;
3. the sequence design objective is aligned with a specific physicochemical target, such as amphipathic helix design, net charge tuning, or hydrophobic moment optimization;
4. the target site is sufficiently tractable for contact-based filtering and structural validation;
5. the study is designed as a reproducible, system-specific computational workflow for peptide discovery and comparative analysis;
6. the project is executed on a GPU-enabled Linux workstation or cluster using the repository-specified conda environments.

Under these conditions, the workflow should be understood as a reproducible and structured method for peptide–target investigation rather than as a universally generalizable predictor across all receptor classes or experimental systems.

---

## Outputs and Interpretation

All outputs are organized under `results/<run_name>/`.

### Key Output Folders

- `library.fasta`: generated peptide library;
- `top_designs/`: top structural candidates and validation outputs;
- `md_simulations/`: molecular-dynamics trajectory and topology files;
- `mmgbsa_results/`: energetic analysis outputs;
- `final_summary/`: aggregate figures and final summary tables.

### Representative Outputs

- candidate-level summary tables;
- per-sequence structural and physicochemical metrics;
- docking quality indicators;
- molecular-dynamics stability traces;
- interaction and energetic summaries;
- final comparative report for downstream analysis.

---

## Recommended Manuscript Reporting Structure

For a publication, the final analysis should report the following:

1. exact target used and receptor metadata;
2. peptide design constraints and selection logic;
3. sequence generation settings and seed values;
4. structure-prediction filtering thresholds;
5. docking protocol and target constraints;
6. molecular-dynamics simulation settings and analysis window;
7. energetic calculations and criteria used;
8. final output directory path and reproducible run identifier.

This allows the computational workflow to be reproduced from the repository in a transparent and auditable manner.

---

## License

This project is provided for research use. Please review the repository license and cite the project appropriately when using it in published work.

---

## Citation

If this workflow contributes to your research, please cite both the software and the relevant manuscript or benchmark study using the project details and software versions corresponding to the experiment.

---

## Notes

This repository is intended for computational biology and protein-design research. Results should be interpreted in conjunction with experimental validation, particularly in the context of peptide binders intended for therapeutic or mechanistic applications.
