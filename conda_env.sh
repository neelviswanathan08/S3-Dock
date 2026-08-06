#!/bin/bash

# Load Conda functions securely
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/.bashrc

# ---------------------------------------------------------------------
# AUTOMATED ENVIRONMENT BUILDER (Flawless First-Time Setup)
# ---------------------------------------------------------------------
echo " Checking and verifying local runtime environments..."

# 1. Boltz Environment
if [ ! -d "envs/boltz_env" ]; then
    echo " Creating local Boltz environment from envs/ directory..."
    conda env create --prefix envs/boltz_env -f envs/boltz_env.yml -y
fi

# 2. HADDOCK3 Environment
if [ ! -d "envs/haddock_env" ]; then
    echo " Creating local HADDOCK3 environment from envs/ directory..."
    conda env create --prefix envs/haddock_env -f envs/haddock_env.yml -y
fi

# 3. OpenMM Physics Environment
if [ ! -d "envs/openmm_env" ]; then
    echo " Creating local OpenMM environment from envs/ directory..."
    conda env create --prefix envs/openmm_env -f envs/openmm_env.yml -y
fi

echo "✅ All environments verified and locked locally."
echo "---------------------------------------------------------------------"
cd scripts/ || { echo "❌ ERROR: scripts/ folder missing!"; exit 1; }

