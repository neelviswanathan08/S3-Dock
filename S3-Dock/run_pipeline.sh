#!/bin/bash

# Find directory where this script sits
REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$REPO_DIR"

# GLOBAL OVERRIDE: Absolutely forbid Python and the OS from buffering logs
export PYTHONUNBUFFERED=1

echo "====================================================================="
echo " S3-DOCK: HIGH-THROUGHPUT DISCOVERY PIPELINE RUNNING"
echo "====================================================================="
echo "Started at: $(date)"
echo "---------------------------------------------------------------------"

if [ ! -d "envs/boltz_env" ] || [ ! -d "envs/haddock_env" ] || [ ! -d "envs/openmm_env" ]; then
    echo "[CRITICAL ERROR] Local runtime environments are missing!"
    exit 1
fi

RUN_NAME=$(grep '^run_folder_name:' config.yaml | sed 's/.*: *//' | sed 's/ *#.*//' | tr -d '""'"'' ")
RUN_DIR="results/${RUN_NAME}"
MAX_LOOPS=$(grep '^max_discovery_loops:' config.yaml | sed 's/.*: *//' | sed 's/ *#.*//' | tr -d '""'"'' ")
if [ -z "$MAX_LOOPS" ]; then MAX_LOOPS=5; fi # Fallback if not in config

echo "[INFO] Active Workspace: ${RUN_DIR}"
echo "---------------------------------------------------------------------"

# --- GLOBAL GUARD: Skip all if final summary already exists ---
if [ -d "${RUN_DIR}/final_summary" ] && ls ${RUN_DIR}/final_summary/*.csv 1> /dev/null 2>&1; then
    echo "[ALL STAGES SKIPPED] Complete final_summary data already exists!"
    echo "====================================================================="
    exit 0
fi

# =====================================================================
# DISCOVERY QUEUE LOOP (Phases 1 & 2)
# =====================================================================
LOOP_COUNT=1
HAS_CANDIDATES=false

while [ $LOOP_COUNT -le $MAX_LOOPS ]; do
    if [ $MAX_LOOPS -gt 1 ]; then
        echo "---------------------------------------------------------------------"
        echo "[DISCOVERY LOOP] Attempt $LOOP_COUNT of $MAX_LOOPS"
        echo "---------------------------------------------------------------------"
    fi

    echo "[Phase 1] Generating Sequence Library..."
    # Smart resume only triggers on the FIRST loop if files already exist
    if [ $LOOP_COUNT -eq 1 ] && [ -f "${RUN_DIR}/library.fasta" ]; then
        echo "   [Smart Resume] Found existing library.fasta. Skipping Phase 1."
    else
        # Pass a random seed to the optimizer on loops > 1 so we get brand new sequences
        SEED_ARG=""
        if [ $LOOP_COUNT -gt 1 ]; then
            SEED_ARG="--seed $RANDOM"
        fi
        
        ./envs/boltz_env/bin/python -u scripts/01_generate_seq.py $SEED_ARG
        if [ $? -ne 0 ]; then echo "[ERROR] Phase 1 failed!"; exit 1; fi
    fi

    echo "[Phase 2] Boltz-2 Structural Folding..."
    if [ $LOOP_COUNT -eq 1 ] && [ -d "${RUN_DIR}/top_designs" ] && ls ${RUN_DIR}/top_designs/*_best.cif 1> /dev/null 2>&1; then
        echo "   [Smart Resume] Found existing folded candidates. Skipping Phase 2."
    else
        ./envs/boltz_env/bin/python -u scripts/02_predict_struct.py
        if [ $? -ne 0 ]; then echo "[ERROR] Phase 2 failed!"; exit 1; fi
    fi

    # --- GATEKEEPER CHECK ---
    if ls ${RUN_DIR}/top_designs/*_best.cif 1> /dev/null 2>&1; then
        HAS_CANDIDATES=true
        break # We found a winner! Break out of the loop and continue to Phase 3.
    else
        echo "[DISCOVERY GATEKEEPER] No sequences survived Phase 2 validation."
        if [ $LOOP_COUNT -lt $MAX_LOOPS ]; then
            echo "[RETRY] Wiping failed candidates and triggering next Discovery Loop..."
            rm -f "${RUN_DIR}/library.fasta"
            rm -rf "${RUN_DIR}"/Seed_* # Clean out the failed Boltz directories
        else
            echo "[STOP] Maximum discovery loops ($MAX_LOOPS) reached without success."
            echo "[SKIP] Gracefully skipping heavy docking and MD (Phases 3, 4, 5)."
        fi
    fi

    LOOP_COUNT=$((LOOP_COUNT + 1))
done

# =====================================================================
# DEEP VALIDATION (Phases 3, 4, 5)
# Only executes if a candidate survived the Discovery Queue
# =====================================================================
if [ "$HAS_CANDIDATES" = true ]; then
    echo "---------------------------------------------------------------------"
    echo "[Phase 3] Unbiased Blind Global Docking..."
    if [ -d "${RUN_DIR}/haddock_runs" ] && ls ${RUN_DIR}/haddock_runs/*/haddock3_output 1> /dev/null 2>&1; then
        echo "   [Smart Resume] Found existing HADDOCK3 outputs. Skipping Phase 3."
    else
        ./envs/haddock_env/bin/python -u scripts/03_blind_dock.py
        if [ $? -ne 0 ]; then echo "[ERROR] Phase 3 failed!"; exit 1; fi
    fi

    echo "[Phase 4] OpenMM Molecular Dynamics Simulation..."
    if [ -d "${RUN_DIR}/md_simulations" ] && [ "$(ls -A ${RUN_DIR}/md_simulations 2>/dev/null)" ]; then
        echo "   [Smart Resume] Found existing MD directories. Skipping Phase 4."
    else
        ./envs/openmm_env/bin/python -u scripts/04_md_simulate.py
        if [ $? -ne 0 ]; then echo "[ERROR] Phase 4 failed!"; exit 1; fi
    fi

    echo "[Phase 5] Calculating Native OpenMM MM-GBSA Free Energy..."
    if [ -d "${RUN_DIR}/mmgbsa_results" ] && ls ${RUN_DIR}/mmgbsa_results/*_mmgbsa.csv 1> /dev/null 2>&1; then
        echo "   [Smart Resume] Found existing Free Energy calculations. Skipping Phase 5."
    else
        ./envs/openmm_env/bin/python -u scripts/05_calculate_energy.py
        if [ $? -ne 0 ]; then echo "[ERROR] Phase 5 failed!"; exit 1; fi
    fi
fi

echo "---------------------------------------------------------------------"
echo "[Phase 6] Compiling Final Master Discovery Report..."
# If Phase 2 failed on all loops, this will cleanly wrap up whatever metadata is left
./envs/boltz_env/bin/python -u scripts/06_compile_summary.py

echo "---------------------------------------------------------------------"
echo "[SUCCESS] PIPELINE COMPLETE! Workspace successfully processed."
echo "Ended at: $(date)"
echo "====================================================================="