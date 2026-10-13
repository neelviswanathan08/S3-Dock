import os
import csv
import numpy as np
from scipy.stats import spearmanr

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------
# This matches the experimentally validated MDM2 PMI dataset ranking 
# (Rank 1 = Strongest affinity / lowest Kd, Rank 15 = Weakest affinity / highest Kd)
EXPERIMENTAL_BENCHMARK = [
    {"mutation": "N8A", "sequence": "TSFAEYWALLSP", "true_rank": 1, "kd_nm": 0.49},
    {"mutation": "P12A", "sequence": "TSFAEYWNLLSA", "true_rank": 2, "kd_nm": 2.1},
    {"mutation": "L9A", "sequence": "TSFAEYWNALSP", "true_rank": 3, "kd_nm": 2.4},
    {"mutation": "WT_PMI", "sequence": "TSFAEYWNLLSP", "true_rank": 4, "kd_nm": 3.2},
    {"mutation": "S11A", "sequence": "TSFAEYWNLLAP", "true_rank": 5, "kd_nm": 3.9},
    {"mutation": "T1A", "sequence": "ASFAEYWNLLSP", "true_rank": 6, "kd_nm": 6.2},
    {"mutation": "trunc_1", "sequence": "TSFAEYWNLL", "true_rank": 7, "kd_nm": 8.6},
    {"mutation": "E5A", "sequence": "TSFAAYWNLLSP", "true_rank": 8, "kd_nm": 21.0},
    {"mutation": "S2A", "sequence": "TAFAEYWNLLSP", "true_rank": 9, "kd_nm": 27.0},
    {"mutation": "trunc_2", "sequence": "SFAEYWNLL", "true_rank": 10, "kd_nm": 170.0},
    {"mutation": "Y6A", "sequence": "TSFAEAWNLLSP", "true_rank": 11, "kd_nm": 610.0},
    {"mutation": "L10A", "sequence": "TSFAEYWNLASP", "true_rank": 12, "kd_nm": 890.0},
    {"mutation": "trunc_3", "sequence": "FAEYWNLLSP", "true_rank": 13, "kd_nm": 6500.0},
    {"mutation": "F3A", "sequence": "TSAAEYWNLLSP", "true_rank": 14, "kd_nm": 38000.0},
    {"mutation": "W7A", "sequence": "TSFAEYANLLSP", "true_rank": 15, "kd_nm": 160000.0},
]

def evaluate_pyrosetta_rankings(results_csv_path):
    print("====================================================", flush=True)
    print(" S3-DOCK: SPEARMAN RANK CORRELATION BENCHMARK EVALUATOR", flush=True)
    print("====================================================", flush=True)

    if not os.path.exists(results_csv_path):
        print(f"❌ [ERROR] Results CSV not found at: {results_csv_path}")
        print("   -> Ensure Phase 2/PyRosetta scoring has generated outputs.")
        return

    # Ingest pipeline scores mapped by sequence
    pipeline_scores = {}
    with open(results_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row.get("Binder_Sequence", "").strip()
            # Look for Rosetta score column (adjust name if your output differs)
            score_val = row.get("Rosetta_dG_kcal_mol") or row.get("Rosetta_Score")
            if seq and score_val and score_val != "N/A":
                try:
                    pipeline_scores[seq] = float(score_val)
                except ValueError:
                    continue

    matched_true_ranks = []
    matched_predicted_scores = []

    print(f"\n[INFO] Aligning experimental benchmark with pipeline results...", flush=True)
    print(f"{'Mut':<8} | {'Sequence':<14} | {'True Rank':<10} | {'Rosetta ΔG (kcal/mol)':<22}")
    print("-" * 62)

    for item in EXPERIMENTAL_BENCHMARK:
        seq = item["sequence"]
        true_rank = item["true_rank"]
        
        # Match by sequence string
        score = pipeline_scores.get(seq)
        if score is not None:
            matched_true_ranks.append(true_rank)
            matched_predicted_scores.append(score)
            print(f"{item['mutation']:<8} | {seq:<14} | {true_rank:<10} | {score:<22.3f}")
        else:
            print(f"{item['mutation']:<8} | {seq:<14} | {true_rank:<10} | {'MISSING IN OUTPUT':<22}")

    if len(matched_true_ranks) < 5:
        print(f"\n❌ [ERROR] Insufficient matched entries ({len(matched_true_ranks)}/15). Check sequence identifiers.")
        return

    # Compute Spearman Rank Correlation
    # Note: Lower (more negative) PyRosetta energy = stronger predicted binding.
    # We negate the scores or let scipy handle directionality. Scipy spearmanr handles monotonic relationships.
    rho, p_value = spearmanr(matched_true_ranks, matched_predicted_scores)

    print("\n----------------------------------------------------", flush=True)
    print(f"📊 STATISTICAL BENCHMARK RESULTS (n = {len(matched_true_ranks)})", flush=True)
    print("----------------------------------------------------", flush=True)
    print(f"   -> Spearman Rank Correlation (ρ) : {rho:.4f}")
    print(f"   -> Statistical Significance (p)  : {p_value:.4e}")
    
    if p_value < 0.05:
        print("   -> 🟢 RESULT: Statistically significant correlation (p < 0.05)!")
    else:
        print("   -> 🟡 RESULT: Correlation not statistically significant at p < 0.05.")
    print("====================================================\n", flush=True)

if __name__ == "__main__":
    # Point directly to your compiled pipeline summary report or phase 2 metrics CSV
    target_csv = "/home/aicd3-2/Desktop/nviswanathan/S3-Dock/results/MAIN/MDM2_Sanity_Check/final_summary/pipeline_summary_report.csv"
    evaluate_pyrosetta_rankings(target_csv)