import os
import sys
import csv
import json
import math
import yaml
import random
import argparse
import peptides

print("===================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 1 - SEQUENCE GENERATION ENGINE", flush=True)
print("====================================================", flush=True)

# ---------------------------------------------------------------------------
# 0. CLI: allow overriding the RNG seed without touching config.yaml
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, default=None,
                     help="Override config.yaml's rng_seed for this run.")
cli_args = parser.parse_args()

# ---------------------------------------------------------------------------
# 1. Setup absolute paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.abspath(os.path.join(SCRIPT_DIR, '..', 'config.yaml'))

with open(CONFIG_PATH, 'r') as file:
    config = yaml.safe_load(file)

run_name = config['run_folder_name']
num_seeds = config['num_seeds']
peptide_len = config['peptide_length']
custom_pattern = config.get('custom_pattern', "")

run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
os.makedirs(run_dir, exist_ok=True)
fasta_bridge_path = os.path.join(run_dir, "library.fasta")
trajectory_path = os.path.join(run_dir, "convergence_log.csv")
metadata_path = os.path.join(run_dir, "run_metadata.json")

# ---------------------------------------------------------------------------
# 1b. REPRODUCIBILITY: seed the RNG and record it
# ---------------------------------------------------------------------------
rng_seed = cli_args.seed if cli_args.seed is not None else config.get('rng_seed', 42)
random.seed(rng_seed)

# ---------------------------------------------------------------------------
# 1c. 🔬 BENCHMARK GATEKEEPER BYPASS
# ---------------------------------------------------------------------------
if config.get('benchmark_mode', False):
    benchmark_seq = config.get('benchmark_sequence', "").strip().upper()
    if not benchmark_seq:
        print("[ERROR] benchmark_mode is True but benchmark_sequence is empty!", flush=True)
        sys.exit(1)
        
    benchmark_count = config.get('benchmark_samples', 1000)
    print(f"[BENCHMARK] Bypassing evolutionary optimizer.", flush=True)
    print(f"[INFO] Cloning explicit target control sequence into {benchmark_count} slots...", flush=True)
    
    with open(fasta_bridge_path, "w") as f:
        for idx in range(1, benchmark_count + 1):
            f.write(f">Seed_{idx}\n{benchmark_seq}\n")
            
    print(f"\n[SUCCESS] STAGE 1 COMPLETE: Fixed benchmark fasta generated at results/{run_name}/library.fasta", flush=True)
    sys.exit(0)

# ---------------------------------------------------------------------------
# 2. STANDARD EVOLUTIONARY OPTIMIZER
# ---------------------------------------------------------------------------
print(f"[INFO] Evolving {num_seeds} highly-optimized candidates for run: {run_name}...", flush=True)
print(f"[INFO] RNG seed: {rng_seed} (reproducible; override with --seed)", flush=True)

if config.get('secondary_structure') == "helix":
    H_AAs = ['A', 'L', 'M', 'F', 'I']
    P_AAs = ['E', 'Q', 'K', 'R', 'H']
    calc_angle = 100
elif config.get('secondary_structure') == "beta_sheet":
    H_AAs = ['V', 'I', 'Y', 'W', 'F', 'C']
    P_AAs = ['T', 'S', 'N', 'D', 'Q']
    calc_angle = config.get('beta_sheet_angle', 180)
else:
    H_AAs = ['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y']
    P_AAs = ['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T']
    calc_angle = 100

ALL_AAs = H_AAs + P_AAs

if not custom_pattern:
    pattern_list = []
    if config.get('secondary_structure') == "helix" and config.get('peptide_property') == "amphipathic":
        for i in range(peptide_len):
            pattern_list.append('H' if i % 7 in [0, 3, 4] else 'P')
    elif config.get('secondary_structure') == "beta_sheet" and config.get('peptide_property') == "amphipathic":
        for i in range(peptide_len):
            pattern_list.append('H' if i % 2 == 0 else 'P')
    elif config.get('peptide_property') == "hydrophobic":
        pattern_list = ['H'] * peptide_len
    elif config.get('peptide_property') == "hydrophilic":
        pattern_list = ['P'] * peptide_len
    else:
        pattern_list = ['X'] * peptide_len
    custom_pattern = "".join(pattern_list)
else:
    if len(custom_pattern) > peptide_len:
        custom_pattern = custom_pattern[:peptide_len]
    elif len(custom_pattern) < peptide_len:
        custom_pattern = (custom_pattern * (peptide_len // len(custom_pattern) + 1))[:peptide_len]

target_charge = config['target_net_charge']
target_boman = config['boman_index_min']
target_moment = config['hydrophobic_moment_target']

charge_tolerance = config.get('charge_tolerance', 0)
charge_weight = config.get('charge_weight', 10)
boman_weight = config.get('boman_weight', 5)
moment_weight = config.get('moment_weight', 5)

max_mutations = config.get('max_mutations', 25000)
plateau_patience = config.get('plateau_patience', 800)
initial_temp = config.get('initial_temp', 1.5)
cooling_rate = config.get('cooling_rate', 0.9995)
dedupe_output = config.get('dedupe_output', True)


def get_fitness(seq):
    pep = peptides.Peptide(seq)
    c = pep.charge(pH=7.4)
    b = pep.boman()
    m = pep.hydrophobic_moment(angle=calc_angle)

    charge_gap = max(0, abs(c - target_charge) - charge_tolerance)
    c_pen = charge_gap * charge_weight
    b_pen = max(0, target_boman - b) * boman_weight
    m_pen = max(0, target_moment - m) * moment_weight

    return (c_pen + b_pen + m_pen), c, b, m


def random_seq_for_pattern(pattern):
    return "".join(
        random.choice(H_AAs) if ch == 'H' else
        random.choice(P_AAs) if ch == 'P' else
        random.choice(ALL_AAs)
        for ch in pattern
    )


valid_sequences = []
run_stats = []          
trajectory_rows = []    

for seed_idx in range(num_seeds):
    current_seq = random_seq_for_pattern(custom_pattern)
    current_fitness, c, b, m = get_fitness(current_seq)
    best_fitness, best_c, best_b, best_m = current_fitness, c, b, m

    attempts = 0
    since_improvement = 0
    temperature = initial_temp

    while current_fitness > 0 and attempts < max_mutations:
        attempts += 1
        temperature *= cooling_rate

        mut_idx = random.randint(0, peptide_len - 1)
        char_type = custom_pattern[mut_idx]
        new_aa = (random.choice(H_AAs) if char_type == 'H' else
                  random.choice(P_AAs) if char_type == 'P' else
                  random.choice(ALL_AAs))
        new_seq = current_seq[:mut_idx] + new_aa + current_seq[mut_idx + 1:]

        n_fitness, n_c, n_b, n_m = get_fitness(new_seq)
        delta = n_fitness - current_fitness

        accept = delta <= 0 or random.random() < math.exp(-delta / max(temperature, 1e-6))
        if accept:
            current_seq = new_seq
            current_fitness = n_fitness
            c, b, m = n_c, n_b, n_m

        if current_fitness < best_fitness:
            best_fitness, best_c, best_b, best_m = current_fitness, c, b, m
            since_improvement = 0
        else:
            since_improvement += 1

        if since_improvement >= plateau_patience:
            perturb_positions = random.sample(range(peptide_len), k=min(3, peptide_len))
            seq_chars = list(current_seq)
            for p in perturb_positions:
                char_type = custom_pattern[p]
                seq_chars[p] = (random.choice(H_AAs) if char_type == 'H' else
                                 random.choice(P_AAs) if char_type == 'P' else
                                 random.choice(ALL_AAs))
            current_seq = "".join(seq_chars)
            current_fitness, c, b, m = get_fitness(current_seq)
            since_improvement = 0

        if attempts % 250 == 0 or current_fitness == 0:
            trajectory_rows.append({
                "seed_idx": seed_idx + 1, "attempt": attempts,
                "fitness": round(current_fitness, 3), "temperature": round(temperature, 4),
                "charge": c, "boman": round(b, 3), "moment": round(m, 3)
            })

    success = current_fitness == 0
    if success:
        valid_sequences.append(current_seq)
        print(f"  [+] Optimized Seed {seed_idx+1}: {current_seq} "
              f"(Charge: {c}, Moment: {m:.2f}, Boman: {b:.2f}) | {attempts} mutations", flush=True)
    else:
        print(f"  [!] Failed to optimize Seed {seed_idx+1} after {max_mutations} mutations. "
              f"(Closest found — Charge: {c}, Boman: {b:.2f}, Moment: {m:.2f})", flush=True)

    run_stats.append({
        "seed_idx": seed_idx + 1, "success": success, "attempts": attempts,
        "final_charge": c, "final_boman": round(b, 3), "final_moment": round(m, 3),
        "final_fitness": round(current_fitness, 3)
    })

if len(valid_sequences) == 0:
    print(f"\n[ERROR] The requested biochemical parameters are mathematically impossible.", flush=True)
    sys.exit(1)

n_before_dedup = len(valid_sequences)
if dedupe_output:
    seen = set()
    deduped = []
    for s in valid_sequences:
        if s not in seen:
            seen.add(s)
            deduped.append(s)
    valid_sequences = deduped
n_after_dedup = len(valid_sequences)


def mean_pairwise_hamming(seqs):
    if len(seqs) < 2:
        return None
    total, count = 0, 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            total += sum(a != b for a, b in zip(seqs[i], seqs[j]))
            count += 1
    return total / count / peptide_len


diversity_score = mean_pairwise_hamming(valid_sequences)

print(f"\n[INFO] Dedup: {n_before_dedup} -> {n_after_dedup} unique sequences.", flush=True)
if diversity_score is not None:
    print(f"[INFO] Mean normalized pairwise Hamming distance: {diversity_score:.3f}", flush=True)

with open(fasta_bridge_path, "w") as f:
    for idx, sequence in enumerate(valid_sequences, 1):
        f.write(f">Seed_{idx}\n{sequence}\n")

if trajectory_rows:
    with open(trajectory_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(trajectory_rows[0].keys()))
        writer.writeheader()
        writer.writerows(trajectory_rows)

metadata = {
    "run_name": run_name,
    "rng_seed": rng_seed,
    "config_used": {
        "peptide_length": peptide_len,
        "custom_pattern": custom_pattern,
        "secondary_structure": config.get('secondary_structure'),
        "peptide_property": config.get('peptide_property'),
        "target_net_charge": target_charge,
        "charge_tolerance": charge_tolerance,
        "boman_index_min": target_boman,
        "hydrophobic_moment_target": target_moment,
        "calc_angle": calc_angle,
        "weights": {"charge": charge_weight, "boman": boman_weight, "moment": moment_weight},
        "max_mutations": max_mutations,
        "plateau_patience": plateau_patience,
        "initial_temp": initial_temp,
        "cooling_rate": cooling_rate,
    },
    "results": {
        "num_seeds_requested": num_seeds,
        "num_successful": sum(1 for r in run_stats if r["success"]),
        "success_rate": sum(1 for r in run_stats if r["success"]) / num_seeds,
        "n_unique_sequences": n_after_dedup,
        "mean_pairwise_hamming_diversity": diversity_score,
    },
    "per_seed_stats": run_stats,
}
with open(metadata_path, "w") as f:
    json.dump(metadata, f, indent=2)

print(f"\n[SUCCESS] STAGE 1 COMPLETE: {n_after_dedup} unique optimized sequences -> results/{run_name}/library.fasta", flush=True)