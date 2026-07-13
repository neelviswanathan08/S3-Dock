import os
import sys
import yaml
import random
import argparse
import peptides

print("====================================================")
print(" S3-DOCK: EVOLUTIONARY SEQUENCE OPTIMIZER (v2)")
print("====================================================")

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
run_dir = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results", run_name))
os.makedirs(run_dir, exist_ok=True)
fasta_bridge_path = os.path.join(run_dir, "library.fasta")

# Check Benchmark Gatekeeper
if config.get('benchmark_mode', False):
    benchmark_seq = config.get('benchmark_sequence', "").strip().upper()
    if not benchmark_seq:
        print("[ERROR] benchmark_mode is True but benchmark_sequence is empty!", flush=True)
        sys.exit(1)
        
    print(f"[BENCHMARK] Bypassing optimizer. Routing explicit control sequence.", flush=True)
    with open(fasta_bridge_path, "w") as f:
        f.write(f">Benchmark_Control\n{benchmark_seq}\n")
    print("[SUCCESS] Phase 1 Complete: Control assigned to library.fasta", flush=True)
    sys.exit(0)

# --- ORIGINAL EVOLUTIONARY OPTIMIZER CONTINUES BELOW ---
import random
import math
import argparse
import peptides

parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, default=None)
cli_args = parser.parse_args()

num_seeds = config['num_seeds']
peptide_len = config['peptide_length']
custom_pattern = config.get('custom_pattern', "")
rng_seed = cli_args.seed if cli_args.seed is not None else config.get('rng_seed', 42)
random.seed(rng_seed)

print(f" Evolving {num_seeds} highly-optimized candidates for run: {run_name}...")
print(f" RNG seed: {rng_seed} (reproducible; override with --seed)")

if config.get('secondary_structure') == "helix":
    H_AAs = ['A', 'L', 'M', 'F', 'I']
    P_AAs = ['E', 'Q', 'K', 'R', 'H']
    calc_angle = 100
elif config.get('secondary_structure') == "beta_sheet":
    H_AAs = ['V', 'I', 'Y', 'W', 'F', 'C']
    P_AAs = ['T', 'S', 'N', 'D', 'P', 'G']
    # NOTE: 180 is the conventional alternating-face angle for extended/beta
    # strands in the Eisenberg moment formalism (100 for helix, 180 for
    # strand). 160 was used previously — confirm which convention you are
    # citing before publishing; both appear in different tools' defaults.
    calc_angle = config.get('beta_sheet_angle', 160)
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
    if len(custom_pattern) > peptide_len: custom_pattern = custom_pattern[:peptide_len]
    elif len(custom_pattern) < peptide_len: custom_pattern = (custom_pattern * (peptide_len // len(custom_pattern) + 1))[:peptide_len]

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
    return "".join(random.choice(H_AAs) if ch == 'H' else random.choice(P_AAs) if ch == 'P' else random.choice(ALL_AAs) for ch in pattern)


# ---------------------------------------------------------------------------
# 5. SIMULATED-ANNEALING HILL CLIMBER
#    - accepts strictly-improving moves always
#    - accepts non-improving moves with probability exp(-delta/T), T decaying
#      each attempt, so early exploration is looser and it tightens over time
#    - if stuck on a plateau for `plateau_patience` steps with no improvement,
#      forces a multi-site perturbation to escape local optima
# ---------------------------------------------------------------------------
valid_sequences = []
run_stats = []          # per-seed summary, for run_metadata.json
trajectory_rows = []    # per-attempt fitness log, for convergence_log.csv

for seed_idx in range(num_seeds):
    current_seq = random_seq_for_pattern(custom_pattern)
    current_fitness, c, b, m = get_fitness(current_seq)
    best_fitness = current_fitness
    attempts, since_improvement = 0, 0
    temperature = initial_temp

    while current_fitness > 0 and attempts < max_mutations:
        attempts += 1
        temperature *= cooling_rate
        mut_idx = random.randint(0, peptide_len - 1)
        char_type = custom_pattern[mut_idx]
        new_aa = random.choice(H_AAs) if char_type == 'H' else random.choice(P_AAs) if char_type == 'P' else random.choice(ALL_AAs)
        new_seq = current_seq[:mut_idx] + new_aa + current_seq[mut_idx + 1:]
        n_fitness, n_c, n_b, n_m = get_fitness(new_seq)
        delta = n_fitness - current_fitness
        if delta <= 0 or random.random() < math.exp(-delta / max(temperature, 1e-6)):
            current_seq, current_fitness, c, b, m = new_seq, n_fitness, n_c, n_b, n_m
        if current_fitness < best_fitness:
            best_fitness = current_fitness
            since_improvement = 0
        else:
            since_improvement += 1

        # Escape plateaus: force a small multi-site jump rather than stalling
        if since_improvement >= plateau_patience:
            perturb_positions = random.sample(range(peptide_len), k=min(3, peptide_len))
            seq_chars = list(current_seq)
            for p in perturb_positions:
                char_type = custom_pattern[p]
                seq_chars[p] = random.choice(H_AAs) if char_type == 'H' else random.choice(P_AAs) if char_type == 'P' else random.choice(ALL_AAs)
            current_seq = "".join(seq_chars)
            current_fitness, c, b, m = get_fitness(current_seq)
            since_improvement = 0

    if current_fitness == 0:
        valid_sequences.append(current_seq)
        print(f"  [+] Optimized Seed {seed_idx+1}: {current_seq} "
              f"(Charge: {c}, Moment: {m:.2f}, Boman: {b:.2f}) | {attempts} mutations")
    else:
        print(f"  [!] Failed to optimize Seed {seed_idx+1} after {max_mutations} mutations. "
              f"(Closest found — Charge: {c}, Boman: {b:.2f}, Moment: {m:.2f})")

    run_stats.append({
        "seed_idx": seed_idx + 1, "success": success, "attempts": attempts,
        "final_charge": c, "final_boman": round(b, 3), "final_moment": round(m, 3),
        "final_fitness": round(current_fitness, 3)
    })

if len(valid_sequences) == 0:
    print(f"\n❌ CRITICAL: The requested biochemical parameters are mathematically "
          f"impossible for a {peptide_len}-mer given the current pattern.")
    print("   Try lowering 'boman_index_min' or 'hydrophobic_moment_target', or "
          "adjusting 'target_net_charge' / 'charge_tolerance' in config.yaml.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# 6. DEDUPLICATION + DIVERSITY METRICS
#    A "library" of optimized candidates should be reported as distinct
#    sequences, with a diversity figure so reviewers can see the optimizer
#    isn't just rediscovering one solution repeatedly.
# ---------------------------------------------------------------------------
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
    return total / count / peptide_len  # normalized 0-1


diversity_score = mean_pairwise_hamming(valid_sequences)

print(f"\n📊 Dedup: {n_before_dedup} → {n_after_dedup} unique sequences.")
if diversity_score is not None:
    print(f"📊 Mean normalized pairwise Hamming distance: {diversity_score:.3f} "
          f"(0 = identical, 1 = maximally different)")

# ---------------------------------------------------------------------------
# 7. Export FASTA (same format/location as before — downstream folding
#    pipeline is unaffected)
# ---------------------------------------------------------------------------
with open(fasta_bridge_path, "w") as f:
    for idx, sequence in enumerate(valid_sequences, 1): f.write(f">Seed_{idx}\n{sequence}\n")

# ---------------------------------------------------------------------------
# 8. Export convergence trajectory (for methods-section figures)
# ---------------------------------------------------------------------------
if trajectory_rows:
    with open(trajectory_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(trajectory_rows[0].keys()))
        writer.writeheader()
        writer.writerows(trajectory_rows)

# ---------------------------------------------------------------------------
# 9. Export run metadata (for reproducibility / supplementary materials)
# ---------------------------------------------------------------------------
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

print(f"\n✅ STAGE 1 COMPLETE: {n_after_dedup} unique optimized sequences → "
      f"results/{run_name}/library.fasta")
print(f"   Convergence log  → results/{run_name}/convergence_log.csv")
print(f"   Run metadata     → results/{run_name}/run_metadata.json")