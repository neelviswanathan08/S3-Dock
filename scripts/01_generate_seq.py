import os
import sys
import yaml

print("====================================================", flush=True)
print("[SYSTEM] S3-DOCK: PHASE 1 - SEQUENCE GENERATION ENGINE", flush=True)
print("====================================================", flush=True)

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

print(f"[INFO] Evolving {num_seeds} optimized candidates for run: {run_name}...", flush=True)

if config.get('secondary_structure') == "helix":
    H_AAs, P_AAs, calc_angle = ['A', 'L', 'M', 'F', 'I'], ['E', 'Q', 'K', 'R', 'H'], 100
elif config.get('secondary_structure') == "beta_sheet":
    H_AAs, P_AAs, calc_angle = ['V', 'I', 'Y', 'W', 'F', 'C'], ['T', 'S', 'N', 'D', 'P', 'G'], config.get('beta_sheet_angle', 160)
else:
    H_AAs, P_AAs, calc_angle = ['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y'], ['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T'], 100

ALL_AAs = H_AAs + P_AAs

if not custom_pattern:
    pattern_list = []
    if config.get('secondary_structure') == "helix" and config.get('peptide_property') == "amphipathic":
        for i in range(peptide_len): pattern_list.append('H' if i % 7 in [0, 3, 4] else 'P')
    elif config.get('secondary_structure') == "beta_sheet" and config.get('peptide_property') == "amphipathic":
        for i in range(peptide_len): pattern_list.append('H' if i % 2 == 0 else 'P')
    elif config.get('peptide_property') == "hydrophobic": pattern_list = ['H'] * peptide_len
    elif config.get('peptide_property') == "hydrophilic": pattern_list = ['P'] * peptide_len
    else: pattern_list = ['X'] * peptide_len
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

valid_sequences = []
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

if len(valid_sequences) == 0:
    print("\n[ERROR] Parameters mathematically impossible.", flush=True)
    sys.exit(1)

if dedupe_output: valid_sequences = list(set(valid_sequences))
with open(fasta_bridge_path, "w") as f:
    for idx, sequence in enumerate(valid_sequences, 1): f.write(f">Seed_{idx}\n{sequence}\n")

print(f"\n[SUCCESS] Phase 1 Complete: {len(valid_sequences)} unique optimized seeds saved.", flush=True)