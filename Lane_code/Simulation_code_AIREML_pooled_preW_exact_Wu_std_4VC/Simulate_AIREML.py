from Function_AIREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score traces tr(V^{-1} K_i) are EXACT (dense Cholesky of V + potri), so
# there is no --nmc, no probe seed, and no truncation --r: W is always the
# exact kernel the phenotype was drawn from.
#
# Print the AI-REML trace (one line per iteration: the four components, rho,
# the exact log-likelihood gain) to STDOUT.  Off by default.
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
rep = args.rep
mode = args.mode
DIR = "/home/ziyanzha/MOM_within_gene/AIREML_pooled_preW_exact_Wu_std_4VC"

# Load the GENOTYPE and build the two standardized designs.  ai_reml forms
# K_a = Z_a Z_a'/m and K_d = Z_d Z_d'/m densely from them on entry.  The
# dominance design is built from the RAW 0/1/2 dosages -- it needs the allele
# frequencies -- by the SAME transformation the Cholesky job used, so K_d is
# identical on both sides.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)
Zd = dominance_design(SNP)

# Load the PRECOMPUTED epistasis kernel the Cholesky job wrote
# (W/W_<mode>_n<n>_m<m>_G<G>.npy plus a JSON sidecar).  It is the
# UNSTANDARDIZED, C-NORMALIZED kernel with the 1/c-hat already inside, so
# nothing is divided here.
#
# c-hat is RECOMPUTED from the genotype (O(nm), one mat-vec per gene) purely as
# a consistency check: load_W_cache refuses the kernel if the sidecar's c-hat
# disagrees, which catches a kernel built from a different genotype file under
# the same name -- or one whose (mode, n, m, G) does not match these arguments.
c_check = pooled_c(split_into_genes(Z, G))
w_npy, w_json = w_cache_paths(f"{DIR}/W", mode, n, m, G)
W, w_meta = load_W_cache(w_npy, w_json, mode, n, m, G, c_check)

# Load phenotype (s2a_s2d_s2gxg_s2e order)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"{DIR}/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Exact AI-REML: V = s2a K_a + s2d K_d + s2gxg W + s2e I, all four dense, V
# Cholesky-factorized once per score/AI evaluation.
#
# Time the estimation only (the AI_REML call) for the per-rep timing record.
# It INCLUDES forming K_a and K_d densely (O(n^2 m), inside ai_reml) and
# EXCLUDES loading W from disk and building it in the Cholesky job (the latter
# is in time/result/W_timing_<FILENAME>.txt).
t_start = time.perf_counter()
if args.verbose:
    print(f"--- rep{rep}: exact AI-REML trace, var(y)={y.var():.6f}, "
          f"psd={w_meta['psd']}, columns s2a s2d s2gxg s2e ---", flush=True)
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, AI = AI_REML(Z, Zd, W, y, iters=iters,
                                                   verbose=args.verbose,
                                                   w_psd=bool(w_meta['psd']))
elapsed = time.perf_counter() - t_start
# How much dense linear algebra that fit cost.  ai_reml zeroes the counters on
# entry, so this snapshot is THIS replicate and nothing else.
op_counts = get_op_counts()

# NO post-fit realized-variance correction: the same c-hat is already INSIDE
# the cached kernel, so s2gxg_hat is read directly on the realized-variance
# scale.  Applying compute_c_pooled here as well would double-count.

# Save result, 4 columns -- the FIT and nothing else:
#
#   (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat)
#    V_a      V_d      V_gamma    V_e
run_tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
output_dir = f"{DIR}/result/{run_tag}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"{DIR}/time/rep_times/{run_tag}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")

# Record this replicate's COST COUNTS, one file per replicate, as labelled
# "<name> <value>" lines in a FIXED order with 0 defaults.  The combine step
# averages all reps into time/result/opcount_<FILENAME>.txt.
OP_KEYS = ("reml_iters", "V_factorizations", "rejected_steps")
op_dir = f"{DIR}/time/op_counts/{run_tag}"
os.makedirs(op_dir, exist_ok=True)
with open(f"{op_dir}/rep{rep}.txt", 'w') as f:
    for key in OP_KEYS:
        f.write(f"{key} {op_counts.get(key, 0)}\n")
print(f"rep{rep}: s-hat=({s2a_hat:.6f}, {s2d_hat:.6f}, {s2gxg_hat:.6f}, "
      f"{s2e_hat:.6f})  SE={np.round(reml_se(AI), 6)}  "
      f"{op_counts.get('reml_iters', 0)} REML iters, "
      f"{op_counts.get('V_factorizations', 0)} Cholesky factorizations of V, "
      f"{op_counts.get('rejected_steps', 0)} rejected steps, "
      f"{elapsed:.3f} s -> {op_dir}/rep{rep}.txt")
