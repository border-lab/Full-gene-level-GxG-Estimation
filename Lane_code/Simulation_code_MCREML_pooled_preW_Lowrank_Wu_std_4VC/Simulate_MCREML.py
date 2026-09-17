from Function_MCREML import *
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
parser.add_argument('--nmc', type=int, default=100)   # FINE phase; coarse is 15
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score traces tr(V^{-1} K_i) are Hutchinson's, with Nmc CG solves per
# REML iteration -- the ONLY estimator; there is no trace_method to pass.
# Which PRECOMPUTED epistasis kernel to fit with -- the one knob this variant
# has, and it must match the Cholesky job's --r:
#   r <= 0  (default)  the EXACT c-normalized W the phenotype was drawn from.
#                      Simulation and estimation share the matrix; no
#                      truncation bias.
#   r >  0             the dense matrix of the _matfree_ sibling's rank-r
#                      operator.  Same deterministic truncation bias as the
#                      sibling at that r, applied as a dense gemm instead.
# It applies to the EPISTASIS kernel ONLY: K_a = Z_a Z_a'/m and K_d = Z_d Z_d'/m
# are applied exactly and matrix-free at every r.
parser.add_argument('--r', type=int, default=R_DEFAULT)
# Print the AI-REML trace (one line per iteration: the four components and
# max|step|) to STDOUT.  Off by default -- with it on, Step 3 must also be
# launched with a real --output, since the pipeline sends stdout to /dev/null.
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
nmc = args.nmc
rep = args.rep
mode = args.mode
r = args.r
DIR = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC"

# Load the GENOTYPE and build the two standardized designs.  K_a and K_d are
# applied MATRIX-FREE from them, exactly as in the _matfree_ sibling: the
# additive GRM as K_a B = Z_a(Z_a'B)/m, the dominance GRM as K_d B = Z_d(Z_d'B)/m.
# The dominance design is built from the RAW 0/1/2 dosages -- it needs the
# allele frequencies -- and is the SAME transformation the Cholesky job used, so
# K_d is identical on both sides.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)
Zd = dominance_design(SNP)

# Load the PRECOMPUTED epistasis kernel -- THE difference from the _matfree_
# sibling, which rebuilds W u from per-gene SVD factors on every CG iteration.
# The Cholesky job wrote it (W/W_<mode>_n<n>_m<m>_G<G>_<rexact|r<r>>.npy plus a
# JSON sidecar); it is the UNSTANDARDIZED, C-NORMALIZED kernel with the 1/c-hat
# already inside, so nothing is divided here.
#
# c-hat is RECOMPUTED from the genotype (the O(nm) third-moment plug-in,
# C_METHOD, one mat-vec per gene) purely as a consistency check: load_W_cache
# refuses the kernel if the sidecar's c-hat disagrees, which is what catches a
# kernel built from a different genotype file under the same name -- or one
# whose (mode, n, m, G, r) does not match these arguments.
c_check = pooled_c(split_into_genes(Z, G))
w_npy, w_json = w_cache_paths(f"{DIR}/W", mode, n, m, G, r)
W, w_meta = load_W_cache(w_npy, w_json, mode, n, m, G, r, c_check)

# Load phenotype (s2a_s2d_s2gxg_s2e order)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"{DIR}/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# NOTHING IS READ BACK FROM THE PHENOTYPE STEP EXCEPT y.  simulate_phenotype
# rescales each component to hit its target exactly (force_realized=True), so
# all four realized variances are the nominal targets in every replicate and
# there is nothing per-replicate to carry.

# Monte-Carlo AI-REML: V = s2a K_a + s2d K_d + s2gxg W + s2e I, W the dense
# precomputed kernel, K_a and K_d matrix-free.  seed=rep fixes the Hutchinson
# probe randomness per replicate (reproducible, yet varied).  w_psd comes from
# the sidecar: the exact kernel is PSD by construction and gets lam_min(W) = 0
# in the feasibility bound, the truncated one has both ends estimated.
#
# Time the estimation only (the MC_REML call) for the per-rep timing record.
# Loading W from disk and building it in the Cholesky job are NOT in this
# number -- the latter is in time/result/W_est_timing_<FILENAME>.txt -- whereas
# the _matfree_ sibling's number INCLUDES its per-gene SVD setup.  Keep that in
# mind before comparing the two timing files.
t_start = time.perf_counter()
if args.verbose:
    print(f"--- rep{rep}: AI-REML trace, var(y)={y.var():.6f}, "
          f"kernel={w_cache_tag(r)} (psd={w_meta['psd']}), "
          f"columns s2a s2d s2gxg s2e ---", flush=True)
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, _ = MC_REML(Z, Zd, W, y, iters=iters,
                                                  Nmc=nmc, seed=rep,
                                                  verbose=args.verbose,
                                                  w_psd=bool(w_meta['psd']))
elapsed = time.perf_counter() - t_start
# How many LINEAR OPERATOR APPLIES that fit actually cost.  mc_reml zeroes the
# counters on entry, so this snapshot is THIS replicate and nothing else.  The
# keys are the sibling's, so the two op-count files compare line for line: the
# COUNTS should be close (same optimizer, same V up to the truncation), and the
# per-apply cost is what differs -- O(n^2) per W column here, O(n m r) there.
op_counts = get_op_counts()

# NO post-fit realized-variance correction: the same c-hat is already INSIDE
# the cached kernel, so s2gxg_hat is read directly on the realized-variance
# scale.  Applying compute_c_pooled here as well would double-count the
# normalization and inflate the epistasis estimate by a factor of c-hat.

# Save result, 4 columns -- the FIT and nothing else:
#
#   (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat)
#    V_a      V_d      V_gamma    V_e
#
# FILENAME carries the kernel tag (rexact / r<r>) so an exact run and a
# truncated run at the same targets land in different files.
run_tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}_{w_cache_tag(r)}"
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

# Record this replicate's OPERATOR-APPLY COUNTS, one file per replicate, as
# labelled "<name> <value>" lines in a FIXED order with 0 defaults.  The combine
# step averages all reps into time/result/opcount_<FILENAME>.txt.
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "K_applies", "K_columns",
           "W_applies", "W_columns")
op_dir = f"{DIR}/time/op_counts/{run_tag}"
os.makedirs(op_dir, exist_ok=True)
with open(f"{op_dir}/rep{rep}.txt", 'w') as f:
    for key in OP_KEYS:
        f.write(f"{key} {op_counts.get(key, 0)}\n")
print(f"operator applies this replicate: V={op_counts.get('V_applies', 0)} "
      f"({op_counts.get('V_columns', 0)} columns), "
      f"K={op_counts.get('K_applies', 0)} "
      f"({op_counts.get('K_columns', 0)} columns), "
      f"W={op_counts.get('W_applies', 0)} "
      f"({op_counts.get('W_columns', 0)} columns); "
      f"{op_counts.get('reml_iters', 0)} REML iters, "
      f"{op_counts.get('cg_iters', 0)} CG iters in "
      f"{op_counts.get('cg_solves', 0)} solves -> {op_dir}/rep{rep}.txt")
