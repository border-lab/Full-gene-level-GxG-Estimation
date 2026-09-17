from Function_MCREML import *
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode


# Read genotype (m SNPs), then split into G contiguous genes inside
# simulate_Cholesky_gxg (several Z, one per gene) for the POOLED kernel.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) for drawing the epistasis effect.
# The dense pooled W is built here because a Cholesky factor needs an explicit
# matrix -- but UNLIKE the _preW pipeline it is NOT cached: estimation rebuilds
# the action of W matrix-free from the genotype, so no n-by-n array is ever
# written to disk or loaded downstream.  w_build_time is the wall-clock cost of
# building the kernel W (build_W_pooled only), tracked separately from the
# per-rep estimation time.
Lgxg, w_build_time = simulate_Cholesky_gxg(SNP, G, s2gxg=s2gxg, s2e=s2e)

# Save the Cholesky factor (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_new/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Record the kernel-build wall-clock time (seconds).  W is built once by this
# single Cholesky job, so there is nothing to average -- write it straight to
# time/result/W_timing_<FILENAME>.txt (parallels estimation's timing_<...>.txt).
# NOTE this is SIMULATION-ONLY cost; the estimator never pays it.
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_new/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time (simulation only): {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")
