from Function_MCREML import *
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) for drawing the epistasis effect,
# plus the epistasis GRM W itself.  W here is the MEAN-CENTRED-kernel GRM,
# h_ab = Z_a.Z_b - mean(Z_a.Z_b), centred but with no 1/sigma_ab scaling.  It
# is deterministic per genotype, so it is PRE-COMPUTED and cached once under
# stored_genotype;
# the estimation step loads that cached dense W and applies it as a plain
# W @ b mat-vec.
# w_build_time is the wall-clock cost of building the kernel W (build_W_batched
# only), tracked separately from the per-rep estimation time.
Lgxg, W, w_build_time = simulate_Cholesky_gxg(SNP, s2gxg=s2gxg, s2e=s2e)

# Save the Cholesky factor (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_gxg_different_kernal/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}"
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Cache the pre-computed GRM W once (genotype-only; reused across variance
# settings and every replicate at estimation time).
#
# NOTE THE `W_ctr_` PREFIX.  MCREML_gxg caches the STANDARDIZED kernel as
# W_<mode>_n<n>_m<m>.npy -- keyed by genotype only, with no kernel in the name
# -- and both pipelines write only `if not os.path.exists`.  Sharing that name
# would mean whichever pipeline ran first wins and the other silently draws its
# phenotype from one kernel while estimating with the other.  Keep them apart.
# Old W_raw_*.npy caches belong to the previous RAW (un-centred) kernel and
# must not be reused for this one.
w_dir = "/home/ziyanzha/MOM_within_gene/stored_genotype"
os.makedirs(w_dir, exist_ok=True)
w_path = f"{w_dir}/W_ctr_{mode}_n{n}_m{m}.npy"
if not os.path.exists(w_path):
    np.save(w_path, W)
    print(f"W saved to: {w_path}")

# Record the kernel-build wall-clock time (seconds).  W is built once by this
# single Cholesky job, so there is nothing to average -- write it straight to
# time/result/W_timing_<FILENAME>.txt (parallels estimation's timing_<...>.txt).
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_gxg_different_kernal/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time: {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")
