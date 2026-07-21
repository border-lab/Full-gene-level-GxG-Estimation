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
# Per-gene share of the m SNPs, e.g. --ratio 0.1,0.2,0.3,0.2,0.2 with G=5 cuts
# the genotype into 5 CONTIGUOUS blocks taking 10%, 20%, 30%, 20%, 20% of the
# SNPs in order.  Omit (or "even") for the equal split.
parser.add_argument('--ratio', type=str, default=None)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
ratio = parse_ratio(args.ratio, G)          # normalized shares (None = even)
split = gene_split_tag(G, args.ratio)       # cache key, from the raw string


# Read genotype (m SNPs), then split into G contiguous genes inside
# simulate_Cholesky_gxg (several Z, one per gene) for the POOLED kernel.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) for drawing the epistasis effect,
# plus the POOLED within-gene epistasis GRM W = (1/P) sum_g H_g H_g' itself.
# W is deterministic per (genotype, G), so it is PRE-COMPUTED and cached once
# under stored_genotype; the estimation step loads that cached dense W and
# applies it as a plain W @ b mat-vec.  w_build_time is the wall-clock cost of
# building the kernel W (build_W_pooled only), tracked separately from the
# per-rep estimation time.
Lgxg, W, w_build_time = simulate_Cholesky_gxg(SNP, G, s2gxg=s2gxg, s2e=s2e,
                                              ratio=ratio)

# Save the Cholesky factor (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_{split}"
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Cache the pre-computed pooled GRM W once (genotype+split only; reused across
# variance settings and every replicate at estimation time).  The split tag
# (G plus the ratio) keeps different gene splits from colliding.
w_dir = "/home/ziyanzha/MOM_within_gene/stored_genotype"
os.makedirs(w_dir, exist_ok=True)
w_path = f"{w_dir}/W_{mode}_n{n}_m{m}_{split}.npy"
if not os.path.exists(w_path):
    np.save(w_path, W)
    print(f"W saved to: {w_path}")

# Record the kernel-build wall-clock time (seconds).  W is built once by this
# single Cholesky job, so there is nothing to average -- write it straight to
# time/result/W_timing_<FILENAME>.txt (parallels estimation's timing_<...>.txt).
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_{split}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time: {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")
