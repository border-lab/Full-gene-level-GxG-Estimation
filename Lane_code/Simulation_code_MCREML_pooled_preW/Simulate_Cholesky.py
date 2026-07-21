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
# 1-based gene indices used AT ESTIMATION ONLY, e.g. --estimate 1,2 fits the
# REML model with a kernel built from the front two genes while the phenotype
# is still simulated from all G.  Omit (or "all") for the correctly specified
# all-genes fit.
parser.add_argument('--estimate', type=str, default=None)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
ratio = parse_ratio(args.ratio, G)          # normalized shares (None = even)
split = gene_split_tag(G, args.ratio)       # cache key, from the raw string
subset = parse_gene_subset(args.estimate, G)    # 0-based genes to fit with
esttag = gene_subset_tag(args.estimate, G)      # "" or "_est1-2"


# Read genotype (m SNPs), then split into G contiguous genes inside
# simulate_Cholesky_gxg (several Z, one per gene) for the POOLED kernel.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W_full) for drawing the epistasis
# effect -- the phenotype is ALWAYS simulated from the full G-gene kernel --
# plus the GRMs themselves.  With --estimate, W_est is a second kernel pooled
# over the selected genes only; that is what REML is fitted with, so the fit is
# deliberately misspecified while the truth stays the full panel.  Both kernels
# are deterministic per (genotype, split, subset), so they are PRE-COMPUTED and
# cached once under stored_genotype and merely loaded at estimation time.
Lgxg, W, W_est, info = simulate_Cholesky_gxg(SNP, G, s2gxg=s2gxg, s2e=s2e,
                                             ratio=ratio, subset=subset)
w_build_time = info["w_build_time"]
print(f"gene sizes: {info['gene_sizes']}")
print(f"estimation genes: {info['subset'] or 'all'}  "
      f"P_est/P_full = {info['pair_share']:.6f}")

# Save the Cholesky factor (per variance target).  It depends only on the FULL
# kernel, so its tag carries no est suffix and one Cholesky/phenotype set is
# shared by every --estimate subset at the same split.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_{split}"
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Cache the pre-computed pooled GRMs once (genotype+split+subset only; reused
# across variance settings and every replicate at estimation time).  The split
# tag (G plus the ratio) keeps different gene splits from colliding, and the
# est suffix keeps different estimation subsets from colliding.
w_dir = "/home/ziyanzha/MOM_within_gene/stored_genotype"
os.makedirs(w_dir, exist_ok=True)
w_path = f"{w_dir}/W_{mode}_n{n}_m{m}_{split}.npy"
if not os.path.exists(w_path):
    np.save(w_path, W)
    print(f"W saved to: {w_path}")

# The estimation kernel (only when a strict subset was requested; with all genes
# the estimation step just loads the full W above).
if W_est is not None:
    w_est_path = f"{w_dir}/W_{mode}_n{n}_m{m}_{split}{esttag}.npy"
    if not os.path.exists(w_est_path):
        np.save(w_est_path, W_est)
        print(f"W_est saved to: {w_est_path}")

# Record the gene split and the subset's share of within-gene pairs.  Because
# the truth is s2gxg * W_full but the fit spans only the subset's genes, the
# estimate is attenuated towards s2gxg * P_est/P_full -- keep that number with
# the run so the results can be read on either scale.
info_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/info"
os.makedirs(info_dir, exist_ok=True)
info_name = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_{split}{esttag}"
with open(f"{info_dir}/split_{info_name}.txt", 'w') as f:
    f.write(f"gene_sizes\t{info['gene_sizes']}\n")
    f.write(f"estimate_genes\t{info['subset'] or 'all'}\n")
    f.write(f"P_full\t{info['P_full']}\n")
    f.write(f"P_est\t{info['P_est']}\n")
    f.write(f"pair_share\t{info['pair_share']:.10f}\n")
    f.write(f"s2gxg_true\t{s2gxg}\n")
    f.write(f"s2gxg_attenuated_target\t{s2gxg * info['pair_share']:.10f}\n")

# Record the kernel-build wall-clock time (seconds).  W is built once by this
# single Cholesky job, so there is nothing to average -- write it straight to
# time/result/W_timing_<FILENAME>.txt (parallels estimation's timing_<...>.txt).
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_{split}{esttag}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time: {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")
