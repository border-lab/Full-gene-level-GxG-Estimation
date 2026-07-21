from Function_MCREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=50)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# Per-gene SNP shares; must match the Cholesky step so the cached W is found.
parser.add_argument('--ratio', type=str, default=None)
# Genes used to build the ESTIMATION kernel; must match the Cholesky step so the
# cached W_est is found.  The phenotype was simulated from all G genes.
parser.add_argument('--estimate', type=str, default=None)
args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode
parse_ratio(args.ratio, G)              # validate: same spec as the Cholesky step
parse_gene_subset(args.estimate, G)     # validate: same spec as the Cholesky step
split = gene_split_tag(G, args.ratio)
esttag = gene_subset_tag(args.estimate, G)      # "" or "_est1-2"

# Load the PRE-COMPUTED pooled epistasis GRM (cached once by the Cholesky step)
# and apply it as a dense W @ b mat-vec inside MC AI-REML -- no genotype needed
# here.  With --estimate this is W_est, pooled over the SELECTED genes only,
# while the phenotype below was simulated from the full G-gene kernel: the fit
# is deliberately misspecified.  esttag is "" for the all-genes fit, so that
# case loads the same full W as before.
W = np.load(f"/home/ziyanzha/MOM_within_gene/stored_genotype/W_{mode}_n{n}_m{m}_{split}{esttag}.npy")

# Load phenotype (s2gxg_s2e order).  Its tag carries NO est suffix -- the
# phenotype depends only on the full kernel, so one set is shared across every
# estimation subset at the same split.
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_{split}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2gxg W + s2e I (dense pre-computed pooled W).
# seed=rep fixes the probe randomness per replicate (reproducible, yet varied).
# Time the estimation only (the MC_REML call) for the per-rep timing record.
t_start = time.perf_counter()
s2gxg_hat, s2e_hat, _ = MC_REML(W, y, iters=iters, Nmc=nmc, seed=rep)
elapsed = time.perf_counter() - t_start

# Save result (s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/result/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_{split}{esttag}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/time/rep_times/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_{split}{esttag}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")
