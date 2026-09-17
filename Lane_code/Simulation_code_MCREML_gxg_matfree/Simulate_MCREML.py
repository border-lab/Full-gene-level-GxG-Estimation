from Function_MCREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=50)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load the genotype and column-standardize it.  MC AI-REML applies the
# epistasis GRM W = (1/p) H H' MATRIX-FREE straight from Z -- rebuilt one
# SNP-pair batch at a time inside every CG pass -- so NO n-by-n W is loaded or
# stored (only the n-by-m genotype).
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)

# Load phenotype (s2gxg_s2e order)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_gxg_matfree/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2gxg W + s2e I, W applied matrix-free from Z.
# seed=rep fixes the probe randomness per replicate (reproducible, yet varied).
# Time the estimation only (the MC_REML call) for the per-rep timing record.
t_start = time.perf_counter()
s2gxg_hat, s2e_hat, _ = MC_REML(Z, y, iters=iters, Nmc=nmc, seed=rep)
elapsed = time.perf_counter() - t_start

# Save result (s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_gxg_matfree/result/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_gxg_matfree/time/rep_times/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")
