from Function_MCREML import *
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=50)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load Z (standardized genotype, same file the MoM / AI-REML pipelines use)
Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

# Load phenotype (s2a_s2e order)
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_additive/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2e{s2e}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2a*K + s2e*I, K = ZZ'/m built implicitly.
# seed=rep fixes the probe randomness per replicate (reproducible, yet varied
# across reps).
s2a_hat, s2e_hat, _ = MC_REML(Z, y, iters=iters, nmc=nmc, seed=rep)

# Save result (s2a_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_additive/result/{mode}_n{n}m{m}_s2a{s2a}_s2e{s2e}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2e_hat})\n")
