from Function_MCREML import *
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=50)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load Za, Zd (additive / dominance designs).  W is applied MATRIX-FREE inside
# MC_REML via the (S, R, T) weight matrices rebuilt from Za -- never stored.
Za = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Za_{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Zd = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Zd_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

# Load phenotype (s2a_s2d_s2gxg_s2e order)
tag = f"{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_three_var/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2a K_a + s2d K_d + s2gxg W + s2e I (W matrix-free).
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, _ = MC_REML(Za, Zd, y, iters=iters, Nmc=nmc, seed=rep)

# Save result (s2a_s2d_s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_three_var/result/{mode}_n{n}m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2gxg_hat},{s2e_hat})\n")
