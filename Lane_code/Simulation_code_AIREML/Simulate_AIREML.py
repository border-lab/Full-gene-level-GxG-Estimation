from Function_AIREML import *
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
parser.add_argument('--iters', type=int, default=12)
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
rep = args.rep
mode = args.mode

# Load Z (standardized genotype, same file the MoM pipeline uses)
Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

# Load raw genotype dosages -- needed to build the dominance GRM
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None).to_numpy()

# Load phenotype (s2a_s2d_s2gxg_s2e order)
y_path = f"/home/ziyanzha/MOM_within_gene/AIREML/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Four-component AI-REML: V = s2a*K + s2d*D + s2gxg*W + s2e*I.
# Fitting the dominance GRM D jointly keeps the additive / GxG estimates
# identifiable when SNPs are in LD (Hivert et al. 2021).
a, d, gxg, e, _ = AI_REML(Z, y, SNP=SNP, iters=iters)

# Save result (s2a_s2d_s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/AIREML/result/{mode}_n{n}m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({a},{d},{gxg},{e})\n")
