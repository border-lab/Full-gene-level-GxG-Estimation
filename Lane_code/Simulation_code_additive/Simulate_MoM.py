from Function_additive import *
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--nmc', type=int, default=100)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
s2a = args.s2a
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load Z
Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

# Load phenotype (s2a_s2gxg_s2e order)
y_path = f"/home/ziyanzha/MOM_within_gene/Whole_model/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

a, gxg, e, _ = MoM_std(Z, y, nmc=nmc)

# Save result (s2a_s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/Whole_model/result/{mode}_n{n}m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({a},{gxg},{e})\n")