from Function_AIREML_pooled import *
import argparse
import pandas as pd
import os

# Pooled-kernel counterpart of Simulate_AIREML.py.
# Fits V = s2a*K + s2gxg*W + s2e*I, where W is the pooled kernel built ONCE and
# stored by Simulate_Cholesky_pooled.py (same W that generated the phenotype).

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--iters', type=int, default=12)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--region_size', type=int, required=True)
parser.add_argument('--seed', type=int, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
s2a = args.s2a
iters = args.iters
rep = args.rep
mode = args.mode
G = args.G
region_size = args.region_size
seed = args.seed

tag = f"G{G}_rs{region_size}_seed{seed}"

# Load Z (standardized genotype, same file the single-gene pipeline uses)
Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

# Load phenotype (pooled, tagged by G/region_size/seed)
y_path = f"/home/ziyanzha/MOM_within_gene/AIREML_pooled/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Load the pooled kernel W built once by Simulate_Cholesky_pooled.py
W_path = f"/home/ziyanzha/MOM_within_gene/AIREML_pooled/W/W_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}.npy"
W = np.load(W_path)

a, gxg, e, _ = AI_REML_pooled(Z, y, W, iters=iters)

# Save result (s2a_s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/AIREML_pooled/result/{mode}_n{n}m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({a},{gxg},{e})\n")
