from Function_sparse import *
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)  
parser.add_argument('--mode', type=str, required=True) 
parser.add_argument('--density', type=float, required=True)
parser.add_argument('--sparseVersion', type=str, required=True) 
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--nmc', type=int, default=100)


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
density = args.density
sparseVersion = args.sparseVersion 
rep= args.rep
nmc = args.nmc

# Load Z
Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()

y_path = f"/home/ziyanzha/MOM_within_gene/Sparse_model/Phenotype/{sparseVersion}y_{mode}_n{n}_m{m}_s2gxg{s2gxg}_s2e{s2e}_density{density}/rep{rep}.csv"

y = pd.read_csv(y_path, header=None).to_numpy().flatten()

gxg, e, _ = MoM_only_M_std(Z, y, nmc=nmc)


output_dir = f"/home/ziyanzha/MOM_within_gene/Sparse_model/result/{mode}_n{n}m{m}_s2gxg{s2gxg}_s2e{s2e}_{sparseVersion}density{density}"

os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"

with open(filename, 'w') as f:
    f.write(f"({gxg},{e})\n")