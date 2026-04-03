from function_whole_model import *
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--nmc', type=int, default=40)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
nmc = args.nmc
rep = args.rep
mode = args.mode


Z = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/Stored_data/genotype/Z_{mode}_n{n}_m{m}.csv", header=None).to_numpy()
y = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/Stored_data/phenotype/STD_y_{mode}_n{n}_m{m}/rep{rep}.csv", header=None).to_numpy().flatten()

gxg, e, _ = MoM_only_M_std(Z, y, nmc=nmc)

# Save result
filename = f"/home/ziyanzha/MOM_within_gene/result/result_328/{mode}_STD_n{n}_m{m}.txt"
with open(filename, 'a') as f:
    f.write(f"({gxg.item()},{e.item()})\n")
