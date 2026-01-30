import argparse
import random
import pandas as pd
from function_whole_model import *

parser = argparse.ArgumentParser(description="A script with arguments")
parser.add_argument('--s2gxg', type=float)
parser.add_argument('--s2e', type=float)
parser.add_argument('--n', type=int)
parser.add_argument('--nmc', type=int)
parser.add_argument('--iter', type=int)

args = parser.parse_args()
s2gxg = args.s2gxg
n = args.n
s2e= args.s2e
nmc= args.nmc
iter = args.iter


raw_file = pd.read_csv("/home/ziyanzha/MOM_within_gene/genotype_matrix_plink/chr1_block1000_pruned_09.raw", sep=r"\s+")
raw_data = raw_file.iloc[:, 6:].to_numpy(dtype=float)
subset = raw_data[:n]


Z, y = simulate_from_raw_only_W(subset, s2gxg, s2e)
m = Z.shape[1]
s2gxg, s2e ,_ = MoM_only_M(Z, y, nmc)
result = (
    float(s2gxg.item()),
    float(s2e.item()),
)


with open(f"/home/ziyanzha/MOM_within_gene/result/raw_method_{n}m{m}.txt", "a") as f:
     f.write(str(result) + '\n')
