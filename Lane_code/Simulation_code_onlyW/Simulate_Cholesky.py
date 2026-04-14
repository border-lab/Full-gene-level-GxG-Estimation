from Function_onlyM_std import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)  
parser.add_argument('--mode', type=str, required=True) 


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

Lgxg = simulate_Cholesky_from_std(SNP , s2gxg, s2e)

# Save Lgxg
save_dir = "/home/ziyanzha/MOM_within_gene/Only_W_model/Cholesky_Lgxg"
save_path = f"{save_dir}/Lgxg_{mode}_n{n}_m{m}_s2gxg{s2gxg}_s2e{s2e}.npy"
np.save(save_path, Lgxg)

