from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2e = args.s2e
mode = args.mode


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

La = simulate_Cholesky_from_std_withadd(SNP, s2a=s2a, s2e=s2e)

# Save La (additive Cholesky factor)
save_dir_La = "/home/ziyanzha/MOM_within_gene/MCREML_additive/Cholesky_La"
os.makedirs(save_dir_La, exist_ok=True)
save_path = f"{save_dir_La}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2e{s2e}.npy"
np.save(save_path, La)
