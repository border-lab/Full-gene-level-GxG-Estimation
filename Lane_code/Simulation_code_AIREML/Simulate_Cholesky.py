from Function_AIREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

Lgxg, La, Ld = simulate_Cholesky_from_std_withadd(SNP, s2a=s2a, s2d=s2d, s2gxg=s2gxg, s2e=s2e)

# Save Lgxg
save_dir = "/home/ziyanzha/MOM_within_gene/AIREML/Cholesky_Lgxg"
os.makedirs(save_dir, exist_ok=True)
save_path = f"{save_dir}/Lgxg_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}.npy"
np.save(save_path, Lgxg)


# Save La
save_dir_La = "/home/ziyanzha/MOM_within_gene/AIREML/Cholesky_La"
os.makedirs(save_dir_La, exist_ok=True)
save_path = f"{save_dir_La}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}.npy"
np.save(save_path, La)


# Save Ld (dominance)
save_dir_Ld = "/home/ziyanzha/MOM_within_gene/AIREML/Cholesky_Ld"
os.makedirs(save_dir_Ld, exist_ok=True)
save_path = f"{save_dir_Ld}/Ld_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}.npy"
np.save(save_path, Ld)
