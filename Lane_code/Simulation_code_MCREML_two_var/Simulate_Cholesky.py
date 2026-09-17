from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2e = args.s2e
mode = args.mode


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Additive + dominance Cholesky factors: La La' = s2a K_a, Ld Ld' = s2d K_d.
La, Ld = simulate_Cholesky_two_var(SNP, s2a=s2a, s2d=s2d, s2e=s2e)

# Save La, Ld (additive / dominance Cholesky factors)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_two_var/Cholesky"
os.makedirs(save_dir, exist_ok=True)
np.save(f"{save_dir}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2e{s2e}.npy", La)
np.save(f"{save_dir}/Ld_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2e{s2e}.npy", Ld)
