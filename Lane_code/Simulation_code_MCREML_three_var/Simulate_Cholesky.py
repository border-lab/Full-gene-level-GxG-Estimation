from Function_MCREML import *
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
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factors: La La'=s2a K_a, Ld Ld'=s2d K_d, Lgxg Lgxg'=s2gxg W.
# (W is built once here only to form Lgxg; estimation never sees it -- it
# applies W matrix-free via compute_WU.)
La, Ld, Lgxg, _ = simulate_Cholesky_three_var(SNP, s2a=s2a, s2d=s2d, s2gxg=s2gxg, s2e=s2e)

# Save the three Cholesky factors
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_three_var/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}"
np.save(f"{save_dir}/La_{tag}.npy", La)
np.save(f"{save_dir}/Ld_{tag}.npy", Ld)
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)
