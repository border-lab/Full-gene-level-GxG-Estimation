from Function_MCREML import *
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
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W) for drawing the epistasis effect.
# The dense W is built once here purely to form Lgxg (a one-off per genotype);
# it is NOT cached, because the estimation step applies W matrix-free straight
# from the genotype.
Lgxg = simulate_Cholesky_gxg(SNP, s2gxg=s2gxg, s2e=s2e)

# Save the Cholesky factor (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_gxg_matfree/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}"
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)
