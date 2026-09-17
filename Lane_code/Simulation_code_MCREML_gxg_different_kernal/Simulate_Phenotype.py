from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
rep = args.rep


# Load the precomputed Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W), W being the
# MEAN-CENTRED-kernel pairwise-epistasis GRM.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_gxg_different_kernal/Cholesky"
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}"
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_gxg + e with exact-variance rescaling.
y = simulate_remove_sampling_err(Lgxg, n, s2gxg=s2gxg, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_gxg_different_kernal/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
