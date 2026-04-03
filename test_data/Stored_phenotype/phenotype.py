from function_whole_model import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--nmc', type=int, default=40)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str)  # Fixed: str not string
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
nmc = args.nmc
rep = args.rep
mode = args.mode

# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/Stored_data/genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()


Lgxg, Le = simulate_Cholesky_from_std(SNP, s2gxg, s2e, stability=1e-8)
Z, y = simulate_from_raw_only_W_remove_sampling_err(SNP, Lgxg, Le, s2gxg, s2e)

# Write Z once (only on first rep)
z_path = f"/home/ziyanzha/MOM_within_gene/Stored_data/genotype/Z_{mode}_n{n}_m{m}.csv"
if rep == 1 and not os.path.exists(z_path):
    pd.DataFrame(Z).to_csv(z_path, index=False, header=False)
    print(f"Saved Z: {z_path}")

output_dir = f"/home/ziyanzha/MOM_within_gene/Stored_data/phenotype/STD_y_{mode}_n{n}_m{m}"
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
