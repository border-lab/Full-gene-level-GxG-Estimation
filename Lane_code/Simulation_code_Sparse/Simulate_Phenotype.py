from Function_sparse import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--density', type=float, required=True)
parser.add_argument('--sparseVersion', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
density = args.density
sparseVersion = args.sparseVersion
rep = args.rep



# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Load Lgxg
lgxg_dir = "/home/ziyanzha/MOM_within_gene/Sparse_model/Cholesky_Lgxg"
lgxg_path = f"{lgxg_dir}/{sparseVersion}Lgxg_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_density{density}.npy"
Lgxg = np.load(lgxg_path)

# Load La
la_dir = "/home/ziyanzha/MOM_within_gene/Sparse_model/Cholesky_La"
la_path = f"{la_dir}/{sparseVersion}La_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_density{density}.npy"
La = np.load(la_path)

Z, y = simulate_remove_sampling_err(SNP, Lgxg, La, s2a, s2gxg, s2e)

# Save phenotype
output_dir = f"/home/ziyanzha/MOM_within_gene/Sparse_model/Phenotype/{sparseVersion}y_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_density{density}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)