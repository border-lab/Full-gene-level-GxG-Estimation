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



args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
density = args.density
sparseVersion = args.sparseVersion 

# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()
if sparseVersion == "SNP":
    La, Lgxg = simulate_Cholesky_SNPsparse(SNP, density, s2a, s2gxg, s2e, stability=1e-10)
elif sparseVersion == "pair":
    La = simulate_Cholesky_additivesparse(SNP, density, s2a, stability=1e-10)
    Lgxg = simulate_Cholesky_pairsparse(SNP, density, s2gxg, stability=1e-10)

# Save Lgxg
save_dir = "/home/ziyanzha/MOM_within_gene/Sparse_model/Cholesky_Lgxg"
save_path = f"{save_dir}/{sparseVersion}Lgxg_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_density{density}.npy"
np.save(save_path, Lgxg)


# Save La
save_dir = "/home/ziyanzha/MOM_within_gene/Sparse_model/Cholesky_La"
save_path = f"{save_dir}/{sparseVersion}La_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_density{density}.npy"
np.save(save_path, La)


