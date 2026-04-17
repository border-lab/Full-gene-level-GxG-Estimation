from Function_additive import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str,required=True)  


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
s2a = args.s2a
mode = args.mode
rep = args.rep


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

save_dir = "/home/ziyanzha/MOM_within_gene/Whole_model/Cholesky_Lgxg"
save_path = f"{save_dir}/Lgxg_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}.npy"
Lgxg = np.load(save_path)

save_dir_s2a = "/home/ziyanzha/MOM_within_gene/Whole_model/Cholesky_La"
save_path_s2a = f"{save_dir_s2a}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}.npy"
La = np.load(save_path_s2a)

Z, y = simulate_remove_sampling_err(SNP, Lgxg, La, s2a,s2gxg,s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/Whole_model/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}"

os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
