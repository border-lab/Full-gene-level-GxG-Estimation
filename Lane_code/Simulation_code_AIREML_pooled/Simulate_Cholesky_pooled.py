from Function_AIREML_pooled import *
import argparse
import os

# Pooled-kernel counterpart of Simulate_Cholesky.py.
# Builds Lgxg from the pooled within-gene AxA kernel over G subsampled
# contiguous regions; (G, region_size, seed) fix the regions and must match
# the values used later in Simulate_AIREML_pooled.py.

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--region_size', type=int, required=True)
parser.add_argument('--seed', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
s2a = args.s2a
mode = args.mode
G = args.G
region_size = args.region_size
seed = args.seed

tag = f"G{G}_rs{region_size}_seed{seed}"


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

Lgxg, La, W = simulate_Cholesky_pooled_withadd(SNP, G, region_size, seed,
                                               s2a=s2a, s2gxg=s2gxg, s2e=s2e)

# Save Lgxg
save_dir = "/home/ziyanzha/MOM_within_gene/AIREML_pooled/Cholesky_Lgxg"
os.makedirs(save_dir, exist_ok=True)
save_path = f"{save_dir}/Lgxg_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}.npy"
np.save(save_path, Lgxg)


# Save La
save_dir_La = "/home/ziyanzha/MOM_within_gene/AIREML_pooled/Cholesky_La"
os.makedirs(save_dir_La, exist_ok=True)
save_path = f"{save_dir_La}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}.npy"

np.save(save_path, La)


# Save the pooled kernel W once -- reused directly by Simulate_AIREML_pooled.py
# so build_pooled_W is never recomputed.
save_dir_W = "/home/ziyanzha/MOM_within_gene/AIREML_pooled/W"
os.makedirs(save_dir_W, exist_ok=True)
save_path = f"{save_dir_W}/W_{mode}_n{n}_m{m}_s2a{s2a}_s2gxg{s2gxg}_s2e{s2e}_{tag}.npy"
np.save(save_path, W)
