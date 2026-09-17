from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True, help='number of gene blocks')
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)
# Per-gene SNP shares; must match the Cholesky step so the tags line up.
parser.add_argument('--ratio', type=str, default=None)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
rep = args.rep
parse_ratio(args.ratio, G)              # validate: same spec as the Cholesky step
rtag = ratio_tag(args.ratio)


# Load the precomputed Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W_full).  Its tag
# carries NO est suffix -- the phenotype depends only on the full kernel, so one
# set is shared across every estimation subset at the same split.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_perG/Cholesky"
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_G{G}_n{n}_m{m}{rtag}"
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_gxg + e with exact-variance rescaling.
y = simulate_remove_sampling_err(Lgxg, n, s2gxg=s2gxg, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_perG/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
