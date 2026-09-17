from Function_AIREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
rep = args.rep
DIR = "/home/ziyanzha/MOM_within_gene/AIREML_pooled_preW_exact_Wu_std_4VC"


# Load the precomputed Cholesky factors: La (La La' = s2a K_a, K_a = Z_a Z_a'/m),
# Ld (Ld Ld' = s2d K_d, K_d = Z_d Z_d'/m) and Lgxg (Lgxg Lgxg' = s2gxg W), W
# being the UNSTANDARDIZED pooled within-gene kernel DIVIDED BY c-hat -- the
# same matrix the Cholesky job cached for AI-REML.  See Simulate_Cholesky.py.
save_dir = f"{DIR}/Cholesky"
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
La = np.load(f"{save_dir}/La_{tag}.npy")
Ld = np.load(f"{save_dir}/Ld_{tag}.npy")
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_a + g_d + g_gxg + e.  ONLY y is written.
#
# simulate_phenotype defaults to force_realized=True, which rescales each drawn
# component by sqrt(target / Var-hat(component)) so its realized variance IS
# its target to machine precision.  That absorbs the O(1/sqrt n) scatter of
# each realized variance and the residual c / c-hat factor on the epistasis
# component, so an estimate's deviation from the nominal target is the
# ESTIMATOR's error alone.  The cost is that y no longer has exactly the
# Gaussian law REML fits -- see simulate_phenotype's docstring.
y = simulate_phenotype(La, Ld, Lgxg, n, s2a=s2a, s2d=s2d, s2gxg=s2gxg,
                       s2e=s2e)

output_dir = f"{DIR}/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
