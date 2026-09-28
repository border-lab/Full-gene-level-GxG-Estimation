from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2e = args.s2e
mode = args.mode
rep = args.rep


# Load the precomputed Cholesky factor La (La La' = s2a K_a, K_a = Z_a Z_a'/m)
# -- see Simulate_Cholesky.py.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/Cholesky"
tag = f"{mode}_s2a{s2a}_s2e{s2e}_n{n}_m{m}"
La = np.load(f"{save_dir}/La_{tag}.npy")

# Phenotype: y = g_a + e.  ONLY y is written.  As in the 4VC parent,
# simulate_phenotype defaults to force_realized=True, which rescales each
# drawn component by sqrt(target / Var-hat(component)) so its realized
# variance IS its target to machine precision -- an estimate's deviation from
# the nominal target is the ESTIMATOR's error alone, at the price of y no
# longer having exactly the Gaussian law REML fits.  See simulate_phenotype's
# docstring.
y = simulate_phenotype(La, n, s2a=s2a, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
