from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2e = args.s2e
mode = args.mode
rep = args.rep


# Load the precomputed Cholesky factors: La (La La' = s2a K_a, K_a = Z_a Z_a'/m)
# and Ld (Ld Ld' = s2d K_d, K_d = Z_d Z_d'/m) -- see Simulate_Cholesky.py.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/Cholesky"
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}_m{m}"
La = np.load(f"{save_dir}/La_{tag}.npy")
Ld = np.load(f"{save_dir}/Ld_{tag}.npy")

# Phenotype: y = g_a + g_d + e.  ONLY y is written.
#
# Both designs are column-standardized, so both genetic components already had
# their NOMINAL target as their EXPECTATION:
#
#   Var-hat(Z_a beta)  = var(g_a)    E[.] = s2a
#   Var-hat(Z_d delta) = var(g_d)    E[.] = s2d
#   Var-hat(e)         = var(e)      E[.] = s2e
#
# AND THEY HIT THOSE NUMBERS EXACTLY: simulate_phenotype defaults to
# force_realized=True (as in the 4VC parent), which rescales each drawn
# component by sqrt(target / Var-hat(component)) so its realized variance IS
# its target to machine precision.  An estimate's deviation from the nominal
# target is therefore the ESTIMATOR's error alone, the three realized
# variances are constants and are not recorded, and the result row is the
# three ESTIMATES and nothing else.
#
# The cost is that y no longer has the law REML fits: dividing by a random
# sample standard deviation is not a Gaussian operation.  See
# simulate_phenotype's docstring.  force_realized=False restores the plain
# draw.
y = simulate_phenotype(La, Ld, n, s2a=s2a, s2d=s2d, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
