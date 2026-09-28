from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
parser.add_argument('--rep', type=int, required=True)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
rep = args.rep


# Load the precomputed Cholesky factor Lgxg (Lgxg Lgxg' = s2gxg W), W being the
# UNSTANDARDIZED pooled within-gene kernel DIVIDED BY c (exact, not truncated)
# -- see Simulate_Cholesky.py.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/Cholesky"
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_gxg + e.  ONLY y is written.
#
# Because the epistasis kernel is c-normalized, both components already had
# their NOMINAL target as their EXPECTATION:
#
#   Var-hat(H gamma)   = var(g_gxg)  E[.] = (c/c-hat) s2gxg
#   Var-hat(e)         = var(e)      E[.] = s2e
#
# AND THEY HIT THOSE NUMBERS EXACTLY: simulate_phenotype defaults to
# force_realized=True (as in the 4VC parent), which rescales each drawn
# component by sqrt(target / Var-hat(component)) so its realized variance IS
# its target to machine precision.  The residual c / c-hat factor is absorbed
# along with the O(1/sqrt n) scatter, so an estimate's deviation from the
# nominal target is the ESTIMATOR's error alone, the two realized variances
# are constants and are not recorded, and the result row is the two ESTIMATES
# and nothing else.
#
# The cost is that y no longer has the law REML fits: dividing by a random
# sample standard deviation is not a Gaussian operation.  See
# simulate_phenotype's docstring.  force_realized=False restores the plain
# draw.
y = simulate_phenotype(Lgxg, n, s2gxg=s2gxg, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)
