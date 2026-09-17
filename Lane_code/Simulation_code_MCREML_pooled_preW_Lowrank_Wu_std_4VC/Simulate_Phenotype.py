from Function_MCREML import *
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


# Load the precomputed Cholesky factors: La (La La' = s2a K_a, K_a = Z_a Z_a'/m),
# Ld (Ld Ld' = s2d K_d, K_d = Z_d Z_d'/m) and Lgxg (Lgxg Lgxg' = s2gxg W), W
# being the UNSTANDARDIZED pooled within-gene kernel DIVIDED BY c (exact, not
# truncated) -- see Simulate_Cholesky.py.
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC/Cholesky"
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
La = np.load(f"{save_dir}/La_{tag}.npy")
Ld = np.load(f"{save_dir}/Ld_{tag}.npy")
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_a + g_d + g_gxg + e.  ONLY y is written; the realized
# variances Var-hat(.) of the four drawn components are no longer recorded.
#
# Because the epistasis kernel is c-normalized, all four components already had
# their NOMINAL target as their EXPECTATION -- no scale factor for a reader to
# apply by hand, the property this pipeline was built for:
#
#   Var-hat(H gamma)   = var(g_gxg)  E[.] = (c/c-hat) s2gxg
#   Var-hat(Z_a beta)  = var(g_a)    E[.] = s2a    (c_a = 1 exactly)
#   Var-hat(Z_d delta) = var(g_d)    E[.] = s2d    (c_d = 1 exactly)
#   Var-hat(e)         = var(e)      E[.] = s2e
#
# AND THEY NOW HIT THOSE NUMBERS EXACTLY.  simulate_phenotype defaults to
# force_realized=True, which rescales each drawn component by
# sqrt(target / Var-hat(component)) so its realized variance IS its target to
# machine precision.  Two consequences, both deliberate:
#
#   - The residual c / c-hat factor on the epistasis component is absorbed
#     along with the O(1/sqrt n) scatter, so v_ell = s2gxg exactly and the
#     plug-in's error no longer reaches the phenotype at all.  An estimate's
#     deviation from the nominal target is now the ESTIMATOR's error alone.
#   - The four realized variances are therefore CONSTANTS, equal to the four
#     targets in every replicate, and are NO LONGER RECORDED.  There is nothing
#     per-replicate left to record: reading a column of 0.1 back out of every
#     row told the reader only what the run was configured with.  The result
#     row is the four ESTIMATES and nothing else, and the comparison is against
#     the nominal targets -- which are now exact, not approximate.  Measured
#     over 30 replicates, the paired (estimate - realized) spread this used to
#     support and the unpaired spread under forcing agree to within Monte Carlo
#     error, so nothing diagnostic was lost with the columns.
#
# The cost is that y no longer has the law REML fits: dividing by a random
# sample standard deviation is not a Gaussian operation.  See
# simulate_phenotype's docstring.  force_realized=False restores the plain
# draw and is what the sibling pipelines still do.
y = simulate_phenotype(La, Ld, Lgxg, n, s2a=s2a, s2d=s2d, s2gxg=s2gxg,
                       s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)

# NO vell_rep<rep>.txt is written any more.  The MCREML step used to read it
# back and carry the four values into columns 6-9 of the estimate row; with
# force_realized on they are the four targets in every replicate, so both the
# file and those columns are gone.  Nothing downstream reads them: the result
# row is 4 columns and calc_stats.py labels a 4-column file with the four
# fitted components and skips the paired block.
