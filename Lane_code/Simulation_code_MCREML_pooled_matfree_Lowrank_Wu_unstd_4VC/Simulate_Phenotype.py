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
# being the UNSTANDARDIZED pooled within-gene kernel (exact, not truncated).
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/Cholesky"
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
La = np.load(f"{save_dir}/La_{tag}.npy")
Ld = np.load(f"{save_dir}/Ld_{tag}.npy")
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_a + g_d + g_gxg + e, plus the REALIZED variance of ALL FOUR
# drawn components (ddof=0):
#
#   v_ell = Var-hat(H gamma)   = var(g_gxg)  E[.] = c * s2gxg  (c != 1)
#   v_a   = Var-hat(Z_a beta)  = var(g_a)    E[.] = s2a        (c_a = 1 exactly)
#   v_d   = Var-hat(Z_d delta) = var(g_d)    E[.] = s2d        (c_d = 1 exactly)
#   v_e   = Var-hat(e)         = var(e)      E[.] = s2e
#
# g_gxg = Lgxg u has exactly the law of H gamma, gamma ~ N(0, s2gxg/P I), so
# v_ell is the note's V_ell for one gamma draw.  The additive, dominance and
# residual realized variances need no c correction, but they ARE still random:
# one draw scatters around its target by O(1/sqrt n).  Recording them lets each
# estimate be compared against what this replicate actually realized, not only
# against the nominal target.
y, v_ell, v_a, v_d, v_e = simulate_phenotype(La, Ld, Lgxg, n, s2a=s2a, s2d=s2d,
                                             s2gxg=s2gxg, s2e=s2e,
                                             return_realized=True)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)

# Record the realized variances NEXT TO the phenotype they belong to, one file
# per replicate, as LABELLED lines (the two-component pipeline wrote a bare
# number; four values need names so the reader cannot silently mis-order them).
# The MCREML step reads them back and carries them into the estimate row as its
# last four columns, so every replicate's realized variances are kept alongside
# the fit they should be compared with -- no separate result/realized_variance/
# tree, and no reduction to mean/std.  Cleaned up with the rest of the phenotype
# directory by the combine job.
with open(f"{output_dir}/vell_rep{rep}.txt", 'w') as f:
    f.write(f"vell_gxg {v_ell}\n")
    f.write(f"vell_a {v_a}\n")
    f.write(f"vell_d {v_d}\n")
    f.write(f"vell_e {v_e}\n")
