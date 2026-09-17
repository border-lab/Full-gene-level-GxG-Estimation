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
# UNSTANDARDIZED pooled within-gene kernel (exact, not truncated).
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/Cholesky"
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
Lgxg = np.load(f"{save_dir}/Lgxg_{tag}.npy")

# Phenotype: y = g_gxg + e, plus the REALIZED variance of this replicate's
# genetic value, V_ell = Var-hat(H gamma) = var(g_gxg) (ddof=0).  g_gxg =
# Lgxg u1 has exactly the law of H gamma, gamma ~ N(0, s2gxg/P I), so this is
# the note's V_ell for one gamma draw -- random, with E[V_ell] = c * s2gxg.
y, v_ell = simulate_phenotype(Lgxg, n, s2gxg=s2gxg, s2e=s2e, return_realized=True)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/Phenotype/y_{tag}"
os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)

# Record V_ell NEXT TO the phenotype it belongs to, one file per replicate.
# The MCREML step reads it back and carries it into the estimate row as the
# fourth column (V_gamma, V_e, V_l, realized_variance), so every replicate's
# realized variance is kept alongside the fit it should be compared with --
# no separate result/realized_variance/ tree, and no reduction to mean/std.
# Cleaned up with the rest of the phenotype directory by the combine job.
with open(f"{output_dir}/vell_rep{rep}.txt", 'w') as f:
    f.write(f"{v_ell}\n")
