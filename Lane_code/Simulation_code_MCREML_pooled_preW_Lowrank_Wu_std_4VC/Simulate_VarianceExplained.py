from Function_MCREML import *
import argparse
import os

# Diagnostic for the ONE knob this pipeline has: the truncation level r.
#
# The low-rank operator replaces one copy of the Hadamard square (K .* K) by
# its rank-r part, so the truncation bias is governed by how much of K_g =
# Z_g Z_g' the leading r eigenvalues carry:
#
#     rho_g(r) = sum_{s<=r} lam_s / sum_s lam_s ,   lam_s = eig(K_g) = sv(Z_g)^2
#
# and sum_s lam_s = tr(Z_g' Z_g) = ||Z_g||_F^2 needs no factorization.  This is
# exactly Function_MCREML.variance_explained.  It depends on the GENOTYPE ONLY
# -- not on s2a, s2gxg, s2e, or any phenotype draw -- so this job takes neither
# and runs once per (mode, n, m, G).  rho_g(r) ~ 1 means r is already big enough
# for that gene (strong LD); rho_g(r) growing ~ r/min(n,m_g) means the spectrum
# is flat (linkage equilibrium) and no affordable r will be enough.
#
# The ADDITIVE kernel is untouched by r -- K_a = Z Z'/m is applied exactly --
# so rho_all below is reported for context (it is the spectrum of the whole Z,
# which is K_a's spectrum up to the 1/m scale), NOT as an accuracy budget.
#
# The c-NORMALIZATION this pipeline applies to the epistasis kernel does NOT
# enter here at all: rho is a RATIO of eigenvalues of the same matrix, and
# dividing W by a positive constant divides numerator and denominator alike.
# So this table is identical to the _unstd_4VC sibling's for the same genotype,
# and the r it recommends is the same r.

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# Truncation levels to report.  Clipped to each matrix's full rank and
# de-duplicated, so the same list can be used for the whole genotype
# (rank <= min(n, m)) and for a single gene (rank <= min(n, m_g)).
parser.add_argument('--r_list', type=str,
                    default='1,2,5,10,20,30,40,50,60,80,100,150,200')
args = parser.parse_args()

n, m, G, mode = args.n, args.m, args.G, args.mode
r_grid = sorted({int(x) for x in args.r_list.split(',') if x.strip()})

# Same genotype and same additive design the estimator uses: the truncation is
# applied to the COLUMN-STANDARDIZED Z_g (only the interaction weights are
# unstandardized in this variant), so rho must be measured on that Z_g too.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv",
                  header=None).to_numpy()
Z = additive_design(SNP)
genes = split_into_genes(Z, G)               # the SAME contiguous blocks

out_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC/result/variance_explained"
os.makedirs(out_dir, exist_ok=True)
fname = f"{mode}_n{n}m{m}_G{G}"

rank_all = min(Z.shape)
rank_g = [min(Zg.shape) for Zg in genes]

lines = []
lines.append(f"# variance explained by the leading r eigenvalues of Z Z'")
lines.append(f"# mode={mode} n={n} m={m} G={G}")
lines.append(f"# full rank: whole Z = {rank_all}; genes = {rank_g}")
lines.append(f"# rho(r) = sum_(s<=r) lam_s / ||Z||_F^2   (1.0 => r is full rank)")
lines.append(f"# rho_all is the spectrum of the ADDITIVE kernel K_a (context only:")
lines.append(f"#   K_a is applied exactly, r truncates the per-gene rho_g below)")
lines.append("# r\trho_all\t" + "\t".join(f"rho_g{g+1}" for g in range(len(genes)))
             + "\trho_gene_mean\trho_gene_min")

for r in r_grid:
    rho_all = variance_explained(Z, r)
    rho_gs = []
    for Zg in genes:
        if Zg.shape[1] < 2:                  # a 1-SNP gene contributes no pair
            rho_gs.append(float('nan'))
            continue
        rho_gs.append(variance_explained(Zg, r))
    finite = [v for v in rho_gs if v == v]
    row = (f"{r}\t{rho_all:.6f}\t" + "\t".join(f"{v:.6f}" for v in rho_gs)
           + f"\t{np.mean(finite):.6f}\t{np.min(finite):.6f}")
    lines.append(row)
    print(f"r={r:4d}  rho_all={rho_all:.6f}  "
          f"gene mean={np.mean(finite):.6f}  gene min={np.min(finite):.6f}", flush=True)

with open(f"{out_dir}/{fname}.txt", 'w') as f:
    f.write("\n".join(lines) + "\n")
print(f"-> {out_dir}/{fname}.txt")
