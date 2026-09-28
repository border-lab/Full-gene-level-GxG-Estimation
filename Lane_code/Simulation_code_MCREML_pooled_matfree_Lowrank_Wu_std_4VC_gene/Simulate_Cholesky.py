from Function_MCREML import *
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)
# Percent of the genotype that is GENES: the first round(m * gene_pct / 100)
# SNPs are cut into G contiguous genes for the epistasis kernel; the rest are
# in no gene.  K_a and K_d use all m SNPs regardless.  e.g. --gene_pct 5.
parser.add_argument('--gene_pct', type=float, required=True)


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode
gene_pct = args.gene_pct


# Read genotype (m SNPs).  The additive kernel K_a = Z_a Z_a'/m and the
# dominance kernel K_d = Z_d Z_d'/m each use all m SNPs at once; the epistasis
# kernel takes only the FIRST gene_pct percent of them and splits those into G
# contiguous genes inside simulate_Cholesky_4vc (several Z, one per gene) for
# the POOLED kernel.  The RAW 0/1/2 dosages are
# passed through: the dominance coding needs the allele frequencies, which the
# standardized designs no longer carry.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factors La (La La' = s2a K_a), Ld (Ld Ld' = s2d K_d) and Lgxg
# (Lgxg Lgxg' = s2gxg W) for drawing the three genetic effects.  W is the
# UNSTANDARDIZED within-gene kernel -- h_ab = Z_a .* Z_b with no centering and
# no 1/sigma_ab scaling -- DIVIDED BY THE REALIZED-VARIANCE FACTOR c-hat:
#
#     W = W_raw / c-hat ,   c = (1/P) sum_{pairs} Var-hat(Z_a .* Z_b) .
#
# That single division is what this pipeline adds to its _unstd_4VC sibling.
# It makes E[Var-hat(g_gxg)] = (c / c-hat) s2gxg, so a run with s2a = s2d =
# s2gxg = 0.1, s2e = 0.7 realizes 0.1 / 0.1 / 0.1 / 0.7 -- the additive and
# dominance designs are column-standardized and hit their targets exactly, and
# the epistasis component now does too up to the plug-in's error.
#
# c-hat IS THE THIRD-MOMENT PLUG-IN (Function_MCREML.pooled_c, the note's O(nm)
# closed form with the skewness read off the genotype).  It is deliberately NOT
# the O(n sum_g m_g^2) exact c: keeping the divisor O(nm) keeps it computable
# from what a real analysis has.
#
# It is also the EXACT kernel, not the rank-r truncation the estimator applies:
# the phenotype must be drawn from the full W so that any gap the estimator
# shows at small r is honestly the operator's truncation bias and nothing else.
# The estimator divides by the SAME c-hat (setup_pooled -> pooled_c on the same
# genotype), so the normalization introduces no new mismatch and the
# plug-in's own error cancels between the two sides.  K_a and
# K_d have no truncation and no normalization on either side, so those two
# components are drawn and fitted with literally the same matrices.
#
# The dense GRMs are built here because a Cholesky factor needs an explicit
# matrix, but UNLIKE the _preW pipeline none is cached: estimation rebuilds the
# action of all three kernels matrix-free from the genotype, so no n-by-n array
# is ever written to disk or loaded downstream.  simulate_Cholesky_4vc builds
# and releases them one at a time to keep the peak near 4 n-by-n arrays -- this
# is the memory-critical job of the pipeline.
La, Ld, Lgxg, w_build_time, _ = simulate_Cholesky_4vc(
    SNP, G, s2a=s2a, s2d=s2d, s2gxg=s2gxg, s2e=s2e, gene_pct=gene_pct)

# Save the Cholesky factors (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}_gp{gene_pct:g}"
print(f"genic SNPs: first {n_gene_snps(m, gene_pct)} of m={m} "
      f"({gene_pct:g}%), in G={G} genes", flush=True)
np.save(f"{save_dir}/La_{tag}.npy", La)
np.save(f"{save_dir}/Ld_{tag}.npy", Ld)
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Record the EPISTASIS kernel-build wall-clock time (seconds).  W is built once
# by this single Cholesky job, so there is nothing to average -- write it
# straight to time/result/W_timing_<FILENAME>.txt (parallels estimation's
# timing_<...>.txt).  NOTE this is SIMULATION-ONLY cost; the estimator never
# pays it.  The additive and dominance GRMs are not timed separately: they are
# a single gemm each and are dwarfed by the O(n^2 m) Hadamard square.
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}_gp{gene_pct:g}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time (simulation only): {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")
