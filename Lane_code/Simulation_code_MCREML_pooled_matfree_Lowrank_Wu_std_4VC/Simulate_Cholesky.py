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


args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
mode = args.mode


# Read genotype (m SNPs).  The additive kernel K_a = Z_a Z_a'/m and the
# dominance kernel K_d = Z_d Z_d'/m each use all m SNPs at once; the epistasis
# kernel splits them into G contiguous genes inside simulate_Cholesky_4vc
# (several Z, one per gene) for the POOLED kernel.  The RAW 0/1/2 dosages are
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
# c-hat IS THE THIRD-MOMENT PLUG-IN (Function_MCREML.C_METHOD = 'moment', the
# note's O(nm) closed form with the skewness read off the genotype rather than
# predicted from the allele frequency under HWE).  It is deliberately NOT the
# O(n sum_g m_g^2) exact c: keeping the divisor O(nm) keeps it computable from
# what a real analysis has, and the residual scale factor c / c-hat is not
# hidden but written out below next to all three routes.
#
# It is also the EXACT kernel, not the rank-r truncation the estimator applies:
# the phenotype must be drawn from the full W so that any gap the estimator
# shows at small r is honestly the operator's truncation bias and nothing else.
# The estimator divides by the SAME c-hat (setup_pooled -> pooled_c on the same
# genotype with the same C_METHOD), so the normalization introduces no new
# mismatch and the plug-in's own error cancels between the two sides.  K_a and
# K_d have no truncation and no normalization on either side, so those two
# components are drawn and fitted with literally the same matrices.
#
# The dense GRMs are built here because a Cholesky factor needs an explicit
# matrix, but UNLIKE the _preW pipeline none is cached: estimation rebuilds the
# action of all three kernels matrix-free from the genotype, so no n-by-n array
# is ever written to disk or loaded downstream.  simulate_Cholesky_4vc builds
# and releases them one at a time to keep the peak near 4 n-by-n arrays -- this
# is the memory-critical job of the pipeline.
La, Ld, Lgxg, w_build_time, c_norm = simulate_Cholesky_4vc(
    SNP, G, s2a=s2a, s2d=s2d, s2gxg=s2gxg, s2e=s2e)

# Save the Cholesky factors (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
np.save(f"{save_dir}/La_{tag}.npy", La)
np.save(f"{save_dir}/Ld_{tag}.npy", Ld)
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Record the EPISTASIS kernel-build wall-clock time (seconds).  W is built once
# by this single Cholesky job, so there is nothing to average -- write it
# straight to time/result/W_timing_<FILENAME>.txt (parallels estimation's
# timing_<...>.txt).  NOTE this is SIMULATION-ONLY cost; the estimator never
# pays it.  It now includes the exact-c computation, one gemm per gene, which
# is negligible against the O(n^2 m) Hadamard square.  The additive and
# dominance GRMs are not timed separately: they are a single gemm each and are
# dwarfed by that square.
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time (simulation only): {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")

# The normalization constant, written out next to the run -- ALL THREE ROUTES,
# so the one that was used can be scored rather than trusted.
#
# c_norm_applied is the divisor that WAS applied to the kernel, taken straight
# from simulate_Cholesky_4vc so this file cannot drift from what was built.  It
# is the THIRD-MOMENT plug-in (C_METHOD = 'moment'), the O(nm) route: the
# note's closed form Var(Z_a Z_b) = 1 + r_ab s_a s_b with s-hat_a =
# mean_t(Z_ta^3) estimated from the standardized genotype, assuming nothing
# about Hardy-Weinberg.
#
# c_hat_hwe is the same closed form with s predicted from the allele frequency
# under HWE -- the note's original estimator.  Its gap against c_hat_moment is
# the HWE departure alone, the two routes differing in nothing else.
#
# c_exact is the literal mean per-pair sample variance, O(n sum_g m_g^2): the
# YARDSTICK, affordable here because this job already pays O(n^2 m), but not
# used as the divisor.  Both plug-ins are scored against it, and the ratio
#
#     c_gxg_after_normalization = c_exact / c_norm_applied
#
# is the realized-variance factor the NORMALIZED kernel still carries -- 1 if
# the plug-in were perfect, and the number the epistasis target is multiplied
# by otherwise.  It is the honest expected_realized_variance_gxg, and it is
# also what verify_lowrank.py measures directly as tr(P_c W)/n.
#
# The additive and dominance designs are column-standardized, so their factors
# are 1 exactly with nothing to compute and nothing to estimate; they are
# written as constants so a reader does not have to remember which components
# carry a correction.
c_hat_moment = compute_c_pooled(SNP, G, method='moment')
c_hat_hwe = compute_c_pooled(SNP, G, method='hwe')
c_exact = compute_c_pooled(SNP, G, method='exact')
c_gxg = c_exact / c_norm            # what the normalized kernel still carries
result_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/result"
os.makedirs(result_dir, exist_ok=True)
with open(f"{result_dir}/c_{fname}.txt", 'w') as f:
    f.write(f"c_method moment\n")
    f.write(f"c_norm_applied {c_norm:.10f}\n")
    f.write(f"c_hat_moment {c_hat_moment:.10f}\n")
    f.write(f"c_hat_hwe {c_hat_hwe:.10f}\n")
    f.write(f"c_exact {c_exact:.10f}\n")
    f.write(f"moment_rel_error {abs(c_hat_moment - c_exact) / c_exact:.10f}\n")
    f.write(f"hwe_rel_error {abs(c_hat_hwe - c_exact) / c_exact:.10f}\n")
    f.write(f"c_gxg_after_normalization {c_gxg:.10f}\n")
    f.write(f"expected_realized_variance_gxg {c_gxg * s2gxg:.10f}\n")
    f.write(f"c_additive 1.0000000000\n")
    f.write(f"expected_realized_variance_additive {s2a:.10f}\n")
    f.write(f"c_dominance 1.0000000000\n")
    f.write(f"expected_realized_variance_dominance {s2d:.10f}\n")
    f.write(f"expected_realized_variance_residual {s2e:.10f}\n")
# c_norm must BE the moment plug-in -- the same deterministic function of the
# same genotype the kernel called.  If it is not, the file above would describe
# a kernel that was never built, so say so loudly rather than write it quietly.
assert c_norm == c_hat_moment, (
    f"the kernel divided by {c_norm!r} but the moment plug-in gives "
    f"{c_hat_moment!r}: C_METHOD and build_W_pooled disagree.")
print(f"c applied to W (third moment) = {c_norm:.6f}; HWE = {c_hat_hwe:.6f}; "
      f"exact = {c_exact:.6f}; after normalization E[V_ell] = "
      f"{c_gxg:.6f} * s2gxg = {c_gxg * s2gxg:.6f} (target {s2gxg:.6f}) "
      f"-> {result_dir}/c_{fname}.txt")
