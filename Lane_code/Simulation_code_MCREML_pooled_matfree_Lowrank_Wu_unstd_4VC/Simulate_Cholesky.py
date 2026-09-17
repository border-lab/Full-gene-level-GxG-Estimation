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
# UNSTANDARDIZED within-gene kernel, h_ab = Z_a .* Z_b with no centering and no
# 1/sigma_ab scaling -- and it is the EXACT kernel, not the rank-r truncation
# the estimator applies: the phenotype must be drawn from the full W so that any
# gap the estimator shows at small r is honestly the operator's truncation bias
# and nothing else.  K_a and K_d have no truncation on either side, so those two
# components are drawn and fitted with literally the same matrices.
#
# The dense GRMs are built here because a Cholesky factor needs an explicit
# matrix, but UNLIKE the _preW pipeline none is cached: estimation rebuilds the
# action of all three kernels matrix-free from the genotype, so no n-by-n array
# is ever written to disk or loaded downstream.  simulate_Cholesky_4vc builds
# and releases them one at a time to keep the peak near 4 n-by-n arrays -- this
# is the memory-critical job of the pipeline.
La, Ld, Lgxg, w_build_time = simulate_Cholesky_4vc(SNP, G, s2a=s2a, s2d=s2d,
                                                   s2gxg=s2gxg, s2e=s2e)

# Save the Cholesky factors (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
np.save(f"{save_dir}/La_{tag}.npy", La)
np.save(f"{save_dir}/Ld_{tag}.npy", Ld)
np.save(f"{save_dir}/Lgxg_{tag}.npy", Lgxg)

# Record the EPISTASIS kernel-build wall-clock time (seconds).  W is built once
# by this single Cholesky job, so there is nothing to average -- write it
# straight to time/result/W_timing_<FILENAME>.txt (parallels estimation's
# timing_<...>.txt).  NOTE this is SIMULATION-ONLY cost; the estimator never
# pays it.  The additive and dominance GRMs are not timed separately: they are
# a single gemm each and are dwarfed by the Hadamard-square build next to them.
w_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/time/result"
os.makedirs(w_time_dir, exist_ok=True)
fname = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
with open(f"{w_time_dir}/W_timing_{fname}.txt", 'w') as f:
    f.write(f"{w_build_time:.4f}\n")
print(f"W build time (simulation only): {w_build_time:.4f} s -> {w_time_dir}/W_timing_{fname}.txt")

# Realized-variance scale factor c (realized_variance.pdf): E[Var-hat(H gamma)]
# = c * s2gxg.  It depends on the genotype only, so it is computed ONCE here.
# Both versions are written so the two error sources the note names can be
# read off directly:  c_hat (the O(nm) HWE closed form the estimator uses)
# vs c_exact (mean sample variance of the P interaction columns, O(n m_g^2)
# per gene -- affordable in this job, which already pays O(n^2 m) for W).
# Their gap is the HWE + moment-sampling error; the scatter of the per-rep
# realized variances around c * s2gxg is the gamma-draw sampling error.
#
# This is an EPISTASIS-ONLY quantity.  The additive and dominance designs are
# column-standardized, so their analogues are c_a = c_d = 1 exactly -- nothing
# to estimate and nothing to write; s2a_hat and s2d_hat are already on the
# realized-variance scale.  They are recorded as constants anyway so that a
# reader of this file does not have to remember which components carry a
# correction and which do not.
c_hat = compute_c_pooled(SNP, G, method='hwe')
c_exact = compute_c_pooled(SNP, G, method='exact')
result_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/result"
os.makedirs(result_dir, exist_ok=True)
with open(f"{result_dir}/c_{fname}.txt", 'w') as f:
    f.write(f"c_hat_hwe {c_hat:.10f}\n")
    f.write(f"c_exact {c_exact:.10f}\n")
    f.write(f"expected_realized_variance_hwe {c_hat * s2gxg:.10f}\n")
    f.write(f"expected_realized_variance_exact {c_exact * s2gxg:.10f}\n")
    f.write(f"c_additive 1.0000000000\n")
    f.write(f"expected_realized_variance_additive {s2a:.10f}\n")
    f.write(f"c_dominance 1.0000000000\n")
    f.write(f"expected_realized_variance_dominance {s2d:.10f}\n")
print(f"c_hat (HWE) = {c_hat:.6f}, c_exact = {c_exact:.6f}; "
      f"E[V_ell] = c * s2gxg = {c_exact * s2gxg:.6f} -> {result_dir}/c_{fname}.txt")
