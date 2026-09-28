from Function_MCREML import *
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2e = args.s2e
mode = args.mode


# Read genotype (m SNPs).  The additive kernel K_a = Z_a Z_a'/m uses all m SNPs
# at once -- it is the ONLY genetic kernel in this pipeline: no dominance, no
# epistasis, so no gene split and no G.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factor La (La La' = s2a K_a) for drawing the additive effect.  The
# design is column-standardized, so c_a = 1 exactly and a run with s2a = 0.1,
# s2e = 0.9 realizes 0.1 / 0.9 in expectation with no scale factor to correct
# -- and, with force_realized on in the Phenotype step, exactly.
#
# The dense GRM is built here because a Cholesky factor needs an explicit
# matrix, but it is not cached: estimation applies K_a matrix-free from the
# genotype, so no n-by-n array is ever loaded downstream except La itself.
La, k_build_time = simulate_Cholesky_additive(SNP, s2a=s2a, s2e=s2e)

# Save the Cholesky factor (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2a{s2a}_s2e{s2e}_n{n}_m{m}"
np.save(f"{save_dir}/La_{tag}.npy", La)

# Record the GRM-build-and-factor wall-clock time (seconds): the gemm and the
# Cholesky.  Built once by this single job, so there is nothing to average --
# write it straight to time/result/K_timing_<FILENAME>.txt (parallels
# estimation's timing_<...>.txt, and the 4VC parent's W_timing).
# SIMULATION-ONLY cost; the estimator never pays it.
k_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/time/result"
os.makedirs(k_time_dir, exist_ok=True)
fname = f"{mode}_s2a{s2a}_s2e{s2e}_n{n}m{m}"
with open(f"{k_time_dir}/K_timing_{fname}.txt", 'w') as f:
    f.write(f"{k_build_time:.4f}\n")
print(f"K_a build + Cholesky time (simulation only): {k_build_time:.4f} s "
      f"-> {k_time_dir}/K_timing_{fname}.txt")
