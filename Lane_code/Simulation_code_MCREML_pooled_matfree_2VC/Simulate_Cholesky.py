from Function_MCREML import *
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2e = args.s2e
mode = args.mode


# Read genotype (m SNPs).  The additive kernel K_a = Z_a Z_a'/m and the
# dominance kernel K_d = Z_d Z_d'/m each use all m SNPs at once; there is no
# gene split because there is no epistasis kernel in this pipeline.  The RAW
# 0/1/2 dosages are passed through: the dominance coding needs the allele
# frequencies, which the standardized designs no longer carry.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Cholesky factors La (La La' = s2a K_a) and Ld (Ld Ld' = s2d K_d) for drawing
# the two genetic effects.  Both designs are column-standardized, so a run with
# s2a = s2d = 0.1, s2e = 0.7 realizes 0.1 / 0.1 / 0.7 in expectation with no
# scale factor to correct -- and, with force_realized on in the Phenotype
# step, exactly.
#
# The dense GRMs are built here because a Cholesky factor needs an explicit
# matrix, but none is cached: estimation rebuilds the action of both kernels
# matrix-free from the genotype, so no n-by-n array is ever written to disk or
# loaded downstream.  simulate_Cholesky_2vc builds and releases them one at a
# time; this is still the memory-critical job of the pipeline.
La, Ld, k_build_time = simulate_Cholesky_2vc(SNP, s2a=s2a, s2d=s2d, s2e=s2e)

# Save the Cholesky factors (per variance target)
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/Cholesky"
os.makedirs(save_dir, exist_ok=True)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}_m{m}"
np.save(f"{save_dir}/La_{tag}.npy", La)
np.save(f"{save_dir}/Ld_{tag}.npy", Ld)

# Record the GRM-build-and-factor wall-clock time (seconds): both gemms and
# both Choleskys.  Built once by this single job, so there is nothing to
# average -- write it straight to time/result/K_timing_<FILENAME>.txt
# (parallels estimation's timing_<...>.txt, and the 4VC parent's W_timing).
# SIMULATION-ONLY cost; the estimator never pays it.
k_time_dir = "/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/time/result"
os.makedirs(k_time_dir, exist_ok=True)
fname = f"{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}m{m}"
with open(f"{k_time_dir}/K_timing_{fname}.txt", 'w') as f:
    f.write(f"{k_build_time:.4f}\n")
print(f"K_a/K_d build + Cholesky time (simulation only): {k_build_time:.4f} s "
      f"-> {k_time_dir}/K_timing_{fname}.txt")
