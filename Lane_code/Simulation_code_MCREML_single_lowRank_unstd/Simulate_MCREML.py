from Function_MCREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--G', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=50)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score traces tr(V^{-1} K_i) are Hutchinson's, with Nmc CG solves per
# REML iteration -- the ONLY estimator now; the SLQ alternative is gone, so
# there is no trace_method to pass.
# Truncation level of the low-rank W u apply -- the ONE knob this variant adds,
# and the only accuracy control there is.  Per gene, the top-r eigenpairs of
# K_g = Z_g Z_g' (from the thin SVD of Z_g) replace one copy of the Hadamard
# square, so the apply is O(n m_g r) and DETERMINISTIC: the error is a
# truncation BIAS set by the discarded eigenvalue tail, not sampling noise --
# small under LD, large under linkage equilibrium (Note_0814, Tables 1-3) --
# and only raising r reduces it.  r >= min(n, m_g) makes the apply EXACT.
# verify_lowrank.py carries an independent exact implementation to check it
# against, and that is the only place an alternative lives.
parser.add_argument('--r', type=int, default=R_DEFAULT)
args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode
r = args.r

# Load the GENOTYPE and column-standardize it.  MC AI-REML applies the pooled
# within-gene epistasis GRM W MATRIX-FREE straight from Z -- the m SNPs are
# split into G contiguous genes inside mc_reml.  NO n-by-n W is loaded or stored
# (only the n-by-m genotype), which is the whole point of this variant: the
# _preW pipeline loads a cached dense W here instead.
#
# The kernel is the UNSTANDARDIZED one (h_ab = Z_a .* Z_b, no centering, no
# 1/sigma_ab scaling), which is exactly the kernel Simulate_Cholesky.py drew the
# phenotype from.  So any bias seen in the results here is the rank-r
# TRUNCATION bias, NOT a kernel mismatch -- and unlike the StochasticWu
# sibling's Monte-Carlo error it is deterministic: every replicate fits the
# same W-hat, and the gap closes only by raising r.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)

# Load phenotype (s2gxg_s2e order)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# This replicate's REALIZED variance V_ell = Var-hat(H gamma), written by the
# Phenotype step next to the y it belongs to.  It is a property of the gamma
# draw, not of the fit, so it is only carried through here -- to end up in the
# same row as the estimates it should be compared against.
vell_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/Phenotype/y_{tag}/vell_rep{rep}.txt"
with open(vell_path) as f:
    v_ell = float(f.read().strip())

# Monte-Carlo AI-REML: V = s2gxg W + s2e I, pooled W applied matrix-free from Z
# by the rank-r truncation.  The genotype-only setup (per-gene SVD factors,
# pair total P) is built once inside mc_reml, outside the REML iteration.
# seed=rep fixes the Hutchinson probe randomness per replicate (reproducible,
# yet varied).  There is NO operator seed: W-hat is a deterministic function of
# (genotype, r), so it is automatically the same matrix for every replicate --
# the property the StochasticWu sibling had to engineer with frozen probes.
# Time the estimation only (the MC_REML call) for the per-rep timing record --
# this INCLUDES the per-gene SVD setup, which is part of the estimator's cost
# here (in the _preW pipeline the equivalent work sits in the Cholesky job).
# NOTE: timings are NOT comparable to runs made on the SLQ path, which skipped
# the Nmc probe solves in most iterations.
t_start = time.perf_counter()
s2gxg_hat, s2e_hat, _ = MC_REML(Z, y, G, iters=iters, Nmc=nmc, seed=rep, r=r)
elapsed = time.perf_counter() - t_start

# Realized-variance correction (realized_variance.pdf).  REML on the raw H
# estimates the variance COMPONENT V_gamma = s2gxg; the variance of the genetic
# value across individuals is V_ell = Var-hat(H gamma) with E[V_ell] = c * V_gamma,
# c = (1/P) sum_{pairs} Var-hat(Z_a .* Z_b).  c is a genotype-only constant,
# computed here by the note's O(nm) HWE closed form (allele frequencies +
# skewness, one mat-vec per gene) -- negligible next to the REML fit, and
# outside the timed block so timing_<...>.txt stays comparable to earlier runs.
# The exact c and c-hat for this genotype are written once by the Cholesky job
# (result/c_<FILENAME>.txt).
c_hat = compute_c_pooled(SNP, G, method='hwe')
Vl_hat = c_hat * s2gxg_hat

# Save result: (s2gxg_hat, s2e_hat, Vl_hat, v_ell) -- V_gamma, V_e, V_l and the
# realized variance of THIS replicate.  The third column is the c-corrected
# ESTIMATE of the realized variance and the fourth is the realized variance
# actually drawn, so the two are paired replicate by replicate; calc_stats.py
# summarises both columns over the reps.
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/result/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2gxg_hat},{s2e_hat},{Vl_hat},{v_ell})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/time/rep_times/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")
