from Function_MCREML import *
import argparse
import pandas as pd
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
# It applies to the EPISTASIS kernel ONLY: K_a = Z_a Z_a'/m and K_d = Z_d Z_d'/m
# are applied exactly at every r, so r never biases those estimates directly (it
# can still move them through the coupling in the AI step).
# verify_lowrank.py carries an independent exact implementation to check it
# against, and that is the only place an alternative lives.
parser.add_argument('--r', type=int, default=R_DEFAULT)
args = parser.parse_args()

m = args.m
n = args.n
G = args.G
s2a = args.s2a
s2d = args.s2d
s2gxg = args.s2gxg
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode
r = args.r

# Load the GENOTYPE and build the two standardized designs.  MC AI-REML applies
# ALL THREE genetic kernels MATRIX-FREE straight from them: the additive GRM as
# K_a B = Z_a(Z_a'B)/m, the dominance GRM as K_d B = Z_d(Z_d'B)/m, and the
# pooled within-gene epistasis GRM from the G contiguous genes Z_a is split into
# inside mc_reml.  NO n-by-n GRM is loaded or stored (only the two n-by-m
# designs), which is the whole point of this variant: the _preW pipeline loads a
# cached dense W here instead.
#
# The dominance design is built from the RAW 0/1/2 dosages -- it needs the
# allele frequencies -- and is the SAME transformation the Cholesky job used, so
# K_d is identical on both sides.  The epistasis kernel is the UNSTANDARDIZED
# one (h_ab = Z_a .* Z_b, no centering, no 1/sigma_ab scaling), also exactly the
# kernel Simulate_Cholesky.py drew the phenotype from.  So any bias seen in the
# results here is the rank-r TRUNCATION bias, NOT a kernel mismatch -- and
# unlike the StochasticWu sibling's Monte-Carlo error it is deterministic: every
# replicate fits the same W-hat, and the gap closes only by raising r.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)
Zd = dominance_design(SNP)

# Load phenotype (s2a_s2d_s2gxg_s2e order)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# This replicate's REALIZED variances -- epistasis V_ell = Var-hat(H gamma),
# additive Var-hat(Z_a beta), dominance Var-hat(Z_d delta) and residual
# Var-hat(e) -- written by the Phenotype step next to the y they belong to, as
# labelled "<name> <value>" lines.  They are properties of the effect draws, not
# of the fit, so they are only carried through here -- to end up in the same row
# as the estimates they should be compared against.  Read by NAME, not by
# position, so the row cannot silently scramble if the Phenotype step's write
# order ever changes.
vell_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/Phenotype/y_{tag}/vell_rep{rep}.txt"
with open(vell_path) as f:
    vell = {k: float(v) for k, v in
            (line.split() for line in f if line.strip())}
v_ell = vell['vell_gxg']
v_a = vell['vell_a']
v_d = vell['vell_d']
v_e = vell['vell_e']

# Monte-Carlo AI-REML: V = s2a K_a + s2d K_d + s2gxg W + s2e I, every GRM
# applied matrix-free from the designs (W by the rank-r truncation, K_a and K_d
# exactly).  The genotype-only setup (per-gene SVD factors, pair total P) is
# built once inside mc_reml, outside the REML iteration.  seed=rep fixes the
# Hutchinson probe randomness per replicate (reproducible, yet varied).  There
# is NO operator seed: W-hat is a deterministic function of (genotype, r), so it
# is automatically the same matrix for every replicate -- the property the
# StochasticWu sibling had to engineer with frozen probes.
# Time the estimation only (the MC_REML call) for the per-rep timing record --
# this INCLUDES the per-gene SVD setup, which is part of the estimator's cost
# here (in the _preW pipeline the equivalent work sits in the Cholesky job).
# NOTE: timings are NOT comparable to the two- or three-component pipelines',
# which solve a smaller system with fewer kernel applies per iteration.
t_start = time.perf_counter()
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, _ = MC_REML(Z, Zd, y, G, iters=iters,
                                                  Nmc=nmc, seed=rep, r=r)
elapsed = time.perf_counter() - t_start

# Realized-variance correction (realized_variance.pdf).  REML on the raw H
# estimates the variance COMPONENT V_gamma = s2gxg; the variance of the genetic
# value across individuals is V_ell = Var-hat(H gamma) with E[V_ell] = c * V_gamma,
# c = (1/P) sum_{pairs} Var-hat(Z_a .* Z_b).  c is a genotype-only constant,
# computed here by the note's O(nm) HWE closed form (allele frequencies +
# skewness, one mat-vec per gene) -- negligible next to the REML fit, and
# outside the timed block so timing_<...>.txt stays comparable to earlier runs.
# The exact c and c-hat for this genotype are written once by the Cholesky job
# (result/c_<FILENAME>.txt).  s2a_hat and s2d_hat get NO such correction: their
# designs are standardized, so their factors are 1 exactly.
c_hat = compute_c_pooled(SNP, G, method='hwe')
Vl_hat = c_hat * s2gxg_hat

# Save result, 9 columns:
#
#   (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, Vl_hat, v_ell, v_a, v_d, v_e)
#    V_a      V_d      V_gamma    V_e      V_l     realized: gxg, add, dom, res
#
# The first four are the FIT (the variance components REML estimated), the
# fifth is the c-corrected ESTIMATE of the realized epistasis variance, and the
# last four are what this replicate's draws actually REALIZED.  So each estimate
# has its draw-level counterpart in the same row:
#
#   V_l  <-> v_ell   (both on the realized-variance scale; c applied to V_gamma)
#   V_a  <-> v_a     (no correction: c_a = 1 exactly)
#   V_d  <-> v_d     (no correction: c_d = 1 exactly)
#   V_e  <-> v_e     (no correction)
#
# V_gamma stays on the raw-H component scale and is the column to compare with
# the nominal S2GXG.  Note v_ell + v_a + v_d + v_e != var(y): the sample
# cross-terms between the four independent draws are O(1/sqrt n), not 0.
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/result/{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2gxg_hat},{s2e_hat},{Vl_hat},"
            f"{v_ell},{v_a},{v_d},{v_e})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC/time/rep_times/{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")
