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
parser.add_argument('--nmc', type=int, default=100)   # FINE phase; coarse is 15
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
# Print the AI-REML trace (one line per iteration: the four components and
# max|step|) to STDOUT.  Off by default -- with it on, Step 3 must also be
# launched with a real --output, since the pipeline sends stdout to /dev/null.
parser.add_argument('--verbose', action='store_true')
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
# K_d is identical on both sides.  The epistasis kernel is the UNSTANDARDIZED,
# C-NORMALIZED one (h_ab = Z_a .* Z_b, no centering and no 1/sigma_ab scaling,
# the whole kernel then divided by c-hat = pooled_c, i.e. the O(nm)
# third-moment plug-in C_METHOD selects), also exactly the kernel
# Simulate_Cholesky.py drew the phenotype from -- the divisor is recomputed
# here by the same deterministic function of the same genotype, so it is the
# same number, bit for bit, and the plug-in's approximation error is common to
# the two sides rather than a mismatch between them.  So any bias seen in the results here is the rank-r TRUNCATION
# bias, NOT a kernel mismatch -- and unlike the StochasticWu sibling's
# Monte-Carlo error it is deterministic: every replicate fits the same W-hat,
# and the gap closes only by raising r.  The normalization does not touch that
# gap: it scales W-hat and W by the same constant.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)
Zd = dominance_design(SNP)

# Load phenotype (s2a_s2d_s2gxg_s2e order)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# NOTHING IS READ BACK FROM THE PHENOTYPE STEP EXCEPT y.  The realized
# variances Var-hat(.) of the four effect draws used to be read from
# Phenotype/y_<TAG>/vell_rep<rep>.txt and carried into columns 6-9 of the row.
# simulate_phenotype now rescales each component to hit its target exactly
# (force_realized=True), so all four are the nominal targets in every replicate
# and there is nothing per-replicate left to carry: the Phenotype step no longer
# writes the file and this step no longer looks for it.

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
if args.verbose:
    print(f"--- rep{rep}: AI-REML trace, var(y)={y.var():.6f}, "
          f"columns s2a s2d s2gxg s2e ---", flush=True)
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, _ = MC_REML(Z, Zd, y, G, iters=iters,
                                                  Nmc=nmc, seed=rep, r=r,
                                                  verbose=args.verbose)
elapsed = time.perf_counter() - t_start
# How many LINEAR OPERATOR APPLIES that fit actually cost.  mc_reml zeroes the
# counters on entry, so this snapshot is THIS replicate and nothing else -- the
# setup applies (K_a U, K_d U, W U for the fixed probes) included, since they
# are part of what estimating this replicate took.  Written below next to the
# wall-clock time and averaged over replicates by the combine job.
op_counts = get_op_counts()

# NO post-fit realized-variance correction -- this is where this pipeline
# differs from its _unstd_4VC sibling.  There, REML on the raw H estimated the
# component V_gamma = s2gxg and the realized variance had to be recovered as
# V_l = c-hat * s2gxg-hat, so the plug-in's error landed on the ESTIMATE and
# only on the estimate.  Here the same c-hat -- the O(nm) third-moment plug-in,
# C_METHOD = 'moment' -- is already INSIDE the kernel, both when the phenotype
# was drawn and when it is fitted, so
#
#     E[Var-hat(g_gxg)] = (c / c-hat) s2gxg   and   V_l = s2gxg-hat ,
#
# with the post-fit correction factor equal to 1 by construction.  Applying
# compute_c_pooled here as well would double-count the normalization and
# inflate the epistasis estimate by a factor of c-hat.
#
# WHERE THE PLUG-IN'S ERROR NOW GOES: NOWHERE.  The residual factor c / c-hat
# used to rescale the SIMULATED epistasis component and the FITTED kernel by the
# same constant, shifting the target a run aims at (not the accuracy of aiming)
# and leaving column 3's target at (c / c-hat) * s2gxg rather than s2gxg.  With
# force_realized=True the simulated component is rescaled to hit s2gxg EXACTLY,
# so that factor is absorbed at the draw and column 3's target is the nominal
# number with nothing to correct for.  The Cholesky job still writes the ratio
# out per run as c_gxg_after_normalization (result/c_<FILENAME>.txt) as a
# property of the genotype, but no column depends on it any more.
#
# NO Vl_hat COLUMN.  It was identically s2gxg_hat in this pipeline -- kept only
# so the row matched the sibling layouts -- and a column that repeats another
# column is a column a reader has to be warned about.  A file from here no
# longer has the same width as a sibling's; read it with THIS directory's
# calc_stats.py.  s2a_hat, s2d_hat and s2e_hat needed no correction either: the
# additive and dominance designs are standardized, so their factors were 1 to
# begin with -- exactly, and with nothing estimated.

# Save result, 4 columns -- the FIT and nothing else:
#
#   (s2a_hat, s2d_hat, s2gxg_hat, s2e_hat)
#    V_a      V_d      V_gamma    V_e
#
# ALL FOUR ARE ON THE REALIZED-VARIANCE SCALE and all four are directly
# comparable to the nominal targets: with S2A = S2D = S2GXG = 0.1 and S2E = 0.7
# the four column means should sit at 0.1 / 0.1 / 0.1 / 0.7, with no rescaling
# and no per-replicate reference to pair against.  That is what forcing the
# realized variances bought -- the targets are exact, so the nominal comparison
# IS the paired comparison.
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/result/{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
run_tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/time/rep_times/{run_tag}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")

# Record this replicate's OPERATOR-APPLY COUNTS, one file per replicate, as
# labelled "<name> <value>" lines, read back by name so the row cannot scramble
# if this list ever grows.  The combine
# step averages all reps into time/result/opcount_<FILENAME>.txt.
#
# WHY NEXT TO THE TIME.  The wall-clock number above is what this machine took;
# these are what the ALGORITHM did, and they are deterministic -- rerun the
# replicate anywhere and they come out identical, because CG's iteration count
# is fixed by (V, rhs, tol).  So a timing difference between two pipelines, or
# two r, or two Nmc, splits cleanly: the counts say how much work was asked
# for, the seconds say how fast the machine did it.  V_columns is the one to
# quote as "the" cost -- every operator here is linear in the column count --
# and V_applies next to it says how well that work was batched (a solve with
# Nmc right-hand sides is ONE apply and Nmc columns).
#
# The keys are written in a FIXED order with 0 defaults, so every rep file has
# the same lines whether or not a counter was ever bumped.
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "K_applies", "K_columns",
           "W_applies", "W_columns")
op_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/time/op_counts/{run_tag}"
os.makedirs(op_dir, exist_ok=True)
with open(f"{op_dir}/rep{rep}.txt", 'w') as f:
    for key in OP_KEYS:
        f.write(f"{key} {op_counts.get(key, 0)}\n")
print(f"operator applies this replicate: V={op_counts.get('V_applies', 0)} "
      f"({op_counts.get('V_columns', 0)} columns), "
      f"K={op_counts.get('K_applies', 0)} "
      f"({op_counts.get('K_columns', 0)} columns), "
      f"W={op_counts.get('W_applies', 0)} "
      f"({op_counts.get('W_columns', 0)} columns); "
      f"{op_counts.get('reml_iters', 0)} REML iters, "
      f"{op_counts.get('cg_iters', 0)} CG iters in "
      f"{op_counts.get('cg_solves', 0)} solves -> {op_dir}/rep{rep}.txt")
