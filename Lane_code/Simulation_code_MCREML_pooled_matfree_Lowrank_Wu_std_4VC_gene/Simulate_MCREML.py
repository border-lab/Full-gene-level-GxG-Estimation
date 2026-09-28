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
# Percent of the genotype that is GENES: the first round(m * gene_pct / 100)
# SNPs are cut into G contiguous genes for the epistasis kernel; the rest are
# in no gene.  K_a and K_d use all m SNPs regardless.  e.g. --gene_pct 5.
parser.add_argument('--gene_pct', type=float, required=True)
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
parser.add_argument('--r', type=int, default=R_DEFAULT)
# The W-hat apply picks its route per apply by column count (Function_MCREML.
# W_MIN): stored-A -- At = [sqrt(lam_s) q_s .* Z_g] built once, two gemm pairs
# -- for wide applies, broadcast from the per-gene SVD factors for narrow ones.
# At is ~r times Z in size and always built (its size is printed); --A_dtype
# float32 halves it.  Memory is bounded by the Slurm --mem request.
parser.add_argument('--A_dtype', choices=('float64', 'float32'),
                    default='float64')
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
gene_pct = args.gene_pct
r = args.r

# Load the GENOTYPE and build the two standardized designs.  MC AI-REML applies
# ALL THREE genetic kernels MATRIX-FREE straight from them: the additive GRM as
# K_a B = Z_a(Z_a'B)/m, the dominance GRM as K_d B = Z_d(Z_d'B)/m, and the
# pooled within-gene epistasis GRM from the G contiguous genes the FIRST
# gene_pct percent of Z_a's columns is split into inside mc_reml (the remaining
# SNPs are in no gene; K_a and K_d still use all m).  NO n-by-n GRM is loaded or stored (only the two n-by-m
# designs), which is the whole point of this variant: the _preW pipeline loads a
# cached dense W here instead.
#
# The dominance design is built from the RAW 0/1/2 dosages -- it needs the
# allele frequencies -- and is the SAME transformation the Cholesky job used, so
# K_d is identical on both sides.  The epistasis kernel is the UNSTANDARDIZED,
# C-NORMALIZED one (h_ab = Z_a .* Z_b, no centering and no 1/sigma_ab scaling,
# the whole kernel then divided by c-hat = pooled_c, the O(nm)
# third-moment plug-in), also exactly the kernel
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
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}_gp{gene_pct:g}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/Phenotype/y_{tag}/rep{rep}.csv"
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
s2a_hat, s2d_hat, s2gxg_hat, s2e_hat, _, info = MC_REML(
    Z, Zd, y, G, iters=iters, Nmc=nmc, seed=rep, r=r, verbose=args.verbose,
    A_dtype=args.A_dtype, gene_pct=gene_pct)
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
# only on the estimate.  Here the same c-hat -- the O(nm) third-moment plug-in
# pooled_c -- is already INSIDE the kernel, both when the phenotype
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
# number with nothing to correct for; no column depends on the ratio.
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
#
# ALL FOUR ARE ON THE REALIZED-VARIANCE SCALE and directly comparable to the
# nominal targets: with S2A = S2D = S2GXG = 0.1 and S2E = 0.7 the four column
# means should sit at 0.1 / 0.1 / 0.1 / 0.7, with no rescaling and no
# per-replicate reference to pair against.  The optimizer's diagnostics
# (converged, rejected steps, components at their bound) go to this job's
# stdout log only.
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/result/{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}_gp{gene_pct:g}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
run_tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}_gp{gene_pct:g}"
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/time/rep_times/{run_tag}"
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
#
# WALL-CLOCK TWINS (keys ending "_sec", seconds on THIS machine).  Same file,
# same averaging, so a run's summary carries both what was asked for and how
# long each kind of work took here.  They NEST -- V_sec contains the K and W
# time spent inside V applies, W_sec contains the W_* pieces, each
# solve_*_sec contains its V applies -- so they are not summed; the
# *_sec_per_col lines (seconds per n-vector) are the comparable unit for the
# three operators.  The W_* pieces split W's time by route (each apply picks
# one by its column count; W_storedA_columns / W_bcast_columns count the
# columns each handled):
#   storedA  W_gemm_A = At(At'U), W_gemm_D = D(D'U); setup_A_sec = the At build
#   bcast    W_bcast = the q_s .* U broadcast, W_gemm1 = Zg'(.), W_gemm2 =
#            Zg(.), W_einsumD = the lam contraction plus the D correction
# Solve groups:
# solve_probe_coarse / _fine (the merged solve V^-1 [y, U], split by phase),
# solve_ai (V^-1 [K_i x]); setup is svd + At + K_i U; phase_coarse /
# phase_fine stamp the switch; cg_negcurv counts CG solves stopped because V
# was not positive definite at the point being evaluated.
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "K_applies", "K_columns",
           "W_applies", "W_columns",
           "W_storedA_columns", "W_bcast_columns",
           "cg_negcurv",
           "V_sec", "K_sec", "W_sec",
           "V_sec_per_col", "K_sec_per_col", "W_sec_per_col",
           "W_gemm_A_sec", "W_gemm_D_sec",
           "W_bcast_sec", "W_gemm1_sec", "W_gemm2_sec", "W_einsumD_sec",
           "solve_probe_coarse_sec", "solve_probe_fine_sec",
           "solve_ai_sec",
           "setup_sec", "setup_svd_sec", "setup_KU_sec",
           "setup_A_sec",
           "phase_coarse_sec", "phase_fine_sec")
# Seconds per column for the three operators: the machine-dependent unit cost
# of one n-vector through each apply.  Nested as above: V includes its K and W.
for _op in ("V", "K", "W"):
    _cols = op_counts.get(f"{_op}_columns", 0)
    op_counts[f"{_op}_sec_per_col"] = (
        op_counts.get(f"{_op}_sec", 0.0) / _cols if _cols else 0.0)
op_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC_gene/time/op_counts/{run_tag}"
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
_g = lambda k: op_counts.get(k, 0.0)
_wsec = _g('W_sec')
_wpct = lambda k: 100.0 * _g(k) / _wsec if _wsec else 0.0
print(f"timing this replicate ({elapsed:.2f} s total): "
      f"setup {_g('setup_sec'):.2f} s "
      f"(svd {_g('setup_svd_sec'):.2f}, At {_g('setup_A_sec'):.2f}, "
      f"K_iU/WU {_g('setup_KU_sec'):.2f}); "
      f"phase coarse {_g('phase_coarse_sec'):.2f} s, fine {_g('phase_fine_sec'):.2f} s; "
      f"solves: [y,U] coarse "
      f"{_g('solve_probe_coarse_sec'):.2f}, [y,U] fine "
      f"{_g('solve_probe_fine_sec'):.2f}, AI {_g('solve_ai_sec'):.2f}; "
      f"cg_negcurv {op_counts.get('cg_negcurv', 0)}")
print(f"optimizer: converged={info['converged']}, iters={info['n_iters']}, "
      f"rejected={info['n_reject']}, switch at iter {info['switch_it']}, "
      f"at_bound={info['at_bound']}")
print(f"seconds per column: V {_g('V_sec_per_col'):.3e}, "
      f"K {_g('K_sec_per_col'):.3e}, W {_g('W_sec_per_col'):.3e}; "
      f"W columns by route: stored-A {op_counts.get('W_storedA_columns', 0)}, "
      f"broadcast {op_counts.get('W_bcast_columns', 0)}; "
      f"W time split: At gemm {_wpct('W_gemm_A_sec'):.1f}%, "
      f"D gemm {_wpct('W_gemm_D_sec'):.1f}%, "
      f"broadcast {_wpct('W_bcast_sec'):.1f}%, "
      f"gemm1 {_wpct('W_gemm1_sec'):.1f}%, "
      f"gemm2 {_wpct('W_gemm2_sec'):.1f}%, "
      f"lam+D {_wpct('W_einsumD_sec'):.1f}%")
