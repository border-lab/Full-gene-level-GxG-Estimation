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
parser.add_argument('--nmc', type=int, default=100)   # FINE phase; coarse is 15
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score trace tr(V^{-1} W) is Hutchinson's, with Nmc CG solves per REML
# iteration.
# Truncation level of the low-rank W u apply -- the ONE knob this variant has,
# and the only accuracy control there is.  Per gene, the top-r eigenpairs of
# K_g = Z_g Z_g' (from the thin SVD of Z_g) replace one copy of the Hadamard
# square, so the apply is O(n m_g r) and DETERMINISTIC: the error is a
# truncation BIAS set by the discarded eigenvalue tail, not sampling noise --
# small under LD, large under linkage equilibrium -- and only raising r
# reduces it.  r >= min(n, m_g) makes the apply EXACT.  With no other kernel
# in the model, r biases s2gxg-hat and s2e-hat alone.
parser.add_argument('--r', type=int, default=R_DEFAULT)
# Print the AI-REML trace (one line per iteration: the two components and
# the predicted gain) to STDOUT.  Off by default -- with it on, Step 3 must
# also be launched with a real --output, since the pipeline sends stdout to
# /dev/null.
parser.add_argument('--verbose', action='store_true')
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

# Load the GENOTYPE and build the standardized design.  MC AI-REML applies the
# pooled within-gene epistasis GRM MATRIX-FREE from the G contiguous genes Z_a
# is split into inside mc_reml.  NO n-by-n GRM is loaded or stored (only the
# n-by-m design).
#
# The epistasis kernel is the UNSTANDARDIZED, C-NORMALIZED one (h_ab =
# Z_a .* Z_b, the whole kernel divided by c-hat = pooled_c, the O(nm)
# third-moment plug-in C_METHOD selects), exactly the kernel
# Simulate_Cholesky.py drew the phenotype from -- the divisor is recomputed
# here by the same deterministic function of the same genotype, so it is the
# same number, bit for bit.  So any bias seen in the results here is the
# rank-r TRUNCATION bias, NOT a kernel mismatch, and it is deterministic:
# every replicate fits the same W-hat, and the gap closes only by raising r.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)

# Load phenotype (s2gxg_s2e order)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2gxg W + s2e I, W applied matrix-free by the
# rank-r truncation.  The genotype-only setup (per-gene SVD factors, pair
# total P, c-hat) is built once inside mc_reml, outside the REML iteration.
# seed=rep fixes the Hutchinson probe randomness per replicate (reproducible,
# yet varied).  There is NO operator seed: W-hat is a deterministic function
# of (genotype, r).
# Time the estimation only (the MC_REML call) for the per-rep timing record --
# this INCLUDES the per-gene SVD setup, which is part of the estimator's cost
# here.
# NOTE: timings are NOT comparable to the 4VC parent's totals, which solve a
# larger system with two more kernel applies per V apply; the comparison is
# made through the per-column seconds below.  Against the _additive sibling,
# W_sec_per_col here versus K_sec_per_col there is the cost of the epistasis
# kernel per column relative to the additive one.
t_start = time.perf_counter()
if args.verbose:
    print(f"--- rep{rep}: AI-REML trace, var(y)={y.var():.6f}, "
          f"columns s2gxg s2e ---", flush=True)
s2gxg_hat, s2e_hat, _ = MC_REML(Z, y, G, iters=iters, Nmc=nmc, seed=rep, r=r,
                                verbose=args.verbose)
elapsed = time.perf_counter() - t_start
# How many LINEAR OPERATOR APPLIES that fit actually cost, and how long each
# kind took.  mc_reml zeroes the counters on entry, so this snapshot is THIS
# replicate and nothing else -- the setup applies (W U for the fixed probes)
# included, since they are part of what estimating this replicate took.
# Written below next to the wall-clock time and averaged over replicates by
# the combine job.
op_counts = get_op_counts()

# NO post-fit realized-variance correction.  The same c-hat is already INSIDE
# the kernel, both when the phenotype was drawn and when it is fitted, so
# V_l = s2gxg-hat with the correction factor equal to 1 by construction.
# Applying compute_c_pooled here would double-count the normalization.

# Save result, 2 columns -- the FIT and nothing else:
#
#   (s2gxg_hat, s2e_hat)
#    V_gamma    V_e
#
# Both are on the realized-variance scale and directly comparable to the
# nominal targets: with S2GXG = 0.1 and S2E = 0.9 the two column means should
# sit at 0.1 / 0.9, with no rescaling and no per-replicate reference to pair
# against (the Phenotype step forced the realized variances onto the targets).
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/result/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
run_tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/time/rep_times/{run_tag}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")

# Record this replicate's OPERATOR-APPLY COUNTS AND TIMERS, one file per
# replicate, as labelled "<name> <value>" lines, read back by name so the row
# cannot scramble if this list ever grows.  The combine step averages all reps
# into time/result/opcount_<FILENAME>.txt.
#
# WHY NEXT TO THE TIME.  The wall-clock number above is what this machine took;
# the counts are what the ALGORITHM did, and they are deterministic -- rerun
# the replicate anywhere and they come out identical, because CG's iteration
# count is fixed by (V, rhs, tol).  So a timing difference between two
# pipelines, or two r, or two Nmc, splits cleanly: the counts say how much
# work was asked for, the seconds say how fast the machine did it.  V_columns
# is the one to quote as "the" cost -- every operator here is linear in the
# column count -- and V_applies next to it says how well that work was
# batched (a solve with Nmc right-hand sides is ONE apply and Nmc columns).
#
# WALL-CLOCK TWINS (keys ending "_sec", seconds on THIS machine).  They NEST --
# V_sec contains the W time spent inside V applies, W_sec contains the four
# W_* pieces, each solve_*_sec contains its V applies -- so they are not
# summed; the *_sec_per_col lines (seconds per n-vector) are the comparable
# unit for the two operators.  The W_* pieces answer how much of a W column is
# gemm (W_gemm1: Zg.T @ Tb, W_gemm2: Zg @ .) versus elementwise (W_bcast: the
# Q*U broadcast, W_einsumD: the einsum and the D correction) -- the
# elementwise share being the part BLAS threads cannot touch.  Solve
# groups: solve_y (V^-1 y), solve_probe_coarse / _fine (V^-1 U, split by
# phase), solve_ai (V^-1 [W x, x]); setup is svd + spectral + W U (the last
# kept under the parent's setup_KU name so opcount files line up);
# phase_coarse / phase_fine stamp the switch; lam_min_exact is the exact
# feasibility fallback.
#
# The keys are written in a FIXED order with 0 defaults, so every rep file has
# the same lines whether or not a counter was ever bumped.  The list is the
# 4VC parent's with the K keys removed (no K_a, no K_d here).
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "W_applies", "W_columns",
           "lam_min_exact",
           "V_sec", "W_sec",
           "V_sec_per_col", "W_sec_per_col",
           "W_bcast_sec", "W_gemm1_sec", "W_gemm2_sec", "W_einsumD_sec",
           "solve_y_sec", "solve_probe_coarse_sec", "solve_probe_fine_sec",
           "solve_ai_sec",
           "setup_sec", "setup_svd_sec", "setup_spectral_sec", "setup_KU_sec",
           "phase_coarse_sec", "phase_fine_sec", "lam_min_exact_sec")
# Seconds per column for the two operators: the machine-dependent unit cost of
# one n-vector through each apply.  Nested as above: V includes its W.
for _op in ("V", "W"):
    _cols = op_counts.get(f"{_op}_columns", 0)
    op_counts[f"{_op}_sec_per_col"] = (
        op_counts.get(f"{_op}_sec", 0.0) / _cols if _cols else 0.0)
op_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/time/op_counts/{run_tag}"
os.makedirs(op_dir, exist_ok=True)
with open(f"{op_dir}/rep{rep}.txt", 'w') as f:
    for key in OP_KEYS:
        f.write(f"{key} {op_counts.get(key, 0)}\n")
print(f"operator applies this replicate: V={op_counts.get('V_applies', 0)} "
      f"({op_counts.get('V_columns', 0)} columns), "
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
      f"(svd {_g('setup_svd_sec'):.2f}, spectral {_g('setup_spectral_sec'):.2f}, "
      f"WU {_g('setup_KU_sec'):.2f}); "
      f"phase coarse {_g('phase_coarse_sec'):.2f} s, fine {_g('phase_fine_sec'):.2f} s; "
      f"solves: y {_g('solve_y_sec'):.2f}, probe coarse "
      f"{_g('solve_probe_coarse_sec'):.2f}, probe fine "
      f"{_g('solve_probe_fine_sec'):.2f}, AI {_g('solve_ai_sec'):.2f}; "
      f"lam_min_exact {_g('lam_min_exact_sec'):.2f} s")
print(f"seconds per column: V {_g('V_sec_per_col'):.3e}, "
      f"W {_g('W_sec_per_col'):.3e}; "
      f"W apply split: bcast {_wpct('W_bcast_sec'):.1f}%, "
      f"gemm1 (Zg.T@Tb) {_wpct('W_gemm1_sec'):.1f}%, "
      f"gemm2 (Zg@.) {_wpct('W_gemm2_sec'):.1f}%, "
      f"einsum+D {_wpct('W_einsumD_sec'):.1f}%")
