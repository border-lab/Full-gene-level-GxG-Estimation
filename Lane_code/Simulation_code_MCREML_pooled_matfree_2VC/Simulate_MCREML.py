from Function_MCREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=100)   # FINE phase; coarse is 15
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score traces tr(V^{-1} K_i) are Hutchinson's, with Nmc CG solves per
# REML iteration.  There is no --G and no --r: this pipeline has no epistasis
# kernel, so nothing is split into genes and nothing is truncated.  Both
# kernels are applied EXACTLY.
# Print the AI-REML trace (one line per iteration: the three components and
# the predicted gain) to STDOUT.  Off by default -- with it on, Step 3 must
# also be launched with a real --output, since the pipeline sends stdout to
# /dev/null.
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load the GENOTYPE and build the two standardized designs.  MC AI-REML applies
# BOTH genetic kernels MATRIX-FREE straight from them: the additive GRM as
# K_a B = Z_a(Z_a'B)/m and the dominance GRM as K_d B = Z_d(Z_d'B)/m.  NO
# n-by-n GRM is loaded or stored (only the two n-by-m designs).
#
# The dominance design is built from the RAW 0/1/2 dosages -- it needs the
# allele frequencies -- and is the SAME transformation the Cholesky job used,
# so K_d is identical on both sides; K_a likewise.  There is nothing else:
# the phenotype was drawn from exactly the two matrices being fitted, so any
# bias seen in the results here is the ESTIMATOR's, not a kernel mismatch.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)
Zd = dominance_design(SNP)

# Load phenotype (s2a_s2d_s2e order)
tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}_m{m}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2a K_a + s2d K_d + s2e I, both GRMs applied
# matrix-free from the designs.  seed=rep fixes the Hutchinson probe
# randomness per replicate (reproducible, yet varied).
# Time the estimation only (the MC_REML call) for the per-rep timing record.
# NOTE: timings are NOT comparable to the 4VC pipelines', which solve a larger
# system with one more (and far more expensive) kernel apply per V apply; that
# comparison is the point of this pipeline, and it is made through the
# per-column seconds below, not through the totals.
t_start = time.perf_counter()
if args.verbose:
    print(f"--- rep{rep}: AI-REML trace, var(y)={y.var():.6f}, "
          f"columns s2a s2d s2e ---", flush=True)
s2a_hat, s2d_hat, s2e_hat, _ = MC_REML(Z, Zd, y, iters=iters, Nmc=nmc,
                                       seed=rep, verbose=args.verbose)
elapsed = time.perf_counter() - t_start
# How many LINEAR OPERATOR APPLIES that fit actually cost, and how long each
# kind took.  mc_reml zeroes the counters on entry, so this snapshot is THIS
# replicate and nothing else -- the setup applies (K_a U, K_d U for the fixed
# probes) included, since they are part of what estimating this replicate
# took.  Written below next to the wall-clock time and averaged over
# replicates by the combine job.
op_counts = get_op_counts()

# Save result, 3 columns -- the FIT and nothing else:
#
#   (s2a_hat, s2d_hat, s2e_hat)
#    V_a      V_d      V_e
#
# All three are directly comparable to the nominal targets: with S2A = S2D =
# 0.1 and S2E = 0.7 the three column means should sit at 0.1 / 0.1 / 0.7, with
# no rescaling and no per-replicate reference to pair against (the Phenotype
# step forced the realized variances onto the targets).
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/result/{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}m{m}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2d_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
run_tag = f"{mode}_s2a{s2a}_s2d{s2d}_s2e{s2e}_n{n}m{m}"
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/time/rep_times/{run_tag}"
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
# pipelines or two Nmc splits cleanly: the counts say how much work was asked
# for, the seconds say how fast the machine did it.  V_columns is the one to
# quote as "the" cost -- every operator here is linear in the column count --
# and V_applies next to it says how well that work was batched (a solve with
# Nmc right-hand sides is ONE apply and Nmc columns).
#
# WALL-CLOCK TWINS (keys ending "_sec", seconds on THIS machine).  They NEST --
# V_sec contains the K time spent inside V applies, each solve_*_sec contains
# its V applies -- so they are not summed; the *_sec_per_col lines (seconds per
# n-vector) are the comparable unit for the two operators, and the number to
# put against the 4VC parent's V_sec_per_col to see what the epistasis kernel
# costs per column.  Solve groups: solve_y (V^-1 y), solve_probe_coarse /
# _fine (V^-1 U, split by phase), solve_ai (V^-1 [K_i x]); setup is spectral +
# K_i U (no svd here); phase_coarse / phase_fine stamp the switch;
# lam_min_exact is the exact feasibility fallback.
#
# The keys are written in a FIXED order with 0 defaults, so every rep file has
# the same lines whether or not a counter was ever bumped.
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "K_applies", "K_columns",
           "lam_min_exact",
           "V_sec", "K_sec",
           "V_sec_per_col", "K_sec_per_col",
           "solve_y_sec", "solve_probe_coarse_sec", "solve_probe_fine_sec",
           "solve_ai_sec",
           "setup_sec", "setup_spectral_sec", "setup_KU_sec",
           "phase_coarse_sec", "phase_fine_sec", "lam_min_exact_sec")
# Seconds per column for the two operators: the machine-dependent unit cost of
# one n-vector through each apply.  Nested as above: V includes its two K's.
for _op in ("V", "K"):
    _cols = op_counts.get(f"{_op}_columns", 0)
    op_counts[f"{_op}_sec_per_col"] = (
        op_counts.get(f"{_op}_sec", 0.0) / _cols if _cols else 0.0)
op_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/time/op_counts/{run_tag}"
os.makedirs(op_dir, exist_ok=True)
with open(f"{op_dir}/rep{rep}.txt", 'w') as f:
    for key in OP_KEYS:
        f.write(f"{key} {op_counts.get(key, 0)}\n")
print(f"operator applies this replicate: V={op_counts.get('V_applies', 0)} "
      f"({op_counts.get('V_columns', 0)} columns), "
      f"K={op_counts.get('K_applies', 0)} "
      f"({op_counts.get('K_columns', 0)} columns); "
      f"{op_counts.get('reml_iters', 0)} REML iters, "
      f"{op_counts.get('cg_iters', 0)} CG iters in "
      f"{op_counts.get('cg_solves', 0)} solves -> {op_dir}/rep{rep}.txt")
_g = lambda k: op_counts.get(k, 0.0)
print(f"timing this replicate ({elapsed:.2f} s total): "
      f"setup {_g('setup_sec'):.2f} s "
      f"(spectral {_g('setup_spectral_sec'):.2f}, K_iU {_g('setup_KU_sec'):.2f}); "
      f"phase coarse {_g('phase_coarse_sec'):.2f} s, fine {_g('phase_fine_sec'):.2f} s; "
      f"solves: y {_g('solve_y_sec'):.2f}, probe coarse "
      f"{_g('solve_probe_coarse_sec'):.2f}, probe fine "
      f"{_g('solve_probe_fine_sec'):.2f}, AI {_g('solve_ai_sec'):.2f}; "
      f"lam_min_exact {_g('lam_min_exact_sec'):.2f} s")
print(f"seconds per column: V {_g('V_sec_per_col'):.3e}, "
      f"K {_g('K_sec_per_col'):.3e}")
