from Function_MCREML import *
import argparse
import pandas as pd
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--iters', type=int, default=30)
parser.add_argument('--nmc', type=int, default=100)   # FINE phase; coarse is 15
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)
# The score traces tr(V^{-1} K_i) are Hutchinson's, with Nmc CG solves per
# REML iteration.  There is no --r and no --A_dtype: the only genetic kernel is
# K_a, applied exactly, so there is no truncation and no stored-A apply.
# Print the AI-REML trace (one line per iteration: the two components and
# max|step|) to STDOUT.  Off by default -- with it on, Step 3 must also be
# launched with a real --output, since the pipeline sends stdout to /dev/null.
parser.add_argument('--verbose', action='store_true')
args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2e = args.s2e
iters = args.iters
nmc = args.nmc
rep = args.rep
mode = args.mode

# Load the GENOTYPE and build the standardized design.  MC AI-REML applies the
# additive GRM MATRIX-FREE straight from it, K_a B = Z_a(Z_a'B)/m.  NO n-by-n
# GRM is loaded or stored (only the n-by-m design), and it is the SAME
# transformation the Cholesky job used, so K_a is identical on both sides:
# any bias seen here is the estimator's, not a kernel mismatch.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)

# Load phenotype (s2a_s2e order)
tag = f"{mode}_s2a{s2a}_s2e{s2e}_n{n}_m{m}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2a K_a + s2e I, K_a applied matrix-free and
# exactly.  seed=rep fixes the Hutchinson probe randomness per replicate
# (reproducible, yet varied).  Time the estimation only (the MC_REML call) for
# the per-rep timing record.
# NOTE: timings are NOT comparable to the three- or four-component pipelines',
# which solve a larger system with more kernel applies per iteration.
t_start = time.perf_counter()
if args.verbose:
    print(f"--- rep{rep}: AI-REML trace, var(y)={y.var():.6f}, "
          f"columns s2a s2e ---", flush=True)
s2a_hat, s2e_hat, _, info = MC_REML(
    Z, y, iters=iters, Nmc=nmc, seed=rep, verbose=args.verbose)
elapsed = time.perf_counter() - t_start
# How many LINEAR OPERATOR APPLIES that fit actually cost.  mc_reml zeroes the
# counters on entry, so this snapshot is THIS replicate and nothing else -- the
# setup apply (K_a U for the fixed probes) included.  Written below next to the
# wall-clock time and averaged over replicates by the combine job.
op_counts = get_op_counts()

# Save result, 2 columns -- the FIT and nothing else:
#
#   (s2a_hat, s2e_hat)
#
# BOTH ARE ON THE REALIZED-VARIANCE SCALE and directly comparable to the
# nominal targets (the design is standardized, c_a = 1 exactly): with
# S2A = 0.1 and S2E = 0.9 the column means should sit at 0.1 / 0.9.  The
# optimizer's diagnostics (converged, rejected steps, components at their
# bound) go to this job's stdout log only.
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/result/{mode}_s2a{s2a}_s2e{s2e}_n{n}m{m}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2a_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
run_tag = f"{mode}_s2a{s2a}_s2e{s2e}_n{n}m{m}"
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/time/rep_times/{run_tag}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")

# Record this replicate's OPERATOR-APPLY COUNTS, one file per replicate, as
# labelled "<name> <value>" lines, read back by name.  The combine step
# averages all reps into time/result/opcount_<FILENAME>.txt.
#
# The counts are what the ALGORITHM did and are deterministic (CG's iteration
# count is fixed by (V, rhs, tol)); the "_sec" keys are wall-clock seconds on
# THIS machine.  The timers NEST -- V_sec contains the K time spent inside V
# applies, each solve_*_sec contains its V applies -- so they are not summed;
# the *_sec_per_col lines (seconds per n-vector) are the comparable unit.
# Solve groups: solve_probe_coarse / _fine (the merged solve V^-1 [y, U], split
# by phase), solve_ai (V^-1 [K_a x, x]); setup is K_a U; phase_coarse /
# phase_fine stamp the switch; cg_negcurv counts CG solves stopped because V
# was not positive definite (a guard; expected 0 here).
#
# The keys are written in a FIXED order with 0 defaults, so every rep file has
# the same lines whether or not a counter was ever bumped.
OP_KEYS = ("reml_iters", "cg_solves", "cg_iters",
           "V_applies", "V_columns",
           "K_applies", "K_columns",
           "cg_negcurv",
           "V_sec", "K_sec",
           "V_sec_per_col", "K_sec_per_col",
           "solve_probe_coarse_sec", "solve_probe_fine_sec",
           "solve_ai_sec",
           "setup_sec", "setup_KU_sec",
           "phase_coarse_sec", "phase_fine_sec")
# Seconds per column for the two operators: the machine-dependent unit cost
# of one n-vector through each apply.  Nested as above: V includes its K.
for _op in ("V", "K"):
    _cols = op_counts.get(f"{_op}_columns", 0)
    op_counts[f"{_op}_sec_per_col"] = (
        op_counts.get(f"{_op}_sec", 0.0) / _cols if _cols else 0.0)
op_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/time/op_counts/{run_tag}"
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
      f"setup {_g('setup_sec'):.2f} s (K_aU {_g('setup_KU_sec'):.2f}); "
      f"phase coarse {_g('phase_coarse_sec'):.2f} s, fine {_g('phase_fine_sec'):.2f} s; "
      f"solves: [y,U] coarse "
      f"{_g('solve_probe_coarse_sec'):.2f}, [y,U] fine "
      f"{_g('solve_probe_fine_sec'):.2f}, AI {_g('solve_ai_sec'):.2f}; "
      f"cg_negcurv {op_counts.get('cg_negcurv', 0)}")
print(f"optimizer: converged={info['converged']}, iters={info['n_iters']}, "
      f"rejected={info['n_reject']}, switch at iter {info['switch_it']}, "
      f"at_bound={info['at_bound']}")
print(f"seconds per column: V {_g('V_sec_per_col'):.3e}, "
      f"K {_g('K_sec_per_col'):.3e}")
