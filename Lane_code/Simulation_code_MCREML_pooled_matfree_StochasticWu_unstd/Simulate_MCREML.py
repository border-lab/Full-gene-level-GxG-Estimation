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
# Score-trace estimator: 'slq' (stochastic Lanczos quadrature, default) runs
# Lanczos on W ONCE and reuses the nodes/weights for every variance setting;
# 'hutchinson' pays Nmc CG solves per REML iteration.  They target the same
# quantity with the same probes and agree to the REML stopping tolerance.
parser.add_argument('--trace_method', type=str, default='slq',
                    choices=['slq', 'hutchinson'])
parser.add_argument('--slq_k', type=int, default=25)
# Frozen Rademacher probes in the O(n m Nw) apply of W u -- the ONE knob this
# variant adds, and the only accuracy control there is: the relative operator
# error is 0.24 sqrt(n/Nw) and FLAT in m, so Nw tracks n, not m.  There is a
# single W-apply configuration (the note's stochastic operator with the diagonal
# estimator); verify_unstd.py carries an independent exact implementation to
# check it against, and that is the only place an alternative lives.
parser.add_argument('--nw', type=int, default=NW_DEFAULT)
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
trace_method = args.trace_method
slq_k = args.slq_k
nw = args.nw

# Load the GENOTYPE and column-standardize it.  MC AI-REML applies the pooled
# within-gene epistasis GRM W MATRIX-FREE straight from Z -- the m SNPs are
# split into G contiguous genes inside mc_reml.  NO n-by-n W is loaded or stored
# (only the n-by-m genotype), which is the whole point of this variant: the
# _preW pipeline loads a cached dense W here instead.
#
# The kernel is the UNSTANDARDIZED one (h_ab = Z_a .* Z_b, no centering, no
# 1/sigma_ab scaling), which is exactly the kernel Simulate_Cholesky.py drew the
# phenotype from.  So unlike the standardized StochasticWu pipeline, any bias
# seen in the results here is Monte Carlo in Nw, NOT a kernel mismatch -- there
# is no deterministic gap left for a larger Nw to fail to close.
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}_n{n}_m{m}.csv", header=None).to_numpy()
Z = additive_design(SNP)

# Load phenotype (s2gxg_s2e order)
tag = f"{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}_m{m}_G{G}"
y_path = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_StochasticWu_unstd/Phenotype/y_{tag}/rep{rep}.csv"
y = pd.read_csv(y_path, header=None).to_numpy().flatten()

# Monte-Carlo AI-REML: V = s2gxg W + s2e I, pooled W applied matrix-free from Z.
# The genotype-only setup (per-gene state, pair total P) is built once inside
# mc_reml, outside the REML iteration.
# seed=rep fixes the TRACE probe randomness per replicate (reproducible, yet
# varied).  The OPERATOR probes are deliberately left at the fixed default seed:
# W-hat must be the same matrix for every replicate, or reps would be fitting
# different kernels and their spread would mix operator noise with sampling
# noise.
# Time the estimation only (the MC_REML call) for the per-rep timing record --
# this INCLUDES the matrix-free setup AND, on the SLQ path, the one-off Lanczos
# run, which is part of the estimator's cost here (in the _preW pipeline the
# equivalent work sits in the Cholesky job).
t_start = time.perf_counter()
s2gxg_hat, s2e_hat, _ = MC_REML(Z, y, G, iters=iters, Nmc=nmc, seed=rep,
                                trace_method=trace_method, slq_k=slq_k, Nw=nw)
elapsed = time.perf_counter() - t_start

# Save result (s2gxg_s2e order)
output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_StochasticWu_unstd/result/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(output_dir, exist_ok=True)
filename = f"{output_dir}/rep{rep}.txt"
with open(filename, 'w') as f:
    f.write(f"({s2gxg_hat},{s2e_hat})\n")

# Record this replicate's estimation wall-clock time (seconds).  The combine
# step averages all reps into time/result/timing_<FILENAME>.txt.
time_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_StochasticWu_unstd/time/rep_times/{mode}_s2gxg{s2gxg}_s2e{s2e}_n{n}m{m}_G{G}"
os.makedirs(time_dir, exist_ok=True)
with open(f"{time_dir}/rep{rep}.txt", 'w') as f:
    f.write(f"{elapsed}\n")
