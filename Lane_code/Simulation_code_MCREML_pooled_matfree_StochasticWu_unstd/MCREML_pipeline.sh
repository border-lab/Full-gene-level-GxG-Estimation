#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline
# (POOLED within-gene pairwise-epistasis, 2-component model).
# UNSTANDARDIZED-KERNEL variant of MCREML_pooled_matfree_StochasticWu.
#   V = s2gxg*W + s2e*I ,
#   W = (1/P) sum_{g=1}^G sum_{a<b in g} h_ab h_ab' ,  h_ab = Z_a .* Z_b ,
#   P = sum_g C(m_g, 2)  (total within-gene pairs; per-PAIR equal weight).
# The m SNPs are split into G contiguous genes (several Z, one per gene); only
# WITHIN-gene pairs enter W.  This is the "Pooled Model" of generative_model.typ.
#
# WHAT IS DIFFERENT FROM StochasticWu.  The interaction column h_ab is RAW --
# no mean-centering, no 1/sigma_ab scaling.  That single change removes the one
# defect the standardized pipeline could not fix: there, the O(n m Nw) operator
# can only reach the centered-unscaled kernel while the simulation drew from the
# standardized one, so the estimator carried a deterministic KERNEL BIAS on top
# of its Monte-Carlo error.  Here simulation and estimation are the same object:
#
#   W u = 1/(2P) [ (K_w .* K_w) u - D (D' u) ] ,   K_w = Z Z' ,  D = Z .* Z
#       ~ 1/(2P) [ 1/Nw ((Z(Z'(u .* (Z(Z'V))))) .* V) 1_Nw - D (D' u) ] ,
#
# with V in R^{n x Nw} i.i.d. Rademacher, frozen once so W-hat is a FIXED matrix
# (CG and Lanczos both require that).  The second line is the ONLY W apply the
# module implements -- the stochastic operator with the diagonal estimator -- so
# NW is the only knob this variant adds.  The exact first line is not a code
# path in the pipeline at all: verify_unstd.py carries its own independent
# implementation of it, purely to check the operator against.
#
# Verified by verify_unstd.py at n=300, m=120, G=4:
#   exact operator vs dense W          rel 5e-16
#   stochastic operator vs dense W     rel 0.243 * sqrt(n/Nw)   (flat in m)
#   mc_reml vs the exact likelihood optimum, per replicate, within trace noise.
#
# TWO PROPERTIES ARE DELIBERATELY GIVEN UP, both consequences of no centering:
# W 1 != 0 (so lambda_min(W) = 0 is an inequality, not an identity -- the SLQ
# conditioning guard is re-derived from PSD-ness) and tr(W) != n.  In exchange
# the spectrum of W is far more spread out than the standardized kernel's, which
# moves W AWAY from I -- the near-collinearity of W and I is exactly what makes
# the 2-component AI matrix near-singular and h^2 weakly identified at large m.
#
# COST.  8 n m_g Nw per right-hand side per gene, against 2 n m_g^2 for an exact
# apply and O(n^2) storage for the dense pre-computed W of _preW.  Nothing
# m-by-m is ever formed.
#
# SIMULATION still forms the dense W (a Cholesky factor needs an explicit
# matrix) but does NOT cache it; dropping the standardization also lets it use
# the same Hadamard identity, so the build is O(n^2 m) here against the
# sibling's O(n^2 m^2) pair loop.  ESTIMATION loads only the n-by-m genotype.
#
# Mirrors the additive / two_var / three_var / preW pipelines: a 4-step SLURM
# chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=1000
M=1000
G=10
S2GXG=0.2
S2E=0.8
MODE=RandomSNP
ITERS=30
NMC=100
# Score-trace estimator.  slq = stochastic Lanczos quadrature: Lanczos on W is
# run once and, because V = s2gxg*W + s2e*I is affine in W, its nodes/weights
# are reused at every variance setting -- so the Nmc probe solves leave the REML
# iteration entirely.  hutchinson = the original Nmc-solves-per-iteration path.
# Identical estimates (agree to the 1e-8 REML stopping tolerance).
TRACE_METHOD=slq
SLQ_K=25
# Frozen Rademacher probes in the W u apply.  Measured relative operator error
# is 0.24 sqrt(N/NW) and FLAT in m, so NW must scale with n, not with m, to hold
# accuracy -- at N=1000 the default NW=200 gives ~0.54, which is why an NW sweep
# is the natural experiment here.  NOTE that sweeping NW as written overwrites
# its own results: FILENAME carries no NW, and step 4 deletes the Cholesky
# factor and phenotypes.  To sweep, add _Nw${NW} to FILENAME and drop the
# LGXG_FILE / PHENO_DIR removals from step 4, so every NW is scored on the SAME
# phenotypes and the spread is operator noise alone.
NW=200
MEM=2G
ARRAY=20
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_StochasticWu_unstd

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg from the dense pooled
# UNSTANDARDIZED GRM W.  Nothing is cached for estimation: W is discarded after
# the factorisation.
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml \
    -p $PARTITION \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=16 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $G $S2GXG $S2E $MODE)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1)
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    --job-name=pheno_mcreml \
    -p $PARTITION \
    --error=$DIR/error/phenotype_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --array=1-300%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $G $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2) -- matrix-free, reads the
# genotype only.  No cached W is needed, so this step has no dependency on any
# n-by-n artifact.  W u comes from the stochastic operator at $NW probes and the
# score traces tr(V^{-1} K_i) from $TRACE_METHOD.
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=mcreml \
    -p $PARTITION \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --array=1-300%$ARRAY \
    --time=12:00:00 \
    $DIR/MCREML.sh $N $M $G $S2GXG $S2E $MODE $ITERS $NMC $TRACE_METHOD $SLQ_K $NW)
echo "MC-AI-REML job: $JOB3  (W apply: stochastic, Nw=$NW; trace estimator: $TRACE_METHOD, slq_k=$SLQ_K)"

# Step 4: Combine results and clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_mcreml \
    -p $PARTITION \
    --error=$DIR/error/combine_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME && python3 $DIR/time/average_time.py $DIR/time/rep_times/$FILENAME $DIR/time/result/timing_$FILENAME.txt && rm -rf $DIR/result/$FILENAME && rm -rf $DIR/time/rep_times/$FILENAME && rm -f $LGXG_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
echo "Average estimation time -> $DIR/time/result/timing_$FILENAME.txt"
echo "Kernel W build time      -> $DIR/time/result/W_timing_$FILENAME.txt  (simulation only)"
