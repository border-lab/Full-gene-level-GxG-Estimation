#!/bin/bash

N=1000
M=1000
G=10
S2A=0.1
S2D=0.1
S2GXG=0.1
S2E=0.7
MODE=ContiguousSNP
ITERS=30
NMC=100                 # FINE phase; mc_reml runs a coarse S=15 phase first
# Which epistasis kernel REML fits, PRECOMPUTED by the Cholesky job:
#   R=0   the EXACT c-normalized W the phenotype is drawn from (no truncation)
#   R>0   the dense matrix of the _matfree_ sibling's rank-R operator
# The phenotype is drawn from the exact W either way.
R=0
# Set to --verbose to print the AI-REML trace (one line per iteration: the
# four components and max|step|) into error/mcreml_<jobid>_<task>.out, one
# file per replicate.  Leave EMPTY to keep the trace off.
VERBOSE=--verbose
# The REML job now HOLDS the dense n-by-n W (8 N^2 bytes), so its memory scales
# with N^2: twice W's size plus 2G of headroom (2G at N=1000, 6G at N=16000).
MEM=$(( (N * N * 16) / 1000000000 + 2 ))G
MEM_COMBINE=2G
MEM_CHOL=16G
MEM_PHENO=8G
ARRAY=50
PARTITION=statgen-gpu

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC

# Kernel tag: "rexact" for R<=0, "r<R>" otherwise.  Must match
# w_cache_tag() in Function_MCREML.py -- it keys the cached W and FILENAME.
if [ "$R" -le 0 ]; then
    RTAG=rexact
else
    RTAG=r${R}
fi

TAG=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}_${RTAG}

mkdir -p $DIR/error

# The Cholesky factors La/Ld/Lgxg_${TAG}.npy and Phenotype/y_${TAG} are named
# from TAG (no kernel tag: they do not depend on R) by the jobs that write them
# and by combine_code.sh, which deletes them.  The precomputed kernel is
# W/W_${MODE}_n${N}_m${M}_G${G}_${RTAG}.npy (+ .json) and is KEPT.

# Step 1: Cholesky (single job) -- builds La from the dense additive GRM
# K_a = Z_a Z_a'/m, Ld from the dense dominance GRM K_d = Z_d Z_d'/m, and Lgxg
# from the dense pooled UNSTANDARDIZED, C-NORMALIZED GRM W = W_raw / c-hat (the
# EXACT kernel).  c-hat is the O(nm) THIRD-MOMENT plug-in
# (Function_MCREML.C_METHOD = 'moment') -- no flag for it, on purpose.  The job
# also writes all three c routes to result/c_$FILENAME.txt.
#
# UNLIKE the _matfree_ sibling it also PRECOMPUTES THE ESTIMATION KERNEL: with
# R=0 the very W it factorized is written to W/ before being freed; with R>0 the
# dense rank-R W-hat is built and written instead.  Overwritten atomically on
# every run.  K_a and K_d are not cached (REML applies them matrix-free).
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml \
    -p $PARTITION \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_CHOL \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE $R)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1).  Draws y = g_a + g_d + g_gxg
# + e from the three factors, each component forced to its target variance.
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    --job-name=pheno_mcreml \
    -p $PARTITION \
    --error=$DIR/error/phenotype_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_PHENO \
    --ntasks=1 \
    --cpus-per-task=1 \
    --array=1-50%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2) -- loads the PRECOMPUTED W
# written by Step 1 (refusing it if its sidecar's (mode, n, m, G, R, c-hat)
# does not match) and applies it as a dense W @ B inside CG; K_a u and K_d u
# from a gemm pair against their design (exact); the score traces
# tr(V^{-1} K_i) from Hutchinson probes ($NMC CG solves per REML iteration).
# s2gxg_hat comes out ON THE REALIZED-VARIANCE SCALE -- no c is applied
# afterwards.
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=mcreml \
    -p $PARTITION \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=$DIR/error/mcreml_%A_%a.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --array=1-50%$ARRAY \
    --time=12:00:00 \
    $DIR/MCREML.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE $ITERS $NMC $R $VERBOSE)
echo "MC-AI-REML job: $JOB3  (4 VC: s2a K_a + s2d K_d + s2gxg W + s2e I, W = W_raw/c c-normalized; W apply: PRECOMPUTED dense, kernel=$RTAG; trace estimator: hutchinson, Nmc=$NMC)"

# Step 4: Combine results and clean up (waits for Step 3).  ALL of it lives in
# combine_code.sh:
#
#   result/$FILENAME/rep*.txt        -> result/$FILENAME.txt
#   time/rep_times/$FILENAME/        -> time/result/timing_$FILENAME.txt
#   time/op_counts/$FILENAME/        -> time/result/opcount_$FILENAME.txt
#
# Each part runs INDEPENDENTLY, each per-rep directory is removed only once its
# own reduction has landed, and the script exits non-zero if any part failed.
# ONLY ONE EXTRA FILE IS NEEDED ON THE CLUSTER for the reductions:
# time/average_time.py.
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_mcreml \
    -p $PARTITION \
    --error=$DIR/error/combine_%j.err \
    --output=$DIR/error/combine_%j.out \
    --nodes=1 \
    --mem=$MEM_COMBINE \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME $TAG")
echo "Combine job: $JOB4"
