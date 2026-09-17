#!/bin/bash

N=32000
M=16000
G=160
S2A=0.1
S2D=0.1
S2GXG=0.1
S2E=0.7
MODE=ContiguousSNP
ITERS=30
NMC=100                 # FINE phase; mc_reml runs a coarse S=15 phase first
R=30
VERBOSE=--verbose
MEM=48G
MEM_CHOL=96G
MEM_PHENO=48G
ARRAY=8                 # concurrent array tasks; each takes 64 cores
PARTITION=statgen-gpu

# Start step: 1 = all, 2 = Cholesky done, 3 = phenotypes done, 4 = combine only.
# A skipped step is assumed FINISHED (no SLURM dependency on it), so check the
# files first.  Usage: bash MCREML_pipeline.sh 2
START=${1:-1}

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC

TAG=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}

mkdir -p $DIR/error

# Cholesky factors and phenotypes are named from TAG by the jobs that write
# them and by combine_code.sh, which deletes them.

# Step 1: Cholesky (single job) -- factors the three dense GRMs: K_a = Z_a Z_a'/m,
# K_d = Z_d Z_d'/m, and the pooled unstandardized W = W_raw/c-hat (EXACT kernel,
# not the truncation).  c-hat is the O(nm) third-moment plug-in, fixed by the
# module constant Function_MCREML.C_METHOD = 'moment' (no flag, so simulation and
# estimation cannot diverge).  All three c routes plus c_exact/c-hat go to
# result/c_$FILENAME.txt.  The GRMs are discarded after factorisation; the
# estimator recomputes c-hat from the genotype.
DEP=""
if [ "$START" -le 1 ]; then
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml \
    -p $PARTITION \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_CHOL \
    --ntasks=1 \
    --cpus-per-task=64 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE)
echo "Cholesky job: $JOB1"
DEP="--dependency=afterok:$JOB1"
fi

# Step 2: Phenotype (array, after Step 1).  Draws y = g_a + g_d + g_gxg + e from
# the factors and records each replicate's realized variances, which Step 3
# carries into the estimate row.
if [ "$START" -le 2 ]; then
JOB2=$(sbatch --parsable $DEP \
    --job-name=pheno_mcreml \
    -p $PARTITION \
    --error=$DIR/error/phenotype_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_PHENO \
    --ntasks=1 \
    --cpus-per-task=64 \
    --array=1-100%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"
DEP="--dependency=afterok:$JOB2"
fi

# Step 3: MC-AI-REML (array, after Step 2) -- matrix-free, genotype only, no
# cached GRM.  W u is the rank-$R truncation of the same c-normalized kernel the
# phenotype was drawn from (same c from the same genotype, so the truncation is
# the only difference); K_a u / K_d u are exact gemm pairs; tr(V^{-1} K_i) comes
# from Hutchinson probes ($NMC CG solves per iteration).  s2gxg_hat is already on
# the realized-variance scale -- no c applied afterwards.
if [ "$START" -le 3 ]; then
JOB3=$(sbatch --parsable $DEP \
    --job-name=mcreml \
    -p $PARTITION \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=$DIR/error/mcreml_%A_%a.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=64 \
    --array=1-100%$ARRAY \
    --time=48:00:00 \
    $DIR/MCREML.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE $ITERS $NMC $R $VERBOSE)
echo "MC-AI-REML job: $JOB3  (4 VC: s2a K_a + s2d K_d + s2gxg W + s2e I, W = W_raw/c c-normalized; W apply: low-rank, r=$R; trace estimator: hutchinson, Nmc=$NMC)"
DEP="--dependency=afterok:$JOB3"
fi

# Step 4: Combine and clean up (after Step 3).  Everything -- combine, both
# reductions, every rm -- lives in combine_code.sh, not an sbatch --wrap chain:
#
#   result/$FILENAME/rep*.txt  -> result/$FILENAME.txt
#   time/rep_times/$FILENAME/  -> time/result/timing_$FILENAME.txt
#   time/op_counts/$FILENAME/  -> time/result/opcount_$FILENAME.txt
#
# The last holds operator-apply counts (applies of V, K_a/K_d and the rank-$R
# epistasis operator per replicate, and to how many columns).  Those are fixed by
# (genotype, r, $NMC, tolerances) and independent of BLAS threads or node load,
# so with the timings they separate work asked for from cluster speed.
#
# The three reductions run independently, each per-rep directory is removed only
# after its own reduction lands, and combine_code.sh exits non-zero on any
# failure (the old `&&` chain silently skipped later steps and the cleanups).
# Only time/average_time.py needs to be on the cluster; the op-count reduction is
# awk inside combine_code.sh, so it needs no python and no conda env.
JOB4=$(sbatch --parsable $DEP \
    --job-name=combine_mcreml \
    -p $PARTITION \
    --error=$DIR/error/combine_%j.err \
    --output=$DIR/error/combine_%j.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME $TAG")
echo "Combine job: $JOB4"
