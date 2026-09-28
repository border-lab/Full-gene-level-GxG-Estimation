#!/bin/bash

N=400
M=4800
S2A=0.1
S2E=0.9
MODE=ContiguousSNP
ITERS=30
NMC=100                 # FINE phase; mc_reml runs a coarse S=15 phase first
VERBOSE=--verbose
MEM=2G
MEM_CHOL=48G
MEM_PHENO=24G
ARRAY=50                 # concurrent array tasks
PARTITION=statgen-gpu
NODE=lanec2-4-1          # every job runs on this node only

# Start step: 1 = all, 2 = Cholesky done, 3 = phenotypes done, 4 = combine only.
# A skipped step is assumed FINISHED (no SLURM dependency on it), so check the
# files first.  Usage: bash MCREML_pipeline.sh 2
START=${1:-1}

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive

# No s2d / s2gxg, no G and no R in the keys: this pipeline has only the
# additive kernel, so there is no gene split and no truncation level.
TAG=${MODE}_s2a${S2A}_s2e${S2E}_n${N}_m${M}
FILENAME=${MODE}_s2a${S2A}_s2e${S2E}_n${N}m${M}

mkdir -p $DIR/error

# Cholesky factor and phenotypes are named from TAG by the jobs that write
# them and by combine_code.sh, which deletes them.

# Step 1: Cholesky (single job) -- factors the dense GRM K_a = Z_a Z_a'/m.  The
# GRM is discarded after factorisation; the estimator applies K_a matrix-free
# from the genotype.
DEP=""
if [ "$START" -le 1 ]; then
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml \
    -p $PARTITION \
    --nodelist=$NODE \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_CHOL \
    --ntasks=1 \
    --cpus-per-task=64 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $S2A $S2E $MODE)
echo "Cholesky job: $JOB1"
DEP="--dependency=afterok:$JOB1"
fi

# Step 2: Phenotype (array, after Step 1).  Draws y = g_a + e from the factor.
if [ "$START" -le 2 ]; then
JOB2=$(sbatch --parsable $DEP \
    --job-name=pheno_mcreml \
    -p $PARTITION \
    --nodelist=$NODE \
    --error=$DIR/error/phenotype_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM_PHENO \
    --ntasks=1 \
    --cpus-per-task=1 \
    --array=1-50%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $S2A $S2E $MODE)
echo "Phenotype job: $JOB2"
DEP="--dependency=afterok:$JOB2"
fi

# Step 3: MC-AI-REML (array, after Step 2) -- matrix-free, genotype only, no
# cached GRM.  K_a u is an exact gemm pair; tr(V^{-1} K_a) comes from
# Hutchinson probes ($NMC CG solves per iteration).
if [ "$START" -le 3 ]; then
JOB3=$(sbatch --parsable $DEP \
    --job-name=mcreml \
    -p $PARTITION \
    --nodelist=$NODE \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=$DIR/error/mcreml_%A_%a.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --array=1-50%$ARRAY \
    --time=48:00:00 \
    $DIR/MCREML.sh $N $M $S2A $S2E $MODE $ITERS $NMC $VERBOSE)
echo "MC-AI-REML job: $JOB3  (2 VC: s2a K_a + s2e I; trace estimator: hutchinson, Nmc=$NMC)"
DEP="--dependency=afterok:$JOB3"
fi

# Step 4: Combine and clean up (after Step 3).  Everything lives in
# combine_code.sh:
#
#   result/$FILENAME/rep*.txt  -> result/$FILENAME.txt
#   time/rep_times/$FILENAME/  -> time/result/timing_$FILENAME.txt
#   time/op_counts/$FILENAME/  -> time/result/opcount_$FILENAME.txt
#
# The three reductions run independently, each per-rep directory is removed only
# after its own reduction lands, and combine_code.sh exits non-zero on any
# failure.  All three reductions are plain shell + awk inside
# combine_code.sh, so the combine job needs no python and no conda env.
JOB4=$(sbatch --parsable $DEP \
    --job-name=combine_mcreml \
    -p $PARTITION \
    --nodelist=$NODE \
    --error=$DIR/error/combine_%j.err \
    --output=$DIR/error/combine_%j.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME $TAG")
echo "Combine job: $JOB4"
