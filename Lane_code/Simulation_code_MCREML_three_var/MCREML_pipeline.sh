#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline
# (additive + dominance + pairwise-epistasis, 4-component model).
#   V = s2a*K_a + s2d*K_d + s2gxg*W + s2e*I ,
#   K_a = Z_a Z_a'/m ,  K_d = Z_d Z_d'/m ,  W = (1/p) sum_{a<b} std(Z_a.Z_b)(...)'.
# SIMULATION forms K_a, K_d, W and their Cholesky factors La, Ld, Lgxg to draw
# correctly-correlated effects, and saves W once (it is deterministic per
# genotype).  ESTIMATION applies additive / dominance via implicit mat-vecs and
# W via the prebuilt DENSE W @ b (O(n^2 c) per CG iteration) + conjugate gradient
# + Hutchinson trace.  See MCREML_three_var.typ.
#
# NOTE ON COST: reusing the prebuilt dense W avoids the O(n m^2) matrix-free
# compute_WU inside every CG iteration (~2500x faster at n = m = 1000).  The
# storage-free compute_WU remains in Function_MCREML.py for the large-m,
# memory-bound regime.
#
# Mirrors the additive / two_var pipelines: a 4-step SLURM dependency chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=1000
M=1000
S2A=0.3
S2D=0.3
S2GXG=0.3
S2E=0.1
MODE=RandomSNP
ITERS=30
NMC=100
MEM=2G
ARRAY=20
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_three_var

TAG=${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}
FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}

mkdir -p $DIR/error

LA_FILE=$DIR/Cholesky/La_${TAG}.npy
LD_FILE=$DIR/Cholesky/Ld_${TAG}.npy
LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job)
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
    $DIR/Cholesky.sh $N $M $S2A $S2D $S2GXG $S2E $MODE)
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
    $DIR/Phenotype.sh $N $M $S2A $S2D $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2)
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
    $DIR/MCREML.sh $N $M $S2A $S2D $S2GXG $S2E $MODE $ITERS $NMC)
echo "MC-AI-REML job: $JOB3"

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
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LA_FILE $LD_FILE $LGXG_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
