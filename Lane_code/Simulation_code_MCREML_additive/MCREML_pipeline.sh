#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline (additive-only, 2-component model).
# Both generation and estimation use additive + noise only:
#   V = s2a*K + s2e*I ,   K = Z Z' / m .
# SIMULATION forms K and its Cholesky factor La (La La' = s2a K) to draw a
# correctly-correlated additive effect; ESTIMATION is MATRIX-FREE (implicit GRM
# mat-vec + conjugate gradient + Hutchinson stochastic trace).  See
# MCREML_additive.typ for the derivation.
# Mirrors the AI-REML pipeline: a 4-step SLURM dependency chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup
# Self-contained model directory (own Cholesky/Phenotype/result).

N=2000
M=1000
S2A=0.5
S2E=0.5
MODE=Random
ITERS=30
NMC=100
MEM=16G
ARRAY=30
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_additive

FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2e${S2E}

mkdir -p $DIR/error

LA_FILE=$DIR/Cholesky_La/La_${MODE}_n${N}_m${M}_s2a${S2A}_s2e${S2E}.npy
PHENO_DIR=$DIR/Phenotype/y_${MODE}_n${N}_m${M}_s2a${S2A}_s2e${S2E}

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
    $DIR/Cholesky.sh $N $M $S2A $S2E $MODE)
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
    $DIR/Phenotype.sh $N $M $S2A $S2E $MODE)
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
    $DIR/MCREML.sh $N $M $S2A $S2E $MODE $ITERS $NMC)
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
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LA_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
