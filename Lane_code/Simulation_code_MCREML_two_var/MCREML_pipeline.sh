#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline (additive + dominance, 3-component model).
# Both generation and estimation use additive + dominance + noise:
#   V = s2a*K_a + s2d*K_d + s2e*I ,   K_a = Z_a Z_a'/m ,   K_d = Z_d Z_d'/m .
# SIMULATION forms K_a, K_d and their Cholesky factors La (La La' = s2a K_a) and
# Ld (Ld Ld' = s2d K_d) to draw correctly-correlated additive / dominance
# effects; ESTIMATION is MATRIX-FREE (implicit GRM mat-vec + conjugate gradient
# + Hutchinson stochastic trace).  See MCREML_two_var.typ for the derivation.
# Mirrors the additive pipeline: a 4-step SLURM dependency chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup
# Self-contained model directory (own Cholesky/Phenotype/result).

N=8000
M=1000
S2A=0.5
S2D=0.25
S2E=0.25
MODE=Contiguous
ITERS=30
NMC=100
MEM=6G
ARRAY=20
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_two_var

FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2d${S2D}_s2e${S2E}

mkdir -p $DIR/error

LA_FILE=$DIR/Cholesky/La_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2e${S2E}.npy
LD_FILE=$DIR/Cholesky/Ld_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2e${S2E}.npy
PHENO_DIR=$DIR/Phenotype/y_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2e${S2E}

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
    $DIR/Cholesky.sh $N $M $S2A $S2D $S2E $MODE)
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
    $DIR/Phenotype.sh $N $M $S2A $S2D $S2E $MODE)
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
    $DIR/MCREML.sh $N $M $S2A $S2D $S2E $MODE $ITERS $NMC)
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
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LA_FILE $LD_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
