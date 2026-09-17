#!/bin/bash
# Exact AI-REML simulation pipeline (4-component model).
# Both generation and estimation use additive + dominance + GxG + noise:
#   V = s2a*K + s2d*D + s2gxg*W + s2e*I
# The dominance GRM keeps the additive/GxG estimates identifiable under LD
# (Hivert et al. 2021).  Set S2D=0.0 to keep dominance as an estimation-only
# nuisance term (phenotype then has no true dominance variance).
# Mirrors Only_W_pipeline / Sparse_pipeline: a 4-step SLURM dependency chain
#   Cholesky (1 job) -> Phenotype (array) -> AI-REML (array) -> combine + cleanup
# Self-contained model directory (own Cholesky/Phenotype/result, like Sparse_model).

N=4000
M=1000
S2A=0.2
S2D=0.0
S2GXG=0.7
S2E=0.1
MODE=Random
ITERS=30
MEM=16G
ARRAY=20

# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/AIREML

FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky_Lgxg/Lgxg_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}.npy
LA_FILE=$DIR/Cholesky_La/La_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}.npy
LD_FILE=$DIR/Cholesky_Ld/Ld_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}.npy
PHENO_DIR=$DIR/Phenotype/y_${MODE}_n${N}_m${M}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}

# Step 1: Cholesky (single job)
JOB1=$(sbatch --parsable \
    --job-name=chol_aireml \
    -p mzhang \
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
    --job-name=pheno_aireml \
    -p mzhang \
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

# Step 3: AI-REML (array job, waits for Step 2)
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=aireml \
    -p mzhang \
    --error=$DIR/error/aireml_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --array=1-300%$ARRAY \
    --time=12:00:00 \
    $DIR/AIREML.sh $N $M $S2A $S2D $S2GXG $S2E $MODE $ITERS)
echo "AI-REML job: $JOB3"

# Step 4: Combine results and clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_aireml \
    -p mzhang \
    --error=$DIR/error/combine_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LGXG_FILE && rm -f $LA_FILE && rm -f $LD_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
