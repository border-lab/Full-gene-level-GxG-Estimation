#!/bin/bash
# Pooled within-gene AxA AI-REML pipeline (3-component: additive + pooled GxG + noise).
# Same 4-step SLURM dependency chain as the single-gene AIREML pipeline:
#   Cholesky (1 job) -> Phenotype (array) -> AI-REML (array) -> combine + cleanup
# The epistatic component is the pooled kernel W = (1/G) sum_g K_g over G
# subsampled contiguous regions (region_size markers each). The (G, region_size,
# seed) triple fixes which regions are drawn; W is built ONCE in the Cholesky
# step, stored, and reused by the AI-REML step (no re-draw at fit time).
# Self-contained model directory (own Cholesky/Phenotype/W/result).

N=4000
M=1000
S2A=0.2
S2GXG=0.7
S2E=0.1
MODE=Random
ITERS=30
MEM=16G
ARRAY=20

# Pooled-kernel knobs:
G=50            # number of subsampled contiguous regions ("genes")
REGION_SIZE=20  # markers per region (m_g)
SEED=42         # fixes which regions are drawn (region placement)

# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/AIREML_pooled

TAG=G${G}_rs${REGION_SIZE}_seed${SEED}
FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${TAG}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky_Lgxg/Lgxg_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${TAG}.npy
LA_FILE=$DIR/Cholesky_La/La_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${TAG}.npy
W_FILE=$DIR/W/W_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${TAG}

# Step 1: Cholesky (single job) -- builds & stores the pooled kernel W once
JOB1=$(sbatch --parsable \
    --job-name=chol_aireml_pooled \
    -p mzhang \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=16 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $S2A $S2GXG $S2E $MODE $G $REGION_SIZE $SEED)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1)
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    --job-name=pheno_aireml_pooled \
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
    $DIR/Phenotype.sh $N $M $S2A $S2GXG $S2E $MODE $G $REGION_SIZE $SEED)
echo "Phenotype job: $JOB2"

# Step 3: AI-REML (array job, waits for Step 2)
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=aireml_pooled \
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
    $DIR/AIREML.sh $N $M $S2A $S2GXG $S2E $MODE $ITERS $G $REGION_SIZE $SEED)
echo "AI-REML job: $JOB3"

# Step 4: Combine results and clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_aireml_pooled \
    -p mzhang \
    --error=$DIR/error/combine_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LGXG_FILE && rm -f $LA_FILE && rm -f $W_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
