#!/bin/bash
# Sparse simulation pipeline

N=4000
M=1000
S2A=0.2
S2GXG=0.7
S2E=0.1
MODE=Random
DENSITY= 1
SPARSE=SNP
NMC=100
MEM=6G
ARRAY=20

DIR=/home/ziyanzha/MOM_within_gene/Sparse_model
FILENAME=${MODE}_n${N}m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_${SPARSE}density${DENSITY}
LGXG_FILE=$DIR/Cholesky_Lgxg/${SPARSE}Lgxg_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_density${DENSITY}.npy
LA_FILE=$DIR/Cholesky_La/${SPARSE}La_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_density${DENSITY}.npy
PHENO_DIR=$DIR/Phenotype/${SPARSE}y_${MODE}_n${N}_m${M}_s2a${S2A}_s2gxg${S2GXG}_s2e${S2E}_density${DENSITY}

# Step 1: Cholesky (single job)
JOB1=$(sbatch --parsable \
    --job-name=cholesky \
    -p mzhang \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=16 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $S2A $S2GXG $S2E $MODE $DENSITY $SPARSE)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1)
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    --job-name=phenotype \
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
    $DIR/Phenotype.sh $N $M $S2A $S2GXG $S2E $MODE $DENSITY $SPARSE)
echo "Phenotype job: $JOB2"

# Step 3: MoM (array job, waits for Step 2)
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=mom \
    -p mzhang \
    --error=$DIR/error/mom_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --array=1-300%$ARRAY \
    --time=12:00:00 \
    $DIR/MoM.sh $N $M $S2A $S2GXG $S2E $MODE $DENSITY $SPARSE $NMC)
echo "MoM job: $JOB3"

# Step 4: Combine results and clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine \
    -p mzhang \
    --error=$DIR/error/combine_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME && rm -rf $DIR/result/$FILENAME && rm -f $LGXG_FILE && rm -f $LA_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"