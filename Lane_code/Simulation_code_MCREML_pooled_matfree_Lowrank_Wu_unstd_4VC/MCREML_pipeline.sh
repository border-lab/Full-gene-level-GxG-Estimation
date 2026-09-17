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
NMC=50
R=30
MEM=2G
MEM_CHOL=64G
MEM_PHENO=24G
ARRAY=50
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC

TAG=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2a${S2A}_s2d${S2D}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}

mkdir -p $DIR/error

LA_FILE=$DIR/Cholesky/La_${TAG}.npy
LD_FILE=$DIR/Cholesky/Ld_${TAG}.npy
LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds La from the dense additive GRM
# K_a = Z_a Z_a'/m, Ld from the dense dominance GRM K_d = Z_d Z_d'/m, and Lgxg
# from the dense pooled UNSTANDARDIZED GRM W (the EXACT kernel, not the
# truncation).  Nothing is cached for estimation: all three GRMs are discarded
# after their factorisations.
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
    $DIR/Cholesky.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1).  Draws y = g_a + g_d + g_gxg
# + e from the three factors, and records each replicate's realized variances
# (epistasis, additive, dominance, residual) next to its phenotype, for Step 3
# to carry into the estimate row.
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
    --array=1-200%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2) -- matrix-free, reads the
# genotype only.  No cached GRM is needed, so this step has no dependency on
# any n-by-n artifact.  W u comes from the rank-$R truncation, K_a u and K_d u
# from a gemm pair against their design (exact), and the score traces
# tr(V^{-1} K_i) from Hutchinson probes ($NMC CG solves per REML iteration).
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
    --cpus-per-task=1 \
    --array=1-200%$ARRAY \
    --time=12:00:00 \
    $DIR/MCREML.sh $N $M $G $S2A $S2D $S2GXG $S2E $MODE $ITERS $NMC $R)
echo "MC-AI-REML job: $JOB3  (4 VC: s2a K_a + s2d K_d + s2gxg W + s2e I; W apply: low-rank, r=$R; trace estimator: hutchinson, Nmc=$NMC)"

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
    --wrap="bash $DIR/combine_code.sh $FILENAME && python3 $DIR/time/average_time.py $DIR/time/rep_times/$FILENAME $DIR/time/result/timing_$FILENAME.txt && rm -rf $DIR/result/$FILENAME && rm -rf $DIR/time/rep_times/$FILENAME && rm -f $LA_FILE && rm -f $LD_FILE && rm -f $LGXG_FILE && rm -rf $PHENO_DIR")
