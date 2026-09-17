#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline -- MATRIX-FREE ESTIMATION variant
# (pairwise-epistasis ONLY, 2-component model).
#   V = s2gxg*W + s2e*I ,   W = (1/p) sum_{a<b} std(Z_a.Z_b)(...)'.
#
# SIMULATION forms the dense W once and its Cholesky factor Lgxg to draw a
# correctly-correlated epistasis effect (same as the dense pipeline -- using the
# O(n^2) W for the one-off simulation is fine).  W is NOT cached.
# ESTIMATION applies W MATRIX-FREE via compute_WU(Z, U, S, R, T) -- the same
# storage-free product used in Simulation_code_MCREML_three_var and the MoM
# onlyW pipeline (O(n m^2 c) per apply, m-by-m weight matrices S,R,T built once).
# No n-by-n W is loaded at estimation time -- only the genotype.
#
# TRADE-OFF vs the pre-computed-dense-W estimator (Simulation_code_MCREML_gxg):
#   + estimation stores only the n-by-m genotype, never the O(n^2) W
#     -> scales to n where a cached dense W no longer fits in memory.
#   - recomputes W U every CG iteration -> ~O(m^2 / n) more arithmetic per mat-vec
#     (the ~2500x slowdown quoted in the dense pipeline at n = m = 1000).
# compute_WU(Z,U,S,R,T) == build_W_batched(Z) @ U to float precision, so results
# match the dense estimator on the same phenotype.
#
# Mirrors the dense / additive pipelines: a 4-step SLURM chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=16000
M=1000
S2GXG=0.2
S2E=0.8
MODE=ContiguousSNP
ITERS=30
NMC=100
MEM=16G
ARRAY=20
PARTITION=mzhang


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_gxg_matfree

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg from the dense W (W not cached)
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml_mf \
    -p $PARTITION \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=16 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $S2GXG $S2E $MODE)
echo "Cholesky job: $JOB1"

# Step 2: Phenotype (array job, waits for Step 1)
JOB2=$(sbatch --parsable \
    --dependency=afterok:$JOB1 \
    --job-name=pheno_mcreml_mf \
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
    $DIR/Phenotype.sh $N $M $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2) -- matrix-free W B in CG
JOB3=$(sbatch --parsable \
    --dependency=afterok:$JOB2 \
    --job-name=mcreml_mf \
    -p $PARTITION \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=16 \
    --array=1-300%$ARRAY \
    --time=12:00:00 \
    $DIR/MCREML.sh $N $M $S2GXG $S2E $MODE $ITERS $NMC)
echo "MC-AI-REML job: $JOB3"

# Step 4: Combine results and clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_mcreml_mf \
    -p $PARTITION \
    --error=$DIR/error/combine_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=00:30:00 \
    --wrap="bash $DIR/combine_code.sh $FILENAME && python3 $DIR/time/average_time.py $DIR/time/rep_times/$FILENAME $DIR/time/result/timing_$FILENAME.txt && rm -rf $DIR/result/$FILENAME && rm -rf $DIR/time/rep_times/$FILENAME && rm -f $LGXG_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
echo "Average estimation time -> $DIR/time/result/timing_$FILENAME.txt"
