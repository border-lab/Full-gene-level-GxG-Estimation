#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline
# (pairwise-epistasis ONLY, 2-component model).
# MEAN-CENTRED-KERNEL variant of MCREML_gxg.
#   V = s2gxg*W + s2e*I ,   W = (1/p) sum_{a<b} h_ab h_ab' ,
#   h_ab = Z_a.Z_b - mean(Z_a.Z_b)   (mean-centred, no 1/sigma_ab scaling),
#   p = m(m-1)/2.
# SIMULATION forms W and its Cholesky factor Lgxg to draw a correctly-correlated
# epistasis effect, and CACHES W once (it is deterministic per genotype).
# ESTIMATION loads the PRE-COMPUTED dense W and applies it as W @ b
# (O(n^2 c) per CG iteration) + conjugate gradient + Hutchinson trace.
# No genotype is needed at estimation time.  See MCREML_gxg.typ.
#
# THE ONLY DIFFERENCE FROM MCREML_gxg IS THE KERNEL.  The interaction column
# h_ab keeps its mean-centring but drops the 1/sigma_ab scaling on BOTH sides
# -- the same W is cached, factorized and estimated with -- so a run here
# differs from an MCREML_gxg run in the kernel alone.  The SNPs are still
# column-standardized; only the pair column changes.  Of the standardized
# kernel's identities, tr(W) = n is given up but W 1 = 0 still holds (the
# columns stay centred).  W stays PSD, so the Cholesky draw is unaffected, and
# the spectrum spreads out, which moves W away from I and helps
# identifiability.
#
# CACHE COLLISION -- the one thing to watch.  MCREML_gxg caches its STANDARDIZED
# W as stored_genotype/W_<mode>_n<n>_m<m>.npy, keyed by genotype only with no
# kernel in the name, and writes it only `if not os.path.exists`.  This pipeline
# uses W_ctr_<mode>_n<n>_m<m>.npy for exactly that reason: sharing the name
# would let whichever pipeline ran first win, and the other would silently draw
# phenotypes from one kernel and estimate with the other.  Stale W_raw_*.npy
# caches belong to the old RAW (un-centred) kernel and must not be reused.
#
# NOTE ON COST: reusing the cached dense W avoids the O(n m^2) matrix-free
# rebuild inside every CG iteration (~2500x faster at n = m = 1000).  The BUILD
# is also much cheaper here than in MCREML_gxg: dropping the 1/sigma_ab scaling
# lets W come from the Hadamard identity (K.*K - D D')/(2p), double-centred, in
# O(n^2 m) instead of the O(n^2 m^2) pair loop -- the mean-centring factors
# through the Gram sum as W <- M W M with M = I - 11'/n, an O(n^2) step -- and
# the gap grows ~linearly in m, so it is the Cholesky job's runtime that
# changes most.
#
# Mirrors the additive / two_var / three_var pipelines: a 4-step SLURM chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=8000
M=1000
S2GXG=0.2
S2E=0.8
MODE=RandomSNP
ITERS=30
NMC=100
MEM=16G
ARRAY=20
PARTITION=statgen-gpu



# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_gxg_different_kernal

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg and caches the dense mean-centred GRM W
JOB1=$(sbatch --parsable \
    --job-name=chol_mcreml \
    -p $PARTITION \
    --error=$DIR/error/cholesky_%j.err \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=48:00:00 \
    $DIR/Cholesky.sh $N $M $S2GXG $S2E $MODE)
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
    --cpus-per-task=1 \
    --array=1-150%$ARRAY \
    --time=12:00:00 \
    $DIR/Phenotype.sh $N $M $S2GXG $S2E $MODE)
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
    --cpus-per-task=1 \
    --array=1-150%$ARRAY \
    --time=12:00:00 \
    $DIR/MCREML.sh $N $M $S2GXG $S2E $MODE $ITERS $NMC)
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
    --wrap="bash $DIR/combine_code.sh $FILENAME && python3 $DIR/time/average_time.py $DIR/time/rep_times/$FILENAME $DIR/time/result/timing_$FILENAME.txt && rm -rf $DIR/result/$FILENAME && rm -rf $DIR/time/rep_times/$FILENAME && rm -f $LGXG_FILE && rm -rf $PHENO_DIR")
echo "Combine job: $JOB4"
echo "Average estimation time -> $DIR/time/result/timing_$FILENAME.txt"
echo "Kernel W build time      -> $DIR/time/result/W_timing_$FILENAME.txt"
