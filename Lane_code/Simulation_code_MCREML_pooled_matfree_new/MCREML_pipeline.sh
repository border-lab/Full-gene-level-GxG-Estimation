#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline
# (POOLED within-gene pairwise-epistasis, 2-component model).
# MATRIX-FREE ESTIMATION variant of Simulation_code_MCREML_pooled_preW.
#   V = s2gxg*W + s2e*I ,
#   W = (1/P) sum_{g=1}^G H_g H_g' = (1/P) sum_g sum_{a<b in g} std(Z_a.Z_b)(...)',
#   P = sum_g C(m_g, 2)  (total within-gene pairs; per-PAIR equal weight).
# The m SNPs are split into G contiguous genes (several Z, one per gene); only
# WITHIN-gene pairs enter W.  This is the "Pooled Model" of generative_model.typ.
#
# The MODEL, PHENOTYPE and ESTIMATOR are identical to _preW; only the way W is
# APPLIED differs:
#
#   _preW    : builds the dense n-by-n W once, caches it, applies W @ b.
#              O(n^2 c) time and O(n^2) storage per CG iteration.
#   matfree  : never forms W.  Per gene, contracts through the m_g-by-m_g
#              weight matrices -- M_g = Z_g' diag(u) Z_g, then the
#              back-contraction -- for 2 c n sum_g m_g^2 time and
#              O(n m + sum_g m_g^2) storage.  See Wu_complexity.typ.
#
# With G equal-sized genes sum_g m_g^2 = m^2/G, so this variant wins whenever
# m^2/G << n: the large-n regime where a dense W will not fit in memory.
#
# SIMULATION still forms the dense W (a Cholesky factor needs an explicit
# matrix) to draw a correctly-correlated epistasis effect -- but it is NOT
# cached, since estimation rebuilds the action from the genotype.  ESTIMATION
# loads only the n-by-m genotype.
#
# Mirrors the additive / two_var / three_var / preW pipelines: a 4-step SLURM
# chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=1000
M=1000
G=10
S2GXG=0.2
S2E=0.8
MODE=RandomSNP
ITERS=30
NMC=100
# Score-trace estimator.  slq = stochastic Lanczos quadrature: Lanczos on W is
# run once and, because V = s2gxg*W + s2e*I is affine in W, its nodes/weights
# are reused at every variance setting -- so the Nmc probe solves leave the REML
# iteration entirely.  hutchinson = the original Nmc-solves-per-iteration path.
# Identical estimates (agree to the 1e-8 REML stopping tolerance).
TRACE_METHOD=slq
SLQ_K=25
MEM=2G
ARRAY=20
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_new

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg from the dense pooled GRM W.
# Unlike _preW nothing is cached for estimation: W is discarded after the
# factorisation.
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
    $DIR/Cholesky.sh $N $M $G $S2GXG $S2E $MODE)
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
    $DIR/Phenotype.sh $N $M $G $S2GXG $S2E $MODE)
echo "Phenotype job: $JOB2"

# Step 3: MC-AI-REML (array job, waits for Step 2) -- matrix-free, reads the
# genotype only.  No cached W is needed, so this step has no dependency on any
# n-by-n artifact.  The score traces tr(V^{-1} K_i) come from $TRACE_METHOD.
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
    $DIR/MCREML.sh $N $M $G $S2GXG $S2E $MODE $ITERS $NMC $TRACE_METHOD $SLQ_K)
echo "MC-AI-REML job: $JOB3  (trace estimator: $TRACE_METHOD, slq_k=$SLQ_K)"

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
echo "Kernel W build time      -> $DIR/time/result/W_timing_$FILENAME.txt  (simulation only)"
