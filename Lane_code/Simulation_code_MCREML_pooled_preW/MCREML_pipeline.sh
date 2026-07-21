#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline
# (POOLED within-gene pairwise-epistasis, 2-component model).
#   V = s2gxg*W + s2e*I ,
#   W = (1/P) sum_{g=1}^G H_g H_g' = (1/P) sum_g sum_{a<b in g} std(Z_a.Z_b)(...)',
#   P = sum_g C(m_g, 2)  (total within-gene pairs; per-PAIR equal weight).
# The m SNPs are split into G contiguous genes -- equally, or in the per-gene
# proportions given by RATIO below (several Z, one per gene); only
# WITHIN-gene pairs enter W.  This is the "Pooled Model" of generative_model.typ,
# distinct from the equal-weight-per-gene W = (1/G) sum_g K_g of
# Simulation_code_MCREML_gxg_1_over_G (they coincide only for equal gene sizes).
#
# SIMULATION forms W and its Cholesky factor Lgxg to draw a correctly-correlated
# epistasis effect, and CACHES W once (it is deterministic per (genotype, G)).
# ESTIMATION loads the PRE-COMPUTED dense W and applies it as W @ b
# (O(n^2 c) per CG iteration) + conjugate gradient + Hutchinson trace.
# No genotype is needed at estimation time.
#
# ESTIMATE below decouples the two: the phenotype is always simulated from the
# FULL G-gene kernel, but the REML fit can use a kernel pooled over a SUBSET of
# genes (ESTIMATE="1,2" = the front two).  That is a misspecified fit, and the
# recovered s2gxg is attenuated by the subset's share of within-gene pairs.
#
# Mirrors the additive / two_var / three_var pipelines: a 4-step SLURM chain
#   Cholesky (1 job) -> Phenotype (array) -> MC-AI-REML (array) -> combine + cleanup

N=1000
M=1000
G=5
# Per-gene share of the M SNPs.  The genes are CONTIGUOUS and in SNP order, so
# G=5 with RATIO="0.1,0.2,0.3,0.2,0.2" and M=1000 cuts the genotype as
#   gene1 = SNPs 1-100 (10%), gene2 = 101-300 (20%), gene3 = 301-600 (30%),
#   gene4 = 601-800 (20%), gene5 = 801-1000 (20%).
# Must have exactly G positive entries (they are rescaled, so they need not sum
# to 1); leftover SNPs from rounding go to the largest fractional remainders.
# Leave EMPTY (RATIO="") for the old equal split.
RATIO="0.1,0.2,0.3,0.2,0.2"
# Genes used to build the ESTIMATION kernel, 1-based, e.g. ESTIMATE="1,2" fits
# REML with a W pooled over the FRONT TWO genes only.  The phenotype is always
# simulated from ALL G genes, so this is a deliberately MISSPECIFIED fit: it
# asks how much of s2gxg a partial gene panel recovers.  Leave EMPTY
# (ESTIMATE="") to fit with all G genes (the correctly specified model).
#
# Because W_est is normalized by its OWN pair total (tr(W_est) = N, i.e. the
# same Pooled Model simply run on the genes you have), the estimate is
# ATTENUATED: expect s2gxg_hat ~ S2GXG * P_est/P_full, where P_est/P_full is the
# subset's share of within-gene pairs.  With RATIO 0.1,0.2,0.3,0.2,0.2 at M=1000
# the gene pair counts are 4950/19900/44850/19900/19900 (P_full=109500), so
# ESTIMATE="1,2" keeps 24850/109500 ~ 0.227 and S2GXG=0.2 should come back near
# 0.045.  Note gene 1 is 10% of the SNPs but only 4.5% of the PAIRS -- pairs
# grow like m_g^2.  The exact share is written to
# info/split_<FILENAME>.txt; multiply by P_full/P_est to rescale to the full
# panel.
ESTIMATE="1,2"
S2GXG=0.2
S2E=0.8
MODE=RandomSNP
ITERS=30
NMC=100
MEM=2G
ARRAY=20
PARTITION=statgen-gpu


# Scripts and data both live under this single model directory:
DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW

# Gene-split tag: "G5" for the even split, "G5_r0.1-0.2-0.3-0.2-0.2" otherwise.
# W and Lgxg are deterministic per (genotype, split), so the split -- not just G
# -- keys every cached file; this must match gene_split_tag() in Function_MCREML.py.
if [ -z "$RATIO" ]; then
    SPLIT=G${G}
else
    SPLIT=G${G}_r$(echo "$RATIO" | tr ', ' '--' | tr -s '-')
fi

# Estimation-subset suffix: "" for the all-genes fit, "_est1-2" otherwise.  Must
# match gene_subset_tag() in Function_MCREML.py.  It keys the estimation kernel
# W_est and the result/time/info dirs -- but NOT the Cholesky factor or the
# phenotype (TAG), which depend only on the full kernel, so one phenotype set is
# reused across every ESTIMATE subset at the same split.
if [ -z "$ESTIMATE" ]; then
    ESTTAG=
else
    ESTTAG=_est$(echo "$ESTIMATE" | tr ', ' '--' | tr -s '-')
fi

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_${SPLIT}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_${SPLIT}${ESTTAG}
echo "Gene split:       $SPLIT"
echo "Estimation genes: ${ESTIMATE:-all}"

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg and caches the dense pooled GRM W
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
    $DIR/Cholesky.sh $N $M $G $S2GXG $S2E $MODE "$RATIO" "$ESTIMATE")
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
    $DIR/Phenotype.sh $N $M $G $S2GXG $S2E $MODE "$RATIO")
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
    $DIR/MCREML.sh $N $M $G $S2GXG $S2E $MODE "$RATIO" "$ESTIMATE" $ITERS $NMC)
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
echo "Gene split / pair share  -> $DIR/info/split_$FILENAME.txt"
