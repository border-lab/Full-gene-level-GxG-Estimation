#!/bin/bash
# Monte-Carlo AI-REML simulation pipeline -- POOLED within-gene epistasis,
# pre-computed dense W (2-component model), EQUAL WEIGHT PER GENE.
#   V = s2gxg*W + s2e*I ,   W = (1/G) sum_{g=1}^G K_g   (pooled within-gene GRM),
#   K_g = (1/p_g) sum_{a<b in g} std(Z_a.Z_b)(...)'  (gene g's own GRM, tr=N).
# The m SNPs are split into G contiguous genes -- equally, or in the per-gene
# proportions given by RATIO below; only WITHIN-gene SNP pairs enter W.  Every
# gene enters with the SAME 1/G weight regardless of size ("The Model in Notes"),
# distinct from the per-PAIR W = (1/P) sum_g H_g H_g' of
# Simulation_code_MCREML_pooled_preW (they coincide only for equal gene sizes).
#
# SIMULATION forms the pooled W and its Cholesky factor Lgxg to draw a correctly-
# correlated epistasis effect, and CACHES W once (deterministic per genotype, G,
# split).  ESTIMATION loads the PRE-COMPUTED dense W and applies it as W @ b +
# CG + Hutchinson trace -- only the kernel W differs.
#
# ESTIMATE below decouples the two: the phenotype is always simulated from the
# FULL G-gene kernel, but the REML fit can use a kernel pooled over a SUBSET of
# genes (ESTIMATE="1,2" = the front two).  That is a misspecified fit, and the
# recovered s2gxg is attenuated by the subset's SHARE OF GENES.
#
# 4-step SLURM chain:
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
# Leave EMPTY (RATIO="") for the old equal split.  NOTE: under this 1/G kernel
# the ratio only changes WHICH SNPs are in each gene -- every gene still weighs
# the same, so it does not change the attenuation target below.
RATIO="0.1,0.2,0.3,0.2,0.2"
# Genes used to build the ESTIMATION kernel, 1-based, e.g. ESTIMATE="1,2" fits
# REML with a W pooled over the FRONT TWO genes only.  The phenotype is always
# simulated from ALL G genes, so this is a deliberately MISSPECIFIED fit: it
# asks how much of s2gxg a partial gene panel recovers.  Leave EMPTY
# (ESTIMATE="") to fit with all G genes (the correctly specified model).
#
# Listing every gene explicitly (ESTIMATE="1,2,3,4,5" at G=5) is also a valid
# correctly-specified run -- W_est then equals W_full and gene_share is 1 -- but
# it is NOT the same as ESTIMATE="": it tags every output _est1-2-3-4-5, so use
# it when you want the all-genes baseline filed alongside the subset runs.
#
# Because W_est averages its OWN genes equally (tr(W_est) = N, i.e. the same 1/G
# model simply run on the genes you have), the estimate is ATTENUATED: expect
# s2gxg_hat ~ S2GXG * G_est/G_full, the subset's SHARE OF GENES.  Unlike the
# per-PAIR preW model, this share is COUNT-based, so it ignores gene sizes: at
# G=5 with ESTIMATE="1,2" it is 2/5 = 0.4 whatever the RATIO, and S2GXG=0.2
# should come back near 0.08.  The exact share is written to
# info/split_<FILENAME>.txt; multiply by G_full/G_est to rescale to the full
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
DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_perG

# Gene-SIZE tag: "" for the even split, "_r0.1-0.2-0.3-0.2-0.2" otherwise.  W and
# Lgxg are deterministic per (genotype, split), so the split -- not just G --
# keys every cached file; this must match ratio_tag() in Function_MCREML.py.
# (The gene COUNT G already sits in the base tag, so only the ratio is added.)
if [ -z "$RATIO" ]; then
    RTAG=
else
    RTAG=_r$(echo "$RATIO" | tr ', ' '--' | tr -s '-')
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

TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_G${G}_n${N}_m${M}${RTAG}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_G${G}_n${N}m${M}${RTAG}${ESTTAG}
echo "Gene split:       G${G}${RTAG}"
echo "Estimation genes: ${ESTIMATE:-all}"

mkdir -p $DIR/error

LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy
PHENO_DIR=$DIR/Phenotype/y_${TAG}

# Step 1: Cholesky (single job) -- builds Lgxg and caches the pooled dense GRM W
JOB1=$(sbatch --parsable \
    --job-name=chol_pooled \
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
    --job-name=pheno_pooled \
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
    --job-name=mcreml_pooled \
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

# Step 4: Combine results + average timing, then clean up (waits for Step 3)
JOB4=$(sbatch --parsable \
    --dependency=afterok:$JOB3 \
    --job-name=combine_pooled \
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
echo "Gene split / gene share  -> $DIR/info/split_$FILENAME.txt"
