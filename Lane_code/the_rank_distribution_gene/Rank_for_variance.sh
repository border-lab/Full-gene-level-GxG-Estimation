#!/bin/bash
# Per-gene rank r explaining THRESHOLD of the variance of K_g = Z_g Z_g' -- the
# data for choosing ONE global r for the rank-r truncation.  The stored genotype
# is split into G contiguous genes exactly as the estimator splits it, and each
# gene gets the smallest r whose top-r eigenvalues reach THRESHOLD of tr(K_g).
# GENOTYPE-ONLY: no Cholesky factor, no phenotype, so it is a single standalone
# pre-step run BEFORE the main MCREML pipeline: its output picks the global r
# passed there as R.  Usage: bash Rank_for_variance.sh

N=32000
M=16000
G=160
MODE=ContiguousSNP
THRESHOLD=0.99
MEM=48G                 # the n-by-m CSV is ~4 GB as int64, plus pandas' parse
PARTITION=statgen-gpu

DIR=/home/ziyanzha/MOM_within_gene/the_rank_distribution_gene
GENO=/home/ziyanzha/MOM_within_gene/stored_genotype/${MODE}_n${N}_m${M}.csv

FILENAME=${MODE}_n${N}m${M}_G${G}_thr${THRESHOLD}
OUT=$DIR/result/${FILENAME}.txt

mkdir -p $DIR/error
mkdir -p $DIR/result

JOB=$(sbatch --parsable \
    --job-name=rank_var \
    -p $PARTITION \
    --error=$DIR/error/rank_var_%j.err \
    --output=$DIR/error/rank_var_%j.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --time=04:00:00 \
    --wrap="source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh && \
conda activate tfenv && \
export MKL_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1} && \
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1} && \
cd $DIR && \
python3 $DIR/rank_for_variance.py \
    --genotype $GENO --G $G --threshold $THRESHOLD --out $OUT")

echo "Rank-for-variance job: $JOB  (n=$N, m=$M, G=$G, mode=$MODE, threshold=$THRESHOLD)"
echo "per-gene r (one number per line) -> $OUT"

