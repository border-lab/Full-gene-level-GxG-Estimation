#!/bin/bash
# Variance explained by the leading r eigenvalues of Z Z' -- the diagnostic for
# the rank-r truncation this pipeline's W apply uses.  GENOTYPE-ONLY: it needs
# no Cholesky factor and no phenotype, so it is a single standalone job, not a
# stage of MCREML_pipeline.sh.

N=2000
M=1000
G=10
MODE=ContiguousSNP
R_LIST=1,2,5,10,20,30,40,50,60,80,100,150,200
MEM=8G
PARTITION=statgen-gpu

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC

FILENAME=${MODE}_n${N}m${M}_G${G}

mkdir -p $DIR/error
mkdir -p $DIR/result/variance_explained

JOB=$(sbatch --parsable \
    --job-name=varexp \
    -p $PARTITION \
    --error=$DIR/error/varexp_%j.err \
    --output=$DIR/error/varexp_%j.out \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=4 \
    --time=04:00:00 \
    --wrap="source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh && \
conda activate tfenv && \
export MKL_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1} && \
export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1} && \
python3 $DIR/Simulate_VarianceExplained.py \
    --n $N --m $M --G $G --mode $MODE --r_list $R_LIST")

echo "Variance-explained job: $JOB  (n=$N, m=$M, G=$G, mode=$MODE, r in {$R_LIST})"
echo "rho(r) table -> $DIR/result/variance_explained/${FILENAME}.txt"
