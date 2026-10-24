#!/bin/bash
#SBATCH -p dbgap
#SBATCH -J subset_100K
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH -o /dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_100K_HM3/logs/chr%a_%A.log

# Usage (head node):  sbatch subset_100K.sh
# Keep a fixed random 100K samples (keep_100K_seed1.txt, made in the notebook)
# and all HM3 SNPs.  One array task per chromosome.
SRC=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3
OUT=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_100K_HM3
PLINK2=/cluster-home/ziyanzha/software/plink2/plink2
c=$SLURM_ARRAY_TASK_ID

$PLINK2 --pfile $SRC/ukb_imp_chr${c}_v3 \
        --keep $OUT/keep_100K_seed1.txt \
        --threads 4 --memory 12000 \
        --make-pgen --out $OUT/ukb_imp_chr${c}_v3