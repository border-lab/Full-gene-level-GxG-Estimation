#!/bin/bash
#SBATCH -p dbgap
#SBATCH -J hm3_bed
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o /dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3/%x_%j.out
set -e
PLINK2=/cluster-home/ziyanzha/software/plink2/plink2
GENO_PATH=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3

for CHROM in $(seq 1 22); do
    # 2-bit hard calls (.bed/.bim/.fam), same samples and SNPs as the .pgen
    $PLINK2 \
        --pfile ${GENO_PATH}/ukb_imp_chr${CHROM}_v3\
        --threads ${SLURM_CPUS_PER_TASK}\
        --memory 12000\
        --make-bed \
        --out ${GENO_PATH}/ukb_imp_chr${CHROM}_v3\
        && rm ${GENO_PATH}/ukb_imp_chr${CHROM}_v3.log
done