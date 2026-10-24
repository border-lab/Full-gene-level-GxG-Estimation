#!/bin/bash
#SBATCH -p dbgap
#SBATCH -J hm3_subset
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o /dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3/%x_%j.out
set -e
PLINK2=/cluster-home/ziyanzha/software/plink2/plink2
INPUT_PATH=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_MAF001_INFO06
OUTPUT_PATH=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3
LIST_PATH=/dbgap/martinzh_ukb16549/ziyanzha/Annotation/hm3_extract

for CHROM in $(seq 1 22); do
    # HM3 SNPs only (same samples, same QC as the input)
    $PLINK2 \
        --pfile ${INPUT_PATH}/ukb_imp_chr${CHROM}_v3\
        --extract ${LIST_PATH}/hm3_chr${CHROM}.txt\
        --threads ${SLURM_CPUS_PER_TASK}\
        --memory 12000\
        --make-pgen \
        --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log
    # Recompute the frq info
    $PLINK2 \
        --pfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
        --freq --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
            && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log
done