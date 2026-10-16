#!/bin/bash
#SBATCH -p dbgap
#SBATCH -J imp_chr22
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -o /dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_MAF001_INFO06/%x_%j.out

# Usage (head node):  sbatch bgen_to_plink.sh
# The output folder must exist BEFORE submitting (Slurm writes the .out there):
#   srun -p dbgap mkdir -p /dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_MAF001_INFO06
# plink2 also prints its log to stdout, so the .out keeps it after the .log is removed.

set -e                 

CHROM=22
PLINK2=/cluster-home/ziyanzha/software/plink2/plink2
OUTPUT_PATH=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_MAF001_INFO06

# MAF>1% INFO>0.6
$PLINK2 \
    --bgen /dbgap/martinzh_ukb16549/data/download_500K/ukb_imp_chr${CHROM}_v3.bgen ref-unknown\
    --sample /dbgap/martinzh_ukb16549/data/download_500K/ukb14048_imp_chr1_v3_s487395.sample\
    --keep /dbgap/martinzh_ukb16549/martinzh/UKBimp_337K_MAF001/unrelated_337K.txt\
    --extract /dbgap/martinzh_ukb16549/martinzh/UKBimp_337K_MAF001/snp_list_chr${CHROM}.MAF_001_INFO_06.txt\
    --rm-dup force-first\
    --maj-ref\
    --geno 0.1\
    --maf 0.01\
    --hwe 1e-50 0\
    --threads 4\
    --memory 28000\
    --make-pgen \
    --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
    && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log

# Recompute the frq info
$PLINK2 \
    --pfile ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
    --freq --out ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3\
        && rm ${OUTPUT_PATH}/ukb_imp_chr${CHROM}_v3.log