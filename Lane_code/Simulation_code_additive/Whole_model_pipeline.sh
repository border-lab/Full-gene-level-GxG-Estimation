#!/bin/bash
#SBATCH --job-name=grm_filter
#SBATCH -p mzhang
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=grm_filter_%j.log
#SBATCH --error=grm_filter_%j.err

# ============================================================
# GRM relatedness filtering on HPC (mzhang partition)
# GRM already generated as: chr1_4000snp_grm
# Filters relatedness > 0.05, extracts unrelated genotype.
# ============================================================

set -e

# ---- CONFIG (edit if needed) ----
WORKDIR=/home/ziyanzha/MOM_within_gene/genotype_matrix_plink
DATA=chr1_4000snp
GRM=${DATA}_grm
CUTOFF=0.05
# ---------------------------------

cd $WORKDIR
echo "=== Job started: $(date) ==="
echo "GRM: $GRM, cutoff: $CUTOFF"

# Check GRM files
ls -lh ${GRM}.grm.id ${GRM}.grm.bin ${GRM}.grm.N.bin
N_TOTAL=$(wc -l < ${GRM}.grm.id)
echo "Total individuals: $N_TOTAL"
echo ""

# Step 1: GCTA grm-cutoff (greedy independent set)
echo "=== Step 1: GCTA grm-cutoff filtering ==="
gcta64 --grm $GRM --grm-cutoff $CUTOFF --make-grm --out ${DATA}_unrelated --thread-num 16

N_KEEP=$(wc -l < ${DATA}_unrelated.grm.id)
echo "Remaining after filtering: $N_KEEP / $N_TOTAL"
echo ""

# Step 2: extract unrelated genotype
echo "=== Step 2: extract unrelated genotype ==="
plink --bfile $DATA \
      --keep ${DATA}_unrelated.grm.id \
      --make-bed \
      --allow-no-sex \
      --out ${DATA}_unrelated_geno
echo "Final: $(wc -l < ${DATA}_unrelated_geno.fam) individuals, $(wc -l < ${DATA}_unrelated_geno.bim) SNPs"
echo ""

# Step 3: convert to CSV for MoM pipeline
echo "=== Step 3: convert to CSV ==="
plink --bfile ${DATA}_unrelated_geno --recode A --allow-no-sex --out tmp_unrel
tail -n +2 tmp_unrel.raw | cut -d' ' -f7- | tr ' ' ',' > ${DATA}_unrelated_geno.csv
rm -f tmp_unrel.raw tmp_unrel.log tmp_unrel.nosex
echo "CSV: $(wc -l < ${DATA}_unrelated_geno.csv) individuals"

echo "=== Job finished: $(date) ==="