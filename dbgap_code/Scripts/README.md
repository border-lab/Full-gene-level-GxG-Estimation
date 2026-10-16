# UKB imputed genotypes (337K unrelated, MAF > 1%, INFO > 0.6)

## bgen_to_plink.sh

Converts the UK Biobank imputed genotypes (bgen) of one chromosome to plink2 format (`.pgen/.pvar/.psam`).

- **Samples:** Martin's 337K unrelated individuals (`unrelated_337K.txt`)
- **SNPs:** INFO > 0.6 (Martin's SNP list), MAF > 1%, missing rate < 10%, HWE p > 1e-50, duplicates removed
- **Output:** `ukb_imp_chr${CHROM}_v3.{pgen,pvar,psam,afreq}` in this folder

Set `CHROM` (and the job name `-J`) at the top of the script, then submit on the head node with `sbatch bgen_to_plink.sh`.
