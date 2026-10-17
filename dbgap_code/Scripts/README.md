# UKB imputed genotypes (337K unrelated, MAF > 1%, INFO > 0.6)

## bgen_to_plink.sh

Converts the UK Biobank imputed genotypes (bgen) of one chromosome to plink2 format (`.pgen/.pvar/.psam`).

- **Samples:** Martin's 337K unrelated individuals (`unrelated_337K.txt`)
- **SNPs:** INFO > 0.6 (Martin's SNP list), MAF > 1%, missing rate < 10%, HWE p > 1e-50, duplicates removed
- **Output:** `ukb_imp_chr${CHROM}_v3.{pgen,pvar,psam,afreq}` in this folder

Set `CHROM` (and the job name `-J`) at the top of the script, then submit on the head node with `sbatch bgen_to_plink.sh`.

## Get_Gene_region.ipynb

Defines the gene regions used by the epistasis (within-gene AxA) kernel: for every gene, its **exons plus TSS ± 1 kb**.

- **Input:** GENCODE v19 annotation (GRCh37), `gencode.v19.annotation.gtf.gz`
- **Genes:** all genes on chr1–22 (all gene types, all transcripts); `gene_type` is kept for later filtering
- **Region of a gene:** union of the exons of all its transcripts, plus ±1 kb around the TSS of every transcript (TSS = start on the + strand, end on the − strand); overlapping intervals within a gene are merged. Regions of different genes may overlap.
- **Output:** `Annotation/gene_regions_exon_tss1kb_gencode_v19.tsv`, one row per merged interval: `chrom start end gene_id gene_name gene_type` (chrom without `chr`; 1-based, inclusive, as in `.pvar`)
- **Size:** 260,655 intervals, about 276 Mb (9.6% of the autosomes; overlaps between genes counted once per gene)
