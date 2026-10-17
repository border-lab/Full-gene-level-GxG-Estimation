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

## Get_Gene_table.ipynb

Assigns SNPs to genes for the epistasis (within-gene AxA) kernel: for each gene, the SNPs of the genotype set that fall in its exons + TSS ± 1 kb regions. Runs one chromosome at a time.

- **Usage:** set `chrom` in the first cell and run all cells; repeat for each chromosome.
- **Input:** `UKBimp_337K_MAF001_INFO06/ukb_imp_chr${CHROM}_v3.pvar` and `Annotation/gene_regions_exon_tss1kb_gencode_v19.tsv` (from `Get_Gene_region.ipynb`)
- **Method:** for each interval, find the SNPs with `start <= POS <= end` (binary search on the sorted `.pvar`); collect them per gene; keep genes with **≥ 2 SNPs**. A SNP may belong to more than one gene.
- **Output:** `UKBimp_337K_MAF001_INFO06/gxg_gene_table/gene_table_chr${CHROM}.tsv`, long format, one row per gene–SNP pair: `gene_id gene_name gene_type snp_idx`. `snp_idx` is the 0-based row in that chromosome's `.pvar` (= column of its genotype matrix), so each table is valid only for its own `.pvar`.
- **Example (chr22):** 1,135 genes, 20,703 gene–SNP pairs, 18,925 distinct SNPs (15% of 126,417); SNPs per gene: median 12, max 147