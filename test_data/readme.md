# Test Data

This directory contains test datasets for the gene-level GxG estimation experiments and the code used to generate them.

## Files

### Genotype Datasets

**Base datasets:**

| File | Description |
|------|-------------|
| `chr1_Contiguous1KSNP.zip` | 1,000 contiguous SNPs from chromosome 1 (high LD) |
| `chr1_Random1KSNP.zip` | 1,000 randomly selected SNPs from chromosome 1 (low LD) |

**Simulation datasets (2 x 2 x 3 design):**

Datasets are named as `{LD level}_{dataset size}_{variant type}.raw`, varying across three LD pruning thresholds, two dataset sizes, and two variant filters. Dataset Size small mean there are 600 SNPs, and large means there are 1800 SNPs. The data file can be seen in google drive: https://drive.google.com/drive/folders/1Qs0kP-gOBEEs4r6v5XV32ddA4NinnstR?dmr=1&ec=wgc-drive-%5Bmodule%5D-goto

| File | LD Pruning | Dataset Size | Variant Filter |
|------|------------|--------------|----------------|
| `highLD_small_all.raw` | r² = 0.9 | Small | All variants |
| `highLD_small_common.raw` | r² = 0.9 | Small | MAF > 0.05 |
| `highLD_large_all.raw` | r² = 0.9 | Large | All variants |
| `highLD_large_common.raw` | r² = 0.9 | Large | MAF > 0.05 |
| `middleLD_small_all.raw` | r² = 0.5 | Small | All variants |
| `middleLD_small_common.raw` | r² = 0.5 | Small | MAF > 0.05 |
| `middleLD_large_all.raw` | r² = 0.5 | Large | All variants |
| `middleLD_large_common.raw` | r² = 0.5 | Large | MAF > 0.05 |
| `lowLD_small_all.raw` | r² = 0.2 | Small | All variants |
| `lowLD_small_common.raw` | r² = 0.2 | Small | MAF > 0.05 |
| `lowLD_large_all.raw` | r² = 0.2 | Large | All variants |
| `lowLD_large_common.raw` | r² = 0.2 | Large | MAF > 0.05 |

### Simulation Results

The `All_simulation_result/` directory contains MoM estimation results from cluster runs, organized by LD level, dataset size, and variant type:

| Folder | LD Pruning | Dataset Size | Variant Filter | SNPs (m) |
|--------|------------|--------------|----------------|----------|
| `highLD_largeDataset_allVariants/` | r² = 0.9 | Large | All | 1800 |
| `highLD_largeDataset_commonVariants/` | r² = 0.9 | Large | MAF > 0.05 | 1800 |
| `highLD_smallDataset_allVariants/` | r² = 0.9 | Small | All | 600 |
| `highLD_smallDataset_commonVariants/` | r² = 0.9 | Small | MAF > 0.05 | 600 |
| `middleLD_largeDataset_allVariants/` | r² = 0.5 | Large | All | 1800 |
| `middleLD_largeDataset_commonVariants/` | r² = 0.5 | Large | MAF > 0.05 | 1800 |
| `middleLD_smallDataset_allVariants/` | r² = 0.5 | Small | All | 600 |
| `middleLD_smallDataset_comonVariants/` | r² = 0.5 | Small | MAF > 0.05 | 600 |
| `lowLD_largeDataset_allVariants/` | r² = 0.2 | Large | All | 1800 |
| `lowLD_largeDataset_commonVariants/` | r² = 0.2 | Large | MAF > 0.05 | 1800 |
| `lowLD_smallDataset_allVariants/` | r² = 0.2 | Small | All | 600 |
| `lowLD_smallDataset_commonVariants/` | r² = 0.2 | Small | MAF > 0.05 | 600 |

Each folder contains result files named `{N}m{M}.txt`, where `N` is the sample size (1000, 2000, 4000, or 8000) and `M` is the number of SNPs. Each line contains a tuple of two estimated variance components.

### Scripts

| File | Description |
|------|-------------|
| `process_code.sh` | Shell script for extracting SNP blocks, LD pruning, and MAF filtering using PLINK |

## How to use

Extract the zip files before running the experiments:

```bash
# Base datasets
unzip chr1_Contiguous1KSNP.zip
unzip chr1_Random1KSNP.zip

# Simulation datasets

Google Drive: https://drive.google.com/drive/folders/1Qs0kP-gOBEEs4r6v5XV32ddA4NinnstR?dmr=1&ec=wgc-drive-%5Bmodule%5D-goto
---