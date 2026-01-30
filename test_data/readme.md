# Test Data

This directory contains test datasets for the gene-level GxG estimation experiments and the code used to generate them.

## Files

### Datasets

| File | Description |
|------|-------------| 
| `chr1_Contiguous1KSNP.zip` | 1,000 contiguous SNPs from chromosome 1 |
| `chr1_Random1KSNP.zip` | 1,000 randomly selected SNPs from chromosome 1 |
| `chr1_Contiguous1KSNP_LDpruned_08.zip` | Contiguous SNPs after LD pruning (r² threshold = 0.8) |
| `chr1_Contiguous1KSNP_LDpruned_09.zip` | Contiguous SNPs after LD pruning (r² threshold = 0.9) |
| `chr1_Contiguous10KSNP_LDpruned_05.zip` | 10K contiguous SNPs, LD pruned (r² = 0.5) |
| `chr1_Contiguous10KSNP_LDpruned_05_maf005.zip` | 10K contiguous SNPs, LD pruned (r² = 0.5), MAF > 0.05 |
| `chr1_block3500_pruned_09.zip` | Block 3500 SNPs, LD pruned (r² = 0.9) |

### Simulation results from cluster

| Folder | Description | LD Pruning | MAF Filter |
|-----------|------|------------|------------|
| `MoM_result_ContiguousSNP_after_pruning_09_SNP_291/` | 291 | r² = 0.9 | - |
| `MoM_result_ContiguousSNP_after_pruning_09_SNP_991/` | 991 | r² = 0.9 | - |
| `MoM_result_ContiguousSNP_after_pruning_05_SNP_861/` | 861 | r² = 0.5 | - |
| `MoM_result_ContiguousSNP_after_pruning_05_SNP_304 _maf_005/` | 304 | r² = 0.5 | > 0.05 |

### SH Scripts

| File | Description |
|------|-------------|
| `process_code.sh` | Shell script for processing and preparing the test data |
| `simulate_code.sh` | Simulation code (in each results directory) |

## How to use

Extract the zip files before running the experiments:

```bash
unzip chr1_Contiguous1KSNP.zip
unzip chr1_Random1KSNP.zip
unzip chr1_Contiguous1KSNP_LDpruned_08.zip
unzip chr1_Contiguous1KSNP_LDpruned_09.zip
unzip chr1_Contiguous10KSNP_LDpruned_05.zip
unzip chr1_Contiguous10KSNP_LDpruned_05_maf005.zip
unzip chr1_block3500_pruned_09.zip