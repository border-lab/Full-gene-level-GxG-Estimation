# Test Data

This directory contains test datasets for the gene-level GxG estimation experiments and the code used to generate them.

## Files

### SNP Datasets

| File | Description |
|------|-------------| 
| `chr1_Contiguous1KSNP.zip` | 1,000 contiguous SNPs from chromosome 1 |
| `chr1_Random1KSNP.zip` | 1,000 randomly selected SNPs from chromosome 1 |
| `chr1_Contiguous1KSNP_LDpruned_08.zip` | Contiguous SNPs after LD pruning (r² threshold = 0.8) |
| `chr1_Contiguous1KSNP_LDpruned_09.zip` | Contiguous SNPs after LD pruning (r² threshold = 0.9) |


### Scripts

| File | Description |
|------|-------------|
| `process_code.sh` | Shell script for processing and preparing the test data |

## Usage

Extract the zip files before running the experiments:

```bash
unzip chr1_Contiguous1KSNP.zip
unzip chr1_Random1KSNP.zip
unzip chr1_Contiguous1KSNP_LDpruned_08.zip
unzip chr1_Contiguous1KSNP_LDpruned_09.zip