#####pick a random 1000 snps  from chr1 chromosome to simulate the low-LD situation
plink --bfile chr1 --thin-count 1000 --make-bed --out chr1_sub1000
plink --bfile chr1_sub1000 --recode A --out chr1_sub1000

#####pick the whole dataset used for simulation
N=$(wc -l < chr1.bim)
START=$(shuf -i 1-$((N-99999)) -n 1)
END=$(( START + 99999 ))
sed -n "${START},${END}p" chr1.bim | awk '{print $2}' > snp_block10k.txt
plink --bfile chr1 \
      --extract snp_block10k.txt \
      --make-bed \
      --out chr1_block10k
plink --bfile chr1_block10k --recode A --out chr1_block10k


####filter out rare varients (MAF > 0.05)
plink --bfile chr1_block10k --maf 0.05 --make-bed --out chr1_block10k_maf005

##### perform LD pruning with a threshold of R2 > 0.9 (common variant)
plink \
  --bfile chr1_block10k_maf005 \
  --indep-pairwise 1000 5 0.9 \
  --out snps_to_keep_09_maf005

plink \
  --bfile chr1_block10k_maf005 \
  --extract snps_to_keep_09_maf005.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_09_maf005

plink \
  --bfile chr1_block10k_pruned_09_maf005 \
  --recode A \
  --out chr1_block10k_pruned_09_maf005

##### perform LD pruning with a threshold of R2 > 0.9 (All variant)

plink \
  --bfile chr1_block10k \
  --indep-pairwise 1000 5 0.9 \
  --out snps_to_keep_09

plink \
  --bfile chr1_block10k \
  --extract snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_09

plink \
  --bfile chr1_block10k_pruned_09 \
  --recode A \
  --out chr1_block10k_pruned_09




##### perform LD pruning with a threshold of R2 > 0.5 (common variant)
plink \
  --bfile chr1_block10k_maf005 \
  --indep-pairwise 1000 5 0.5 \
  --out snps_to_keep_05_maf005

plink \
  --bfile chr1_block10k_maf005 \
  --extract snps_to_keep_05_maf005.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_05_maf005

plink \
  --bfile chr1_block10k_pruned_05_maf005 \
  --recode A \
  --out chr1_block10k_pruned_05_maf005



##### perform LD pruning with a threshold of R2 > 0.5 (All variant)

plink \
  --bfile chr1_block10k \
  --indep-pairwise 1000 5 0.5 \
  --out snps_to_keep_05

plink \
  --bfile chr1_block10k \
  --extract snps_to_keep_05.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_05

plink \
  --bfile chr1_block10k_pruned_05 \
  --recode A \
  --out chr1_block10k_pruned_05



##### perform LD pruning with a threshold of R2 > 0.2 (common variant)
plink \
  --bfile chr1_block10k_maf005 \
  --indep-pairwise 1000 5 0.2 \
  --out snps_to_keep_02_maf005

plink \
  --bfile chr1_block10k_maf005 \
  --extract snps_to_keep_02_maf005.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_02_maf005

plink \
  --bfile chr1_block10k_pruned_02_maf005 \
  --recode A \
  --out chr1_block10k_pruned_02_maf005



##### perform LD pruning with a threshold of R2 > 0.2 (All variant)
plink \
  --bfile chr1_block10k \
  --indep-pairwise 1000 5 0.2 \
  --out snps_to_keep_02

plink \
  --bfile chr1_block10k \
  --extract snps_to_keep_02.prune.in \
  --make-bed \
  --out chr1_block10k_pruned_02

plink \
  --bfile chr1_block10k_pruned_02 \
  --recode A \
  --out chr1_block10k_pruned_02



#### choose the dataset with R2 > 0.9 for common variant and all variant to simulate the high-LD situation
cut -d' ' -f1-606 chr1_block10k_pruned_09.raw > highLD_small_all.raw




