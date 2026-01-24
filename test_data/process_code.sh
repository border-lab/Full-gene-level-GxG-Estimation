#####pick a continues 1000 snps  from chr1 chromosome to simulate the high-LD situation
N=$(wc -l < chr1.bim)
START=$(shuf -i 1-$((N-999)) -n 1)
END=$(( START + 999 ))
sed -n "${START},${END}p" chr1.bim | awk '{print $2}' > snp_block1000.txt
plink --bfile chr1 \
      --extract snp_block1000.txt \
      --make-bed \
      --out chr1_block1000

plink --bfile chr1_block1000 --recode A --out chr1_block1000

#####pick a random 1000 snps  from chr1 chromosome to simulate the low-LD situation
plink --bfile chr1 --thin-count 1000 --make-bed --out chr1_sub1000
plink --bfile chr1_sub1000 --recode A --out chr1_sub1000


##### perform LD pruning with a threshold of R2 > 0.8
plink \
  --bfile chr1_sub1000 \
  --indep-pairwise 50 5 0.8 \
  --out snps_to_keep_08

plink \
  --bfile chr1_sub1000 \
  --extract snps_to_keep_08.prune.in \
  --make-bed \
  --out chr1_sub1000_pruned_08

plink \
  --bfile chr1_sub1000_pruned_08 \
  --recode A \
  --out chr1_sub1000_pruned_08

##### perform LD pruning with a threshold of R2 > 0.9
plink \
  --bfile chr1_sub1000 \
  --indep-pairwise 50 5 0.9 \
  --out snps_to_keep_09

plink \
  --bfile chr1_sub1000 \
  --extract snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_sub1000_pruned_09

plink \
  --bfile chr1_sub1000_pruned_09 \
  --recode A \
  --out chr1_sub1000_pruned_09
