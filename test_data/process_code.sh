#####pick a continues 1000 snps  from chr1 chromosome to simulate the high-LD situation
N=$(wc -l < chr1.bim)
START=$(shuf -i 1-$((N-39999)) -n 1)
END=$(( START + 39999 ))
sed -n "${START},${END}p" chr1.bim | awk '{print $2}' > snp_block400000.txt
plink --bfile chr1 \
      --extract snp_block400000.txt \
      --make-bed \
      --out chr1_block400000

plink --bfile chr1_block400000 --recode A --out chr1_block400000

#####pick a random 1000 snps  from chr1 chromosome to simulate the low-LD situation
plink --bfile chr1 --thin-count 1000 --make-bed --out chr1_sub1000
plink --bfile chr1_sub1000 --recode A --out chr1_sub1000


##### perform LD pruning with a threshold of R2 > 0.8
plink \
  --bfile chr1_block1000 \
  --indep-pairwise 1000 5 0.8 \
  --out snps_to_keep_08

plink \
  --bfile chr1_block1000 \
  --extract snps_to_keep_08.prune.in \
  --make-bed \
  --out chr1_block1000_pruned_08

plink \
  --bfile chr1_block1000_pruned_08 \
  --recode A \
  --out chr1_block1000_pruned_08

##### perform LD pruning with a threshold of R2 > 0.9
plink \
  --bfile chr1_block1000 \
  --indep-pairwise 1000 5 0.9 \
  --out snps_to_keep_09

plink \
  --bfile chr1_block1000 \
  --extract snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block1000_pruned_09

plink \
  --bfile chr1_block1000_pruned_09 \
  --recode A \
  --out chr1_block1000_pruned_09


 ##### perform LD pruning for 3500 SNP with a threshold of R2 > 0.9
plink \
  --bfile chr1_block3500 \
  --indep-pairwise 3500 5 0.9 \
  --out 3500snps_to_keep_09
plink \
  --bfile chr1_block3500 \
  --extract 3500snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block3500_pruned_09

plink \
  --bfile chr1_block3500_pruned_09 \
  --recode A \
  --out chr1_block3500_pruned_09


  ##### perform LD pruning for 10K SNP with a threshold of R2 > 0.5
plink \
  --bfile chr1_block10000 \
  --indep-pairwise 10000 5 0.5 \
  --out 10ksnps_to_keep_05

plink \
  --bfile chr1_block10000 \
  --extract 10ksnps_to_keep_05.prune.in \
  --make-bed \
  --out chr1_block10000_pruned_05

plink \
  --bfile chr1_block10000_pruned_05 \
  --recode A \
  --out chr1_block10000_pruned_05

##### filter SNPs with MAF < 0.05
plink --bfile chr1_block10000_pruned_05 --maf 0.05 --make-bed --out chr1_block10000_pruned_05_maf005

plink \
  --bfile chr1_block10000_pruned_05_maf005 \
  --recode A \
  --out chr1_block10000_pruned_05_maf005


###### 2 by 2 by 3

N=$(wc -l < chr1.bim)
START=$(shuf -i 1-$((N-9999)) -n 1)
END=$(( START + 9999 ))
sed -n "${START},${END}p" chr1.bim | awk '{print $2}' > snp_block10000.txt
plink --bfile chr1 \
      --extract snp_block10000.txt \
      --make-bed \
      --out chr1_block10000

plink --bfile chr1_block10000 --recode A --out chr1_block10000


#######for LD > 0.9 pruning test
 ##### R2 > 0.9， small sample size, all maf.
plink \
  --bfile chr1_block1500 \
  --indep-pairwise 1500 5 0.9 \
  --out 1500snps_to_keep_09
plink \
  --bfile chr1_block1500 \
  --extract 1500snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block1500_pruned_09

plink \
  --bfile chr1_block1500_pruned_09 \
  --recode A \
  --out chr1_block1500_pruned_09


 ##### R2 > 0.9， large sample size, all maf.

 plink \
  --bfile chr1_block5000 \
  --indep-pairwise 5000 5 0.9 \
  --out 5000snps_to_keep_09
plink \
  --bfile chr1_block5000 \
  --extract 5000snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block5000_pruned_09

plink \
  --bfile chr1_block5000_pruned_09 \
  --recode A \
  --out chr1_block5000_pruned_09


#####  ##### R2 > 0.9， small sample size, common variant.
plink --bfile chr1_block5000_pruned_09 --maf 0.05 --make-bed --out chr1_block5000_pruned_09_maf005

plink \
  --bfile chr1_block5000_pruned_09_maf005 \
  --recode A \
  --out chr1_block5000_pruned_09_maf005


#####  ##### R2 > 0.9， large sample size, common variant.
plink \
  --bfile chr1_block10000 \
  --indep-pairwise 10000 5 0.9 \
  --out 10000snps_to_keep_09
plink \
  --bfile chr1_block10000 \
  --extract 10000snps_to_keep_09.prune.in \
  --make-bed \
  --out chr1_block10000_pruned_09

plink \
  --bfile chr1_block10000_pruned_09 \
  --recode A \
  --out chr1_block10000_pruned_09


plink --bfile chr1_block10000_pruned_09 --maf 0.05 --make-bed --out chr1_block10000_pruned_09_maf005

plink \
  --bfile chr1_block10000_pruned_09_maf005 \
  --recode A \
  --out chr1_block10000_pruned_09_maf005


#######for LD > 0.5 pruning test

 ##### R2 > 0.5， small sample size, all variant.
plink \
  --bfile chr1_block4000 \
  --indep-pairwise 4000 5 0.5 \
  --out 4000snps_to_keep_05
plink \
  --bfile chr1_block4000 \
  --extract 4000snps_to_keep_05.prune.in \
  --make-bed \
  --out chr1_block4000_pruned_05

plink \
  --bfile chr1_block4000_pruned_05 \
  --recode A \
  --out chr1_block4000_pruned_05

##### R2 > 0.5， large sample size, all variant.
plink \
  --bfile chr1_block15000 \
  --indep-pairwise 15000 5 0.5 \
  --out 15000snps_to_keep_05

plink \
  --bfile chr1_block15000 \
  --extract 15000snps_to_keep_05.prune.in \
  --make-bed \
  --out chr1_block15000_pruned_05

plink \
  --bfile chr1_block15000_pruned_05 \
  --recode A \
  --out chr1_block15000_pruned_05

##### R2 > 0.5， small sample size, common variant.

plink --bfile chr1_block15000_pruned_05 --maf 0.05 --make-bed --out chr1_block15000_pruned_05_maf005

plink \
  --bfile chr1_block15000_pruned_05_maf005 \
  --recode A \
  --out chr1_block15000_pruned_05_maf005

##### R2 > 0.5， large sample size, common variant.

plink \
  --bfile chr1_block20000 \
  --indep-pairwise 20000 5 0.5 \
  --out 20000snps_to_keep_05

plink \
  --bfile chr1_block20000 \
  --extract 20000snps_to_keep_05.prune.in \
  --make-bed \
  --out chr1_block20000_pruned_05

plink \
  --bfile chr1_block20000_pruned_05 \
  --recode A \
  --out chr1_block20000_pruned_05


plink --bfile chr1_block20000_pruned_05 --maf 0.05 --make-bed --out chr1_block20000_pruned_05_maf005

plink \
  --bfile chr1_block20000_pruned_05_maf005 \
  --recode A \
  --out chr1_block20000_pruned_05_maf005

  #######for LD > 0.2 pruning test

   ##### R2 > 0.2， small sample size, all variant.
plink \
  --bfile chr1_block10000 \
  --indep-pairwise 10000 5 0.2 \
  --out 10000snps_to_keep_02

plink \
  --bfile chr1_block10000 \
  --extract 10000snps_to_keep_02.prune.in \
  --make-bed \
  --out chr1_block10000_pruned_02

plink \
  --bfile chr1_block10000_pruned_02 \
  --recode A \
  --out chr1_block10000_pruned_02


     ##### R2 > 0.2， large sample size, all variant.
plink \
  --bfile chr1_block20000 \
  --indep-pairwise 20000 5 0.2 \
  --out 20000snps_to_keep_02
  
plink \
  --bfile chr1_block20000 \
  --extract 20000snps_to_keep_02.prune.in \
  --make-bed \
  --out chr1_block20000_pruned_02

plink \
  --bfile chr1_block20000_pruned_02 \
  --recode A \
  --out chr1_block20000_pruned_02


##### R2 > 0.2， small sample size, common variant.
plink --bfile chr1_block20000_pruned_02 --maf 0.05 --make-bed --out chr1_block20000_pruned_02_maf005

plink \
  --bfile chr1_block20000_pruned_02_maf005 \
  --recode A \
  --out chr1_block20000_pruned_02_maf005


  ##### R2 > 0.2， large sample size, common variant.

plink --bfile chr1_first30000 --recode A --out chr1_first30000
plink --bfile chr1_last30000 --recode A --out chr1_last30000



plink \
  --bfile chr1_last30000 \
  --indep-pairwise 30000 5 0.2 \
  --out last_30000snps_to_keep_02

plink \
  --bfile chr1_last30000 \
  --extract last_30000snps_to_keep_02.prune.in \
  --make-bed \
  --out chr1_last30000_pruned_02

plink \
  --bfile chr1_last30000_pruned_02 \
  --recode A \
  --out chr1_last30000_pruned_02


plink --bfile chr1_last30000_pruned_02 --maf 0.05 --make-bed --out chr1_last30000_pruned_02_maf005

plink \
  --bfile chr1_last30000_pruned_02_maf005 \
  --recode A \
  --out chr1_last30000_pruned_02_maf005