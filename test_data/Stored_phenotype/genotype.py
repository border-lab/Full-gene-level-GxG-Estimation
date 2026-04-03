import pandas as pd
import numpy as np 



Contiguous= pd.read_csv("/home/ziyanzha/MOM_within_gene/genotype_matrix_plink/chr1_n50k_m5k.raw", sep=r'\s+')
real_data = Contiguous.iloc[:, 6:].to_numpy()

m = 1000
n = 1000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n1000_m1000.csv", index=False, header=False)

n = 2000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n2000_m1000.csv", index=False, header=False)

n = 4000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n4000_m1000.csv", index=False, header=False)

n = 8000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n8000_m1000.csv", index=False, header=False)

n = 16000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n16000_m1000.csv", index=False, header=False)

n = 32000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n32000_m1000.csv", index=False, header=False)


m = 2000
n = 1000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n1000_m2000.csv", index=False, header=False)

n = 2000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n2000_m2000.csv", index=False, header=False)

n = 4000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n4000_m2000.csv", index=False, header=False)

n = 8000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n8000_m2000.csv", index=False, header=False)

n = 16000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n16000_m2000.csv", index=False, header=False)

n = 32000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n32000_m2000.csv", index=False, header=False)


m = 4000
n = 1000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n1000_m4000.csv", index=False, header=False)

n = 2000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n2000_m4000.csv", index=False, header=False)

n = 4000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n4000_m4000.csv", index=False, header=False)

n = 8000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n8000_m4000.csv", index=False, header=False)

n = 16000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n16000_m4000.csv", index=False, header=False)

n = 32000
Z_1000 = real_data[:n, :m]
pd.DataFrame(Z_1000).to_csv("/home/ziyanzha/MOM_within_gene/Stored_data/genotype/ContiguousSNP_n32000_m4000.csv", index=False, header=False)

