#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/Sparse_model/error/sim_%j.err
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

#infinitesimal sparse sparse_pair
#Contiguous  Random

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_Cholesky.py --n 32000 --m 1000 --s2gxg 0.9 --s2e 0.1 --mode Random --density 0.05 --sparseVersion SNP

