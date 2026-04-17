#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/Sparse_model/error/sim_%j.err
#SBATCH --nodes=1
#SBATCH --mem=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

#infinitesimal sparse pair
#Contiguous  Random

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_Cholesky.py --n 1000 --m 1000 --s2a 0.2 --s2gxg 0.7 --s2e 0.1 --mode Random --density 1 --sparseVersion pair

