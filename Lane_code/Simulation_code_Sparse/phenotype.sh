#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/Sparse_model/error/sim_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-300%10
#SBATCH --time=12:00:00

#infinitesimal sparse sparse_pair
#Contiguous  Random

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_Phenotype.py --n 1000 --m 1000 --s2gxg 0.9 --s2e 0.1 --mode Random --density 0.02 --sparseVersion SNP --rep $SLURM_ARRAY_TASK_ID

