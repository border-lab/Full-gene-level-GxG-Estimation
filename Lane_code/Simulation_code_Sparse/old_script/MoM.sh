#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/Sparse_model/error/sim_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-300%5
#SBATCH --time=12:00:00


#infinitesimal sparse pair
#Contiguous  Random

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=1 
export OMP_NUM_THREADS=1 

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_MoM.py --n 1000 --m 1000 --s2a 0.2 --s2gxg 0.7 --s2e 0.1 --mode Random --density 1 --nmc 100 --sparseVersion SNP --rep $SLURM_ARRAY_TASK_ID

