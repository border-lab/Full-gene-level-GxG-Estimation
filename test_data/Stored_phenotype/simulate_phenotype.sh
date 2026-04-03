#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/error/sim_%A_%a.err
#SBATCH --output=/home/ziyanzha/MOM_within_gene/outputInf/sim_%A_%a.out
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-300%33
#SBATCH --time=12:00:00

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Stored_data/phenotype.py --n 16000 --m 4000 --s2gxg 0.9 --s2e 0.1 --rep $SLURM_ARRAY_TASK_ID --mode Contiguous

#python3 /home/ziyanzha/MOM_within_gene/Stored_data/phenotype_sparse.py --n 16000 --m 4000 --s2gxg 0.9 --s2e 0.02 --density 0.02 --rep $SLURM_ARRAY_TASK_ID --mode Random
