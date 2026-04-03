#!/bin/bash
#SBATCH --job-name=sim
#SBATCH -p mzhang
#SBATCH --error=/home/ziyanzha/MOM_within_gene/error/sim_%A_%a.err
#SBATCH --output=/home/ziyanzha/MOM_within_gene/outputInf/sim_%A_%a.out
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-300%5
#SBATCH --time=12:00:00

source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=1 
export OMP_NUM_THREADS=1 

python3 /home/ziyanzha/MOM_within_gene/run_MoM_std.py --n 8000 --m 4000 --s2gxg 0.9 --s2e 0.1 --rep $SLURM_ARRAY_TASK_ID --mode Contiguous
