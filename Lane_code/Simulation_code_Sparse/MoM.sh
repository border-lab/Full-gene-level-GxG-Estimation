#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_MoM.py \
    --n $1 --m $2 --s2a $3 --s2gxg $4 --s2e $5 \
    --mode $6 --density $7 --sparseVersion $8 \
    --nmc $9 --rep $SLURM_ARRAY_TASK_ID