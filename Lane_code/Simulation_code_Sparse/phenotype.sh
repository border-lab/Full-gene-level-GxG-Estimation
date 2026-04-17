#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Sparse_model/Simulate_Phenotype.py \
    --n $1 --m $2 --s2a $3 --s2gxg $4 --s2e $5 \
    --mode $6 --density $7 --sparseVersion $8 \
    --rep $SLURM_ARRAY_TASK_ID