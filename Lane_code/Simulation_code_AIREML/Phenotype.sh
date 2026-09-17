#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/AIREML/Simulate_Phenotype.py \
    --n $1 --m $2 --s2a $3 --s2d $4 --s2gxg $5 --s2e $6 --mode $7 \
    --rep $SLURM_ARRAY_TASK_ID
