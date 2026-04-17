#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

python3 /home/ziyanzha/MOM_within_gene/Only_W_model/Simulate_Cholesky.py \
    --n $1 --m $2 --s2gxg $3 --s2e $4 --mode $5