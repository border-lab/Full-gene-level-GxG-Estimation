#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Builds the three Cholesky factors and caches the exact W under W/ for REML.
python3 /home/ziyanzha/MOM_within_gene/AIREML_pooled_preW_exact_Wu_std_4VC/Simulate_Cholesky.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8
