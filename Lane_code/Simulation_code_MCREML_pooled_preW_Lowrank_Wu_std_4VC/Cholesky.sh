#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# $9 is R: which estimation kernel to precompute and cache under W/ --
# 0 (default) = the exact W the phenotype is drawn from, r > 0 = its dense
# rank-r truncation.
python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC/Simulate_Cholesky.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8 \
    --r ${9:-0}
