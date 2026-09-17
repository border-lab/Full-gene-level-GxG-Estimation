#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}


python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd/Simulate_Phenotype.py \
    --n $1 --m $2 --G $3 --s2gxg $4 --s2e $5 --mode $6 \
    --rep $SLURM_ARRAY_TASK_ID
