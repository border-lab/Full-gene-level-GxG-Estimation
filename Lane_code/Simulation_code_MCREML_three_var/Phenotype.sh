#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}


python3 /home/ziyanzha/MOM_within_gene/MCREML_three_var/Simulate_Phenotype.py \
    --n $1 --m $2 --s2a $3 --s2d $4 --s2gxg $5 --s2e $6 --mode $7 \
    --rep $SLURM_ARRAY_TASK_ID
