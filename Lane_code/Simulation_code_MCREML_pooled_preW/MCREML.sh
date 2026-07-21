#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the epistasis mat-vec (W B) inside every CG solve,
# so let BLAS use the cores allocated to this task.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2gxg $4 --s2e $5 --mode $6 --ratio "$7" \
    --iters $8 --nmc $9 \
    --rep $SLURM_ARRAY_TASK_ID
