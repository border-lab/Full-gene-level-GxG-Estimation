#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the epistasis mat-vec (W B) inside every CG solve,
# so let BLAS use the cores allocated to this task.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_three_var/Simulate_MCREML.py \
    --n $1 --m $2 --s2a $3 --s2d $4 --s2gxg $5 --s2e $6 --mode $7 \
    --iters $8 --nmc $9 \
    --rep $SLURM_ARRAY_TASK_ID
