#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the thin genotype mat-vecs Z (Z' b) inside every CG
# solve, so let BLAS use the cores allocated to this task.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_additive/Simulate_MCREML.py \
    --n $1 --m $2 --s2a $3 --s2e $4 --mode $5 \
    --iters $6 --nmc $7 \
    --rep $SLURM_ARRAY_TASK_ID
