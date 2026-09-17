#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# Exact AI-REML is dominated by one dense O(n^3) Cholesky + potri of V per
# score/AI evaluation, so let BLAS use the cores allocated to this task
# (CPUS_REML in AIREML_pipeline.sh).
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# ${10} is the optional --verbose flag.
python3 /home/ziyanzha/MOM_within_gene/AIREML_pooled_preW_exact_Wu_std_4VC/Simulate_AIREML.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8 \
    --iters $9 \
    --rep $SLURM_ARRAY_TASK_ID ${10:-}
