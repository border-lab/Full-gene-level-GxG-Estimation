#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE kernel applies inside every CG
# solve.  Here every apply is one gemm pair against the n-by-m additive design
# and nothing else, so let BLAS use the cores allocated to this task.
#
# Arguments: 1-7 are positional (n m s2a s2e mode iters nmc); everything from
# 8 on is forwarded to Simulate_MCREML.py verbatim, which is how
# MCREML_pipeline.sh passes --verbose (possibly empty).
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_additive/Simulate_MCREML.py \
    --n $1 --m $2 --s2a $3 --s2e $4 --mode $5 \
    --iters $6 --nmc $7 \
    --rep $SLURM_ARRAY_TASK_ID "${@:8}"
