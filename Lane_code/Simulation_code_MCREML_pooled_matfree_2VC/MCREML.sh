#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE kernel applies inside every CG
# solve.  Here every apply is a gemm pair against an n-by-m design (K_a and
# K_d, one each per V apply) and nothing else, so let BLAS use the cores
# allocated to this task.  For a single-thread profile of the per-column
# costs, launch with --cpus-per-task=1.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_2VC/Simulate_MCREML.py \
    --n $1 --m $2 --s2a $3 --s2d $4 --s2e $5 --mode $6 \
    --iters $7 --nmc $8 \
    --rep $SLURM_ARRAY_TASK_ID ${9:-}
