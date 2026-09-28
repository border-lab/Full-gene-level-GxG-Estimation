#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE kernel applies inside every CG
# solve.  For the low-rank epistasis operator that apply is pure gemm against
# Z_g -- Z_g'(q_s .* u) then the back-contraction, r terms per gene with no
# m-by-m weight matrix anywhere -- and it is the ONLY kernel in this model, so
# let BLAS use the cores allocated to this task.  For a single-thread profile
# of the per-column costs, launch with --cpus-per-task=1.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2gxg $4 --s2e $5 --mode $6 \
    --iters $7 --nmc $8 \
    --r ${9:-20} \
    --rep $SLURM_ARRAY_TASK_ID ${10:-}
