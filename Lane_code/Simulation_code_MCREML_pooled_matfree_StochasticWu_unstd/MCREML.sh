#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE epistasis apply inside every CG
# solve.  For the UNSTANDARDIZED kernel that apply is pure gemm against Z_g --
# Z_g'(u .* K_g V) then the back-contraction, with no m-by-m weight matrix
# anywhere -- so let BLAS use the cores allocated to this task.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_StochasticWu_unstd/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2gxg $4 --s2e $5 --mode $6 \
    --iters $7 --nmc $8 \
    --trace_method ${9:-slq} --slq_k ${10:-25} \
    --nw ${11:-200} \
    --rep $SLURM_ARRAY_TASK_ID
