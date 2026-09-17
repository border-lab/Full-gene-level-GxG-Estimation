#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the kernel applies inside every CG solve.  Here the
# epistasis apply is ONE dense n-by-n gemm against the PRECOMPUTED W (loaded
# once per replicate from W/, written by the Cholesky job), and the additive and
# dominance terms are one more gemm pair each against their n-by-m design, so
# let BLAS use the cores allocated to this task.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# ${11} is R: which cached kernel to load -- 0 (default) = the exact W,
# r > 0 = the dense rank-r W-hat.  It must match the Cholesky job's R.
python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_Lowrank_Wu_std_4VC/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8 \
    --iters $9 --nmc ${10} \
    --r ${11:-0} \
    --rep $SLURM_ARRAY_TASK_ID ${12:-}
