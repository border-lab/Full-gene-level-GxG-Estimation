#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE kernel applies inside every CG
# solve.  For the low-rank epistasis operator that apply is pure gemm against
# Z_g -- Z_g'(q_s .* u) then the back-contraction, r terms per gene with no
# m-by-m weight matrix anywhere -- and the additive and dominance terms are one
# more gemm pair each against their n-by-m design, so let BLAS use the cores
# allocated to this task.  The c-normalization adds one Gram matrix per gene in
# the setup phase (exact c, computed once outside the REML loop) and a scalar
# division inside the apply -- neither changes what this job needs from BLAS.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8 \
    --iters $9 --nmc ${10} \
    --r ${11:-20} \
    --rep $SLURM_ARRAY_TASK_ID ${12:-}
