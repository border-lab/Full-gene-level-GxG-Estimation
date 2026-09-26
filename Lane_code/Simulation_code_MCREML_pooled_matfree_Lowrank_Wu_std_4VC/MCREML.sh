#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# MC-AI-REML is dominated by the MATRIX-FREE kernel applies inside every CG
# solve.  The low-rank epistasis apply is gemm either way -- two pairs against
# the stored At and D (--w_route storedA), or per gene against Z_g
# (--w_route bcast) -- and the additive and dominance terms are one more gemm
# pair each against their n-by-m design, so let BLAS use the cores allocated
# to this task.
#
# Arguments: 1-11 are positional (n m G s2a s2d s2gxg s2e mode iters nmc r);
# everything from 12 on is forwarded to Simulate_MCREML.py verbatim, which is
# how MCREML_pipeline.sh passes --verbose (possibly empty) and
# --w_route $W_ROUTE.  Named flags there, not positions, so an empty
# VERBOSE cannot shift the route into the wrong slot.
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC/Simulate_MCREML.py \
    --n $1 --m $2 --G $3 --s2a $4 --s2d $5 --s2gxg $6 --s2e $7 --mode $8 \
    --iters $9 --nmc ${10} \
    --r ${11:-20} \
    --rep $SLURM_ARRAY_TASK_ID "${@:12}"
