#!/bin/bash
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# AI-REML is dominated by the dense V inverse each iteration, so let BLAS
# use the cores allocated to this task (unlike the single-threaded MoM job).
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

python3 /home/ziyanzha/MOM_within_gene/AIREML_pooled/Simulate_AIREML_pooled.py \
    --n $1 --m $2 --s2a $3 --s2gxg $4 --s2e $5 --mode $6 \
    --iters $7 --G $8 --region_size $9 --seed ${10} \
    --rep $SLURM_ARRAY_TASK_ID
