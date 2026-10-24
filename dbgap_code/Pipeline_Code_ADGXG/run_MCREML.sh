#!/bin/bash
# Job script: fit MC-REML on ONE real phenotype.  Submitted by MCREML_pipeline.sh, e.g.
#   sbatch run_MCREML.sh --trait body_BMIz --chroms 22 --components A --verbose
#   sbatch run_MCREML.sh --trait body_BMIz --chroms 22 --components A D GxG --gene_dir <dir> --verbose
# All arguments are passed unchanged to run_MCREML.py.
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv

# The kernel products (A, D, GxG) are parallelized by their own threads (--threads);
# BLAS runs 1 thread per worker so the workers do not oversubscribe the cores.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1      # write the log as it goes

CODE=/cluster-home/ziyanzha/Pipeline_Code_update
cd $CODE                       # so Python finds Function_input / Function_estimation
python3 $CODE/run_MCREML.py --threads ${SLURM_CPUS_PER_TASK:-1} "$@"