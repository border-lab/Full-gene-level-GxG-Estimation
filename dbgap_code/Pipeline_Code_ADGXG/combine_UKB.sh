#!/bin/bash
# Fit MC-AI-REML on ONE real phenotype.  Submit with sbatch, e.g.
#   sbatch run_MCREML.sh N M PGEN GENE_TABLE PHEN_DIR TRAIT ITERS NMC R [extra flags]
#   PGEN       plink2 prefix (without .pgen / .psam)
#   GENE_TABLE gene table (idx_start / idx_stop index the rows of PGEN.pvar)
#   PHEN_DIR   folder with <trait>.phen files (FID, IID, <trait>)
#   TRAIT      trait name, e.g. body_BMIz
# Anything after the 9th argument (e.g. --verbose) is passed to run_MCREML.py.
source /home/ziyanzha/miniforge3/etc/profile.d/conda.sh
conda activate tfenv
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

if [ $# -lt 9 ]; then
    echo "Usage: $0 N M PGEN GENE_TABLE PHEN_DIR TRAIT ITERS NMC R [extra flags]" >&2
    exit 1
fi

CODE=/cluster-home/ziyanzha/Simulation_Code/MCREML_real
TRAIT=$6

python3 $CODE/run_MCREML.py \
    --n $1 --m $2 --pgen $3 --gene_table $4 \
    --phen $5/$TRAIT.phen --trait $TRAIT \
    --iters $7 --nmc $8 --r $9 "${@:10}"