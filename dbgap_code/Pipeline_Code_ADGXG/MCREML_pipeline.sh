#!/bin/bash

# ------------------------------------------------------------- How to Use

# MC-REML on ONE real phenotype, dbgap / UKB layout.  Run on the HEAD NODE:
#
#     bash MCREML_pipeline.sh
#
# One job: read int8 genotype (HM3 SNPs) + phenotype + covariates,
# fit the variance components in COMPONENTS (A, D, GxG) plus s2e, with every
# kernel product computed on the fly (multithreaded),
# write one-row TSV $DIR/<TRAIT>_<chr label>_<components>.tsv
set -u

# ------------------------------------------------------------- parameters
TRAIT=body_BMIz          # name of $PHEN_DIR/<TRAIT>.phen
CHROMS="22"              # all chromosomes: CHROMS="$(seq 1 22)"
COMPONENTS="A D GxG"     # subset of A D GxG, in this order (E is always included)
NMC=100                  # Hutchinson probes, fine phase
NMC_COARSE=15            # Hutchinson probes, coarse phase
SEED=1
GXG_RANK=9999            # rank r of the GxG kernel; 9999 = exact

VERBOSE=--verbose        # set to "" to drop the AI-REML trace

# -------------------------------------------------------------- resources
# int8 G = n x m bytes.  100K samples: chr22 ~1.7 GB, chr1-22 ~117 GB
# (whole genome: MEM=180G, CPUS=32, NODE=compute-0-30, TIME=7-00:00:00)
# CPUS is also the number of threads for the kernel products.
PARTITION=dbgap
CPUS=8
MEM=32G
TIME=24:00:00
NODE=""                  # pin a node, e.g. NODE=compute-0-30; "" = any node

# ------------------------------------------------------------------- data
# Compute-node paths.  GENO_DIR holds ukb_imp_chr{c}_v3.{pgen,pvar,psam};
# its samples are the samples of the fit (100K subset or full 337K).
# GENE_DIR holds gene_table_chr{c}.tsv (GxG gene sets, built on the 337K .pvar).
GENO_DIR=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_100K_HM3
PHEN_DIR=/dbgap/martinzh_ukb16549/ziyanzha/Phenotype/phen
GENE_DIR=/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3/gxg_gene_table

# ------------------------------------------------------------------ paths
# CODE is the HEAD-NODE path: sbatch reads run_MCREML.sh here at submit
# time.  DIR exists on compute nodes only.
CODE=/home/ziyanzha/Pipeline_Code_update
DIR=/dbgap/martinzh_ukb16549/ziyanzha/results/MCREML_real/$(basename $GENO_DIR)
LOG=$DIR/logs

# chr label, must match the name run_MCREML.py writes: chr22 or chr1-22
C=($CHROMS)
if [ ${#C[@]} -eq 1 ]; then CHR=chr${C[0]}; else CHR=chr${C[0]}-${C[-1]}; fi
COMP=${COMPONENTS// /}    # e.g. ADGxG; matches run_MCREML.py when COMPONENTS is in A D GxG order
RUN=${TRAIT}_${CHR}_${COMP}


# The log directory lives on the compute nodes and must exist before the job
# starts, or Slurm cannot open the job's output file.
if srun -p $PARTITION --chdir=/tmp mkdir -p $LOG; then
    :
else
    echo "ERROR: could not create $LOG on a compute node" >&2
    exit 1
fi


# ------------------------------------------------------------------ Submit the job
echo "run: $RUN"
echo "genotype: $GENO_DIR  ($CHR)"
echo "phenotype: $PHEN_DIR/$TRAIT.phen"
echo "components: $COMPONENTS"
case " $COMPONENTS " in *" GxG "*) echo "gene tables: $GENE_DIR  (rank $GXG_RANK)";; esac

JOB=$(sbatch --parsable \
    -p $PARTITION --nodes=1 --ntasks=1 --chdir=$DIR \
    --job-name=reml_$TRAIT \
    --output=$LOG/reml_${RUN}_%j.log \
    --mem=$MEM \
    --cpus-per-task=$CPUS \
    --time=$TIME \
    ${NODE:+-w $NODE} \
    $CODE/run_MCREML.sh \
        --trait $TRAIT --chroms $CHROMS --components $COMPONENTS \
        --geno_dir $GENO_DIR --phen_dir $PHEN_DIR --gene_dir $GENE_DIR --out_dir $DIR \
        --gxg_rank $GXG_RANK \
        --nmc $NMC --nmc_coarse $NMC_COARSE --seed $SEED $VERBOSE)
if [ -z "$JOB" ]; then echo "ERROR: submit failed" >&2; exit 1; fi
echo "MC-REML job: $JOB  (Nmc=$NMC, $CPUS CPUs, $MEM${NODE:+, $NODE})"

echo
echo "watch:   squeue -u $USER"
echo "log:     $LOG/reml_${RUN}_${JOB%%;*}.log   (compute nodes only)"
echo "result:  $DIR/$RUN.tsv"