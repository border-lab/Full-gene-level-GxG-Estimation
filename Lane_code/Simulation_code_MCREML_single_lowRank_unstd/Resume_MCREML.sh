#!/bin/bash
# Resume an interrupted Step-3 (MC-AI-REML) array job.
#
# Re-runs ONLY the replicates whose result/<FILENAME>/rep{r}.txt is missing or
# empty.  Safe to run repeatedly: seed=rep and the low-rank operator is
# deterministic, so a resumed replicate is bit-identical to what the killed job
# would have written -- finishing a run in two batches biases nothing.
#
# It does NOT touch the Cholesky factor or the phenotypes, and it does NOT run
# the combine job (combine DELETES Lgxg, the phenotype dir and the per-rep
# result dir -- you want that only once all reps are in).  See the end of this
# script for the combine command to run afterwards.

N=16000
M=1000
G=10
S2GXG=0.2
S2E=0.8
MODE=ContiguousSNP
ITERS=30
NMC=50
R=100
NREPS=200
MEM=2G          # use whatever the ORIGINAL run used, not the default
ARRAY=40
PARTITION=statgen-gpu

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd

# The two tags differ by ONE underscore -- result/ uses n{n}m{m}, while
# Cholesky/ and Phenotype/ use n{n}_m{m}.  Getting this wrong is the easy
# mistake here, so both are spelled out.
TAG=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}_m${M}_G${G}
FILENAME=${MODE}_s2gxg${S2GXG}_s2e${S2E}_n${N}m${M}_G${G}

RESULT_DIR=$DIR/result/$FILENAME
PHENO_DIR=$DIR/Phenotype/y_${TAG}
LGXG_FILE=$DIR/Cholesky/Lgxg_${TAG}.npy

mkdir -p $DIR/error

# --- which replicates are still missing? -----------------------------------
# -s (non-empty) rather than -f: a replicate killed mid-write leaves a 0-byte
# file, which must count as missing.
missing=()
for r in $(seq 1 $NREPS); do
    [ -s "$RESULT_DIR/rep${r}.txt" ] || missing+=($r)
done

if [ ${#missing[@]} -eq 0 ]; then
    echo "All $NREPS replicates present in $RESULT_DIR -- nothing to resume."
    echo "Run the combine step (see bottom of this script)."
    exit 0
fi

echo "Present: $(( NREPS - ${#missing[@]} ))/$NREPS      Missing: ${#missing[@]}"

# --- are the inputs those replicates need still on disk? -------------------
# Step 3 reads the PHENOTYPE only.  Lgxg matters solely if a phenotype has to
# be redrawn -- and at n=$N rebuilding it is the expensive job you want to
# avoid, so report its state rather than silently depending on it.
no_pheno=()
for r in "${missing[@]}"; do
    if [ ! -s "$PHENO_DIR/rep${r}.csv" ] || [ ! -s "$PHENO_DIR/vell_rep${r}.txt" ]; then
        no_pheno+=($r)
    fi
done

if [ ${#no_pheno[@]} -ne 0 ]; then
    echo "ERROR: ${#no_pheno[@]} of the missing replicates have no phenotype in"
    echo "       $PHENO_DIR"
    echo "       reps: ${no_pheno[*]}"
    if [ -s "$LGXG_FILE" ]; then
        echo "  Lgxg IS present, so redraw just those phenotypes first:"
        echo "    sbatch -p $PARTITION --mem=$MEM --time=12:00:00 \\"
        echo "      --array=$(IFS=,; echo "${no_pheno[*]}")%$ARRAY \\"
        echo "      --error=$DIR/error/phenotype_%A.err --output=/dev/null \\"
        echo "      $DIR/Phenotype.sh $N $M $G $S2GXG $S2E $MODE"
        echo "  then re-run this script."
    else
        echo "  Lgxg is ALSO gone ($LGXG_FILE)."
        echo "  The combine job must have run and cleaned up.  Rebuilding it at"
        echo "  n=$N means redoing the O(n^2 m) W build and the O(n^3) Cholesky"
        echo "  -- decide deliberately before doing that."
    fi
    exit 1
fi

# --- submit exactly the missing indices ------------------------------------
# sbatch takes an explicit index list, so this is the same array job as Step 3
# with a sparse index set; the throttle (%$ARRAY) is kept.
ARRAY_SPEC=$(IFS=,; echo "${missing[*]}")%$ARRAY
echo "Submitting --array=$ARRAY_SPEC"

JOB=$(sbatch --parsable \
    --job-name=mcreml_resume \
    -p $PARTITION \
    --error=$DIR/error/mcreml_%A.err \
    --open-mode=append \
    --output=/dev/null \
    --nodes=1 \
    --mem=$MEM \
    --ntasks=1 \
    --cpus-per-task=1 \
    --time=12:00:00 \
    --array=$ARRAY_SPEC \
    $DIR/MCREML.sh $N $M $G $S2GXG $S2E $MODE $ITERS $NMC $R)

echo "Resume job: $JOB   (${#missing[@]} replicates, r=$R, Nmc=$NMC)"
echo
echo "When it finishes, re-run this script to confirm 0 missing, THEN combine:"
echo "  bash $DIR/combine_code.sh $FILENAME && \\"
echo "  python3 $DIR/time/average_time.py $DIR/time/rep_times/$FILENAME \\"
echo "          $DIR/time/result/timing_$FILENAME.txt && \\"
echo "  rm -rf $RESULT_DIR $DIR/time/rep_times/$FILENAME $PHENO_DIR && \\"
echo "  rm -f $LGXG_FILE"
