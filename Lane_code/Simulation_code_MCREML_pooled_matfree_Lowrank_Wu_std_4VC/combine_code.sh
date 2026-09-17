#!/bin/bash
#
# Step 4 of MCREML_pipeline.sh: combine the per-replicate outputs, reduce the
# per-replicate diagnostics, and clean up.
#
#     combine_code.sh <FILENAME> <TAG>
#
# FILENAME is the result/timing key  <MODE>_s2a..._n<N>m<M>_G<G>  and TAG the
# Cholesky/Phenotype key  <MODE>_s2a..._n<N>_m<M>_G<G>  (the two differ only in
# the underscore before m).  Both are built once in MCREML_pipeline.sh and
# passed in, so this script never has to reconstruct them.
#
# WHY THIS IS A SCRIPT AND NOT AN sbatch --wrap CHAIN.  It used to be one long
# `a && b && c && rm && rm && ...` string.  That couples everything to
# everything: any single failure -- a reduction script not yet copied to the
# cluster, an empty input directory -- silently skipped every later step,
# including the cleanups, which is how a run ends with the per-rep result
# directory still on disk and no explanation.  Here each step is attempted
# INDEPENDENTLY and reports for itself, and the exit status at the bottom is
# non-zero if any of them failed, so the job shows FAILED and error/combine_*
# .err says which part.
#
# WHAT IS DELETED, AND WHEN.  A per-replicate directory is removed ONLY IF the
# thing that consumes it actually produced output.  Losing 200 rep files to a
# tidy-up that ran after its own reduction failed is not recoverable without
# rerunning the whole array, so:
#
#     result/<FILENAME>/         removed iff result/<FILENAME>.txt is non-empty
#     time/rep_times/<FILENAME>/ removed iff time/result/timing_<FILENAME>.txt
#                                exists
#     time/op_counts/<FILENAME>/ removed iff time/result/opcount_<FILENAME>.txt
#                                exists
#     Cholesky factors, Phenotype/y_<TAG>/
#                                removed iff the estimates combined -- they are
#                                the big artifacts and are regenerable, but not
#                                worth discarding while the run is still broken
#
set -u

if [ $# -lt 2 ]; then
    echo "Usage: $0 <FILENAME> <TAG>" >&2
    exit 1
fi

filename=$1
tag=$2

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_std_4VC

status=0

# ---------------------------------------------------------------- estimates
# One row per replicate, FOUR columns,
#
#   (V_a, V_d, V_gamma, V_e)
#
# the four variance components REML fitted and nothing else.  All four are on
# the realized-variance scale (the kernel carries the 1/c-hat, so no post-fit
# correction is applied), so the column means are directly comparable to the
# nominal S2A / S2D / S2GXG / S2E.
#
# THE ROW USED TO BE NINE COLUMNS: these four, then V_l (identically V_gamma
# here, kept only to match the sibling layouts), then the realized variances
# Var-hat(.) of that replicate's four effect draws.  The Phenotype step now
# rescales each drawn component to hit its target exactly, which made those
# four realized columns constants equal to the targets, so they and V_l are
# gone.  A result file from THIS pipeline is therefore narrower than a
# sibling's -- read it with this directory's calc_stats.py.
REP_DIR=$DIR/result/$filename
OUT=$DIR/result/${filename}.txt
combined=0

n_rep=$(ls -1 $REP_DIR/rep*.txt 2>/dev/null | wc -l)
if [ "$n_rep" -gt 0 ]; then
    cat $REP_DIR/rep*.txt > $OUT
    if [ -s "$OUT" ]; then
        combined=1
        echo "combined $n_rep replicate rows -> $OUT"
    else
        echo "ERROR: $OUT came out empty from $n_rep rep files" >&2
        status=1
    fi
else
    echo "ERROR: no rep*.txt under $REP_DIR -- nothing to combine" >&2
    status=1
fi

# ------------------------------------------------------- timing reduction
# Mean wall-clock estimation time over replicates -> time/result/timing_*.txt.
TIME_REP=$DIR/time/rep_times/$filename
TIME_OUT=$DIR/time/result/timing_${filename}.txt
if python3 $DIR/time/average_time.py $TIME_REP $TIME_OUT; then
    :
else
    echo "ERROR: average_time.py failed on $TIME_REP" >&2
    status=1
fi

# --------------------------------------------- operator-apply reduction
# time/op_counts/<FILENAME>/rep*.txt  ->  ONE file, time/result/opcount_*.txt:
# mean / std / min / max over replicates of how many applies of V, of K_a/K_d
# and of the rank-r epistasis operator one estimation costs, and to how many
# columns each.  This is the machine-independent companion to the timing above
# -- fixed by (genotype, r, Nmc, tolerances), unmoved by BLAS threads or node
# load -- and the per-rep directory is removed once it exists.
#
# DONE HERE, IN awk, ON PURPOSE.  It used to be a separate python helper
# (time/average_op_counts.py), which meant the reduction silently did not
# happen on a cluster where that one extra file had not been copied up yet --
# leaving the op_counts directory sitting there with no summary and no obvious
# reason.  This script has to be deployed for ANY of the combine to work, so
# folding the arithmetic into it removes that failure mode entirely: there is
# now nothing to forget.  awk is in coreutils, so it also needs no python and
# no conda, which this job does not load.
#
# Keys are emitted in FIRST-SEEN order, which is the fixed order
# Simulate_MCREML.py writes them in, so the summary reads the same every run
# and picks up new counters automatically if that list ever grows.  std is
# ddof=0, matching the realized-variance summaries elsewhere.
OP_REP=$DIR/time/op_counts/$filename
OP_OUT=$DIR/time/result/opcount_${filename}.txt
n_op=$(ls -1 $OP_REP/rep*.txt 2>/dev/null | wc -l)
if [ "$n_op" -gt 0 ]; then
    mkdir -p $DIR/time/result
    awk '
        FNR == 1 { nfile++ }
        NF == 2 {
            k = $1; v = $2 + 0
            if (!(k in cnt)) keys[++nk] = k
            cnt[k]++; sum[k] += v; sq[k] += v * v
            if (cnt[k] == 1 || v < min[k]) min[k] = v
            if (cnt[k] == 1 || v > max[k]) max[k] = v
        }
        END {
            for (i = 1; i <= nk; i++) {
                k = keys[i]; n = cnt[k]
                mean = sum[k] / n
                var = sq[k] / n - mean * mean
                if (var < 0) var = 0            # rounding, not a real negative
                printf "%s_mean %.4f\n", k, mean
                printf "%s_std %.4f\n",  k, sqrt(var)
                printf "%s_min %.0f\n",  k, min[k]
                printf "%s_max %.0f\n",  k, max[k]
            }
            printf "n_reps %d\n", nfile
        }
    ' $OP_REP/rep*.txt > $OP_OUT
    if [ -s "$OP_OUT" ]; then
        echo "averaged $n_op replicates' operator counts -> $OP_OUT"
    else
        echo "ERROR: $OP_OUT came out empty from $n_op rep files" >&2
        status=1
    fi
else
    echo "ERROR: no rep*.txt under $OP_REP -- no operator counts to average." \
         "Is Simulate_MCREML.py on the cluster the version that writes them?" >&2
    status=1
fi

# ------------------------------------------------------------------ cleanup
# Each per-rep tree goes only when its own reduction has landed.
if [ "$combined" -eq 1 ]; then
    rm -rf $REP_DIR
    rm -f $DIR/Cholesky/La_${tag}.npy
    rm -f $DIR/Cholesky/Ld_${tag}.npy
    rm -f $DIR/Cholesky/Lgxg_${tag}.npy
    rm -rf $DIR/Phenotype/y_${tag}
    echo "cleaned: $REP_DIR, the three Cholesky factors, Phenotype/y_${tag}"
else
    echo "KEPT $REP_DIR and the Cholesky/Phenotype artifacts: the estimates" \
         "did not combine" >&2
fi

if [ -f "$TIME_OUT" ]; then
    rm -rf $TIME_REP
    echo "cleaned: $TIME_REP"
else
    echo "KEPT $TIME_REP: no $TIME_OUT was written" >&2
fi

if [ -f "$OP_OUT" ]; then
    rm -rf $OP_REP
    echo "cleaned: $OP_REP"
else
    echo "KEPT $OP_REP: no $OP_OUT was written" >&2
fi

exit $status
