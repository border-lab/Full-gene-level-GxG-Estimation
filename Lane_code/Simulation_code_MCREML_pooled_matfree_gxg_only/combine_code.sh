#!/bin/bash
#
# Step 4 of MCREML_pipeline.sh: combine the per-replicate outputs, reduce the
# per-replicate diagnostics, and clean up.
#
#     combine_code.sh <FILENAME> <TAG>
#
# FILENAME is the result/timing key  <MODE>_s2gxg..._n<N>m<M>_G<G>  and TAG the
# Cholesky/Phenotype key  <MODE>_s2gxg..._n<N>_m<M>_G<G>  (the two differ only
# in the underscore before m).  Both are built once in MCREML_pipeline.sh and
# passed in, so this script never has to reconstruct them.
#
# Each step is attempted INDEPENDENTLY and reports for itself, and the exit
# status at the bottom is non-zero if any of them failed, so the job shows
# FAILED and error/combine_*.err says which part.
#
# WHAT IS DELETED, AND WHEN.  A per-replicate directory is removed ONLY IF the
# thing that consumes it actually produced output:
#
#     result/<FILENAME>/         removed iff result/<FILENAME>.txt is non-empty
#     time/rep_times/<FILENAME>/ removed iff time/result/timing_<FILENAME>.txt
#                                exists
#     time/op_counts/<FILENAME>/ removed iff time/result/opcount_<FILENAME>.txt
#                                exists
#     Cholesky factor, Phenotype/y_<TAG>/
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

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_gxg_only

status=0

# ---------------------------------------------------------------- estimates
# One row per replicate, TWO columns,
#
#   (V_gamma, V_e)
#
# the two variance components REML fitted and nothing else.  Both are on the
# realized-variance scale (the kernel carries the 1/c-hat, so no post-fit
# correction is applied), so the column means are directly comparable to the
# nominal S2GXG / S2E.  Read the combined file with this directory's
# calc_stats.py.
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
# Mean wall-clock estimation time over replicates -> time/result/timing_*.txt
# (one number, seconds, 4 decimals).
#
# IN awk, NOT python.  This job is launched with sbatch --wrap and never
# activates conda, so python3 is whatever the bare compute node has -- often
# nothing.  Every reduction in this script needs only coreutils and awk.
TIME_REP=$DIR/time/rep_times/$filename
TIME_OUT=$DIR/time/result/timing_${filename}.txt
n_time=$(ls -1 $TIME_REP/rep*.txt 2>/dev/null | wc -l)
if [ "$n_time" -gt 0 ]; then
    mkdir -p $DIR/time/result
    awk '
        NF >= 1 && $1 ~ /^[0-9.eE+-]+$/ { sum += $1 + 0; n++ }
        END {
            if (n == 0) exit 1
            printf "%.4f\n", sum / n
        }
    ' $TIME_REP/rep*.txt > $TIME_OUT
    if [ -s "$TIME_OUT" ]; then
        echo "averaged $n_time replicates' estimation times -> $TIME_OUT" \
             "($(cat $TIME_OUT) s)"
    else
        echo "ERROR: $TIME_OUT came out empty from $n_time rep files" \
             "(are the rep*.txt one number each?)" >&2
        rm -f $TIME_OUT
        status=1
    fi
else
    echo "ERROR: no rep*.txt under $TIME_REP -- no estimation times to average." \
         "Did Step 3 finish, or was time/rep_times/ already cleaned?" >&2
    status=1
fi

# --------------------------------------------- operator-apply reduction
# time/op_counts/<FILENAME>/rep*.txt  ->  ONE file, time/result/opcount_*.txt:
# mean / std / min / max over replicates of how many applies of V and of the
# rank-r epistasis operator one estimation costs, to how many columns each,
# and of the wall-clock timers written next to them.  The counts are the
# machine-independent companion to the timing above -- fixed by (genotype, r,
# Nmc, tolerances), unmoved by BLAS threads or node load.
#
# DONE HERE, IN awk, ON PURPOSE: this script has to be deployed for ANY of the
# combine to work, so folding the arithmetic into it means there is nothing
# else to forget to copy up, and it needs no python and no conda.
#
# Keys are emitted in FIRST-SEEN order, which is the fixed order
# Simulate_MCREML.py writes them in, so the summary reads the same every run
# and picks up new counters automatically if that list ever grows.  std is
# ddof=0, matching the summaries elsewhere.  Keys ending "_sec" are the
# wall-clock timers (seconds on the node that ran the rep); same reduction,
# decimals kept.
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
                # counts: 4 decimals on mean/std, integers on min/max.  The
                # "_sec" timers (per-column ones are ~1e-4) get 6 throughout.
                if (k ~ /_sec/) { fm = "%.6f"; fx = "%.6f" }
                else            { fm = "%.4f"; fx = "%.0f" }
                printf "%s_mean " fm "\n", k, mean
                printf "%s_std "  fm "\n", k, sqrt(var)
                printf "%s_min "  fx "\n", k, min[k]
                printf "%s_max "  fx "\n", k, max[k]
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
    rm -f $DIR/Cholesky/Lgxg_${tag}.npy
    rm -rf $DIR/Phenotype/y_${tag}
    echo "cleaned: $REP_DIR, the Cholesky factor, Phenotype/y_${tag}"
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
