#!/bin/bash
#
# Step 4 of AIREML_pipeline.sh: combine the per-replicate outputs, reduce the
# per-replicate diagnostics, and clean up.
#
#     combine_code.sh <FILENAME> <TAG>
#
# FILENAME is the result/timing key  <MODE>_s2a..._n<N>m<M>_G<G>  and TAG the
# Cholesky/Phenotype key  <MODE>_s2a..._n<N>_m<M>_G<G>  (the two differ only in
# the underscore before m).  Both are built once in AIREML_pipeline.sh and
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
#     Cholesky factors, Phenotype/y_<TAG>/
#                                removed iff the estimates combined
#
set -u

if [ $# -lt 2 ]; then
    echo "Usage: $0 <FILENAME> <TAG>" >&2
    exit 1
fi

filename=$1
tag=$2

DIR=/home/ziyanzha/MOM_within_gene/AIREML_pooled_preW_exact_Wu_std_4VC

status=0

# ---------------------------------------------------------------- estimates
# One row per replicate, FOUR columns,
#
#   (V_a, V_d, V_gamma, V_e)
#
# the four variance components exact AI-REML fitted and nothing else.  All four
# are on the realized-variance scale (the kernel carries the 1/c-hat), so the
# column means are directly comparable to the nominal S2A / S2D / S2GXG / S2E.
# Read the file with this directory's calc_stats.py.
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

# ------------------------------------------------------ cost-count reduction
# time/op_counts/<FILENAME>/rep*.txt  ->  ONE file, time/result/opcount_*.txt:
# mean / std / min / max over replicates of reml_iters, V_factorizations (each
# a dense O(n^3) Cholesky + potri of V) and rejected_steps.  This is the
# machine-independent companion to the timing above -- fixed by (genotype, y,
# tolerances), unmoved by BLAS threads or node load.
#
# Done in awk so the reduction needs no extra python helper on the cluster.
# Keys are emitted in FIRST-SEEN order, which is the fixed order
# Simulate_AIREML.py writes them in.  std is ddof=0.
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
        echo "averaged $n_op replicates' cost counts -> $OP_OUT"
    else
        echo "ERROR: $OP_OUT came out empty from $n_op rep files" >&2
        status=1
    fi
else
    echo "ERROR: no rep*.txt under $OP_REP -- no cost counts to average." >&2
    status=1
fi

# ------------------------------------------------------------------ cleanup
# Each per-rep tree goes only when its own reduction has landed.
#
# The PRECOMPUTED kernel W/W_<mode>_n<n>_m<m>_G<G>.npy is NOT removed: it does
# not depend on the variance targets, so other runs on the same genotype use
# the same matrix (their Cholesky jobs rewrite it identically).  It is 8 n^2
# bytes -- delete it by hand when a genotype is retired.
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
