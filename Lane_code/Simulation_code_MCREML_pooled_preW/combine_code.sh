#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1
RESULT_DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/result/${filename}

# Fail loudly rather than writing an empty .txt.  `cat missing/rep*.txt > out`
# still CREATES out (the shell truncates it before cat runs), so a tag mismatch
# between the pipeline and the Python steps used to surface as a silent empty
# file sitting next to a perfectly full result directory under a slightly
# different name.  Refuse to combine unless the directory really has reps, and
# point at what is actually on disk so the mismatch is obvious.
if [ ! -d "$RESULT_DIR" ]; then
    echo "ERROR: no result directory $RESULT_DIR" >&2
    echo "Existing result directories:" >&2
    ls -1d /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/result/*/ 2>/dev/null >&2
    exit 1
fi

nrep=$(ls -1 ${RESULT_DIR}/rep*.txt 2>/dev/null | wc -l)
if [ "$nrep" -eq 0 ]; then
    echo "ERROR: $RESULT_DIR contains no rep*.txt" >&2
    exit 1
fi

cat ${RESULT_DIR}/rep*.txt > /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW/result/${filename}.txt
echo "Combined $nrep reps -> result/${filename}.txt"
