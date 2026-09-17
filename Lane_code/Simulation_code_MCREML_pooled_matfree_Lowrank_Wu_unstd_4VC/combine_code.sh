#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd_4VC

# Estimates: one row per replicate,
#
#   (V_a, V_d, V_gamma, V_e, V_l, realized_gxg, realized_a, realized_d, realized_e)
#
# the four variance components REML fitted, the c-corrected estimate
# V_l = c_hat*s2gxg_hat, and the realized variances Var-hat(.) of that
# replicate's four effect draws (carried through from the Phenotype step by
# Simulate_MCREML.py).  Nothing else to combine: the realized variances live in
# this file now, one row per replicate, instead of being reduced to a mean/std
# under result/realized_variance/.
cat $DIR/result/${filename}/rep*.txt > $DIR/result/${filename}.txt
