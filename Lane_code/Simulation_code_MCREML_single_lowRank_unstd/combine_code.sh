#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

DIR=/home/ziyanzha/MOM_within_gene/MCREML_pooled_matfree_Lowrank_Wu_unstd

# Estimates: one "(V_gamma,V_e,V_l,realized_variance)" row per replicate --
# the two variance components, the c-corrected estimate V_l = c_hat*s2gxg_hat,
# and the realized variance Var-hat(H gamma) of that replicate's gamma draw
# (carried through from the Phenotype step by Simulate_MCREML.py).  Nothing
# else to combine: the realized variances live in this file now, one per row,
# instead of being reduced to a mean/std under result/realized_variance/.
cat $DIR/result/${filename}/rep*.txt > $DIR/result/${filename}.txt
