#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

cat /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_perG/result/${filename}/rep*.txt > /home/ziyanzha/MOM_within_gene/MCREML_pooled_preW_perG/result/${filename}.txt
