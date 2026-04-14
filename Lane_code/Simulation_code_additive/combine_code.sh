#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

cat /home/ziyanzha/MOM_within_gene/Whole_model/result/${filename}/rep*.txt > /home/ziyanzha/MOM_within_gene/Whole_model/result/${filename}.txt

python3 /home/ziyanzha/MOM_within_gene/Whole_model/result/calc_stats.py /home/ziyanzha/MOM_within_gene/Whole_model/result/${filename}.txt