#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

cat /home/ziyanzha/MOM_within_gene/AIREML/result/${filename}/rep*.txt > /home/ziyanzha/MOM_within_gene/AIREML/result/${filename}.txt
