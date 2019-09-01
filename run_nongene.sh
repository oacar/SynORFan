#!/usr/bin/env bash

for i in *;
do
name="$(cut -d'_' -f1 <<<$i)";
#python ../../bioconductor.py -p $i -n $name -y ../orf_genomic_all.fasta
python ../../bioconductor.py -p $i -n $i -y $i/$i"_sequence.fa" -a
#python ../../analysis.py -p $i -n $name
done;
