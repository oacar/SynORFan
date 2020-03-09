#!/usr/bin/env bash

for i in *;
do
name="$(cut -d'_' -f1,2,3,4 <<<$i)";
#echo $name
#python ../../bioconductor.py -p $i -n $name -y ../orf_genomic_all.fasta
python ../../../bioconductor.py -p $i -n $name -y $i/$name"_sequence.fa"

#python ../../analysis.py -p $i -n $name
done;
