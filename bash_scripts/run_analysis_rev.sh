#!/usr/bin/env bash

#for i in *;
ls | tac | while read i;
do
name="$(cut -d'_' -f1 <<<$i)";
#muscle_file_name = ""
if ls $i/*muscle* 1> /dev/null 2>&1;
then
#    echo "file exists"
    python ../../../bioconductor.py -p $i -n $name -y ../../orf_genomic_all.fasta -m
else
    python ../../../bioconductor.py -p $i -n $name -y ../../orf_genomic_all.fasta -ap
fi
#python ../../analysis.py -p $i -n $name
done;
