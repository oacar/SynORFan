#!/usr/bin/env bash

for i in *;
do
name="$(cut -d'_' -f1 <<<$i)";
#muscle_file_name = ""
if ls $i/*subalignment* 1> /dev/null 2>&1;
then
#    echo "file exists"
#    python ../../../bioconductor.py -p $i -n $name -y ../../orf_genomic_all.fasta -m
    continue
else
    python ../../../bioconductor.py -p $i -n $name -y ../../orf_genomic_all.fasta -ap
fi

#python ../../analysis.py -p $i -n $name
done;
