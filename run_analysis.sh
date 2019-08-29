#!/usr/bin/env bash

for i in *;
do
name="$(cut -d'_' -f1 <<<$i)";
python ../../analysis.py -p $i -n $name
done;
