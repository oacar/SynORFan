# SynORFan
This package is designed to find overlapping open reading frames in a multiple sequence alignment(MSA) given a reference alignment. Two main scripts are bioconductor.py and analysis.py which are both designed to be used as CLI programs. 

```
python bioconductor.py --help
```

will give you the options and input files needed to use this program.

```
usage: bioconductor.py [-h] -p PATH -n ORF_NAME [-a] -y YEAST [-m] [-ap]

optional arguments:
  -h, --help   show this help message and exit
  -p PATH      Directory path for alignment and output folder
  -n ORF_NAME  ORF name for output names
  -a           Is the sequence is annotated?
  -y YEAST     Fasta file containing dna sequence for annotated yeast genes
```


## Example Usage:
```
python bioconductor.py -p input_folder/ -n YBR196C-A -y orf_genomic_all.fasta -ap -a
```

## Requirements:
Python requirements are in `requirements.txt` file however you also need mafft to be on your system path and a tmp folder on your Home folder.(i.e. $HOME/tmp/ should be available) 
-y argument needs `orf_genomic_all.fasta` file which can be downloaded from SGD for yeast to get the sequence if -a is specified or the ORF sequence can be given directly to -y. 
