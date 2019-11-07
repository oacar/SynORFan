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
  -m           is the input alignment is already aligned?
  -ap          Use only pairwise alignments for faster analysis.
```


## Example Usage:
```
python bioconductor.py -p input_folder/ -n YBR196C-A -y orf_genomic_all.fasta -ap  
