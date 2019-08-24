import pandas as pd
import numpy as np
from Bio import AlignIO
import os
from Bio.SeqRecord import SeqRecord
from bioconductor import find_best_overlap_id
from bioconductor import count_identical_chars

def pairwise_analyze(path, id):
    try:
        overlap_aa_file = [s for s in os.listdir(path) if 'AATranslation_overlap_' + str(id) in s][0]
        overlap_aa = AlignIO.read(path + '/' + overlap_aa_file, 'fasta')
        overlap_dna_file = [s for s in os.listdir(path) if 'subalignment_overlap_' + str(id) in s][0]
        overlap_dna = AlignIO.read(path + '/' + overlap_dna_file, 'fasta')

        union_aa_file = [s for s in os.listdir(path) if 'subalignment_' + str(id) in s][0]
        union_aa = AlignIO.read(path + '/' + union_aa_file, 'fasta')

        union_dna_file = [s for s in os.listdir(path) if 'AATranslation_' + str(id) in s][0]
        union_dna = AlignIO.read(path + '/' + union_dna_file, 'fasta')
        orf_aa_file = [s for s in os.listdir(path) if 'orf_aa_'+str(id) in s][0]
        orf_aa = AlignIO.read(path + '/' + orf_aa_file, 'fasta')
    except (FileNotFoundError, TypeError) as error:
        raise error
    else:
        # columns = ['A']
        df = pd.DataFrame()
        length = len(orf_aa[0].seq)
        common_aa = count_identical_chars(overlap_aa[0].seq,overlap_aa[1].seq)

    return 0


def main():
    path = 'data/ybr_deneme/'
    pairwise_folders = next(os.walk(path))[1]
    for folder in pairwise_folders:
        for i in range(int(len(os.listdir('data/ybr_deneme/Sjur')) / 4)):
            pairwise_analyze(path + '/' + folder + '/', i)


if __name__ == '__main__':
    main()
