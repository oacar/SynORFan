import pandas as pd
import argparse
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
        orf_aa_file = [s for s in os.listdir(path) if 'orf_aa_' + str(id) in s][0]
        orf_aa = AlignIO.read(path + '/' + orf_aa_file, 'fasta')
    except (FileNotFoundError, TypeError) as error:
        return pd.DataFrame()
    else:
        # columns = ['A']
        df = pd.DataFrame()
        length = len(orf_aa[0].seq)
        if len(overlap_aa) == 1 or len(overlap_aa) == 0:
            print(path + '/' + overlap_aa_file + ' has no or alignment ')
            return df
        common_aa = count_identical_chars(overlap_aa[0].seq, overlap_aa[1].seq)
        df[orf_aa[0].id + "_length_" + str(id)] = [length]
        df[orf_aa[0].id + "_common_aa_" + str(id)] = [common_aa]

    return df


def subalignment_analysis(path, ref_id='Scer'):
    sub_aa_filename = [s for s in os.listdir(path) if '_AATranslation' in s][0]
    sub_aa = AlignIO.read(path + '/' + sub_aa_filename, 'fasta')
    df = pd.DataFrame()
    for i, record in enumerate(sub_aa):
        if (record.id == ref_id):
            continue
        df[record.id + '_start_codon'] = [record.seq.ungap('-')[0]]
        df[record.id + '_stop_codon'] = [record.seq.ungap('-')[-1]]

    return df


def main(path, orf_name):
    # path = 'data/pgs/YDR169C-A'
    # orf_name = 'YDR169C-A'
    df = pd.DataFrame()
    df['orf_name'] = [orf_name]
    df = df.join(subalignment_analysis(path))
    pairwise_folders = next(os.walk(path))[1]
    for folder in pairwise_folders:
        best_id = find_best_overlap_id(path + '/' + folder)
        if best_id is None:
            continue
        df[folder + '_best_id'] = [best_id]
        sub_split = list(set([u.split('_')[-1] for u in os.listdir(path + '/' + folder)]))
        list_of_ids = [int(i.split('.')[0]) for i in sub_split]

        if len(list_of_ids) == 0:
            continue
        for i in list_of_ids:
            df = df.join(pairwise_analyze(path + '/' + folder + '/', i))

    # return df
    df.to_csv(path + '/' + orf_name + '_data.csv')


#    subalignment_analysis(path)


if __name__ == '__main__':
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-p', action="store", dest='path', help='Directory path for alignment and output folder',
    #                    required=True)
    #parser.add_argument('-n', action="store", dest='orf_name', help='ORF name for output names', required=True)
    #res = parser.parse_args()
    #path = res.path
    #orf_name = res.orf_name

    #main(path, orf_name)
    main('data/pgs/YOR314W/','YOR314W')
