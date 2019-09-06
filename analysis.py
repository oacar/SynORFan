import pandas as pd
import argparse
from Bio import AlignIO
import os
from Bio import SeqIO


def count_identical_chars(seq1, seq2):
    """
    Calculate identical characters between two sequences
    :param seq1: string or Bio.Seq object
    :param seq2: string or Bio.Seq object
    :return: int, number of identical characters
    """
    count = 0
    for i in range(len(seq1)):
        if seq1[i] == '-' and seq2[i] == '-':
            continue
        elif seq1[i] == seq2[i]:
            count += 1
    return count


def find_best_overlap_id(path):
    """
    This function read overlapping orf files that were written by find_homologs function and decides which is the best
    overlapping orf pair with highest number of identical amino acid pairs
    :param path: string, path which contains the overlapping protein files
    :return: identifier for best overlapping pair
    """
    try:
        file_name = [s for s in os.listdir(path) if 'AATranslation_overlap' in s]
    except FileNotFoundError:
        print("%s is not found" % path)
    else:
        best = None
        best_count = 0
        # best_length = 0
        if len(file_name) == 0:
            return None
        else:
            for i in range(len(file_name)):
                aln = AlignIO.read(path + '/' + file_name[i], 'fasta')
                id = file_name[i].split('_')[-1].split('.')[0]
                if len(aln) == 1 or len(aln) == 0:
                    print(path + '/' + file_name[i] + ' has no alignment ')
                    continue
                count = count_identical_chars(aln[0].seq, aln[1].seq)
                if count > best_count:
                    best_count = count
                    best = id
                # elif count == best_count &
            return best


def pairwise_analyze(path, id):
    try:
        overlap_aa_file = [s for s in os.listdir(path) if 'AATranslation_overlap_' + str(id) + '.fa' in s][0]
        overlap_aa = AlignIO.read(path + '/' + overlap_aa_file, 'fasta')
        overlap_dna_file = [s for s in os.listdir(path) if 'subalignment_overlap_' + str(id) + '.fa' in s][0]
        overlap_dna = AlignIO.read(path + '/' + overlap_dna_file, 'fasta')

        union_aa_file = [s for s in os.listdir(path) if 'subalignment_' + str(id) + '.fa' in s][0]
        union_aa = AlignIO.read(path + '/' + union_aa_file, 'fasta')

        union_dna_file = [s for s in os.listdir(path) if 'AATranslation_' + str(id) + '.fa' in s][0]
        union_dna = AlignIO.read(path + '/' + union_dna_file, 'fasta')
        orf_aa_file = [s for s in os.listdir(path) if 'orf_aa_' + str(id) + '.fa' in s][0]
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
        if record.id == ref_id:
            continue
        df[record.id + '_start_codon'] = [record.seq.ungap('-')[0]]
        df[record.id + '_stop_codon'] = [record.seq.ungap('-')[-1]]

    return df


def subalignment_dna_id(path, ref_id='Scer'):
    sub_dna_filename = [s for s in os.listdir(path) if '_subalignment.fa' in s][0]
    sub_dna = AlignIO.read(path + '/' + sub_dna_filename, 'fasta')
    df = pd.DataFrame()
    ref_seq_id = [i for i, rec in enumerate(sub_dna) if rec.id == ref_id][0]
    for i, record in enumerate(sub_dna):
        if i == ref_seq_id:
            continue
        identical_chars = count_identical_chars(record.seq, sub_dna[ref_seq_id].seq)
        df[record.id + '_dna_identity'] = [identical_chars / len(sub_dna[ref_seq_id].seq.ungap('-'))]

    return df


def main(path, orf_name, yeast_fname, is_annotated):
    # path = 'data/pgs/YDR169C-A'
    # orf_name = 'YDR169C-A'
    df = pd.DataFrame()
    df['orf_name'] = [orf_name]
    df = df.join(subalignment_analysis(path))
    df = df.join(subalignment_dna_id(path))
    if is_annotated:
        yeast = SeqIO.parse(yeast_fname, 'fasta')
        for record in yeast:
            if record.id == orf_name:
                orf_seq = record.seq
        if orf_seq is None:
            print(orf_name + ' is not found in ' + yeast_fname)
            return 0
    else:
        yeast = SeqIO.parse(yeast_fname, 'fasta')
        for record in yeast:
            orf_seq = record.seq
    df['orf_length'] = [len(orf_seq)]
    pairwise_folders = next(os.walk(path))[1]
    for folder in pairwise_folders:
        best_id = find_best_overlap_id(path + '/' + folder)
        if best_id is None:
            continue
        # df[folder + '_best_id'] = [best_id]
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
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", dest='path', help='Directory path for alignment and output folder',
                        required=True)
    parser.add_argument('-n', action="store", dest='orf_name', help='ORF name for output names', required=True)
    parser.add_argument('-a', action='store_false', dest='is_annotated', help='Is the sequence is annotated?',
                        default=True)
    parser.add_argument('-y', action='store', dest='yeast',
                        help='Fasta file containing dna sequence for annotated yeast genes', required=True)
    res = parser.parse_args()
    path = res.path
    orf_name = res.orf_name
    yeast_fname = res.yeast
    is_annotated = res.is_annotated
    main(path, orf_name, yeast_fname, is_annotated)
    # main('data/pgs/YOR314W/', 'YOR314W')
