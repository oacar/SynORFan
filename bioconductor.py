import os
import re
import subprocess
import sys

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord


def muscle_align(input_seq):
    """
    This function uses muscle to align the input_seq whether it is dna or protein sequence
    using stdin and stdout. It does not write any output to file
    :param input_seq:
    :return: muscle aligned input_seq
    """
    muscle_cline = MuscleCommandline()
    child = subprocess.Popen(str(muscle_cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform != 'win32'))

    SeqIO.write(input_seq, child.stdin, 'fasta')
    child.stdin.close()

    aligned_seq = AlignIO.read(child.stdout, 'fasta')
    return (aligned_seq)


def translate_alignment(aln):
    """
    This function takes a dna sequence alignment (Bio.Align object)
    for every sequence in this alignment, removes '-' characters and translate.
    Stop codon is represented by 'X'
    Then this alignment is being aligned using muscle_align
    :param aln: Bio.Align object
    :return: aligned aa_translation as Bio.Align object
    """
    aa_recs = []
    for rec in aln:
        # print(rec.id)
        aa_recs.append(SeqRecord(rec.seq.ungap('-').translate(stop_symbol='X'), id=str(rec.id), description=''))
    aa_translation = muscle_align(aa_recs)  # AlignIO.read(child.stdout,'fasta')
    return aa_translation


def get_subalignment(input_seq, ref_seq_id, smorfSeq, outputDirectory, orf_name, is_aligned=False):
    """
    This function takes a series of sequences, aligns them and finds a sequence on the reference sequences.
    Then extracts subalignment and translate it. Saves these files.
    :param input_seq: Bio.Align object, Multiple Sequence alignment file
    :param ref_seq_id: int, the id of the reference sequence. This might change if is_aligned is False
    :param smorfSeq: the sequence that will be searched for in reference
    :param outputDirectory: string, the directory to write the subalignment and aa translation files
    :param orf_name: string, an identifier for outputs
    :param is_aligned: boolean, shows whether input_seq is already aligned or not
    :return: integers, start and end positions of smorfSeq on the alignment.
    """
    input_seqs = SeqIO.parse(input_seq, 'fasta')
    if is_aligned:
        align = input_seqs
    else:
        align = muscle_align(input_seqs)
        ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
    align_file = open(outputDirectory + '/' + orf_name + '_alignment_muscle.fa', 'w')
    AlignIO.write(align, align_file, 'fasta')

    # smorf = Seq('ATGTCCCGTG')

    re_smorf = re.compile("[-]*".join(smorfSeq))

    res = re_smorf.search(str(align[ref_seq_id].seq))
    start = res.start()
    end = res.end()

    subalign_seq = align[:, start:end + 100]

    subalign = muscle_align(subalign_seq)

    aa_translation = translate_alignment(subalign)
    subalign_file = open(outputDirectory + '/' + orf_name + '_subalignment.fa', 'w')
    AlignIO.write(subalign, subalign_file, 'fasta')

    aa_file = open(outputDirectory + '/' + orf_name + '_AATranslation.fa', 'w')
    AlignIO.write(aa_translation, aa_file, 'fasta')
    return (start, end)


def map_aln_to_seq(gapped):
    """
    This function maps positions of regular sequence to aligned, gap containing sequence
    :param gapped: Bio.Seq object. sequence which contains gaps as '-' characters
    :return: dictionary mapping un-gapped nucleotide positions to gapped positions
    """
    ungapped = gapped.ungap('-')
    j = 0
    map_ungapped_gapped = {}
    for i in range(len(ungapped)):
        if ungapped[i] == gapped[j]:
            map_ungapped_gapped[i] = j
            j += 1
        else:
            while gapped[j] == '-':
                j += 1
            map_ungapped_gapped[i] = j
            j += 1

    return map_ungapped_gapped


def find_orfs(seq, startCodons=['ATG'], stopCodons=['TAA', 'TAG', 'TGA']):
    """
    This function finds open reading frames in seq
    :param seq: Bio.Seq object or string, without gaps
    :param startCodons: as default 'ATG' is considered as start codon
    :param stopCodons: list of 'TAA','TAG','TGA' are considered as stop codons
    :return: np array containing start and stop positions of open reading frames
    """
    start_pos = [m.start() for m in re.finditer('|'.join(startCodons), str(seq))]
    stop_pos = [m.start() for m in re.finditer('|'.join(stopCodons), str(seq))]
    ms = np.array(np.meshgrid(start_pos, stop_pos)).T.reshape(-1, 2)
    orf_start_stop = ms[((ms[:, 1] - ms[:, 0]) % 3 == 0) & (ms[:, 0] < ms[:, 1])]
    return orf_start_stop


def do_ranges_overlap(l1, l2):
    """
    This function returns a boolean if 2 ranges given by [start,stop] overlaps
    :param l1: range 1, [start1, stop1]
    :param l2: range 2, [start2, stop2]
    :return: True if ranges overlap, False otherwise
    """
    if l1[0] <= l2[1] and l2[0] <= l1[1]:
        return True
    else:
        return False


def find_overlapping_orf_ranges(gapped, ref_start, ref_stop):
    """
    This function finds open reading frames which has an overlap with the range [ref_start, ref_stop]
    :param gapped: Bio.Seq object with gapped sequence
    :param ref_start: int, reference range start
    :param ref_stop: int, reference range stop
    :return: np.array if there are overlapping orfs, None otherwise
    """
    map_gapped = map_aln_to_seq(gapped)

    orf_ranges = find_orfs(gapped.ungap('-'))

    starts = [map_gapped[k] for k in orf_ranges[:, 0]]
    ends = [map_gapped[k] for k in orf_ranges[:, 1]]
    orfstartstop = np.array(list(zip(starts, ends)))

    l1 = [ref_start, ref_stop]
    oss = np.unique(orfstartstop, axis=0)
    c = []
    for i in range(len(oss)):
        l2 = oss[i, :]
        if do_ranges_overlap(l1, l2):
            c.append(i)

    if len(c) != 0:
        return np.unique(oss[c], axis=0)
    else:
        return None


def exclude_shorter_orfs(gapped, oor):
    """
    Name is a little less than what this function actually does.
    This function checks whether the orfs on gapped sequence given by oor np.array files contain
    a full orf. If there are stop codons in between, they are removed.
    Then takes longest open reading frame from the orfs that have stop codon at the same position.
    :param gapped: Bio.Seq object
    :param oor: np.array returned by find_overlapping_orf_ranges function
    :return: returns np.array containing orf ranges determined by above criteria
    """
    good_oor = []
    for i in range(len(oor)):
        start_ = oor[i, 0]
        end_ = oor[i, 1]
        sub_gapped = gapped[start_:(end_ + 3)]
        tr_seq = sub_gapped.seq.ungap('-').translate()
        if tr_seq.find('*') == len(tr_seq) - 1:
            good_oor.append(i)

    selected = oor[good_oor]
    # print(selected)
    # longer_oor = np.array([])
    xx = []
    yy = []
    for i in np.unique(selected[:, 1]):
        # print(np.array([max(selected[:, 0][selected[:, 1] == i]), i]))
        yy.append(i)
        xx.append(min(selected[:, 0][selected[:, 1] == i]))
        # res = np.array([, i])
        # np.append(res, longer_oor)
    longer_oor = np.stack((xx, yy), axis=1)
    return longer_oor


def get_frame_mapping(gapped, start):
    """
    This function finds translation frame positions with respect to start
    :param gapped: Bio.Seq object
    :param start: the ATG position that is the start codon of the ORF of interest
    :return: dictionary containing frame correspondence of each position
    """
    frame_mapping = {start: 0}
    for i in range(start - 1, -1, -1):
        if gapped[i] == '-':
            frame_mapping[i] = frame_mapping.get(i + 1)
        else:
            frame_mapping[i] = (frame_mapping.get(i + 1) - 1) % 3
    for i in range(start + 1, len(gapped), 1):
        if gapped[i] == '-':
            frame_mapping[i] = frame_mapping.get(i - 1)
        else:
            frame_mapping[i] = (frame_mapping.get(i - 1) + 1) % 3
    return frame_mapping


def get_correct_frame(gapped, range, start):
    frame_mapping = get_frame_mapping(gapped, start)

    range_start = range[0]
    frame = frame_mapping.get(range_start)
    while frame != 0:
        range_start = range_start - 1
        frame = frame_mapping.get(range_start)
    return range_start


def union_range(range1, range2):
    """

    :param range1:
    :param range2:
    :return: list containing the start and stop position of union of range1 and range2
    """
    start = min(range1[0], range2[0])
    end = max(range1[1], range2[1])

    return [start, end]


def intersect_range(range1, range2):
    """

    :param range1:
    :param range2:
    :return: list containing the start and stop position of intersection of range1 and range2
    """
    start = max(range1[0], range2[0])
    end = min(range1[1], range2[1])

    return [start, end]


def write_pairwise(aln, filename):
    """
    wrapper function for writing aln to filename
    :param aln: Bio.Align object
    :param filename: string for the filename of the output
    :return: None
    """
    aln_file = open(filename, 'w')
    AlignIO.write(aln, aln_file, 'fasta')
    aln_file.close()


def find_homologs(align, ref_seq_id, ref_range, orf_name, out_path='./'):
    """
    This function uses a multiple sequence alignment and finds orfs in all species except ref_seq_id
    which has overlaps with ref_range. Then it saves the pairwise alignment file that contains union and intersection of
    two orfs, one on the ref_seq_id and one on the comparison species.
    :param align: Bio.Align object with MSA
    :param ref_seq_id: int, the id for reference species
    :param ref_range: [start,stop] of ORF on the reference species
    :param orf_name: identifier for output writing
    :param out_path: output path for output writing
    :return: None
    """
    ref_seq = align[ref_seq_id, :].seq
    for i in range(len(align)):
        if i == ref_seq_id:
            continue

        pairwise_aln = MultipleSeqAlignment([align[i], align[ref_seq_id]])
        spec_name = align[i].id
        path = out_path + '/' + spec_name + '/'
        if os.path.isdir(path) is False:
            try:
                os.mkdir(path)
            except OSError:
                print("Creation of the directory %s failed" % path)
            else:
                print("Successfully created the directory %s " % path)

        oor = exclude_shorter_orfs(align[i, :],
                                   find_overlapping_orf_ranges(align[i, :].seq, ref_range[0], ref_range[1]))
        if oor is None:
            continue
        for u, itr in enumerate(oor):
            un_range = union_range(itr, ref_range)
            if itr[0] > ref_range[0]:
                un_range[0] = get_correct_frame(align[ref_seq_id, :].seq, un_range, ref_range[0])
            else:
                un_range[0] = get_correct_frame(align[i, :].seq, un_range, itr[0])

            int_range = intersect_range(itr, ref_range)
            if itr[0] < ref_range[0]:
                int_range[0] = get_correct_frame(align[ref_seq_id, :].seq, int_range, ref_range[0])
            else:
                int_range[0] = get_correct_frame(align[i, :].seq, int_range, itr[0])

            int_aln = muscle_align(pairwise_aln[:, int_range[0]:int_range[1]])
            write_pairwise(int_aln, path + orf_name + '_subalignment_overlap_' + str(u) + '.fa')

            int_trans = translate_alignment(int_aln)
            write_pairwise(int_trans, path + orf_name + '_AATranslation_overlap_' + str(u) + '.fa')

            uni_aln = muscle_align(pairwise_aln[:, un_range[0]:un_range[1]])
            write_pairwise(uni_aln, path + orf_name + '_subalignment_' + str(u) + '.fa')
            uni_trans = translate_alignment(uni_aln)
            write_pairwise(uni_trans, path + orf_name + '_AATranslation_' + str(u) + '.fa')
            orf_file = open(path + orf_name + '_orf_aa_' + str(u) + '.fa', 'w')
            SeqIO.write(SeqRecord(align[i][itr[0]:itr[1] + 3].seq.ungap('-').translate(stop_symbol = 'X'), id=align[i].id, description=''), orf_file,
                        'fasta')
            orf_file.close()


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
        if len(file_name) == 0 or len(file_name) == 1:
            return 0
        else:
            for i in range(len(file_name)):
                aln = AlignIO.read(path + '/' + file_name[i], 'fasta')
                count = count_identical_chars(aln[0].seq, aln[1].seq)
                if count > best_count:
                    best_count = count
                    best = i
            return best


def main():
    sub = AlignIO.read('data/YBR196C-A/YBR196C-A_subalignment.fa', 'fasta')

    start = 2754
    end = 2918
    align = AlignIO.read('data/ybr_deneme/YBR196C-A_alignment_muscle.fa', 'fasta')
    try:
        ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
    except IndexError:
        print('Reference sequence name is not in the alignment')
    # find_best_overlap_id('data/ybr_deneme/Spar')
    start, end = get_subalignment('data/YBR196C-A/YBR196C-A_alignment.fa', ref_seq_id, str(sub[0, :].seq.ungap('-')),
                                  'data/ybr_deneme', 'YBR196C-A')
    find_homologs(align=align, ref_seq_id=ref_seq_id, ref_range=[start, end], orf_name='YBR196C-A',
                  out_path='data/ybr_deneme')


if __name__ == '__main__':
    main()
