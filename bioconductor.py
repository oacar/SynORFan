import os
import argparse
import re
import subprocess
import sys

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from io import StringIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline, MafftCommandline
from Bio.Application import ApplicationError
from Bio.SeqRecord import SeqRecord
import analysis


def mafft_align(input_seq):
    """
    this uses mafft for alignment
    :param input_seq:
    :return:
    """
    tf = open('tmp.fa', 'w')
    SeqIO.write(input_seq, tf, 'fasta')
    mafft_cline = MafftCommandline(input='tmp.fa')
    tf.close()
    stdout, stderr = mafft_cline()
    if os.path.isfile('tmp.fa'):
        os.remove('tmp.fa')
    if len(stdout) == 0:
        raise ValueError(stderr)
    aligned_seq = AlignIO.read(StringIO(stdout), 'fasta')
    return aligned_seq


def muscle_align(input_seq):
    """
    This function uses muscle to align the input_seq whether it is dna or protein sequence
    using stdin and stdout. It does not write any output to file
    :param input_seq:
    :return: muscle aligned input_seq
    """
    muscle_cline = MuscleCommandline(maxiters=2, diags=True)
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
    Then this alignment is being aligned using mafft_align
    :param aln: Bio.Align object
    :return: aligned aa_translation as Bio.Align object
    """
    aa_recs = []
    for rec in aln:
        # print(rec.id)
        aa_recs.append(SeqRecord(rec.seq.ungap('-').translate(stop_symbol='X'), id=str(rec.id), description=''))
    aa_translation = mafft_align(aa_recs)  # AlignIO.read(child.stdout,'fasta')
    return aa_translation


def get_subalignment(input_seq, smorfSeq, outputDirectory, orf_name, is_aligned=False):
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
    if isinstance(input_seq, str):
        if is_aligned:
            input_seqs = AlignIO.read(input_seq, 'fasta')
            align = input_seqs
            ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]

        else:
            input_seqs = SeqIO.parse(input_seq, 'fasta')
            align = mafft_align(input_seqs)
            ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
            align_file = open(outputDirectory + '/' + orf_name + '_alignment_muscle.fa', 'w')
            AlignIO.write(align, align_file, 'fasta')
    else:
        input_seqs = input_seq
        align = mafft_align(input_seqs)
        ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
        other_name = [rec.id for rec in align if rec.id != 'Scer'][0]
        align_file = open(outputDirectory + '/' + orf_name + '_alignment_' + other_name + '.fa', 'w')
        AlignIO.write(align, align_file, 'fasta')
    # smorf = Seq('ATGTCCCGTG')

    re_smorf = re.compile("[-]*".join(smorfSeq), re.IGNORECASE)

    res = re_smorf.search(str(align[ref_seq_id].seq))
    if res is None:
        print(orf_name + ' is not found in alignment')
        raise ValueError
    start = res.start()
    end = res.end()
    sub_extended = True
    if isinstance(input_seq, str):
        subalign_seq = align[:, start:end]

        subalign = mafft_align(subalign_seq)

        aa_translation = translate_alignment(subalign)

        subalign_file = open(outputDirectory + '/' + orf_name + '_subalignment.fa', 'w')
        AlignIO.write(subalign, subalign_file, 'fasta')

        aa_file = open(outputDirectory + '/' + orf_name + '_AATranslation.fa', 'w')
        AlignIO.write(aa_translation, aa_file, 'fasta')
        if sub_extended:
            subalign_seq = align[:, max(start - 1000, 0):min(end + 1000, len(align[0]))]

            subalign = mafft_align(subalign_seq)

            subalign_file = open(outputDirectory + '/' + orf_name + '_subalignment_extended.fa', 'w')
            AlignIO.write(subalign, subalign_file, 'fasta')
            res = re_smorf.search(str(subalign_seq[ref_seq_id].seq))
            if res is None:
                print(orf_name + ' is not found in alignment')
                raise ValueError
            start = res.start()
            end = res.end()
    else:
        subalign_seq = align[:, start:end]

        subalign = mafft_align(subalign_seq)

        aa_translation = translate_alignment(subalign)

        subalign_file = open(outputDirectory + '/' + orf_name + '_subalignment_' + other_name + '.fa', 'w')
        AlignIO.write(subalign, subalign_file, 'fasta')

        aa_file = open(outputDirectory + '/' + orf_name + '_AATranslation_' + other_name + '.fa', 'w')
        AlignIO.write(aa_translation, aa_file, 'fasta')
        if sub_extended:
            subalign_seq = align[:, max(start - 1000, 0):min(end + 2000, len(align[0]))]

            subalign = mafft_align(subalign_seq)

            subalign_file = open(outputDirectory + '/' + orf_name + '_subalignment_extended_'+other_name+'.fa', 'w')
            AlignIO.write(subalign, subalign_file, 'fasta')
            res = re_smorf.search(str(subalign_seq[ref_seq_id].seq))
            if res is None:
                print(orf_name + ' is not found in alignment')
                raise ValueError
            start = res.start()
            end = res.end()
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
    start_pos = [m.start() for m in re.finditer('|'.join(startCodons), str(seq), re.IGNORECASE)]
    stop_pos = [m.start() for m in re.finditer('|'.join(stopCodons), str(seq), re.IGNORECASE)]
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

    orf_range = find_orfs(gapped.ungap('-'))
    if len(orf_range)==0:
        return None
    orf_ranges = exclude_shorter_orfs(gapped, orf_range)
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
    stopcodons = ['TAA', 'TAG', 'TGA']
    if oor is None:
        return None
    starts = np.unique(oor[:, 0])
    selected = np.empty([len(starts), 2])

    for i in range(len(starts)):
        start_ = starts[i]
        ranges = oor[oor[:, 0] == start_, :]
        widths = np.abs(ranges[:, 0] - ranges[:, 1])
        shorthest = ranges[np.where(widths == np.min(widths))]
        # end_ = oor[i, 1]
        # sub_gapped = gapped[start_:end_]
        # sub_ungapped = sub_gapped.seq.ungap('-')
        # stop_pos = [m.start() for m in re.finditer('|'.join(stopcodons), str(sub_ungapped), re.IGNORECASE)]
        # stop_pos_same_frame = [s for s in stop_pos if s % 3 == 0]
        # tr_seq = sub_gapped.seq.ungap('-').translate()
        # if len(stop_pos_same_frame) == 0:
        selected[i, :] = shorthest

    # selected = oor[good_oor]
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
    range_start_forward = range_start
    range_start_backward = range_start
    frame_forward = frame_mapping.get(range_start_forward)
    frame_backward = frame_mapping.get(range_start_backward)

    while frame_forward!=0 and frame_backward!=0:
        range_start_forward = range_start_forward+1
        range_start_backward = range_start_backward-1
        frame_forward = frame_mapping.get(range_start_forward)
        frame_backward = frame_mapping.get(range_start_backward)
        if range_start_backward == -1 and range_start_forward == len(gapped):
            raise ValueError
    # if range_start == 0:
    #     while frame != 0:
    #         range_start = range_start + 1
    #         frame = frame_mapping.get(range_start)
    #
    # else:
    #     while frame != 0:
    #         range_start = range_start - 1
    #         frame = frame_mapping.get(range_start)
    #         if range_start == -1:
    #             raise ValueError
    if frame_forward == 0:
        return range_start_forward
    elif frame_backward == 0:
        return range_start_backward
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
        # new_start = max(0, ref_range[0] - 1000)
        # new_stop = min(len(align[i]), ref_range[1] + 1000)

        aln_sub = pairwise_aln  # [:, new_start:new_stop]

        new_ref_start = ref_range[0]  # min(1000, ref_range[0])
        new_ref_stop = ref_range[1]  # new_ref_start + ref_range[1] - ref_range[0] + 1
        oor = find_overlapping_orf_ranges(aln_sub[0].seq,
                                          new_ref_start, new_ref_stop)
        new_ref_range = [new_ref_start, new_ref_stop]
        if oor is None:
            continue
        for u, itr in enumerate(oor):

            try:
                un_range = union_range(itr, new_ref_range)
                # correct for translation frames so that when translated reference orf has original aa sequence
                # create 2 range values one for the the original range calculated by union_range or intersect_range
                # and the other for frame correction

                if itr[0] < new_ref_range[0]:
                    new_un_range_start = get_correct_frame(align[ref_seq_id, :].seq, un_range, new_ref_range[0])
                    un_range = np.append([un_range], [[new_un_range_start, un_range[1]]], axis=0)
                elif itr[0] > new_ref_range[0]:
                    new_un_range_start = get_correct_frame(align[i, :].seq, un_range, itr[0])
                    un_range = np.append([[new_un_range_start, un_range[1]]], [un_range], axis=0)
                else:
                    un_range = np.append([un_range], [un_range], axis=0)

                int_range = intersect_range(itr, new_ref_range)
                if itr[0] > new_ref_range[0]:
                    new_int_range_start = get_correct_frame(align[ref_seq_id, :].seq, int_range, new_ref_range[0])
                    int_range = np.append([int_range], [[new_int_range_start, int_range[1]]], axis=0)
                elif itr[0] < new_ref_range[0]:
                    new_int_range_start = get_correct_frame(align[i, :].seq, int_range, itr[0])
                    int_range = np.append([[new_int_range_start, int_range[1]]], [int_range], axis=0)
                else:
                    int_range = np.append([int_range], [int_range], axis=0)
                int_aln_seqrecords = [align[0, int_range[0][0]:int_range[0][1]],
                                      align[1, int_range[1][0]:int_range[1][1]]]
                int_aln = mafft_align(int_aln_seqrecords)
                write_pairwise(int_aln, path + orf_name + '_subalignment_overlap_' + str(u) + '.fa')

                int_trans = translate_alignment(int_aln)
                write_pairwise(int_trans, path + orf_name + '_AATranslation_overlap_' + str(u) + '.fa')

                uni_aln_seqrecords = [align[0, un_range[0][0]:un_range[0][1]],
                                      align[1, un_range[1][0]:un_range[1][1]]]
                uni_aln = mafft_align(uni_aln_seqrecords)
                write_pairwise(uni_aln, path + orf_name + '_subalignment_' + str(u) + '.fa')
                uni_trans = translate_alignment(uni_aln)
                write_pairwise(uni_trans, path + orf_name + '_AATranslation_' + str(u) + '.fa')
                orf_file = open(path + orf_name + '_orf_aa_' + str(u) + '.fa', 'w')
                if len(int_trans) == 1 or len(int_aln) == 1 or len(uni_trans) == 1 or len(uni_aln) == 1:
                    raise ValueError
                SeqIO.write(
                    SeqRecord(align[i][itr[0]:itr[1]].seq.ungap('-').translate(stop_symbol='X'), id=align[i].id,
                              description=''), orf_file,
                    'fasta')
                orf_file.close()
            except (ApplicationError, ValueError) as error:
                print(orf_name + '_' + spec_name + '_' + str(u) + ' gave error. Removing its files')
                if os.path.exists(path + orf_name + '_subalignment_overlap_' + str(u) + '.fa'):
                    os.remove(path + orf_name + '_subalignment_overlap_' + str(u) + '.fa')
                if os.path.exists(path + orf_name + '_AATranslation_overlap_' + str(u) + '.fa'):
                    os.remove(path + orf_name + '_AATranslation_overlap_' + str(u) + '.fa')
                if os.path.exists(path + orf_name + '_subalignment_' + str(u) + '.fa'):
                    os.remove(path + orf_name + '_subalignment_' + str(u) + '.fa')
                if os.path.exists(path + orf_name + '_AATranslation_' + str(u) + '.fa'):
                    os.remove(path + orf_name + '_AATranslation_' + str(u) + '.fa')
                if os.path.exists(path + orf_name + '_orf_aa_' + str(u) + '.fa'):
                    os.remove(path + orf_name + '_orf_aa_' + str(u) + '.fa')


def main(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise):
    print(orf_name)
    #    start = 2754
    #    end = 2918

    #        path = 'data/pgs/YLL059C_2/'
    #        orf_name = 'YLL059C'
    #        yeast_fname = 'data/orf_genomic_all.fasta'
    #        is_annotated = True
    if is_aligned:
        filename = [s for s in os.listdir(path) if 'muscle.fa' in s][0]
    else:
        filename = [s for s in os.listdir(path) if '_alignment.fa' in s][0]
    # mcl = MuscleCommandline(input='data/ybr_deneme/YBR196C-A_alignment.fa',out = 'data/ybr_deneme/YBR196C-A_alignment_muscle.fa')
    # find_best_overlap_id('data/ybr_deneme/Spar')
    aln = SeqIO.parse(path + '/' + filename,'fasta')
    maxlen = 0
    for rec in aln:
        l = len(rec.seq)
        if l > maxlen:
            maxlen = l
    #    if maxlen>100000:
    #        return 0
    orf_seq = None
    if is_annotated:
        yeast = SeqIO.parse(yeast_fname, 'fasta')
        for record in yeast:
            if record.id == orf_name:
                orf_seq = record.seq
        if orf_seq is None:
            print(orf_name + ' is not found in ' + yeast_fname)
            return (0)
    else:
        yeast = SeqIO.parse(yeast_fname, 'fasta')
        for record in yeast:
            orf_seq = record.seq
    if align_pairwise:
        msa_file = list(SeqIO.parse(path + '/' + filename, 'fasta'))
        ref_seq_record = [rec for rec in msa_file if rec.id == 'Scer'][0]

        for record in msa_file:
            if record.id == 'Scer' or len(record.seq) == 0:
                continue

            start, end = get_subalignment([ref_seq_record, record], str(orf_seq),
                                          path, orf_name, is_aligned=is_aligned)
            aln_file_name = [s for s in os.listdir(path) if '_subalignment_extended_' + record.id in s][0]
            align = AlignIO.read(path + '/' + aln_file_name, 'fasta')
            try:
                ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
            except IndexError:
                print('Reference sequence name is not in the alignment')

            find_homologs(align=align, ref_seq_id=ref_seq_id, ref_range=[start, end], orf_name=orf_name,
                          out_path=path)
    else:
        start, end = get_subalignment(path + '/' + filename, str(orf_seq),
                                      path, orf_name, is_aligned=is_aligned)
        aln_file_name = [s for s in os.listdir(path) if '_alignment_muscle' in s][0]
        align = AlignIO.read(path + '/' + aln_file_name, 'fasta')
        try:
            ref_seq_id = [i for i, rec in enumerate(align) if rec.id == 'Scer'][0]
        except IndexError:
            print('Reference sequence name is not in the alignment')

        find_homologs(align=align, ref_seq_id=ref_seq_id, ref_range=[start, end], orf_name=orf_name,
                      out_path=path)
        # ss = []
    analysis.main(path, orf_name, yeast_fname, is_annotated, align_pairwise)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", dest='path', help='Directory path for alignment and output folder',
                        required=True)
    parser.add_argument('-n', action="store", dest='orf_name', help='ORF name for output names', required=True)
    parser.add_argument('-a', action='store_false', dest='is_annotated', help='Is the sequence is annotated?',
                        default=True)
    parser.add_argument('-y', action='store', dest='yeast',
                        help='Fasta file containing dna sequence for annotated yeast genes', required=True)
    parser.add_argument('-m', action='store_true', dest='is_aligned', help='is the input alignment is already aligned?')
    parser.add_argument('-ap', action='store_true', dest='align_pairwise', help='Use only pairwise alignments for faster analysis.')
    #    parser.add_argument('')
    res = parser.parse_args()
    path = res.path
    orf_name = res.orf_name
    yeast_fname = res.yeast
    is_annotated = res.is_annotated
    is_aligned = res.is_aligned
    align_pairwise = res.align_pairwise
    main(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise)

    #main('data/YPR204W/', 'YPR204W', 'data/orf_genomic_all.fasta', True, False, True)
    #main('data/0_10090_10398_0', '0_10090_10398_0', 'data/0_10090_10398_0/0_10090_10398_0_sequence.fa', False, False, True)

    #main('data/YBR196C-A/', 'YBR196C-A', 'data/orf_genomic_all.fasta', True, False, True)
    # analysis.main('data/YAL067C/', 'YAL067C','data/orf_genomic_all.fasta',True,True)
