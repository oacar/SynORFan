from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
import re
import numpy as np
import subprocess
import sys


def muscle_align(input_seq):
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
    return(aligned_seq)

def getSubalignment(input_seq, smorfSeq,outputDirectory,orf_id):
    input_seqs = SeqIO.parse(input_seq, 'fasta')
    align = muscle_align(input_seqs)
    align_file = open(outputDirectory + '/' + orf_id + '_alignment_muscle.fa', 'w')
    AlignIO.write(align, align_file, 'fasta')
    scer_id = [i for i, rec in enumerate(align) if rec.id=='Scer'][0]

    #smorf = Seq('ATGTCCCGTG')

    re_smorf = re.compile("[-]*".join(smorfSeq))

    res = re_smorf.search(str(align[scer_id].seq))
    start = res.start()
    end = res.end()

    subalign_seq = align[:,start:end+100]


    subalign = muscle_align(subalign_seq)

    aa_recs = []
    for rec in subalign:
        #print(rec.id)
        aa_recs.append(SeqRecord(rec.seq.ungap('-').translate(), id=str(rec.id)))

    aa_translation = muscle_align(aa_recs)#AlignIO.read(child.stdout,'fasta')
    subalign_file = open(outputDirectory+'/'+orf_id+'_subalignment.fa','w')
    AlignIO.write(subalign,subalign_file,'fasta')

    aa_file = open(outputDirectory+'/'+orf_id+'_AATranslation.fa','w')
    AlignIO.write(aa_translation,aa_file,'fasta')
    return (start,end)
sub = AlignIO.read('data/YBR196C-A/YBR196C-A_subalignment.fa', 'fasta')

start, end = getSubalignment('data/YBR196C-A/YBR196C-A_alignment.fa',str(sub[0,:].seq.ungap('-')),'data/ybr_deneme','YBR196C-A')

start = 2754
end = 2918
seqId = 3
align = AlignIO.read('data/ybr_deneme/YBR196C-A_alignment_muscle.fa','fasta')
def map_aln_to_seq(gapped):
    ungapped = gapped.ungap('-')
    j=0
    map_ungapped_gapped = {}
    for i in range(len(ungapped)):
        if ungapped[i] == gapped[j]:
            map_ungapped_gapped[i]=j
            j+=1
        else:
            while gapped[j] == '-':
                j+=1
            map_ungapped_gapped[i] = j
            j+=1

    return map_ungapped_gapped


def findORFs(seq, startCodons=['ATG'], stopCodons=['TAA','TAG','TGA']):
    start_pos = [m.start() for m in re.finditer('|'.join(startCodons), str(seq))]
    stop_pos = [m.start() for m in re.finditer('|'.join(stopCodons), str(seq))]
    ms = np.array(np.meshgrid(start_pos,stop_pos)).T.reshape(-1,2)
    orf_start_stop = ms[((ms[:,1]-ms[:,0])%3==0 )& (ms[:,0]<ms[:,1])]
    return orf_start_stop

def isRangesOverlap(l1 , l2):
    if l1[0]<=l2[1] and l2[0]<=l1[1]:
        return True
    else:
        return False

def findOverlappingOrfRanges(gapped, ref_start,ref_stop):
    map_gapped = map_aln_to_seq(gapped)

    orf_ranges = findORFs(gapped.ungap('-'))

    starts = [map_gapped[k] for k in orf_ranges[:, 0]]
    ends = [map_gapped[k] for k in orf_ranges[:, 1]]
    orfstartstop = np.array(list(zip(starts, ends)))

    l1 = [ref_start,ref_stop]
    oss = np.unique(orfstartstop , axis = 0)
    c = []
    for i in range(len(oss)):
        l2 = oss[i, :]
        if isRangesOverlap(l1, l2):
            c.append(i)

    if len(c)!=0:
        return np.unique(oss[c],axis=0)
    else:
        return None


oor = findOverlappingOrfRanges(align[seqId, :].seq, start,end+100)


def excludeShorterOrfs(gapped,oor):
    good_oor = []
    for i in range(len(oor)):
        start_ = oor[i, 0]
        end_ = oor[i, 1]
        sub_gapped = gapped[start_:(end_ + 3)]
        tr_seq = sub_gapped.seq.ungap('-').translate()
        if tr_seq.find('*') == len(tr_seq) - 1:
            good_oor.append(i)

    selected = oor[good_oor]
    #print(selected)
    #longer_oor = np.array([])
    xx = []
    yy = []
    for i in np.unique(selected[:, 1]):
        #print(np.array([max(selected[:, 0][selected[:, 1] == i]), i]))
        yy.append(i)
        xx.append(min(selected[:, 0][selected[:, 1] == i]))
        #res = np.array([, i])
        #np.append(res, longer_oor)
    longer_oor = np.stack((xx,yy),axis = 1)
    return longer_oor


good_oor = excludeShorterOrfs(align[seqId,:],oor)
for i in range(2):
    print(align[seqId,good_oor[i,0]:good_oor[i,1]+3])

min_pro_len = 10
#
# def find_orfs_with_trans(seq, min_protein_length):
#     answer = []
#     seq_len = len(seq)
#     for strand, nuc in [(+1, seq)]:
#         for frame in range(3):
#             trans = str(nuc[frame:].translate())
#             trans_len = len(trans)
#             aa_start_pos = [m.start() for m in re.finditer('M', str(trans))]
#             aa_end = 0
#             aa_start = aa_start_pos[0]
#             while aa_start < trans_len:
#                 aa_end = trans.find("*", aa_start)
#                 if aa_end == -1:
#                     aa_end = trans_len
#                 if aa_end-aa_start >= min_protein_length:
#                     if strand == 1:
#                         start = frame+aa_start*3
#                         end = min(seq_len,frame+aa_end*3+3)
#                     answer.append((start, end, strand,
#                                    trans[aa_start:aa_end]))
#                 aa_start = min(itr for itr in aa_start_pos if itr>=aa_end+1 & itr is not None)
#     answer.sort()
#     return answer
#
# orf_list = find_orfs_with_trans(align[0,:].seq.ungap('-'), min_pro_len)
# for start, end, strand, pro in orf_list:
#     if pro[0]=='M':
#         print("%s...%s - length %i, strand %i, %i:%i" \
#               % (pro[:30], pro[-3:], len(pro), strand, start, end))
#
# min_val = min(i for i in test_list if i > k)
