from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord

from Bio.Alphabet import IUPAC
import re

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


input_seqs = SeqIO.parse('data/YBR196C-A/YBR196C-A_alignment.fa', 'fasta')
subalign = AlignIO.read('data/YBR196C-A/YBR196C-A_subalignment.fa', 'fasta')




align = muscle_align(input_seqs)#AlignIO.read(child.stdout,'fasta')
#print(align)

scer_id = [i for i, rec in enumerate(align) if rec.id=='Scer'][0]

smorf = Seq('ATGTCCCGTG')

re_smorf = re.compile("[-]*".join(smorf))

res = re_smorf.search(str(align[3].seq))
start = res.start()
end = res.end()

subalign_seq = align[:,start:end+100]


subalign = muscle_align(subalign_seq)

aa_recs = []
for rec in subalign:
    print(rec.id)
    aa_recs.append(SeqRecord(rec.seq.ungap('-').translate(), id=str(rec.id)))

aa_translation = muscle_align(aa_recs)#AlignIO.read(child.stdout,'fasta')

