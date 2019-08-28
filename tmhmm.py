import os
import re
import subprocess
import sys
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_tmhmm(input_seq):
    """
    this function runs TMHMM command line program and return TM helix start and stop if they exists
    :param input_seq: SeqRecord or AlignIO object
    :return: numpy array with start and stop positions
    """
    tmhmm_cline = 'tmhmm'
    res = np.empty((0, 2))
    child = subprocess.Popen(str(tmhmm_cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform != 'win32'))

    SeqIO.write(input_seq, child.stdin, 'fasta')
    child.stdin.close()
    # print(child.stdout.read())o

    while True:
        line = child.stdout.readline()
        if not line:
            break
        if "TMhelix" in line:
            ll = re.split('\s+', line)
            res = np.append(res, [[int(ll[-3]), int(ll[-2])]], axis=0)

    child.stdout.close()
    child.stderr.close()
    child.wait()
    return res


#seq = AlignIO.read("data/YBR196C-A/Spar/YBR196C-A_AATranslation_Spar_best.fa", 'fasta')
#res = run_tmhmm(SeqRecord(Seq('MRRRRRMRRRR'), id='s'))
#print(res)
