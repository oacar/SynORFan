import pandas as pd
import argparse
from Bio import AlignIO
import os
from Bio import SeqIO
import json
from tqdm import tqdm


from multiprocessing import Pool
class Orf(object):
    def __init__(self, name, sequence, species='Scer'):
        self.name = name
        self.sequence = sequence
        self.species = species
        self.pairwise = {'union_aa': {}, 'union_dna':{}, 'overlap_aa':{},'overlap_dna':{}}
        self.msa = {}
    def add_pairwise(self, pw_obj):
        if pw_obj.species2 not in self.pairwise[pw_obj._type]:
            self.pairwise[pw_obj._type][pw_obj.species2] = {}
            self.pairwise[pw_obj._type][pw_obj.species2][pw_obj._id] = pw_obj.__dict__
        else:
            self.pairwise[pw_obj._type][pw_obj.species2][pw_obj._id] = pw_obj.__dict__

    def add_msa(self, msa_obj):
        self.msa[msa_obj._type] = msa_obj.__dict__
        
class Pairwise(object):
    def __init__(self, name, species1, seq1, species2, seq2, _type ='dna', _id=1):
        self.name = name
        self.species1 = species1
        self.seq1 = seq1
        self.species2 = species2
        self.seq2 = seq2
        self._type=_type
        self._id = _id
class Msa(object):
    def __init__(self, name, biomsa, _type='dna'):
        self.name = name
        order=['Scer','Spar', 'Smik','Sjur','Skud','Suva','Seub','Sarb']
        order = [o for o in order if o in [s.id for s in biomsa]]
        self.species = order
        self.sequences = [str(s.seq) for s in biomsa]
        res = {} 
        for key in [s.id for s in biomsa]: 
            for value in self.sequences: 
                res[key] = value 
                self.sequences.remove(value) 
                break
        result_list = [res[x] for x in  self.species]
        self.sequences=result_list
        
        self._type = _type

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

def add_pairwise_from_dict(p, name, _type, _id):
    aln = p[_type]
    pw = Pairwise(name,aln[0].id , str(aln[0].seq), aln[1].id, str(aln[1].seq), _type, _id)
    return pw

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


def pairwise_read(path, id):
    try:
        overlap_aa_file = [s for s in os.listdir(path) if 'AATranslation_overlap_' + str(id) + '.fa' in s][0]
        overlap_aa = AlignIO.read(path + '/' + overlap_aa_file, 'fasta')
        overlap_dna_file = [s for s in os.listdir(path) if 'subalignment_overlap_' + str(id) + '.fa' in s][0]
        overlap_dna = AlignIO.read(path + '/' + overlap_dna_file, 'fasta')

        union_dna_file = [s for s in os.listdir(path) if 'subalignment_' + str(id) + '.fa' in s][0]
        union_aa_file = [s for s in os.listdir(path) if 'AATranslation_' + str(id) + '.fa' in s][0]

        union_aa = AlignIO.read(path + '/' + union_aa_file, 'fasta')

        union_dna = AlignIO.read(path + '/' + union_dna_file, 'fasta')
        orf_aa_file = [s for s in os.listdir(path) if 'orf_aa_' + str(id) + '.fa' in s][0]
        orf_aa = AlignIO.read(path + '/' + orf_aa_file, 'fasta')
    except (FileNotFoundError, TypeError) as error:
        return None
    else:
        # columns = ['A']
   #     df = pd.DataFrame()
    #    length = len(orf_aa[0].seq)
        if len(overlap_aa) == 1 or len(overlap_aa) == 0:
            print(path + '/' + overlap_aa_file + ' has no or alignment ')
            return None
#        common_aa = count_identical_chars(overlap_aa[0].seq, overlap_aa[1].seq)
        return {'union_aa': union_aa, 'overlap_aa':overlap_aa, 'union_dna':union_dna, 'overlap_dna':overlap_dna}

 #       df[orf_aa[0].id + "_length_" + str(id)] = [length]
  #      df[orf_aa[0].id + "_common_aa_" + str(id)] = [common_aa]



def read_subalignment(path, align_pairwise, ref_id='Scer', aa=True):
    if aa:
        search_str = '_AATranslation'
    else:
        search_str = '_subalignment'
    if align_pairwise:
        sub_dna_filenames = [s for s in os.listdir(path) if search_str in s and 'extended' not in s]
        seqs = []
        for fname in sub_dna_filenames:
            sub_dna = AlignIO.read(path + '/' + fname, 'fasta')
            ref_seq_id = [i for i, rec in enumerate(sub_dna) if rec.id == ref_id][0]
            for i, record in enumerate(sub_dna):
                if i == ref_seq_id or len(record.seq.ungap('-')) == 0:
                    continue
                seqs.append(record)
                #identical_chars = count_identical_chars(record.seq, sub_dna[ref_seq_id].seq)
                #df[record.id + '_dna_identity'] = [identical_chars / len(sub_dna[ref_seq_id].seq.ungap('-'))]
        seqs.append(sub_dna[ref_seq_id])
        from bioconductor import align_sequences

        seqs_aligned = align_sequences(seqs,algorithm='mafft')


    
    else:
        sub_dna_filename = [s for s in os.listdir(path) if search_str+'.fa' in s][0]
        sub_dna = AlignIO.read(path + '/' + sub_dna_filename, 'fasta')
        #df = pd.DataFrame()
        #ref_seq_id = [i for i, rec in enumerate(sub_dna) if rec.id == ref_id][0]
        seqs_aligned= sub_dna
        #for i, record in enumerate(sub_dna):
         #   if i == ref_seq_id or len(record.seq.ungap('-')) == 0:
          #      continue
            #identical_chars = count_identical_chars(record.seq, sub_dna[ref_seq_id].seq)
            #df[record.id + '_dna_identity'] = [identical_chars / len(sub_dna[ref_seq_id].seq.ungap('-'))]


    return seqs_aligned
    

# def order_seqs(seqs, order=['Scer','Spar', 'Smik','Sjur','Skud','Seub','Sbay','Sarb']):
#     seqs_ordered=[]
#     names = [i.id for i in seqs]
#     order = order[order in names]
#     for i in seqs:
#         if

def parallel_func(args):
    path = args[0]
    f = args[1]
    yeast_fname=args[2]
    is_annotated=args[3]
    align_pairwise=args[4]
    orf_name = f.split('_')[0]
    return {f:readOrfData(path+'/'+f,orf_name,yeast_fname,is_annotated,align_pairwise)}

def main(path, yeast_fname, is_annotated, align_pairwise,output,parallel):
    df_dict = {}
    count = 0 
    #parallel=False
   # [u.split('_')[-1] for u in os.listdir(path + '/' + folder)])
    if parallel:
        fls = os.listdir(path)
 
        with Pool(12) as p:
                fls = sorted(os.listdir(path))
                max_ = len(fls)
                args = [(path, f, yeast_fname, is_annotated, align_pairwise) for f in fls]
                with tqdm(total=max_) as pbar:
                    for i, _ in enumerate(p.imap_unordered(parallel_func, args)):
                        df_dict.update(_)
                        pbar.update()    
    else:
        for f in os.listdir(path):
            orf_name = f.split('_')[0]
            try:
                df_dict[f] = readOrfData(path+'/'+f, orf_name, yeast_fname, is_annotated, align_pairwise)
            except:
                print(f)
            count = count + 1
            if count%100==0:
                print(count)
    import json
    with open(output, 'w') as outfile:
        json.dump(df_dict
              , outfile)

def readOrfData(path, orf_name, yeast_fname, is_annotated, align_pairwise):
    # path = 'data/pgs/YDR169C-A'
    # orf_name = 'YDR169C-A'
    # df = pd.DataFrame()
    # df['orf_name'] = [orf_name]
   # subalignment_analysis(path, align_pairwise))
    #subalignment_dna_id(path, align_pairwise))

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
    orf = Orf(orf_name, str(orf_seq))

    #df['orf_length'] = [len(orf_seq)]
    pairwise_folders = next(os.walk(path))[1]
    for folder in pairwise_folders:
    #    best_id = find_best_overlap_id(path + '/' + folder)
     #   if best_id is None:
      #      continue
        # df[folder + '_best_id'] = [best_id]
        sub_split = list(set([u.split('_')[-1] for u in os.listdir(path + '/' + folder)]))
        list_of_ids = [int(i.split('.')[0]) for i in sub_split]

        if len(list_of_ids) == 0:
            continue
        for i in list_of_ids:
            p = pairwise_read((path + '/' + folder + '/'),i)
           # print(path)
            union_aa = add_pairwise_from_dict(p, orf_name, 'union_aa', i)
            union_dna = add_pairwise_from_dict(p, orf_name, 'union_dna', i)
            overlap_aa = add_pairwise_from_dict(p, orf_name, 'overlap_aa', i)
            overlap_dna = add_pairwise_from_dict(p, orf_name, 'overlap_dna', i)
            orf.add_pairwise(union_aa)
            orf.add_pairwise(union_dna)
            orf.add_pairwise(overlap_aa)
            orf.add_pairwise(overlap_dna)
            #df = df.join(pairwise_analyze(path + '/' + folder + '/', i))
    msa_dna = Msa(orf_name, read_subalignment(path, align_pairwise, aa=False), _type = 'dna')
    msa_aa = Msa(orf_name, read_subalignment(path, align_pairwise, aa=True), _type='aa')
    orf.add_msa(msa_dna)
    orf.add_msa(msa_aa)

    return(orf.__dict__)
    # return df
    #df.to_csv(path + '/' + orf_name + '_data.csv')


#    subalignment_analysis(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", dest='path', help='Directory path for results folder',
                        required=True)
   # parser.add_argument('-n', action="store", dest='orf_name', help='ORF name for output names', required=True)
    parser.add_argument('-a', action='store_false', dest='is_annotated', help='Is the sequence is annotated?',
                        default=False)
    parser.add_argument('-y', action='store', dest='yeast',
                        help='Fasta file containing dna sequence for annotated yeast genes', required=True)
    parser.add_argument('-ap', action='store_true', dest='align_pairwise')
    parser.add_argument('-o', action='store', dest='output', required=True, help='Output json file name')
    parser.add_argument('-pa', action='store_true', dest='parallel',help='if given, pool of workers will be used for speed')
    res = parser.parse_args()
    path = res.path
    #orf_name = res.orf_name
    yeast_fname = res.yeast
    is_annotated = res.is_annotated
    align_pairwise = res.align_pairwise
    main(path, yeast_fname, is_annotated, align_pairwise,res.output, res.parallel)

    # main('data/pgs/YOR314W/', 'YOR314W')
