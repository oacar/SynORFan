import Bio
import os
import numpy as np
import readFiles
from tqdm import tqdm
import pandas as pd 
import argparse
from multiprocessing import Pool


def combine_sequences(seqs):
    from itertools import chain
    seqs = list(chain.from_iterable(seqs))
    seqs_scer = [u for u in seqs if u.id=='Scer']
    seqs_else = [u for u in seqs if u.id!='Scer']
    seqs_else.append(seqs_scer[0])
    from bioconductor import align_sequences

    seqs_aligned = align_sequences(seqs_else,algorithm='mafft')
    return seqs_aligned

def parallel_msa_read(args):
    df_ = args[0]
    orf_name = args[1]
    dir_path = args[2]
    seqs_dna  = []
    seqs_aa = []
    fls = df_.loc[df_['orf_name']==orf_name].iloc[:,0]
    #orf_name = orf_names_unique[i]
    orf = readFiles.Orf(orf_name, 'dummy_seq')

    for j in fls:
        try:
            path = f'{dir_path}/{j}/'
            s_dna = (readFiles.read_subalignment(path, True, ref_id='Scer',aa = False))
            seqs_dna.append([s_ for s_ in s_dna])
            s_aa = (readFiles.read_subalignment(path, True, ref_id='Scer',aa = True))
            seqs_aa.append([s_ for s_ in s_aa])
        except:
            print(f"error for {j}")
    orf.add_msa(readFiles.Msa(orf_name, combine_sequences(seqs_aa),_type='aa'))
    orf.add_msa(readFiles.Msa(orf_name, combine_sequences(seqs_dna),_type='dna'))
    #orf_dict[orf.__dict__['name']] = orf.__dict__
    return {orf.__dict__['name']: orf.__dict__}
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action="store", dest='path', help='Directory path for results folder',
                        required=True)
    parser.add_argument('-o', action='store', dest='output', required=True, help='Output json file name')
    parser.add_argument('-pa', action='store_true', dest='parallel',help='if given, pool of workers will be used for speed')

    #path = 'data/python_analysis/gene_all_1009/'
    res = parser.parse_args()
    fls = sorted(os.listdir(res.path))

    orf_names = [u.split('_')[0] for u in fls]

    d = dict(zip(fls, [u.split('_')[0] for u in fls]))
    orf_names_unique = np.unique(orf_names)


    df_ = pd.DataFrame(data={'fls':fls,'orf_name':orf_names})
    #df
    orf_dict = {}

    if res.parallel:
        with Pool(12) as p:
                #fls = sorted(os.listdir(path))
                max_ = len(orf_names_unique)
                args = [(df_, orf_name, res.path) for orf_name in orf_names_unique]
                with tqdm(total=max_) as pbar:
                    for i, _ in enumerate(p.imap_unordered(parallel_msa_read, args)):
                        orf_dict.update(_)
                        pbar.update()    
        
    else:
        for i in tqdm(range(len(orf_names_unique))):#range():
            seqs_dna  = []
            seqs_aa = []
            fls = df_.loc[df_['orf_name']==orf_names_unique[i]].iloc[:,0]
            orf_name = orf_names_unique[i]
            orf = readFiles.Orf(orf_name, 'dummy_seq')

            for j in fls:
                try:
                    path = f'{res.path}/{j}/'
                    s_dna = (readFiles.read_subalignment(path, True, ref_id='Scer',aa = False))
                    seqs_dna.append([s_ for s_ in s_dna])
                    s_aa = (readFiles.read_subalignment(path, True, ref_id='Scer',aa = True))
                    seqs_aa.append([s_ for s_ in s_aa])
                except:
                    print(f"error for {j}")
            orf.add_msa(readFiles.Msa(orf_name, combine_sequences(seqs_aa),_type='aa'))
            orf.add_msa(readFiles.Msa(orf_name, combine_sequences(seqs_dna),_type='dna'))
            orf_dict[orf.__dict__['name']] = orf.__dict__
            
    import json
    with open(res.output, 'w') as outfile:
        json.dump(orf_dict
            , outfile)

if __name__ == '__main__':
    main()