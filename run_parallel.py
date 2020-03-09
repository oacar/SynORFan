import numpy as np
from bioconductor import use_pairwise_extended
import os
import re
import argparse
from tqdm import tqdm


from multiprocessing import Pool

# def f(x):
#     return x*x

#if __name__ == '__main__':

def parallel_func(i):
    p='data/python_analysis/gene_all_0228_2/'
    path=p+'/'+i
    orf_name = i.split('_')[0]
    yeast_fname = 'data/orf_genomic_all.fasta'
    is_annotated= True
    is_aligned= False
    align_pairwise = True
    algorithm = 'mafft'
    if len([s for s in os.listdir(path) if 'extended' in s])!=0:
        try:
            use_pairwise_extended(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise, algorithm=algorithm)
        except:
            print(f"{path} has extended file but gave error")
    else:
        print(f"{path} don't have extended files")
        
def main(p):
    #p = 'data/python_analysis/gene_all_0228/'
    fls = sorted(os.listdir(p))
    for i in tqdm(fls):
        path=p+'/'+i
        orf_name = i.split('_')[0]
        yeast_fname = 'data/orf_genomic_all.fasta'
        is_annotated= True
        is_aligned= False
        align_pairwise = True
        algorithm = 'mafft'
        if len([s for s in os.listdir(path) if 'extended' in s])!=0:
            try:
                use_pairwise_extended(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise, algorithm=algorithm)
            except:
                print(f"{path} has extended file but gave error")
        else:
            print(f"{path} don't have extended files")

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-p', action="store", dest='path', help='Directory path for alignment and output folder',
    #                     required=True)
    # res = parser.parse_args()
    path = 'data/python_analysis/gene_all_0228_2/'#res.path
    with Pool(12) as p:
        fls = sorted(os.listdir(path))
        max_ = len(fls)
        with tqdm(total=max_) as pbar:
            for i, _ in enumerate(p.imap_unordered(parallel_func, fls)):
                pbar.update()
        #p.map(parallel_func, fls)

#print(algorithm)
# if False:
#     path = 'tmp/YIL089W'
#     orf_name = 'YIL089W'
#     yeast_fname = 'data/orf_genomic_all.fasta'
#     is_annotated= True
#     is_aligned= False
#     align_pairwise = True
#     algorithm = 'mafft'
#     #main(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise, algorithm=algorithm)
#     use_pairwise_extended(path, orf_name, yeast_fname, is_annotated, is_aligned, align_pairwise, algorithm=algorithm)
