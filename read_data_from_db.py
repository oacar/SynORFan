import json
import pandas as pd 
from tqdm import tqdm
import argparse
from Bio import SeqIO
from Bio import AlignIO
from sqlalchemy import create_engine
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import numpy as np
def writeMsa(df,path,orf_name):
    for t in ['aa','dna']:
        subset = df.loc[df['_type']==t]
        sr = [SeqRecord(Seq(subset.iloc[i,:].loc['sequence']),id=subset.iloc[i,:].loc['species'],description="") for i in range(subset.shape[0])]
        if t == 'aa':
            SeqIO.write(sr, f"{path}/{orf_name}_AATranslation.fa", "fasta")
        else:
            SeqIO.write(sr, f"{path}/{orf_name}_subalignment.fa", "fasta")
        #print(sr)
    return 0

def writePairwise(df,path,orf_name):
    species = df.species2.unique()
    for sp in species:
        for t in ['union_aa','union_dna','overlap_aa','overlap_dna']:
            subset = df.loc[(df['type']==t) & (df['species2']==sp)]
            #print((subset.iloc[0]['sequence1']))
            sr = [SeqRecord(Seq(subset.iloc[0]['sequence1']),id='Scer',description=''),SeqRecord(Seq(subset.iloc[0]['sequence2']),id=sp,description='')]
            if t=='union_dna':
                SeqIO.write(sr, f"{path}/{sp}/{orf_name}_subalignment.fa", "fasta")
            elif t == 'union_aa':
                SeqIO.write(sr, f"{path}/{sp}/{orf_name}_AATranslation.fa", "fasta")
            elif t == 'overlap_dna':
                SeqIO.write(sr, f"{path}/{sp}/{orf_name}_subalignment_overlap.fa", "fasta")
            elif t == 'overlap_aa':
                SeqIO.write(sr, f"{path}/{sp}/{orf_name}_AATranslation_overlap.fa", "fasta")
                #SeqIO.write(sr, f"{path}_{sp}_{t}.fa", "fasta")
            
    return 0
def main(path,orf_name):
    #engine = create_engine('mysql+mysqldb://root:dktgp275@127.0.0.1:3315/deneme_synal')
    engine_paris = create_engine('mysql+mysqldb://oma21:dktgp2750@127.0.0.1:3306/omer')
    print(orf_name)
    msa_query = f"""select * from msa where orf_name='{orf_name}'"""
    pairwise_query=f"""select * from pairwise_best where orf_name='{orf_name}'"""
    pairwise_data = pd.read_sql_query(pairwise_query, engine_paris)
    msa_data = pd.read_sql_query(msa_query, engine_paris)

    #print(msa_data)
    #print(pairwise_data)
    species = np.union1d(pairwise_data.species2.unique(),msa_data.species.unique())
    for s in species:
        if os.path.isdir(path+'/'+s) is False:
            try:
                os.mkdir(path+'/'+s)
            except OSError:
                print("Creation of the directory %s failed" % path)
            else:
                print("Successfully created the directory %s " % path)
    writeMsa(msa_data,path,orf_name)
    writePairwise(pairwise_data,path,orf_name)
    engine_paris.dispose()
    return 0


if __name__=='__main__':
    main('read_db_deneme','YGR146C-A')