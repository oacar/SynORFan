import json
import pandas as pd 
from tqdm import tqdm
import argparse

def read_json_file(fname):
    with open(fname) as f:
        u = json.load(f)
    return u


def get_msa_data_to_csv(data_ex):
    msa_df = pd.DataFrame(columns=['orf_name', 'species','sequence','_type'])
    for _type in data_ex['msa']:
        d = data_ex['msa'][_type]
        name = d['name']
        seqs = d['sequences']
        species = d['species']
        for i in range(len(seqs)):
            msa_df = msa_df.append(pd.DataFrame(columns=['orf_name', 'species','sequence','_type'],
                        data= [[name,species[i], seqs[i] ,_type ]]))
    return msa_df


def get_pairwise_data_to_csv(data_ex):
    pairwise_df = pd.DataFrame(columns=['orf_name', 'species1', 'sequence1','species2','sequence2','type', 'id'])
    for group in data_ex['pairwise']:
        for species in data_ex['pairwise'][group]:
            for _id in data_ex['pairwise'][group][species]:
                d=data_ex['pairwise'][group][species][_id]
                name = d['name']
                species1 = d['species1']
                seq1= d['seq1']
                species2 = d['species2']
                seq2=d['seq2']
                _type = d['_type']
                pairwise_df = pairwise_df.append(
                    pd.DataFrame(
                        columns = ['orf_name', 'species1', 'sequence1','species2','sequence2','type', 'id'],
                        data= [[name, species1, seq1, species2, seq2 , _type, _id]]))
    return pairwise_df


def write_to_file(df,fname):
    df.to_csv(fname,index=False)


def write_to_database(df, engine,table_name):
    df.to_sql(table_name,engine,if_exists='replace',index=False)


def main(datafile, pairwise,output):#, db_username='root', db_password='dktgp275', db_ip='127.0.0.1',db_port=3315,db_name='deneme_synal'):
    json_file = read_json_file(datafile)
    if pairwise:
        orf = pd.DataFrame(columns = ['orf_name', 'sequence', 'species'])
        pairwise_df = pd.DataFrame(columns=['orf_name', 'species1', 'sequence1','species2','sequence2','type', 'id'])
        for i in tqdm(json_file):
            try:
                data_ex = json_file[i]
                df = pd.DataFrame(columns = ['orf_name', 'sequence', 'species'], 
                            data = [[data_ex['name'], data_ex['sequence'], data_ex['species']]])
                orf = orf.append(df)
                pairwise_df = pairwise_df.append(get_pairwise_data_to_csv(data_ex))
                #msa_df = msa_df.append(get_msa_data_to_csv(data_ex))
            except:
                print(i)
    #        continue
    
        write_to_file(orf,f'{output}_orf.csv')
        write_to_file(pairwise_df,f'{output}_pairwise.csv')
    else:
        msa_df = pd.DataFrame(columns=['orf_name', 'species','sequence','_type'])

        for i in tqdm(json_file):
            try:
                data_ex = json_file[i]
                # df = pd.DataFrame(columns = ['orf_name', 'sequence', 'species'], 
                #             data = [[data_ex['name'], data_ex['sequence'], data_ex['species']]])
                #orf = orf.append(df)
                #pairwise_df = pairwise_df.append(get_pairwise_data_to_csv(data_ex))
                msa_df = msa_df.append(get_msa_data_to_csv(data_ex))
            except:
                print(i)

    
        write_to_file(msa_df,f'{output}_msa.csv')

    #from sqlalchemy import create_engine
    #engine = create_engine(f'mysql+mysqldb://{db_username}:{db_password}@{db_ip}:{db_port}/{db_name}')
    #write_to_database(orf,engine,'orf')
    #write_to_database(pairwise_df,engine,'pairwise')
    #write_to_database(msa_df,engine,'msa')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', action="store", dest='datafile', help='json file for results',
                        required=True)
   # parser.add_argument('-n', action="store", dest='orf_name', help='ORF name for output names', required=True)
    parser.add_argument('-p', action='store_true', dest='pairwise', help='If given, saves pairwise and orf csv\'s',
                        default=False)
    parser.add_argument('-o', action="store", dest='output', help='json file output prefix',
                        required=True)
    res = parser.parse_args()
    main(res.datafile,res.pairwise, res.output)    
