import pandas as pd
import numpy as np

import os

DIR='/Users/benceszalai/Downloads/why-pathway-methods-work/results/'+\
    'benchmark/scores/'

k=0
for setname1 in ['BIOCARTA','CGP','KEGG','REACTOME','dorothea_AB']:
    for setname2 in ['BIOCARTA','CGP','KEGG','REACTOME','dorothea_AB']:
        if setname1<setname2:
            for set_type in ['minus','minus_rev','overlap']:
                k+=1
                print(k)
                if set_type=='minus_rev':
                    fname='_'.join([setname2,'minus',setname1])
                else:
                    fname='_'.join([setname1,set_type,setname2])
                progeny=pd.read_csv(DIR+'progeny/'+fname+'_0.csv',
                                    sep=',',header=0,index_col=0)
                tcga=pd.read_csv(DIR+'tcga/'+fname+'_0.csv',
                                    sep=',',header=0,index_col=0)
                os.remove(DIR+'progeny/'+fname+'_0.csv')
                os.remove(DIR+'tcga/'+fname+'_0.csv')
                i=0
                while True:
                    i+=1
                    try:
                        temp=pd.read_csv(DIR+'progeny/'+fname+'_'+str(i)+'.csv',
                                        sep=',',header=0,index_col=0)
                    except:
                        break
                    os.remove(DIR+'progeny/'+fname+'_'+str(i)+'.csv')
                    progeny=pd.concat([progeny,temp])
                    temp=pd.read_csv(DIR+'tcga/'+fname+'_'+str(i)+'.csv',
                                        sep=',',header=0,index_col=0)
                    os.remove(DIR+'tcga/'+fname+'_'+str(i)+'.csv')
                    tcga=pd.concat([tcga,temp])
                progeny.to_csv('../results/benchmark/scores/progeny/'\
                                +fname+'.csv',sep=',')
                tcga.to_csv('../results/benchmark/scores/tcga/'\
                                +fname+'.csv',sep=',')