import pandas as pd
import numpy as np

import os

DIR='../results/benchmark/scores/'

for geneset in ['KEGG','BIOCARTA','REACTOME']:
    for set_type in ['minus','minus_rev','overlap']:
        if set_type=='minus_rev':
            fname='_'.join(['dorothea_AB','minus',geneset])
        else:
            fname='_'.join([geneset,set_type,'dorothea_AB'])
        progeny=pd.read_csv(DIR+'progeny/overlap/'+fname+'_0.csv',
                                        sep=',',header=0,index_col=0)
        tcga=pd.read_csv(DIR+'tcga/overlap/'+fname+'_0.csv',
                                        sep=',',header=0,index_col=0)
        os.remove(DIR+'progeny/overlap/'+fname+'_0.csv')
        os.remove(DIR+'tcga/overlap/'+fname+'_0.csv')
        i=0
        while True:
            i+=1
            try:
                temp=pd.read_csv(DIR+'progeny/overlap/'+\
                                    fname+'_'+str(i)+'.csv',
                                    sep=',',header=0,index_col=0)
            except:
                break
            os.remove(DIR+'progeny/overlap/'+fname+'_'+str(i)+'.csv')
            progeny=pd.concat([progeny,temp])
            temp=pd.read_csv(DIR+'tcga/overlap/'+fname+'_'+str(i)+'.csv',
                                sep=',',header=0,index_col=0)
            os.remove(DIR+'tcga/overlap/'+fname+'_'+str(i)+'.csv')
            tcga=pd.concat([tcga,temp])
        progeny.to_csv(DIR+'progeny/overlap/'+fname+'.csv',sep=',')
        tcga.to_csv(DIR+'tcga/overlap/'+fname+'.csv',sep=',')
        