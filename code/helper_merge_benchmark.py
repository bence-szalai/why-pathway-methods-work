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
        progeny=pd.read_csv(DIR+'progeny/overlap/'+fname+'_0_abs.csv',
                                        sep=',',header=0,index_col=0)
        tcga=pd.read_csv(DIR+'tcga/overlap/'+fname+'_0_abs.csv',
                                        sep=',',header=0,index_col=0)
        os.remove(DIR+'progeny/overlap/'+fname+'_0_abs.csv')
        os.remove(DIR+'tcga/overlap/'+fname+'_0_abs.csv')
        i=0
        while True:
            i+=1
            try:
                temp=pd.read_csv(DIR+'progeny/overlap/'+\
                                    fname+'_'+str(i)+'_abs.csv',
                                    sep=',',header=0,index_col=0)
            except:
                break
            os.remove(DIR+'progeny/overlap/'+fname+'_'+str(i)+'_abs.csv')
            progeny=pd.concat([progeny,temp])
            temp=pd.read_csv(DIR+'tcga/overlap/'+fname+'_'+str(i)+'_abs.csv',
                                sep=',',header=0,index_col=0)
            os.remove(DIR+'tcga/overlap/'+fname+'_'+str(i)+'_abs.csv')
            tcga=pd.concat([tcga,temp])
        progeny.to_csv(DIR+'progeny/overlap/'+fname+'_abs.csv',sep=',')
        tcga.to_csv(DIR+'tcga/overlap/'+fname+'_abs.csv',sep=',')
        