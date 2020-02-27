import pandas as pd
import numpy as np
import os

gdsc=pd.read_csv('../results/benchmark/datasets/gdsc_meta.csv',
                        sep=',',header=0,index_col=0)
gdsc.columns=gdsc.columns.astype(int)
cell_anno=pd.read_excel('../data/gdsc/Cell_Lines_Details.xlsx',skipfooter=1)
cosmics=cell_anno['COSMIC identifier'][cell_anno['Cancer Type\n(matching TCGA label)']=='BRCA'].values    
cosmics=list(set(cosmics)&set(gdsc.columns))
gdsc=gdsc[cosmics]                 
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/gdsc/single/')]

try:
    methods.remove('.DS_S')
except:
    pass
for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/gdsc/single/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    scores.index=pd.Series(scores.index).apply(lambda x:x[1:]).astype(int).values
    cosmics=list(set(gdsc.columns) & set(scores.index))
    
   
  
    results=pd.concat([gdsc[cosmics].T,scores.loc[cosmics]],1)
    results=results.corr().loc[gdsc.index,scores.columns]
        
    results.abs().to_csv('../results/benchmark/z_scores/gdsc/single/'+method+'.csv',
                    sep=',')
        
        
        
    
    

    
    



