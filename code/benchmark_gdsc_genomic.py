import pandas as pd
import numpy as np
import os

from sklearn.metrics import roc_auc_score as ROCAUC

def my_rocauc(x):
    return ROCAUC(y_tr,x)

gdsc=pd.read_csv('../data/gdsc/PANCAN_simple_MOBEM.rdata.tsv',
                sep='\t',header=0,index_col=0)
gdsc.columns=gdsc.columns.astype(int)
scores=pd.read_csv('../results/benchmark/scores/gdsc/single/BEST_dorothea_AB_filtered.csv',
                sep=',',header=0,index_col=0).T
scores.index=pd.Series(scores.index).apply(lambda x:x[1:]).astype(int).values
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/gdsc/single/')]
cosmics=list(set(gdsc.columns)&set(scores.index))
gdsc=gdsc[cosmics]
features=gdsc.index[gdsc.sum(1)>10]
gdsc=gdsc.loc[features]
fil=pd.Series(gdsc.index).apply(lambda x: '_mut' in x).values
gdsc=gdsc[fil]

try:
    methods.remove('.DS_S')
except:
    pass
    
for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/gdsc/single/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    scores.index=pd.Series(scores.index).apply(lambda x:x[1:]).astype(int).values
    scores=scores.loc[gdsc.columns]
    results=pd.DataFrame(index=gdsc.index,
                            columns=scores.columns)
    for pathway in results.index:
        y_tr=gdsc.loc[pathway]
        results.loc[pathway,:]=scores.apply(my_rocauc,0)
    results.to_csv('../results/benchmark/rocaucs/gdsc/single/%s.csv' % method,
                    sep=',')
    
results=pd.DataFrame(index=gdsc.index,columns=range(1000))
for i in range(1000):
    print(i)
    y_pr=np.random.uniform(size=gdsc.shape[1])
    for pathway in results.index:
        y_tr=gdsc.loc[pathway]
        results.loc[pathway,i]=ROCAUC(y_tr,y_pr)
results.to_csv('../results/benchmark/rocaucs/gdsc/random_dist.csv',
                    sep=',')
