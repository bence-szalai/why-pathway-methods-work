import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score as ROCAUC

def my_rocauc(x):
    return ROCAUC(y_tr,x)

progeny=pd.read_csv('../results/benchmark/datasets/progeny_meta.csv',
                        sep=',',header=0,index_col=0)
                        
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/progeny/overlap/')]
methods=[x for x in methods if 'random' in x]
try:
    methods.remove('.DS_S')
except:
    pass

for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/progeny/overlap/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    results=pd.DataFrame(index=list(set(progeny['pathway'])),
                            columns=scores.columns)
    for pathway in results.index:
        y_tr=(progeny['pathway']==pathway)*1
        results.loc[pathway,:]=scores.apply(my_rocauc,0)
    results.to_csv('../results/benchmark/rocaucs/progeny/overlap/%s.csv' % method,
                    sep=',')