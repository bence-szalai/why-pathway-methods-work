import pandas as pd
import numpy as np
import os

from sklearn.metrics import roc_auc_score as ROCAUC

def my_rocauc(x):
    return ROCAUC(y_tr,x)

gdsc=pd.read_csv('../results/benchmark/datasets/gdsc_mut.csv',
                sep=',',header=0,index_col=0)
gdsc.columns=gdsc.columns.astype(int)


methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/gdsc/overlap/')]
methods=[x for x in methods if 'random' in x]
try:
    methods.remove('.DS_S')
except:
    pass
                                                         
for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/gdsc/overlap/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    scores.index=pd.Series(scores.index).apply(lambda x:x[1:]).astype(int).values
    scores=scores.loc[gdsc.columns]
    results=pd.DataFrame(index=gdsc.index,
                            columns=scores.columns)
    for pathway in results.index:
        print(pathway)
        y_tr=gdsc.loc[pathway]
        results.loc[pathway,:]=scores.apply(my_rocauc,0)
    results.to_csv('../results/benchmark/rocaucs/gdsc/overlap/%s.csv' % method,
                    sep=',')