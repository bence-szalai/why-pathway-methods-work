import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score as ROCAUC

progeny=pd.read_csv('../results/benchmark/datasets/progeny_meta.csv',
                        sep=',',header=0,index_col=0)
                        
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/progeny/')]
try:
    methods.remove('.DS_S')
except:
    pass

for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/progeny/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    results=pd.DataFrame(index=list(set(progeny['pathway'])),
                            columns=scores.columns)
    for pathway in results.index:
        y_tr=(progeny['pathway']==pathway)*1
        for score in results.columns:
            y_pr=scores[score]
            results.loc[pathway,score]=ROCAUC(y_tr,y_pr)
    results.to_csv('../results/benchmark/progeny/rocaucs/%s.csv' % method,
                    sep=',')

results=pd.DataFrame(index=list(set(progeny['pathway'])),
                    columns=range(1000))
for i in range(1000):
    y_pr=np.random.uniform(size=len(progeny))
    for pathway in results.index:
        y_tr=(progeny['pathway']==pathway)*1
        results.loc[pathway,i]=ROCAUC(y_tr,y_pr)
results.to_csv('../results/benchmark/progeny/rocaucs/random_dist.csv',
                    sep=',')