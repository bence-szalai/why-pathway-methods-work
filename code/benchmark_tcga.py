import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score as ROCAUC

tcga=pd.read_csv('../results/benchmark/tcga/raw/tcga_meta.csv',
                        sep=',',header=0,index_col=0)
                        
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/tcga/')]
for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/tcga/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    results=pd.DataFrame(index=list(set(tcga['TCGA'])),
                            columns=scores.columns)
    for tissue in results.index:
        indexes=tcga.index[tcga['TCGA']==tissue]
        y_tr=tcga.loc[indexes,'Tumor']
        for score in results.columns:
            y_pr=scores.loc[indexes,score]
            results.loc[tissue,score]=ROCAUC(y_tr,y_pr)
    results.to_csv('../results/benchmark/tcga/rocaucs/%s.csv' % method,
                    sep=',')

results=pd.DataFrame(index=list(set(tcga['TCGA'])),
                    columns=range(1000))
for i in range(1000):
    scores=pd.Series(np.random.uniform(size=len(tcga)),index=tcga.index)
    for tissue in results.index:
        indexes=tcga.index[tcga['TCGA']==tissue]
        y_tr=tcga.loc[indexes,'Tumor']
        y_pr=scores[indexes]
        results.loc[tissue,i]=ROCAUC(y_tr,y_pr)
results.to_csv('../results/benchmark/tcga/rocaucs/random_dist.csv',
                    sep=',')