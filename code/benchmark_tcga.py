import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score as ROCAUC

def my_rocauc(x):
    return ROCAUC(y_tr,x)

tcga=pd.read_csv('../results/benchmark/datasets/tcga_meta.csv',
                        sep=',',header=0,index_col=0)
new_index=[]
for sample in tcga.index:
    new_index.append(sample.replace('-','.'))
tcga.index=new_index
                        
methods=[x[:-4] for x in os.listdir('../results/benchmark/scores/tcga/overlap/')]

try:
    methods.remove('.DS_S')
except:
    pass
    
for method in methods:
    print(method)
    scores=pd.read_csv('../results/benchmark/scores/tcga/overlap/%s.csv' % method,
                    sep=',',header=0,index_col=0).T
    results=pd.DataFrame(index=list(set(tcga['TCGA'])),
                            columns=scores.columns)
    for tissue in results.index:
        indexes=tcga.index[tcga['TCGA']==tissue]
        y_tr=tcga.loc[indexes,'Tumor']
        results.loc[tissue,:]=scores.loc[indexes].apply(my_rocauc,0)
    results.to_csv('../results/benchmark/rocaucs/tcga/overlap/%s.csv' % method,
                    sep=',')

# results=pd.DataFrame(index=list(set(tcga['TCGA'])),
#                     columns=range(1000))
# for i in range(1000):
#     scores=pd.Series(np.random.uniform(size=len(tcga)),index=tcga.index)
#     for tissue in results.index:
#         indexes=tcga.index[tcga['TCGA']==tissue]
#         y_tr=tcga.loc[indexes,'Tumor']
#         y_pr=scores[indexes]
#         results.loc[tissue,i]=ROCAUC(y_tr,y_pr)
# results.to_csv('../results/benchmark/rocaucs/tcga/random_dist.csv',
#                     sep=',')