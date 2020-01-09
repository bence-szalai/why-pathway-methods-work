import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numpy.linalg import inv

def remove_self_edges(x):
    y=x.copy()
    y[y.name]=0
    return y


AM=pd.read_csv('../data/omnipath/AM_GC.csv',sep=',',header=0,index_col=0)
AM=AM.apply(remove_self_edges,axis=0)

prot_ids=AM.index
AM=AM.values

n = AM.shape[0]
degree = AM.sum(axis=1)
p = AM / degree
for k in [4,5,6,7]:
    c = np.eye(n)
    for i in range(k):
        c = np.dot(c,p) + np.eye(n)
    DM=squareform(pdist(c,metric='cityblock'))
    DM=pd.DataFrame(DM,index=prot_ids,columns=prot_ids)
    DM.to_csv('../data/omnipath/DSD_%i.csv' % k,sep=',')
pi = degree / degree.sum()
DM=squareform(pdist(inv(np.eye(n) - p - pi.T),metric='cityblock'))
DM=pd.DataFrame(DM,index=prot_ids,columns=prot_ids)
DM.to_csv('../data/omnipath/DSD_conv.csv',sep=',')