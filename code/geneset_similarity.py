import pandas as pd
import numpy as np
import pickle as pckl
import os

def get_composational_similarity(dname1,setname1,dname2,setname2,
                                sim_type='jaccard'):
    """get composational similarity between gene sets
    sim_type: {jaccard,overlap}"""
    fin=open('../results/genesets/%s/dicts/%s.pkl' % (dname1,setname1),
            'rb')
    data1=pckl.load(fin)
    fin.close()
    
    fin=open('../results/genesets/%s/dicts/%s.pkl' % (dname2,setname2),
            'rb')
    data2=pckl.load(fin)
    fin.close()
    
    results=pd.DataFrame(0,index=data1.keys(),columns=data2.keys())
    
    for gs1 in data1.keys():
        for gs2 in data2.keys():
            if sim_type=='jaccard':
                num=len(data1[gs1]&data2[gs2])
                denom=len(data1[gs1]|data2[gs2])
                results.loc[gs1,gs2]=num/denom
            elif sim_type=='overlap':
                num=len(data1[gs1]&data2[gs2])
                denom=np.min([len(data1[gs1]),len(data2[gs2])])
                results.loc[gs1,gs2]=num/denom
            elif sim_type=='first':
                num=len(data1[gs1]&data2[gs2])
                denom=len(data1[gs1])
                results.loc[gs1,gs2]=num/denom
            elif sim_type=='second':
                num=len(data1[gs1]&data2[gs2])
                denom=len(data1[gs2])
                results.loc[gs1,gs2]=num/denom
    results.to_csv('../results/similarity/%s_%s_%s.csv' \
                                        % (setname1,setname2,sim_type),
                    sep=',')

snames=['dorothea_AB','KEGG','REACTOME','CGP','BIOCARTA']    
            
for set1 in snames:
    for set2 in snames:
        if set1<set2:
            for sim_type in ['jaccard','overlap','first','second']:
                try:
                    get_composational_similarity('single',set1,'single',set2,
                                sim_type=sim_type)
                except:
                    print('Problem',set1,set2)