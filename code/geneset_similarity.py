import pandas as pd
import numpy as np
import pickle as pckl

def get_composational_similarity(dorothea,msigdb,sim_type='jaccard'):
    """get composational similarity between gene sets
    sim_type: {jaccard,overlap}"""
    fin=open('../results/genesets/dorothea/%s_set_gene.pkl' % dorothea,'rb')
    data1=pckl.load(fin)
    fin.close()
    
    fin=open('../results/genesets/msigdb/%s_set_gene.pkl' % msigdb,'rb')
    data2=pckl.load(fin)
    fin.close()
    
    results=pd.DataFrame(0,index=data1.keys(),columns=data2.keys())
    
    for tf in data1.keys():
        for gs in data2.keys():
            if sim_type=='jaccard':
                num=len(set(data1[tf])&set(data2[gs]))
                denom=len(set(data1[tf])|set(data2[gs]))
                results.loc[tf,gs]=num/denom
            elif sim_type=='overlap':
                num=len(set(data1[tf])&set(data2[gs]))
                denom=np.min([len(data1[tf]),len(data2[gs])])
                results.loc[tf,gs]=num/denom
    results.to_csv('../results/similarity/%s_%s_%s.csv' \
                                        % (dorothea,msigdb,sim_type),
                    sep=',')
    return results
    
for dorothea in ['dorothea_A','dorothea_B','dorothea_C','dorothea_D',
                'dorothea_E','dorothea_BEST']:
    for msigdb in ['BIOCARTA','CGP','KEGG','REACTOME']:
        for sim_type in ['jaccard','overlap']:
            get_composational_similarity(dorothea,msigdb,sim_type)
        