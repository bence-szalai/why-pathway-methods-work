import pandas as pd
import numpy as np

def get_composational_similarity(setname1,setname2):
    """get composational similarity between gene sets
    sim_type: {jaccard,overlap}"""

    data1=pd.read_csv('../results/genesets/single/csvs/%s.csv' % setname1,
                    sep=',',header=0,index_col=0)
    data1['Indicator']=1
    data1=data1.pivot(index='Gene',columns='Set',values='Indicator')


    data2=pd.read_csv('../results/genesets/single/csvs/%s.csv' % setname2,
                    sep=',',header=0,index_col=0)
    data2['Indicator']=1
    data2=data2.pivot(index='Gene',columns='Set',values='Indicator')

    all_genes=list(set(data1.index) | set(data2.index))

    not_genes1 = list(set(all_genes) - set(data1.index))
    not_genes2 = list(set(all_genes) - set(data2.index))

    data1 = pd.concat([data1,
                    pd.DataFrame(0, index=not_genes1, columns=data1.columns)])
    data2 = pd.concat([data2,
                    pd.DataFrame(0, index=not_genes2, columns=data2.columns)])

    data1 = data1.loc[all_genes]
    data2 = data2.loc[all_genes]

    #fil=pd.isnull(data1); data1[fil]=0;
    data1=data1.astype(int)
    #fil=pd.isnull(data2); data2[fil]=0;
    data2=data2.astype(int)

    intersection=pd.DataFrame(np.dot(data1.T,data2),
                            index=data1.columns,columns=data2.columns)
    union=data1.sum().values.reshape((-1,1))+\
            data2.sum().values.reshape((1,-1))-\
            intersection.values
    union=pd.DataFrame(union,index=data1.columns,columns=data2.columns)

    jaccard=intersection/union
    overlap1=(intersection.T/data1.sum()).T
    overlap2=intersection/data2.sum()
    overlap=overlap1.copy()
    fil=overlap2>overlap1
    overlap[fil]=overlap2[fil]

    jaccard.to_csv('../results/similarity/%s_%s_jaccard.csv' %\
                    (setname1,setname2),sep=',')
    overlap.to_csv('../results/similarity/%s_%s_overlap.csv' %\
                    (setname1,setname2),sep=',')
    overlap1.to_csv('../results/similarity/%s_%s_first.csv' %\
                    (setname1,setname2),sep=',')
    overlap2.to_csv('../results/similarity/%s_%s_second.csv' %\
                    (setname1,setname2),sep=',')
