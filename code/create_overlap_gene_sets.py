import pandas as pd
import numpy as np

def create_overlap(setname1,setname2):
    geneset1=pd.read_csv('../results/genesets/single/csvs/%s.csv' % setname1,
                    sep=',',header=0,index_col=0)
    geneset1['Indicator']=1
    geneset2=pd.read_csv('../results/genesets/single/csvs/%s.csv' % setname2,
                        sep=',',header=0,index_col=0)
    geneset2['Indicator']=1
    geneset1=geneset1.pivot(index='Gene',columns='Set',values='Indicator')
    geneset2=geneset2.pivot(index='Gene',columns='Set',values='Indicator')
    geneset1[geneset1.isnull()]=0; geneset1=geneset1.astype(int)
    geneset2[geneset2.isnull()]=0; geneset2=geneset2.astype(int)
    all_genes=list(set(geneset1.index)&set(geneset2.index))
    geneset1=geneset1.loc[all_genes]; geneset2=geneset2.loc[all_genes]
    
    tensor_o=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_o+=geneset1.values.reshape((len(all_genes),geneset1.shape[1],1))
    tensor_o*=geneset2.values.reshape((len(all_genes),1,geneset2.shape[1]))
    tensor_o=tensor_o.reshape((1,-1))[0]
    
    tensor_1m2=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_1m2+=geneset1.values.reshape((len(all_genes),geneset1.shape[1],1))
    tensor_1m2-=geneset2.values.reshape((len(all_genes),1,geneset2.shape[1]))
    tensor_1m2=tensor_1m2.reshape((1,-1))[0]
    
    tensor_2m1=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_2m1-=geneset1.values.reshape((len(all_genes),geneset1.shape[1],1))
    tensor_2m1+=geneset2.values.reshape((len(all_genes),1,geneset2.shape[1]))
    tensor_2m1=tensor_2m1.reshape((1,-1))[0]
    
    tensor_set1=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_set1+=np.array(range(geneset1.shape[1])).reshape((1,-1,1))
    tensor_set1=tensor_set1.reshape((1,-1))[0]
    
    tensor_set2=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_set2+=np.array(range(geneset2.shape[1])).reshape((1,1,-1))
    tensor_set2=tensor_set2.reshape((1,-1))[0]
    
    tensor_gene=np.zeros((len(all_genes),geneset1.shape[1],geneset2.shape[1]))
    tensor_gene+=np.array(range(geneset2.shape[0])).reshape((-1,1,1))
    tensor_gene=tensor_gene.reshape((1,-1))[0]
    
    results_o=pd.DataFrame(index=range(len(tensor_o)))
    results_o['Gene']=tensor_gene
    results_o['Set1']=tensor_set1
    results_o['Set2']=tensor_set2
    results_o['Indicator']=tensor_o
    results_o=results_o.astype(int)
    fil=results_o['Indicator']==1
    results_o=results_o[fil]
    results_o['Gene']=geneset1.index[results_o['Gene'].values]
    results_o['Set1']=geneset1.columns[results_o['Set1'].values]
    results_o['Set2']=geneset2.columns[results_o['Set2'].values]
    results_o['Set']=results_o['Set1'] + '*' + results_o['Set2']
    results_o=results_o[['Set','Gene']]
    
    
    
    
    
    
    
               