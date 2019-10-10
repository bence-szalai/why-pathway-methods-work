import pandas as pd
import numpy as np

th=4.5

bname='progeny'
dname1='REACTOME'
dname2='dorothea_AB'

data=pd.read_csv('../results/benchmark/datasets/%s_data.csv' % bname,
                    sep=',',header=0,index_col=0)

activity1=pd.read_csv('../results/benchmark/z_scores/%s/%s.csv'\
                    % (bname,dname1),sep=',',header=0,index_col=0)

activity2=pd.read_csv('../results/benchmark/z_scores/%s/%s.csv'\
                    % (bname,dname2),sep=',',header=0,index_col=0)
                    
activity1=(activity1.abs()>th)*1
activity2=(activity2.abs()>th)*1

overlap_activity=pd.DataFrame(np.dot(activity1.T,activity2),
                    index=activity1.columns,columns=activity2.columns)
overlap_activity[overlap_activity>0]=1

geneset1=pd.read_csv('../results/genesets/single/csvs/%s.csv' % dname1,
                    sep=',',header=0,index_col=0)
geneset1['Indicator']=1
geneset2=pd.read_csv('../results/genesets/single/csvs/%s.csv' % dname2,
                    sep=',',header=0,index_col=0)
geneset2['Indicator']=1

geneset1=geneset1.pivot(index='Gene',columns='Set',values='Indicator')
geneset2=geneset2.pivot(index='Gene',columns='Set',values='Indicator')

all_genes=list(set(geneset1.index)&set(geneset2.index)\
                    &set(data.index)) #we need only common

geneset1=geneset1.loc[all_genes]
geneset1[geneset1.isnull()]=0; geneset1=geneset1.astype(int)
geneset2=geneset2.loc[all_genes]
geneset2[geneset2.isnull()]=0; geneset2=geneset2.astype(int)

overlap_comp=pd.DataFrame(np.dot(geneset1.T,geneset2),
                    index=geneset1.columns,columns=geneset2.columns)
fil=overlap_comp>=4
overlap_comp=pd.DataFrame(0,index=overlap_comp.index,
                    columns=overlap_comp.columns)
overlap_comp[fil]=1

indexes1=list(set(overlap_activity.index)&set(overlap_comp.index))
indexes2=list(set(overlap_activity.columns)&set(overlap_comp.columns))

overlap=overlap_activity.loc[indexes1,indexes2] * overlap_comp.loc[indexes1,
                                                                indexes2]







