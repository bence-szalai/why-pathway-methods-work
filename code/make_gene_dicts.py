import pandas as pd
import pickle as pckl

def make_gene_dict_dorothea(confidence):
    """makes gene - geneset and geneset - gene dicts for dorothea regulons
    saves them as pickles.
    confidence: {'A','B','C','D','E','BEST'}"""
    data=pd.read_csv('../results/genesets/dorothea/raw/%s_viperRegulon.csv' \
                                                                % confidence,
                    sep=',',header=0,index_col=0)
    data=pd.Series(data.index)
    fil=~(data.apply(lambda x: 'likelihood' in x))
    data=data[fil]
    data.index=range(len(data))
    
    results=pd.DataFrame(index=data.index,columns=['Set','Gene'])
    results['Set']=data.apply(lambda x:x.split(' - ')[0])
    results['Gene']=data.apply(lambda x:x.split('.')[-1])
    results=results.drop_duplicates()
    #set - gene dict
    set_gene=results.groupby('Set')['Gene'].apply(list).to_dict()
    fout=open('../results/genesets/dorothea/dorothea_%s_set_gene.pkl' \
                                                            % confidence
            ,'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
    #gene - set dict
    gene_set=results.groupby('Gene')['Set'].apply(list).to_dict()
    fout=open('../results/genesets/dorothea/dorothea_%s_gene_set.pkl' \
                                                            % confidence
            ,'wb')
    pckl.dump(gene_set,fout)
    fout.close()
    
def make_gene_dict_msigdb(setname):
    """makes gene - geneset and geneset - gene dicts for msigdb
    saves them as pickles.
    setname: {'BIOCARTA','CGP','KEGG','REACTOME'}"""
    results=pd.read_csv('../results/genesets/msigdb/raw/%s.csv' % setname,
                    sep=',',header=0,index_col=0)
    results=results.drop_duplicates()
    results.columns=['Set','Gene']
    #set - gene dict
    set_gene=results.groupby('Set')['Gene'].apply(list).to_dict()
    fout=open('../results/genesets/msigdb/%s_set_gene.pkl' % setname,'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
    #gene - set dict
    gene_set=results.groupby('Gene')['Set'].apply(list).to_dict()
    fout=open('../results/genesets/msigdb/%s_gene_set.pkl' % setname,'wb')
    pckl.dump(gene_set,fout)
    fout.close()