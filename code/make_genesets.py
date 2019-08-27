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
    results['Gene']=data.apply(lambda x:x.split('tfmode.')[-1])
    results=results.drop_duplicates()
    #set - gene dict
    set_gene=results.groupby('Set')['Gene'].apply(list).to_dict()
    fout=open('../results/genesets/dorothea/dicts/dorothea_%s_set_gene.pkl' \
                                                            % confidence, 'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
    #gene - set dict
    gene_set=results.groupby('Gene')['Set'].apply(list).to_dict()
    fout=open('../results/genesets/dorothea/dicts/dorothea_%s_gene_set.pkl' \
                                                            % confidence, 'wb')
    pckl.dump(gene_set,fout)
    fout.close()
    
def make_multi_level_dorothea():
    """makes multi level dorothea gene sets"""
    fin=open('../results/genesets/dorothea/dicts/dorothea_A_set_gene.pkl','rb')
    dorothea=pckl.load(fin)
    fin.close()
    fname='dorothea_A'
    for confidence in ['B','C','D','E']:
        fname+=confidence
        fin=open('../results/genesets/dorothea/dicts/dorothea_%s_set_gene.pkl' \
                                                            % confidence, 'rb')
        genesets=pckl.load(fin)
        fin.close()
        for tf in genesets.keys():
            try: dorothea[tf]+=genesets[tf]
            except KeyError: dorothea[tf]=genesets[tf]
        fout=open('../results/genesets/dorothea/dicts/%s_set_gene.pkl' \
                                                            % fname, 'wb')
        pckl.dump(dorothea,fout)
        fout.close()
        
    fin=open('../results/genesets/dorothea/dicts/dorothea_A_gene_set.pkl','rb')
    dorothea=pckl.load(fin)
    fin.close()
    fname='dorothea_A'
    for confidence in ['B','C','D','E']:
        fname+=confidence
        fin=open('../results/genesets/dorothea/dicts/dorothea_%s_gene_set.pkl' \
                                                            % confidence, 'rb')
        genesets=pckl.load(fin)
        fin.close()
        for gene in genesets.keys():
            try: dorothea[gene]+=genesets[gene]
            except KeyError: dorothea[gene]=genesets[gene]
        fout=open('../results/genesets/dorothea/dicts/%s_gene_set.pkl' \
                                                            % fname, 'wb')
        pckl.dump(dorothea,fout)
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
    fout=open('../results/genesets/msigdb/dicts/%s_set_gene.pkl' % setname,'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
    #gene - set dict
    gene_set=results.groupby('Gene')['Set'].apply(list).to_dict()
    fout=open('../results/genesets/msigdb/dicts/dicts%s_gene_set.pkl' % setname,
                'wb')
    pckl.dump(gene_set,fout)
    fout.close()
    
def make_csv_geneset(dname,fname):
    """ makese csv from dict, for viper"""
    fin=open('../results/genesets/%s/dicts/%s.pkl' % (dname,fname),'br')
    geneset=pckl.load(fin)
    fin.close()
    l=max([len(geneset[x]) for x in geneset.keys()])
    results=pd.DataFrame('NONE',index=geneset.keys(),columns=range(l))
    for setname in geneset:
        results.loc[setname,range(len(geneset[setname]))]=geneset[setname]
    results.to_csv('../results/genesets/%s/csvs/%s.csv' % (dname,fname),sep=',')
    

def make_overlap_gene_sets(dname1,setname1,dname2,setname2):
    fin=open('../results/genesets/%s/dicts/%s_set_gene.pkl' % (dname1,setname1),
                'br')
    geneset1=pckl.load(fin)
    fin.close()
    fin=open('../results/genesets/%s/dicts/%s_set_gene.pkl' % (dname2,setname2),
                'br')
    geneset2=pckl.load(fin)
    fin.close()
    overlap={}
    set1_minus_set2={}
    set2_minus_set1={}
    for set1 in geneset1:
        for set2 in geneset2:
            overlap[set1+'*'+set2]=list(set(geneset1[set1]) & \
                                            set(geneset2[set2]))
            set1_minus_set2[set1+'*'+set2]=list(set(geneset1[set1]) - \
                                            set(geneset2[set2]))
            set2_minus_set1[set1+'*'+set2]=list(set(geneset2[set2]) - \
                                            set(geneset1[set1]))
    empty=[]
    for geneset in overlap:
        if len(overlap[geneset])<4:
            empty.append(geneset)
        elif len(set1_minus_set2[geneset])<4:
            empty.append(geneset)
        elif len(set2_minus_set1[geneset])<4:
            empty.append(geneset)
    for geneset in empty:
        del overlap[geneset]
        del set1_minus_set2[geneset]
        del set2_minus_set1[geneset]
    fout=open('../results/genesets/overlaps/dicts/%s_overlap_%s.pkl' % (setname1,
                                                                setname2),
                'wb')
    pckl.dump(overlap,fout)
    fout=open('../results/genesets/overlaps/dicts/%s_minus_%s.pkl' % (setname1,
                                                                setname2),
                'wb')
    pckl.dump(set1_minus_set2,fout)
    fout=open('../results/genesets/overlaps/dicts/%s_minus_%s.pkl' % (setname2,
                                                                setname1),
                'wb')
    pckl.dump(set2_minus_set1,fout)
    
       
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    