import pandas as pd
import numpy as np
import pickle as pckl

def make_gene_dict_dorothea(confidence):
    """makes gene - geneset and geneset - gene dicts for dorothea regulons
    saves them as pickles.
    confidence: {'A','B','C','D','E','BEST'}"""
    data=pd.read_csv('../results/genesets/single/raw/%s_viperRegulon.csv' \
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
    results.to_csv('../results/genesets/single/csvs/dorothea_%s.csv' 
                                                            % confidence,
                                                            sep=',')
    
    set_gene=results.groupby('Set')['Gene'].apply(list).to_dict()
    for s in set_gene:
        set_gene[s]=set(set_gene[s])
    fout=open('../results/genesets/single/dicts/dorothea_%s.pkl' \
                                                            % confidence, 'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
def make_multi_level_dorothea():
    """makes multi level dorothea gene sets"""
    fin=open('../results/genesets/single/dicts/dorothea_A.pkl','rb')
    dorothea=pckl.load(fin)
    fin.close()
    fname='dorothea_A'
    for confidence in ['B','C','D','E']:
        fname+=confidence
        fin=open('../results/genesets/single/dicts/dorothea_%s.pkl' \
                                                            % confidence, 'rb')
        genesets=pckl.load(fin)
        fin.close()
        for tf in genesets.keys():
            try: dorothea[tf]=dorothea[tf] | genesets[tf]
            except KeyError: dorothea[tf]=genesets[tf]
        fout=open('../results/genesets/single/dicts/%s.pkl' \
                                                            % fname, 'wb')
        pckl.dump(dorothea,fout)
        fout.close()
        sets=[]
        genes=[]
        for s in dorothea:
            sets+=[s]*len(dorothea[s])
            genes+=dorothea[s]
        results=pd.DataFrame(index=range(len(sets)),columns=['Set','Gene'])
        results['Set']=sets
        results['Gene']=genes
        results.to_csv('../results/genesets/single/csvs/%s.csv' \
                                                            % fname,sep=',')
    
def make_gene_dict_msigdb(setname):
    """makes gene - geneset and geneset - gene dicts for msigdb
    saves them as pickles.
    setname: {'BIOCARTA','CGP','KEGG','REACTOME'}"""
    results=pd.read_csv('../results/genesets/single/raw/%s.csv' % setname,
                    sep=',',header=0,index_col=0)
    results.columns=['Set','Gene']
    results=results.drop_duplicates()
    results.to_csv('../results/genesets/single/csvs/%s.csv' % setname,sep=',')
    set_gene=results.groupby('Set')['Gene'].apply(list).to_dict()
    for s in set_gene:
        set_gene[s]=set(set_gene[s])
    fout=open('../results/genesets/single/dicts/%s.pkl' % setname,'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
def make_overlap_gene_sets(setname1,setname2):
    fin=open('../results/genesets/single/dicts/%s.pkl' % setname1, 'br')
    geneset1=pckl.load(fin)
    fin.close()
    fin=open('../results/genesets/single/dicts/%s.pkl' % setname2, 'br')
    geneset2=pckl.load(fin)
    fin.close()
    overlap={}
    set1_minus_set2={}
    set2_minus_set1={}
    for set1 in geneset1:
        for set2 in geneset2:
            overlap[set1+'*'+set2]=geneset1[set1] & geneset2[set2]
            set1_minus_set2[set1+'*'+set2]=geneset1[set1] - geneset2[set2]
            set2_minus_set1[set1+'*'+set2]=geneset2[set2] - geneset1[set1]
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
    fnames=['%s_overlap_%s.pkl' % (setname1,setname2),
            '%s_minus_%s.pkl' % (setname1,setname2),
            '%s_minus_%s.pkl' % (setname2,setname1)]
    genesets=[overlap,set1_minus_set2,set2_minus_set1]    
    for i in range(3):
        fout=open('../results/genesets/overlap/dicts/'+fnames[i],'wb')
        pckl.dump(genesets[i],fout)
        fout.close()
    fnames=['%s_overlap_%s' % (setname1,setname2),
            '%s_minus_%s' % (setname1,setname2),
            '%s_minus_%s' % (setname2,setname1)]
    for i in range(3):
        for j in range(len(genesets[i])//1000+1):
            sets=[]
            genes=[]
            for s in list(genesets[i].keys())[j*1000:(j+1)*1000]:
                sets+=[s]*len(genesets[i][s])
                genes+=genesets[i][s]
            results=pd.DataFrame(index=range(len(sets)),columns=['Set','Gene'])
            results['Set']=sets
            results['Gene']=genes
            results.to_csv('../results/genesets/overlap/csvs/%s_%i.csv' \
                                        % (fnames[i],j),sep=',')
                                                        
def make_random_geneset(fnames=['KEGG'],s=19890904,n=200):
    np.random.seed(s)
    genes=set([])
    lens=[]
    for fname in fnames:
        fin=open('../results/genesets/single/dicts/%s.pkl' % fname,'rb')
        geneset=pckl.load(fin)
        fin.close()
        for setname in geneset:
            genes=genes | geneset[setname]
            lens.append(len(geneset[setname]))
    genes=list(genes)
    random_set={}
    for i in range(n):
        l=np.random.choice(lens,1)[0]
        random_set['Geneset'+str(i)]=set(np.random.choice(genes,l,False))
    fname='_'.join(fnames+[str(n),str(s)])
    fout=open('../results/genesets/random/dicts/%s.pkl' % fname,'wb')
    pckl.dump(random_set,fout)
    fout.close()
    
    sets=[]
    genes=[]
    for s in random_set:
        sets+=[s]*len(random_set[s])
        genes+=random_set[s]
    results=pd.DataFrame(index=range(len(sets)),columns=['Set','Gene'])
    results['Set']=sets
    results['Gene']=genes
    results.to_csv('../results/genesets/random/csvs/%s.csv' \
                                                        % fname,sep=',')
                                    
for setname1 in ['BIOCARTA','CGP','KEGG','REACTOME','dorothea_AB']:
    for setname2 in ['BIOCARTA','CGP','KEGG','REACTOME','dorothea_AB']:
        if setname1<setname2:    
            try:
                make_overlap_gene_sets(setname1,setname2)
            except:
                print('Problem',setname1,setname2)
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    