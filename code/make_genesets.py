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
    
def make_different_level_BEST_dorothea():
    data=pd.read_csv('../results/genesets/single/raw/BEST_viperRegulon.csv',
                    sep=',',header=0,index_col=0)
    data=pd.Series(data.index)
    fil=~(data.apply(lambda x: 'likelihood' in x))
    data=data[fil]
    data.index=range(len(data))
        
    results=pd.DataFrame(index=data.index,columns=['Set','Gene','Confidence'])
    results['Set']=data.apply(lambda x:x.split('_')[0])
    results['Confidence']=data.apply(lambda x:x.split('_')[1])\
                    .apply(lambda x:x.split(' - ')[0])
    results['Gene']=data.apply(lambda x:x.split('tfmode.')[-1])
    results=results.drop_duplicates()
    for c in ['A','B','C','D','E']:
        fil=results['Confidence']==c
        results_conf=results[fil][['Set','Gene']].copy()
        results_conf.to_csv('../results/genesets/single/csvs/BEST_dorothea_%s.csv' 
                                                            % c, sep=',')
        set_gene=results_conf.groupby('Set')['Gene'].apply(list).to_dict()
        for s in set_gene:
            set_gene[s]=set(set_gene[s])
        fout=open('../results/genesets/single/dicts/BEST_dorothea_%s.pkl' \
                                                            % c, 'wb')
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
def make_multi_level_dorothea_BEST():
    data_dict={}
    for c in ['A','B','C','D']:
        data_dict[c]=pd.read_csv('../results/genesets/single/csvs/BEST_dorothea_%s.csv'
                    % c,sep=',',header=0,index_col=0)
    data_AB=pd.concat([data_dict['A'],data_dict['B']])
    data_CD=pd.concat([data_dict['C'],data_dict['D']])
    data_AB.to_csv('../results/genesets/single/csvs/BEST_dorothea_AB.csv',sep=',')
    data_CD.to_csv('../results/genesets/single/csvs/BEST_dorothea_CD.csv',sep=',')
    set_gene=data_AB.groupby('Set')['Gene'].apply(list).to_dict()
    for s in set_gene:
        set_gene[s]=set(set_gene[s])
    fout=open('../results/genesets/single/dicts/BEST_dorothea_AB.pkl', 'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    set_gene=data_CD.groupby('Set')['Gene'].apply(list).to_dict()
    for s in set_gene:
        set_gene[s]=set(set_gene[s])
    fout=open('../results/genesets/single/dicts/BEST_dorothea_CD.pkl', 'wb')
    pckl.dump(set_gene,fout)
    fout.close()
    
    
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
    
def make_pickle_from_genese(setname):
    data=pd.read_csv('../results/genesets/single/csvs/%s.csv' % setname,
                    sep=',',header=0,index_col=0)
    data_dict={}
    for geneset in data['Set'].unique():
        fil=data['Set']==geneset
        genes=set(data[fil]['Gene'].unique())
        data_dict[geneset]=genes
    fout=open('../results/genesets/single/dicts/%s.pkl' % setname, 'wb')
    pckl.dump(data_dict,fout)
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
                                    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    