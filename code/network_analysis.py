import pandas as pd
import numpy as np

data=pd.read_csv('../data/interactions.txt',sep='\t',header=0,index_col=None)

anno=pd.read_csv('../data/uniprot_symbol.csv',sep=',',header=0,index_col=0)

fil=~pd.isnull(anno['hgnc_symbol'])
anno=anno[fil]
anno=anno.drop_duplicates('uniprotswissprot')
anno=anno.drop_duplicates('hgnc_symbol')
anno.index=anno['uniprotswissrot']
anno=anno