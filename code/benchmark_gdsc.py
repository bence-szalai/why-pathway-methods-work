import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy.stats import pearsonr as pcor

response=pd.read_excel('../data/gdsc/v17.3_fitted_dose_response.xlsx')
response=response[['DRUG_ID','COSMIC_ID','LN_IC50']]
response['DRUG_ID']=response['DRUG_ID'].astype(str)
response['COSMIC_ID']=response['COSMIC_ID'].astype(str)

anno=pd.read_csv('../results/benchmark/gdsc/raw/gdsc_cell.csv',sep=',',
                header=0,index_col=0)
anno['COSMIC']=anno['COSMIC'].astype(str)
anno.index=anno.index.astype(str)
                
scores=pd.read_csv('../results/benchmark/scores/gdsc/KEGG.csv',sep=',',
                header=0,index_col=0)

cosmics=list(set(anno['COSMIC'])&set(response['COSMIC_ID'])&set(scores.columns))                
fil=np.in1d(anno['COSMIC'],cosmics)
anno=anno[fil]
fil=np.in1d(response['COSMIC_ID'],cosmics)
response=response[fil]
scores=scores[cosmics]
scores=scores.T

results=pd.DataFrame(index=list(set(response['DRUG_ID'])),
                    columns=scores.columns)
for drug in results.index:
    print(drug)
    fil=response['DRUG_ID']==drug
    response_drug=response[fil]
    response_drug['TCGA']=anno.loc[response_drug['COSMIC_ID'].values,
                                                        'TCGA'].values
    response_drug['MSI']=anno.loc[response_drug['COSMIC_ID'].values,
                                                        'MSI'].values
    model_1=smf.ols('LN_IC50 ~ TCGA + MSI',data=response_drug).fit()
    resid_1=model_1.resid
    for score in results.columns:
        response_drug['Score']=scores.loc[response_drug['COSMIC_ID'].values,
                                                        score].values
        model_2=smf.ols('Score ~ TCGA + MSI',data=response_drug).fit()
        resid_2=model_2.resid
        results.loc[drug,score]=pcor(resid_1,resid_2)[0]
        
        
                




