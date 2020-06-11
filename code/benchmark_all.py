import pandas as pd
import numpy as np
import os

from sklearn.metrics import roc_auc_score


def my_rocauc(x, y_tr):
    return roc_auc_score(y_tr, x)


np.random.seed(19890904)
for dname in ['progeny', 'gdsc']:
    meta = pd.read_csv('../results/benchmark/datasets/%s_meta.csv' % dname,
                       sep=',', header=0, index_col=0)
    if dname == 'progeny':
        meta.index = pd.Series(meta.index).apply(lambda x:x.replace('-', '.')).values
        meta.index = pd.Series(meta.index).apply(lambda x: x.replace(' ', '.')).values
        meta['index'] = meta.index
        meta['val'] = 1
        meta = meta.pivot_table(values='val', index='index', columns='pathway')
        meta[meta.isna()] = 0
        meta = meta.astype(int)
    else:
        meta = meta.T
    for set_group in ['single', 'overlap']:
        score_list = os.listdir('../results/benchmark/scores/%s/%s/' % (dname, set_group))
        score_list = [x for x in score_list if x[-4:]=='.csv']
        for score_name in score_list:
            print(dname, set_group, score_name)
            score = pd.read_csv('../results/benchmark/scores/%s/%s/%s' % (dname, set_group,score_name),
                                sep=',', header=0, index_col=0)
            if dname == 'gdsc':
                score.columns = pd.Series(score.columns).apply(lambda x: x[1:]).values
            score = score[meta.index].T
            results = pd.DataFrame(index=score.columns, columns=meta.columns)
            for group in results.columns:
                y_tr = meta[group]
                results[group] = score.apply(my_rocauc, y_tr=y_tr)
            results.to_csv('../results/benchmark/rocaucs/%s/%s/%s' % (dname, set_group, score_name), sep=',')
    results = pd.DataFrame(index=range(1000), columns=meta.columns)
    score = pd.DataFrame(np.random.uniform(size=(len(score.index), 1000)), index=score.index, columns=range(1000))
    for group in results.columns:
        y_tr = meta[group]
        results[group] = score.apply(my_rocauc, y_tr=y_tr)
    results.to_csv('../results/benchmark/rocaucs/%s/random.csv' % dname, sep=',')




