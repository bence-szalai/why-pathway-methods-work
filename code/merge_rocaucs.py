#%%
import pandas as pd
import numpy as np
import os

#%%
set_names = ['BEST_dorothea_AB', 'BEST_dorothea_CD', 'KEGG', 'REACTOME', 'BIOCARTA', 'CGP']
set_groups = ['filtered', 'random_uniform', 'random_dist']
overlaps = ['overlap', 'minus']

#%%
for benchmark in ['progeny', 'gdsc']:
    all_files = os.listdir('../results/benchmark/rocaucs/%s/overlap/' % benchmark)
    for sn1 in set_names:
        for sn2 in set_names:
            for sg in set_groups:
                if sg == 'filtered':
                    for ovrlp in overlaps:
                        for ovrlp in overlaps:
                            fnames = []
                            prefix = '_'.join([sn1, sg, ovrlp, sn2, sg])
                            for f in all_files:
                                if (prefix in f) & ('abs' in f):
                                    fnames.append(f)
                            if len(fnames) > 0:
                                data = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                   fnames[0]),
                                                   sep=',', header=0, index_col=0)
                                for fn in fnames[1:]:
                                    data_tmp = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                           fn),
                                                           sep=',', header=0, index_col=0)
                                    data = pd.concat([data, data_tmp], 0)
                                data.to_csv('../results/benchmark/rocaucs/%s/overlap_merged/%s_abs.csv' % (benchmark,
                                                                                                       prefix), sep=',')
                            all_files = list(set(all_files) - set(fnames))
                            fnames = []
                            prefix = '_'.join([sn1, sg, ovrlp, sn2, sg])
                            for f in all_files:
                                if (prefix in f) & ('abs' not in f):
                                    fnames.append(f)
                            if len(fnames) > 0:
                                data = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                   fnames[0]),
                                                   sep=',', header=0, index_col=0)
                                for fn in fnames[1:]:
                                    data_tmp = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                           fn),
                                                           sep=',', header=0, index_col=0)
                                    data = pd.concat([data, data_tmp], 0)
                                data.to_csv('../results/benchmark/rocaucs/%s/overlap_merged/%s.csv' % (benchmark,
                                                                                                       prefix), sep=',')
                            all_files = list(set(all_files) - set(fnames))

                else:
                    for i in range(5):
                        for ovrlp in overlaps:
                            fnames = []
                            prefix = '_'.join([sn1, sg, str(i), ovrlp, sn2, sg, str(i)])
                            for f in all_files:
                                if (prefix in f) & ('abs' in f):
                                    fnames.append(f)
                            if len(fnames) > 0:
                                data = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                   fnames[0]),
                                                   sep=',', header=0, index_col=0)
                                for fn in fnames[1:]:
                                    data_tmp = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                           fn),
                                                           sep=',', header=0, index_col=0)
                                    data = pd.concat([data, data_tmp], 0)
                                data.to_csv('../results/benchmark/rocaucs/%s/overlap_merged/%s_abs.csv' % (benchmark,
                                                                                                       prefix), sep=',')
                            all_files = list(set(all_files) - set(fnames))
                            fnames = []
                            prefix = '_'.join([sn1, sg, str(i), ovrlp, sn2, sg, str(i)])
                            for f in all_files:
                                if (prefix in f) & ('abs' not in f):
                                    fnames.append(f)
                            if len(fnames) > 0:
                                data = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                   fnames[0]),
                                                   sep=',', header=0, index_col=0)
                                for fn in fnames[1:]:
                                    data_tmp = pd.read_csv('../results/benchmark/rocaucs/%s/overlap/%s' % (benchmark,
                                                                                                           fn),
                                                           sep=',', header=0, index_col=0)
                                    data = pd.concat([data, data_tmp], 0)
                                data.to_csv('../results/benchmark/rocaucs/%s/overlap_merged/%s.csv' % (benchmark,
                                                                                                       prefix), sep=',')
                            all_files = list(set(all_files) - set(fnames))
    print(all_files)
