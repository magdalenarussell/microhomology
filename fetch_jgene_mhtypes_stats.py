import numpy as np
from collections import Counter
from tqdm import tqdm
import pandas as pd
import sys

selection = sys.argv[1]


all_counts = []

fnames = pd.read_csv('fnames',header=None)
fnames = list(fnames[0])

for ptname in tqdm(fnames):

    data = pd.read_csv(f'processed_DJ/{ptname}_withmh.csv')
    data['new_j_trim'] = data['j_trim']-pd.Series([min(max(0,i),2) for i in data['j_pnucs']])
    data['new_d1_trim'] = data['d1_trim']-pd.Series([min(max(0,i),2) for i in data['d1_pnucs']])

    if selection == 'productive':
        data = data[data['productive']]
    elif selection == 'non-productive':
        data = data[~data['productive']]
    elif selection == 'no-insertions':
        data = data[data['dj_insert']==0]


    d = 'TRBD1*01'
    mh_trimming_groups = {}
    for j in list(set(data['j_gene'])):

        temp = data[data['j_gene']==j]
        temp = temp[temp['d_gene']==d]
        t = temp[temp['new_j_trim']>=0]
        if Counter(t['mismatch_mh'])[0]/np.sum(list(Counter(t['mismatch_mh']).values()))>0.1:
            mh_trimming_groups[j]='type1'
        else:
            mh_trimming_groups[j]='type2'
    data['j_mh_group'] = data['j_gene'].replace(mh_trimming_groups)
    toplot = pd.DataFrame(data.groupby('j_mh_group').count()['j_gene'])
    if len(all_counts)==0:
        all_counts = pd.DataFrame(data.groupby('j_mh_group').count()['j_gene'])
        all_counts.columns = [f'{ptname}']
    else:
        toplot = pd.DataFrame(data.groupby('j_mh_group').count()['j_gene'])
        toplot.columns = [f'{ptname}']
        all_counts = all_counts.merge(toplot, left_index=True, right_index=True)

all_counts.to_csv(f'allpt_counts_jgeneMH_{selection}.csv')
