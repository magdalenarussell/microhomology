import numpy as np
import pandas as pd
from collections import Counter
import os
from tqdm import tqdm


pt_names = pd.read_csv('fnames',header=None)
pt_names = list(pt_names[0])
nb_samples = len(pt_names)

result_dj = np.zeros((nb_samples, 8))
result_vd = np.zeros((nb_samples, 8))
found_names = []
for pt_ix,pt_name in tqdm(enumerate(pt_names)):
    fname = f'{pt_name}_withmh.csv'
    if fname in os.listdir('processed_both'):
        data = pd.read_csv(f'processed_both/{pt_name}_withmh.csv')
        ## for the DJ junction
        cnt = Counter(data['mismatch_mh_dj'])
        total = np.sum(list(cnt.values()))
        found_names.append(pt_name)
        for i in range(8):
            result_dj[pt_ix,i] = cnt[i]/total

        ## for the VD junction
        cnt = Counter(data['mismatch_mh_vd'])
        total = np.sum(list(cnt.values()))
        found_names.append(pt_name)
        for i in range(8):
            result_vd[pt_ix,i] = cnt[i]/total

result_dj = pd.DataFrame(result_dj[:len(found_names),:])
result_dj.columns = np.arange(8)
result_dj.index = found_names

result_vd = pd.DataFrame(result_vd[:len(found_names),:])
result_vd.columns = np.arange(8)
result_vd.index = found_names

result_dj.to_csv('DJ_junctionMH_freq.csv')
result_vd.to_csv('VD_junctionMH_freq.csv')
