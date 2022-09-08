import numpy as np
import pandas as pd
from collections import Counter
import os


pt_names = pd.read_csv('fnames',header=None)
pt_names = list(pt_names[0])
nb_samples = len(pt_names)

result = np.zeros((nb_samples, 8))
found_names = []
for pt_ix,pt_name in enumerate(pt_names):
    fname = f'{pt_name}_withmh.csv'
    if fname in os.listdir('processed_DJ'):
        data = pd.read_csv(f'processed_DJ/{pt_name}_withmh.csv')
        cnt = Counter(data['mismatch_mh'])
        total = np.sum(list(cnt.values()))
        found_names.append(pt_name)
        print(pt_ix)
        for i in range(8):
            result[pt_ix,i] = cnt[i]/total

result = pd.DataFrame(result[:len(found_names),:])
result.columns = np.arange(8)
result.index = found_names

result.to_csv('DJ_junctionMH_freq.csv')
