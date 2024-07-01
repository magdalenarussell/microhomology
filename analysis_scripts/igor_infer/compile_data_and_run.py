import numpy as np
import pandas as pd
import pygor3 as p3
import sys
import os
import subprocess
import logging
logging.getLogger().setLevel(logging.INFO)

file = sys.argv[1]
script_wd = sys.argv[2]
fname = os.path.basename(file)
data = pd.read_csv(file,sep='\t')
print (data.shape)

### From here on, we use the new stitched sequences to infer the IGoR model.
ptname = fname.split('.')[0]
igor_working_directory = f'./{ptname}_igor_infer'
os.mkdir(igor_working_directory)
logging.info('Running IGoR to obtain CDR3 nucleotide sequences.')
igor_infile = os.path.join(igor_working_directory, fname + '_align.csv')
data['V1'].to_csv(igor_infile, sep=';',
                  index_label='seq_id', header=['sequence'])
igor_command = f'bash {script_wd}/analysis_scripts/igor_infer/run_igor.sh {igor_working_directory} {igor_infile} batchname '
subprocess.call(igor_command.split(' '))
