import sys
sys.path.append('/home/mrussel2/microhomology/jax_scripts/')
import os
import glob
import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import importlib
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import ConditionalLogisticRegressor, ConditionalLogisticRegressionPredictor
from config import MOD_OUTPUT_PATH, MOD_PROJECT_PATH
import variable_configuration

ANNOTATION_TYPE = sys.argv[1]
PARAM_GROUP = sys.argv[2]
LEFT_NUC_MOTIF_COUNT = int(sys.argv[3])
RIGHT_NUC_MOTIF_COUNT = int(sys.argv[4])
MODEL_TYPE = sys.argv[5]
L2 = sys.argv[6]
L2 = (L2.lower() == 'true')
L2reg = float(sys.argv[7])
PROP = sys.argv[8]
NCPU = int(sys.argv[9])

# initialize parallelized pandas
pandarallel.initialize(nb_workers=NCPU, progress_bar=True)

# set global variables
params = variable_configuration.global_paramaters(MOD_OUTPUT_PATH,
                                                  MOD_PROJECT_PATH,
                                                  ANNOTATION_TYPE,
                                                  PARAM_GROUP,
                                                  LEFT_NUC_MOTIF_COUNT,
                                                  RIGHT_NUC_MOTIF_COUNT,
                                                  MODEL_TYPE)

# set model type specific parameters
model_params = variable_configuration.model_specific_parameters(PARAM_GROUP,
                                                                MODEL_TYPE,
                                                                LEFT_NUC_MOTIF_COUNT,
                                                                RIGHT_NUC_MOTIF_COUNT)
model_params = model_params.process_model_parameters()
print('loaded parameters')

# read in data
processed_data_path = params.R_subsampling_processed_data_path(PROP)
pattern = f"{processed_data_path}/*.tsv"
files_list = glob.glob(pattern)

all_coefs = pd.DataFrame()

for file in files_list:
    index = files_list.index(file)
    print(f"starting experiment {index} of {len(files_list)}")
    processed_data = pd.read_csv(file, sep = '\t')
    print('read in data')

    # initialize model 
    model = ConditionalLogisticRegressor(training_df = processed_data,
                                         variable_colnames = model_params.variable_colnames,
                                         count_colname = model_params.count_colname,
                                         group_colname = model_params.group_colname,
                                         repeat_obs_colname = model_params.repeat_obs_colname,
                                         choice_colname = model_params.choice_colname)
    print('initialized model')

    # train model
    model = model.train_model(l2=L2, l2reg_value=L2reg, maxiter=10000, tolerance=1e-8)
    print('trained model')

    # get coefficients
    coefs = model.get_coefficients_df()
    coefs['iteration'] = index
    all_coefs = pd.concat([all_coefs, coefs], ignore_index=True)
    print('finished processing model coefficients')

coefs_filename = params.subsampling_coefs_path(PROP, L2)
all_coefs.to_csv(coefs_filename, sep='\t', index=False)

# removing all temporary files
for file in files_list:
    os.remove(file)
