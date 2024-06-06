import sys
sys.path.append('/home/mrussel2/microhomology/jax_scripts/')
import os
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
from jax_twostep_model_classes import TwoStepDataTransformer, TwoStepConditionalLogisticRegressor, TwoStepConditionalLogisticRegressionPredictor
from jax_twostep_model_em_classes import TwoStepDataTransformerEM, TwoStepConditionalLogisticRegressorEM
from config import MOD_OUTPUT_PATH, MOD_PROJECT_PATH
import variable_configuration

ANNOTATION_TYPE = sys.argv[1]
PARAM_GROUP = sys.argv[2]
LEFT_NUC_MOTIF_COUNT = int(sys.argv[3])
RIGHT_NUC_MOTIF_COUNT = int(sys.argv[4])
MODEL_TYPE = sys.argv[5]
L2 = sys.argv[6]
L2 = (L2.lower() == 'true')
NCPU = int(sys.argv[7])

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

params.sample_annotation = True

assert params.sample_annotation == True, "Annotations must be sampled for each observed sequence"

# set model type specific parameters
model_params = variable_configuration.model_specific_parameters(PARAM_GROUP,
                                                                MODEL_TYPE,
                                                                LEFT_NUC_MOTIF_COUNT,
                                                                RIGHT_NUC_MOTIF_COUNT)
model_params = model_params.process_model_parameters()
print('loaded parameters')

# read in data
processed_data_filename = params.R_processed_data_path()
processed_data = pd.read_csv(processed_data_filename, sep = '\t')
print('read in data')

# initialize model 
model = TwoStepConditionalLogisticRegressor(training_df = processed_data,
                                           variable_colnames = model_params.variable_colnames,
                                           choice1_variable_colnames = model_params.choice1_variable_colnames,
                                           choice2_variable_colnames = model_params.choice2_variable_colnames,
                                           count_colname = model_params.count_colname,
                                           group_colname = model_params.group_colname,
                                           repeat_obs_colname = model_params.repeat_obs_colname,
                                           choice_colname = model_params.choice_colname,
                                           choice2_colname = model_params.choice2_colname,
                                           params = params)
print('initialized model')


# train model
model = model.train_model(l2=L2, maxiter=250, tolerance=1e-8)
print('trained model')

# make predictions on full training dataset
# predictor = TwoStepConditionalLogisticRegressionPredictor(model=model,
#                                                       variable_colnames = model_params.variable_colnames,
#                                                       choice1_variable_colnames = model_params.choice1_variable_colnames,
#                                                       choice2_variable_colnames = model_params.choice2_variable_colnames,
#                                                       count_colname = model_params.count_colname,
#                                                       group_colname = model_params.group_colname,
#                                                       repeat_obs_colname = model_params.repeat_obs_colname,
#                                                       choice_colname = model_params.choice_colname,
#                                                       choice2_colname = model_params.choice2_colname,
#                                                       params = params)

# write predictions and coefficients
# training_pred = predictor.predict(new_df=processed_data)
# predictions_filename = params.predictions_data_path(L2)
# training_pred.to_csv(predictions_filename, sep='\t', index=False)

if MODEL_TYPE != 'null':
    coefs = model.get_coefficients_df()
    coefs['training_error'] = float(model.training_info.state.error)
    coefs_filename = params.trained_coefs_path(L2)
    coefs.to_csv(coefs_filename, sep='\t', index=False)

print('finished processing model predictions')

# save trained model
model_filename = params.model_output_path(L2)
model.save_model(model_filename)
