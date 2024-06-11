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
ANNOTATION_TYPE_VALIDATION = sys.argv[8]
CALCULATE_EXPECTED_LOSS = sys.argv[9]
CALCULATE_EXPECTED_LOSS = (CALCULATE_EXPECTED_LOSS.lower() == 'true')

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

assert params.sample_annotation == True, "Annotations must be sampled for each observed sequence"

# read in data
processed_data_filename = params.R_processed_data_path(annotation = ANNOTATION_TYPE_VALIDATION)
processed_data = pd.read_csv(processed_data_filename, sep = '\t')
print('read in data')

if ANNOTATION_TYPE_VALIDATION == ANNOTATION_TYPE:
    train_df=processed_data
else:
    train_df=None

# validate model
model_filename = params.model_output_path(L2)
evaluator = TwoStepConditionalLogisticRegressionEvaluator(model_filename,
                                                   params,
                                                   training_df=train_df,
                                                   validation_df=processed_data)

result = evaluator.compile_evaluation_results_df(True,
                                                 CALCULATE_EXPECTED_LOSS)

result['validation_data_type'] = ANNOTATION_TYPE_VALIDATION

# write predictions and coefficients
path = params.model_eval_results_path(L2)

if os.path.isfile(path):
    # Read the file as a Pandas DataFrame
    df = pd.read_csv(path)
    result = pd.concat([df, result], axis = 0)

result.to_csv(path, sep='\t', index=False)
print('finished validating model')

# make predictions on validation dataset
if ANNOTATION_TYPE_VALIDATION != ANNOTATION_TYPE:
    predictor = TwoStepConditionalLogisticRegressionPredictor(model=evaluator.model,
                                                       variable_colnames = model_params.variable_colnames,
                                                       count_colname = model_params.count_colname,
                                                       group_colname = model_params.group_colname,
                                                       repeat_obs_colname = model_params.repeat_obs_colname,
                                                       choice_colname = model_params.choice_colname,
                                                       params = params)


    # write predictions and coefficients
    training_pred = predictor.predict(new_df=processed_data)
    predictions_filename = params.validation_predictions_data_path(L2)
    training_pred.to_csv(predictions_filename, sep='\t', index=False)
    print('finished making predictions with validation data')

print('done')
