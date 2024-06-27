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
from jax_twostep_model_em_classes import TwoStepDataTransformerEM, TwoStepConditionalLogisticRegressorEM, TwoStepConditionalLogisticRegressionPredictorEM, TwoStepConditionalLogisticRegressionEvaluatorEM
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
PARAM_GROUP_VALIDATION = sys.argv[9]

# initialize parallelized pandas
pandarallel.initialize(nb_workers=NCPU, progress_bar=True)

# set global variables
trained_params = variable_configuration.global_paramaters(MOD_OUTPUT_PATH,
                                                  MOD_PROJECT_PATH,
                                                  ANNOTATION_TYPE,
                                                  PARAM_GROUP,
                                                  LEFT_NUC_MOTIF_COUNT,
                                                  RIGHT_NUC_MOTIF_COUNT,
                                                  MODEL_TYPE)


val_params = variable_configuration.global_paramaters(MOD_OUTPUT_PATH,
                                                  MOD_PROJECT_PATH,
                                                  ANNOTATION_TYPE_VALIDATION,
                                                  PARAM_GROUP_VALIDATION,
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

assert trained_params.sample_annotation == False, "Annotations must not be sampled"

# read in data
processed_data_filename = trained_params.R_processed_data_path(annotation = ANNOTATION_TYPE_VALIDATION, param_group = PARAM_GROUP_VALIDATION)
processed_data = pd.read_csv(processed_data_filename, sep = '\t')
print('read in data')

if ANNOTATION_TYPE_VALIDATION == ANNOTATION_TYPE:
    train_df=processed_data
else:
    train_df=None

# validate model
model_filename = trained_params.model_output_path(L2)
evaluator = TwoStepConditionalLogisticRegressionEvaluatorEM(model_filename,
                                                   val_params,
                                                   training_df=train_df,
                                                   validation_df=processed_data)

result = evaluator.compile_evaluation_results_df(ANNOTATION_TYPE,
                                                 trained_params.productivity,
                                                 True)

result['validation_data_type'] = ANNOTATION_TYPE_VALIDATION

# write predictions and coefficients
path = trained_params.model_eval_results_path(L2)

if os.path.isfile(path):
    # Read the file as a Pandas DataFrame
    df = pd.read_csv(path, sep='\t')
    result = pd.concat([df, result], axis = 0)

result.to_csv(path, sep='\t', index=False)
print('finished validating model')

# make predictions on validation dataset
if ANNOTATION_TYPE_VALIDATION != ANNOTATION_TYPE:
    predictor = TwoStepConditionalLogisticRegressionPredictorEM(model=evaluator.model,
                                                       variable_colnames = model_params.variable_colnames,
                                                       choice1_variable_colnames = model_params.choice1_variable_colnames,
                                                       choice2_variable_colnames = model_params.choice2_variable_colnames,

                                                       count_colname = model_params.count_colname,
                                                       group_colname = model_params.group_colname,
                                                       repeat_obs_colname = model_params.repeat_obs_colname,
                                                       choice_colname = model_params.choice_colname,
                                                       choice2_colname = model_params.choice2_colname,
                                                       params = trained_params)

    # write predictions and coefficients
    training_pred = predictor.predict(new_df=processed_data)
    predictions_filename = trained_params.validation_predictions_data_path(L2, ANNOTATION_TYPE_VALIDATION)
    training_pred.to_csv(predictions_filename, sep='\t', index=False)
    print('finished making predictions with validation data')

print('done')
