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
import matplotlib.pyplot as plt
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_twostep_model_em_classes import TwoStepConditionalLogisticRegressionEvaluatorEM
from jax_twostep_model_classes import TwoStepConditionalLogisticRegressionEvaluator
from config import MOD_OUTPUT_PATH, MOD_PROJECT_PATH
import variable_configuration

ANNOTATION_TYPE_ROOT = 'igor_mh_sim_alpha_from_motif_two-side-base-count-beyond_MHprob'
PARAM_GROUP = 'both_v-j_trim_ligation-mh'
LEFT_NUC_MOTIF_COUNT = 1
RIGHT_NUC_MOTIF_COUNT = 2
MODEL_TYPE = 'motif_two-side-base-count-beyond_mh-config-count_ligation-mh'
L2 = False
NCPU = 2

results = {}

for mh in ['0', '1']:
    for trim_mh in ['0', '1']:
        MHkey = 'ligMH param = ' + mh + ', trimMH param = ' + trim_mh
        results[MHkey] = {}

        ANNOTATION_TYPE = ANNOTATION_TYPE_ROOT + mh + '_trimMHprob' + trim_mh

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

        model_filename = params.model_output_path(L2)
        evaluator = TwoStepConditionalLogisticRegressionEvaluatorEM(model_filename,
                                                           params,
                                                           None)

        results[MHkey]['EM_losses'] = evaluator.model.EM_iteration_losses
        results[MHkey]['EM_errors'] = evaluator.model.EM_iteration_errors

        params.sample_annotation = True


# Use Set2 colormap for color palette
colors = plt.cm.Set2(range(len(results)))

# Iterate over keys of the dictionary
for i, (key, value) in enumerate(results.items()):
    # Extract EM_losses
    EM_losses = value['EM_losses']
    # Plot EM_losses against iteration (index)
    plt.plot(range(len(EM_losses)), EM_losses, label=key, color=colors[i])
    # Plot points
    plt.scatter(range(len(EM_losses)), EM_losses, color=colors[i], zorder=5)  # zorder to plot points on top

# Add labels and legend
plt.xlabel('EM iteration')
plt.ylabel('Loss')
plt.legend(title='Simulator MH involvement')

# Add grid lines
plt.grid(True)

file_name = MOD_PROJECT_PATH + '/plotting_scripts/manuscript_figs/mh_simulator_experiments/compare_EM_noEM_simulator/both_motif/convergence.pdf'

plt.savefig(file_name)
