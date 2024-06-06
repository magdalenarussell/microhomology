import sys
sys.path.append('/home/mrussel2/microhomology/mechanistic-trimming/jax_scripts/')
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

# Parameters for data generation
N = 100  # Number of groups
K = 4 # Number of choice 1 options
J = 3 # Number of choice 2 options
beta = np.array([0.51, -0.89, 0.34])  # Known coefficients for predictors (first two for choice 1 and last one for choice 2)

# Generate predictors
X1 = np.random.randn(N, K, len(beta) - 1)
X2 = np.random.randn(N, K, J, 1)
mask = np.ones((N, K, J, 1))

# remove some choices
# Specify the fraction of elements to set to zero
fraction_zeros = 0.2

# Calculate the number of elements to set to zero
total_elements = mask.size
num_zeros = int(fraction_zeros * total_elements)

# Randomly choose indices to set to zero
indices = np.random.choice(total_elements, num_zeros, replace=False)

# Convert linear indices to multi-dimensional indices and set those positions to zero
np.put(mask, indices, 0)
np.put(mask, np.array([i for i in range(0, 24)]), 0)


# Generate outcome based on the predictors and known coefficients
model = TwoStepConditionalLogisticRegressor(training_df = None,
                                           variable_colnames = ['b1', 'b2', 'a1'],
                                           choice1_variable_colnames = ['b1', 'b2'],
                                           choice2_variable_colnames = ['a1'],
                                           count_colname = None,
                                           group_colname = None,
                                           repeat_obs_colname = None,
                                           choice_colname = None,
                                           choice2_colname = None,
                                           params = None)

print('initialized model')

# get probabilities
probs = np.array(model.get_joint_prob(X1, X2, mask, beta))

# Initialize an array to count selections
selection_counts = np.zeros_like(probs)

# Number of simulations
num_simulations = 1000
num_matrices, rows, cols = probs.shape

for simulation in range(num_simulations):
    for matrix_index in range(num_matrices):
        # Flatten the matrix to make the entire matrix work as a single probability distribution
        flat_probs = probs[matrix_index].flatten()

        if np.isnan(flat_probs.sum()):
            continue

        # Ensure probabilities sum to 1 (if needed due to precision issues)
        flat_probs /= flat_probs.sum()

        # Sample a single position from the flattened probability distribution
        sampled_flat_index = np.random.choice(np.arange(rows * cols), p=flat_probs)

        # Convert the flat index back to a 2D index
        sampled_row, sampled_col = divmod(sampled_flat_index, cols)

        # Increment the count for the selected position
        selection_counts[matrix_index, sampled_row, sampled_col] += 1

# get observation proportions
counts_matrix = selection_counts/num_simulations

# set variables
model.choice1_variable_matrix = X1
model.choice2_variable_matrix = X2
model.mask_matrix = mask
model.counts_matrix = counts_matrix

# train model
model = model.train_model(l2=False, maxiter=100, tolerance=1e-8)

# compare coefficients
print("Data generating parameters: ", beta)
print("Inferred model coefficients: ", np.array(round(model.coefs.flatten(), 2)))
print("Equivalent?", all(beta == np.array(round(model.coefs.flatten(), 2))))
