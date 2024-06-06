import pandas as pd
import numpy as np
import pytest
from jax_model_classes import DataTransformer  # Replace 'your_module' with the actual module name

# Define some sample data for testing
data = {
    'Variable1': ['A', 'B', 'A', 'C'],
    'Variable2': [1, 2, 3, 4],
    'Count': [10, 20, 30, 40],
    'Group': ['Group1', 'Group1', 'Group2', 'Group2'],
    'RepeatObs': ['Obs1', 'Obs1', 'Obs1', 'Obs1'],
    'Choice': ['Choice1', 'Choice2', 'Choice1', 'Choice2']
}


df = pd.DataFrame(data)

# Instantiate the DataTransformer class for testing
data_transformer = DataTransformer(df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

def test_get_random_coefs():
    coefs = data_transformer.get_random_coefs()
    assert coefs.shape == (len(data_transformer.variable_colnames), 1)

def test_preprocess_data():
    preprocessed_df = data_transformer.preprocess_data()
    assert 'Variable1_A' in preprocessed_df.columns
    assert 'Group_int_RepeatObs_int' in preprocessed_df.columns
    assert 'Choice_int' in preprocessed_df.columns

def test_get_matrices():
    var_matrix, counts_matrix, nonrep_grp_matrix = data_transformer.get_matrices()
    df = data_transformer.training_df
    assert var_matrix.shape == (len(df['Group_int_RepeatObs_int'].unique()), len(df['Choice'].unique()), len(data_transformer.variable_colnames))
    assert counts_matrix.shape == (len(df['Group_int_RepeatObs_int'].unique()), len(df['Choice'].unique()), 1)
    assert nonrep_grp_matrix.shape == (len(df['Group_int_RepeatObs_int'].unique()), 1)

def test_get_coefficients():
    coefs = np.array([[0.5], [0.3], [0.2]])
    data_transformer.coefs = coefs
    coef_dict = data_transformer.get_coefficients()
    assert coef_dict['Variable1_A'] == 0.3
    assert 'Variable1_C' in coef_dict
    assert coef_dict['Variable2'] == 0.5

