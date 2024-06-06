import pandas as pd
import numpy as np
import pytest
from jax_model_classes import DataPreprocessor

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

# Instantiate the DataPreprocessor class for testing
data_preprocessor = DataPreprocessor(df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

def test_get_mapping_dict():
    mapping_dict = data_preprocessor.get_mapping_dict(df, 'Variable1')
    assert mapping_dict == {'A': 0, 'B': 1, 'C': 2}

def test_transform_categorical_response_vars():
    transformed_df = data_preprocessor.transform_categorical_response_vars(df, 'group_colname', 'group_colname')
    assert 'Group_int' in transformed_df.columns

def test_get_contrast_matrix():
    contrast_matrix = data_preprocessor.get_contrast_matrix(df, 'Variable1')
    assert 'Variable1_A' in contrast_matrix.columns
    assert 'Variable1_B' in contrast_matrix.columns

def test_get_dropped_contrast_var():
    dropped_vars, missing_var = data_preprocessor.get_dropped_contrast_var(df, 'Variable1')
    assert 'Variable1_C' not in dropped_vars
    assert missing_var == 'Variable1_C'

def test_transform_categorical_vars():
    transformed_df = data_preprocessor.transform_categorical_vars(df, 'Variable1')
    assert 'Variable1_A' in transformed_df.columns
    assert 'Variable1_A' in data_preprocessor.variable_colnames
    assert 'Variable1' not in data_preprocessor.variable_colnames

def test_expand_multivariable():
    data_preprocessor.group_colname = ['Group', 'RepeatObs']
    expanded_df = data_preprocessor.expand_multivariable(df, 'group_colname', 'new_group_colname')
    assert 'Group_RepeatObs' in expanded_df.columns
    assert 'Group' not in [data_preprocessor.new_group_colname]

