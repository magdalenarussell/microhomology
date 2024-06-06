import pandas as pd
import numpy as np
import pytest
from jax_model_classes import ConditionalLogisticRegressor  # Replace 'your_module' with the actual module name

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

# Instantiate the ConditionalLogisticRegressor class for testing
logistic_regressor = ConditionalLogisticRegressor(df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

def test_get_random_coefs():
    coefs = logistic_regressor.get_random_coefs()
    assert coefs.shape == (len(logistic_regressor.variable_colnames), 1)

def test_train_model():
    trained_model = logistic_regressor.train_model(l2=True, maxiter=1000, stepsize=0.001)
    assert trained_model.coefs is not None
    assert (trained_model.coefs != trained_model.initial_coefs).all()
    assert trained_model.l2reg != 0
    assert trained_model.maxiter == 1000
    assert trained_model.stepsize == 0.001

def test_grid_search_cv():
    results = logistic_regressor.grid_search_cv(2, [0.1, 0.5, 1.0])
    assert len(results) == len([0.1, 0.5, 1.0])
    assert 'l2reg' in results.columns
    assert 'mean_CV_loss' in results.columns

def test_get_l2reg():
    l2reg = logistic_regressor.get_l2reg()
    assert isinstance(l2reg, float)

