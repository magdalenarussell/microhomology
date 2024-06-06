import pandas as pd
import numpy as np
import pytest
from jax_model_classes import ConditionalLogisticRegressionPredictor  # Import your class

# Create a fixture to set up a sample model and data for testing
@pytest.fixture
def sample_data_and_model():
    # Create a sample DataFrame for testing
    data = {
        'Variable1': ['A', 'B', 'A', 'C'],
        'Variable2': [1, 2, 3, 4],
        'Count': [10, 20, 30, 40],
        'Group': ['Group1', 'Group1', 'Group2', 'Group2'],
        'RepeatObs': ['Obs1', 'Obs1', 'Obs1', 'Obs1'],
        'Choice': ['Choice1', 'Choice2', 'Choice1', 'Choice2']
    }

    df = pd.DataFrame(data)

    # Create a sample model for testing
    from jax_model_classes import ConditionalLogisticRegressor  # Import your model class
    model = ConditionalLogisticRegressor(df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')
    model.train_model()

    return df, model

# Test the predict method of ConditionalLogisticRegressionPredictor
def test_predict(sample_data_and_model):
    df, model = sample_data_and_model

    # Create a predictor instance
    predictor = ConditionalLogisticRegressionPredictor(model, df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

    # Make predictions
    predictions = predictor.predict()

    # Ensure that predictions is a DataFrame
    assert isinstance(predictions, pd.DataFrame)

    # Check if the predicted probabilities are within the valid range (0 to 1)
    assert (predictions['predicted_prob'] >= 0).all()
    assert (predictions['predicted_prob'] <= 1).all()

# Test the compute_loss method of ConditionalLogisticRegressionPredictor
def test_compute_loss(sample_data_and_model):
    df, model = sample_data_and_model

    # Create a predictor instance
    predictor = ConditionalLogisticRegressionPredictor(model, df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

    # Compute the loss
    loss = predictor.compute_loss()

    assert isinstance(loss, float)

# Test the get_coefficients method of ConditionalLogisticRegressionPredictor
def test_get_coefficients(sample_data_and_model):
    df, model = sample_data_and_model

    # Create a predictor instance
    predictor = ConditionalLogisticRegressionPredictor(model, df, ['Variable1', 'Variable2'], 'Count', 'Group', 'RepeatObs', 'Choice')

    # Get the coefficients as a dictionary
    coefficients = predictor.get_coefficients()

    # Ensure that coefficients is a dictionary
    assert isinstance(coefficients, dict)

    # Check if the values are numbers
    assert all(isinstance(value, (int, float)) for value in coefficients.values())

# Run pytest to execute the tests
if __name__ == "__main__":
    pytest.main()

