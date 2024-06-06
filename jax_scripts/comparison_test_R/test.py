import sys
sys.path.append('/home/mrussel2/microhomology/jax_scripts/')
import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import ConditionalLogisticRegressor, ConditionalLogisticRegressionPredictor

ld = pd.read_csv('comparison_test_R/long_test_data.tsv', sep = '\t')

# initialize model 
model = ConditionalLogisticRegressor(training_df = ld,
                                     variable_colnames=['height','width', 'rating'],
                                     count_colname='count',
                                     group_colname='color',
                                     repeat_obs_colname='garden',
                                     choice_colname='species')

# train model
model = model.train_model(l2=False)

# example of getting results
predictor = ConditionalLogisticRegressionPredictor(model=model,
                                                   new_df=ld,
                                                   variable_colnames=['height','width', 'rating'],
                                                   count_colname='count',
                                                   group_colname='color',
                                                   repeat_obs_colname='garden',
                                                   choice_colname='species')

y_hat = predictor.predict()
print(y_hat)
loss = predictor.compute_loss()
print(loss)
c = predictor.get_coefficients()
print(c)

