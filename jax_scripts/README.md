# Model training scripts

This directory contains files for training and validating all models: 

The main scripts are contained in the following files:

* Model training scripts:
    * [Train trimming model](fit_model.py): Note that this model does not infer ligation related signals and does not infer ligation configuration
    * [Train trimming and ligation model without EM](fit_twostep_model.py): only valid when using sequences without N-insertions where the ligation configuration can be determined using the germline sequences; training without EM is only valid when the true trimming and ligation configuration is known (e.g. simulated data)
    * [Train trimming and ligation model with EM](fit_twostep_model_EM.py): only valid when using sequences without N-insertions where the ligation configuration can be determined using the germline sequences; training with EM is what should be done with real data

* Model class files:
    * [Trimming model classes](jax_model_classes.py)
    * [Trimming and ligation model classes](jax_twostep_model_classes.py)
    * [Trimming and ligation model with EM classes](jax_twostep_model_em_classes.py)
    * [Model variable configuration classes](variable_configuration.py)

* Model validation scripts: can be used for making predictions and evaluating loss on an unseen dataset using a pre-trained model
    * [Validate using a trimming model](validate_model.py)
    * [Validate using a trimming and ligation model (trained without EM)](validate_twostep_model.py)
    * [Validate using a trimming and ligation model (trained with EM)](validate_twostep_model_EM.py)

* Model prediction scripts: can be used for making predictions on the training dataset using a pre-trained model (most of the time you should just use the validation scripts)
    * [Predict using a trimming model](predict.py)
    * [Predict using a trimming and ligation model (trained without EM)](predict_twostep.py)
    * [Predict using a trimming and ligation model (trained with EM)](predict_twostep_EM.py)

Note: all data should be processed prior to model training. See [data processing directory](../scripts) for more detail.

Argument specific functions are located in the other directories housed here:

* [annotation configuration files](annotation_configs) which are specific for the type of sequence annotation software and locus sequenced for the data  
* [parameter group files](param_group_configs) which are specific to the modeling parameters being used
* [model type files](model_type_configs) which are specific to the model being used
