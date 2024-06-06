# Analysis scripts

This directory contains files for running the all of the analyses: 

The main scripts are contained in the following files:

* [Data compilation functions](data_compilation_functions.R) which are sourced by nearly all top-level analysis scripts
* [Model training script](fit_model.R) which is sourced by top level model training script (e.g. ../fit_model.sh)
* [Coefficient significance script](evaluate_coeffs.R) which is sourced by top level script to quantify the significance of model coefficients (e.g. ../evaluate_coef_significance.sh)
* [Model evaluation script](evaluate_models.R) which is sourced by top level model evaluation script (e.g. ../evaluate_model.sh)
* [Model validation script](validate_model.R) which is sourced by top level model validation script (e.g. ../validate_model.sh)

Additional analysis scripts are contained in the [analysis scripts](analysis_scripts) directory. See [README](analysis_scripts/README.md) for more details.
IGoR annotation specific scripts are contained within the [igor scripts](igor_annotation_scripts) directory. See [README](igor_annotation_scripts/README.md) for more details.
General [model evaluation](model_evaluation_functions.R) and [model training](model_fitting_functions.R) functions are also located here.

Argument specific functions are located in the other directories housed here:

* [annotation_specific_functions](annotation_specific_functions) which are specific for the type of sequence annotation used for the data  
* [data grouping functions](data_grouping_functions) which are specific to how the data should be grouped 
* [gene specific functions](gene_specific_functions) which are specific to whether V-genes or J-genes are being modeled
* [model evaluation type functions](model_evaluation_type_functions) which are specific to the type of model evaluation being used
* [model formula functions](model_formula_functions) which are specific to the model being used
* [model group functions](model_group_functions) which are specific to whether the model is trained using all individuals, are separately for each individual
* [motif class functions](motif_class_functions) which are specific to the specified hairpin-opening position
* [sampling procedure functions](sampling_procedure_functions) which are specific to how the weighted likelihood function is formulated

