# Data processing scripts

This directory contains files for running the all of the data processing and supplementary analyses: 

The main scripts are contained in the following files:

* [Data compilation functions](data_compilation_functions.R) which are sourced by nearly all top-level data processing scripts
* [Addition MH related data processing functions](mh_functions.R) which are sourced by nearly all top-level data processing scripts
* [Data processing script](process_data_for_model_fitting.R) which actually processes the data for model training
* [Bootstrapped data processing script](process_bootstrap_datasets.R) which generates bootstrapped datasets for supplementary analyses
* [Coefficient significance script](run_coef_signif_test.R) which runs coefficient significance tests using inferred coefficients from bootstrapped dataset models (requires training models beforehand)
* [Likelihood ratio test script](run_lrt.R) which runs the likelihood ratio test using pre-trained models

Model training scripts are contained in the [jax_scripts](../jax_scripts) directory. See [README](../jax_scripts/README.md) for more details.

MH simulator scripts are contained in the [mh simulation scripts](../mh_simulation_scripts) directory. See [README](../mh_simulation_scripts/README.md) for more details.

IGoR annotation specific scripts are contained within the [igor scripts](../igor_annotation_scripts) directory. See [README](../igor_annotation_scripts/README.md) for more details.

Argument specific functions are located in the other directories housed here:

* [annotation_specific_functions](annotation_specific_functions) which are specific for the type of sequence annotation software and locus sequenced for the data  
* [data type functions](data_type_functions) which are specific to the annotation software used 
* [locus specific functions](locus_specific_functions) which are specific to the gene locus being modeled
* [gene specific functions](gene_specific_functions) which are specific to whether V-genes or J-genes are being modeled
* [gene count specific functions](gene_count_specific_functions) which are specific to whether a locus with two junctions or one junction is being modeled
* [parameter group functions](param_groups) which are specific to the modeling parameters being used
* [model formula functions](model_formula_functions) which are specific to the model being used
* [motif class functions](motif_class_functions) which are specific to the specified hairpin-opening position
* [sampling procedure functions](sampling_procedure_functions) which are specific to how the weighted likelihood function is formulated
