# microhomology
The goal of this project is to use model-based statistical inference to identify the extent of microhomology involvement in nucleotide trimming and ligation processes during V(D)J recomination of adaptive immune receptor loci.


# Install
Everything R 4.1.3 and/or Python 3.11.5 based. Python and R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
conda env create -f environment.yml
conda activate microhomology_jax
```

You will also need to install IGoR if you wish to annotate sequences using IGoR (Marcou et.al Nature Communications 2018)



# Requirements: 
Most of these analyses can be run on any machine.
However, some of the data preparation steps, such as sequence annotation using IGoR (Marcou et.al Nature Communications 2018), are computationally intensive and require a cluster to run efficiently.
This [sequence annotation script](scripts/annotate_with_igor.sh) is written specifically for a cluster set up to use the Slurm job scheduler. 
Some minor modifications could allow this step to be run locally or using a different cluster workload manager. 


# Analysis outline: 

__Table of Contents:__

* [Model training](#model-training)
    * [Summary of parameter groups and model types](#parameter-group-and-model-type-options)
* [Quantifying coefficient significance](#quantifying-coefficient-significance)
* [Model comparison](#model-comparison)
* [Model validation](#model-validation)
    * [Summary of model validation options](#model-validation-options)
* [Plot results](#plot-results)

## Model training

0. Download your training dataset of interest; We are using a dataset consisting of sequences published in these two publications: [here](https://doi.org/10.1016/j.molimm.2020.09.003) and [here](https://doi.org/10.1016/j.jaut.2021.102616)
    1. We downloaded processed repertoire data from these publications from the Adaptive ImmuneACCESS portal: [here](https://doi.org/10.21417/NH2020MI) and [here](https://doi.org/10.21417/NH2021JA)
1. Annotate these sequences using IGoR. You can run the [annotation 
script](igor_annotation_scripts/annotate_with_igor.sh) to do this. As described above, this script is written specifically for a cluster set up to use the Slurm job scheduler. This script takes five arguments:
    1. the raw file path (wherever you downloaded the files in step 0 above)
    2. a temporary directory location (wherever you want the intermediate annotation files to be stored)
    3. a boolean value representing whether the sequences were downloaded from the Adaptive ImmuneACCESS portal (if `TRUE`, the sequences require some extra processing)
    4. the locus of the sequences (e.g. `alpha` or `beta`)
    5. an output directory location (wherever you want the final annotation files to be stored--you will eventually save this location as `TCR_REPERTOIRE_DATA_igor_alpha` within the [config](config/config.R) file in step 3)
    6. the number of possible rearrangement scenarios you want to sample from for each sequence (we used 10 scenarios)
    7. the model parameters you would like to use (e.g. we use the `default` argument) 
    8. the model marginals you would like to use (e.g. we use the `default` argument) 
    9. the number of CPUs you want to use
    10. the cluster partition that you want to submit to
2. Download TRA gene name (or whatever locus you'd like to use) and germline sequences from IMGT (we downloaded these data in September 2023); save these data to a file--you will eventually save the location of this file as `WHOLE_NUCSEQS_alpha` within the [config](config/config.R) file in step 3 
3. Edit [config](config/config.R) file to be project and/or computer specific. See the [README](config/README.md) for more details.
4. Pre-process data for model fitting using the [data processing script](process_model_data.sh). This script takes 6 arguments, and can be run locally:
    1. the annotation type (e.g. `igor_alpha` for IGoR-annotated TRA sequences)
    2. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-and-model-type-options)
    3. the number of CPUs you want to use
    4. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    6. the desired model type; the options are described in the [model options section](#parameter-group-and-model-type-options)
4. Train model using the script suggested in the [model options section](#parameter-group-and-model-type-options) (e.g. the model training script will depend on the model type you are wanting to train). All of the modeling scripts take 7 arguments, and can be run locally or on a cluster: 
    1. the annotation type (e.g. `igor_alpha` for IGoR-annotated TRA sequences)
    2. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-options)
    3. the number of CPUs you want to use
    4. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    6. the desired model type; the options are described in the [model options section](#parameter-group-and-model-type-options)
    7. boolean value indicating whether you want to use L2 regularization for the MH-related terms in your model (e.g. choices are `True` or `False`)   

This trained model will be stored in the [trained_models](trained_models/) directory. 

### Parameter group and model type options

Here is a summary of some of the parameter group options:

| Parameter group argument             | Description                                                                                     | Model training script                              | Model type option (long name, same as manuscript)             | Trimming-related features parameterized                                                                                                                                                                                                                        | Ligation-related features parameterized                | Data processing model argument                            | Model training model argument                                     |
|--------------------------------------|-------------------------------------------------------------------------------------------------|----------------------------------------------------|---------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------|-----------------------------------------------------------|-------------------------------------------------------------------|
| `motif_two-side-base-count-beyond`   | Nonproductive sequences used to infer joint V-J trimming probabilities, independent of ligation | [fit_trim_model.sh](fit_trim_model.sh)             | null                                                          | None                                                                                                                                                                                                                                                           | NA                                                     | `null`                                                    | `null`                                                            |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond                            | 1x2 trimming motif and the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site                                                                                                      | NA                                                     | `motif_two-side-base-count-beyond`                        | `motif_two-side-base-count-beyond`                                |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond + interior-mh-count        | 1x2 trimming motif and the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site, as well as the number of non-contiguous MH nucleotides within overlapping regions between sequences | NA                                                     | `motif_two-side-base-count-beyond_interior-mh-count`      | `motif_two-side-base-count-beyond_interior-mh-count`              |
| `nonproductive_v-j_trim_ligation-mh` | Nonproductive sequences used to infer joint V-J trimming and ligation probabilities             | [fit_twostep_model_EM.sh](fit_twostep_model_EM.sh) | null                                                          | None                                                                                                                                                                                                                                                           | None                                                   | `null`                                                    | `twostep_null`                                                    |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond                            | 1x2 trimming motif and the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site                                                                                                      | None                                                   | `motif_two-side-base-count-beyond`                        | `twostep_motif_two-side-base-count-beyond`                        |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond + average-mh               | 1x2 trimming motif, the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site, and the average number of contiguous MH nucleotides across all possible ligation configurations        | None                                                   | `motif_two-side-base-count-beyond_average-mh`             | `twostep_motif_two-side-base-count-beyond_average-mh`             |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond + ligation-mh              | 1x2 trimming motif and the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site                                                                                                      | Number of MH nucleotides within ligation configuration | `motif_two-side-base-count-beyond_ligation-mh`            | `twostep_motif_two-side-base-count-beyond_ligation-mh`            |
|                                      |                                                                                                 |                                                    | motif + two-side-base-count-beyond + average-mh + ligation-mh | 1x2 trimming motif, the count of AT and GC nucleotides on either side of the trimming site, for both the V-gene trimming site and J-gene trimming site, and the average number of contiguous MH nucleotides across all possible ligation configurations        | Number of MH nucleotides within ligation configuration | `motif_two-side-base-count-beyond_average-mh_ligation-mh` | `twostep_motif_two-side-base-count-beyond_average-mh_ligation-mh` |



## Quantifying coefficient significance

If you would like to measure the significance of each model coefficient, you can follow the following steps:

1. Train the model of interest
2. Determine the relevant bootstrapping file which will depend on the parameter group and model type you are using. If you are using the [fit_trim_model.sh](fit_trim_model.sh) model training script, you will want to use the [bootstrap_trim_models.sh](bootstrap_trim_models.sh) script; Likewise, if you are using the [fit_twostep_model_EM.sh](fit_twostep_model_EM.sh) model training script, you will want to use the [bootstrap_twostep_models_EM.sh](bootstrap_twostep_models_EM.sh) script.
3. Run the bootstrapping script. These scripts take 7 arguments, and can be run on a cluster: 
    1. the annotation type (e.g. `igor_alpha` for IGoR-annotated TRA sequences)
    2. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-and-model-type-options)
    3. the number of CPUs you want to use
    4. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    6. the desired model type; the options are described in the [model options section](#parameter-group-and-model-type-options)
    7. boolean value indicating whether you want to use L2 regularization for the MH-related terms in your model (e.g. choices are `True` or `False`)   
4. Once the bootstrapped models have finished being trained, you can run the [significance testing script](scripts/run_coef_signif_test.R) to evaluate the significance of the model coefficients. This script takes 7 arguments and can be run locally:
    1. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-and-model-type-options) 
    2. the number of CPUs you want to use
    3. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    4. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. the model type of the model of interest; the options are described in the [model options section](#parameter-group-and-model-type-options)
    6. boolean value indicating whether you want to use L2 regularization for the MH-related terms in your model (e.g. choices are `True` or `False`)   
    7. the annotation type (e.g. `igor_alpha` for IGoR-annotated TRA sequences)

This analysis will save a file containing the results in a directory called `bootstraps` located in the indicated `OUTPUT_PATH` as specified in the [config](config) files


## Model comparison

If you would like to use a likelihood ratio test to compare two nested models, you can follow these steps:
1. Train the two models of interest using the same dataset (e.g. same input arguments to model training script, except for model type); Make sure the two models are nested!
2. Run the [LRT script](scripts/run_lrt.R) which takes 9 arguments and can be ran locally:
    1. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-and-model-type-options) 
    2. the number of CPUs you want to use
    3. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    4. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. the model type of the full model; the options are described in the [model options section](#parameter-group-and-model-type-options)
    6. the model type of the alternative model (which should contain a subset of the parameters that the full model contains)
    7. boolean value indicating whether you want to use L2 regularization for the MH-related terms in your model (e.g. choices are `True` or `False`)   
    8. the annotation type (e.g. `igor_alpha` for IGoR-annotated TRA sequences)
    9. the number of degrees of freedom (e.g. should be 1 if only one parameter is different between the two models being analyzed)

The resulting output file will be located in a directory within the `OUTPUT_PATH` as specified in the [config](config) file

## Model validation

If you would like to use the model to make predictions on a new data set and/or validate the model on a testing data set, you can follow these steps:
    
1. Download the processed testing data set, and make sure it contains the same column names as the training data set 
2. Depending on the locus of the testing data, download the appropriate gene name and germline sequences from IMGT (we downloaded these data in April 2023); save these data to a file--you will save the location of this file in the next step
3. Edit [config](config/config.R) file to be project and/or computer specific (be sure to add the path to the file from the last step)
2. Run the [model validation script](validate_model.sh) or the [twostep model validation script](validate_twostep_model_EM.sh), depending on the parameter group and model type you are using (e.g. see [parameter group and model type options](#parameter-group-and-model-type-options) section for more details). This script takes 8 arguments:
    1. the annotation type of the data used for model training (e.g. `igor_alpha` for IGoR-annotated TRA sequences)
    2. the parameter group which specifies sequence productivity and the prediction task; the options are described in the [parameter group options section](#parameter-group-and-model-type-options)
    3. the number of CPUs you want to use
    4. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    5. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    6. the desired model type; the options are described in the [model options section](#parameter-group-and-model-type-options)
    7. boolean value indicating whether you want to use L2 regularization for the MH-related terms in your model (e.g. choices are `True` or `False`)
    8. the annotation type of the data you want to use for validation (e.g. see the [model validation options summary](#model-validation-options) for details)

All output files will be located in a directory within the `OUTPUT_PATH` as specified in the [config](config) file

### Model validation options

Summary of the model validation data sets used in our analysis:

| Locus/Junction name  | Code argument           | Data used in the manuscript                                                                                                                                                                                                                                 |
|-------------|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _TRA locus/VJ junction_ | `validation_data_alpha` | [here](https://doi.org/10.1016/j.jaut.2021.102616) |
| _TRG locus/VJ junction_ | `validation_data_gamma` | [here](https://clients.adaptivebiotech.com/pub/TCRB-TCRG-comparison)  |


## Plot results

Plot [figures](plotting_scripts/manuscript_plots) from the manuscript.

Also, have a look at the plotting [README](plotting_scripts/README.md) for more details.



# About the analysis

With this analysis, we want to quantify the extent of microhomology involvement for nucleotide trimming and ligation during V(D)J recombination.
See the manuscript for specific model and methods details: 

TODO: add manuscript details

The following packages were especially helpful in our analyses:

Python:

* `jax`
* `jaxopt`
* `pandas`

R:

- `data.table` (Dowle and Srinivasan, 2021)
- `tidyverse` (Wickham et. al, 2019) 
- `doParallel` (Corporation and Steve Weston, 2020)
- `cowplot` (Wilke, 2020)
- `Biostrings` (Pagès et. al, 2021)
