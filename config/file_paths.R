python_output_path <- function(model_type=NULL){
    if (is.null(model_type)){
        model_type = MODEL_TYPE
    }
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT,'_', model_type))
    return(path)
}

get_model_coef_file_path <- function(L2, model_type=NULL, sample_annotation=SAMPLE_ANNOT){
    if (is.null(model_type)){
        model_type = MODEL_TYPE
    }
    path = python_output_path(model_type)
    if (sample_annotation==TRUE){
        file_name = paste0(path, '/trained_coefs_L2', L2, '.tsv')
    } else {
        file_name = paste0(path, '/trained_coefs_all_annotations_L2', L2, '.tsv')
    }
    return(file_name)
}

get_model_predictions_file_path <- function(L2, model_type=NULL, sample_annotation=SAMPLE_ANNOT){
    if (is.null(model_type)){
        model_type = MODEL_TYPE
    }
    path = python_output_path(model_type)
    if (sample_annotation==TRUE){
        file_name = paste0(path, '/predicted_dist_data_L2', L2, '.tsv')
    } else {
        file_name = paste0(path, '/predicted_dist_data_all_annotations_L2', L2, '.tsv')
    }
    return(file_name)
}

get_validation_predictions_file_path <- function(L2, validation_annotation, sample_annotation=SAMPLE_ANNOT){
    path = python_output_path(MODEL_TYPE)
    if (sample_annotation==TRUE){
        file_name = paste0(path, '/validation_prediction/', validation_annotation, '_predicted_dist_data_L2', L2, '.tsv')
    } else {
        file_name = paste0(path, '/validation_prediction/', validation_annotation, '_predicted_dist_data_all_annotations_L2', L2, '.tsv')
    }
    return(file_name)
}

get_model_eval_file_path <- function(L2, model_type=NULL, sample_annotation=SAMPLE_ANNOT){
    if (is.null(model_type)){
        model_type = MODEL_TYPE
    }
    path = python_output_path(model_type)
    if (sample_annotation==TRUE){
        file_name = paste0(path, '/model_evaluation_results_L2', L2, '.tsv')
    } else {
        file_name = paste0(path, '/model_evaluation_results_all_annotations_L2', L2, '.tsv')
    }
    return(file_name)
}

get_LRT_file_path <- function(L2, sample_annotation=SAMPLE_ANNOT){
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))

    if (sample_annotation==TRUE){
        file_name = paste0(path, '/model_lrt_results_L2', L2, '.tsv')
    } else {
        file_name = paste0(path, '/model_lrt_results_all_annotations_L2', L2, '.tsv')
    }
    return(file_name)
}

get_subsample_file_path <- function(L2, model_type=NULL){
    if (is.null(model_type)){
        model_type = MODEL_TYPE
    }
    path = python_output_path(model_type)
    path = file.path(path, 'subsampling_experiment')
    return(path)
}

get_manuscript_path <- function(){
    path = file.path(MOD_PROJECT_PATH, 'plots', ANNOTATION_TYPE, PARAM_GROUP, paste0(TRIM_TYPE, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'manuscript')
    dir.create(path, recursive = TRUE)
    return(path)
}
