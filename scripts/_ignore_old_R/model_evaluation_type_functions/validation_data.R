HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

evaluate_loss <- function(validation_data, model, trim_type = TRIM_TYPE, gene_type = GENE_NAME) {
    if (MODEL_TYPE == 'null'){
        parameter_count = 0
    } else {
        parameter_count = length(model$coefficients) 
    }
    # compute conditional logistic loss value for training data 
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', LOSS_GENE_WEIGHT, '.R'), local = TRUE)
    validation_data = calculate_subject_gene_weight(validation_data, gene_type = gene_type, trim_type = trim_type)
    log_loss = calculate_cond_expected_log_loss(model, validation_data)
    return(list(loss = log_loss, model_parameter_count = parameter_count, held_out_cluster_number = NA, held_out_genes = VALIDATION_TYPE))
}
