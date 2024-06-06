HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

evaluate_loss <- function(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME) {
    # fit model to the training data 
    if (MODEL_TYPE != 'null') {
        model = fit_model(motif_data, trim_type= trim_type)
        parameter_count = length(model$coefficients) 
    } else {
        model = 'null'
        parameter_count = 0
    } 
    # adjust weighting to reflect actual p_gene_given_subject weight for log loss calculation
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/p_gene_pooled.R'), local = TRUE)
    motif_data = calculate_subject_gene_weight(motif_data,gene_type = gene_type, trim_type = trim_type)
    # compute conditional logistic loss value for training data 
    log_loss = calculate_cond_expected_log_loss(model, motif_data)
    return(list(loss = log_loss, model_parameter_count = parameter_count, held_out_cluster_number = NA, held_out_genes = NA))
}
