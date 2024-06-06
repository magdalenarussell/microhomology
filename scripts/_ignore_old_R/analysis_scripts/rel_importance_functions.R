#TODO need to generalize to both gene motifs/base counts
get_model_feature_scores <- function(validation_data, pwm, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    base_count_pwm = get_all_base_count_pwm_score(validation_data, pwm, trim_type)
    setnames(base_count_pwm, 'pwm_score', 'base_count_score')
    motif_pwm =  get_all_pwm_score(validation_data, pwm, NULL, trim_type = trim_type, gene_type = gene_type) 
    setnames(motif_pwm, 'pwm_score', 'motif_score')
    validation_data = merge(validation_data, base_count_pwm, by = c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC')))
    positions = get_positions(trim_type)
    validation_data = merge(validation_data, motif_pwm, by = c(paste0(gene_type, '_group'), trim_type, paste0(trim_type, '_motif'), positions))
    return(validation_data)
}

fit_rel_importance_model <- function(motif_data, formula = NULL, gene_type = GENE_NAME){
    if (is.null(formula)){
        formula = formula(paste0('cbind(weighted_observation, interaction(', gene_type, '_group ,subject)) ~ base_count_score + motif_score'))
    }
    model = mclogit(formula, 
                    data = motif_data)
    return(model)
}

get_per_run_rel_importance_model_path <- function(){
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'rel_importance', VALIDATION_TYPE, VALIDATION_TRIM_TYPE, VALIDATION_PRODUCTIVITY) 
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    return(path)
}

get_per_run_rel_importance_model_file_name <- function(){
    path = get_per_run_rel_importance_model_path()
    name = file.path(path, paste0('rel_importance_model_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

compile_rel_importance_result <- function(model_coefs){
    result = c(motif_length_5_end = LEFT_NUC_MOTIF_COUNT, motif_length_3_end = RIGHT_NUC_MOTIF_COUNT, motif_type = MOTIF_TYPE, gene_weight_type = GENE_WEIGHT_TYPE, upper_bound = UPPER_TRIM_BOUND, lower_bound = LOWER_TRIM_BOUND, model_type = MODEL_TYPE, terminal_melting_5_end_length = LEFT_SIDE_TERMINAL_MELT_LENGTH) 
    result = c(result, model_coefs)
    return(data.table(t(result)))
}
 
write_rel_importance_result_dt <- function(model_coefs, loss_results, file_name){
    result = compile_rel_importance_result(model_coefs)
    result$loss = loss_results$loss
    fwrite(result, file_name, sep = '\t')
    return(result)
}

compile_rel_importance_results <- function(){
    files = fs::dir_ls(path = get_per_run_rel_importance_model_path())
    require(parallel)
    cluster = makeCluster(NCPU)
    files_dt = parLapply(cluster, files, function(x){
                        data.table::fread(x)})
    stopCluster(cluster)
    rbound = rbindlist(files_dt, fill = TRUE)
    return(rbound)
}


