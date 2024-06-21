stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)
stopifnot(TRIM_TYPE == 'trim_ligation-mh')

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/mh.R'))

get_parameter_vector <- function(trims, genes){
    return(c('mh_config_count', 'ligation_mh'))
}

process_single_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    together = get_mh_config_count(together, trim_type = trim_type, gene_type = gene_type)
    return(together)
}
