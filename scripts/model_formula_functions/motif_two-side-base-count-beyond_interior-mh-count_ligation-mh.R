stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)
stopifnot(TRIM_TYPE == 'v-j_trim_ligation-mh')

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/mh.R'))

get_parameter_vector <- function(trims, genes){
    trims = trims[trims %like% 'trim']
    motif_positions = c()
    for (i in seq(length(trims))){
        motif_positions = c(motif_positions, get_positions(trims[i]))
    }

    mh_vars = get_all_mh_count_variables(overlap_vector = c(1, 2, 3, 4), pos = c('mid'))

    left_vars = c()
    right_vars = c()
    for (i in seq(length(trims))){
        left_vars = c(left_vars, get_all_base_variables('5end', trims[i]))
        right_vars = c(right_vars, get_all_base_variables('3end', trims[i]))
    }

    return(c(motif_positions, left_vars, right_vars, mh_vars, 'ligation_mh'))
}

process_single_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE, whole_nucseq = whole_nucseq, gene_type = gene_type, trim_type = trim_type)
    together = process_for_mh(together, whole_nucseq = whole_nucseq, overlap_vector = c(1, 2, 3, 4), trim_type = TRIM_TYPE, gene_type = GENE_NAME, prop = FALSE, positions = c('mid'))
    return(together)
}
