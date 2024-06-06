stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)
stopifnot(TRIM_TYPE %like% 'v-j_trim')

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/mh.R'))

get_parameter_vector <- function(trims, genes){
    trims = trims[trims %like% 'trim']
    motif_positions = c()
    for (i in seq(length(trims))){
        motif_positions = c(motif_positions, get_positions(trims[i]))
    }

    mh_vars = get_all_mh_prop_variables(overlap_vector = c(0, 1, 2, 3, 4))

    left_vars = c()
    right_vars = c()
    for (i in seq(length(trims))){
        left_vars = c(left_vars, paste0(get_all_base_variables('5end', trims[i])))
        right_vars = c(right_vars, paste0(get_all_base_variables('3end', trims[i]), '_prop'))
    }

    return(c(motif_positions, left_vars, right_vars, mh_vars, 'v_length', 'j_length'))
}

process_single_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE, whole_nucseq = whole_nucseq, gene_type = gene_type, trim_type = trim_type)
    together = process_for_mh(together, whole_nucseq = whole_nucseq, overlap_vector = c(0, 1, 2, 3, 4), trim_type = trim_type, gene_type = gene_type)
    base_count_cols = colnames(together)[colnames(together) %like% '3end_base_count']
    base_count_cols = base_count_cols[base_count_cols %like% trim_type]
    base_count_col_types = unique(substring(base_count_cols, 1, 22))

    for (col in base_count_col_types){
        gc_temp = paste0(col, '_GC')
        at_temp = paste0(col, '_AT')
        gc_prop_temp = paste0(col, '_GC_prop')
        at_prop_temp = paste0(col, '_AT_prop')
        together[((get(gc_temp) + get(at_temp)) > 0), paste(gc_prop_temp) := get(gc_temp)/(get(gc_temp) + get(at_temp))]
        together[((get(gc_temp) + get(at_temp)) > 0), paste(at_prop_temp) := get(at_temp)/(get(gc_temp) + get(at_temp))]
        together[((get(gc_temp) + get(at_temp)) == 0), paste(at_prop_temp) := 0]
        together[((get(gc_temp) + get(at_temp)) == 0), paste(gc_prop_temp) := 0]

        cols = colnames(together)[!(colnames(together) %in% c(gc_temp, at_temp))]
        together = together[, ..cols]
    }
    new_name = paste0(substring(gene_type, 1, 1), '_length')
    together[, paste(new_name) := get(trim_type) + 2]
    return(together)
}
