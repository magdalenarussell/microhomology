calculate_subject_gene_weight <- function(compiled_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = FALSE, sample_annotation = SAMPLE_ANNOT){
    genes = get_gene_order(gene_type)
    trims = get_trim_vars(trim_type)

    if (sample_annotation == FALSE){
        compiled_data[, sub_index := 1:.N, by = index]
        subset = compiled_data[sub_index == 1]
        subset[, total_tcr := sum(count, na.rm = TRUE)]
        col = c(paste0(genes))
        subset[, paste0(gene_type, '_count') := sum(count, na.rm = TRUE), by = col]
        cols = c(col, 'total_tcr', paste0(gene_type, '_count'))
        subset_subset = unique(subset[, ..cols])
        compiled_data = merge(compiled_data, subset_subset)
    } else{
        compiled_data[, total_tcr := sum(count, na.rm = TRUE)]
        col = c(paste0(genes))
        compiled_data[, paste0(gene_type, '_count') := sum(count, na.rm = TRUE), by = col]
    }
    compiled_data[[paste0('p_', gene_type)]] = compiled_data[[paste0(gene_type, '_count')]]/compiled_data$total_tcr
    compiled_data[[paste0('p_', trim_type, '_given_', gene_type)]] = compiled_data$count/compiled_data[[paste0(gene_type, '_count')]]
    compiled_data$weighted_observation = compiled_data[[paste0('p_', trim_type, '_given_', gene_type)]]*compiled_data[[paste0('p_', gene_type)]]
    compiled_data$gene_weight_type = 'p_gene_pooled' 
    return(compiled_data)
}
