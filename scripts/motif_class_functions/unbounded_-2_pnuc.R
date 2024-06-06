PNUC_COUNT <<- -2

get_all_nuc_contexts <- function(tcr_dataframe, subject_id, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    motif_data = general_get_all_nuc_contexts(tcr_dataframe, subject_id, gene_type = gene_type, trim_type = trim_type)
    return(motif_data)
}
