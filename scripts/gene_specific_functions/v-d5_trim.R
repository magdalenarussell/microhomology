stopifnot(LOCUS %like% 'beta')

source(paste0(MOD_PROJECT_PATH, '/scripts/gene_specific_functions/junction_specific_functions/vd.R'))

filter_motif_data_for_possible_sites <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES){
    return(motif_data)
}

adjust_trimming_sites_for_ligation_mh <- function(motif_data, sample_annotation, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    return(motif_data)
}

get_all_possible_sites <- function(gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    return(NULL)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    unobserved = data.table()
    return(unobserved)
}

get_all_annotations <- function(tcr_dataframe, configs = read_frames_data(), trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    return(tcr_dataframe)
}
