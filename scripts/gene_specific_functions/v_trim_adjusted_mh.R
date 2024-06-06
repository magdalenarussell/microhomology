stopifnot(LOCUS == 'alpha')
stopifnot(INSERTIONS == 'zero')

get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[substring(gene, 4, 4) == 'J', sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    whole_nucseq[substring(gene, 4, 4) == 'J', j_gene_sequence := sequence]
    whole_nucseq[substring(gene, 4, 4) == 'V', v_gene_sequence := sequence]
    return(whole_nucseq[, -c('sequence')])
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/single.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))

filter_motif_data_for_possible_sites <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE, only_nonprod_sites = ONLY_NONPROD_SITES){
    return(motif_data)
}

get_all_possible_sites <- function(gene_type = GENE_NAME){
    return(NULL)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    unobserved = data.table()
    return(unobserved)
}

adjust_trimming_sites_for_ligation_mh <- function(tcr_dataframe, sample_annotation){
    mh_dt = get_possible_mh(tcr_dataframe, keep_gene_seqs = FALSE)
    mh_dt = count_mh_bordering_trim(mh_dt)
    mh_dt = reassign_trimming_sites_with_mh(mh_dt)
    mh_dt[, v_trim_adjusted_mh := adjusted_v_trim]
    mh_dt[, j_trim_adjusted_mh := adjusted_j_trim]
    return(mh_dt)
}

get_all_annotations <- function(tcr_dataframe, configs = read_frames_data()){
    return(tcr_dataframe)
}
