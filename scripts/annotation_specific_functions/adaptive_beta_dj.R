LOCUS <<- 'beta_dj'
DATA_TYPE <<- 'adaptive'
CONVERT_FREQ <<- FALSE
source(paste0(MOD_PROJECT_PATH, '/scripts/locus_specific_functions/', LOCUS, '.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/data_type_functions/', DATA_TYPE, '.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    split_name = str_split(tcr_repertoire_file_path, "/")[[1]]

    file_name = split_name[length(split_name)]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', LOCUS)))[, -c('name')]
    return(whole_nucseq[, c('gene', 'sequence')])
}
