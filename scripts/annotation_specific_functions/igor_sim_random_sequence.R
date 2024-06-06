LOCUS <<- 'alpha'
DATA_TYPE <<- 'igor'
source(paste0(MOD_PROJECT_PATH, '/scripts/locus_specific_functions/', LOCUS, '.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/data_type_functions/', DATA_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/annotation_specific_functions/frame_data_functions/adjusted_frames.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    split_name = str_split(tcr_repertoire_file_path, "/")[[1]]

    file_name = split_name[length(split_name)]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    # localID = str_split(file_root_name, "_")[[1]][1]
    return(file_root_name)
}

generate_random_dna <- function(length_list) {
    bases <- c("A", "C", "G", "T")
    seqs = c() 
    for (length in length_list){
        random_sequence = sample(bases, size = length, replace = TRUE)
        seq = paste(random_sequence, collapse = "")
        seqs = c(seqs, seq)
    }
  return(seqs)
}

get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', LOCUS)))[, -c('name')]
    whole_nucseq[, length := nchar(sequence)]
    cols = c('gene', 'length')
    random_seqs = whole_nucseq[, ..cols] 
    random_seqs[, sequence := generate_random_dna(length)]    
    return(random_seqs[, c('gene', 'sequence')])
}
