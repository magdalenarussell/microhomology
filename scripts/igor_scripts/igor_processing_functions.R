extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = tail(str_split(tcr_repertoire_file_path, "/")[[1]], n = 1)
    localID = str_split(file_name, ".tsv")[[1]][1]
    stopifnot(!(is.na(localID)))
    return(localID)
}

get_raw_cdr3_seqs <- function(tcr_repertoire_file){
    data = fread(tcr_repertoire_file)
    if (!('rearrangement' %in% colnames(data))){
        if ('cdr3_nucseq' %in% colnames(data)){
            setnames(data, 'cdr3_nucseq', 'rearrangement')
        }
    }
    stopifnot('rearrangement' %in% colnames(data))
    cdr3s = data$rearrangement
    v_index = data$v_index
    return(list(cdr3s = cdr3s, v_index = v_index))
}

read_igor_output <- function(output_location){
    output = fread(file.path(output_location, 'foo_output/best_scenarios_counts.csv'))
}

sample_sequences <- function(output){
    set.seed(55)
    sampled = foreach(seq = unique(output$seq_index), .combine = rbind) %do% {
        subset = output[seq_index == seq]
        subset[sample(nrow(subset), 1, prob = subset$scenario_proba_cond_seq)]
    }
    return(sampled)
}

create_final_output_location <- function(output_directory, final_output_directory){
    subject = tail(str_split(output_directory, '/')[[1]], 1)
    name = paste0(subject, '_igor_sampled_annotations.tsv')
    dir.create(final_output_directory, recursive = TRUE, showWarnings = FALSE)
    return(file.path(final_output_directory, name))
}

get_gene_from_long_igor_column <- function(long_gene_column, sep){
    split = str_split(long_gene_column, sep)
    genes = sapply(split, function(x) x[x %like% 'TRA'])
    return(genes)
}
