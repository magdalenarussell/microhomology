extract_subject_ID <- function(tcr_repertoire_file_path){
    require(stringr)
    file_name = tail(str_split(tcr_repertoire_file_path, "/")[[1]], n = 1)
    localID = str_split(file_name, ".tsv")[[1]][1]
    stopifnot(!(is.na(localID)))
    return(localID)
}

convert_adaptive_gene_names_to_imgt <- function(data, gene_col){
    data = copy(data) 
    stopifnot((data[[gene_col]][1] %like% 'TCR') | (data[[gene_col]][1] == 'unknown')) 
    mapping = fread('https://raw.githubusercontent.com/kmayerb/tcrdist3/master/tcrdist/db/adaptive_imgt_mapping.csv')[species == 'human']
    gene_type = toupper(substring(gene_col, 1, 1))
    data[, c(gene_col, "allele") := tstrsplit(get(gene_col), "\\*0")]
    tog = merge(data, mapping, by.x = gene_col, by.y = 'adaptive', all.x = TRUE, allow.cartesian = TRUE)
    setnames(tog, 'imgt', paste0(tolower(gene_type), '_gene'))
    return(tog[, -c('allele', 'species')])
}

get_whole_v_nucseqs_and_anchors <- function(){
    require(stringr)
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', LOCUS, '_only')))[gene %like% 'V']

    whole_nucseq[, alignment := sapply(names, function(x) str_split(x, '\\|')[[1]][13])]
    whole_nucseq[, offset := sapply(alignment, function(x) str_split(x, '\\+')[[1]][2])]
    whole_nucseq[, offset := sapply(offset, function(x) str_split(x, '\\=')[[1]][1])]
    whole_nucseq[, anchor := 309 - as.numeric(offset)]
    whole_nucseq[, anchor_to_tail_sequence := substring(sequence, anchor + 1)]
    whole_nucseq[, begin_to_anchor_sequence := substring(sequence, 1, anchor)]
    return(whole_nucseq[, c('gene', 'sequence', 'anchor', 'anchor_to_tail_sequence', 'begin_to_anchor_sequence')])
}

find_subseq <- function(rearr, anchor_to_tail_seq){
    mismatch = 0
    start = NULL
    while (is.null(start)){
        match_index = as.data.frame(matchPattern(anchor_to_tail_seq, DNAString(rearr), max.mismatch = mismatch))
        if (nrow(match_index) == 0){
            start = NULL 
            mismatch = mismatch + 1
        } else if (nrow(match_index) > 1){
            match_index = match_index[1,]
            start = match_index$start 
        } else {
            start = match_index$start 
        }
    } 
    return(data.table('rearrangement' = as.character(rearr), 'start' = as.numeric(start), 'mismatch' = as.numeric(mismatch)))
}

stitch_V_seq_to_cdr3 <- function(data){
    require(Biostrings)
    require(parallel)

    cols = c('rearrangement', 'anchor_to_tail_sequence')
    subset = unique(data[, ..cols])
    match_indices = data.table(t(mcmapply(function(x, y) find_subseq(x, y), subset$rearrangement, subset$anchor_to_tail_sequence, mc.cores = NCPU)))
    match_indices = data.table(rearrangement = unlist(match_indices$rearrangement),
                               start = unlist(match_indices$start),
                               mismatch = unlist(match_indices$start))
    data = merge(data, match_indices, by = 'rearrangement')
    data[, old_rearrangement := rearrangement]
    data[, rearrangement := paste0(begin_to_anchor_sequence, substring(rearrangement, start))]
    return(data)
}

get_raw_cdr3_seqs <- function(tcr_repertoire_file, adaptive){
    data = fread(tcr_repertoire_file)
    if (isTRUE(adaptive)){
        data = convert_adaptive_gene_names_to_imgt(data, 'v_resolved')
        # stitch V-gene sequence to the beginning of adaptive cdr3
        whole_nucseq = get_whole_v_nucseqs_and_anchors()
        data = merge(data, whole_nucseq, by.x = 'v_gene', by.y = 'gene')
        data = stitch_V_seq_to_cdr3(data)
    }

    if (!('rearrangement' %in% colnames(data))){
        if ('cdr3_nucseq' %in% colnames(data)){
            setnames(data, 'cdr3_nucseq', 'rearrangement')
        }
        if ('sequence' %in% colnames(data)){
            setnames(data, 'sequence', 'rearrangement')
        }
        if ('CDR3nt' %in% colnames(data)){
            setnames(data, 'CDR3nt', 'rearrangement')
        }
    }
    stopifnot('rearrangement' %in% colnames(data))
    cdr3s = data$rearrangement
    if ('v_index' %in% colnames(data)){
        v_index = data$v_index
    } else {
        v_index = NULL
    }
    if ('count' %in% colnames(data)){
        cols = c('rearrangement', 'count')
        if (isTRUE(adaptive)){
            cols = c(cols, 'v_gene')
        }
        subset = data[, ..cols]
    } else if ('templates' %in% colnames(data)){
        cols = c('rearrangement', 'templates')
        if (isTRUE(adaptive)){
            cols = c(cols, 'v_gene')
        }
        subset = data[, ..cols]
    } else {
        subset = NULL
    }
    return(list(cdr3s = cdr3s, v_index = v_index, counts = subset))
}

get_raw_cdr3_seqs_first_subsample <- function(tcr_repertoire_file, sample_size, adaptive){
    data = fread(tcr_repertoire_file)
    if (isTRUE(adaptive)){
        data = data[frame_type == 'Out' | frame_type == 'Stop']
        data = data[v_resolved != 'unknown']
        inds = sample(seq(nrow(data)), sample_size, replace = TRUE, prob = data$templates)
    } else {
        data = data[productive == 'nonproductive']
        inds = sample(seq(nrow(data)), sample_size, replace = FALSE)
    }
    data = data[inds,]
    if (isTRUE(adaptive)){
        data = convert_adaptive_gene_names_to_imgt(data, 'v_resolved')
        # stitch V-gene sequence to the beginning of adaptive cdr3
        whole_nucseq = get_whole_v_nucseqs_and_anchors()
        data = merge(data, whole_nucseq, by.x = 'v_gene', by.y = 'gene')
        data = stitch_V_seq_to_cdr3(data)
    }

    if (!('rearrangement' %in% colnames(data))){
        if ('cdr3_nucseq' %in% colnames(data)){
            setnames(data, 'cdr3_nucseq', 'rearrangement')
        }
        if ('sequence' %in% colnames(data)){
            setnames(data, 'sequence', 'rearrangement')
        }
    }
    stopifnot('rearrangement' %in% colnames(data))
    cdr3s = data$rearrangement
    if ('v_index' %in% colnames(data)){
        v_index = data$v_index
    } else {
        v_index = NULL
    }
    if ('count' %in% colnames(data)){
        cols = c('rearrangement', 'count')
        if (isTRUE(adaptive)){
            cols = c(cols, 'v_gene')
        }
        subset = data[, ..cols]
    } else if ('templates' %in% colnames(data)){
        cols = c('rearrangement', 'templates')
        if (isTRUE(adaptive)){
            cols = c(cols, 'v_gene')
        }
        subset = data[, ..cols]
    } else {
        subset = NULL
    }
    return(list(cdr3s = cdr3s, v_index = v_index, counts = subset))
}



read_igor_output <- function(output_location){
    output = fread(file.path(output_location, 'foo_output/best_scenarios_counts.csv'))
}

sample_sequences <- function(output){
    set.seed(55)
    sampled = foreach(seq = unique(output$seq_index), .combine = rbind) %do% {
        subset = output[seq_index == seq]
        if (nrow(subset) > 0){
            subset[sample(nrow(subset), 1, prob = subset$scenario_proba_cond_seq)]
        }
    }
    return(sampled)
}

create_final_output_location <- function(output_directory, final_output_directory){
    require(stringr)
    subject = tail(str_split(output_directory, '/')[[1]], 1)
    name = paste0(subject, '_igor_sampled_annotations.tsv')
    dir.create(final_output_directory, recursive = TRUE, showWarnings = FALSE)
    return(file.path(final_output_directory, name))
}

get_gene_from_long_igor_column <- function(long_gene_column, sep){
    require(stringr)
    split = str_split(long_gene_column, sep)
    genes = sapply(split, function(x) x[x %like% 'TRA'])
    return(genes)
}
