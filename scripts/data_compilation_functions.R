source(paste0(MOD_PROJECT_PATH, '/scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/annotation_specific_functions/', ANNOTATION_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/gene_specific_functions/', TRIM_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_formula_functions/', MODEL_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/mh_functions.R'))

REQUIRED_COMMON_NUCS_5 <<- UPPER_TRIM_BOUND + 10 

get_gene_order <- function(gene_type){
    type = strsplit(gene_type, '_')[[1]][1]
    if (type %like% '-'){
        type = strsplit(type, '-')[[1]]
    }
    final = paste0(type, '_gene')
    return(final)
}

get_trim_order <- function(trim_type){
    type = strsplit(trim_type, '_')[[1]][1]
    if (type %like% '-'){
        type = strsplit(type, '-')[[1]]
    }
    final = paste0(type, '_trim')
    if (trim_type %like% 'adjusted_mh'){
        final = paste0(final, '_adjusted_mh')
    }
    return(final)
}

get_trim_vars <- function(trim_type){
    trim_sites = get_trim_order(trim_type)
    if (trim_type %like% 'ligation-mh'){
        trim_sites = c(trim_sites, 'ligation_mh')
    }
    return(trim_sites)
}

get_subject_motif_output_location <- function(sample_annotation=SAMPLE_ANNOT){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'motif_data')
    if (isFALSE(sample_annotation)){
        output_location = file.path(output_location, 'all_annotations')
    }
    return(output_location)
}

get_nuc_context <- function(whole_gene_seqs, trim_lengths){
     # Ensure input vectors have the same length
     stopifnot(length(whole_gene_seqs) == length(trim_lengths))

     # Initialize lists for storing nucleotide contexts
     left_nuc_list = c()
     right_nuc_list = c()
     for (index in seq(1, length(whole_gene_seqs))){
        # Convert gene sequence to DNAString
        whole_gene_seq = DNAString(whole_gene_seqs[index])

        # Handle positive and negative PNUC_COUNT values
        if (PNUC_COUNT > 0){
            possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, PNUC_COUNT)
        } else if (PNUC_COUNT < 0){
            possible_pnucs_5_to_3 = DNAString() 
            whole_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq) + PNUC_COUNT)             
        } else {
            possible_pnucs_5_to_3 = DNAString() 
        }

        # Concatenate gene sequence and possible PNUCs
        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        # Trim gene sequence based on specified length
        trim_length = trim_lengths[index]
        trimmed_gene_seq = substring(whole_gene_and_pnucs, 1, nchar(whole_gene_seq)-trim_length)
        trimmed_length = nchar(trimmed_gene_seq)

        # Extract left nucleotide context
        left_nucs = substring(trimmed_gene_seq, trimmed_length - (REQUIRED_COMMON_NUCS_5 - 1), trimmed_length) 

        # Extract right nucleotide context
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + PNUC_COUNT) + 1, nchar(whole_gene_and_pnucs))

        # Adjust right nucleotide context based on REQUIRED_COMMON_NUCS_5
        if (nchar(seq_right_of_trim) < REQUIRED_COMMON_NUCS_5){
            missing_nucs = DNAString(strrep('-', REQUIRED_COMMON_NUCS_5 - nchar(seq_right_of_trim)))
            seq_right_of_trim = c(unlist(seq_right_of_trim), unlist(missing_nucs))
        } else if (nchar(seq_right_of_trim) > REQUIRED_COMMON_NUCS_5){
            seq_right_of_trim = substring(seq_right_of_trim, 1, REQUIRED_COMMON_NUCS_5)
        }

        # Add nucleotide contexts to their respective lists
        left_nuc_list = c(left_nuc_list, as.character(left_nucs))
        right_nuc_list = c(right_nuc_list, as.character(seq_right_of_trim))
    }
    return(list(left_nucs = left_nuc_list, right_nucs = right_nuc_list))
}

apply_get_nuc_context <- function(df, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    # Calculate nucleotide context
    for (i in seq(length(trims))){
        u_cols = c(paste0(genes[i], '_whole_seq'), paste0(genes[i]), trims[i])
        subset = unique(df[, ..u_cols])
        subset = subset[, c(paste0(trims[i], '_left_nucs'), paste0(trims[i], '_right_nucs')):= get_nuc_context(get(paste0(genes[i], '_whole_seq')), get(trims[i]))] 
        df = merge(df, subset, by = u_cols)
    }
    return(df)
}


get_unobserved_nuc_context <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    # Retrieve trim and gene orders based on type
    trims = get_trim_order(trim_type)
    trim_vars = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    # Filter out specified columns to get observed data
    remove_cols = c('count', paste0(genes, 'whole_seq'), paste0(trims, '_left_nucs'), paste0(trims, '_right_nucs'), 'index')
    cols = colnames(tcr_dataframe)[!(colnames(tcr_dataframe) %in% remove_cols)]
    tcr_dataframe_observed = tcr_dataframe[,.N, by = cols]
    
    # Extract unique observations excluding trim variables
    non_trim_cols = cols[!(cols %in% c(trim_vars, 'N'))]
    unique_obs = unique(tcr_dataframe_observed[,..non_trim_cols])

    # Generate desired observations with varying trim lengths
    trim_lengths = seq(LOWER_TRIM_BOUND,UPPER_TRIM_BOUND)

    desired_obs = unique_obs
    for (i in seq(length(trims))){
        desired_obs = desired_obs %>%
            mutate(!!trims[i]:= list(trim_lengths)) %>%
            unnest(cols = trims[i]) %>%
            as.data.table()
    }

    # Determine unobserved data by comparison
    unobserved_cols = cols[cols != 'ligation_mh']
    unobserved = desired_obs[!tcr_dataframe_observed, on = unobserved_cols]

    unobserved = inner_unobserved_nucleotide_context_function(unobserved, trim_type, gene_type)
    return(unobserved)
}

inner_unobserved_nucleotide_context_function <- function(unobserved, trim_type=TRIM_TYPE, gene_type=GENE_NAME){
    # Retrieve trim and gene orders based on type
    trims = get_trim_order(trim_type)
    trim_vars = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    # Calculate unobserved nucleotide context
    unobserved = apply_get_nuc_context(unobserved, trim_type, gene_type)

    # Initialize count and ligation mismatch (if applicable) to zero
    unobserved$count = 0
    if (!('ligation_mh' %in% colnames(unobserved)) & ('ligation_mh' %in% trim_vars)){
        unobserved$ligation_mh = 0
    }
    return(unobserved)
}

adaptive_data_filtering <- function(data){
    if ('adaptive_v_gene_call' %in% colnames(data)){
        data = data[v_gene == adaptive_v_gene_call]
        if (LOCUS == 'alpha'){
            data = data[get(JOINING_INSERT) <= 15]
        }
    }
    return(data)
}

filter_by_productivity <- function(data){
    # only filter when we are looking at all sequences, regardless of insertion status
    # when we are looking at ligation MH (e.g. with zero insertion sequences)
    # within VJ junctions, we will do productivity filtering later
    stopifnot(PRODUCTIVITY %in% c('productive', 'nonproductive', 'both'))
    if (PRODUCTIVITY == 'productive' & !(grepl('ligation', PARAM_GROUP))){
        if (all(unique(data$productive) %in% c(TRUE, FALSE))){
            data = data[productive == TRUE]
        } else {
            data = data[productive == 'productive']
        }
    } else if (PRODUCTIVITY == 'nonproductive' & !(grepl('ligation', PARAM_GROUP))){
        if (all(unique(data$productive) %in% c(TRUE, FALSE))){
            data = data[productive == FALSE]
        } else {
            data = data[productive == 'nonproductive']
        }
    }
    return(data)
}

filter_by_insertions <- function(tcr_dataframe, insertions = INSERTIONS){
    # Filter the dataframe based on the 'insertions' criteria
    if (insertions == 'nonzero'){
        tcr_dataframe = tcr_dataframe[get(JOINING_INSERT) != 0]
    } else if (insertions == 'zero'){
        tcr_dataframe = tcr_dataframe[get(JOINING_INSERT) == 0]
    }
    return(tcr_dataframe)
}

condense_tcr_data <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE, insertions = INSERTIONS){
    stopifnot(insertions %in% c('nonzero', 'zero', 'all'))

    # Retrieve trim variables and gene orders based on input types
    trim_vars = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    # Define columns to be retained from the input dataframe and filter
    cols = c(paste0(genes, '_sequence'), paste0(genes), paste(trim_vars), JOINING_INSERT)
    if ('count' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'count')
    }
    if ('index' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'index', 'processed_sequence')
    }

    tcr_dataframe = tcr_dataframe[,..cols]

    # Filter the dataframe based on the 'insertions' criteria
    tcr_dataframe = filter_by_insertions(tcr_dataframe, insertions)

    # Rename gene sequence columns for clarity
    setnames(tcr_dataframe, paste0(genes, '_sequence'), paste0(genes, '_whole_seq'))

    # Condense data by gene, trim, and other factors
    cols = c(paste0(genes, '_whole_seq'), trim_vars, paste0(genes))
    if ('index' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'index', 'processed_sequence')
    }
    if ('count' %in% colnames(tcr_dataframe)){
        # Summarize data by count if 'count' column exists
        condensed_tcr = tcr_dataframe[, sum(count), by = cols]
        setnames(condensed_tcr, 'V1', 'count')
    } else {
        # Otherwise, use the number of rows (.N) as a count
        condensed_tcr = tcr_dataframe[, .N, by = cols]
        setnames(condensed_tcr, 'N', 'count')
    }
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    trim_vars = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    cols = c(paste0(genes), trim_vars, paste0(trims, '_left_nucs'), paste0(trims, '_right_nucs'))
    if ('index' %in% colnames(condensed_tcr_dataframe)){
        cols = c(cols, 'index', 'processed_sequence')
    }
    if ('frame_type' %in% colnames(condensed_tcr_dataframe)){
        cols = c(cols, 'frame_type', 'frame_stop')
    }
    summed = condensed_tcr_dataframe[, sum(count, na.rm = TRUE), by = cols]
    setnames(summed, 'V1', 'count')
    return(summed)
}

general_get_all_nuc_contexts <- function(tcr_dataframe, subject_id, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    # Retrieve trim and gene orders based on input types
    trims = get_trim_order(trim_type)
    trim_vars = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    # Filter data based on trim bounds
    if (length(trims) > 1){
        tcr_dataframe = tcr_dataframe[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND & get(trims[2]) >= LOWER_TRIM_BOUND & get(trims[2]) <= UPPER_TRIM_BOUND]
    } else {
         tcr_dataframe = tcr_dataframe[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND]
    }

    # Return empty dataframe if no rows match criteria
    if (nrow(tcr_dataframe) == 0){
        return(tcr_dataframe)
    } else {
        # Condense data by gene and trim
        tcr_dataframe = condense_tcr_data(tcr_dataframe, gene_type = gene_type, trim_type = trim_type)

        # Extract motifs for observed trims
        motif_dataframe = tcr_dataframe
        motif_dataframe = apply_get_nuc_context(motif_dataframe, trim_type)

        # Include motifs for unobserved sequences
        if (length(genes) == 1){
            desired_count = length(unique(motif_dataframe[[paste0(genes)]]))*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)
        } else if (length(genes) == 2){
            desired_count = length(unique(interaction(motif_dataframe[[paste0(genes[1])]], motif_dataframe[[paste0(genes[2])]])))*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)
        }
 
        if ((nrow(motif_dataframe) < desired_count) & !grepl('ligation-mh', PARAM_GROUP)){
            unobserved = get_unobserved_nuc_context(motif_dataframe, gene_type = gene_type, trim_type = trim_type)
            together = rbind(motif_dataframe, unobserved)
        } else {
            together = motif_dataframe
        }

        # Re-condense observations by group
        recondensed = sum_trim_observations(together, gene_type = gene_type, trim_type = trim_type)

        # Mark observations as observed or unobserved
        recondensed[count == 0, paste0(trim_type, '_observed') := FALSE]
        recondensed[count != 0, paste0(trim_type, '_observed') := TRUE]

        # Add gene type and subject information
        recondensed$gene_type = gene_type 
        if (!is.null(subject_id)){
            recondensed$subject = subject_id
        }
        return(recondensed)
    }
}

convert_data_to_motifs <- function(compiled_data, left_window_size = LEFT_NUC_MOTIF_COUNT, right_window_size = RIGHT_NUC_MOTIF_COUNT, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    for (i in seq(length(trims))){
        left =paste0(trims[i], '_left_nucs') 
        right = paste0(trims[i], '_right_nucs')
        compiled_data[, paste0(trims[i], '_left_motif') := substring(get(left), REQUIRED_COMMON_NUCS_5-(left_window_size-1), REQUIRED_COMMON_NUCS_5)]
        compiled_data[, paste0(trims[i], '_right_motif') := substring(get(right), 1, right_window_size)]
        compiled_data[, paste0(trims[i], '_motif') := paste0(get(left), get(right))] 
    }
    return(compiled_data)
}

process_data_for_model_fit <- function(group_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    processed = group_motif_data
    for (i in seq(length(trims))){
        processed = process_single_data_for_model_fit(processed, gene_type = genes[i], trim_type = trims[i])
    }
    return(processed)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    for (i in seq(length(trims))){
        positions = get_positions(trims[i])
        for (side in c('5end', '3end')){
            subset = positions[positions %like% side]
            lr_side = ifelse(side == '5end', 'left', 'right')
            if (length(subset) > 0){
                aggregated_subject_data[, paste0(subset) := tstrsplit(get(paste0(trims[i], '_', lr_side, '_motif')), "")]
            } 
        }
    }
    return(aggregated_subject_data)
}

compile_data_for_subject <- function(file_path=NULL, dataset=NULL, write = TRUE, gene_type = GENE_NAME, trim_type = TRIM_TYPE, sample_annotation=SAMPLE_ANNOT, only_nonprod_sites = ONLY_NONPROD_SITES){
    # Validate that either file_path or dataset is provided, but not both
    stopifnot(!(is.null(file_path)) | !(is.null(dataset)))
    stopifnot(!(!(is.null(file_path)) & !(is.null(dataset))))

    # Load data from file or dataset
    if (!(is.null(file_path))){
        temp_data = fread(file_path)
        subject_id = extract_subject_ID(file_path)
    } else if (!(is.null(dataset))){
        temp_data = dataset
        subject_id = unique(temp_data$sample_name)
    }

    # Reformat and filter data
    temp_data = reformat_data(temp_data)
    if (gene_type == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }

    # Filter data by productivity and other adaptive data factors
    temp_data = filter_by_productivity(temp_data)    
    temp_data = adaptive_data_filtering(temp_data)

    # Prepare output location
    output_location = get_subject_motif_output_location(sample_annotation) 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)

    # filter by insertions
    temp_data = filter_by_insertions(temp_data)

    # Adjust trimming sites if necessary based on MH ligation
    TRIMMING_LIGATION_REANNOTATED <<- FALSE
    temp_data = adjust_trimming_sites_for_ligation_mh(temp_data, sample_annotation, trim_type, gene_type)

    # Get oriented full sequences and group genes by common features
    together = get_oriented_full_sequences(temp_data, gene_type = gene_type)

    # Get NT context for each scenario
    if ('subject' %in% colnames(together)){
        motif_data = data.table()
        for (subject_i in unique(together$subject)){
            temp = get_all_nuc_contexts(together[subject == subject_i], subject_id = subject_i, gene_type = gene_type, trim_type = trim_type)
            motif_data = rbind(motif_data, temp)
        }
    } else {
        motif_data = get_all_nuc_contexts(together, subject_id, gene_type = gene_type, trim_type = trim_type)
    }

    # Filter motif data for possible sites
    motif_data = filter_motif_data_for_possible_sites(motif_data, gene_type = gene_type, trim_type = trim_type, only_nonprod_sites = only_nonprod_sites)

    # Write or return motif data
    if (isTRUE(write)){
        fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
    } else {
        return(motif_data)
    }
}

compile_all_data <- function(directory, gene_type = GENE_NAME, trim_type = TRIM_TYPE, sample_annotation = SAMPLE_ANNOT, only_nonprod_sites = ONLY_NONPROD_SITES){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        compile_data_for_subject(file, gene_type = gene_type, trim_type = trim_type, sample_annotation=sample_annotation, only_nonprod_sites = only_nonprod_sites)
        print(paste0(file))
    }
    stopImplicitCluster()
}

frame_data_path <- function(){
    output_location = file.path(MOD_OUTPUT_PATH, 'meta_data', CHAIN_TYPE, SUB_JUNCTION_TYPE)
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    filename = file.path(output_location, 'frame_data.tsv')
    return(filename)
} 

read_frames_data <- function(){
    file_name = frame_data_path() 
    if (!file.exists(file_name)){
        frame_data = get_frames_data()
        fwrite(frame_data, file_name, sep = '\t')
    } else {
        frame_data = fread(file_name)
    }
    return(frame_data)
}

get_cdr3_indices <- function(frame_data){
    v_type = paste0(CHAIN_TYPE, 'V')
    j_type = paste0(CHAIN_TYPE, 'J')

    vframe = frame_data[id %like% v_type]
    vframe[, c("extra", "cdr1", "cdr2", 'cdr3') := tstrsplit(cdr_columns, ";", fixed = TRUE)]
    vframe[, c("cdr3_start", "CDR3_end") := tstrsplit(cdr3, "-", fixed = TRUE)]
    vframe[, seq_start := substring(aligned_protseq, 1, cdr3_start)]
    vframe[, dot_count := stringr::str_count(seq_start, "\\.")]
    vframe[, cys_index := as.numeric(cdr3_start) - dot_count]
    vframe[, cys_character := substring(seq_start, nchar(seq_start))]
    vframe[cys_character == '.', cys_character := ""]
    setnames(vframe, 'id', 'v_gene')

    jframe = frame_data[id %like% j_type]
    jframe[, c("cdr3_start", "CDR3_end") := tstrsplit(cdr_columns, "-", fixed = TRUE)]
    jframe[, seq_start := substring(aligned_protseq, 1, CDR3_end)]
    jframe[, seq_end := substring(aligned_protseq, CDR3_end)]
    jframe[, dot_count := stringr::str_count(seq_start, "\\.")]
    jframe[, dot_count_end := stringr::str_count(seq_end, "\\.")]
    jframe[, phe_index := as.numeric(CDR3_end) - dot_count]
    jframe[, phe_from_end_index := nchar(seq_end) - 1 - dot_count_end]
    jframe[, phe_character := substring(seq_end, 1, 1)]
    setnames(jframe, 'id', 'j_gene')

    return(list(v = vframe[, c('v_gene', 'cys_index', 'cys_character')], j = jframe[, c('j_gene', 'phe_index', 'phe_from_end_index', 'phe_character')]))
}

get_cdr3_frame <- function(adjusted_data){
    # check for matching cysteine and phenylalanine/tryptophan
    adjusted_data[(substring(processed_cdr3, nchar(processed_cdr3)) == phe_character), frame_type := 'In']
    adjusted_data[!(substring(processed_cdr3, nchar(processed_cdr3)) == phe_character), frame_type := 'Out']
    return(adjusted_data)
}

get_stop_codons <- function(adjusted_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME, keep_extra_cols = c()){
    trims = get_trim_vars(trim_type)
    genes = get_gene_order(gene_type)

    require(Biostrings)
    
    stops = grepl("\\*", adjusted_data$processed_cdr3)
    adjusted_data[, frame_stop := stops]

    # Select relevant columns for output
    cols2 = c(genes, trims, 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'frame_stop', keep_extra_cols)
    adjusted_data = unique(adjusted_data[, ..cols2])
    return(adjusted_data)
}

processed_data_path <- function(sample_annotation=SAMPLE_ANNOT){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_', MODEL_TYPE))
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    if (sample_annotation){
        name = 'processed_data.tsv' 
    } else {
        name = 'processed_data_all_annotations.tsv'
    }
    filename = file.path(output_location, name)
    return(filename)
}

bootstrap_data_path <- function(iteration, sample_annotation=SAMPLE_ANNOT){
    filename = processed_data_path(sample_annotation)
    bn = basename(filename)
    bn = str_replace(bn, '.tsv', paste0('_bootstrap_', iteration, '.tsv'))
    path = file.path(dirname(filename), 'bootstraps')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    boot_name = file.path(path, bn)
    return(boot_name)
}


subsampling_processed_data_path <- function(prop, iter){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_', MODEL_TYPE), 'temp_subsampling_exp', paste0('prop', prop))
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    filename = file.path(output_location, paste0('processed_data_', iter, '.tsv'))
    return(filename)
}

subset_processed_data <- function(data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trim_vars = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    params = get_parameter_vector(trim_vars, genes)
    other = c(paste0(genes), trim_vars, 'weighted_observation', 'count', 'total_tcr')
    cols = c(other, params)
    if ('index' %in% colnames(data)){
        cols = c('index', cols)
    }
    if ('frame_type' %in% colnames(data)){
        cols = c(cols, 'frame_type', 'frame_stop')
    }
    return(data[, ..cols])
}

get_positions <- function(trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    if (LEFT_NUC_MOTIF_COUNT > 0){
        left = c()
        for (i in seq(length(trims))){
            left = c(left, paste0(trims[i], '_motif_5end_pos', seq(LEFT_NUC_MOTIF_COUNT, 1)))
        }
    } else {
        left = c()
    } 
    if (RIGHT_NUC_MOTIF_COUNT > 0){
        right = c()
        for (i in seq(length(trims))){
            right = c(right, paste0(trims[i], '_motif_3end_pos', seq(1, RIGHT_NUC_MOTIF_COUNT)))
        }
    } else {
        right = c()
    }
    positions = c(left, right)
    return(positions)
}

aggregate_all_subject_data <- function(directory = get_subject_motif_output_location(sample_annotation), gene_type = GENE_NAME, trim_type = TRIM_TYPE, sample_annotation = SAMPLE_ANNOT, only_nonprod_sites = ONLY_NONPROD_SITES){
    stopifnot(LEFT_NUC_MOTIF_COUNT <= 10)
    stopifnot(LEFT_SIDE_TERMINAL_MELT_LENGTH <= 10 | is.na(LEFT_SIDE_TERMINAL_MELT_LENGTH))
    desired_file_count = length(list.files(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE))))
    if (!dir.exists(directory) | !(length(list.files(directory)) == desired_file_count)) {
        print('compiling motif data, first')
        compile_all_data(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE)), gene_type = gene_type, trim_type = trim_type, sample_annotation = sample_annotation) 
    }  
    
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        file_data
    }
    
    weighted_together = inner_aggregation_processing(together, gene_type, trim_type, only_nonprod_sites, sample_annotation)
    stopImplicitCluster()
    return(weighted_together)
}

inner_aggregation_processing <- function(together, gene_type, trim_type, only_nonprod_sites = ONLY_NONPROD_SITES, sample_annotation = SAMPLE_ANNOT){
    genes = get_gene_order(gene_type)
    if ('index' %in% colnames(together)){
        cols = c(genes, 'processed_sequence')
        # recondense across all individuals
        together[, index := .GRP, by = cols]
    } 
    together = sum_trim_observations(together, gene_type, trim_type) 
    if (only_nonprod_sites == TRUE & grepl('ligation-mh', PARAM_GROUP)){
        if (grepl('nonprod', PARAM_GROUP)){
            together[frame_type == 'In' & frame_stop == FALSE, count := NA]
        } else if (grepl('prod', PARAM_GROUP)){
            together[frame_type == 'Out' | frame_stop == TRUE, count := NA]
        }
    }

    if (MODEL_TYPE %like% 'dna_shape') {
        together = convert_data_to_motifs(together, left_window_size = LEFT_NUC_MOTIF_COUNT + 2, right_window_size = RIGHT_NUC_MOTIF_COUNT + 2)
        processed_motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
        motif_data = convert_data_to_motifs(processed_motif_data)
    } else {
        together = convert_data_to_motifs(together)
        motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
    }

    cols = colnames(motif_data)[!(colnames(motif_data) %like% 'left_nucs')]
    cols = cols[!(cols %like% 'right_nucs')]

    motif_data = motif_data[, ..cols]
    together_pos = split_motif_column_by_motif_position(motif_data, trim_type = trim_type) 
    weighted_together = calculate_subject_gene_weight(together_pos, gene_type = gene_type, trim_type = trim_type, only_nonprod_sites, sample_annotation)
    return(weighted_together)
}

get_coef_pvalues <- function(bootstrap_results, original_model_results){
    sd_coeff = bootstrap_results[, sd(coefficient), by = .(parameter, base)]
    setnames(sd_coeff, 'V1', 'sd')

    together = merge(original_model_results, sd_coeff)

    together[, zstat := coefficient/sd]
    together[, pvalue := 2*pnorm(-abs(zstat))]
    together$iterations = max(bootstrap_results$iteration)
    return(together)
}

subsample <- function(motif_data, prop, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    # sample size is the number of gene combos
    if (length(genes) == 1){
        motif_data$cluster = motif_data[[paste0(gene_type)]]
    } else if (length(genes) == 2){
        motif_data$cluster = interaction(motif_data[[paste0(genes[1])]], motif_data[[paste0(genes[2])]])
    }

    # sample proportion of sequences for each individual
    size = ceiling(unique(motif_data$total_tcr)*prop)
    motif_data[, subsample_total_tcr := size]
    vars = c(colnames(motif_data)[colnames(motif_data) %like% 'mh_'],
             colnames(motif_data)[colnames(motif_data) %like% 'base_count'],
             colnames(motif_data)[colnames(motif_data) %like% 'motif'])

    cols = c(paste0(genes), trims, vars, 'count', 'subsample_total_tcr', 'cluster')
    motif_data[, row := seq(1, .N)]
    subset = motif_data[motif_data[, sample(.I, size, replace = TRUE, prob = count)]]
    subset[, count := .N, by = .(row)]
    subset_final = unique(subset[, ..cols])

    # fill in unobserved seq cases
    subset_orig_small = motif_data[, ..cols][, -c('count')]
    subset_final_small = subset_final[, ..cols][, -c('count')]
    unsampled = fsetdiff(subset_orig_small, subset_final_small) 
    unsampled$count = 0
    subset_final = rbind(subset_final, unsampled)

    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(subset_final, gene_type = gene_type, trim_type = trim_type)
    return(sample_data)
}
