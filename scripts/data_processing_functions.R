read_all_data <- function(directory){
    require(vroom)
    all_data = data.table()
    files = list.files(directory, pattern = 'sv', full.names = TRUE)
    if (length(files) > 1){
        iter = ceiling(length(files)/1000)
        for (it in seq(iter)){
            start = (it-1)*1000 + 1
            end = min(it*1000, length(files))
            subset = files[start:end]
            data = vroom(subset)
            data = as.data.table(data)
            all_data = rbind(all_data, data)
        }
    } else {
        all_data = fread(files)
    }
    return(all_data)
}

convert_adaptive_gene_names_to_imgt <- function(data, gene_col){
    if (data[[gene_col]][1] %like% 'TCR'){
        mapping = fread('https://raw.githubusercontent.com/kmayerb/tcrdist3/master/tcrdist/db/adaptive_imgt_mapping.csv')[species == 'human']
        gene_type = toupper(substring(gene_col, 1, 1))
        data[, c(gene_col, "allele") := tstrsplit(get(gene_col), "\\*0")]
        tog = merge(data, mapping, by.x = gene_col, by.y = 'adaptive', all.x = TRUE, allow.cartesian = TRUE)
        setnames(tog, 'imgt', paste0(tolower(gene_type), '_gene'))
        return(tog[, -c('allele', 'species')])
    } else {
        return(data)
    }
}

convert_frame_type_to_productivity <- function(data){
    data[frame_type %in% c('Out', 'Stop'), productive := 'nonproductive']
    data[frame_type == 'In', productive := 'productive']
    return(data[, -c('frame_type')])
}

convert_trims <- function(data, trim_col){
    trim_type = substring(trim_col, 1, 1)
    if (trim_type == 'd'){
        trim_type = substring(trim_col, 1, 2)
    }
    data[get(trim_col) == 'no data', paste(trim_col) := 0]
    setnames(data, trim_col, paste0(trim_type, '_trim'))
    return(data)
}

convert_adaptive_style_to_imgt <- function(data){
    adaptive_cols = c('v_resolved', 'd_resolved', 'j_resolved')
    if (any(adaptive_cols %in% colnames(data))){
        # convert genes 
        for (col in adaptive_cols){
            data = convert_adaptive_gene_names_to_imgt(data, col)
        } 

        # convert productivity
        data = convert_frame_type_to_productivity(data)

        # convert trims
        for (trim in c('v_', 'd3_', 'd5_', 'j_')){
            data = convert_trims(data, paste0(trim, 'deletions'))
            data[[paste0(trim, 'trim')]] = as.numeric(data[[paste0(trim, 'trim')]])
        }

        data[n1_insertions == 'no data', n1_insertions := 0]
        data$n1_insertions = as.numeric(data$n1_insertions)
        data = create_processing_vars(data)
        return(data)
    } else {
        return(data)
    }
}

filter_data <- function(data, filter_frequency = FALSE){
    # filter by productivity
    subset = data[productive == PRODUCTIVITY]
    
    # reomve observations containing a d_gene
    if ('d_gene' %in% colnames(subset)){
        subset = subset[is.na(d_gene)]
    }

    # remove observations with missing genes
    for (gene in c('v_gene', 'j_gene')){
        subset = subset[!is.na(get(gene))] 
    }

    cols = c('sample_name', 'v_trim', 'j_trim', 'v_gene', 'j_gene', 'n1_insertions', 'zero_process', 'zero_insert', 'productive')
    subset = subset[, ..cols]

    if (TRIM_TYPE %like% 'trim'){
        # bound to reasonable trims
        subset = subset[get(TRIM_TYPE) <= 24]
        subset[get(TRIM_TYPE) < 0 , paste(TRIM_TYPE) := 0]
    }

    # remove delta
    subset = subset[!(get(GENE_NAME) %like% 'TRD')]
    subset = subset[!(get(JOINING_GENE) %like% 'TRD')]

    subset[, total_count := .N, by = .(sample_name, productive)]
    if (filter_frequency == TRUE){
        subset = filter_by_freq(subset, output_uncondensed = TRUE)
    }
    return(subset)
}

get_unobserved_observations <- function(data){
    # get unobserved observations
    cols = c('v_gene', 'j_gene', 'sample_name', TRIM_TYPE, 'total_count')
    no_trim_cols = cols[!(cols %in% c('sample_name', TRIM_TYPE, 'total_count'))]
    count_pair = unique(data[, c('sample_name', 'total_count')])
    # get unique genes
    unique_obs = unique(data[,..no_trim_cols])

    # get desired trim lengths, genes
    trim_lengths = seq(0,max(data[[TRIM_TYPE]]))
    desired_obs = unique_obs %>%
            mutate(trim_length = list(trim_lengths)) %>%
            mutate(sample_name = list(unique(data$sample_name))) %>%
            unnest(cols = c(trim_length)) %>%
            unnest(cols = c(sample_name)) %>%
            as.data.table()
    setnames(desired_obs, 'trim_length', TRIM_TYPE)
    # get unobserved subset
    desired_obs = merge(desired_obs, count_pair)
    unobserved = desired_obs[!data, on = cols]
    unobserved$count = 0
    unobserved$productive = PRODUCTIVITY
    tog = rbind(data, unobserved)
    return(tog)
}

condense_data <- function(data, filter = TRUE){
    subset = filter_by_freq(data, filter)
    
    # get average trimming profiles
    cols = c(TRIM_TYPE, 'v_gene', 'j_gene', 'productive', 'avg_paired_freq', 'subj_count')

    subset[, p_trim_given_gene := count/paired_gene_count]
    avg = subset[, mean(p_trim_given_gene), by = cols]     
    setnames(avg, 'V1', 'avg_prob')
    avg[avg_paired_freq < 5e-04, low_freq := TRUE]
    avg[avg_paired_freq >= 5e-04, low_freq := FALSE]
    return(avg)
}

filter_by_freq <- function(data, filter = TRUE, output_uncondensed = FALSE){
    orig_data = data
    # get counts
    cols = c('sample_name', TRIM_TYPE, 'v_gene', 'j_gene', 'productive', 'total_count')
    data = data[, .N, by = cols]
    setnames(data, 'N', 'count')

    # get unobserved data
    if (all(is.numeric(data[[TRIM_TYPE]]))){
        data = get_unobserved_observations(data)
    }

    # get frequencies
    data[, v_gene_count := sum(count), by = .(sample_name, v_gene, productive)]
    data[, j_gene_count := sum(count), by = .(sample_name, j_gene, productive)]
    data[, paired_gene_count := sum(count), by = .(sample_name, v_gene, j_gene, productive)]
    data[, paired_freq := paired_gene_count/total_count]

    # filter by frequencies
    max_trims = length(unique((data[[TRIM_TYPE]])))
    if (isTRUE(filter)){
        subset = data[paired_freq > 5e-04]
        subset[, subj_count := .N/max_trims, by = .(v_gene, j_gene, productive)]
        thresh = ceiling(0.75*length(unique(subset$sample_name)))
        subset = subset[subj_count >= thresh]
    } else {
        subset = data
        subset[, subj_count := .N/max_trims, by = .(v_gene, j_gene, productive)]
    }
    subset[, avg_paired_freq := mean(paired_freq), by = .(v_gene, j_gene, productive)]
    if (output_uncondensed == TRUE){
        subset[, temp := paste0(v_gene, j_gene)]
        orig_data[, temp := paste0(v_gene, j_gene)]
        pairs = unique(subset$temp)
        subset = orig_data[temp %in% pairs]
    }
    return(subset)
}

get_original_mean_curve <- function(original_condensed_data, new_condensed_data, gene_type = GENE_NAME){
    cols0 = c('v_trim', gene_type, 'productive')
    original_condensed_data[, mean_p_trim_given_pair := mean(avg_prob), by = cols0]
    cols = c('mean_p_trim_given_pair', gene_type, 'v_trim', 'productive')
    new_condensed_data = merge(new_condensed_data, unique(original_condensed_data[, ..cols])) 
    return(new_condensed_data)
}

calculate_sum_of_abs_diffs <- function(condensed_data, original_condensed_data, repeated = FALSE, gene_type = GENE_NAME){
    condensed_data = get_original_mean_curve(original_condensed_data, condensed_data, gene_type)
    condensed_data[, abs_diff := abs(mean_p_trim_given_pair - avg_prob)]
    if (isFALSE(repeated)){
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count)]
    } else {
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count, draw)]
    }
    setnames(temp, 'V1', 'sum_abs_diff')
    return(temp)
}

bootstrap_diffs <- function(uncondensed_data, filter = TRUE){
    uncondensed_data[, total := .N, by = .(sample_name)]
    sample = uncondensed_data[uncondensed_data[, sample(.I, .N, replace = TRUE), by = sample_name]$V1]
    condensed = condense_data(sample, filter) 
    return(list(condensed = condensed, uncondensed = sample))
}

fit_multinom_models <- function(uncondensed_data, joining_gene = JOINING_GENE){
    require(nnet)
    if (length(unique(uncondensed_data$sample_name)) > 1){
        null_form = as.formula(paste0(TRIM_TYPE, '~ sample_name'))
        form = as.formula(paste0(TRIM_TYPE, '~ sample_name + ', joining_gene))
    } else {
        null_form = as.formula(paste0('as.factor(', TRIM_TYPE, ') ~ 1'))
        form = as.formula(paste0('as.factor(',TRIM_TYPE, ') ~ as.factor(', joining_gene, ')'))
    }
    null = multinom(null_form, data = uncondensed_data, maxit = 1000, MaxNWts=5000)
    varying = multinom(form, data = uncondensed_data, maxit = 1000, MaxNWts=5000)
    return(list(null = null, model = varying))
}

get_pval <- function(null_model, varying_model){
    lrt = anova(null_model, varying_model)
    result = data.table(p = lrt[['Pr(Chi)']][2], LRstat = lrt[['LR stat.']][2], df = lrt[['   Df']][2], null_lik = logLik(null_model), varying_lik = logLik(varying_model))
    result[, p:= pchisq(LRstat, df = df, lower.tail = FALSE)]
    return(result)
}

get_output_path <- function(){
    path = file.path(OUTPUT_PATH, 'results', PRODUCTIVITY, LOCUS, paste0(TRIM_TYPE, '_joining_', JOINING_GENE))
    dir.create(path, recursive = TRUE)
    return(path)
}

get_multinom_file_name <- function(type = 'lrt_pvalue'){
    path = get_output_path()
    name = file.path(path, paste0(type, '.tsv'))
    return(name)
}

convert_model_coefs_to_dt <- function(model){
    coefs = coef(model)
    coefs_df = as.data.frame(coefs)
    coefs_df$trim_length = rownames(coefs_df)
    
    pivot = coefs_df %>%
        pivot_longer(!trim_length, names_to = 'coef', values_to = 'value') %>%
        as.data.table()
    return(pivot)
}

get_predicted_probs <- function(model, uncondensed_data){
    interaction = interaction(unique(uncondensed_data$sample_name), unique(uncondensed_data[[JOINING_GENE]]))
    df = data.table(temp = levels(interaction))
    df[, sample_name := sapply(temp, function(x) str_split(x, '\\.')[[1]][1])]
    df[, paste(JOINING_GENE) := sapply(temp, function(x) str_split(x, '\\.')[[1]][2])]

    pred = as.data.table(predict(model, newdata = df, type = 'probs'))
    pred = cbind(df, pred)

    long_cols = c('sample_name', JOINING_GENE)

    pred_long = pred[, -c('temp')] %>%
        pivot_longer(!any_of(long_cols), names_to = TRIM_TYPE, values_to = 'prob')%>%
        as.data.table()
    
    cols = c(TRIM_TYPE, JOINING_GENE)
    avg = pred_long[, mean(prob), by = cols]
    setnames(avg, 'V1', 'avg_prob')
    avg[[GENE_NAME]] = unique(uncondensed_data[[GENE_NAME]])
    return(avg)
}

create_processing_vars <- function(rep_data){
    # create var for zero processing
    rep_data[v_trim == 0 & j_trim == 0 & n1_insertions == 0, zero_process := TRUE]
    rep_data[!(v_trim == 0 & j_trim == 0 & n1_insertions == 0), zero_process := FALSE]

    # create var for zero insertions
    rep_data[n1_insertions == 0, zero_insert := TRUE]
    rep_data[!(n1_insertions == 0), zero_insert := FALSE]
    
    return(rep_data)
}
