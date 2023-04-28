source(paste0(PROJECT_PATH, '/scripts/processing_functions/data_type_functions/', DATA_TYPE, '.R'))

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

process_data <- function(data, filter_frequency = FALSE){
    data = reformat_data(data)
    filtered = filter_data(data, filter_frequency)
    return(filtered)
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

    cols = c('v_trim', 'j_trim', 'v_gene', 'j_gene', 'vj_insert', 'productive', 'v_pnuc', 'j_pnuc')
    cols = cols[cols %in% colnames(subset)]
    subset = subset[, ..cols]

    if (TRIM_TYPE %like% 'trim'){
        # bound to reasonable trims
        subset = subset[get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
        subset = subset[get(TRIM_TYPE) >= LOWER_TRIM_BOUND]
    }

    # remove delta
    subset = subset[!(get(GENE_NAME) %like% 'TRD')]
    subset = subset[!(get(JOINING_GENE) %like% 'TRD')]

    subset[, total_count := .N, by = .(productive)]
    if (filter_frequency == TRUE){
        subset = filter_by_freq(subset, output_uncondensed = TRUE)
    }
    return(subset)
}

get_unobserved_observations <- function(data){
    # get unobserved observations
    cols = c('v_gene', 'j_gene', TRIM_TYPE, 'total_count')
    no_trim_cols = cols[!(cols %in% c(TRIM_TYPE, 'total_count'))]
    count_pair = unique(data[, c('total_count')])
    # get unique genes
    unique_obs = unique(data[,..no_trim_cols])

    # get desired trim lengths, genes
    trim_lengths = seq(0,max(data[[TRIM_TYPE]]))
    desired_obs = unique_obs %>%
            mutate(trim_length = list(trim_lengths)) %>%
            unnest(cols = c(trim_length)) %>%
            as.data.table()
    setnames(desired_obs, 'trim_length', TRIM_TYPE)
    # get unobserved subset
    if (nrow(count_pair) == 1){
        desired_obs$total_count = count_pair$total_count
    } else {
        desired_obs = merge(desired_obs, count_pair)
    }
    unobserved = desired_obs[!data, on = cols]
    unobserved$count = 0
    unobserved$productive = PRODUCTIVITY
    tog = rbind(data, unobserved)
    return(tog)
}

condense_data <- function(data, filter = TRUE){
    subset = filter_by_freq(data, filter)
    
    # get average trimming profiles
    cols = c(TRIM_TYPE, 'v_gene', 'j_gene', 'productive')

    subset[, p_trim_given_gene := count/paired_gene_count]
    avg = subset[, mean(p_trim_given_gene), by = cols]     
    setnames(avg, 'V1', 'avg_prob')
    return(avg)
}

filter_by_freq <- function(data, filter = TRUE, output_uncondensed = FALSE){
    orig_data = data
    # get counts
    cols = c(TRIM_TYPE, 'v_gene', 'j_gene', 'productive', 'total_count')
    data = data[, .N, by = cols]
    setnames(data, 'N', 'count')

    # get unobserved data
    if (all(is.numeric(data[[TRIM_TYPE]]))){
        data = get_unobserved_observations(data)
    }

    # get frequencies
    data[, v_gene_count := sum(count), by = .(v_gene, productive)]
    data[, j_gene_count := sum(count), by = .(j_gene, productive)]
    data[, paired_gene_count := sum(count), by = .(v_gene, j_gene, productive)]
    data[, paired_freq := paired_gene_count/total_count]

    # filter by frequencies
    max_trims = length(unique((data[[TRIM_TYPE]])))
    if (isTRUE(filter)){
        subset = data[paired_gene_count > 5000]
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
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive)]
    } else {
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive, draw)]
    }
    setnames(temp, 'V1', 'sum_abs_diff')
    return(temp)
}

bootstrap_diffs <- function(uncondensed_data, filter = TRUE){
    sample = uncondensed_data[uncondensed_data[, sample(.I, .N, replace = TRUE)]$V1]
    condensed = condense_data(sample, filter) 
    return(list(condensed = condensed, uncondensed = sample))
}

fit_multinom_models <- function(uncondensed_data, joining_gene = JOINING_GENE){
    require(nnet)
    null_form = as.formula(paste0('as.factor(', TRIM_TYPE, ') ~ 1'))
    form = as.formula(paste0('as.factor(',TRIM_TYPE, ') ~ as.factor(', joining_gene, ')'))

    if (length(unique(uncondensed_data[[TRIM_TYPE]])) > 2){
        null = multinom(null_form, data = uncondensed_data, maxit = 1000, MaxNWts=5000)
        varying = multinom(form, data = uncondensed_data, maxit = 1000, MaxNWts=5000)
    } else {
        null = glm(null_form, data = uncondensed_data, family = 'binomial')
        varying = glm(form, data = uncondensed_data, family = 'binomial')
    }
    return(list(null = null, model = varying))
}

get_pval <- function(null_model, varying_model){
    if (any(class(varying_model) %like% 'glm')){
        lrt = anova(null_model, varying_model, test = 'LRT')
    } else {
        lrt = anova(null_model, varying_model)
    }
    if ('LR stat.' %in% names(lrt)){
        result = data.table(p = lrt[['Pr(Chi)']][2], LRstat = lrt[['LR stat.']][2], df = lrt[['   Df']][2], null_lik = logLik(null_model), varying_lik = logLik(varying_model))
    } else {
        result = data.table(p = lrt[['Pr(>Chi)']][2], LRstat = lrt[['Deviance']][2], df = lrt[['Df']][2], null_lik = logLik(null_model), varying_lik = logLik(varying_model))
    }
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
    df = data.table(temp= unique(uncondensed_data[[JOINING_GENE]]))
    setnames(df, 'temp', JOINING_GENE)

    if (any(class(model) %like% 'glm')){
        pred = as.data.table(predict(model, newdata = df, type = 'response'))
    } else {
        pred = as.data.table(predict(model, newdata = df, type = 'probs'))
    }

    pred = cbind(df, pred)

    long_cols = c(JOINING_GENE)

    pred_long = pred %>%
        pivot_longer(!any_of(long_cols), names_to = TRIM_TYPE, values_to = 'prob')%>%
        as.data.table()
    
    cols = c(TRIM_TYPE, JOINING_GENE)
    avg = pred_long[, mean(prob), by = cols]
    setnames(avg, 'V1', 'avg_prob')
    avg[[GENE_NAME]] = unique(uncondensed_data[[GENE_NAME]])
    return(avg)
}

get_frames_for_pair <- function(v_seq, v_frame, v_trim, j_seq, j_frame, j_trim){
    # get processed V-gene sequence
    adjusted_v_seq = substring(v_seq, v_frame) 
    if (v_trim < 0){
        adjusted_v_seq = paste0(adjusted_v_seq, paste0(rep('X', -1*v_trim), collapse = ''), collapse = '')
    } else {
        adjusted_v_seq = substring(adjusted_v_seq, 1, nchar(adjusted_v_seq) - v_trim)
    }

    # get processed J-gene sequence
    j_extra = (nchar(j_seq) - (j_frame - 1))%%3
    adjusted_j_seq = substring(j_seq, 1, nchar(j_seq) - j_extra)
    if (j_trim < 0){
        adjusted_j_seq = paste0(paste0(rep('X', -1*j_trim), collapse = ''), adjusted_j_seq, collapse = '')
    } else {
        adjusted_j_seq = substring(adjusted_j_seq, j_trim + 1)
    }

    tog = paste0(adjusted_v_seq, adjusted_j_seq, collapse = '')
    frame = nchar(tog) %% 3
    if (frame == 0){
        frame_type = 'In'
    } else {
        frame_type = 'Out'
    }
    return(frame_type)
}

get_all_frames <- function(){
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv')[organism == 'human' & chain == substring(LOCUS, 3, 3)]

    v = frames[region == 'V'][, c('id', 'frame', 'nucseq')]
    colnames(v) = c('v_gene', 'v_frame', 'v_seq')
    j = frames[region == 'J'][, c('id', 'frame', 'nucseq')]
    colnames(j) = c('j_gene', 'j_frame', 'j_seq')

    v$dummy = 1
    j$dummy = 1

    pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

    trims = data.table(v_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), each = (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND) + 1), j_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND) + 1))
    trims$dummy = 1

    pairs = merge(pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]

    pairs[, frame_type := apply(.SD, 1, function(x) get_frames_for_pair(x[1], x[2], x[3], x[4], x[5], x[6])), .SDcols = c("v_seq", "v_frame", "v_trim", "j_seq", "j_frame", "j_trim")]

    return(pairs[, -c('v_seq', 'j_seq')])
}

get_frames_data <- function(){
    path = file.path(PROJECT_PATH, 'data')
    dir.create(path)
    file_name = file.path(path, paste0(LOCUS, '_paired_frame_data.tsv'))

    if (!file.exists(file_name)) {
        frame_data = get_all_frames()
        fwrite(frames_data, file_name, sep = '\t')
    } else {
        frame_data = fread(file_name)
    }
    return(frame_data)
}
