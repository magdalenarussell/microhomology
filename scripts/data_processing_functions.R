get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_alpha')))[, -c('name')]
    return(whole_nucseq[, c('gene', 'sequence')])
}

get_similarity_distance <- function(nt_count, gene_type = 'J', whole_nucseqs){
    # TODO this function only works for J-genes
    require(DECIPHER)
    require(Biostrings)
    tr = whole_nucseqs[gene %like% gene_type]
    tr[, short:= substring(sequence, 1, nt_count)]
    dists = DistanceMatrix(DNAStringSet(tr$short))

    tr[, dist := rowMeans(dists)]
    tr[, dist_nt_count := nt_count]
    return(tr)
}

read_all_data <- function(directory){
    require(vroom)
    files = list.files(directory, pattern = 'tsv', full.names = TRUE)
    data = vroom(files)
    data = as.data.table(data)
    return(data)
}

convert_adaptive_gene_names_to_imgt <- function(data, gene_col){
    mapping = fread('https://raw.githubusercontent.com/kmayerb/tcrdist3/master/tcrdist/db/adaptive_imgt_mapping.csv')[species == 'human']
    gene_type = toupper(substring(gene_col, 1, 1))
    data[, c(gene_col, "allele") := tstrsplit(get(gene_col), "\\*0")]
    tog = merge(data, mapping, by.x = gene_col, by.y = 'adaptive', all.x = TRUE, allow.cartesian = TRUE)
    setnames(tog, 'imgt', paste0(tolower(gene_type), '_gene'))
    return(tog[, -c('allele', 'species')])
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
    # convert genes 
    for (col in c('v_resolved', 'd_resolved', 'j_resolved')){
        data = convert_adaptive_gene_names_to_imgt(data, col)
    } 

    # convert productivity
    data = convert_frame_type_to_productivity(data)

    # convert trims
    for (trim in c('v_', 'd3_', 'd5_', 'j_')){
        data = convert_trims(data, paste0(trim, 'deletions'))
        data[[paste0(trim, 'trim')]] = as.numeric(data[[paste0(trim, 'trim')]])
    }
    return(data)
}

filter_data <- function(data){
    # filter by productivity
    subset = data[productive == PRODUCTIVITY]
    
    # reomve observations containing a d_gene
    subset = subset[is.na(d_gene)]

    # remove observations with missing genes
    for (gene in c('v_gene', 'j_gene')){
        subset = subset[!is.na(get(gene))] 
    }

    cols = c('sample_name', 'v_trim', 'j_trim', 'v_gene', 'j_gene', 'productive')
    subset = subset[, ..cols]

    # bound to reasonable trims
    subset = subset[get(TRIM_TYPE) <= 24]

    # remove delta
    subset = subset[!(get(GENE_NAME) %like% 'TRD')]
    return(subset)
}

get_unobserved_observations <- function(data){
    # get unobserved observations
    cols = c('v_gene', 'j_gene', 'sample_name', TRIM_TYPE)
    no_trim_cols = cols[!(cols %in% c('sample_name', TRIM_TYPE))]
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
    unobserved = desired_obs[!data, on = cols]
    unobserved$count = 0
    unobserved$productive = PRODUCTIVITY
    tog = rbind(data, unobserved)
    return(tog)
}

condense_data <- function(data){
    subset = filter_by_freq(data)
    # get average trimming profiles
    cols = c(TRIM_TYPE, 'v_gene', 'j_gene', 'productive', 'avg_paired_freq', 'subj_count')

    subset[, p_trim_given_gene := count/paired_gene_count]
    avg = subset[, mean(p_trim_given_gene), by = cols]     
    setnames(avg, 'V1', 'p_trim_given_pair')
    return(avg)
}

filter_by_freq <- function(data){
    # get counts
    cols = c('sample_name', TRIM_TYPE, 'v_gene', 'j_gene', 'productive')
    data = data[, .N, by = cols]
    setnames(data, 'N', 'count')

    # get unobserved data
    data = get_unobserved_observations(data)

    # get frequencies
    data[, total_count := sum(count), by = .(sample_name, productive)]
    data[, v_gene_count := sum(count), by = .(sample_name, v_gene, productive)]
    data[, j_gene_count := sum(count), by = .(sample_name, j_gene, productive)]
    data[, paired_gene_count := sum(count), by = .(sample_name, v_gene, j_gene, productive)]
    data[, paired_freq := paired_gene_count/total_count]

    # filter by frequencies
    subset = data[paired_freq > 5e-04]
    subset[, subj_count := .N/25, by = .(v_gene, j_gene, productive)]
    subset = subset[subj_count >= 3]
    subset[, avg_paired_freq := mean(paired_freq), by = .(v_gene, j_gene, productive)]
    return(subset)
}

get_original_mean_curve <- function(original_condensed_data, new_condensed_data){
    original_condensed_data[, mean_p_trim_given_pair := mean(p_trim_given_pair), by = .(v_trim, v_gene, productive)]
    cols = c('mean_p_trim_given_pair', 'v_gene', 'v_trim', 'productive')
    new_condensed_data = merge(new_condensed_data, unique(original_condensed_data[, ..cols])) 
    return(new_condensed_data)
}

calculate_mean_of_abs_diffs <- function(condensed_data, original_condensed_data, repeated = FALSE){
    condensed_data = get_original_mean_curve(original_condensed_data, condensed_data)
    condensed_data[, abs_diff := abs(mean_p_trim_given_pair - p_trim_given_pair)]
    if (isFALSE(repeated)){
        temp = condensed_data[, mean(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count)]
    } else {
        temp = condensed_data[, mean(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count, draw)]
    }
    setnames(temp, 'V1', 'mean_abs_diff')
    return(temp)
}


calculate_sum_of_abs_diffs <- function(condensed_data, original_condensed_data, repeated = FALSE){
    condensed_data = get_original_mean_curve(original_condensed_data, condensed_data)
    condensed_data[, abs_diff := abs(mean_p_trim_given_pair - p_trim_given_pair)]
    if (isFALSE(repeated)){
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count)]
    } else {
        temp = condensed_data[, sum(abs_diff), by = .(v_gene, j_gene, productive, avg_paired_freq, subj_count, draw)]
    }
    setnames(temp, 'V1', 'sum_abs_diff')
    return(temp)
}

bootstrap_diffs <- function(set_of_joining_genes, empirical_data){
    sample = sample(set_of_joining_genes, length(set_of_joining_genes), replace = TRUE)
    random_dt = data.table()
    random_joins = data.table()
    for (rep in seq(1, length(set_of_joining_genes))){
        rep_sample = sample[rep]
        temp = empirical_data[get(JOINING_GENE) == rep_sample][, -c('abs_diff')]
        temp$draw = rep
        join_temp = whole_nucseqs[gene == rep_sample]
        random_joins = rbind(random_joins, join_temp)
        random_dt = rbind(random_dt, temp)
    }

    return(list(random_whole_nucseqs = random_joins, random = random_dt))
}

run_diff_bootstrap <- function(set_of_joining_genes, empirical_data, boot_count){
    results = data.table()
    results_mean = data.table()
    for (boot in seq(boot_count)){
        bootstrap = bootstrap_diffs(set_of_joining_genes, empirical_data)
        random_dt = bootstrap$random
        random_joins = bootstrap$random_whole_nucseqs
        random_diverged = calculate_sum_of_abs_diffs(random_dt, random_dt, repeated = TRUE)
        random_mean_diverged = calculate_mean_of_abs_diffs(random_dt, empirical_data, repeated = TRUE)
        # random_dists = get_similarity_distance(nt_count = 10, gene_type = toupper(substring(JOINING_GENE, 1, 1)), random_joins)
        # random_diverged = merge(random_diverged, unique(random_dists), by.x = JOINING_GENE, by.y = 'gene', all.x = TRUE)
        random_diverged$bootstrap = boot
        random_mean_diverged$bootstrap = boot

        results = rbind(results, random_diverged)
        results_mean = rbind(results_mean, random_mean_diverged)

        print(paste0('bootstrap ', boot, ' complete'))
    }
    return(list(sum = results, mean = results_mean))
}

run_t_test <- function(boot_results, mean_divergence, sample_size){
    df = sample_size - 1
    sd = sd(boot_results$mean_divergence) 
    tstat = (mean_divergence)/(sd/sqrt(sample_size)) 
    # one-tail t test
    p = pt(tstat, df = df, lower.tail = FALSE)
    result = data.table(p = p, mean = mean_divergence, sd = sd, t_statistic = tstat, bootstraps = nrow(boot_results), total_joining_gene_count = sample_size)
    result[[GENE_NAME]] = unique(boot_results[[GENE_NAME]])
    return(result)
}
