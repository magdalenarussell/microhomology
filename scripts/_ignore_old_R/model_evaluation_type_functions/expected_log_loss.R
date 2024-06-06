HELD_OUT_FRACTION <<- 0.4
REPETITIONS <<- 20

get_unique_genes_col <- function(motif_data, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    if (length(genes) == 1){
        motif_data$unique_gene = motif_data[[paste0(gene_type, '_group')]] 
    } else if (length(genes) == 2){
        motif_data$unique_gene = interaction(motif_data[[paste0(genes[1], '_group')]], motif_data[[paste0(genes[2], '_group')]]) 
    }
    return(motif_data)
}

generate_hold_out_sample <- function(motif_data, sample_size, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    genes = get_gene_order(gene_type)
    motif_data = get_unique_genes_col(motif_data, gene_type)
    unique_g = unique(motif_data$unique_gene)

    stopifnot(sample_size < length(unique_g))

    sampled_vector = c()
    all_genes = unique_g
    while (length(sampled_vector) < sample_size){
        # sample size is the number of gene, subject combos
        sample = sample(unique_g, 1, replace = FALSE)
        if (genes > 1){
            expanded_sample_v = unique(motif_data[unique_gene == sample]$v_gene_group)
            expanded_sample_j = unique(motif_data[unique_gene == sample]$j_gene_group)
            expanded_sample = unique(motif_data[v_gene_group == expanded_sample_v | j_gene_group == expanded_sample_j]$unique_gene)
            sampled_vector = unique(c(sampled_vector,as.character(expanded_sample)))
        } else {
            sampled_vector = c(sampled_vector, sample)
        }
        unique_g = unique_g[unique_g != sample]
    }

    unsampled_vector = setdiff(all_genes, sampled_vector)

    motif_data_subset = motif_data[unique_gene %in% unsampled_vector]
    sample_data = motif_data[unique_gene %in% sampled_vector]

    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    motif_data_subset = calculate_subject_gene_weight(motif_data_subset, gene_type = gene_type, trim_type = trim_type)
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/p_gene_pooled.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(sample_data, gene_type = gene_type, trim_type = trim_type)
    return(list(sample = sample_data[, -c('unique_gene')], motif_data_subset = motif_data_subset[, -c('unique_gene')]))
}

get_hold_out_sample_probability <- function(sample_size, motif_data, gene_type = GENE_NAME){
    motif_data = get_unique_genes_col(motif_data, gene_type)
    unique_g = unique(motif_data$unique_gene)

    total_genes = length(unique_g)

    prob = (1/total_genes)^sample_size
    return(prob)
}

evaluate_loss <- function(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME) {
    set.seed(66)
    motif_data = get_unique_genes_col(motif_data, gene_type)
    unique_g = unique(motif_data$unique_gene)
    gene_count = length(unique_g)
    held_out_gene_count = round(HELD_OUT_FRACTION*gene_count)
    log_loss_vector = c()
    sample_prob_vector = c()
    parameter_count_vector = c()
    for (rep in 1:REPETITIONS){
        # Generate a held out sample and motif data subset
        sample_data = generate_hold_out_sample(motif_data, sample_size = held_out_gene_count, gene_type = gene_type, trim_type = trim_type) 
        motif_data_subset = sample_data$motif_data_subset
        sample = sample_data$sample

        # Fit model to the motif_data_subset
        if (MODEL_TYPE != 'null'){
            model = fit_model(motif_data_subset, trim_type = trim_type)
            parameter_count_vector = c(parameter_count_vector, length(model$coefficients)) 
        } else {
            model = 'null'
            parameter_count_vector = c(parameter_count_vector, 0)
        }

        # Compute conditional logistic loss value for held out sample using model
        log_loss = calculate_cond_expected_log_loss(model, sample)
        log_loss_vector = c(log_loss_vector, log_loss)

        # Compute probability of held out sample
        prob = get_hold_out_sample_probability(held_out_gene_count, motif_data)
        sample_prob_vector = c(sample_prob_vector, prob)
    }

    expected_log_loss = sum((1/REPETITIONS) * log_loss_vector)

    return(list(loss = expected_log_loss, model_parameter_count = unique(parameter_count_vector), held_out_cluster_number = NA, held_out_genes = 'averaged', vect = log_loss_vector))
}
