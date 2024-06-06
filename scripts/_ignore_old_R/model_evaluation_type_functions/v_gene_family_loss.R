HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

source(paste0(MOD_PROJECT_PATH, '/scripts/model_evaluation_type_functions/evaluation_type_classes/gene_family.R'))

evaluate_loss <- function(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME) {
    gene_families = get_gene_families(cluster_count = 3, combine_by_terminal = TRUE, full_sequence = FALSE, align = FALSE, gene_type = 'v_gene')$cluster_data
    cluster_counts = gene_families[, .N, by = clusters_grouped]
    largest_cluster = cluster_counts[N == max(N)]$clusters_grouped
    big_enough_clusters = cluster_counts[N > 1]$clusters_grouped
    held_out_clusters = unique(gene_families[clusters_grouped != largest_cluster]$clusters_grouped)
    held_out_clusters = held_out_clusters[held_out_clusters %in% big_enough_clusters]
    log_loss_vector = c()
    parameter_count_vector = c()
    genes = c()
    clusters = c()
    for (cluster_group in c(as.list(held_out_clusters), list(held_out_clusters))) {
        held_out_genes = unique(gene_families[clusters_grouped %in% cluster_group][['v_gene_group']])
        if (length(unique(substring(held_out_genes, 1, 6))) == 1) {
            next
        }
        # Generate a held out sample and motif data subset
        sample_data = generate_hold_out_sample(motif_data, held_out_genes, gene_type = 'v_gene', trim_type = trim_type) 
        motif_data_subset = sample_data$motif_data_subset
        sample = sample_data$sample

        if (nrow(motif_data_subset[weighted_observation != 0]) < 50000){
            print('skipping cluster: number of non-zero observations is too small for model fitting')
            next
        }

        # Fit model to the motif_data_subset
        if (MODEL_TYPE != 'null'){
            model = fit_model(motif_data_subset, trim_type= trim_type)
            parameter_count_vector = c(parameter_count_vector, length(model$coefficients)) 
        } else {
            model = 'null'
            parameter_count_vector = c(parameter_count_vector, 0)
        }

        # Compute conditional logistic loss value for held out sample using model
        log_loss = calculate_cond_expected_log_loss(model, sample)
        log_loss_vector = c(log_loss_vector, log_loss)
       
        held_out_genes = paste(held_out_genes, collapse = ', ')
        genes = c(genes, held_out_genes)
        cluster_string = paste(cluster_group, collapse = ', ')
        clusters = c(clusters, cluster_string)
    }
    return(list(loss = log_loss_vector, model_parameter_count = parameter_count_vector, held_out_cluster_number = clusters, held_out_genes = genes))
}
