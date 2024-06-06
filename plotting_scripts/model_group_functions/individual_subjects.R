get_predicted_dist_figure_file_name <- function(gene_name, subgroup){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(subgroup, '_', new_gene_name, '_dist.pdf')
    return(filename)
}

get_coef_heatmap_file_name <- function(type, subgroup){
    filename = paste0(subgroup, '_', type, '_heatmap.pdf')
    return(filename)
}

plot_predicted_trimming_dists <- function(data, gene_name){
    file_path = get_predicted_dist_figure_file_path()
    for (indiv in unique(data$subject)){
        subset_data = data[subject == indiv]
        file_name = get_predicted_dist_figure_file_name(gene_name, indiv)
        complete_path = file.path(file_path, file_name)
        plot_predicted_trimming_dists_single_group(subset_data, gene_name, complete_path)
        print(paste0('finished plotting for ', indiv))
    }
}

plot_model_coefficient_heatmap <- function(model_coef_matrix, with_values = FALSE){
    file_path = get_coef_heatmap_file_path()
    limits = range(sapply(model_coef_matrix[,-c('base', 'model_group')], range, na.rm = TRUE))/log(10)
    for (indiv in unique(model_coef_matrix$model_group)){
        subset_data = model_coef_matrix[model_group == indiv]
        if (MODEL_TYPE %like% 'two_side_terminal_melting'){
            plot_melting_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('melting', indiv)), with_values = with_values, write_plot = write_plot, limits = limits)
        } 
        if (MODEL_TYPE %like% 'distance'){
            plot_distance_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('distance', indiv)), with_values = with_values, write_plot = write_plot, limits = limits)
        }
        if (MODEL_TYPE %like% 'motif'){
            plot_model_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('motif', indiv)), with_values = with_values, write_plot = write_plot, limits = limits)
        }
        if (MODEL_TYPE %like% 'dna_shape'){
            plot_shape_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('dna_shape', indiv)), with_values = with_values, write_plot = write_plot, limits = limits)
        }
        if (MODEL_TYPE %like% 'count'){
            plot_base_count_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('count', indiv)), with_values = with_values, write_plot = write_plot, limits = limits)
        }
        print(paste0('finished plotting for ', indiv))
    }
}
