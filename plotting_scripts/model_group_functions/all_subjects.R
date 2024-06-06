get_predicted_dist_figure_file_name <- function(gene_name, subgroup = NULL){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(new_gene_name, '_dist.pdf')
    return(filename)
}

get_coef_heatmap_file_name <- function(type, subgroup = NULL){
    filename = paste0(type, '_heatmap.pdf')
    return(filename)
}

plot_predicted_trimming_dists <- function(data, gene_name, ylim = NULL, color = 'blue', motif_highlight = "CTT", motif_highlight_color = 'black', seq_text = 7){
    file_path = get_predicted_dist_figure_file_path()
    file_name = get_predicted_dist_figure_file_name(gene_name)
    complete_path = file.path(file_path, file_name)
    plot_predicted_trimming_dists_single_group(data, gene_name, complete_path, ylim, color = color, motif_highlight_color = motif_highlight_color, motif_highlight = motif_highlight)
}

plot_model_coefficient_heatmap <- function(model_coef_matrix, with_values = FALSE, write_plot = TRUE, melt_limits = NULL, motif_limits = NULL, dist_limits = NULL, shape_limits = NULL, count_limits = NULL){
    file_path = get_coef_heatmap_file_path()
    if (MODEL_TYPE %like% 'two_side_terminal_melting'){
        plot_melting_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('melting')), with_values = with_values, write_plot = write_plot, limits = melt_limits)
    } 
    if (MODEL_TYPE %like% 'distance'){
        plot_distance_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('distance')), with_values = with_values, write_plot = write_plot, limits = dist_limits)
    }
    if (MODEL_TYPE %like% 'motif'){
        plot_model_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('motif')), with_values = with_values, write_plot = write_plot, limits = motif_limits)
    }
    if (MODEL_TYPE %like% 'dna_shape'){
        plot_shape_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('dna_shape')), with_values = with_values, write_plot = write_plot, limits = shape_limits)
    }
    if (MODEL_TYPE %like% 'count'){
        plot_base_count_coefficient_heatmap_single_group(model_coef_matrix = model_coef_matrix, file_name = file.path(file_path, get_coef_heatmap_file_name('count')), with_values = with_values, write_plot = write_plot, limits = count_limits)
    }
}
