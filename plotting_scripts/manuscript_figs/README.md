# Manuscript figures

This directory contains a script to produce each manuscript figure.
The scripts map to manuscript figures as follows:

* [Figure 2A](motif_base_count_average-mh_ligation-mh_model/nonproductive_v-j_trim_ligation-mh/motif_base_count_mh_model_indiv_plots.R): figures containing model parameter values in heatmap format
* [Figure 2C](motif_base_count_average-mh_ligation-mh_model/nonproductive_v-j_trim_ligation-mh/rel_size.R): figure showing model parameter values ordered by relative size
* [Figure 3](validation_loss/dropout.R): figure showing variable dropout experiments to show changes in log loss and MAE
* [Figure 4](analysis_annotation_ranking/compare_annotation_ranking.R): figures showing effects of including MH terms on annotation probabilities and rankings

The remaining scripts produce supplementary figures:

- [Figure S1](motif_base_count_average-mh_ligation-mh_model/nonproductive_v-j_trim_ligation-mh/convergence.py): figure illustrating EM convergence
- [Figure S2](ligation_mh_choices_germline_tra/choice_options.R): heatmap showing distribution of complementary sequence regions across trimming configurations
- [Figure S3](ligation_mh_choices_germline_tra/choice_options.R): scatterplot comparing average number of ligation configurations with the average number of microhomologous nucleotides per trimming configuration per gene pair
- [Figure S4](ligation_mh_choices_germline_tra/choice_options.R): heatmap showing distribution of ligation configuration options across trimming configurations
- [Figure S5](mh_simulator_experiments/nonproductive_v-j_trim_ligation-mh/motif_base_count_mh_model.R): heatmaps of model parameter values for microhomology simulation experiments
- [Figure S6](motif_base_count_interior-mh-count_model/motif_base_count_mh_model.R): heatmaps of model parameter values for model trained on sequences with N-insertions which lack germline-dependent microhomology-mediated ligation
- [Figure S7]: comparison of the number of microhomologous nucleotides in the configuration and the difference model probabilities using the MH model versus no MH model 
