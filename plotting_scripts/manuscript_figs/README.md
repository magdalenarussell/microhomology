# Manuscript figures

This directory contains a script to produce each manuscript figure.
The scripts map to manuscript figures as follows:

- [Figure 2](training_loss/multi_parameter_loss_plot.R): figure containing loss values for many different models across various held-out data sets
- [Figure 3](motif_base_count_model/motif_base_count_model.R): figure containing trimming distributions and coefficient heatmaps for the _motif + two-side base-count-beyond_ model
- [Figure 4](motif_exploration/most_improved.R): figure comparing the _two-side base-count_ model to the _motif + two-side base-count-beyond_ model
- [Figure 5](snp_interaction/snp-interaction_coefs_motif_two-side-base-count-beyond.R): figure containing heatmap of Artemis SNP, model parameter interaction coefficients
- [Figure 6](validation_loss/validation.R): figure containing loss values for many validation data sets (nonproductive sequences only)

The remaining scripts produce supplementary figures:

- [Figure S1](supp_figs/murugan_replication.R): figure replicating previous model from Murugan et. al PNAS 2012
- [Figure S2](supp_figs/single_parameter_coefs.R): figure showing model coefficients for the _two-side base-count_ model
- [Figure S3](supp_figs/gallery.R): figure showing inferred distributions for all TRB V-genes with the _motif + two-side base-count-beyond_ model
- [Figure S4](supp_figs/motif_base_count_motif_size.R): figure comparing losses by motif size for the _motif + two-side base-count-beyond_ model
- [Figure S5](supp_figs/motif_base_count_ss_model.R): figure showing trimming distributions and coefficient heatmaps for the _motif + two-side base-count-beyond_ model when single-stranded nucleotides are allowed to be counted in model parameters
- [Figure S6](supp_figs/hairpin.R): figure exploring the effects of different hairing-opening position assumptions on the performance of the _motif + two-side base-count-beyond_ model
- [Figure S7](supp_figs/motif_base_count_model_vj.R): figure containing trimming distributions and coefficient heatmaps for the _motif + two-side base-count-beyond_ model when the model was trained using J-gene sequences from the training data set
- [Figure S8](supp_figs/snp-interaction_coefs_motif_linear-distance_two-side-base-count-beyond.R):figure containing heatmap of Artemis SNP, model parameter interaction coefficients when a distance term was also included in the model
- [Figure S9](supp_figs/parsimony_annotation.R): model coefficient heatmap for the _motif + two-side base-count-beyond_ model when using a parsimony-based annotation method for the V-gene training data set (instead of the IGoR-based annotation method)
- [Figure S10](supp_figs/validation.R): figure containing loss values for many validation data sets (productive sequences only)
- [Figure S11](supp_figs/motif_frequency.R): figure showing the frequency of each 3-nucleotide motif within the TRB and IGH germline
- [Figure S12](supp_figs/rel_importance.R): figure showing the relative weights of the motif and two-side-base-count terms across each testing data set
- [Figure S16](supp_figs/gene_trees.R): figure showing hamming-distance based trees for V-gene sequences
