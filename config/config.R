# (required) the following paths indicate the location of the project
# (MOD_PROJECT_PATH) and the location of where files should be written (MOD_OUTPUT PATH)
MOD_OUTPUT_PATH <<- '/fh/fast/matsen_e/shared/tcr-gwas/microhomology'
MOD_PROJECT_PATH <<- '/home/mrussel2/microhomology'
ROOT_PATH <<- '/home/mrussel2/microhomology' 
PROJECT_PATH <<- ROOT_PATH

#TODO update these
# (optional) Parsed TCR repertoire data (data located within the `emerson_parsed_TCRB.tgz` file available at https://doi.org/10.5281/zenodo.5719520)
# TCR_REPERTOIRE_DATA_parsimony = paste0(MOD_PROJECT_PATH, '/_ignore/emerson_stats/')
# (required) Training data set; change the path to these data after IGoR processing
# TCR_REPERTOIRE_DATA_igor = paste0(MOD_PROJECT_PATH, '/_ignore/emerson_igor_stats/')

TCR_REPERTOIRE_DATA_adaptive_alpha = paste0(MOD_PROJECT_PATH, '/_ignore/heikkila_alpha/adaptive/')
TCR_REPERTOIRE_DATA_igor_alpha = paste0(MOD_PROJECT_PATH, '/_ignore/heikkila_alpha/igor/')
TCR_REPERTOIRE_DATA_igor_sim_alpha = paste0(MOD_PROJECT_PATH, '/_ignore/igor_sim_alpha/')


# (required) TRB gene names and germline sequences from IMGT
WHOLE_NUCSEQS_beta = paste0(MOD_PROJECT_PATH,'/_ignore/igor_imgt_genes.tsv')

# (optional) Other gene name and germline sequences for other loci
#TODO remove
# WHOLE_NUCSEQS_adaptive_beta = paste0(MOD_PROJECT_PATH,'/_ignore/tcrb_processed_geneseq.tsv')
WHOLE_NUCSEQS_alpha = paste0(MOD_PROJECT_PATH,'/_ignore/imgt_tcrad_sequences.tsv')
WHOLE_NUCSEQS_gamma = paste0(MOD_PROJECT_PATH, '/_ignore/imgt_tcrg_sequences.tsv')
WHOLE_NUCSEQS_delta = paste0(MOD_PROJECT_PATH,'/_ignore/imgt_tcrad_sequences.tsv') 
WHOLE_NUCSEQS_igh = paste0(MOD_PROJECT_PATH,'/_ignore/imgt_igh_sequences.tsv')

# (optional) SNP genotype file (download https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001918.v1.p1)
SNP_GDS_FILE <<- paste0(MOD_PROJECT_PATH, "/_ignore/HSCT_comb_geno_combined_v03_tcr.gds")

# (optional) File to map scanIDs to localID (data located within the `gwas_id_mapping.tsv` file available at https://doi.org/10.5281/zenodo.5719520)
ID_MAPPING_FILE <<- paste0(MOD_PROJECT_PATH, '/_ignore/gwas_id_mapping.tsv')
