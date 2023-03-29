# (required) the following paths indicate the location of the project (PROJECT_PATH) and the location of where files should be written (OUTPUT PATH)
OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/microhomology'
PROJECT_PATH = paste0('/home/', Sys.getenv('USER'), '/microhomology')

# (optional) Parsed TCR repertoire data 
TRA_REPERTOIRE_DATA = paste0(OUTPUT_PATH, '/alpha/thymus_alpha/adaptive')

# (optional) Other gene name and germline sequences for other loci
WHOLE_NUCSEQS_alpha = paste0(PROJECT_PATH,'/_ignore/imgt_tcra_sequences.tsv')
