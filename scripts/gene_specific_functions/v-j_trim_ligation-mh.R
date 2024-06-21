stopifnot(LOCUS %in% c('alpha', 'gamma'))
stopifnot(INSERTIONS == 'zero')

source(paste0(MOD_PROJECT_PATH, '/scripts/gene_specific_functions/junction_specific_functions/vj.R'))
