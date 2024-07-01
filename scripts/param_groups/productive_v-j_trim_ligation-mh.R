TRIM_TYPE <<- 'v-j_trim_ligation-mh'
PRODUCTIVITY <<- 'productive'
MOTIF_TYPE <<- 'unbounded'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 3), '_gene')
UPPER_TRIM_BOUND <<- 14 
LOWER_TRIM_BOUND <<- -2
INSERTIONS <<- 'zero'
MODEL_GROUP <<- 'all_subjects'
GENE_WEIGHT_TYPE <<- 'p_gene_pooled'
SAMPLE_ANNOT <<- FALSE
ONLY_NONPROD_SITES <<- TRUE
