TRIM_TYPE <<- 'v_trim'
PRODUCTIVITY <<- 'nonproductive'
MOTIF_TYPE <<- 'unbounded'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
UPPER_TRIM_BOUND <<- 14 
LOWER_TRIM_BOUND <<- -2
INSERTIONS <<- 'all'
MODEL_GROUP <<- 'all_subjects'
GENE_WEIGHT_TYPE <<- 'p_gene_pooled'
SAMPLE_ANNOT <<- TRUE
ONLY_NONPROD_SITES <<- TRUE
