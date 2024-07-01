library(data.table)
source('config/config.R')

args = commandArgs(trailingOnly=TRUE)

input_data = args[1]
sample_count = args[2]
output_dir = args[3]

if (sample_count != 'all'){
    sample_count = as.numeric(sample_count)
}

df = fread(input_data)
df = df[frame_type == 'Out' | frame_stop == TRUE]

# get processed sequences
frame_file = file.path(MOD_OUTPUT_PATH, 'meta_data', 'TRA', 'VJ', 'frame_data.tsv')
frames = fread(frame_file)

df_tog = merge(df, frames, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh'))

df_subset = unique(df_tog[, c('index', 'processed_sequence', 'count')])

df_expand_indices = df_subset[!is.na(count) & count >0][,rep(.I, count)]
df_expand = df_subset[!is.na(count) & count >0][df_expand_indices]

if (sample_count == 'all'){
    sampled_dt = df_expand
} else {
    sampled_indices = sample(1:nrow(df_expand), sample_count, replace = FALSE)
    sampled_dt = df_expand[sampled_indices]
}

filename = file.path(output_dir, paste0('igor_sampled_sequences_', sample_count, '.txt'))
fwrite(data.table(sampled_dt$processed_sequence), filename)

cat(filename)
