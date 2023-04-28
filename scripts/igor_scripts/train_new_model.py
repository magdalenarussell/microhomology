import pygor3 as p3
import sys
import pandas as pd
import glob

input_directory = sys.argv[1]
final_output_directory = sys.argv[2]
sample_count = sys.argv[3]

files = glob.glob(input_directory + "/*.tsv")

# Read each TSV file into DataFrame
# This creates a list of dataframes
df_list = (pd.read_csv(file, sep = '\t') for file in files)

# Concatenate all DataFrames
all_seqs = pd.concat(df_list, ignore_index=True)

# sample sequences to use for model fitting
nonprod = all_seqs[all_seqs['productive'] == False]
sampled_nonprod = nonprod.sample(n = int(sample_count))
sampled_nonprod_seqs = sampled_nonprod['sequence']

# Get Genomic Germline templates from IMGT
imgt_species = 'Homo+sapiens'
imgt_chain = 'TRA'
hb_genomic_dict = p3.imgt.download_ref_genome(imgt_species, imgt_chain, dropna=True)

# create new model
hb_mdl_0 = p3.IgorModel.make_default_from_Dataframe_dict(hb_genomic_dict)

hb_mdl_new, df_likelihoods = p3.infer(sampled_nonprod_seqs, hb_mdl_0, N_iter=20, return_likelihoods=True)

marg_name = final_output_directory + '/model_' + str(sample_count) + '_marginals.txt'
parms_name = final_output_directory + '/model_' + str(sample_count) + '_parms.txt'

hb_mdl_new.write_model(fln_model_parms = parms_name, fln_model_marginals = marg_name)
