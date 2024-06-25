import pygor3 as p3
import sys
import pandas as pd
import glob

input_directory = sys.argv[1]
final_output_directory = sys.argv[2]
sample_count = sys.argv[3]
LOCUS = sys.argv[4]
assert LOCUS == 'alpha'

files = glob.glob(input_directory + "/*.txt")

# Read each TSV file into DataFrame
# This creates a list of dataframes
df_list = (pd.read_csv(file, sep = '\t', header = 0) for file in files)

# Concatenate all DataFrames
all_seqs = pd.concat(df_list, ignore_index=True)

# Get Genomic Germline templates from IMGT
imgt_species = 'Homo+sapiens'
imgt_chain = 'TRA'
hb_genomic = pd.read_csv('/home/mrussel2/microhomology/_ignore/imgt_tcra_sequences.tsv', sep = '\t')
# read anchors
v_anchors = pd.read_csv("https://raw.githubusercontent.com/statbiophys/SONIA/master/sonia/default_models/human_T_alpha/V_gene_CDR3_anchors.csv")
j_anchors = pd.read_csv("https://raw.githubusercontent.com/statbiophys/SONIA/master/sonia/default_models/human_T_alpha/J_gene_CDR3_anchors.csv")

# merge
hb_genomic_v = hb_genomic.drop(columns = ['names']).merge(v_anchors, on = 'gene').reset_index(drop = True)
hb_genomic_j = hb_genomic.drop(columns = ['names']).merge(j_anchors, on = 'gene').reset_index(drop = True)

# reformat
hb_genomic_v.columns = ['value', 'name', 'anchor_index', 'gfunction']
hb_genomic_j.columns = ['value', 'name', 'anchor_index', 'gfunction']
hb_genomic_dict = {'V':hb_genomic_v, 'J':hb_genomic_j}
# hb_genomic_dict = p3.imgt.download_ref_genome(imgt_species, imgt_chain, dropna=True)

# create new model
hb_mdl_0 = p3.IgorModel.make_default_from_Dataframe_dict(hb_genomic_dict)

# get alignments
df_functionality, df_CDR3 = p3.naive_align(all_seqs['V1'], hb_mdl_0)

hb_mdl_new, df_likelihoods = p3.infer(all_seqs, hb_mdl_0, N_iter=10, return_likelihoods=True)

marg_name = final_output_directory + '/model_' + str(sample_count) + '_marginals.txt'
parms_name = final_output_directory + '/model_' + str(sample_count) + '_parms.txt'

hb_mdl_new.write_model(fln_model_parms = parms_name, fln_model_marginals = marg_name)
