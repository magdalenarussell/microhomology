import pygor3 as p3
import sys
import pandas as pd
import numpy as np

MOD_OUTPUT_PATH=sys.argv[1]

mdl_hb = p3.get_default_IgorModel("human", "tcr_alpha")

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_v_trim_params.tsv'

vchoice = {'v_gene':mdl_hb['v_choice'].lbl__v_choice.values, 'v_gene_prob':mdl_hb['v_choice'].values}
v_df = pd.DataFrame(vchoice)

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_vchoice_params.tsv'
v_df.to_csv(path, sep='\t', index=False)

vchoice_reshaped = mdl_hb['v_choice'].values[:, np.newaxis]
jprobs = vchoice_reshaped*mdl_hb['j_choice'].values
jprobs = jprobs.sum(axis=0)

jchoice = {'j_gene':mdl_hb['j_choice'].lbl__j_choice.values, 'j_gene_prob':jprobs}
j_df = pd.DataFrame(jchoice)
prob_sum = j_df.j_gene_prob.sum()
j_df['j_gene_prob'] = j_df.j_gene_prob/prob_sum

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_jchoice_params.tsv'
j_df.to_csv(path, sep='\t', index=False)

