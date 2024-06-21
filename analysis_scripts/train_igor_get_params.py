import pygor3 as p3
import sys
import pandas as pd
import numpy as np

MOD_OUTPUT_PATH=sys.argv[1]
print(MOD_OUTPUT_PATH)

test = pd.read_csv('_ignore/heikkila_alpha/adaptive/Thymus_1.tsv', sep = '\t')
hb_genomic_dict = p3.imgt.download_ref_genome('Homo+sapiens', 'TRA', dropna=True)
mdl_hb = p3.get_default_IgorModel("human", "tcr_alpha")

vgenes = pd.DataFrame()

for i in range(len(mdl_hb['v_3_del'])):
    item = mdl_hb['v_3_del'][i]
    temp_data = {'v_gene': item.lbl__v_choice.values,'v_trim': item.lbl__v_3_del.values,'v_trim_prob': item.values}
    temp = pd.DataFrame(temp_data)
    vgenes = pd.concat([vgenes, temp], ignore_index = True)

jgenes = pd.DataFrame()

for i in range(len(mdl_hb['j_5_del'])):
    item = mdl_hb['j_5_del'][i]
    temp_data = {'j_gene': item.lbl__j_choice.values,'j_trim': item.lbl__j_5_del.values,'j_trim_prob': item.values}
    temp = pd.DataFrame(temp_data)
    jgenes = pd.concat([jgenes, temp], ignore_index = True)

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_j_trim_params.tsv'
jgenes.to_csv(path, sep='\t', index=False)

path=MOD_OUTPUT_PATH + '/meta_data/igor_alpha_v_trim_params.tsv'
vgenes.to_csv(path, sep='\t', index=False)

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
