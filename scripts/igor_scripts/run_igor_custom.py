import pygor3 as p3
import sys
import pandas as pd
import glob

OUTPUT_DIR = sys.argv[1]
FINAL_OUTPUT_DIR = sys.argv[2]
ANNOTATION_COUNT = int(sys.argv[3])
NCPU = int(sys.argv[4])
PARMS_PATH = sys.argv[5]
MARG_PATH = sys.argv[6]

# load nonproductive sequences and model
seqs = pd.read_csv(OUTPUT_DIR + '/cdr3_seqs.txt', sep = '\t', header = None)
new_model = p3.IgorModel.load_from_txt(PARMS_PATH, MARG_PATH)
print('finished loading sequences and model')

# align and annotate sequences
df_scens = p3.evaluate(seqs, new_model, N_scenarios=ANNOTATION_COUNT, igor_wd=OUTPUT_DIR, batch_clean=True, igor_evaluate_dict_opts = {'-threads':{'value':NCPU, 'active':True}}, igor_align_dict_opts={'--all':{'value':None, 'active':True}, '-threads':{'value':NCPU, 'active':True}})
print('finished aligning and annotating sequences')

# sample sequences
df_scens['seq_index'] = df_scens.index
df_scens.index.name = None
sampled = df_scens.groupby('seq_index').sample(n= 1, weights = df_scens['scenario_proba_cond_seq'])
print('finished sampling annotations')

# process outputs
## retrieve actual junction feature quantities
junction_variables = {'v_trim':'v_3_del', 'j_trim':'j_5_del', 'vj_insert':'vj_ins'}

for var in junction_variables.keys():
    sampled[var] = new_model.get_realization_value_from_df_scenarios(sampled, junction_variables[var])

## retrieve actual gene feature values
def get_vgene(scenario):
    v_choice = new_model.realization(scenario, 'v_choice')
    return v_choice

def get_jgene(scenario):
    j_choice = new_model.realization(scenario, 'j_choice')
    return j_choice

sampled['v_gene'] = new_model.get_observable_from_df_scenarios(get_vgene, df_scenarios = sampled)

sampled['j_gene'] = new_model.get_observable_from_df_scenarios(get_jgene, df_scenarios = sampled)
print('finished processing outputs')

# reformat
sampled = sampled.rename(columns = {'vj_dinucl':'vj_insert_nucs'})
cols = ['seq_index', 'scenario_rank', 'scenario_proba_cond_seq', 'v_gene', 'j_gene', 'v_trim', 'j_trim', 'vj_insert', 'vj_insert_nucs']
sampled = sampled[cols]

sampled['v_gene'] = sampled['v_gene'].astype(str).str.split(pat=';').str[0]
sampled['j_gene'] = sampled['j_gene'].astype(str).str.split(pat=';').str[0]
sampled['productive'] = 'FALSE'
print('finished reformatting outputs')

# write results
subject = OUTPUT_DIR.split('/')[-1]
name = FINAL_OUTPUT_DIR + '/' + subject + '_igor_sampled_annotations.tsv'
sampled.to_csv(name, sep = '\t', index = False)
print('finished writing results')
