import pygor3 as p3
import sys
import pandas as pd

output_directory = sys.argv[1]
input_directory = sys.argv[2]
iteration = sys.argv[3]

# seq_count = int(sys.argv[2])

# get default model
default_human_model = p3.get_default_IgorModel("human", "tcr_beta")

# simulate sequences
path = input_directory + '/bar_generated/generated_realizations_werr.csv'
df_scenarios = default_human_model.get_dataframe_from_fln_generated_realizations_werr(path)
# df_scenarios = p3.generate(Nseqs=seq_count, mdl=default_human_model, return_scenarios=True)

# retrieve actual junction feature quantities
junction_variables = {'v_trim':'v_3_del', 'd0_trim':'d_5_del', 'd1_trim':'d_3_del', 'j_trim':'j_5_del', 'dj_insert':'dj_ins', 'vd_insert':'vd_ins'}

for var in junction_variables.keys():
    df_scenarios[var] = default_human_model.get_realization_value_from_df_scenarios(df_scenarios, junction_variables[var])

# retrieve actual gene feature values
def get_vgene(scenario):
    v_choice = default_human_model.realization(scenario, 'v_choice')
    return v_choice

def get_dgene(scenario):
    d_choice = default_human_model.realization(scenario, 'd_gene')
    return d_choice

def get_jgene(scenario):
    j_choice = default_human_model.realization(scenario, 'j_choice')
    return j_choice

df_scenarios['v_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_vgene, df_scenarios = df_scenarios)
df_scenarios['v_gene'] = df_scenarios['v_gene_call'].astype(str).str.split('|').str[1]

df_scenarios['d_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_dgene, df_scenarios = df_scenarios)
df_scenarios['d_gene'] = df_scenarios['d_gene_call'].astype(str).str.split(';').str[0]

df_scenarios['j_gene_call'] = default_human_model.get_observable_from_df_scenarios(get_jgene, df_scenarios = df_scenarios)
df_scenarios['j_gene'] = df_scenarios['j_gene_call'].astype(str).str.split('|').str[1]

vd_nucs = default_human_model.get_df_realizations_dinucl(df_scenarios, 'vd_dinucl')
df_scenarios['vd_insert_nucs'] = vd_nucs.value.apply(lambda x: "".join(x))

dj_nucs = default_human_model.get_df_realizations_dinucl(df_scenarios, 'dj_dinucl')
df_scenarios['dj_insert_nucs'] = dj_nucs.value.apply(lambda x: "".join(x))

df_scenarios['productive'] = 'nonproductive'
df_scenarios['sample_name'] = 'igor_simulation'
cols = ['sample_name', 'v_gene', 'd_gene', 'j_gene', 'productive'] + list(junction_variables.keys())

file_name = output_directory + "/parsed_sampled_output_TRB_" + iteration + ".csv"

df_scenarios[cols].to_csv(file_name)
