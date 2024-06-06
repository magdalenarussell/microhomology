get_igor_gene_usage_params <- function(){
    # get igor params
    j_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_jchoice_params.tsv')
    v_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_vchoice_params.tsv')

    jg = fread(j_choice_path)
    vg = fread(v_choice_path)
    return(list(v_choice = vg, j_choice = jg))
}
