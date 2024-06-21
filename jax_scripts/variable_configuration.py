import os
import importlib

class global_paramaters():
    def __init__(self, root_path, project_path, annotation_type, param_group, left_motif_size, right_motif_size, model_type):
        self.root_path = root_path
        self.project_path = project_path
        self.annotation_type = annotation_type
        self.param_group = param_group
        self.left_nuc_motif_count = left_motif_size
        self.right_nuc_motif_count = right_motif_size
        self.model_type = model_type.replace('twostep_', '')
        self.param_config = self.import_param_config()
        self.trim_type = getattr(self.param_config, "TRIM_TYPE")
        self.productivity = getattr(self.param_config, "PRODUCTIVITY")
        self.motif_type = getattr(self.param_config, "MOTIF_TYPE")
        self.gene_name = getattr(self.param_config, "GENE_NAME")
        self.upper_trim_bound = getattr(self.param_config, "UPPER_TRIM_BOUND")
        self.lower_trim_bound = getattr(self.param_config, "LOWER_TRIM_BOUND")
        self.insertions = getattr(self.param_config, "INSERTIONS")
        self.model_group = getattr(self.param_config, "MODEL_GROUP")
        self.gene_weight_type = getattr(self.param_config, "GENE_WEIGHT_TYPE")
        self.annotation_config = self.import_annotation_config()
        self.chain_type = getattr(self.annotation_config, "CHAIN_TYPE")
        self.junction_type = getattr(self.annotation_config, "JUNCTION_TYPE")
        self.sub_junction_type = getattr(self.annotation_config, "SUB_JUNCTION_TYPE")
        self.trimming_ligation_reannotated = getattr(self.annotation_config, "TRIMMING_LIGATION_REANNOTATED")
        self.sample_annotation = getattr(self.param_config, 'SAMPLE_ANNOT')
        self.only_nonprod_sites = getattr(self.param_config, 'ONLY_NONPROD_SITES')

    def import_param_config(self):
        param_config = importlib.import_module(f"param_group_configs.{self.param_group}")
        return(param_config)

    def import_annotation_config(self):
        if '_from_' in self.annotation_type:
            root_annotation = self.annotation_type.split('_from_')[0]
        else:
            root_annotation = self.annotation_type

        annotation_config = importlib.import_module(f"annotation_configs.{root_annotation}")
        return(annotation_config)

    def R_processed_data_path(self, sample_annotation = None, annotation = None):
        if annotation == None:
            annotation = self.annotation_type
        if sample_annotation == None:
            sample_annotation = self.sample_annotation

        path = self.root_path + '/' + annotation + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        if sample_annotation:
            file_name = path + '/processed_data.tsv'
        else:
            file_name = path + '/processed_data_all_annotations.tsv'
        return(file_name)

    def R_bootstrap_data_path(self, iteration, sample_annotation = None, annotation = None):
        if annotation == None:
            annotation = self.annotation_type
        if sample_annotation == None:
            sample_annotation = self.sample_annotation

        path = self.root_path + '/' + annotation + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/bootstraps'
        if sample_annotation:
            file_name = path + '/processed_data_bootstrap_' + str(iteration) + '.tsv'
        else:
            file_name = path + '/processed_data_all_annotations_bootstrap_' + str(iteration) + '.tsv'
        return(file_name)


    def R_input_domain_data_path(self):
        path = self.root_path + '/meta_data/' + self.chain_type + '/' + self.sub_junction_type
        file_name = path + '/frame_data.tsv'
        return(file_name)

    def R_subsampling_processed_data_path(self, prop, annotation = None):
        if annotation == None:
            annotation = self.annotation_type
        path = self.root_path + '/' + annotation + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type +'/temp_subsampling_exp/prop' + prop
        return(path)

    def model_output_path(self, l2):
        path = self.project_path + '/trained_models/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        os.makedirs(path, exist_ok=True)
        if self.sample_annotation:
            file_name = path + '/trained_model_L2' + str(l2) + '.pkl'
        else:
            file_name = path + '/trained_model_all_annotations_L2' + str(l2)+ '.pkl'
        return(file_name)

    def predictions_data_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        if self.sample_annotation:
            file_name = path + '/predicted_dist_data_L2' + str(l2) + '.tsv'
        else:
            file_name = path + '/predicted_dist_data_all_annotations_L2' + str(l2) + '.tsv'
        return(file_name)

    def validation_predictions_data_path(self, l2, validation_annotation):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/' + str(validation_annotation)
        os.makedirs(path, exist_ok=True)
        if self.sample_annotation:
            file_name = path + '/predicted_dist_data_L2' + str(l2) + '.tsv'
        else:
            file_name = path + '/predicted_dist_data_all_annotations_L2' + str(l2) + '.tsv'
        return(file_name)

    def trained_coefs_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        if self.sample_annotation:
            file_name = path + '/trained_coefs_L2' + str(l2) + '.tsv'
        else:
            file_name = path + '/trained_coefs_all_annotations_L2' + str(l2) + '.tsv'
        return(file_name)

    def trained_bootstrap_coefs_path(self, iteration, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/bootstraps'
        if self.sample_annotation:
            file_name = path + '/trained_coefs_L2' + str(l2) + '_bootstrap_' + str(iteration) + '.tsv'
        else:
            file_name = path + '/trained_coefs_all_annotations_L2' + str(l2) + '_bootstrap_' + str(iteration) + '.tsv'
        return(file_name)


    def subsampling_coefs_path(self, prop, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type + '/subsampling_experiment'
        os.makedirs(path, exist_ok=True)
        if self.sample_annotation:
            file_name = path + '/data_prop' + str(prop) + '_coefs_L2' + str(l2) + '.tsv'
        else:
            file_name = path + '/data_prop' + str(prop) + '_coefs_all_annotations_L2' + str(l2) + '.tsv'
        return(file_name)

    def model_eval_results_path(self, l2):
        path = self.root_path + '/' + self.annotation_type + '/' + self.param_group + '/' + self.motif_type + '_motif_trims_bounded_' + str(self.lower_trim_bound) + '_' + str(self.upper_trim_bound) + '/' + str(self.left_nuc_motif_count) + '_' + str(self.right_nuc_motif_count) + '_' + self.model_type
        if self.sample_annotation:
            file_name = path + '/model_evaluation_results_L2' + str(l2) + '.tsv'
        else:
            file_name = path + '/model_evaluation_results_all_annotations_L2' + str(l2) + '.tsv'
        return(file_name)


class model_specific_parameters():
    def __init__(self, param_group, model_type, left_motif_size, right_motif_size):
        self.param_group = param_group
        self.model_type = model_type
        self.left_nuc_motif_count = left_motif_size
        self.right_nuc_motif_count = right_motif_size
        self.param_config = self.import_param_config()
        self.model_type_config = self.import_model_type_config()
        self.variable_colnames = getattr(self.model_type_config, "VARIABLE_COLNAMES")
        self.choice1_variable_colnames = getattr(self.model_type_config, "CHOICE1_VARIABLE_COLNAMES")
        self.choice2_variable_colnames = getattr(self.model_type_config, "CHOICE2_VARIABLE_COLNAMES")
        self.count_colname = getattr(self.model_type_config, "COUNT_COLNAME")
        self.group_colname = getattr(self.param_config, "GENE_NAME")
        self.repeat_obs_colname = getattr(self.param_config, "REPEAT_OBS_COLNAME")
        self.choice_colname = getattr(self.param_config, "TRIM_TYPE")
        self.choice2_colname = None
        self.twostep = getattr(self.model_type_config, "TWOSTEP")

    def import_param_config(self):
        param_config = importlib.import_module(f"param_group_configs.{self.param_group}")
        return(param_config)

    def import_model_type_config(self):
        model_config = importlib.import_module(f"model_type_configs.{self.model_type}")
        return(model_config)

    def process_choice_colnames(self):
        trims = self.choice_colname.replace('_ligation-mh', '')
        trims = trims.split('_')[0]
        trims = trims.split('-')
        if 'adjusted_mh' in self.choice_colname:
            trims = [item + '_trim_adjusted_mh' for item in trims]
        else:
            trims = [item + '_trim' for item in trims]

        if self.twostep:
            self.choice_colname = trims
            self.choice2_colname = 'ligation_mh'
        else:
            if 'ligation-mh' in self.choice_colname:
                trims.append('ligation_mh')
            self.choice_colname = trims
        return(self.choice_colname, self.choice2_colname)

    def process_group_colnames(self):
        self.group_colname = self.group_colname.split('_')[0]
        self.group_colname = self.group_colname.split('-')
        self.group_colname = [item + '_gene' for item in self.group_colname]
        return(self.group_colname)

    def get_all_mh_variables(self, overlap_list=[0,1,2,3,4], prop=True, pos=['up', 'mid', 'down']):
        choice1 = False
        choice2 = False

        # Define the variables to check based on the value of 'prop'
        if prop:
            variables_to_check = ['mh', 'interior_mh']
        else:
            variables_to_check = ['mh_count', 'interior_mh_count']

        # Loop through each variable in the list of variables to check
        for variable in variables_to_check:
            # Check if the variable is in the main variable_colnames list
            if variable in self.variable_colnames:
                # If the variable is found, remove it from the relevant lists
                self.variable_colnames.remove(variable)

                # Now, check if the variable is in the choice1 or choice variable lists
                if variable in self.choice1_variable_colnames:
                    choice1 = True
                    self.choice1_variable_colnames.remove(variable)
                if variable in self.choice2_variable_colnames:
                    choice2 = True
                    self.choice2_variable_colnames.remove(variable)

                # Once the variable is found and processed, break out of the loop
                # since no need to check the next variable in variables_to_check
                break

        mh_vars = []
        for o in overlap_list:
            for p in pos:
                if o == 0 and p == 'mid':
                    continue
                if prop is True:
                    var = f'mh_prop_{p}_overlap_{o}'
                else:
                    var = f'mh_count_{p}_overlap_{o}'
                mh_vars.append(var)
        self.variable_colnames = self.variable_colnames + mh_vars
        if choice1:
            self.choice1_variable_colnames = self.choice1_variable_colnames + mh_vars
        if choice2:
            self.choice2_variable_colnames = self.choice2_variable_colnames + mh_vars
        return(self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames)

    def get_all_base_variables(self, side):
        assert side + '_base_count' in self.variable_colnames, "base_count is not a variable colname"
        variable= side + '_base_count'
        self.variable_colnames.remove(variable)
        choice1=False
        choice2=False
        if variable in self.choice1_variable_colnames:
            choice1 = True
            self.choice1_variable_colnames.remove(variable)
        if variable in self.choice2_variable_colnames:
            choice2 = True
            self.choice2_variable_colnames.remove(variable)

        bases = ['AT', 'GC']
        trims = [t for t in self.choice_colname if 'trim' in t]
        base_vars = [trim + '_' + side + '_base_count_' + base for trim in trims for base in bases if base + side != 'AT5end']
        self.variable_colnames = self.variable_colnames + base_vars
        if choice1:
            self.choice1_variable_colnames = self.choice1_variable_colnames + base_vars
        if choice2:
            self.choice2_variable_colnames = self.choice2_variable_colnames + base_vars
        return(self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames)

    def get_all_motif_variables(self):
        assert 'motif' in self.variable_colnames, "motif is not a variable colname"
        assert self.left_nuc_motif_count > 0, "left motif size must be greater than zero"
        assert self.right_nuc_motif_count > 0, "right motif size must be greater than zero"
        self.variable_colnames.remove('motif')
        choice1=False
        choice2=False
        if 'motif' in self.choice1_variable_colnames:
            choice1 = True
            self.choice1_variable_colnames.remove('motif')
        if 'motif' in self.choice2_variable_colnames:
            choice2 = True
            self.choice2_variable_colnames.remove('motif')

        trims = [t for t in self.choice_colname if 'trim' in t]
        motif_vars_5 = [trim + '_motif_5end_pos' + str(pos) for trim in trims for pos in range(self.left_nuc_motif_count, 0, -1)]
        motif_vars_3 = [trim + '_motif_3end_pos' + str(pos) for trim in trims for pos in range(1, self.right_nuc_motif_count+1)]
        self.variable_colnames = self.variable_colnames + motif_vars_5 + motif_vars_3
        if choice1:
            self.choice1_variable_colnames = self.choice1_variable_colnames + motif_vars_5 + motif_vars_3
        if choice2:
            self.choice2_variable_colnames = self.choice2_variable_colnames + motif_vars_5 + motif_vars_3
        return(self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames)

    def get_ligation_mh_variables(self):
        assert 'ligation_mh' in self.variable_colnames, "ligation_mh is not a varaible colname"
        return(self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames)

    def get_mh_config_count_variables(self):
        assert 'mh_config_count' in self.variable_colnames, "mh_config_count is not a varaible colname"
        return(self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames)

    def process_model_parameters(self):
        self.choice_colname, self.choice2_colname = self.process_choice_colnames()
        self.group_colname = self.process_group_colnames()
        if 'mh' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_mh_variables()
        if 'interior_mh' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_mh_variables(overlap_list=[1, 2, 3, 4], pos = ['mid'])
        if 'mh_count' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_mh_variables(prop=False)
        if 'interior_mh_count' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_mh_variables(overlap_list=[1, 2, 3, 4], prop=False, pos = ['mid'])
        if '5end_base_count' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_base_variables('5end')

        if '3end_base_count' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_base_variables('3end')

        if 'motif' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_all_motif_variables()

        if 'ligation_mh' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_ligation_mh_variables()

        if 'mh_config_count' in self.variable_colnames:
            self.variable_colnames, self.choice1_variable_colnames, self.choice2_variable_colnames = self.get_mh_config_count_variables()

        return(self)
