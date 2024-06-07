import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import dill
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_twostep_model_classes import TwoStepDataTransformer, TwoStepConditionalLogisticRegressor

class TwoStepDataTransformerEM(TwoStepDataTransformer):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params)
        self.choice2_colname = choice2_colname
        self.choice1_variable_colnames = choice1_variable_colnames
        self.choice2_variable_colnames = choice2_variable_colnames
        self.original_choice1_variable_colnames = choice1_variable_colnames
        self.original_choice2_variable_colnames = choice2_variable_colnames
        self.original_choice2_colname = choice2_colname
        self.input_choice1_variable_colnames = choice1_variable_colnames
        self.input_choice2_variable_colnames = choice2_variable_colnames
        self.input_choice2_colname = choice2_colname
        self.coefs = None

    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares and returns data matrices necessary for the two-step modeling approach from the preprocessed DataFrame.

        Args:
            df (pd.DataFrame): The DataFrame from which matrices are to be generated, typically after preprocessing.
            pretrain (bool, optional): Indicates if the matrix preparation is for pretraining purposes. Defaults to True.
            replace_object (str, optional): Name of the attribute to replace with the preprocessed DataFrame. Defaults to None.
            return_df (bool, optional): Flag indicating whether to return the preprocessed DataFrame along with the matrices. Defaults to False.

        Returns:
            tuple: Contains matrices for choice1 variables, choice2 variables, counts, non-repeat groups, and all_site_mask, along with the optionally returned preprocessed DataFrame.
        """
        df = self.preprocess_data(df, pretrain)
        df['seq_index'] = df['index']

        # Fill in missing counts with zero
        df[self.count_colname] = df[self.count_colname].fillna(0).astype(float)

        # Get matrix shapes
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])
        choice2s = pd.unique(df[self.choice2_colname])
        indices = pd.unique(df.seq_index)

        choice1_final_shape = (len(groups), len(choices), len(self.choice1_variable_colnames))
        choice1_int_shape = (len(groups), len(self.choice1_variable_colnames), len(choices))
        choice2_final_shape = (len(groups), len(choices), len(choice2s), len(self.choice2_variable_colnames))
        choice2_int_shape = (len(groups), len(self.choice2_variable_colnames), len(choices), len(choice2s))
        counts_shape = (len(groups), len(choices), len(choice2s), 1)
        indices_shape = (len(indices), len(groups), len(choices), len(choice2s), 1)
        prod_shape = (len(groups), len(choices), 1)

        # map variables to indices
        choice1_var_mapping = {key: value for key, value in enumerate(self.choice1_variable_colnames)}
        choice1_var_mapping_df = pd.DataFrame(list(choice1_var_mapping.items()), columns=['VarIndex', 'VarName'])
        choice2_var_mapping = {key: value for key, value in enumerate(self.choice2_variable_colnames)}
        choice2_var_mapping_df = pd.DataFrame(list(choice2_var_mapping.items()), columns=['VarIndex', 'VarName'])


        # fill in df with all combos
        all_combinations = pd.MultiIndex.from_product([choices, choice2s],
                                                      names=[self.choice_colname, self.choice2_colname]).to_frame(index=False)
        df_full = pd.merge(all_combinations, df, on=[self.choice_colname, self.choice2_colname], how='left').fillna(0)
        df_full.loc[df_full.prod_domain_indicator == 0, 'seq_index'] = -1.0

        # get counts matrix
        pivot_counts = df_full.pivot_table(index=[self.group_colname],
                                      columns=[self.choice_colname, self.choice2_colname],
                                      values=self.count_colname,
                                      fill_value=0)
        counts_mat = jnp.array(pivot_counts).reshape(counts_shape)

        # get variable matrix
        choice1_pivot_vars = df.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=choice1_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice1_pivot_vars = choice1_pivot_vars.reindex(choice1_var_mapping_df.VarName.tolist(), axis=1, level = 0)

        choice2_pivot_vars = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=choice2_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice2_pivot_vars = choice2_pivot_vars.reindex(choice2_var_mapping_df.VarName.tolist(), axis=1, level = 0)
        assert len(self.choice2_variable_colnames) == 1, "need to verify correct choice2 matrix"

        choice1_mat = jnp.array(choice1_pivot_vars).reshape(choice1_int_shape)
        choice1_mat = choice1_mat.transpose((0, 2, 1))
        assert choice1_mat.shape == choice1_final_shape, "choice1 variable matrix is the incorrect dimension"

        choice2_mat = jnp.array(choice2_pivot_vars).reshape(choice2_int_shape)
        choice2_mat = choice2_mat.transpose((0, 2, 3, 1))
        assert choice2_mat.shape == choice2_final_shape, "choice2 variable matrix is the incorrect dimension"

        # Assuming groups is a list of unique group identifiers
        group_indices = {group: idx for idx, group in enumerate(groups)}
        nonrepeat_groups_mat = jnp.array([group_indices[group] for group in df[self.group_colname].unique()])

        # get mask matrix (1 for valid entries or 0 for nonvalid)
        pivot_mask = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=['all_site_indicator'],
                                    fill_value=0)

        pivot_mask_prod = df_full.pivot_table(index=[self.group_colname],
                                              columns=[self.choice_colname, self.choice2_colname],
                                              values=['prod_domain_indicator'],
                                              fill_value=0)

        all_site_mask_mat = jnp.array(pivot_mask).reshape(counts_shape)
        prod_mask_mat = jnp.array(pivot_mask_prod).reshape(counts_shape)

        # get index matrix
        pivot_index = df_full.pivot_table(index=[self.group_colname],
                            columns=[self.choice_colname, self.choice2_colname],
                            values=['seq_index'],
                            fill_value=-1.0)
        index_mat = jnp.array(pivot_index).reshape(counts_shape)


        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, all_site_mask_mat, index_mat, prod_mask_mat, df
        else:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, all_site_mask_mat, index_mat, prod_mask_mat



class TwoStepConditionalLogisticRegressorEM(TwoStepDataTransformerEM, TwoStepConditionalLogisticRegressor):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params, l2kfold=10):
        super().__init__(training_df=None, variable_colnames=variable_colnames, choice1_variable_colnames=choice1_variable_colnames, choice2_variable_colnames=choice2_variable_colnames, count_colname=count_colname, group_colname=group_colname, repeat_obs_colname=repeat_obs_colname, choice_colname=choice_colname, choice2_colname=choice2_colname, params=params)
        self.training_df = training_df
        if training_df is not None:
            self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix, self.all_site_mask_matrix, self.index_matrix, self.prod_mask_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.original_choice1_variable_colnames = choice1_variable_colnames
        self.original_choice2_variable_colnames = choice2_variable_colnames
        self.original_choice2_colname = choice2_colname
        self.input_choice1_variable_colnames = choice1_variable_colnames
        self.input_choice2_variable_colnames = choice2_variable_colnames
        self.input_choice2_colname = choice2_colname
        self.initial_coefs = self.get_ones_coef_array()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-6
        self.step = 0.1
        self.maxiterEM = 100
        self.toleranceEM = 1e-8
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    def get_index_weights(self, choice1_variables, choice2_variables, all_site_mask, index_matrix, prod_mask_matrix, weighting_coefs):
        probs = self.get_joint_prob(choice1_variables, choice2_variables, all_site_mask, prod_mask_matrix, weighting_coefs)
        flatten_probs = probs.flatten()
        flatten_index = index_matrix.flatten().astype(int)

        total_probabilities = jnp.bincount(flatten_index, weights=flatten_probs)

        marginalized_probs = flatten_probs / total_probabilities[flatten_index]

        marginalized_probs_reshape = marginalized_probs.reshape(index_matrix.shape)
        return(marginalized_probs_reshape)

    # Get cross-entropy loss
    def cross_entropy(self):
        return None

    def l2regularization(self, coefs, size, l2reg):
        """
        Computes the L2 regularization term for the model coefficients.

        Args:
            coefs (ndarray): Coefficients of the model.
            size (int): The sample size.
            l2reg (float): L2 regularization strength.

        Returns:
            float: The L2 regularization term.
        """
        def calculate_coef_sum(choice_var_colnames, coefs):
            var_list = [var for var in choice_var_colnames if 'base_count' not in var or 'interaction' in var]
            var_list = [var for var in var_list if 'motif' not in var]

            mh_list = [var in var_list for var in choice_var_colnames]
            mh_jnp = jnp.array(mh_list).reshape(-1, 1)

            coef_subset = coefs * mh_jnp

            c = jnp.nansum(coef_subset**2)
            return c

        c = 0
        # add coefs from choice1 variables
        c += calculate_coef_sum(self.choice1_variable_colnames + self.choice2_variable_colnames, coefs)
        return(0.5*(1/size)*l2reg*c)

    def loss_fn(self, coefs, marg_weights, choice1_variables, choice2_variables, counts, all_site_mask, prod_mask_matrix, l2reg=0):
        """
        Computes the total loss for the two-step conditional logistic regression model, including the cross-entropy loss

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            choice1_variables (ndarray): Data matrix for the first choice variables.
            choice2_variables (ndarray): Data matrix for the second choice variables.
            counts (ndarray): Counts of choices.
            all_site_mask (ndarray): Data matrix indicating valid data entries.

        Returns:
            float: Total loss, including cross-entropy
        """
        choice_probs = self.get_joint_prob(choice1_variables, choice2_variables, all_site_mask, prod_mask_matrix, coefs)

        if choice_probs.ndim == 1:
            choice_probs = choice_probs.reshape((choice_probs.shape[0], 1))

        size = jnp.count_nonzero(prod_mask_matrix).item()

        counts_reshape = jnp.squeeze(counts)
        marg_weights_reshape = jnp.squeeze(marg_weights)

        loss = -jnp.nansum(jnp.log(jnp.where(choice_probs==0, 1, choice_probs)) * counts_reshape * marg_weights_reshape) + self.l2regularization(coefs, size, l2reg)
        return loss

    def m_step(self, choice1_variable_matrix, choice2_variable_matrix, counts_matrix, all_site_mask_matrix, prod_mask_matrix, l2reg, maxiter, tol, step, initial_coefs, marg_weights):
        """
        Fits the two-step conditional logistic regression model to the data using an optimization algorithm.

        Args:
            choice1_variable_matrix (ndarray), choice2_variable_matrix (ndarray): Data matrices for the first and second choice variables.
            counts_matrix (ndarray): Counts of choices.
            all_site_mask_matrix (ndarray): Mask indicating valid data entries.
            maxiter (int): Maximum number of optimization iterations.
            tol (float): Convergence tolerance.
            step (float): Step size for the optimization.
            initial_coefs (ndarray): Initial coefficients for the model fitting.

        Returns:
            The result of the optimization process, including the optimized coefficients.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.BFGS(fun=self.loss_fn, maxiter=maxiter, tol=tol, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         marg_weights=marg_weights,
                         choice1_variables=choice1_variable_matrix,
                         choice2_variables=choice2_variable_matrix,
                         counts=counts_matrix,
                         all_site_mask=all_site_mask_matrix,
                         prod_mask_matrix=prod_mask_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg, starting_coefs):
        """
        Computes the cross-validation loss for the given L2 regularization strength across the specified number of folds.

        Args:
            fold_count (int): The number of folds to use for cross-validation.
            l2reg (float): The L2 regularization strength to use when fitting the model on each fold.

        Returns:
            list of float: A list containing the cross-validation loss for each fold. This provides an estimate of
                           the model's performance on unseen data when trained with the specified regularization strength.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.choice1_variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_c1_data, train_c2_data, val_c1_data, val_c2_data = self.choice1_variable_matrix[train_index], self.choice2_variable_matrix[train_index], self.choice1_variable_matrix[val_index], self.choice2_variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_all_site_mask, val_all_site_mask = self.all_site_mask_matrix[train_index], self.all_site_mask_matrix[val_index]
            train_pp, val_pp = self.prod_mask_matrix[train_index], self.prod_mask_matrix[val_index]
            train_index, val_index = self.index_matrix[train_index], self.index_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            marg_weights = self.get_index_weights(train_c1_data, train_c2_data, train_all_site_mask, train_index, train_pp, starting_coefs)

            model = self.m_step(train_c1_data, train_c2_data, train_counts, train_all_site_mask, train_pp, l2reg, self.maxiter, self.tolerance, self.step, starting_coefs, marg_weights)

            # Compute the loss on the validation data
            val_marg_weights = self.get_index_weights(val_c1_data, val_c2_data, val_all_site_mask, val_index, val_pp, model.params)

            loss = self.loss_fn(model.params,
                                val_marg_weights,
                                val_c1_data,
                                val_c2_data,
                                val_counts,
                                val_all_site_mask,
                                val_pp)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, starting_coefs, l2kfold, l2reg_values=[10**i for i in range(-6, 2)] + [0]):
        """
        Perform grid search cross-validation for hyperparameter tuning.

        Args:
            l2kfold (int): Number of folds for cross-validation.
            l2reg_values (list): List of L2 regularization strengths to search.

        Returns:
            pd.DataFrame: DataFrame with hyperparameter values and corresponding mean cross-validation loss.
        """
        results = pd.DataFrame({'l2reg':[], 'mean_CV_loss':[]})

        for l2reg in l2reg_values:
            # Perform cross-validation
            scores = self.cv_loss(l2kfold, l2reg, starting_coefs)

            # Calculate the mean cross-validation score
            mean_score = np.mean(scores)

            temp = pd.DataFrame({'l2reg':[l2reg], 'mean_CV_loss':[mean_score]})
            results = pd.concat([results, temp], ignore_index=True)

        return(results)

    def get_l2reg(self, starting_coefs):
        """
        Determine the optimal L2 regularization strength.

        Returns:
            float: Optimal L2 regularization strength.
        """
        self.l2kfold = min(len(np.unique(self.nonrepeat_grp_matrix)), 10)
        grid = self.grid_search_cv(starting_coefs, self.l2kfold)
        self.l2reg_grid = grid
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, l2reg_value=None, maxiter=None, tolerance=None, step=None, maxiterEM=None, toleranceEM=None):
        assert self.counts_matrix is not None, "counts column is needed"

        if maxiter is None:
            maxiter = self.maxiter
        else:
            self.maxiter = maxiter

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        if step is None:
            step = self.step
        else:
            self.step = step

        if l2 is False:
            self.l2reg = 0
        else:
            if l2reg_value is None:
                self.l2reg = self.get_l2reg(self.initial_coefs)
            else:
                self.l2reg = l2reg_value

        if maxiterEM is None:
            maxiterEM = self.maxiterEM
        else:
            self.maxiterEM = maxiterEM

        if toleranceEM is None:
            toleranceEM = self.toleranceEM
        else:
            self.toleranceEM = toleranceEM

        coefs = self.initial_coefs
        iteration = 0
        EMlossTol = 100
        last_loss = 100
        EM_training_info_log = {}

        while iteration < maxiterEM and EMlossTol > toleranceEM:
            iteration += 1

            marg_weights = self.get_index_weights(self.choice1_variable_matrix, self.choice2_variable_matrix, self.all_site_mask_matrix, self.index_matrix, self.prod_mask_matrix, coefs)

            res = self.m_step(self.choice1_variable_matrix,
                           self.choice2_variable_matrix,
                           self.counts_matrix,
                           self.all_site_mask_matrix,
                           self.prod_mask_matrix,
                           self.l2reg,
                           maxiter,
                           tolerance,
                           step,
                           coefs,
                           marg_weights)

            coefs = res.params
            EM_training_info_log[iteration] = res

            EMlossTol = abs(last_loss - float(res.state.value))
            last_loss = float(res.state.value)

        self.training_info_EM_log = EM_training_info_log
        self.EM_iteration_losses = [float(EM_training_info_log[i].state.value) for i in range(1, iteration)]
        self.EM_iteration_errors = [float(EM_training_info_log[i].state.error) for i in range(1, iteration)]
        self.EM_iterations = iteration
        self.EM_loss_tolerance = EMlossTol
        self.coefs = res.params
        self.training_info = res
        self.cov_matrix = self.get_cov_matrix(self.training_info)
        self.standard_errors = self.get_errors(self.training_info)
        return self

    def get_cov_matrix(self, training_info):
        hess_mat = training_info.state.H
        cov_matrix = jnp.linalg.inv(hess_mat)
        return cov_matrix

    def get_errors(self, training_info):
        cov = self.get_cov_matrix(training_info)
        standard_errors = np.sqrt(np.diag(cov))
        return standard_errors

    def save_model(self, file_path):
        """
        Saves the trained model to a specified file path for later use or analysis. This method serializes the model instance, including its coefficients and configuration.

        Args:
            file_path (str): The path where the model should be saved.

        Note:
            Before saving, large attributes such as the training and validation data matrices are set to None to minimize the file size.
        """
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.choice1_variable_matrix = None
        self.choice2_variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        self.all_site_mask_matrix = None
        self.index_matrix = None
        self.training_info_EM_log = None
        self.prod_mask_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)


class TwoStepConditionalLogisticRegressionPredictorEM(TwoStepDataTransformerEM):
    def __init__(self, model, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(training_df=None, variable_colnames=variable_colnames, choice1_variable_colnames=choice1_variable_colnames, choice2_variable_colnames=choice2_variable_colnames, count_colname=count_colname, group_colname=group_colname, repeat_obs_colname=repeat_obs_colname, choice_colname=choice_colname, choice2_colname=choice2_colname, params=params)
        self.original_choice1_variable_colnames = choice1_variable_colnames
        self.original_choice2_variable_colnames = choice2_variable_colnames
        self.original_choice2_colname = choice2_colname
        self.input_choice1_variable_colnames = choice1_variable_colnames
        self.input_choice2_variable_colnames = choice2_variable_colnames
        self.input_choice2_colname = choice2_colname
        self.model = model
        if not isinstance(model, TwoStepConditionalLogisticRegressorEM):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressorEM object")

    # get probability for input parameters given coefficients
    def predict(self, new_df):
        """
        Make predictions using the trained model.

        Returns:
            pd.DataFrame: A DataFrame containing predicted probabilities merged with the original data.
        Raises:
            ValueError: If the number of variable columns in the input DataFrame doesn't match the model's coefficients.
        """
        original_new_df = new_df
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not choice1_variable_matrix.shape[-1] + choice2_variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_joint_prob(choice1_variable_matrix, choice2_variable_matrix, all_site_mask_matrix, prod_mask_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice1_cols = self.get_mapping_dict(new_df, self.choice_colname)
        choice2_cols = self.get_mapping_dict(new_df, self.choice2_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice1_cols), len(choice2_cols)))

        # Convert the 3D array to a 2D array where each "page" is flattened out
        # and then stack them vertically
        array_2d = probs.reshape(-1, probs.shape[2])

        # Create a MultiIndex representing the first two dimensions (depth and rows)
        # np.repeat and np.tile help to properly align the indices with the reshaped array
        index = pd.MultiIndex.from_tuples([(i, j) for i in range(probs.shape[0]) for j in range(probs.shape[1])], names=[self.group_colname, self.choice_colname])

        # Create the DataFrame
        prob_df = pd.DataFrame(array_2d, index=index, columns=[i for i in range(probs.shape[2])])
        prob_df = prob_df.reset_index()

        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname, self.choice_colname],
                            var_name=self.choice2_colname,
                            value_name='predicted_prob')

        # Invert the dictionary
        choice1_mapping = {v: k for k, v in choice1_cols.items()}
        choice2_mapping = {v: k for k, v in choice2_cols.items()}
        group_mapping = {v: k for k, v in group_cols.items()}

        # Use the map function to replace values in 'Column1'
        melted_df[self.group_colname] = melted_df[self.group_colname].map(group_mapping)
        melted_df[self.choice_colname] = melted_df[self.choice_colname].map(choice1_mapping)
        melted_df[self.choice2_colname] = melted_df[self.choice2_colname].map(choice2_mapping)

        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname, self.choice2_colname],
                             how='inner')

        # fill in out of domain values
        common_cols = list(set(original_new_df) & set(merged_df))
        final_df = pd.merge(merged_df, original_new_df, how='outer', on=common_cols)

        # Fill missing values with a specified value
        final_df.fillna(0, inplace=True)
        return(final_df)

    def compute_loss(self, new_df):
        """
        Compute the loss for the given data.

        Returns:
            float: The computed loss.
        """
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        marg_weights = self.model.get_index_weights(choice1_variable_matrix, choice2_variable_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix, self.model.coefs)

        loss = self.model.loss_fn(self.model.coefs,
                                  marg_weights,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  all_site_mask_matrix,
                                  prod_mask_matrix)

        return(float(loss))

class TwoStepConditionalLogisticRegressionEvaluatorEM(TwoStepDataTransformerEM):
    """
    Evaluates a trained Two-Step Conditional Logistic Regression model by calculating log loss on both training and validation datasets. This class extends the TwoStepDataTransformer for consistent data preparation before evaluation and supports loading a saved model for evaluation purposes.

    Args:
        model_path (str): Path to the saved trained model file.
        params (dict): Parameters used for data transformation and evaluation.
        training_df (pd.DataFrame, optional): The DataFrame containing the training data.
        validation_df (pd.DataFrame, optional): The DataFrame containing the validation data for evaluation.

    Attributes:
        model (TwoStepConditionalLogisticRegressor): The trained Two-Step Conditional Logistic Regression model.
        log_loss (float, optional): Log loss on the training data, initialized as None.
        expected_log_loss (float, optional): Expected log loss calculated through cross-validation, initialized as None.
        validation_df (pd.DataFrame): DataFrame containing validation data.

    Methods:
        calculate_log_loss: Calculates log loss on the training dataset.
        calculate_expected_log_loss: Estimates log loss on the training dataset using cross-validation.
        calculate_validation_log_loss: Calculates log loss on the validation dataset.
        compile_evaluation_results_df: Compiles evaluation results into a DataFrame.
    """
    def __init__(self, model_path, params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, TwoStepConditionalLogisticRegressorEM):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressorEM object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.input_choice1_variable_colnames, self.model.input_choice2_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, self.model.input_choice2_colname, params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        """
        Loads a trained Two-Step Conditional Logistic Regression model from a specified file path.

        Args:
            file_path (str): The file path to the saved model.

        Returns:
            TwoStepConditionalLogisticRegressor: The loaded trained model.
        """
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        """
        Calculates the log loss on the training dataset using the loaded model's coefficients.

        Returns:
            float: The log loss value for the training dataset.
        """
        assert self.model.training_df is not None, 'No input training dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        marg_weights = self.model.get_index_weights(choice1_variable_matrix, choice2_variable_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix, self.model.coefs)

        loss = self.model.loss_fn(self.model.coefs,
                                  marg_weights,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  all_site_mask_matrix,
                                  prod_mask_matrix)

        return(loss)

    def calculate_validation_log_loss(self):
        """
        Calculates the log loss on the validation dataset using the loaded model's coefficients.

        Returns:
            float: The log loss value for the validation dataset.
        """
        assert self.validation_df is not None, 'No input validation dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        marg_weights = self.model.get_index_weights(choice1_variable_matrix, choice2_variable_matrix, all_site_mask_matrix, index_matrix, prod_mask_matrix, self.model.coefs)

        loss = self.model.loss_fn(self.model.coefs,
                                  marg_weights,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  all_site_mask_matrix,
                                  prod_mask_matrix)
        return(loss)

    def compile_evaluation_results_df(self, calculate_validation_loss = False):
        """
        Compiles the evaluation results, including training log loss, expected log loss, and validation log loss, into a DataFrame for easy comparison and analysis.

        Args:
            calculate_validation_loss (bool, optional): If True, calculates and includes the log loss on the validation dataset in the results. Default is False.

        Returns:
            pd.DataFrame: A DataFrame containing the compiled evaluation results and the model parameters used during training and evaluation.
        """
        result = {'training_annotation_type':[self.params.annotation_type],
                  'productivity':[self.params.productivity],
                  'motif_length_5_end':[self.params.left_nuc_motif_count],
                  'motif_length_3_end':[self.params.right_nuc_motif_count],
                  'motif_type':[self.params.motif_type],
                  'gene_weight_type':[self.params.gene_weight_type],
                  'upper_bound':[self.params.upper_trim_bound],
                  'lower_bound':[self.params.lower_trim_bound],
                  'insertion_bound':[self.params.insertions],
                  'model_type':[self.params.model_type],
                  'base_count_5end_length':[10],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        final = pd.DataFrame()

        if calculate_validation_loss is True:
            val_loss = self.calculate_validation_log_loss()
            val = results_df.copy()
            val['loss type'] = 'Log loss on validation data'
            val['log_loss'] = val_loss

            final = pd.concat([final, val], axis = 0)

        return(final)



