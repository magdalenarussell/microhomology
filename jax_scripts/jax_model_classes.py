import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import os
os.environ["JAX_PLATFORM_NAME"] = "cpu"
jax.config.update('jax_default_device', jax.devices('cpu')[0])
jax.config.update("jax_disable_jit", True)
import patsy
import dill
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold

class DataPreprocessor():
    """
    A class for preprocessing training data for conditional logit models,
    specifically designed to handle transformations and encoding of categorical
    variables, expansion of multivariable columns, and filtering based on various criteria.

    Attributes:
        training_df (pd.DataFrame): DataFrame containing the training data.
        variable_colnames (list): List of column names to be used as variables in the model.
        choice1_variable_colnames (list): List of column names representing the first choice variables, if applicable.
        choice2_variable_colnames (list): List of column names representing the second choice variables, if applicable.
        count_colname (str): Name of the column containing count data.
        group_colname (str): Name of the column used to group the data.
        repeat_obs_colname (str): Name of the column indicating repeated observations.
        choice_colname (str): Name of the column containing choice data.
        params (dict): Dictionary containing parameters for preprocessing.
    """
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        self.training_df = training_df
        self.variable_colnames = variable_colnames
        self.choice1_variable_colnames = choice1_variable_colnames
        self.choice2_variable_colnames = choice2_variable_colnames
        self.count_colname = count_colname
        self.group_colname = group_colname
        self.repeat_obs_colname = repeat_obs_colname
        self.choice_colname = choice_colname
        self.params = params

    def get_mapping_dict(self, training_df, col):
        """
        Generates a mapping dictionary that maps each unique value in a specified column to a unique integer.

        Parameters:
            training_df (pd.DataFrame): DataFrame containing the data.
            col (str): Name of the column for which to generate the mapping.

        Returns:
            dict: A dictionary where keys are unique values from `col` and values are unique integers.
        """
        # Create a mapping dictionary for unique values in a column
        unique = sorted(list(pd.unique(training_df[col])))
        mapping = {key: index for index, key in enumerate(unique)}
        return mapping

    def transform_categorical_response_vars(self, training_df, col, new_col):
        """
        Transform a categorical variable into integers.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the column to be transformed.
            new_col (str): The name of the new integer-encoded column.

        Returns:
            pd.DataFrame: The DataFrame with the integer-encoded column added.
        """
        # Transform a categorical variable into integers
        col = getattr(self, col)
        if training_df[col].iloc[1] != int:
            new_col_name = col + '_int'
            setattr(self, new_col, new_col_name)
            var_mapping = self.get_mapping_dict(training_df, col)
            training_df[new_col_name] = training_df[col].map(var_mapping)
        return training_df

    def get_contrast_matrix(self, training_df, col, pretrain=True):
        """
        Get the contrast matrix for categorical variables.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the categorical column.

        Returns:
            pd.DataFrame: The contrast matrix for the categorical column.
        """
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            # create sum contrasts matrix
            unique = sorted(list(pd.unique(training_df[col])))
            if not pretrain:
                unique = ['-', 'A', 'C', 'G', 'T']

            contrast = Sum().code_without_intercept(unique).matrix

            # Create an identity matrix with the number of categories
            identity_matrix = np.eye(len(unique))

            # Apply the contrast matrix to the identity matrix
            encoded_values = np.dot(identity_matrix, contrast)

            # Create a new DataFrame with the encoded values
            encoded_data = pd.DataFrame(encoded_values, columns=[f"{col}_{i}" for i in unique[:-1]])
            encoded_data[col] = unique
            nonbase_var = col + '_-'
            if nonbase_var in encoded_data.columns:
                encoded_data[nonbase_var] = 0.0
                encoded_data.loc[encoded_data[col] == '-', encoded_data.columns != col] = 0.0
        else:
            encoded_data = pd.DataFrame({col:unique})
        return(encoded_data)

    def get_dropped_contrast_var(self, training_df, original_col):
        """
        Get the names of dropped contrast variables after transformation.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            original_col (str): The name of the original categorical column.

        Returns:
            tuple: A tuple containing two lists - the names of contrast variables and the name of the dropped contrast variable.
        """
        unique = sorted(list(pd.unique(training_df[original_col])))
        contrast_vars=[f"{original_col}_{i}" for i in unique]
        return(contrast_vars[:-1], contrast_vars[-1])

    def transform_categorical_vars(self, training_df, col, pretrain=True):
        """
        Transforms categorical variables in the DataFrame into contrast columns for model training.

        Parameters:
            training_df (pd.DataFrame): DataFrame containing the data.
            col (str): Name of the categorical column to be transformed.
            pretrain (bool, optional): If set to True, performs transformations assuming pretraining. Defaults to True.

        Returns:
            pd.DataFrame: Updated DataFrame with categorical variables transformed into contrast columns.
        """
        assert col in self.variable_colnames, "Input column name is not a variable name"
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            contrast_df = self.get_contrast_matrix(training_df, col, pretrain)
            training_df = pd.merge(training_df, contrast_df, on=col, how='inner')
            new_cols = [x for x in list(contrast_df.columns) if x != col]
            self.variable_colnames = [x for x in self.variable_colnames if x != col] + new_cols
            if self.choice1_variable_colnames is not None:
                if col in self.choice1_variable_colnames:
                    self.choice1_variable_colnames = [x for x in self.choice1_variable_colnames if x != col] + new_cols

            if self.choice2_variable_colnames is not None:
                if col in self.choice2_variable_colnames:
                    self.choice2_variable_colnames = [x for x in self.choice2_variable_colnames if x != col] + new_cols

        return(training_df)

    def expand_multivariable(self, training_df, col, new_col):
        """
        Expand multi-variable columns into a single string column.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the multi-variable column.
            new_col (str): The name of the new string column.

        Returns:
            pd.DataFrame: The DataFrame with the multi-variable column expanded into a string column.
        """
        # Expand multi-variable columns into a single string column
        if type(getattr(self, col)) == list:
            if len(getattr(self, col)) > 1:
                new_col_name = '_'.join(getattr(self, col))
                setattr(self, new_col, new_col_name)
                training_df[new_col_name] = training_df.parallel_apply(lambda row: '_'.join(row[getattr(self, col)].astype(str)), axis=1)
            else:
                new_col_name = getattr(self, col)[0]
            setattr(self, col, new_col_name)
        return training_df

    def check_within_set_variance(self, training_df, pretrain=True):
        """
        Checks and removes variables with insufficient within-set variance, modifying the class attribute lists of variable names.

        Parameters:
            training_df (pd.DataFrame): DataFrame containing the data.
            pretrain (bool, optional): Determines if the check is performed under the assumption of pretraining. Defaults to True.
        """
        if pretrain:
            for col in self.variable_colnames:
                var_counts = training_df.groupby([self.group_colname, col]).size().reset_index(name = 'N')
                unique_var_counts = var_counts.groupby([self.group_colname]).size().reset_index(name = 'N')
                if 1 in unique_var_counts.N.unique():
                    self.variable_colnames = [item for item in self.variable_colnames if item != col]
                    if self.choice1_variable_colnames is not None:
                        self.choice1_variable_colnames = [item for item in self.choice1_variable_colnames if item != col]
                    if self.choice2_variable_colnames is not None:
                        self.choice2_variable_colnames = [item for item in self.choice2_variable_colnames if item != col]
                    print('removing ' + col + ' from model due to insufficient within-choice set variance')
        else:
            # remove missing NT variables
            self.variable_colnames = [var for var in self.variable_colnames if '_-' not in var]
            if self.choice1_variable_colnames is not None:
                self.choice1_variable_colnames = [var for var in self.choice1_variable_colnames if '_-' not in var]
            if self.choice2_variable_colnames is not None:
                self.choice2_variable_colnames = [var for var in self.choice2_variable_colnames if '_-' not in var]


    def remove_zero_set_counts(self, training_df):
        """
        Removes groups from the DataFrame that have zero counts in the specified count column.

        Parameters:
            training_df (pd.DataFrame): DataFrame containing the data.

        Returns:
            pd.DataFrame: Updated DataFrame with groups having zero counts removed.
        """
        count_sums = training_df.groupby([self.group_colname])[self.count_colname].sum()
        if 0 in count_sums.unique():
            zeros = count_sums[count_sums == 0]
            zero_groups = list(zeros.index)
            training_df = training_df[~training_df[self.group_colname].isin(zero_groups)]
            print('removing ' + str(zero_groups) + ' groups from model due to zero counts')
        return(training_df)

    def filter_input_domain_space(self, df):
        """
        Filters the input DataFrame based on a predefined domain space, ensuring data compatibility.

        Parameters:
            df (pd.DataFrame): DataFrame containing the data to be filtered.

        Returns:
            pd.DataFrame: Filtered DataFrame that fits within the specified domain space.
        """
        domain_file = self.params.R_input_domain_data_path()

        domain_data = pd.read_csv(domain_file, sep = '\t')

        cols = self.input_group_colname + self.input_choice_colname
        if self.input_choice_colname != None:
            cols = cols + [self.input_choice2_colname]
            cols = cols + ['frame_type', 'frame_stop']

        # filter input df
        filtered_df = pd.merge(domain_data, df, how='inner', on=cols)

        return filtered_df


class DataTransformer(DataPreprocessor):
    """
    A class for transforming and preprocessing data for modeling.
    This class introduces methods for initializing model coefficients randomly,
    preprocessing data (including categorical variable transformation and expansion of multivariable
    columns), and preparing data matrices for model input.

    Attributes inherited from DataPreprocessor are complemented with:
    - original_variable_colnames (list): Original variable column names before transformation.
    - original_group_colname (str): Original group identifier column name.
    - original_choice_colname (str): Original choice identifier column name.
    - coefs (ndarray): Model coefficients, initialized as None and can be set using get_random_coefs or externally.

    Args:
        training_df (pd.DataFrame): The DataFrame containing training data.
        variable_colnames (list of str): List of variable column names to be used in the model.
        count_colname (str): Column name for count data, used in count-based models.
        group_colname (str): Column name for group identifiers.
        repeat_obs_colname (str): Column name for identifiers of repeated observations.
        choice_colname (str): Column name for choice identifiers.
        params (dict): Additional parameters for data preprocessing and transformation.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        super().__init__(training_df, variable_colnames, None, None, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname
        self.coefs = None

    def get_random_coefs(self):
        """
        Generates random coefficients for model initialization.

        Returns:
            ndarray: Randomly generated coefficients as a NumPy array.
        """
        # Set random coefficients for model initialization
        key = jax.random.PRNGKey(123)
        coefs = jax.random.normal(key, shape=(len(self.variable_colnames), 1))
        return coefs

    def preprocess_data(self, df, pretrain=True):
        """
        Applies a series of preprocessing steps to the provided DataFrame, preparing it for model training.
        The preprocessing steps include filtering domain space, expanding multivariable columns,
        removing entries with zero counts, and transforming categorical variables into numerical
        or encoded forms as required.

        Args:
            df (pd.DataFrame): The DataFrame to be preprocessed.
            pretrain (bool, optional): Flag indicating whether preprocessing is for pretraining purposes. Defaults to True.

        Returns:
            pd.DataFrame: The preprocessed DataFrame, ready for modeling.
        """
        # filter for possible sites
        if 'ligation_mh' in self.input_choice_colname:
            df = self.filter_input_domain_space(df)


        # Transform group column lists into strings
        df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # remove zero counts
        if pretrain:
            df = self.remove_zero_set_counts(df)

        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")

        # Transform choice column lists into strings
        df = self.expand_multivariable(df, "original_choice_colname", "choice_colname")
        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")
        df = self.transform_categorical_response_vars(df, "original_choice_colname", "choice_colname")

        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            df = self.transform_categorical_vars(df, col, pretrain)

        # check for within-choice-set variance
        self.check_within_set_variance(df, pretrain)

        return(df)

    def reset_weighted_observations(self, counts_mat):
        """
        Normalizes the counts matrix so that the sum of all counts equals 1, effectively resetting the observations
        to a weighted distribution. This function is useful for adjusting the counts data for modeling purposes,
        ensuring that the data represents probabilities or normalized weights rather than absolute counts.

        Args:
            counts_mat (ndarray): A NumPy array containing counts data. This matrix should have dimensions
                                  corresponding to the number of groups by the number of choices.

        Returns:
            ndarray: A NumPy array of the same shape as `counts_mat`, with values normalized so that
                     the total sum of the matrix equals 1.
        """
        mat_sum = jnp.nansum(counts_mat)
        reset_mat = counts_mat/mat_sum
        return(reset_mat)

    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares and returns data matrices necessary for modeling from the preprocessed DataFrame.

        The method prepares matrices for variables, counts, non-repeat groups, and a mask indicating valid data points,
        which are essential for training count-based and choice-based models.

        Args:
            df (pd.DataFrame): The DataFrame from which matrices are to be generated, typically after preprocessing.
            pretrain (bool, optional): Indicates if the matrix preparation is for pretraining purposes. Defaults to True.
            replace_object (str, optional): Name of the attribute to replace with the preprocessed DataFrame. Defaults to None.
            return_df (bool, optional): Flag indicating whether to return the preprocessed DataFrame along with the matrices. Defaults to False.

        Returns:
            tuple: Contains the variable matrix, counts matrix, non-repeat groups matrix, and mask matrix as NumPy arrays.
                   Optionally returns the preprocessed DataFrame if `return_df` is True.
        """
        df = self.preprocess_data(df, pretrain)

        # Fill in missing counts with zero
        df[self.count_colname] = df[self.count_colname].fillna(0).astype(float)

        # Get matrix shapes
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])

        final_shape = (len(groups), len(choices), len(self.variable_colnames))
        int_shape = (len(groups), len(self.variable_colnames), len(choices))
        counts_shape = (len(groups), len(choices), 1)

        # map variables to indices
        var_mapping = {key: value for key, value in enumerate(self.variable_colnames)}
        var_mapping_df = pd.DataFrame(list(var_mapping.items()), columns=['VarIndex', 'VarName'])

        # get counts matrix
        pivot_counts = df.pivot_table(index=[self.group_colname],
                                      columns=[self.choice_colname],
                                      values=self.count_colname,
                                      fill_value=0)
        counts_mat = jnp.array(pivot_counts).reshape(counts_shape)

        # get variable matrix
        pivot_vars = df.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        pivot_vars = pivot_vars.reindex(var_mapping_df.VarName.tolist(), axis=1, level = 0)

        mat = jnp.array(pivot_vars).reshape(int_shape)
        mat = mat.transpose((0, 2, 1))
        assert mat.shape == final_shape, "Variable matrix is the incorrect dimension"

        # Assuming groups is a list of unique group identifiers
        group_indices = {group: idx for idx, group in enumerate(groups)}
        nonrepeat_groups_mat = jnp.array([group_indices[group] for group in df[self.group_colname].unique()])

        # get mask matrix (1 for valid entries or 0 for nonvalid)
        df['indicator'] = 1.0
        pivot_mask = df.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=['indicator'],
                                    fill_value=0)
        mask_mat = jnp.array(pivot_mask).reshape(counts_shape)

        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return mat, counts_mat, nonrepeat_groups_mat, mask_mat, df
        else:
            return mat, counts_mat, nonrepeat_groups_mat, mask_mat

    def get_coefficients(self, coefs=None):
        """
        Retrieves the model coefficients, optionally using an external set of coefficients. This method formats
        the coefficients into a structured dictionary, including both the coefficient values and their standard errors,
        if available. It also adjusts for dropped contrast variables by computing the implied value of the base category
        from the contrasted variables.

        Args:
            coefs (ndarray, optional): An array of coefficients to be used instead of the instance's attribute.
                                       If not provided, the method will use the `self.coefs` attribute. Defaults to None.

        Returns:
            list of dict: A list where each element is a dictionary containing the 'coefficient' name, its 'value',
                          and its 'error'. The list includes one dictionary for each model variable, including
                          inferred coefficients for base categories of categorical variables transformed into contrasts.

        Notes:
            - The method assumes that `self.coefs` and `self.standard_errors` are properly initialized and match
              the dimensions of the model variables.
            - For categorical variables transformed into contrasts, the base category's coefficient is computed
              as the negative sum of the contrast coefficients, and its error is set to None since it's derived
              indirectly.
        """
        if coefs is None:
            coefs = self.coefs
        coefs = coefs.flatten().tolist()
        errors = self.standard_errors.flatten().tolist()
        if self.choice1_variable_colnames is not None:
            choice_sum = self.choice1_variable_colnames
        else:
            choice_sum = []
        if self.choice2_variable_colnames is not None:
            choice_sum = choice_sum + self.choice2_variable_colnames

        if choice_sum != []:
            self.variable_colnames = self.choice1_variable_colnames + self.choice2_variable_colnames
        d = {n: {'value':c, 'error':e} for n, c, e in zip(self.variable_colnames, coefs, errors)}
        for col in list(set(self.original_variable_colnames) - set(self.variable_colnames)):
            contrast_vars, missing_var = self.get_dropped_contrast_var(self.training_df, col)
            for element in contrast_vars:
                if element not in self.variable_colnames:
                    contrast_vars.remove(element)
            d[missing_var] = {'value':-1*sum(d[var]['value'] for var in contrast_vars), 'error':None}

        data = [{'coefficient': name, 'value': values['value'], 'error': values['error']} for name, values in d.items()]
        return data

    def get_coefficients_df(self):
        """
        Constructs and returns a DataFrame of model coefficients, including additional information like base, position, and side
        derived from the coefficient names.

        This method simplifies the interpretation of coefficients by parsing their names and organizing them into a structured format.

        Returns:
            pd.DataFrame: A DataFrame containing detailed information on each model coefficient, including its value and error if available.
        """
        df = pd.DataFrame(self.get_coefficients())
        df.reset_index(drop=True, inplace=True)

        # get bases
        df['base'] = None
        df.loc[df.coefficient.str.contains('motif'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count') & df.coefficient.str.contains('prop'), 'base'] = df.coefficient.str.split('_').str[-2]


        # get positions
        df['position'] = None
        df.loc[df.coefficient.str.contains('motif'), 'position'] = df.coefficient.str.split('_').str[-2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]
        df.loc[df.coefficient.str.contains('mh_count'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]


        # get side
        df['side'] = None
        df.loc[df.coefficient.str.contains('motif'), 'side'] = df.coefficient.str.split('_').str[3]
        df.loc[df.coefficient.str.contains('base_count'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_count'), 'side'] = df.coefficient.str.split('_').str[2]

        # fix trimming specific
        df['trim_type'] = None
        df.loc[df.coefficient.str.contains('_trim'), 'trim_type'] = df.coefficient.str.split('_trim').str[0] + '_trim'
        df.loc[df.coefficient.str.contains('_length'), 'trim_type'] = df.coefficient.str.split('_length').str[0].str.split('_').str[-1] + '_trim'

        # simplify coefficients
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_prop'), 'coefficient'] = 'mh_prop_length_interaction'
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_count'), 'coefficient'] = 'mh_count_length_interaction'
        df.loc[df.coefficient.str.contains('motif'), 'coefficient'] = 'motif'
        df.loc[df.coefficient.str.contains('base_count'), 'coefficient'] = 'base_count'
        df.loc[df.coefficient.str.contains('mh_prop') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_prop'
        df.loc[df.coefficient.str.contains('mh_count') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_count'

        return(df)


class ConditionalLogisticRegressor(DataTransformer):
    """
    Implements a conditional logistic regression model, extending the DataTransformer class
    with functionalities specific to logistic regression, including fitting the model using
    gradient descent, cross-validation for hyperparameter tuning, and model evaluation.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list of str): List of variable column names.
        count_colname (str): Name of the column representing count data.
        group_colname (str): Name of the column representing group identifiers.
        repeat_obs_colname (str): Name of the column representing repeated observation identifiers.
        choice_colname (str): Name of the column representing choice identifiers.
        params (dict): Dictionary of parameters for data preprocessing and transformation.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning. Default is 10.

    Attributes:
        Inherits all attributes from DataTransformer.
        l2reg (float): L2 regularization strength. Default is 0.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning.
        l2reg_grid (pd.DataFrame): DataFrame containing the grid search results for L2 regularization tuning.
        maxiter (int): Maximum number of iterations for the optimization algorithm. Default is 1000.
        tolerance (float): Tolerance for convergence in the optimization algorithm. Default is 1e-6.
        step (float): Step size for the gradient descent optimization. Default is 0.1.
        coefs (ndarray): Fitted coefficients of the logistic regression model.
        training_info: Stores information returned by the optimizer during model fitting.
        cov_matrix (ndarray): Covariance matrix of the model coefficients, computed from the Hessian matrix.
        standard_errors (ndarray): Standard errors of the model coefficients, derived from the covariance matrix.

    Methods:
        fit: Fits the conditional logistic regression model using gradient descent.
        get_prob: Computes choice probabilities for given variables and coefficients.
        cross_entropy: Calculates the cross-entropy loss.
        l2regularization: Computes the L2 regularization term.
        loss_fn: Computes the loss function for optimization, including cross-entropy and L2 regularization.
        grid_search_cv: Performs grid search cross-validation for hyperparameter tuning of L2 regularization.
        get_l2reg: Determines the optimal L2 regularization strength based on cross-validation.
        train_model: Trains the conditional logistic regression model, optionally with L2 regularization.
        get_hessian: Computes the Hessian matrix of the loss function with respect to the coefficients.
        get_cov_matrix: Computes the covariance matrix of the model coefficients from the Hessian matrix.
        get_errors: Computes the standard errors of the model coefficients from the covariance matrix.
        save_model: Saves the model to a file for later use.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params, l2kfold=10):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix, self.mask_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-6
        self.step = 0.1
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    # Get probability for input parameters given coefficients
    def get_prob(self, variables, mask, coefs=None):
        """
        Compute choice probabilities for given variables and coefficients.

        Args:
            variables (ndarray): Data matrix of variables.
            mask (ndarray): Data matrix of variable mask.
            coefs (ndarray, optional): Coefficients for the logistic regression model. If not provided, the trained coefficients will be used.

        Returns:
            ndarray: Choice probabilities for each choice given the variables.
        """
        if coefs is None:
            coefs = self.coefs
        assert coefs is not None, "Need to train the model before making predictions!"
        # Compute the logits for each choice
        cov = jnp.dot(variables, coefs)
        reshape = jnp.squeeze(cov)
        reshaped_mask = jnp.squeeze(mask)

        # Calculate the probability of the observed choices
        # Dimensions of this matrix are groups x choices
        # replace missing choices with -INF so that they will not count towards probability
        probs = jax.nn.softmax(jnp.where(reshaped_mask, reshape, -jnp.inf))

        return probs

    # Get cross-entropy loss
    def cross_entropy(self, probs, counts):
        """
        Calculate the cross-entropy loss.

        Args:
            probs (ndarray): Choice probabilities.
            counts (ndarray): Counts of choices.

        Returns:
            float: Cross-entropy loss.
        """
        counts_reshape = jnp.squeeze(counts)
        loss = -jnp.sum(jnp.log(jnp.where(probs==0, 1, probs)) * counts_reshape)
        return loss

    def l2regularization(self, coefs, size, l2reg):
        """
        Compute L2 regularization term for only non-base-count and non-motif coefficients.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            size (int): Size of the training data
            l2reg (float): L2 regularization strength.

        Returns:
            float: L2 regularization term.
        """
        var_list = [var for var in self.variable_colnames if 'base_count' not in var or 'interaction' in var]
        var_list = [var for var in var_list if 'motif' not in var]

        mh_list = [var in var_list for var in self.variable_colnames]

        if any(mh_list):
            binary_mh_list = [jnp.array(b, dtype=int) for b in mh_list]
            binary_mh_jnp = jnp.array(binary_mh_list).reshape(-1, 1)
            coef_subset = coefs * binary_mh_jnp
            c = jnp.nansum(coef_subset**2)
        else:
            c = 0
        return(0.5*(1/size)*l2reg*c)

    # Compute the loss function
    def loss_fn(self, coefs, variables, counts, mask, l2reg=0):
        """
        Compute the loss function for optimization.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            variables (ndarray): Data matrix of variables.
            counts (ndarray): Counts of choices.
            mask (ndarray): Data matrix of mask.
            l2reg (float): L2 regularization strength.

        Returns:
            float: Total loss, including cross-entropy and L2 regularization.
        """
        probs = self.get_prob(variables, mask, coefs)
        if probs.ndim == 1:
            probs = probs.reshape((probs.shape[0], 1))
        size = jnp.count_nonzero(mask).item()

        loss = self.cross_entropy(probs, counts) + self.l2regularization(coefs, size, l2reg)
        return loss

    def fit(self, variable_matrix, counts_matrix, mask_matrix, l2reg, maxiter, tol, step, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            mask_matrix (ndarray): Data matrix of masks.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            tol (float): tolerance for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.BFGS(fun=self.loss_fn, maxiter=maxiter, tol=tol, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         variables=variable_matrix,
                         counts=counts_matrix,
                         mask=mask_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg):
        """
        Computes the cross-validation loss for a given L2 regularization strength across a specified number of folds. 
        This method splits the data into training and validation sets according to the fold count, trains the model on 
        each training set while applying the specified L2 regularization strength, and then calculates and returns the 
        loss on each validation set.

        Args:
            fold_count (int): The number of folds to use for cross-validation.
            l2reg (float): The L2 regularization strength to use when training the model.

        Returns:
            list of float: A list containing the cross-validation loss for each fold. The length of the list will 
                           be equal to `fold_count`. Each element in the list represents the loss computed on the 
                           validation set of the corresponding fold.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_data, val_data = self.variable_matrix[train_index], self.variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_mask, val_mask = self.mask_matrix[train_index], self.mask_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            model = self.fit(train_data, train_counts, train_mask, l2reg, self.maxiter, self.tolerance, self.step, self.initial_coefs)

            # Compute the loss on the validation data
            loss = self.loss_fn(model.params,
                                val_data,
                                val_counts,
                                val_mask)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-5, 5)] + [0]):
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
            scores = self.cv_loss(l2kfold, l2reg)

            # Calculate the mean cross-validation score
            mean_score = np.mean(scores)

            temp = pd.DataFrame({'l2reg':[l2reg], 'mean_CV_loss':[mean_score]})
            results = pd.concat([results, temp], ignore_index=True)

        return(results)

    def get_l2reg(self):
        """
        Determine the optimal L2 regularization strength.

        Returns:
            float: Optimal L2 regularization strength.
        """
        self.l2kfold = min(len(np.unique(self.nonrepeat_grp_matrix)), 10)
        grid = self.grid_search_cv(self.l2kfold)
        self.l2reg_grid = grid
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, l2reg_value=None, maxiter=None, tolerance=None, step=None):
        """
        Trains the conditional logistic regression model with optional L2 regularization. This method allows for the model
        to be trained with or without L2 regularization, and various training parameters can be specified to control the
        optimization process.

        Args:
            l2 (bool): Indicates whether L2 regularization should be applied during model training. If `True`, L2
                       regularization is applied using either the specified `l2reg_value` or an optimal value determined
                       via cross-validation. Default is `False`, indicating no regularization.
            l2reg_value (float, optional): Specifies the L2 regularization strength to be used if L2 regularization is
                                           enabled. If `None` and `l2` is `True`, the method performs cross-validation to
                                           determine the optimal regularization strength. Default is `None`.
            maxiter (int, optional): The maximum number of iterations for the optimization algorithm. If not provided,
                                     the model uses a default value.
            tolerance (float, optional): The convergence tolerance for the optimization algorithm. Optimization stops when
                                         the loss value change is less than this tolerance between iterations. If not
                                         provided, the model uses a default value.
            step (float, optional): The step size for the gradient descent optimization algorithm. Only applicable if
                                    a gradient descent-based optimizer is used. If not provided, the model uses a default value.

        Returns:
            self: Returns an instance of `ConditionalLogisticRegressor` after training. The trained coefficients can
                  be accessed through the `coefs` attribute of the instance.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        if l2 is False:
            self.l2reg = 0
        else:
            if l2reg_value is None:
                self.l2reg = self.get_l2reg()
            else:
                self.l2reg = l2reg_value

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

        res = self.fit(self.variable_matrix,
                       self.counts_matrix,
                       self.mask_matrix,
                       self.l2reg,
                       maxiter,
                       tolerance,
                       step,
                       self.initial_coefs)

        self.coefs = res.params
        self.training_info = res
        self.maxiter = maxiter
        self.tolerance = tolerance
        self.step = step
        self.cov_matrix = self.get_cov_matrix(self.coefs, self.variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        self.standard_errors = self.get_errors(self.coefs, self.variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        return self

    def get_hessian(self, coefs, variables, counts, mask, l2reg=0):
        """
        Computes the Hessian matrix of the loss function with respect to the model coefficients. This matrix is
        crucial for understanding the curvature of the loss function around the optimal coefficients and is used
        to compute the covariance matrix of the coefficients.

        Args:
            coefs (ndarray): A NumPy array containing the current coefficients of the model.
            variables (ndarray): A NumPy array containing the predictor variables.
            counts (ndarray): A NumPy array containing the counts of each choice.
            mask (ndarray): A NumPy array indicating valid choices with a boolean value.
            l2reg (float, optional): The strength of L2 regularization. Default is 0.

        Returns:
            ndarray: The Hessian matrix of the loss function evaluated at the given coefficients.
        """
        # Wrapper function
        def wrapper_loss_fn(coefs):
            return self.loss_fn(coefs, variables, counts, mask, l2reg)

        hessian_fn = jax.hessian(wrapper_loss_fn, argnums=0)
        hessian_matrix = hessian_fn(coefs.reshape(-1))
        return hessian_matrix

    def get_cov_matrix(self, coefs, variables, counts, mask, l2reg=0):
        """
        Computes the covariance matrix of the model coefficients, derived from the inverse of the Hessian matrix.

        Args:
            coefs (ndarray): A NumPy array containing the current coefficients of the model.
            variables (ndarray): A NumPy array containing the predictor variables.
            counts (ndarray): A NumPy array containing the counts of each choice.
            mask (ndarray): A NumPy array indicating valid choices with a boolean value.
            l2reg (float, optional): The strength of L2 regularization. Default is 0.

        Returns:
            ndarray: The covariance matrix of the model coefficients.
        """
        hess_mat = self.get_hessian(coefs, variables, counts, mask, l2reg)
        cov_matrix = jnp.linalg.inv(hess_mat)
        return cov_matrix

    def get_errors(self, coefs, variables, counts, mask, l2reg=0):
        """
        Computes the standard errors of the model coefficients, which are derived from the square root of the diagonal
        elements of the covariance matrix. These standard errors are essential for assessing the statistical significance
        of the coefficients.

        Args:
            coefs (ndarray): A NumPy array containing the current coefficients of the model.
            variables (ndarray): A NumPy array containing the predictor variables.
            counts (ndarray): A NumPy array containing the counts of each choice.
            mask (ndarray): A NumPy array indicating valid choices with a boolean value.
            l2reg (float, optional): The strength of L2 regularization. Default is 0.

        Returns:
            ndarray: An array containing the standard errors of the model coefficients.
        """
        cov = self.get_cov_matrix(coefs, variables, counts, mask, l2reg)
        standard_errors = np.sqrt(np.diag(cov))
        return standard_errors

    def save_model(self, file_path):
        """
        Saves the current model to a file. This includes the model coefficients and configuration but excludes
        the training data and matrices to minimize file size.

        Args:
            file_path (str): The path to the file where the model should be saved.

        Note:
            - Before saving, the method removes large attributes such as the training data and matrices from the model
              to prevent excessive file sizes. These attributes must be reloaded or recomputed for further model use.
              be overwritten.
        """
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)


class ConditionalLogisticRegressionPredictor(DataTransformer):
    """
    Facilitates making predictions and computing loss using a pre-trained Conditional Logistic Regression model
    on new or unseen data. This class handles the preprocessing of new data to match the format expected by the
    model and computes probabilities and loss values using the model's coefficients.

    Inherits from:
        DataTransformer: Inherits data transformation capabilities for consistent preprocessing.

    Args:
        model (ConditionalLogisticRegressor): A trained ConditionalLogisticRegressor model instance.
        variable_colnames (list of str): List of variable column names as used in the training phase.
        count_colname (str): The column name for the count variable, used in the training phase.
        group_colname (list of str): List of column names representing group identifiers, as used in the training phase.
        choice_colname (list of str): List of column names representing choice identifiers, as used in the training phase.
        params (dict): Additional parameters for data preprocessing and transformation.

    Attributes:
        model (ConditionalLogisticRegressor): Stores the trained logistic regression model.
        original_variable_colnames, original_group_colname, original_choice_colname:
            Store the original column names for variables, groups, and choices, facilitating data preparation for prediction.
        variable_matrix, counts_matrix, nonrepeat_grp_matrix:
            Numpy arrays prepared for prediction, storing variables, counts, and non-repeating group identifiers respectively.

    Methods:
        predict(new_df): Uses the trained model to make predictions on new data.
        compute_loss(new_df): Computes the loss on new data using the model's loss function.
    """
    def __init__(self, model, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, training_params, validation_params):
        super().__init__(None, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, validation_params)
        self.model = model
        if not isinstance(model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname

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
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_prob(variable_matrix, mask_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice_cols = self.get_mapping_dict(new_df, self.choice_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice_cols)))
        prob_df = pd.DataFrame(probs)
        prob_df.columns = list(choice_cols.keys())
        prob_df[self.group_colname] = list(group_cols.keys())
        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname],
                            var_name=self.choice_colname,
                            value_name='predicted_prob')
        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname],
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
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(float(loss))


class ConditionalLogisticRegressionEvaluator(DataTransformer):
    """
    A class dedicated to evaluating a trained Conditional Logistic Regression model. It provides functionalities
    to load a trained model, calculate log loss on training data, expected log loss through cross-validation,
    and log loss on validation data. It compiles evaluation results into a DataFrame for easy analysis and comparison.

    Inherits from:
        DataTransformer: For consistent data preprocessing and transformation.

    Args:
        model_path (str): Path to the file containing the saved trained model.
        params (dict): Dictionary of parameters used for data preprocessing and model evaluation.
        training_df (pd.DataFrame, optional): DataFrame containing the training data. Default is None.
        validation_df (pd.DataFrame, optional): DataFrame containing the validation data. Default is None.

    Attributes:
        model (ConditionalLogisticRegressor): The loaded trained logistic regression model.
        log_loss (float): Log loss calculated on the training data.
        expected_log_loss (float): Expected log loss calculated through cross-validation on the training data.
        validation_df (pd.DataFrame): DataFrame containing the validation data.
    """
    def __init__(self, model_path, training_params, validation_params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.training_params = training_params
        self.validation_params = validation_params
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, validation_params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        """
        Loads the trained Conditional Logistic Regression model from a specified file path.

        Args:
            file_path (str): The path to the file containing the trained model.

        Returns:
            ConditionalLogisticRegressor: The loaded model.
        """
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        """
        Calculates the log loss on the training dataset using the loaded model.

        Returns:
            float: The log loss value on the training dataset.
        """
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def calculate_expected_log_loss(self, fold_count=20):
        """
        Calculates the expected log loss on the training dataset using cross-validation.

        Args:
            fold_count (int, optional): The number of folds to use for cross-validation. Default is 20.

        Returns:
            float: The expected log loss value averaged across all folds.
        """
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the expected loss on the training data
        self.model.variable_matrix = variable_matrix
        self.model.counts_matrix = counts_matrix
        self.model.nonrepeat_grp_matrix = nonrepeat_grp_matrix

        expected = self.model.cv_loss(fold_count, self.model.l2reg)
        e_loss = sum((1/fold_count) * np.array(expected))
        return(e_loss)

    def calculate_validation_log_loss(self):
        """
        Calculates the log loss on the validation dataset using the loaded model.

        Returns:
            float: The log loss value on the validation dataset.
        """
        assert self.validation_df is not None, 'No input validation dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def compile_evaluation_results_df(self, calculate_validation_loss = False, calculate_expected_loss=False):
        """
        Compiles the model evaluation results, including log loss on training data, expected log loss through
        cross-validation, and log loss on validation data, into a DataFrame.

        Args:
            calculate_validation_loss (bool, optional): Whether to calculate and include log loss on validation data.
                                                         Default is False.
            calculate_expected_loss (bool, optional): Whether to calculate and include expected log loss through
                                                      cross-validation. Default is False.

        Returns:
            pd.DataFrame: A DataFrame containing the compiled evaluation results and model parameters.
        """
        result = {'training_annotation_type': [self.training_params.annotation_type],
                  'productivity':[self.training_params.productivity],
                  'motif_length_5_end':[self.training_params.left_nuc_motif_count],
                  'motif_length_3_end':[self.training_params.right_nuc_motif_count],
                  'motif_type':[self.training_params.motif_type],
                  'gene_weight_type':[self.training_params.gene_weight_type],
                  'upper_bound':[self.training_params.upper_trim_bound],
                  'lower_bound':[self.training_params.lower_trim_bound],
                  'insertion_bound':[self.training_params.insertions],
                  'model_type':[self.training_params.model_type],
                  'base_count_5end_length':[10],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        final = pd.DataFrame()

        if calculate_expected_loss is True:
            self.expected_log_loss = self.calculate_expected_log_loss()

            # add expected loss result
            e = results_df.copy()
            e['loss_type'] = 'Expected log loss across training data'
            e['log_loss'] = self.expected_log_loss
            e['validation_annotation_type'] = [self.validation_params.annotation_type]
            e['validation_productivity'] = [self.validation_params.productivity]

            final = pd.concat([final, e], axis = 0)

        if calculate_validation_loss is True:
            val_loss = self.calculate_validation_log_loss()
            val = results_df.copy()
            val['loss type'] = 'Log loss on validation data'
            val['log_loss'] = val_loss
            val['validation_annotation_type'] = [self.validation_params.annotation_type]
            val['validation_productivity'] = [self.validation_params.productivity]

            final = pd.concat([final, val], axis = 0)

        return(final)
