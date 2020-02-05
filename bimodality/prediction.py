"""
...

"""

###############################################################################
# Notes

# If the count of bootstrap shuffles is too small, then it is likely not to
# have any values equal to or greater than the real value, in which case the
# probability is zero.

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import copy
import random
import itertools
import functools
import multiprocessing
import datetime
import time


# Relevant

import numpy
import pandas
import scipy.stats
#from sklearn.linear_model import LinearRegression
import statsmodels

# Custom

import distribution
import modality
import plot
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_selection = os.path.join(dock, "selection", "tight")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation_gencode.pickle"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )
    path_genes_selection = os.path.join(
        path_selection, "genes.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_persons_properties": data_persons_properties,
        "genes_selection": genes_selection,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


##########
# Regression


def define_regression_variables():
    """
    Defines a list of variables' names for regression analysis.

    arguments:

    raises:

    returns:
        (dict<list<str>>): names of independent variables for regression

    """

    independence = [
        "female",
        "age",
        "body",
        "hardiness",
        "season_sequence",

        "genotype_1",
        "genotype_2",
        "genotype_3",
        "genotype_4",
        "genotype_5",
        "genotype_6",
        "genotype_7",
        "genotype_8",
        "genotype_9",
        "genotype_10",

        "delay",

        "tissues_1",
        "tissues_2",
        "tissues_3",
        #"tissues_4", # <- omit
        #"tissues_5", # <- omit
        #"tissues_6", # <- omit
        #"tissues_7", # <- omit
        #"tissues_8", # <- omit
        #"tissues_9", # <- omit
        #"tissues_10", # <- omit

        "facilities_1",
        "facilities_2",
        #"facilities_3", # <- omit

        "batches_isolation_1",
        "batches_isolation_2",
        "batches_isolation_3",
        #"batches_isolation_4", # <- omit
        #"batches_isolation_5", # <- omit
        #"batches_isolation_6", # <- omit
        #"batches_isolation_7", # <- omit
        #"batches_isolation_8", # <- omit
        #"batches_isolation_9", # <- omit
        #"batches_isolation_10", # <- omit

        "batches_analysis_1",
        "batches_analysis_2",
        "batches_analysis_3",
        #"batches_analysis_4", # <- omit
        #"batches_analysis_5", # <- omit
        #"batches_analysis_6", # <- omit
        #"batches_analysis_7", # <- omit
        #"batches_analysis_8", # <- omit
        #"batches_analysis_9", # <- omit
        #"batches_analysis_10", # <- omit
    ]
    hypothesis = [
        "female",
        "age",
        "body",
        "hardiness",
        "season_sequence",
    ]

    # Compile information.
    information = dict()
    information["independence"] = independence
    information["hypothesis"] = hypothesis

    # Return information.
    return information


def organize_dependent_independent_variables(
    variables=None,
    data_persons_properties=None,
    data_signals_genes_persons=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        variables (list<str>): names of independent variables for regression
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        data_signals_genes_persons (object): Pandas data frame of pan-tissue
            signals across genes and persons

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables

    """

    # Remove all columns from persons properties except the covariates
    # Copy data.
    data_persons_properties = data_persons_properties.copy(deep=True)
    data_signals_genes_persons = data_signals_genes_persons.copy(
        deep=True
    )
    # Organize data.
    data_persons_properties = data_persons_properties.loc[
        :, data_persons_properties.columns.isin(variables)
    ]

    # Join signals on persons' properties.
    data_variables = data_persons_properties.join(
        data_signals_genes_persons,
        how="left",
        on="person"
    )

    # Return information.
    return data_variables


def regress_signal_ordinary_residuals(
    dependence=None,
    independence=None,
    proportion=None,
    data=None,
):
    """
    Regresses a quantitative continuous dependent variable against multiple
    independent variables and returns relevant parameters and statistics.

    Data format must have observations across rows and dependent and
    independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first dimension vector of
    observations and for each observation a second dimension vector of scalar
    values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    Description of formats for SKLearn...

    Format of dependent variable is an array of observations, and each
    observation is an array of features.
    [
        [1.3],
        [1.5],
        [1.2],
        ...
    ]

    arguments:
        dependence (str): name of dependent variable
        independence (list<str>): names of independent variables
        proportion (float): minimal proportion of total observations that must
            have nonmissing values
        data (object): Pandas data frame of values and associations between
            dependent and independent variables

    raises:

    returns:
        (dict): information about regression

    """

    if False:
        # Regress by SciKit Learn.
        regression = LinearRegression(
            fit_intercept=True,
            normalize=False,
            copy_X=True,
            n_jobs=None
        ).fit(
            values_independence,
            values_dependence,
            sample_weight=None
        )
        score = regression.score(
            values_independence,
            values_dependence,
            sample_weight=None
        )

    # Organize data.
    data = data.copy(deep=True)
    # Determine minimal count of observations.
    count = data.shape[0]
    threshold = proportion * count
    # Remove observations with any missing values.
    columns = copy.deepcopy(independence)
    columns.insert(0, dependence)
    data = data.loc[
        :, data.columns.isin(columns)
    ]
    data.dropna(
        axis="index",
        how="any",
        #subset=[dependence],
        inplace=True,
    )
    data = data[[*columns]]

    # Note
    # It is very important to keep track of the order of independent variables
    # in order to match parameters and probabilities to the correct variables.

    # Determine whether data have sufficient observations for regression.
    if data.shape[0] > threshold:
        # Extract values of dependent and independent variables.
        values_dependence = data[dependence].to_numpy()
        values_independence = data.loc[ :, independence].to_numpy()
        # Introduce constant value for intercept.
        values_independence_intercept = statsmodels.api.add_constant(
            values_independence,
            prepend=True,
        )
        # Define model.
        model = statsmodels.api.OLS(
            values_dependence,
            values_independence_intercept,
            missing="drop",
        )
        report = model.fit()
        #utility.print_terminal_partition(level=2)
        #print(data)
        #print(independence)
        #print(values_independence)
        #print(report.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(report))
        #print(report.params)
        #print(report.pvalues)

        # Compile information.
        counter = 1
        parameters = dict()
        probabilities = dict()
        parameters["intercept_parameter"] = report.params[0]
        probabilities["intercept_probability"] = report.pvalues[0]
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameters[parameter] = report.params[counter]
            probability = str(variable + ("_probability"))
            probabilities[probability] = report.pvalues[counter]
            counter += 1
            pass
        information = {
            "freedom": report.df_model,
            "observations": report.nobs,
            "r_square": report.rsquared,
            "r_square_adjust": report.rsquared_adj,
            "log_likelihood": report.llf,
            "akaike": report.aic,
            "bayes": report.bic,
        }
        information.update(parameters)
        information.update(probabilities)
    else:
        # Compile information.
        #probabilities = list(
        #    itertools.repeat(float("nan"), len(values_independence_intercept))
        #)
        parameters = dict()
        probabilities = dict()
        parameters["intercept_parameter"] = float("nan")
        probabilities["intercept_probability"] = float("nan")
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameters[parameter] = float("nan")
            probability = str(variable + ("_probability"))
            probabilities[probability] = float("nan")
            pass
        information = {
            "freedom": float("nan"),
            "observations": float("nan"),
            "r_square": float("nan"),
            "r_square_adjust": float("nan"),
            "log_likelihood": float("nan"),
            "akaike": float("nan"),
            "bayes": float("nan"),
        }
        information.update(parameters)
        information.update(probabilities)
    # Return information.
    return information


def regress_cases(
    cases=None,
    variables=None,
    data_variables=None,
):
    """
    Drives iterative regression across cases.

    arguments:
        cases (list<str>): names of independent cases for which to regress
            variables
        variables (list<str>): names of independent variables for regression
        data_variables (object): Pandas data frame of dependent and independent
            variables

    raises:

    returns:
        (object): Pandas data frame of parameters and statistics from
            regressions across cases

    """

    # Regress.
    records = list()
    for case in cases:
        #utility.print_terminal_partition(level=1)
        #print(case)
        report = regress_signal_ordinary_residuals(
            dependence=case,
            independence=variables,
            proportion=0.5,
            data=data_variables,
        )
        report["case"] = case
        records.append(report)
        pass
    # Organize data.
    data_regression = utility.convert_records_to_dataframe(
        records=records
    )
    data_regression.set_index(
        ["case"],
        append=False,
        drop=True,
        inplace=True
    )
    columns = list()
    columns.append("observations")
    columns.append("freedom")
    columns.append("r_square")
    columns.append("r_square_adjust")
    columns.append("log_likelihood")
    columns.append("akaike")
    columns.append("bayes")
    columns.append("intercept_parameter")
    columns.append("intercept_probability")
    for variable in variables:
        columns.append(str(variable + ("_parameter")))
        columns.append(str(variable + ("_probability")))
        pass
    data_regression = data_regression[[*columns]]
    # Return information.
    return data_regression


def calculate_regression_discoveries(
    variables=None,
    threshold=None,
    data=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression variables for
            which to calculate false discovery rates
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (object): Pandas data frame of parameters and statistics from
            regressions across cases

    """

    def calculate_discoveries(
        values=None,
        threshold=None,
        ):
        report = statsmodels.stats.multitest.multipletests(
            values,
            alpha=threshold,
            method="fdr_bh",
            is_sorted=False,
        )
        discoveries = report[1]
        return discoveries


    # Copy data.
    data_discovery = data.copy(deep=True)
    # Iterate on relevant variables.
    for variable in variables:
        probability = str(variable + ("_probability"))
        discovery = str(variable + ("_discovery"))
        # Calculate false discovery rate from probability.
        data_discovery[discovery] = calculate_discoveries(
                values=data_discovery[probability].to_numpy(),
                threshold=threshold,
        )
        pass
    # Return information.
    return data_discovery


def scale_data_columns(
    columns=None,
    data=None,
):
    """
    Calculates the standard (z-score) of values within specific columns.

    arguments:
        columns (list<str>): names of columns for which to scale values
        data (object): Pandas data frame

    raises:

    returns:
        (object): Pandas data frame

    """

    # Calculate standard score.
    # This method inserts missing values if the standard deviation is zero.
    data = data.copy(deep=True)
    data_scale = data.apply(
        lambda column: (
            ((column - column.mean()) / column.std())
        ) if column.name in columns else column,
        axis="index",
    )
    return data_scale


def select_scale_regression_parameters(
    variables=None,
    data_regression_genes=None,
):
    """
    Selects and scales regression parameters.

    arguments:
        variables (list<str>): names of independent regression variables to
            select and scale
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (object): Pandas data frame of parameters and statistics from
            regressions across cases

    """

    # Select relevant variables.
    columns = list()
    for variable in variables:
        columns.append(str(variable + ("_parameter")))
        columns.append(str(variable + ("_probability")))
        columns.append(str(variable + ("_discovery")))
        pass
    data_selection = data_regression_genes.copy(deep=True)
    data_selection = data_selection.loc[
        :, data_selection.columns.isin(columns)
    ]
    # Scale regression parameters across genes.
    parameters = list()
    for variable in variables:
        parameters.append(str(variable + ("_parameter")))
        pass
    data_scale = scale_data_columns(
        columns=parameters,
        data=data_selection,
    )
    # Summarize standardization.
    data_mean = data_scale.aggregate(
        lambda x: x.mean(),
        axis="index"
    )
    utility.print_terminal_partition(level=2)
    print("Mean")
    print(data_mean)
    data_deviation = data_scale.aggregate(
        lambda x: x.std(),
        axis="index"
    )
    utility.print_terminal_partition(level=2)
    print("Standard deviation")
    print(data_deviation)

    # Return information.
    return data_scale


def separate_regression_gene_sets(
    genes_selection=None,
    genes_unimodal=None,
    genes_multimodal=None,
    data_regression_genes=None,
):
    """
    Drives iterative regression across cases.

    arguments:
        genes_selection (list<str>): identifiers of genes
        genes_unimodal (list<str>): identifiers of genes
        genes_multimodal (list<str>): identifiers of genes
        data_regression_genes (object): Pandas data frame of parameters and
            statistics from regressions across cases

    raises:

    returns:
        (dict): regression information across sets of genes

    """

    # Select regression data for specific genes.
    data_selection = data_regression_genes.copy(deep=True)
    data_selection = data_selection.loc[
        data_selection.index.isin(genes_selection), :
    ]
    data_unimodal = data_regression_genes.copy(deep=True)
    data_unimodal = data_unimodal.loc[
        data_unimodal.index.isin(genes_unimodal), :
    ]
    data_multimodal = data_regression_genes.copy(deep=True)
    data_multimodal = data_multimodal.loc[
        data_multimodal.index.isin(genes_multimodal), :
    ]

    # Compile information.
    information = dict()
    information["data_regression_genes_selection_scale"] = data_selection
    information["data_regression_genes_unimodal_scale"] = data_unimodal
    information["data_regression_genes_multimodal_scale"] = data_multimodal
    # Return information.
    return information


def summarize_regression(
    instances=None,
    variables=None,
    probabilities=None,
    type=None,
    threshold_discovery=None,
    data=None,
):
    """
    Summarizes a regression model's parameters across multiple independent
    instances.

    arguments:
        instances (int): count of independent instances
        variables (list<str>): names of independent variables in regression
        probabilities (bool): whether to summarize probabilities for variables
            across instances
        type (str): type of report for probabilities of variables across
            instances, either aggregate false discovery rate or percentage of
            instances within discovery threshold
        threshold_discovery (float): minimal false discovery rate value for
            instances
        data (object): Pandas data frame of information about a regression
            model across multiple independent observations

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    # Organize data.
    data = data.copy(deep=True)

    # Report mean model descriptors.
    descriptors = [
        "freedom",
        "observations",
        "r_square",
        "r_square_adjust",
        "log_likelihood",
        "akaike",
        "bayes",
    ]
    data_descriptors = data.copy(deep=True)
    data_descriptors = data_descriptors.loc[
        :, data_descriptors.columns.isin(descriptors)
    ]
    data_descriptors = data_descriptors[[*descriptors]]
    data_descriptors_mean = data_descriptors.aggregate(statistics.mean)
    utility.print_terminal_partition(level=3)
    print("Model descriptors...")
    print(data_descriptors_mean)

    # Report counts of observations with probabilities below threshold for each
    # covariate.
    if probabilities:
        utility.print_terminal_partition(level=3)
        print("Independent instances with threshold probabilities...")

        for variable in variables:
            # Calculate false discovery rates from probabilities.
            probabilities = data[variable].to_numpy()
            report = statsmodels.stats.multitest.multipletests(
                probabilities,
                alpha=0.05,
                method="fdr_bh",
                is_sorted=False,
            )
            discoveries = report[1]

            # Prepare report.
            if type == "discovery":
                # Aggregation of false discovery rates across independent
                # instances probably is not reasonable.
                discovery = scipy.stats.combine_pvalues(
                    discoveries,
                    method="stouffer",
                    weights=None,
                )[1]
                print(variable + ": " + str(round(discovery, 1)))
            elif type == "percentage":
                count = (discoveries <= threshold_discovery).sum()
                percentage = round((count / instances) * 100, 1)
                print(variable + ": " + str(percentage) + "%")
            pass
    pass



##########
# Tissue-tissue correlation


def calculate_pairwise_correlations(
    data=None,
    key=None,
):
    """
    Calculates the pairwise correlation coefficients for multiple pairs of
    features.

    Data has observations across rows and features across columns.

    index key           feature_1 feature_2 feature_3 feature_4 feature_5
    0     observation_1 ...       ...       ...       ...       ...
    1     observation_2 ...       ...       ...       ...       ...
    2     observation_3 ...       ...       ...       ...       ...
    3     observation_4 ...       ...       ...       ...       ...
    4     observation_5 ...       ...       ...       ...       ...

    arguments:
        key (str): name of column to match observations
        data (object): Pandas data frame of features across observations

    raises:

    returns:
        (dict<dict<dict<float>>>): matrix of counts of common observations and
            correlation coefficients from pairwise comparison of features

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_copy.set_index(
        keys=[key],
        drop=True,
        append=False,
        inplace=True,
    )
    features = data_copy.columns.values.tolist()
    # Determine pairs of features available for each observation.
    # Use sequence-specific permutations to simplify assembly of matrix.
    pairs = list(itertools.permutations(features, 2))

    # Collect counts and correlations across pairs of features.
    correlations = dict()
    # Iterate on pairs of features.
    for pair in pairs:
        # Select data for pair of features.
        data_pair = data_copy.loc[:, list(pair)]
        # Remove observations with missing values for either feature.
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        # Determine observations for which pair of features has matching
        # values.
        count = data_pair.shape[0]
        # Calculate correlation.
        if count > 1:
            correlation = scipy.stats.pearsonr(
                data_pair[pair[0]].values,
                data_pair[pair[1]].values,
            )[0]
            pass
        else:
            correlation = float("nan")
            pass

        # Collect information.
        if not pair[0] in correlations:
            correlations[pair[0]] = dict()
            pass
        correlations[pair[0]][pair[1]] = dict()
        correlations[pair[0]][pair[1]]["count"] = count
        correlations[pair[0]][pair[1]]["correlation"] = correlation

    # Return information.
    return correlations


def calculate_mean_correlation(
    observation,
    correlations=None,
):
    """
    Calculates the mean of pairwise correlation coefficients for multiple pairs
    of features from a single observation.

    arguments:
        observation (object): Pandas series of features from a single
            observation
        correlations (dict<dict<dict<float>>>): matrix of counts of common
            observations and correlation coefficients from pairwise comparison
            of features

    raises:

    returns:
        (float): mean of correlation coefficients

    """

    # Determine pairs of features for which observation has valid values.
    # Remove features with missing values for the observation.
    observation.dropna(
        inplace=True,
    )
    # Extract names of valid features.
    features = observation.index.values.tolist()
    # Determine pairs of features.
    pairs = list(itertools.combinations(features, 2))
    #print(features)
    #print(len(features))
    #print(pairs)
    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on pairs of features.
    for pair in pairs:
        count = correlations[pair[0]][pair[1]]["count"]
        correlation = correlations[pair[0]][pair[1]]["correlation"]
        # Collection information.
        record = dict()
        record["count"] = count
        record["correlation"] = correlation
        records.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    # Apply Fisher's z transformation to correlation coefficients.
    # Fisher's z transformation is equivalent to the inverse hyperbolic
    # tangent.
    if False:
        data["transformation"] = data["correlation"].map(
            lambda x: math.atanh(x),
        )
    data["transformation"] = data["correlation"].apply(
        lambda x: numpy.arctanh(x),
    )

    # Calculate mean across transformations of correlations.
    size = data.shape[0]
    data["numerator"] = data.aggregate(
        lambda row: ((row["count"] - 3) * row["transformation"]),
        axis="columns",
    )
    data["denominator"] = data.aggregate(
        #lambda row: (row["count"] - (3 * size)), # incorrect
        #lambda row: ((row["count"] - 3) * size), # incorrect
        lambda row: ((row["count"] - 3)),
        axis="columns",
    )
    numerator = data["numerator"].aggregate("sum")
    denominator = data["denominator"].aggregate("sum")
    mean_transformation_weight = (numerator / denominator)
    mean_transformation_simple = data["transformation"].aggregate("mean")

    # Apply inverse Fisher's z transformation to mean.
    # Inverse Fisher's z transformation is equivalent to the hyperbolic
    # tangent.
    mean_weight = math.tanh(mean_transformation_weight)
    mean_simple = math.tanh(mean_transformation_simple)

    return mean_weight


def calculate_mean_pairwise_correlations(
    correlations=None,
    data=None,
    key=None,
):
    """
    Calculates the mean pairwise correlation coefficients for multiple pairs of
    features across observations.

    Data has observations across rows and features across columns.

    index key           feature_1 feature_2 feature_3 feature_4 feature_5
    0     observation_1 ...       ...       ...       ...       ...
    1     observation_2 ...       ...       ...       ...       ...
    2     observation_3 ...       ...       ...       ...       ...
    3     observation_4 ...       ...       ...       ...       ...
    4     observation_5 ...       ...       ...       ...       ...

    arguments:
        correlations (dict<dict<dict<float>>>): matrix of counts of common
            observations and correlation coefficients from pairwise comparison
            of features
        key (str): name of column to match observations
        data (object): Pandas data frame of features across observations

    raises:

    returns:
        (object): Pandas data frame of mean correlations across observations

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_copy.set_index(
        keys=[key],
        drop=True,
        append=False,
        inplace=True,
    )
    features = data_copy.columns.values.tolist()

    # Iterate on observations.
    data_copy["correlation"] = data_copy.apply(
        calculate_mean_correlation,
        axis="columns",
        raw=False,
        correlations=correlations,
    )

    # Organize data.
    data_copy.drop(
        labels=features,
        axis="columns",
        inplace=True
    )

    # Return data.
    return data_copy


def calculate_mean_tissue_pairwise_correlations(
    data_gene_persons_tissues_signals=None
):
    """
    Calculates the mean pairwise correlation coefficients between available
    tissues across persons.

    Source data has persons across rows and tissues across columns.

    index    tissue_1 tissue_2 tissue_3 tissue_4 tissue_5
    person_1 ...      ...      ...      ...      ...
    person_2 ...      ...      ...      ...      ...
    person_3 ...      ...      ...      ...      ...
    person_4 ...      ...      ...      ...      ...
    person_5 ...      ...      ...      ...      ...

    Product data.

    index    correlation_mean
    person_1 ...
    person_2 ...
    person_3 ...
    person_4 ...
    person_5 ...

    Method
    1. determine pairs of available tissues
    2. determine counts of usable persons for each pairwise tissue comparison
    3. calculate correlation coefficients between pairs of tissues across
        persons
    4. organize counts and correlation coefficients within a tissue-tissue
        adjacency matrix
    5. calculate mean correlation coefficient for each person
    5.1. determine pairs of tissues available for each person
    5.2. collect counts and correlation coefficients for each pair of tissues
    5.3. apply Fisher's z transformation to correlation coefficients
    5.4. calculate mean Fisher's z transformation with consideration of counts
    5.5. apply inverse Fisher's z transformation

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of mean tissue-tissue correlation
            coefficients across persons

    """

    # Organize data.
    tissues = data_gene_persons_tissues_signals.columns.values.tolist()
    data_index = data_gene_persons_tissues_signals.reset_index(
        level="person",
        drop=False,
        inplace=False
    )
    # Calculate pairwise correlations between pairs of tissues.
    correlations = calculate_pairwise_correlations(
        key="person",
        data=data_index,
    )

    #utility.print_terminal_partition(level=2)
    #print(correlations["adipose"]["artery"])
    #print(correlations["artery"]["adipose"])
    #print(correlations["lung"]["heart"])
    #print(correlations["heart"]["lung"])
    #utility.print_terminal_partition(level=2)

    # Calculate mean tissue-tissue correlations across persons.
    data_gene_mean_correlations = calculate_mean_pairwise_correlations(
        correlations=correlations,
        key="person",
        data=data_index,
    )

    # Return information.
    return data_gene_mean_correlations


########## Keep these functions
# Analysis of Variance (ANOVA)


def select_values_by_binary_category(
    category=None,
    value=None,
    data=None
):
    """
    Selects a gene's aggregate scores by binary category.

    The dependent variable of interest must have column name "value".

    arguments:
        category (str): name of a category
        value (int): value of a category
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (list<float>): values for a group

    """

    # Select values of gene's scores by attribute category.
    data_selection = data.loc[
        data[category] == value, :
    ]
    return data_selection["value"].values.tolist()


def evaluate_variance_by_binary_groups(
    groups=None,
    threshold=None,
    data=None
):
    """
    Evaluates the variance between values in binary groups by one-way analysis
    of variance.

    The probability is sometimes a missing value.
    This missing value does not seem to indicate missing values in the groups.
    Potential problems
        1. a specific group lacks values
        2. multiple groups have the same values

    arguments:
        groups (list<str>): names of binary groups
        threshold (int): minimal count of values in a group for inclusion
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (dict<float>): values of results from the analysis of variance

    """

    # Collect values in each group.
    groups_values = list()
    for group in groups:
        values = select_values_by_binary_category(
            category=group,
            value=1,
            data=data,
        )
        if len(values) > threshold:
            groups_values.append(values)
    # Analyze variance.
    # Use splat operator to extract arguments from list.
    # Alternative: scipy.stats.f_oneway(*[df[col] for col in df])
    # The probability is sometimes a missing value.
    # Even when this probability is missing, there do not seem to be
    # irregularities in the values.
    if (len(groups_values) > 1):
        report = scipy.stats.f_oneway(*groups_values)
        # Compile information.
        information = {
            "statistic": report[0],
            "probability": report[1],
        }
    else:
        # Compile information.
        information = {
            "statistic": float("nan"),
            "probability": float("nan"),
        }
    # Return information.
    return information


# Write.


def write_product(
    dock=None,
    information=None
):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of a single gene
        dock (str): path to root or dock directory for source and product
            directories and files
        information (object): information to write to file

    raises:

    returns:

    """

    # Specify directories and files.
    path_prediction = os.path.join(dock, "prediction")
    utility.create_directory(path_prediction)

    path_data_regression_genes = os.path.join(
        path_prediction, "data_regression_genes.pickle"
    )
    path_data_regression_genes_text = os.path.join(
        path_prediction, "data_regression_genes.tsv"
    )

    path_data_regression_genes_selection_scale = os.path.join(
        path_prediction, "data_regression_genes_selection_scale.pickle"
    )
    path_data_regression_genes_unimodal_scale = os.path.join(
        path_prediction, "data_regression_genes_unimodal_scale.pickle"
    )
    path_data_regression_genes_multimodal_scale = os.path.join(
        path_prediction, "data_regression_genes_multimodal_scale.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_regression_genes"],
        path_data_regression_genes
    )
    information["data_regression_genes"].to_csv(
        path_or_buf=path_data_regression_genes_text,
        sep="\t",
        header=True,
        index=True,
    )

    pandas.to_pickle(
        information["data_regression_genes_selection_scale"],
        path_data_regression_genes_selection_scale
    )
    pandas.to_pickle(
        information["data_regression_genes_unimodal_scale"],
        path_data_regression_genes_unimodal_scale
    )
    pandas.to_pickle(
        information["data_regression_genes_multimodal_scale"],
        path_data_regression_genes_multimodal_scale
    )

    pass


###############################################################################
# Procedure


def execute_procedure(
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Remove previous files to avoid version or batch confusion.
    path_prediction = os.path.join(dock, "prediction")
    utility.remove_directory(path=path_prediction)

    # Read source information from file.
    source = read_source(
        dock=dock,
    )
    # Define variables for regression.
    variables = define_regression_variables()
    utility.print_terminal_partition(level=3)
    print(
        "Count of independent variables: " +
        str(len(variables["independence"]))
    )
    print(
        "Count of hypothesis (non technical) variables: " +
        str(len(variables["hypothesis"]))
    )

    # Organize dependent and independent variables for regression analysis.
    data_variables = organize_dependent_independent_variables(
        variables=variables["independence"],
        data_persons_properties=source["data_persons_properties"],
        data_signals_genes_persons=source["data_signals_genes_persons"],
    )

    # Regress each gene's signal across persons.
    # Iterate on genes.
    #genes_iteration = random.sample(source["genes_selection"], 1000)
    genes_iteration = source["genes_selection"]#[0:10]
    data_regression_genes = regress_cases(
        cases=genes_iteration,
        variables=variables["independence"],
        data_variables=data_variables,
    )
    utility.print_terminal_partition(level=2)
    print("data_regression_genes")
    print(data_regression_genes)

    # Review.
    # 4 February 2020
    # Model statistics match in the individual regression's summary and in the
    # compilation across cases.
    # Independent variables' parameters and probabilities match in the
    # individual regression's summary and in the compilation across cases.

    # Adjust probabilities for multiple hypothesis tests.
    # Calculate false discovery rate by Benjamini-Hochberg's method.
    data_regression_genes_discovery = calculate_regression_discoveries(
        variables=variables["independence"],
        threshold=0.05,
        data=data_regression_genes,
    )
    utility.print_terminal_partition(level=2)
    print("data_regression_genes_discovery")
    print(data_regression_genes_discovery)

    # Select regression parameters and probabilities that are biologically
    # relevant.
    # Scale regression parameters across genes.
    # The purpose of this scaling is to simplify comparison in chart.
    data_regression_genes_scale = select_scale_regression_parameters(
        variables=variables["hypothesis"],
        data_regression_genes=data_regression_genes,
    )
    utility.print_terminal_partition(level=2)
    print("data_regression_genes_scale")
    print(data_regression_genes_scale)

    # Separate regression information by sets of genes.
    bin = separate_regression_gene_sets(
        genes_selection=source["genes_selection"],
        genes_unimodal=source["genes_unimodal"],
        genes_multimodal=source["genes_multimodal"],
        data_regression_genes=data_regression_genes_scale,
    )
    utility.print_terminal_partition(level=2)
    print("Regression genes, original: " + str(data_regression_genes.shape[0]))
    print(
        "Regression genes, selection: " +
        str(bin["data_regression_genes_selection_scale"].shape[0])
    )
    print(
        "Regression genes, unimodal: " +
        str(bin["data_regression_genes_unimodal_scale"].shape[0])
    )
    print(
        "Regression genes, multimodal: " +
        str(bin["data_regression_genes_multimodal_scale"].shape[0])
    )

    if False:

        # Prepare summary of regression model across all genes.
        utility.print_terminal_partition(level=2)
        summarize_regression(
            instances=len(genes_iteration),
            variables=variables,
            probabilities=True,
            type="percentage", # "discovery" or "percentage"
            threshold_discovery=0.05,
            data=data_regression,
        )

    # Compile information.
    information = dict()
    information["data_regression_genes"] = data_regression_genes_discovery
    information["data_regression_genes_selection_scale"] = (
        bin["data_regression_genes_selection_scale"]
    )
    information["data_regression_genes_unimodal_scale"] = (
        bin["data_regression_genes_unimodal_scale"]
    )
    information["data_regression_genes_multimodal_scale"] = (
        bin["data_regression_genes_multimodal_scale"]
    )
    # Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
