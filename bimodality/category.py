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
import statsmodels.api

# Custom

import distribution
import metric
import plot
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source_initial(
    source_genes=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        source_genes (str): name of directory from which to obtain genes list
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    # Read information from file.
    if source_genes == "split":
        path_source = os.path.join(dock, "split")
        path_genes = os.path.join(path_source, "genes.txt")
        genes = utility.read_file_text_list(path_genes)
    elif source_genes == "rank":
        path_source = os.path.join(dock, "rank")
        path_genes = os.path.join(path_source, "genes_consensus.pickle")
        with open(path_genes, "rb") as file_source:
            genes = pickle.load(file_source)

    # Compile and return information.
    return {
        "genes": genes,
    }


def read_source(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of a single gene
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )

    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_signal = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )

    data_gene_samples_signals = pandas.read_pickle(path_signal)

    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Gene's persons and signals


def organize_persons_attributes(
    persons=None,
    data_samples_tissues_persons=None,
):
    """
    Organizes persons' attributes.

    arguments:
        persons (list<str>): list of persons of interest
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    # Copy data.
    data_persons = data_samples_tissues_persons.copy(deep=True)

    # Select persons and attributes of interest.
    data_persons.reset_index(
        level="sample",
        inplace=True
    )
    data_persons = data_persons.loc[
        data_persons["person"].isin(persons),
    ]
    data_persons.drop(
        labels=["sample", "tissue_major", "tissue_minor"],
        axis="columns",
        inplace=True
    )
    data_persons.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_persons.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    # Return data.
    return data_persons


def organize_persons_signals_groups(
    data_samples_tissues_persons=None,
    data_gene_persons_tissues=None,
    data_gene_persons_signals=None,
    data_gene_mean_correlations=None,
):
    """
    Reads and organizes source information from file

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons
        data_gene_mean_correlations (object): Pandas data frame of mean
            tissue-tissue correlation coefficients across persons

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate pan-tissue signals
            across persons with persons's attributes and groups

    """

    # Select persons of interest.
    persons = data_gene_persons_signals.index.to_list()

    # Organize information about persons.
    data_persons_attributes = organize_persons_attributes(
        persons=persons,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )

    # Organize information about persons' tissues.

    # Join
    data_correlations = data_gene_persons_signals.join(
        data_gene_mean_correlations,
        how="left",
        on="person"
    )
    data_attributes = data_correlations.join(
        data_persons_attributes,
        how="left",
        on="person"
    )
    data_tissues = data_attributes.join(
        data_gene_persons_tissues,
        how="left",
        on="person"
    )

    # Return information
    return data_tissues


def translate_groups_binary(
    tissues=None,
    data_persons_signals_groups=None,
):
    """
    Reads and organizes source information from file

    arguments:
        tissues (list<str>): identifiers of tissues
        data_persons_signals_groups (object): Pandas data frame of a gene's
            aggregate pan-tissue signals across persons with persons's
            attributes and groups

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate pan-tissue signals
            across persons with persons's attributes and groups

    """

    # Define functions
    def translate_sex_female(sex):
        if sex == "female":
            return 1
        elif sex == "male":
            return 0
    def translate_sex_male(sex):
        if sex == "male":
            return 1
        elif sex == "female":
            return 0
    def translate_tissue(tissue):
        if tissue:
            return 1
        else:
            return 0
    def translate_season(season):
        if season == "winter":
            return 0
        elif season == "spring":
            return 1
        elif season == "summer":
            return 2
        elif season == "fall":
            return 3
        else:
            return float("nan")

    # Copy data.
    data = data_persons_signals_groups.copy(deep=True)

    # Translate sex.
    data["female"] = (
        data["sex"].apply(translate_sex_female)
    )
    data["male"] = (
        data["sex"].apply(translate_sex_male)
    )
    # Translate season.
    data["season"] = (
        data["season"].apply(translate_season)
    )
    # Translate tissues.
    for tissue in tissues:
        data[tissue] = data[tissue].apply(translate_tissue)

    # Remove unnecessary columns.
    data.drop(
        labels=["sex"],
        axis="columns",
        inplace=True
    )

    # Return information
    return data


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


##########
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


##########
# Regression


def evaluate_correlation_by_linearity_separate(
    dependence=None,
    independence=None,
    threshold=None,
    data=None
):
    """
    Evaluates the correlation between a dependent and independent variable by
    linear regression.

    Returns a single probability value. Combines p-values for multiple
    independent variables.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable is a matrix, a vector of observations with
    each observation a vector of scalar values.
    StatsModels requires an intercept.
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
        independence (list<str>): name of independent variable
        threshold (float): minimal count of valid observations to perform
            regression
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (float): probability value from coefficients of the regression

    """

    if False:
        # Define regression.
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
    # Preserve observation association (bindings) between values of dependent
    # and independent variables.
    # Remove observations with any missing values.
    columns = copy.deepcopy(independence)
    columns.append(dependence)
    data = data.loc[ :, columns]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Determine whether data have sufficient observations for analysis.
    if data.shape[0] > threshold:
        # Access values of dependent and independent variables.
        values_dependence = data[dependence].values
        values_independence = data.loc[ :, independence].values
        #values_independence = statsmodels.api.add_constant(values_independence)
        # Define model.
        model = statsmodels.api.OLS(values_dependence, values_independence)
        report = model.fit()
        #print(report.summary())
        #print(dir(report))
        probabilities = report.pvalues
        if len(independence) > 1:
            probability = scipy.stats.combine_pvalues(
                probabilities,
                method="stouffer",
                weights=None,
            )[1]
        else:
            probability = probabilities[0]
    else:
        probability = float("nan")
    return probability


# TODO:
# introduce a constant for the intercept
def evaluate_correlation_by_linearity_together(
    dependence=None,
    independence=None,
    threshold=None,
    intercept=None,
    data=None
):
    """
    Evaluates the correlation between a dependent and independent variable by
    linear regression.

    Returns a list of probability values for all independent variables.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable is a matrix, a vector of observations with
    each observation a vector of scalar values.
    StatsModels requires an intercept.
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
        independence (list<str>): name of independent variable
        threshold (float): minimal count of valid observations to perform
            regression
        intercept (bool): whether to introduce a constant for intercept
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (float): probability value from coefficients of the regression

    """

    if False:
        # Define regression.
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
    # Preserve observation association (bindings) between values of dependent
    # and independent variables.
    # Remove observations with any missing values.
    columns = copy.deepcopy(independence)
    columns.append(dependence)
    data = data.loc[ :, columns]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Determine whether data have sufficient observations for analysis.
    if data.shape[0] > threshold:
        # Access values of dependent and independent variables.
        values_dependence = data[dependence].values
        values_independence = data.loc[ :, independence].values

        # https://stackoverflow.com/questions/38836465/how-to-get-the-regression-intercept-using-statsmodels-api

        values_independence_intercept = statsmodels.api.add_constant(
            values_independence,
            prepend=False,
        )

        #values_independence = statsmodels.api.add_constant(values_independence)
        # Define model.
        model = statsmodels.api.OLS(
            values_dependence, values_independence_intercept
        )
        report = model.fit()
        #print(report.summary())
        #print(dir(report))
        probabilities = report.pvalues
    else:
        probabilities = list(
            itertools.repeat(float("nan"), len(values_independence_intercept))
        )
    return probabilities


##########
# Organization and evaluation


# This function calls ANOVA and regression procedures.
# TODO: include pairwise Pearson correlation p-values too? <- maybe same as regression p-values
def evaluate_gene_distribution_covariance(
    gene=None,
    sexes=None,
    tissues=None,
    data_persons_signals_groups=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        sexes (list<str>): identifiers of sexes
        tissues (list<str>): identifiers of tissues
        data_persons_signals_groups (object): Pandas data frame of a gene's
            aggregate pan-tissue signals across persons with persons's
            attributes and groups

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Evaluate covariance by categories.
    # Function return probability values from ANOVA by sex and tissue
    # categories.
    category = evaluate_gene_distribution_category(
        gene=gene,
        sexes=sexes,
        tissues=tissues,
        threshold_proportion=0.25,
        data_persons_signals_groups=data_persons_signals_groups,
    )

    # Evaluate covariance by correlation.

    # Evaluate covariance by regression.
    regression = evaluate_gene_distribution_regression_together(
        gene=gene,
        dependence="value",
        tissues=tissues,
        threshold_proportion=0.25,
        combination="minimum",
        data_persons_signals_groups=data_persons_signals_groups,
    )

    if False:
        regression = evaluate_gene_distribution_regression(
            gene=gene,
            dependence="correlation",
            tissues=tissues,
            threshold_proportion=0.25,
            data_persons_signals_groups=data_persons_signals_groups,
        )

    # Compile information.
    information = {
        "category_tissue": category["tissue"],
        "category_sex": category["sex"],
        "regression_ancestry": regression["ancestry"],
        "regression_sex": regression["sex"],
        "regression_age": regression["age"],
        "regression_body": regression["body"],
        "regression_season": regression["season"],
        "regression_hardiness": regression["hardiness"],
        "regression_delay": regression["delay"],
        "regression_correlation": regression["correlation"],
        "regression_tissue": regression["tissue"],
    }

    # Return information.
    return information


def evaluate_gene_distribution_category(
    gene=None,
    sexes=None,
    tissues=None,
    threshold_proportion=None,
    data_persons_signals_groups=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        sexes (list<str>): identifiers of sexes
        tissues (list<str>): identifiers of tissues
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a group in analyses
        data_persons_signals_groups (object): Pandas data frame of a gene's
            aggregate pan-tissue signals across persons with persons's
            attributes and groups

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Determine minimal count of persons for inclusion of a group in
    # analyses.
    count_persons = data_persons_signals_groups.shape[0]
    threshold_count = threshold_proportion * count_persons

    # Analysis of Variance (ANOVA)

    # Signals versus tissues
    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    report = evaluate_variance_by_binary_groups(
        groups=tissues,
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )
    tissue = report["probability"]

    # Signals versus sex
    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    report = evaluate_variance_by_binary_groups(
        groups=sexes,
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )
    sex = report["probability"]

    # Compile information.
    information = {
        "tissue": tissue,
        "sex": sex,
    }
    # Return information.
    return information


def evaluate_gene_distribution_regression_together(
    gene=None,
    dependence=None,
    tissues=None,
    threshold_proportion=None,
    combination=None,
    data_persons_signals_groups=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        dependence (str): identifier of dependent variable
        tissues (list<str>): identifiers of tissues
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a group in analyses
        combination (str): method by which to combine probabilities of multiple
            covariates, either "combination" or "minimum"
        data_persons_signals_groups (object): Pandas data frame of a gene's
            aggregate pan-tissue signals across persons with persons's
            attributes and groups

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Determine minimal count of persons for inclusion of a group in
    # analyses.
    count_persons = data_persons_signals_groups.shape[0]
    threshold_count = threshold_proportion * count_persons

    # Regression

    # Signals versus mean tissue-tissue correlation
    if dependence == "value":
        probabilities = evaluate_correlation_by_linearity_together(
            dependence=dependence,
            independence=[
                "ancestry_1", "ancestry_2", "ancestry_3", # x1, x2, x3,
                "female", # x4
                "age", # x5
                "body", # x6
                "season", # x7
                "hardiness", # x8
                "delay", # x9
                "correlation", # x10
                *tissues, # x11 ...
            ],
            threshold=threshold_count,
            intercept=True,
            data=data_persons_signals_groups,
        )
    elif dependence == "correlation":
        probabilities = evaluate_correlation_by_linearity(
            dependence=dependence,
            independence=["value"],
            threshold=threshold_count,
            data=data_persons_signals_groups,
        )

    # Determine combination probabilities.
    if combination == "combination":
        # Combine probabilities by Stouffer's method.
        ancestry = scipy.stats.combine_pvalues(
            probabilities[0:3],
            method="stouffer",
            weights=None,
        )[1]
        tissue = scipy.stats.combine_pvalues(
            probabilities[10:-1],
            method="stouffer",
            weights=None,
        )[1]
    elif combination == "minimum":
        # Take the minimal probability.
        ancestry = min(probabilities[0:3])
        tissue = min(probabilities[10:-1])

    # Compile information.
    information = {
        "ancestry": ancestry,
        "sex": probabilities[3],
        "age": probabilities[4],
        "body": probabilities[5],
        "season": probabilities[6],
        "hardiness": probabilities[7],
        "delay": probabilities[8],
        "correlation": probabilities[9],
        "tissue": tissue,
    }
    # Return information.
    return information


def evaluate_gene_distribution_regression_separate(
    gene=None,
    dependence=None,
    tissues=None,
    threshold_proportion=None,
    data_persons_signals_groups=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        dependence (str): identifier of dependent variable
        tissues (list<str>): identifiers of tissues
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a group in analyses
        data_persons_signals_groups (object): Pandas data frame of a gene's
            aggregate pan-tissue signals across persons with persons's
            attributes and groups

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Determine minimal count of persons for inclusion of a group in
    # analyses.
    count_persons = data_persons_signals_groups.shape[0]
    threshold_count = threshold_proportion * count_persons

    # Regression

    # Signals versus mean tissue-tissue correlation
    if dependence == "value":
        cross = evaluate_correlation_by_linearity(
            dependence=dependence,
            independence=["correlation"],
            threshold=threshold_count,
            data=data_persons_signals_groups,
        )
    elif dependence == "correlation":
        cross = evaluate_correlation_by_linearity(
            dependence=dependence,
            independence=["value"],
            threshold=threshold_count,
            data=data_persons_signals_groups,
        )

    # Sex
    sex = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["female"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Age
    age = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["age"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Body mass index
    body = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["body"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Ancestry # TODO: I think I'll need to combine the p-values?
    ancestry = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["ancestry_1", "ancestry_2", "ancestry_3"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Tissue
    tissue = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=tissues,
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Hardiness
    hardiness = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["hardiness"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Season
    season = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["season"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Delay
    delay = evaluate_correlation_by_linearity(
        dependence=dependence,
        independence=["delay"],
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )

    # Compile information.
    information = {
        "cross": cross,
        "sex": sex,
        "age": age,
        "body": body,
        "ancestry": ancestry,
        "tissue": tissue,
        "hardiness": hardiness,
        "season": season,
        "delay": delay,
    }
    # Return information.
    return information


# This function organizes each gene's signals across persons and groups.
# TODO: include information from reports from organization (max variation) and restriction procedures...
# TODO: consider including some of those values in regression...
def evaluate_gene_categories_method(
    gene=None,
    method=None,
    data_samples_tissues_persons=None,
    observation=None,
    data_gene_annotation=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Access information.
    report_organization = observation["organization"]["report_gene"]
    report_restriction = observation[method]["report_restriction"]
    report_aggregation = observation[method]["report_aggregation"]
    data_gene_persons_tissues = report_restriction["data_gene_persons_tissues"]
    data_gene_persons_signals = report_aggregation["data_gene_persons_signals"]
    persons = report_aggregation["persons"]
    tissues_mean = report_restriction["tissues_mean"]
    tissues_median = report_restriction["tissues_median"]

    name = data_gene_annotation.loc[gene, "gene_name"]

    # Determine derivative information about gene's signals across persons.
    # Calculate tissue-tissue correlations across persons.
    # Use the data before imputation to grasp correlations across real values.
    data_gene_mean_correlations = (
        calculate_mean_tissue_pairwise_correlations(
            data_gene_persons_tissues_signals=(
                report_restriction["data_gene_persons_tissues_signals"]
            )
        )
    )
    # TODO:
    # Calculate principal components on binary values for tissue selection.
    # Observations are across persons.
    # Features are across tissues available for each person.
    # Values could be signals for each tissue, or they could be binary
    # variables to indicate availability for each person.

    # Organize information about gene's signals across persons and groups.
    data_persons_signals_groups = organize_persons_signals_groups(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=data_gene_persons_tissues,
        data_gene_persons_signals=data_gene_persons_signals,
        data_gene_mean_correlations=data_gene_mean_correlations,
    )

    # Translate categorical variables for analysis.
    tissues = data_gene_persons_tissues.columns.values.tolist()
    data_persons_signals_groups_binary = translate_groups_binary(
        tissues=tissues,
        data_persons_signals_groups=data_persons_signals_groups,
    )

    #print(data_persons_signals_groups_binary)

    # Evaluate gene's distribution by covariates.
    report = evaluate_gene_distribution_covariance(
        gene=gene,
        sexes=["male", "female"],
        tissues=tissues,
        data_persons_signals_groups=data_persons_signals_groups_binary,
    )
    # Include more information.
    report.update(
        gene=gene,
        name=name,
        persons=persons,
        tissues_mean=tissues_mean,
        tissues_median=tissues_median,
    )

    # Compile information.
    information = {
        "report": report,
        "data_persons_signals_groups": data_persons_signals_groups,
    }

    # Return information.
    return information


def evaluate_gene_categories(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
    data_gene_annotation=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Determine gene's distributions of aggregate tissue scores across persons.
    observation = distribution.determine_gene_distributions(
        gene=gene,
        modality=False,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    imputation = evaluate_gene_categories_method(
        gene=gene,
        method="imputation",
        data_samples_tissues_persons=data_samples_tissues_persons,
        observation=observation,
        data_gene_annotation=data_gene_annotation,
    )

    availability = evaluate_gene_categories_method(
        gene=gene,
        method="availability",
        data_samples_tissues_persons=data_samples_tissues_persons,
        observation=observation,
        data_gene_annotation=data_gene_annotation,
    )

    # Compile information.
    information = {
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


# This function drives iteration across genes. It is now obsolete.
def evaluate_genes_categories(
    genes=None,
    data_samples_tissues_persons=None,
    genes_samples_signals=None,
    data_gene_annotation=None,
):
    """
    Evaluates a gene's distribution of signals across persons by multiple
    factors.

    arguments:
        genes (list<str>): identifiers of genes
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        genes_samples_signals (dict<object>): collection of Pandas data frames
            of genes' signals across samples
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (dict): information about gene's distribution and categories

    """

    # Initiate collections.
    records_availability = list()
    records_imputation = list()
    # Iterate on genes.
    for gene in genes:
        collection = evaluate_gene_categories(
            gene=gene,
            data_samples_tissues_persons=data_samples_tissues_persons,
            genes_samples_signals=genes_samples_signals,
            data_gene_annotation=data_gene_annotation,
        )
        records_availability.append(collection["availability"])
        records_imputation.append(collection["imputation"])
    # Organize data.
    data_availability = utility.convert_records_to_dataframe(
        records=records_availability
    )
    data_availability.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )

    data_imputation = utility.convert_records_to_dataframe(
        records=records_imputation
    )
    if False:
        data_imputation.set_index(
            ["gene"],
            append=False,
            drop=True,
            inplace=True
        )

    # Compile information.
    information = {
        "imputation": data_imputation,
        "availability": data_availability,
    }

    # Return information.
    return information


# Write.


def initialize_directories(
    dock=None,
):
    """
    Initializes directories.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Create directories.
    path_category = os.path.join(dock, "category")
    path_collection = os.path.join(path_category, "collection")

    utility.remove_directory(path=path_category)

    utility.create_directory(path_category)
    utility.create_directory(path_collection)

    pass


def write_product(
    gene=None,
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
    path_category = os.path.join(dock, "category")
    path_collection = os.path.join(path_category, "collection")
    path_gene = os.path.join(path_collection, gene)
    utility.create_directory(path_gene)

    path_imputation = os.path.join(path_gene, "report_imputation.pickle")
    path_availability = os.path.join(path_gene, "report_availability.pickle")

    # Write information to file.
    with open(path_imputation, "wb") as file_product:
        pickle.dump(information["report_imputation"], file_product)
    with open(path_availability, "wb") as file_product:
        pickle.dump(information["report_availability"], file_product)

    pass


def read_collect_genes_reports(
    genes=None,
    method=None,
    dock=None,
):
    """
    Collects information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<dict<dict>>): information about genes

    """

    # Specify directories and files.
    path_category = os.path.join(dock, "category")
    path_collection = os.path.join(path_category, "collection")

    # Report process.
    utility.print_terminal_partition(level=1)
    print(
        "Reading and compiling information about genes' distributions."
    )
    # Check contents of directory.
    print("Check that directories exist for all genes.")
    directories = os.listdir(path_collection)
    match = utility.compare_lists_by_mutual_inclusion(
        list_one=genes, list_two=directories
    )
    print("Genes and directories match: " + str(match))

    # Collect information about genes.
    records = list()
    # Iterate on directories for genes.
    for gene in directories:
        # Report progress.
        #print(gene)
        # Specify directories and files.
        path_gene = os.path.join(path_collection, gene)
        path_report = os.path.join(
            path_gene, ("report_" + method + ".pickle")
        )
        # Read information from file.
        with open(path_report, "rb") as file_source:
            record = pickle.load(file_source)
        # Create entry for gene.
        records.append(record)

    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    data.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )

    # Report genes with distributions dependent on each covariate.

    # Specify directories and files.
    path_data = os.path.join(path_category, ("data_" + method + ".tsv"))

    # Write information to file.
    data.to_csv(
        path_or_buf=path_data,
        sep="\t",
        header=True,
        index=False,
    )

    pass




###############################################################################
# Procedure


# TODO: iterate on a list of genes from the rank procedure in a parallel way
# TODO: determine the distribution for each gene's signals
# TODO: evaluate the gene's signal distribution by ANOVA and regression... include attributes of patients
# - sex, age, race, ancestry, time to harvest, hardiness scale, etc...
# construct a mixed effects model


def execute_procedure(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
    data_gene_annotation=None,
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a single gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        data_gene_annotation (object): Pandas data frame of genes' annotations
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Evaluate genes' categories and distributions.
    # Report.
    collection = evaluate_gene_categories(
        gene=gene,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
        data_gene_annotation=data_gene_annotation,
    )

    # Compile information.
    information = {
        "report_imputation": collection["imputation"]["report"],
        "report_availability": collection["availability"]["report"],
    }
    #Write product information to file.
    write_product(gene=gene, dock=dock, information=information)

    pass


def execute_procedure_local(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    initialize_directories(dock=dock)

    # Read source information from file.
    source = read_source_initial(
        source_genes="split",
        dock=dock
    )
    print("count of genes: " + str(len(source["genes"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "shuffles", and
    # "dock" variables.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
        dock=dock,
    )

    # Initialize multiprocessing pool.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=8)

    # Iterate on genes.
    genes_check = [
        "ENSG00000231925", # TAPBP ... tapasin binding protein
        #"ENSG00000147050", # KDM6A ... X-linked analog to UTY
        #"ENSG00000183878", # UTY ... Y-linked analog to KDM6A
    ]
    #report = pool.map(execute_procedure_gene, genes_check)
    #report = pool.map(execute_procedure_gene, source["genes"][0:100])
    report = pool.map(execute_procedure_gene, source["genes"])

    # Pause procedure.
    time.sleep(5.0)

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # Collect genes.
    read_collect_genes_reports(
        genes=source["genes"],
        method="imputation",
        dock=dock,
    )
    read_collect_genes_reports(
        genes=source["genes"],
        method="availability",
        dock=dock,
    )

    # Report date and time.
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))

    pass


def execute_procedure_local_sub(
    gene=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a single gene
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        data_gene_annotation=source["data_gene_annotation"],
        dock=dock,
    )

    # Report contents of directory.
    path_category = os.path.join(dock, "category")
    path_collection = os.path.join(path_category, "collection")
    files = os.listdir(path_collection)
    count = len(files)
    if (count % 10 == 0):
        print("complete genes: " + str(len(files)))

    pass


if (__name__ == "__main__"):
    execute_procedure()
