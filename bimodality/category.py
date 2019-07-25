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

# Relevant

import numpy
import pandas
import scipy.stats
from sklearn.linear_model import LinearRegression

# Custom

import pipe
import metric
import plot
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(
    gene=None,
    method=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of a gene
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
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
    path_pipe = os.path.join(dock, "pipe")
    path_gene = os.path.join(path_pipe, gene)
    path_method = os.path.join(path_gene, method)

    path_report_restriction = os.path.join(
        path_method, "report_restriction.pickle"
    )
    path_report_distribution = os.path.join(
        path_method, "report_distribution.pickle"
    )

    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    with open(path_report_restriction, "rb") as file_source:
        report_restriction = pickle.load(file_source)
    with open(path_report_distribution, "rb") as file_source:
        report_distribution = pickle.load(file_source)

    # Access information.
    data_persons_tissues = report_restriction["data_persons_tissues"]
    data_gene_persons_signals = (
        report_distribution["data_gene_persons_signals"]
    )

    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_persons_tissues": data_persons_tissues,
        "data_gene_persons_signals": data_gene_persons_signals,
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
    def translate_age(age):
        ages_text = [
            "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89",
            "90-99", "100-109",
        ]
        ages_number = [
            25, 35, 45, 55, 65, 75, 85, 95, 105
        ]
        # Match indices.
        if age in ages_text:
            return ages_number[ages_text.index(age)]
        else:
            return float("nan")

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

    # Organize persons' sex and age.
    data_persons["female"] = (
        data_persons["sex"].apply(translate_sex_female)
    )
    data_persons["male"] = (
        data_persons["sex"].apply(translate_sex_male)
    )
    data_persons["age"] = (
        data_persons["age_range"].apply(translate_age)
    )
    data_persons.drop(
        labels=["sex", "age_range"],
        axis="columns",
        inplace=True
    )
    # Return data.
    return data_persons


def organize_persons_tissues(
    data_gene_persons_tissues=None,
):
    """
    Organizes persons' tissues.

    arguments:
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons

    raises:

    returns:
        (object): Pandas data frame of a gene's selection of tissues across
            persons

    """

    # Define functions
    def translate_tissue(tissue):
        if tissue:
            return 1
        else:
            return 0

    # Copy data.
    data_persons_tissues = data_gene_persons_tissues.copy(deep=True)

    # Organize persons' tissues.
    data_persons_tissues = (
        data_persons_tissues.applymap(translate_tissue)
    )
    # Return data.
    return data_persons_tissues


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
    data_gene_persons_tissues_binary = organize_persons_tissues(
        data_gene_persons_tissues=data_gene_persons_tissues
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
        data_gene_persons_tissues_binary,
        how="left",
        on="person"
    )

    # Return information
    return data_tissues


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


def evaluate_correlation_by_linearity(
    dependence=None,
    independence=None,
    data=None
):
    """
    Evaluates the correlation between a dependent and independent variable by
    linear regression.

    arguments:
        dependence (str): name of dependent variable
        independence (list<str>): name of independent variable
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (dict<float>): values of results from the analysis of variance

    """

    # Access values of dependent and independent variables.
    values_dependence = data[dependence].values.reshape(-1, 1)
    if len(independence) > 1:
        data_independence = data.loc[ :, independence]
        values_independence = data_independence.values
        pass
    else:
        values_independence = data[independence[0]].values.reshape(-1, 1)

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
    return score


# Write.


def write_product(dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_analysis = os.path.join(dock, "analysis")
    utility.confirm_path_directory(path_analysis)
    path_probabilities = os.path.join(
        path_analysis, "data_genes_probabilities.pickle"
    )
    path_summary = os.path.join(
        path_analysis, "data_summary_genes.pickle"
    )
    path_summary_text = os.path.join(
        path_analysis, "data_summary_genes.txt"
    )
    path_summary_text_alternative = os.path.join(
        path_analysis, "data_summary_genes_alternative.txt"
    )
    path_data_rank = os.path.join(
        path_analysis, "data_rank_genes.txt"
    )
    path_rank_ensembl = os.path.join(
        path_analysis, "genes_ranks_ensembl.txt"
    )
    path_rank_hugo = os.path.join(
        path_analysis, "genes_ranks_hugo.txt"
    )
    path_genes_ensembl = os.path.join(
        path_analysis, "genes_ensembl.txt"
    )
    path_genes_hugo = os.path.join(
        path_analysis, "genes_hugo.txt"
    )
    path_genes_report = os.path.join(
        path_analysis, "genes_report.txt"
    )
    path_reports = os.path.join(
        path_analysis, "reports.pickle"
    )
    # Write information to file.
    information["data_genes_probabilities"].to_pickle(path_summary)
    information["data_summary_genes"].to_pickle(path_summary)
    information["data_summary_genes"].to_csv(
        path_or_buf=path_summary_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_summary_genes"].reset_index(
        level="identifier", inplace=True
    )
    summary_genes = utility.convert_dataframe_to_records(
        data=information["data_summary_genes"]
    )
    utility.write_file_text_table(
        information=summary_genes,
        path_file=path_summary_text_alternative,
        names=summary_genes[0].keys(),
        delimiter="\t"
    )
    information["data_rank_genes"].to_csv(
        path_or_buf=path_data_rank,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_rank_genes_ensembl"].to_csv(
        path_or_buf=path_rank_ensembl,
        sep="\t",
        header=False,
        index=False,
    )
    information["data_rank_genes_hugo"].to_csv(
        path_or_buf=path_rank_hugo,
        sep="\t",
        header=False,
        index=False,
    )
    utility.write_file_text_list(
        information=information["genes_ensembl"],
        path_file=path_genes_ensembl
    )
    utility.write_file_text_list(
        information=information["genes_hugo"],
        path_file=path_genes_hugo
    )
    utility.write_file_text_list(
        information=information["genes_report"],
        path_file=path_genes_report
    )
    with open(path_reports, "wb") as file_product:
        pickle.dump(information["reports"], file_product)

    pass


def write_product_reports(dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_analysis = os.path.join(dock, "analysis")
    utility.confirm_path_directory(path_analysis)
    path_reports = os.path.join(path_analysis, "reports")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_reports)
    utility.confirm_path_directory(path_reports)
    # Iterate on reports.
    for gene in information["genes_report"]:
        # Access information.
        name = information["reports"][gene]["name"]
        # Specify directories and files.
        path_report = os.path.join(
            path_reports, (name + "_reports.pickle")
        )
        # Write information to file.

        pass

    pass



# TODO: write_product_reports()


###############################################################################
# Procedure


def execute_procedure(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Specify genes of interest.
    gene = "ENSG00000231925" # TAPBP
    #gene = "ENSG00000183878" # UTY ... Y-linked analog to KDM6A
    #gene = "ENSG00000147050" # KDM6A ... X-linked analog to UTY

    # Read source information from file.
    source = read_source(
        gene=gene,
        method="availability",
        dock=dock
    )

    # Organize persons' signals and properties.

    if True:
        data_independence = data_persons_tissues_signals.drop(
            labels=["value", "age", "sex"],
            axis="columns",
            inplace=False
        )
        pass
    if False:
        data_independence = data_persons_tissues_signals.loc[
            :, "age"
        ]
        pass
    #print(data_independence.values.reshape(-1, 1))
    data_dependence = data_persons_tissues_signals.loc[
        :, "value"
    ]
    #print(data_dependence.values.reshape(-1, 1))

    # Regression...
    regression = LinearRegression(
        fit_intercept=True,
        normalize=False,
        copy_X=True,
        n_jobs=None
    ).fit(
        data_independence.values,#.reshape(-1, 1),
        data_dependence.values.reshape(-1, 1),
        sample_weight=None
    )
    score = regression.score(
        data_independence.values,#.reshape(-1, 1),
        data_dependence.values.reshape(-1, 1),
        sample_weight=None
    )
    print(score)




    if False:

        # Compile information.
        information = {
            "data_genes_probabilities": data_genes_probabilities,
            "data_summary_genes": data_summary_genes,
            "data_rank_genes": ranks["data_rank_genes"],
            "data_rank_genes_ensembl": ranks["data_rank_genes_ensembl"],
            "data_rank_genes_hugo": ranks["data_rank_genes_hugo"],
            "genes_ensembl": ranks["genes_ensembl"],
            "genes_hugo": ranks["genes_hugo"],
            "genes_report": genes_report,
            "reports": reports,
        }
        #Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
