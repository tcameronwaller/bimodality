"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import functools
import multiprocessing
import datetime
import gc
import time
import random

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

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
    path_selection = os.path.join(dock, "selection")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_split_genes = os.path.join(path_split, "genes.txt")

    # Use genes' bimodality measures from distribution procedure (not those
    # from the probability procedure) because here we want the raw values
    # before standardization.
    path_distribution = os.path.join(dock, "distribution")
    path_genes_scores = os.path.join(
        path_distribution, "genes_scores.pickle"
    )
    path_scores = os.path.join(
        path_distribution, "scores.pickle"
    )
    path_data_report = os.path.join(
        path_distribution, "data_gene_report.pickle"
    )
    path_data_report_text = os.path.join(
        path_distribution, "data_gene_report.tsv"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    genes_split = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_split_genes,
    )
    genes_distribution = utility.extract_subdirectory_names(
        path=path_distribution
    )
    with open(path_genes_scores, "rb") as file_source:
        genes_scores = pickle.load(file_source)
    with open(path_scores, "rb") as file_source:
        scores = pickle.load(file_source)
    data_distribution_report = pandas.read_pickle(path_data_report)

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_split": genes_split,
        "genes_distribution": genes_distribution,
        "genes_scores": genes_scores,
        "scores": scores,
        "data_distribution_report": data_distribution_report,
    }


def determine_bimodality_measures_thresholds(
    measures=None,
    scores=None,
    factors=None,
):
    """
    Determines values of thresholds for genes' measures of bimodality.

    arguments:
        measures (list<str>): measures of bimodality
        scores (dict<dict>): information about genes' measures of bimodality
        factors (dict<float>): factors for each measure of bimodality

    raises:

    returns:
        (dict<float>): values of thresholds for genes' measures of bimodality

    """

    # Collect thresholds for measures of bimodality.
    thresholds = dict()
    for measure in measures:
        # Calculate selection threshold.
        thresholds[measure] = (
            scores[measure]["mean"] +
            (factors[measure] * scores[measure]["deviation"])
        )
    pass
    return thresholds


def summarize_bimodality_measures_thresholds(
    measures=None,
    scores=None,
    thresholds=None,
):
    """
    Summarizes values of thresholds for genes' measures of bimodality.

    arguments:
        measures (list<str>): measures of bimodality
        scores (dict<dict>): information about genes' measures of bimodality
        thresholds (dict<float>): values of thresholds for genes' measures of
            bimodality

    raises:

    returns:


    """

    utility.print_terminal_partition(level=2)
    for measure in measures:
        utility.print_terminal_partition(level=3)
        print("measure: " + str(measure))
        print("mean: " + str(scores[measure]["mean"]))
        print("deviation: " + str(scores[measure]["deviation"]))
        print("threshold: " + str(thresholds[measure]))
    pass
    utility.print_terminal_partition(level=2)


# Filtration


def copy_split_minimal_gene_data(
    measure=None,
    data_genes_distributions=None,
):
    """
    Copy and split information about genes.

    arguments:
        measure (str): name of a measure of bimodality
        data_genes_distributions (object): Pandas data frame of information
            about genes and their measures of bimodality

    raises:

    returns:
        (object): Pandas data frame of genes' identifiers and values of a
            measure of bimodality

    """

    # Select relevant information.
    columns = list()
    columns.append("identifier")
    columns.append(measure)
    data_copy = data_genes_distributions.copy(deep=True)
    data = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    # Return information.
    return data


def filter_genes_by_bimodality_threshold(
    data=None,
    measure=None,
    threshold=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality measures with permutation probabilities below threshold.

    arguments:
        data (object): Pandas data frame of information about genes and their
            measures of bimodality
        measure (str): name of a measure of bimodality
        threshold (float): value of a threshold for a measure of bimodality

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities

    """

    # Determine whether values pass threshold.
    # Only consider the original modality measures for this threshold.
    data_threshold = data.loc[(data[measure] > threshold), :]
    # Return information.
    return data_threshold


# TODO: this function is useful...
def filter_genes_probabilities_threshold_old(
    data=None,
    threshold=None,
    count=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality measures with permutation probabilities below threshold.

    arguments:
        data (object): Pandas data frame of probabilities from genes' scores
            and permutations
        threshold (float): maximal probability
        count (int): minimal count of probabilities that must pass threshold

    raises:

    returns:
        (object): Pandas data frame of probabilities from genes' scores and
            distributions

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Count how many of the gene's modality probabilities are below threshold
    # Keep genes with at least 2 measures below threshold

    # Determine whether values pass threshold.
    # Only consider the original modality measures for this threshold.
    data_selection = data.loc[ :, ["coefficient", "dip", "mixture"]]
    data_threshold = (data_selection <= threshold)
    # This aggregation operation produces a series.
    # Aggregate columns to determine which rows to keep in the next step.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=2),
        axis="columns",
    )

    # Select rows and columns with appropriate values.
    data_pass = data.loc[data_count, : ]
    # Return information.
    return data_pass


def filter_genes_by_bimodality_thresholds(
    measures=None,
    thresholds=None,
    data_genes_distributions=None,
):
    """
    Copy and split information about genes.

    arguments:
        measures (list<str>): measures of bimodality
        thresholds (dict<float>): values of thresholds for measures of
            bimodality
        data_genes_distributions (object): Pandas data frame of information
            about genes and their measures of bimodality

    raises:

    returns:
        (dict<list<str>>): identifiers of genes that pass filtration by
            thresholds on each measure of bimodality

    """

    utility.print_terminal_partition(level=1)
    print(
        "count of genes filtered by probabilities of each bimodality " +
        "measurement"
    )

    # Collect genes from filtration by each measurement of bimodality.
    entries = dict()
    for measure in measures:
        # Copy minimal genes' data for each measure of bimodality.
        data_measure = copy_split_minimal_gene_data(
            measure=measure,
            data_genes_distributions=data_genes_distributions,
        )

        # Filter genes by threshold on each measure's probabilities.
        data_filter = filter_genes_by_bimodality_threshold(
            data=data_measure,
            measure=measure,
            threshold=thresholds[measure],
        )

        # Extract genes' identifiers.
        genes = data_filter["identifier"].tolist()
        utility.print_terminal_partition(level=3)
        print(measure + ": " + str(len(genes)))

        # Compile information.
        entries[measure] = genes

    # Return information.
    return entries


def select_genes_by_identifier_rank(
    data_genes=None,
    genes=None,
    rank=None,
):
    """
    Prepares ranks of genes for subsequent functional analyses.

    arguments:
        data_genes (object): Pandas data frame of information about genes
        genes (list<str>): identifiers of genes
        rank (str): property to use for ranks

    raises:

    returns:
        (object): Pandas data frame of information about genes

    """

    # Copy data.
    data_copy = data_genes.copy(deep=True)
    # Select data for genes of interest.
    data = data_copy.loc[
        data_copy["identifier"].isin(genes), :
    ]
    # Rank genes.
    data.sort_values(
        by=[rank],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Organize information.
    # Return information.
    return data


##########
# Product


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
    path_candidacy = os.path.join(dock, "candidacy")
    utility.create_directory(path_candidacy)
    path_measures_thresholds = os.path.join(
        path_candidacy, "measures_thresholds.pickle"
    )
    path_genes_sets = os.path.join(
        path_candidacy, "sets_genes_measures.pickle"
    )
    path_genes_candidacy = os.path.join(
        path_candidacy, "genes_candidacy.pickle"
    )
    path_data_genes_candidacy = os.path.join(
        path_candidacy, "data_genes_candidacy.pickle"
    )
    path_data_genes_candidacy_text = os.path.join(
        path_candidacy, "data_genes_candidacy.tsv"
    )
    # Write information to file.
    with open(path_measures_thresholds, "wb") as file_product:
        pickle.dump(information["measures_thresholds"], file_product)
    with open(path_genes_sets, "wb") as file_product:
        pickle.dump(information["sets_genes_measures"], file_product)
    with open(path_genes_candidacy, "wb") as file_product:
        pickle.dump(
            information["genes_candidacy"], file_product
        )
    information["data_genes_candidacy"].to_pickle(
        path=path_data_genes_candidacy
    )
    information["data_genes_candidacy"].to_csv(
        path_or_buf=path_data_genes_candidacy_text,
        sep="\t",
        header=True,
        index=True,
    )

    pass



###############################################################################
# Procedure


def execute_procedure(
    dock=None
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
    path_candidacy = os.path.join(dock, "candidacy")
    utility.remove_directory(path=path_candidacy)

    # Read source information from file.
    source = read_source(dock=dock)

    print(source["data_distribution_report"])

    # Set measures of bimodality.
    measures = ["dip", "mixture", "coefficient"]
    # Set factors for each measure of bimodality.
    # Set factors to select 50-100 genes by each measure of bimodality.
    # Values of bimodality coefficient greater than 0.555 are multimodal.
    factors = dict()
    factors["dip"] = 2 # (mean + (3 * sigma) = 0.021): 55
    factors["mixture"] = 3 # (mean + (5 * sigma) = 0.273): 62
    factors["coefficient"] = 3 # (mean + (4 * sigma) = 0.574): 68
    # Calculate thresholds for each measure of bimodality.
    measures_thresholds = determine_bimodality_measures_thresholds(
        measures=measures,
        scores=source["scores"],
        factors=factors,
    )
    summarize_bimodality_measures_thresholds(
        measures=measures,
        scores=source["scores"],
        thresholds=measures_thresholds,
    )

    # Filter genes by thresholds on measures of bimodality.
    genes_measures = filter_genes_by_bimodality_thresholds(
        measures=measures,
        thresholds=measures_thresholds,
        data_genes_distributions=source["data_distribution_report"],
    )

    # Select genes that pass filters by multiple measures of bimodality.
    genes_candidacy = utility.select_elements_by_sets(
        names=measures,
        sets=genes_measures,
        count=2,
    )
    utility.print_terminal_partition(level=2)
    print(
        "selection of genes by multiple measurements: " +
        str(len(genes_candidacy))
    )
    utility.print_terminal_partition(level=2)

    # Select and rank genes by probabilities.
    data_genes_candidacy = select_genes_by_identifier_rank(
        data_genes=source["data_distribution_report"],
        genes=genes_candidacy,
        rank="combination",
    )
    #print(data_genes_candidacy)
    #print(data_genes_candidacy.iloc[0:10, 0:13])

    # Compile information.
    information = {
        "measures_thresholds": measures_thresholds,
        "sets_genes_measures": genes_measures,
        "genes_candidacy": genes_candidacy,
        "data_genes_candidacy": data_genes_candidacy,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
