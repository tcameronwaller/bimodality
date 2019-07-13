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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import rank
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(dock=None):
    """
    Reads and organizes source information from file

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    r
    aises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_collection = os.path.join(dock, "collection")
    path_scores_availability = os.path.join(
        path_collection, "genes_scores_shuffles_availability.pickle"
    )
    path_scores_imputation = os.path.join(
        path_collection, "genes_scores_shuffles_imputation.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    with open(path_scores_availability, "rb") as file_source:
        genes_scores_shuffles_availability = pickle.load(file_source)
    with open(path_scores_imputation, "rb") as file_source:
        genes_scores_shuffles_imputation = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "genes_scores_shuffles_availability": (
            genes_scores_shuffles_availability
        ),
        "genes_scores_shuffles_imputation": (
            genes_scores_shuffles_imputation
        ),
    }


# Combination scores by mean
# This method is old and less accurate than combination of p-values by Stouffer

def filter_genes_scores_variance(
    genes_scores_distributions=None,
    threshold=None
):
    """
    Filter genes by variance across shuffles of their scores.

    arguments:
        genes_scores_distributions (dict<dict<dict>>): information about genes'
            scores and distributions
        threshold (float): value by which to filter variance

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and distributions

    """

    # Presumably these genes have inadequate real measurements.
    del genes_scores_distributions["ENSG00000183795"]
    del genes_scores_distributions["ENSG00000196866"]
    del genes_scores_distributions["ENSG00000249811"]
    del genes_scores_distributions["ENSG00000250913"]
    del genes_scores_distributions["ENSG00000233136"]
    del genes_scores_distributions["ENSG00000172352"]
    del genes_scores_distributions["ENSG00000233803"]

    if False:

        # Copy information.
        # Copying this information is expensive.
        #information = copy.deepcopy(genes_scores_distributions)
        # Extract identifiers of genes with scores.
        genes = list(genes_scores_distributions.keys())
        # Iterate on genes.
        for gene in genes:
            entry = genes_scores_distributions[gene]
            # Determine whether the gene's shuffle scores have adequate variance.
            # Consider using the isclose method of the math package.
            # Comparison to zero is problematic.
            variance_coefficient = (
                statistics.variance(entry["distributions"]["coefficient"])
            )
            variance_dip = (
                statistics.variance(entry["distributions"]["dip"])
            )
            variance_mixture = (
                statistics.variance(entry["distributions"]["mixture"])
            )
            if (
                (variance_coefficient < threshold) or
                (variance_dip < threshold) or
                (variance_mixture < threshold)
            ):
                print("Inadequate variance in shuffle scores!")
                print("Deleting gene:" + gene)
                del genes_scores_distributions[gene]
                pass
    # Report.
    # Extract identifiers of genes with scores.
    count = len(list(genes_scores_distributions.keys()))
    print("count of genes after filter by variance: " + str(count))
    # Return information.
    return genes_scores_distributions


def calculate_combination_scores(
    genes_scores_distributions=None
):
    """
    Calculates a score from combination of the dip statistic and the bimodality
    coefficient.

    Data structure.
    - genes_scores_distributions (dict)
    -- gene (dict)
    --- scores (dict)
    ---- combination (float)
    ---- coefficient (float)
    ---- dip (float)
    --- distributions (dict)
    ---- combination (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)

    arguments:
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:
        (dict<dict<dict>>): information about genes scores and distributions

    """

    # Report process.
    utility.print_terminal_partition(level=1)
    print("Combining scores for genes.")

    # Collect information about genes.
    #information = copy.deepcopy(genes_scores_distributions)
    # Extract identifiers of genes with scores.
    genes = list(genes_scores_distributions.keys())
    # Iterate on genes.
    for gene in genes:
        # Report progress.
        #print(gene)
        # Access information about gene's scores and distributions.
        gene_scores_distributions = genes_scores_distributions[gene]
        # Convert gene's scores and distributions to standard, z-score, space.
        standard = calculate_standard_scores_distributions(
            scores_distributions=gene_scores_distributions
        )
        # Calculate combination scores.
        combination = statistics.mean([
            standard["scores"]["coefficient"],
            standard["scores"]["dip"],
            standard["scores"]["mixture"],
        ])
        combinations = list()
        for index in range(len(standard["distributions"]["coefficient"])):
            coefficient = standard["distributions"]["coefficient"][index]
            dip = standard["distributions"]["dip"][index]
            mixture = standard["distributions"]["mixture"][index]
            value = statistics.mean([coefficient, dip, mixture])
            combinations.append(value)
            pass
        # Compile information.
        genes_scores_distributions[gene]["scores"]["combination"] = combination
        genes_scores_distributions[gene]["distributions"]["combination"] = (
            combinations
        )
    # Return information.
    return genes_scores_distributions


def calculate_standard_score(
    value=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard score, z-score, of a value.

    arguments:
        value (float): value to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (float): value in standard score space

    """

    return ((value - mean) / deviation)


def calculate_standard_scores(
    values=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard scores, z-scores, of values.

    arguments:
        value (list<float>): values to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (list<float>): value in standard score space

    """

    values_standard = list()
    for value in values:
        value_standard = calculate_standard_score(
            value=value,
            mean=mean,
            deviation=deviation
        )
        values_standard.append(value_standard)
    return values_standard


def calculate_standard_scores_distributions(
    scores_distributions=None
):
    """
    Calculates a gene's scores and distributions in standard, z-score, space.

    Data structure.
    - gene (dict)
    -- scores (dict)
    --- coefficient (float)
    --- dip (float)
    -- distributions (dict)
    --- coefficient (list)
    ---- value (float)
    --- dip (list)
    ---- value (float)

    arguments:
        scores_distributions (dict<dict>): information about a gene's scores
            and distributions

    raises:

    returns:
        (dict<dict>): information about a gene's scores and distributions

    """

    # Access information about gene's scores and distributions.
    coefficient = scores_distributions["scores"]["coefficient"]
    dip = scores_distributions["scores"]["dip"]
    mixture = scores_distributions["scores"]["mixture"]
    coefficients = (
        scores_distributions["distributions"]["coefficient"]
    )
    dips = scores_distributions["distributions"]["dip"]
    mixtures = scores_distributions["distributions"]["mixture"]
    # Convert gene's scores and distributions to standard, z-score, space.
    mean_coefficient = statistics.mean(coefficients)
    deviation_coefficient = statistics.stdev(coefficients)
    mean_dip = statistics.mean(dips)
    deviation_dip = statistics.stdev(dips)
    mean_mixture = statistics.mean(mixtures)
    deviation_mixture = statistics.stdev(mixtures)
    # Calculate standard scores.
    coefficient_standard = calculate_standard_score(
        value=coefficient,
        mean=mean_coefficient,
        deviation=deviation_coefficient
    )
    dip_standard = calculate_standard_score(
        value=dip,
        mean=mean_dip,
        deviation=deviation_dip
    )
    mixture_standard = calculate_standard_score(
        value=mixture,
        mean=mean_mixture,
        deviation=deviation_mixture
    )
    coefficients_standard = calculate_standard_scores(
        values=coefficients,
        mean=mean_coefficient,
        deviation=deviation_coefficient
    )
    dips_standard = calculate_standard_scores(
        values=dips,
        mean=mean_dip,
        deviation=deviation_dip
    )
    mixtures_standard = calculate_standard_scores(
        values=mixtures,
        mean=mean_mixture,
        deviation=deviation_mixture
    )
    # Compile information.
    information = dict()
    information["scores"] = dict()
    information["scores"]["coefficient"] = coefficient_standard
    information["scores"]["dip"] = dip_standard
    information["scores"]["mixture"] = mixture_standard
    information["distributions"] = dict()
    information["distributions"]["coefficient"] = coefficients_standard
    information["distributions"]["dip"] = dips_standard
    information["distributions"]["mixture"] = mixtures_standard
    # Return information.
    return information


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
    path_combination = os.path.join(dock, "combination")
    utility.confirm_path_directory(path_combination)
    path_genes = os.path.join(
        path_combination, "genes.txt"
    )
    # Write information to file.
    # List of genes needs to be easy to read in Bash.
    utility.write_file_text_list(
        information=information["genes"], path_file=path_genes
    )
    pass


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

    # The purpose of this procedure is to compare and integrate genes'
    # distributions and rankings from both the "imputation" and "availability"
    # methods for selection of tissues and persons.

    # Read source information from file.
    source = read_source(dock=dock)

    # Calculate genes' probabilities of bimodality.
    data_genes_probabilities_availability = rank.calculate_probabilities_genes(
        genes_scores_shuffles=(
            source["genes_scores_shuffles_availability"]
        ),
    )
    data_genes_probabilities_imputation = rank.calculate_probabilities_genes(
        genes_scores_shuffles=(
            source["genes_scores_shuffles_imputation"]
        ),
    )
    print(data_genes_probabilities_availability)
    print(data_genes_probabilities_imputation)

    # Rank genes by probabilities.
    data_ranks_availability = rank.rank_genes(
        data_genes_probabilities=data_genes_probabilities_availability,
        rank="combination",
        method="threshold",
        threshold=0.02,
        count=1000,
    )
    data_ranks_imputation = rank.rank_genes(
        data_genes_probabilities=data_genes_probabilities_imputation,
        rank="combination",
        method="threshold",
        threshold=0.02,
        count=1000,
    )
    print(data_ranks_availability)
    print(data_ranks_imputation)

    # Extract identifiers of genes.
    genes_availability = data_ranks_availability.index.to_list()
    genes_imputation = data_ranks_imputation.index.to_list()

    # Combine consensus genes.
    genes_combination = utility.filter_unique_union_elements(
        list_one=genes_availability,
        list_two=genes_imputation,
    )
    print("Count of consensus genes: " + str(len(genes_combination)))

    # Compile information.
    information = {
        "genes": genes_combination,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
