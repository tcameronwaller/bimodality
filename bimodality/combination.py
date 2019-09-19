"""
...

"""

###############################################################################
# Notes

# If the count of bootstrap permutations is too small, then it is likely not to
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
        path_collection, "genes_scores_permutations_availability.pickle"
    )
    path_scores_imputation = os.path.join(
        path_collection, "genes_scores_permutations_imputation.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    with open(path_scores_availability, "rb") as file_source:
        genes_scores_permutations_availability = pickle.load(file_source)
    with open(path_scores_imputation, "rb") as file_source:
        genes_scores_permutations_imputation = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "scores_permutations_availability": (
            genes_scores_permutations_availability
        ),
        "scores_permutations_imputation": (
            genes_scores_permutations_imputation
        ),
    }


def calculate_combination_scores(
    genes_scores_permutations=None
):
    """
    Calculates a score from combination of the dip statistic and the bimodality
    coefficient.

    Data structure.
    - genes_scores_permutations (dict)
    -- gene (dict)
    --- scores (dict)
    ---- combination (float)
    ---- coefficient (float)
    ---- dip (float)
    ---- mixture (float)
    --- permutations (dict)
    ---- combination (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)

    arguments:
        genes_scores_permutations (dict<dict<dict>>): information about genes

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Report process.
    utility.print_terminal_partition(level=1)
    print("Combining scores for genes.")

    # Collect information about genes.
    #information = copy.deepcopy(genes_scores_distributions)
    # Extract identifiers of genes with scores.
    genes = list(genes_scores_permutations.keys())
    # Iterate on genes.
    for gene in genes:
        # Report progress.
        #print(gene)
        # Access information about gene's scores and distributions.
        scores_permutations = genes_scores_permutations[gene]
        # Convert gene's scores and distributions to standard, z-score, space.
        standard = calculate_standard_scores_permutations(
            scores_permutations=scores_permutations
        )
        # Calculate combination scores.
        combination = statistics.mean([
            standard["scores"]["coefficient"],
            standard["scores"]["dip"],
            standard["scores"]["mixture"],
        ])
        combinations = list()
        for index in range(len(standard["permutations"]["coefficient"])):
            coefficient = standard["permutations"]["coefficient"][index]
            dip = standard["permutations"]["dip"][index]
            mixture = standard["permutations"]["mixture"][index]
            value = statistics.mean([coefficient, dip, mixture])
            combinations.append(value)
            pass
        # Compile information.
        genes_scores_permutations[gene]["scores"]["combination"] = combination
        genes_scores_permutations[gene]["permutations"]["combination"] = (
            combinations
        )
    # Return information.
    return genes_scores_permutations


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


def calculate_standard_scores_permutations(
    scores_permutations=None
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
        scores_permutations (dict<dict>): information about a gene's scores
            and permutations

    raises:

    returns:
        (dict<dict>): information about a gene's scores and permutations

    """

    # Access information about gene's scores and distributions.
    coefficient = scores_permutations["scores"]["coefficient"]
    dip = scores_permutations["scores"]["dip"]
    mixture = scores_permutations["scores"]["mixture"]
    coefficients = (
        scores_permutations["permutations"]["coefficient"]
    )
    dips = scores_permutations["permutations"]["dip"]
    mixtures = scores_permutations["permutations"]["mixture"]
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
    information["permutations"] = dict()
    information["permutations"]["coefficient"] = coefficients_standard
    information["permutations"]["dip"] = dips_standard
    information["permutations"]["mixture"] = mixtures_standard
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
    utility.create_directory(path_combination)
    path_scores_permutations_imputation = os.path.join(
        path_combination, "genes_scores_permutations_imputation.pickle"
    )
    path_scores_permutations_availability = os.path.join(
        path_combination, "genes_scores_permutations_availability.pickle"
    )
    # Write information to file.
    with open(path_scores_permutations_imputation, "wb") as file_product:
        pickle.dump(
            information["genes_scores_permutations_imputation"], file_product
        )
    with open(path_scores_permutations_availability, "wb") as file_product:
        pickle.dump(
            information["genes_scores_permutations_availability"], file_product
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

    # Calculate combination scores.
    # Transform coefficient, dip, and mixture scores to standard z-score space.
    # Calculate the combination score as the mean of the other three scores in
    # standard space.
    # Also calculate shuffle distributions of these scores in order to
    # calculate probabilities.
    genes_scores_permutations_availability = calculate_combination_scores(
        genes_scores_permutations=source["scores_permutations_availability"]
    )
    genes_scores_permutations_imputation = calculate_combination_scores(
        genes_scores_permutations=source["scores_permutations_imputation"]
    )

    print("done with calculation of combination scores...")

    # Compile information.
    information = {
        "genes_scores_permutations_imputation": (
            genes_scores_permutations_imputation
        ),
        "genes_scores_permutations_availability": (
            genes_scores_permutations_availability
        ),
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
