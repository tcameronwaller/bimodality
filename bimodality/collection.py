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

import metric
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

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_pipe = os.path.join(dock, "pipe")
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    genes_scores_distributions = read_collect_genes_scores_distributions(
        path=path_pipe
    )
    # Compile and return information.
    return {
        "genes": genes,
        "genes_scores_distributions": genes_scores_distributions
    }


def read_collect_genes_scores_distributions(path=None):
    """
    Collects information about genes.

    Data structure.
    - genes_scores_distributions (dict)
    -- gene (dict)
    --- scores (dict)
    ---- coefficient (float)
    ---- dip (float)
    --- distributions (dict)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)

    arguments:
        path (str): path to directory

    raises:

    returns:
        (dict<dict<dict>>): information about genes

    """

    # Report process.
    utility.print_terminal_partition(level=1)
    print(
        "Reading and compiling information about genes' scores and " +
        "distributions."
    )
    # Collect information about genes.
    genes_scores_distributions = dict()
    # Iterate on directories for genes.
    directories = os.listdir(path)
    for directory in directories:
        # Report progress.
        print(directory)
        # Create entry for gene.
        genes_scores_distributions[directory] = dict()
        # Specify directories and files.
        path_directory = os.path.join(path, directory)
        path_scores = os.path.join(path_directory, "scores.pickle")
        path_distributions = os.path.join(
            path_directory, "distributions.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            scores = pickle.load(file_source)
        with open(path_distributions, "rb") as file_source:
            distributions = pickle.load(file_source)
        # Compile information.
        genes_scores_distributions[directory]["scores"] = scores
        genes_scores_distributions[directory]["distributions"] = distributions
    # Return information.
    return genes_scores_distributions


def check_genes_scores_distributions(
    genes=None,
    genes_scores_distributions=None
):
    """
    Checks and summarizes information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:


    """

    #print(genes_scores_distributions)
    #print(genes_scores_distributions["ENSG00000005483"])

    # Extract identifiers of genes with scores.
    genes_scores = list(genes_scores_distributions.keys())
    print("count of genes: " + str(len(genes)))
    print("count of genes with scores: " + str(len(genes_scores)))
    gene_entries = genes_scores_distributions["ENSG00000005483"]
    print(
        "count of shuffles for bimodality coefficient: " +
        str(len(
            gene_entries["distributions"]["coefficient"]
        ))
    )
    print(
        "count of shuffles for dip statistic: " +
        str(len(
            gene_entries["distributions"]["dip"]
        ))
    )
    print(
        "count of shuffles for mixture model: " +
        str(len(
            gene_entries["distributions"]["mixture"]
        ))
    )
    if False:
        # Check whether all genes have scores.
        # Iterate on genes.
        for gene in genes_scores[5000:]:
            entry = genes_scores_distributions[gene]
            if not (
                len(entry["distributions"]["coefficient"]) ==
                len(entry["distributions"]["dip"])
            ):
                print(
                    "***Error***: difference in counts of shuffle values for " +
                    "scores!"
                )
                print(gene)
                pass
    pass


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
    if False:
        del genes_scores_distributions["ENSG00000183795"]
        del genes_scores_distributions["ENSG00000196866"]
        del genes_scores_distributions["ENSG00000249811"]
        del genes_scores_distributions["ENSG00000250913"]
        del genes_scores_distributions["ENSG00000233136"]
        del genes_scores_distributions["ENSG00000172352"]
        del genes_scores_distributions["ENSG00000233803"]

    # Copy information.
    information = copy.deepcopy(genes_scores_distributions)
    # Extract identifiers of genes with scores.
    genes = list(information.keys())
    # Iterate on genes.
    for gene in genes:
        entry = information[gene]
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
            del information[gene]
            pass
    # Report.
    # Extract identifiers of genes with scores.
    count = len(list(information.keys()))
    print("count of genes after filter by variance: " + str(count))
    # Return information.
    return information


# Score


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
    information = genes_scores_distributions
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
        information[gene]["scores"]["combination"] = combination
        information[gene]["distributions"]["combination"] = combinations
    # Return information.
    return information


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
    path_collection = os.path.join(dock, "collection")
    utility.confirm_path_directory(path_collection)
    path_distributions = os.path.join(
        path_collection, "genes_scores_distributions.pickle"
    )
    # Write information to file.
    with open(path_distributions, "wb") as file_product:
        pickle.dump(
            information["genes_scores_distributions"], file_product
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

    # Read source information from file.
    source = read_source(dock=dock)

    # Check and summarize information about genes.
    check_genes_scores_distributions(
        genes=source["genes"],
        genes_scores_distributions=source["genes_scores_distributions"]
    )

    # Filter genes by variance across shuffles of their scores.
    # Presumably these genes have inadequate real measurements.
    # 18804 genes before filter (TCW, 2019-04-22)
    # 18797 genes after filter (TCW, 2019-04-22).
    genes_scores_distributions = filter_genes_scores_variance(
        genes_scores_distributions=source["genes_scores_distributions"],
        threshold=1E-25
    )

    # Combine modality scores.
    # This procedure requires approximately 1.5 hours to run in a single
    # process on a 3.9 Gigahertz processor.
    genes_scores_distributions_complete = calculate_combination_scores(
        genes_scores_distributions=genes_scores_distributions
    )
    #print(genes_scores_distributions_complete["ENSG00000183878"])

    # Compile information.
    information = {
        "genes_scores_distributions": genes_scores_distributions_complete,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
