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
import gc
import functools

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import assembly
import distribution
import permutation
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality

# Scrap before separation...

def collect_organize_genes_scores_permutations_old(
    genes=None,
    path_distribution=None,
    path_permutation=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        path_distribution (str): path to distribution directory
        path_permutation (str): path to permutation directory

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Collect scores from distribution procedure.
    # Collect scores from permutation procedure.
    genes_scores_permutations = read_collect_genes_scores_permutations(
        genes=genes,
        path_distribution=path_distribution,
        path_permutation=path_permutation,
    )

    # Check and summarize information about genes.
    utility.print_terminal_partition(level=3)
    print("original scores and permutations")
    check_genes_scores_permutations(
        genes=genes,
        genes_scores_permutations=genes_scores_permutations
    )

    # Standardize bimodality measures.
    # Calculate standard scores and permutations using the same values of
    # mean and standard deviation in order to allow comparison between scores
    # and permutations.
    genes_scores_permutations_standard = standardize_scores_permutations(
        genes_scores_permutations=genes_scores_permutations,
    )

    # Check and summarize information about genes.
    utility.print_terminal_partition(level=3)
    print("standardization complete")

    # Calculate combination scores.
    # Combination scores are scaled means of the primary bimodality measures.
    # Combination scores are useful to rank genes.
    genes_scores_permutations_combination = combine_scores_permutations(
        genes_scores_permutations=genes_scores_permutations_standard,
    )

    # Check and summarize information about genes.
    utility.print_terminal_partition(level=3)
    print("combination complete")

    # Return information.
    return genes_scores_permutations_combination

def read_collect_genes_scores_permutations_old(
    genes=None,
    path_distribution=None,
    path_permutation=None,
):
    """
    Collects information about genes.

    Data structure.
    - genes_scores_permutations (dict)
    -- gene (dict)
    --- scores (dict)
    ---- dip (float)
    ---- mixture (float)
    ---- coefficient (float)
    --- permutations (dict)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)

    arguments:
        genes (list<str>): identifiers of genes
        path_distribution (str): path to distribution directory
        path_permutation (str): path to permutation directory

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Report process.
    utility.print_terminal_partition(level=1)
    print(
        "Reading and compiling information about genes' scores and " +
        "permutations."
    )

    # Check contents of directory.
    utility.print_terminal_partition(level=3)
    print("Check that directories exist for all genes to collect.")
    genes_distribution = utility.extract_subdirectory_names(
        path=path_distribution
    )
    genes_permutation = utility.extract_subdirectory_names(
        path=path_permutation
    )
    match_distribution = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes
    )
    print(
        "Genes and distribution directories match: " +
        str(match_distribution)
    )
    match_permutation = utility.compare_lists_by_inclusion(
        list_one=genes_permutation,
        list_two=genes
    )
    print(
        "Genes and permutation directories match: " +
        str(match_permutation)
    )

    # Collect genes' scores and permutations.
    utility.print_terminal_partition(level=3)
    print("Collect genes' scores and permutations.")
    # Collect information about genes.
    genes_scores_permutations = dict()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_distribution_gene = os.path.join(path_distribution, gene)
        path_permutation_gene = os.path.join(path_permutation, gene)
        path_scores = os.path.join(
            path_distribution_gene, "scores.pickle"
        )
        path_permutations = os.path.join(
            path_permutation_gene, "permutations.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            scores = pickle.load(file_source)
        with open(path_permutations, "rb") as file_source:
            permutations = pickle.load(file_source)
        # Create entry for gene.
        genes_scores_permutations[gene] = dict()
        # Compile information.
        genes_scores_permutations[gene]["scores"] = scores
        genes_scores_permutations[gene]["permutations"] = permutations

    # Return information.
    return genes_scores_permutations

def standardize_scores_permutations_old(
    genes_scores_permutations=None
):
    """
    Calculates standard values of scores and permutations of the primary
    bimodality measures.

    There are two reasonable ways to standard bimodality measures.
    1. standardize relative to each gene
    - use mean and standard deviation from that gene
    2. standardize relative to each measure across all genes
    - use mean and standard deviation for each measure across all genes
    - this standardization is most appropriate in order to compare genes

    Data structure.
    - genes_scores_permutations (dict)
    -- gene (dict)
    --- scores (dict)
    ---- dip (float)
    ---- mixture (float)
    ---- coefficient (float)
    --- permutations (dict)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)

    arguments:
        genes_scores_permutations (dict<dict<dict>>): information about genes'
            scores and permutations

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Calculate values of mean and standard deviation from scores and
    # permutations for each primary bimodality measure across all genes.
    # This standardization allows subsequent comparison between genes.
    dip = calculate_mean_deviation_scores_permutations(
        measure="dip",
        genes_scores_permutations=genes_scores_permutations,
    )
    mixture = calculate_mean_deviation_scores_permutations(
        measure="mixture",
        genes_scores_permutations=genes_scores_permutations,
    )
    coefficient = calculate_mean_deviation_scores_permutations(
        measure="coefficient",
        genes_scores_permutations=genes_scores_permutations,
    )

    # Standardize scores and permutations.
    # Use values of mean and standard deviation for each primary bimodality
    # measure across all genes.
    # Collecting standard values in a new collection would require a lot of
    # memory.
    # Rather, change values to standard values in the original collection.
    # Iterate on genes.
    for gene in genes_scores_permutations:
        # Access gene's original scores and permutations.
        entry = genes_scores_permutations[gene]
        # Calculate standard scores.
        entry["scores"]["dip"] = utility.calculate_standard_score(
            value=entry["scores"]["dip"],
            mean=dip["mean"],
            deviation=dip["deviation"],
        )
        entry["scores"]["mixture"] = utility.calculate_standard_score(
            value=entry["scores"]["mixture"],
            mean=mixture["mean"],
            deviation=mixture["deviation"],
        )
        entry["scores"]["coefficient"] = utility.calculate_standard_score(
            value=entry["scores"]["coefficient"],
            mean=coefficient["mean"],
            deviation=coefficient["deviation"],
        )
        entry["permutations"]["dip"] = utility.calculate_standard_scores(
            values=entry["permutations"]["dip"],
            mean=dip["mean"],
            deviation=dip["deviation"],
        )
        entry["permutations"]["mixture"] = utility.calculate_standard_scores(
            values=entry["permutations"]["mixture"],
            mean=mixture["mean"],
            deviation=mixture["deviation"],
        )
        entry["permutations"]["coefficient"] = (
                utility.calculate_standard_scores(
                values=entry["permutations"]["coefficient"],
                mean=coefficient["mean"],
                deviation=coefficient["deviation"],
            )
        )

        # Compile information.
        genes_scores_permutations[gene] = entry
    # Return information.
    return genes_scores_permutations

def combine_scores_permutations_old(
    genes_scores_permutations=None
):
    """
    Calculates combination scores and permutations as the mean of standard
    primary bimodality measures.

    Data structure.
    - genes_scores_permutations (dict)
    -- gene (dict)
    --- scores (dict)
    ---- dip (float)
    ---- mixture (float)
    ---- coefficient (float)
    ---- combination (float)
    --- permutations (dict)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)
    ---- combination (list)
    ----- value (float)

    arguments:
        genes_scores_permutations (dict<dict<dict>>): information about genes'
            scores and permutations

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Extract identifiers of genes with scores.
    genes = list(genes_scores_permutations.keys())
    # Iterate on genes.
    for gene in genes:
        # Access information about gene's scores and distributions.
        entry = genes_scores_permutations[gene]
        # Calculate combination scores.
        combination = statistics.mean([
            entry["scores"]["coefficient"],
            entry["scores"]["dip"],
            entry["scores"]["mixture"],
        ])
        # Calculate combination permutations.
        combinations = list()
        for index in range(len(entry["permutations"]["coefficient"])):
            dip = entry["permutations"]["dip"][index]
            mixture = entry["permutations"]["mixture"][index]
            coefficient = entry["permutations"]["coefficient"][index]
            value = statistics.mean([dip, mixture, coefficient])
            combinations.append(value)
            pass
        # Compile information.
        entry["scores"]["combination"] = combination
        entry["permutations"]["combination"] = combinations
        genes_scores_permutations[gene].update(entry)
    # Return information.
    return genes_scores_permutations


###############################




##########
# Initialization


def initialize_paths(dock=None):
    """
    Initialize directories for procedure's product files.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Define paths to directories.
    path_distribution = os.path.join(dock, "distribution")
    path_distribution_genes = os.path.join(path_distribution, "genes")

    path_permutation = os.path.join(dock, "permutation_all")
    path_permutation_genes = os.path.join(path_permutation, "genes")

    # Collect information.
    paths = dict()
    paths["distribution"] = path_distribution
    paths["distribution_genes"] = path_distribution_genes
    paths["permutation"] = path_permutation
    paths["permutation_genes"] = path_permutation_genes

    # Return information.
    return paths


# Source


def read_source(
    paths=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        paths (dict<str>): collection of paths to directories for procedure's
            files
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
    path_selection_genes = os.path.join(
        path_selection, "genes_selection.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_selection_genes, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    genes_distribution = utility.extract_subdirectory_names(
        path=paths["distribution_genes"]
    )
    genes_permutation = utility.extract_subdirectory_names(
        path=paths["permutation_genes"]
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
        "genes_distribution": genes_distribution,
        "genes_permutation": genes_permutation,
    }


# Collection


def check_genes(
    genes=None,
    genes_distribution=None,
    genes_permutation=None,
):
    """
    Checks genes from split, distribution, and permutation batch procedures.

    arguments:
        genes (list<str>): identifiers of genes
        genes_distribution (list<str>): identifiers of genes
        genes_permutation (list<str>): identifiers of genes

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print(
        "Compare lists of genes from selection, distribution, and " +
        "permutation procedures."
    )

    print("count of genes from selection (master): " + str(len(genes)))
    print("count of genes from distribution: " + str(len(genes_distribution)))
    print("count of genes from permutation: " + str(len(genes_permutation)))

    utility.print_terminal_partition(level=3)

    print(
        "Check whether distributions and permutations are available for all " +
        "genes."
    )
    match_distribution = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes
    )
    print(
        "Genes and distribution directories match: " +
        str(match_distribution)
    )
    match_permutation = utility.compare_lists_by_inclusion(
        list_one=genes_permutation,
        list_two=genes
    )
    print(
        "Genes and permutation directories match: " +
        str(match_permutation)
    )

    utility.print_terminal_partition(level=3)

    # Returns True if all elements in list_two are in list_one.
    match = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes_permutation,
    )
    print("permutation genes match distribution genes: " + str(match))

    utility.print_terminal_partition(level=2)

    pass


def collect_organize_genes_scores_permutations(
    genes=None,
    paths=None,
):
    """
    Collects and organizes information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Collect scores from distribution procedure.
    # Collect scores from permutation procedure.
    collection = read_collect_genes_scores_permutations(
        genes=genes,
        paths=paths,
    )
    utility.print_terminal_partition(level=3)
    print("collection complete")

    # Check and summarize information about genes.
    if False:
        utility.print_terminal_partition(level=3)
        print("original scores and permutations")
        check_genes_scores_permutations(
            genes=genes,
            genes_scores_permutations=genes_scores_permutations
        )

    # Convert permutations.
    genes_permutations = convert_genes_permutations(
        genes_permutations=collection["genes_permutations"]
    )
    utility.print_terminal_partition(level=3)
    print("conversion complete")

    # Compile information.
    information = dict()
    information["genes_scores"] = collection["genes_scores"]
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


def read_collect_genes_scores_permutations(
    genes=None,
    paths=None,
):
    """
    Collects information about genes.

    Data structure.

    - genes_scores (dict)
    -- gene (dict)
    --- dip (float)
    --- mixture (float)
    --- coefficient (float)

    - genes_permutations (dict)
    -- gene (dict)
    --- dip (list)
    ---- value (float)
    --- mixture (list)
    ---- value (float)
    --- coefficient (list)
    ---- value (float)

    arguments:
        genes (list<str>): identifiers of genes
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (dict): information about genes' scores and permutations

    """

    # Collect genes' scores and permutations.
    utility.print_terminal_partition(level=2)
    print("Collect genes' scores and permutations.")
    # Collect information about genes.
    genes_scores = dict()
    genes_permutations = dict()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_distribution_gene = os.path.join(
            paths["distribution_genes"], gene
        )
        path_permutation_gene = os.path.join(
            paths["permutation_genes"], gene
        )
        path_scores = os.path.join(
            path_distribution_gene, "scores.pickle"
        )
        path_permutations = os.path.join(
            path_permutation_gene, "permutations.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            scores = pickle.load(file_source)
        with open(path_permutations, "rb") as file_source:
            permutations = pickle.load(file_source)
        # Create entry for gene.
        # Compile information.
        genes_scores[gene] = scores
        genes_permutations[gene] = permutations

    # Compile information.
    information = dict()
    information["genes_scores"] = genes_scores
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


def check_genes_scores_permutations(
    genes=None,
    genes_scores_permutations=None
):
    """
    Checks and summarizes information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        genes_scores_permutations (dict<dict<dict>>): information about genes'
            scores and permutations

    raises:

    returns:


    """

    # Extract identifiers of genes with scores.
    utility.print_terminal_partition(level=2)
    genes_scores = list(genes_scores_permutations.keys())
    print("count of genes: " + str(len(genes)))
    print("count of genes with scores: " + str(len(genes_scores)))
    gene_check = random.sample(genes, 1)[0]
    print("example gene: " + str(gene_check))
    entry = genes_scores_permutations[gene_check]
    print(
        "count of permutations for bimodality coefficient: " +
        str(len(entry["permutations"]["coefficient"]))
    )
    print(
        "count of shuffles for dip statistic: " +
        str(len(entry["permutations"]["dip"]))
    )
    print(
        "count of shuffles for mixture model: " +
        str(len(entry["permutations"]["mixture"]))
    )

    # Check whether all genes have scores.
    utility.print_terminal_partition(level=3)
    print("check for genes with null scores")
    # Iterate on genes.
    null_genes = list()
    for gene in genes_scores:
        entry = genes_scores_permutations[gene]
        # Check whether gene has valid scores.
        scores = entry["scores"]
        if (
            math.isnan(scores["coefficient"]) or
            math.isnan(scores["dip"]) or
            math.isnan(scores["mixture"])
        ):
            null_genes.append(gene)
        # Check whether gene has valid permutations.
        if (
            (
                len(entry["permutations"]["coefficient"]) !=
                len(entry["permutations"]["dip"])
            ) or
            (
                len(entry["permutations"]["coefficient"]) !=
                len(entry["permutations"]["mixture"])
            ) or
            (
                len(entry["permutations"]["dip"]) !=
                len(entry["permutations"]["mixture"])
            )
        ):
            print(
                "***Error***: difference in counts of permutation values " +
                "for scores!"
            )
            print(gene)
            pass
    print("count of genes with null scores: " + str(len(null_genes)))

    # Check standardization.
    utility.print_terminal_partition(level=3)
    print("check standardization of scores and permutations")
    print("... mean is not zero and standard deviation is not one because")
    print("standardization is across all genes")
    print("example gene: " + str(gene_check))
    entry = genes_scores_permutations[gene_check]
    for measure in ["dip", "mixture", "coefficient"]:
        print("----------")
        print(measure + ": " + str(entry["scores"][measure]))
        mean = statistics.mean(entry["permutations"][measure])
        deviation = statistics.stdev(entry["permutations"][measure])
        print("mean: " + str(mean))
        print("standard deviation: " + str(deviation))

    pass


# Standardization... now obsolete
# As modality measures have little correlation, the z-score combination of
# these measures is not reasonable.


def calculate_mean_deviation_scores_permutations(
    measure=None,
    genes_scores=None,
    genes_permutations=None,
):
    """
    Calculates the mean and standard deviation of scores and permutations for a
    single bimodality measure across all genes.

    arguments:
        measure (str): name of a single bimodality measure
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict<float>): values of mean and standard deviation

    """

    # Collect values of each primary bimodality measure.
    values = list()
    # Iterate on genes.
    for gene in genes_scores:
        # Access information about gene's scores and distributions.
        scores = genes_scores[gene]
        permutations = genes_permutations[gene]
        # Collect values of scores.
        values.append(scores[measure])
        # Collect values of permutations.
        values.extend(permutations[measure])

    # Calculate mean and standard deviation of the primary bimodality measure's
    # scores and permutations across all genes.
    #print(str(len(values)))
    mean = statistics.mean(values)
    deviation = statistics.stdev(values)
    # Compile information.
    information = dict()
    information["mean"] = mean
    information["deviation"] = deviation
    # Return information.
    return information


def standardize_scores_permutations(
    genes_scores=None,
    genes_permutations=None,
):
    """
    Calculates standard values of scores and permutations of the primary
    bimodality measures.

    There are two reasonable ways to standard bimodality measures.
    1. standardize relative to each gene
    - use mean and standard deviation from that gene
    2. standardize relative to each measure across all genes
    - use mean and standard deviation for each measure across all genes
    - this standardization is most appropriate in order to compare genes

    Data structure.
    - genes_scores (dict)
    -- gene (dict)
    --- dip (float)
    --- mixture (float)
    --- coefficient (float)

    - genes_permutations (dict)
    -- gene (dict)
    --- dip (list)
    ---- value (float)
    --- mixture (list)
    ---- value (float)
    --- coefficient (list)
    ---- value (float)

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict): information about genes' scores and permutations

    """

    # Calculate values of mean and standard deviation from scores and
    # permutations for each primary bimodality measure across all genes.
    # This standardization allows subsequent comparison between genes.
    dip = calculate_mean_deviation_scores_permutations(
        measure="dip",
        genes_scores=genes_scores,
        genes_permutations=genes_permutations
    )
    mixture = calculate_mean_deviation_scores_permutations(
        measure="mixture",
        genes_scores=genes_scores,
        genes_permutations=genes_permutations
    )
    coefficient = calculate_mean_deviation_scores_permutations(
        measure="coefficient",
        genes_scores=genes_scores,
        genes_permutations=genes_permutations
    )

    # Standardize scores and permutations.
    # Use values of mean and standard deviation for each primary bimodality
    # measure across all genes.
    # Collecting standard values in a new collection would require a lot of
    # memory.
    # Rather, change values to standard values in the original collection.
    # Iterate on genes.
    for gene in genes_scores:
        # Access gene's original scores and permutations.
        scores = genes_scores[gene]
        permutations = genes_permutations[gene]
        # Calculate standard scores.
        scores["dip"] = utility.calculate_standard_score(
            value=scores["dip"],
            mean=dip["mean"],
            deviation=dip["deviation"],
        )
        scores["mixture"] = utility.calculate_standard_score(
            value=scores["mixture"],
            mean=mixture["mean"],
            deviation=mixture["deviation"],
        )
        scores["coefficient"] = utility.calculate_standard_score(
            value=scores["coefficient"],
            mean=coefficient["mean"],
            deviation=coefficient["deviation"],
        )
        permutations["dip"] = utility.calculate_standard_scores(
            values=permutations["dip"],
            mean=dip["mean"],
            deviation=dip["deviation"],
        )
        permutations["mixture"] = utility.calculate_standard_scores(
            values=permutations["mixture"],
            mean=mixture["mean"],
            deviation=mixture["deviation"],
        )
        permutations["coefficient"] = (
                utility.calculate_standard_scores(
                values=permutations["coefficient"],
                mean=coefficient["mean"],
                deviation=coefficient["deviation"],
            )
        )

        # Compile information.
        genes_scores[gene] = scores
        genes_permutations[gene] = permutations
    # Compile information.
    information = dict()
    information["genes_scores"] = genes_scores
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


def combine_scores_permutations(
    genes_scores=None,
    genes_permutations=None
):
    """
    Calculates combination scores and permutations as the mean of standard
    primary bimodality measures.

    Data structure.
    - genes_scores_permutations (dict)
    -- gene (dict)
    --- scores (dict)
    ---- dip (float)
    ---- mixture (float)
    ---- coefficient (float)
    ---- combination (float)
    --- permutations (dict)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)
    ---- coefficient (list)
    ----- value (float)
    ---- combination (list)
    ----- value (float)

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict): information about genes' scores and permutations

    """

    # Iterate on genes.
    for gene in genes_scores:
        # Access information about gene's scores and distributions.
        scores = genes_scores[gene]
        permutations = genes_permutations[gene]
        # Calculate combination scores.
        combination = statistics.mean([
            scores["coefficient"],
            scores["dip"],
            scores["mixture"],
        ])
        # Calculate combination permutations.
        combinations = list()
        for index in range(len(permutations["coefficient"])):
            dip = permutations["dip"][index]
            mixture = permutations["mixture"][index]
            coefficient = permutations["coefficient"][index]
            value = statistics.mean([dip, mixture, coefficient])
            combinations.append(value)
            pass
        # Compile information.
        scores["combination"] = combination
        permutations["combination"] = combinations
        genes_scores[gene] = scores
        genes_permutations[gene] = permutations
    # Compile information.
    information = dict()
    information["genes_scores"] = genes_scores
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


def standardize_combine_scores_permutations(
    genes_scores=None,
    genes_permutations=None,
):
    """
    Standardizes and combines scores and their permutations.

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict<dict<dict>>): information about genes' scores and permutations

    """

    # Standardize bimodality measures.
    # Calculate standard scores and permutations using the same values of
    # mean and standard deviation in order to allow comparison between scores
    # and permutations.
    standardization = standardize_scores_permutations(
        genes_scores=genes_scores,
        genes_permutations=genes_permutations,
    )
    utility.print_terminal_partition(level=3)
    print("standardization complete")

    # Calculate combination scores.
    # Combination scores are scaled means of the primary bimodality measures.
    # Combination scores are useful to rank genes.
    combination = combine_scores_permutations(
        genes_scores=standardization["genes_scores"],
        genes_permutations=standardization["genes_permutations"],
    )
    utility.print_terminal_partition(level=3)
    print("combination complete")

    # Compile information.
    information = dict()
    information["genes_scores"] = combination["genes_scores"]
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


# Conversion


def convert_genes_permutations(
    genes_permutations=None
):
    """
    Converts genes' permutations from lists to NumPy arrays.

    arguments:
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict<array>): information about genes' permutations

    """

    # Iterate on genes.
    for gene in genes_permutations:
        # Access information about gene's scores and distributions.
        permutations = genes_permutations[gene]
        permutations["dips"] = numpy.asarray(
            permutations["dips"]
        )
        permutations["mixtures"] = numpy.asarray(
            permutations["mixtures"]
        )
        permutations["coefficients"] = numpy.asarray(
            permutations["coefficients"]
        )
        # Compile information.
        genes_permutations[gene] = permutations
    # Return information.
    return genes_permutations


# Probability


def calculate_probability_equal_greater(
    value=None,
    distribution=None,
):
    """
    Calculates from a distribution the probability of obtaining a value equal
    to or greater than a specific value.

    arguments:
        value (float): value
        distribution (list<float>): distribution of values

    raises:

    returns:
        (float): probability of obtaining from the distribution a value equal
            to or greater than the specific value

    """

    count_total = len(distribution)
    count_matches = 0
    for value_distribution in distribution:
        if value_distribution >= value:
            count_matches += 1
    probability = count_matches / count_total
    return probability


def calculate_organize_probabilities_genes(
    genes_scores=None,
    genes_permutations=None,
    data_gene_annotation=None,
):
    """
    Calculates probabilities (p-values) from genes' scores and random
    distributions.

    Data structure.
    - genes_p_values (dict)
    -- gene (dict)
    ---- dip (float)
    ---- mixture (float)
    ---- coefficient (float)
    ---- combination (float)

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' probabilities from permutations

    """

    # Collect information about genes.
    records = list()
    # Iterate on genes.
    for gene in genes_scores:
        name = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
        # Access information about gene's scores and distributions.
        scores = genes_scores[gene]
        permutations = genes_permutations[gene]
        # Scores
        dip = scores["dip"]
        mixture = scores["mixture"]
        coefficient = scores["coefficient"]
        # Permutations
        dips = permutations["dips"].tolist()
        mixtures = permutations["mixtures"].tolist()
        coefficients = permutations["coefficients"].tolist()
        # Calculate p-values.
        # These scores of bimodality indicate greater bimodality as values
        # increase.
        # The scores are unidirectional, so the hypothesis is unidirectional.
        # The p-value is the probability of obtaining by random chance a value
        # equal to or greater than the actual score.
        probability_dip = calculate_probability_equal_greater(
            value=dip,
            distribution=dips,
        )
        probability_mixture = calculate_probability_equal_greater(
            value=mixture,
            distribution=mixtures,
        )
        probability_coefficient = calculate_probability_equal_greater(
            value=coefficient,
            distribution=coefficients,
        )
        # Compile information.
        record = dict()
        record["gene"] = gene
        record["name"] = name
        record["probability_dip"] = probability_dip
        record["probability_mixture"] = probability_mixture
        record["probability_coefficient"] = probability_coefficient
        records.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records,
    )
    data.set_index(
        "gene",
        drop=True,
        inplace=True,
    )
    # Return information.
    return data


def calculate_discoveries_genes(
    threshold=None,
    data_genes_probabilities=None,
):
    """
    Calculates false discovery rates from genes' probabilities.

    arguments:
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        data_genes_probabilities (object): Pandas data frame of genes'
            probabilities from permutations

    raises:

    returns:
        (object): Pandas data frame of genes' false discovery rates from
            permutations

    """

    # Copy data.
    data_copy = data_genes_probabilities.copy(deep=True)
    # Calculate false discovery rates.
    data_coefficient = (
        utility.calculate_false_discovery_rate(
            threshold=threshold,
            probability="probability_coefficient",
            discovery="discovery_coefficient",
            significance="significance_coefficient",
            data_probabilities=data_copy,
        )
    )
    data_mixture = (
        utility.calculate_false_discovery_rate(
            threshold=threshold,
            probability="probability_mixture",
            discovery="discovery_mixture",
            significance="significance_mixture",
            data_probabilities=data_coefficient,
        )
    )
    data_dip = (
        utility.calculate_false_discovery_rate(
            threshold=threshold,
            probability="probability_dip",
            discovery="discovery_dip",
            significance="significance_dip",
            data_probabilities=data_mixture,
        )
    )
    # Return information.
    return data_dip


# Selection


def select_genes_measure_discoveries(
    significance=None,
    data_genes_discoveries=None,
):
    """
    Selects genes with significant false discovery rates for a single measure
    of modality.

    arguments:
        significance (str): name of column of significance
        data_genes_discoveries (object): Pandas data frame of genes'
            false discovery rates for measures of modality

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Copy genes' heritabilities.
    data_copy = data_genes_discoveries.copy(deep=True)

    # Set threshold.
    data_significance = data_copy.loc[data_copy[significance]]

    # Extract identifiers of genes.
    genes = data_significance.index.to_list()

    # Return information.
    return genes


def select_genes_measures_discoveries(
    data_genes_discoveries=None,
):
    """
    Selects genes by their false discovery rates for measures of modality.

    arguments:
        data_genes_discoveries (object): Pandas data frame of genes'
            false discovery rates for measures of modality

    raises:

    returns:
        (dict<list<str>>): sets of genes with significant false discovery rates for each
            measure of modality
    """

    # Select genes for each measure of modality.
    genes_coefficient = select_genes_measure_discoveries(
        significance="significance_coefficient",
        data_genes_discoveries=data_genes_discoveries,
    )
    genes_mixture = select_genes_measure_discoveries(
        significance="significance_mixture",
        data_genes_discoveries=data_genes_discoveries,
    )
    genes_dip = select_genes_measure_discoveries(
        significance="significance_dip",
        data_genes_discoveries=data_genes_discoveries,
    )

    # Compile information.
    information = dict()
    information["coefficient"] = genes_coefficient
    information["mixture"] = genes_mixture
    information["dip"] = genes_dip
    # Return information.
    return information


# Sets


def organize_genes_sets(
    genes_selection=None,
    genes_unimodal=None,
    genes_multimodal=None,
    genes_probability=None,
):
    """
    Organize sets of genes.

    arguments:
        genes_selection (list<str>): identifiers of genes from selection
        genes_unimodal (list<str>): identifiers of genes with unimodal
            distributions in their pan-tissue signals across persons
        genes_multimodal (list<str>): identifiers of genes with multimodal
            distributions in their pan-tissue signals across persons
        genes_probability (list<str>): identifiers of genes with significant
            distributions in their pan-tissue signals across persons

    raises:

    returns:
        (dict): sets of genes
    """

    # Organize sets of genes.
    # Selection
    sets_genes_selection = dict()
    sets_genes_selection["selection"] = genes_selection
    sets_genes_selection["probability"] = genes_probability
    # Unimodal
    sets_genes_unimodal = dict()
    sets_genes_unimodal["unimodal"] = genes_unimodal
    sets_genes_unimodal["probability"] = genes_probability
    # Multimodal
    sets_genes_multimodal = dict()
    sets_genes_multimodal["multimodal"] = genes_multimodal
    sets_genes_multimodal["probability"] = genes_probability

    # Compile information.
    information = dict()
    information["selection"] = sets_genes_selection
    information["unimodal"] = sets_genes_unimodal
    information["multimodal"] = sets_genes_multimodal
    # Return information.
    return information


# Permutation


def permute_collect_gene_signals_all(
    gene=None,
    data_gene_persons_tissues_signals=None,
    shuffles=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues
        shuffles (list<list<list<int>>>): matrices of shuffle indices

    raises:

    returns:
        (array<float>): gene's signals across persons and permutations

    """

    def permute_collect_signals(
        collection=None,
        shuffle_indices=None,
        data_gene_signals=None,
    ):
        # Shuffle gene's signals.
        data_shuffle = permutation.permute_gene_signals(
            data_gene_signals=data_gene_signals,
            shuffle_indices=shuffle_indices
        )
        # Distribution
        bins = distribution.prepare_describe_distribution(
            data_gene_persons_tissues_signals=data_shuffle,
        )
        # Extract signals.
        #signals = bins["bin_aggregation"]["values"]
        data = bins["bin_aggregation"]["data_gene_persons_signals"].copy(
            deep=True
        )
        data.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        signals = data["value"].to_numpy()
        # Collect signals.
        collection.append(signals)
        return collection

    # Initialize collections of measures.
    collection = list()

    # Iterate on permutations.
    signals_permutations = functools.reduce(
        lambda collection, shuffle_indices: permute_collect_signals(
            collection=collection,
            shuffle_indices=shuffle_indices,
            data_gene_signals=data_gene_persons_tissues_signals,
        ),
        shuffles, # Iterable
        collection, # Initializer
    )
    # Collapse signals to a single dimension.
    signals = numpy.concatenate(
        signals_permutations,
        axis=None,
    )
    # Return information.
    return signals


def permute_collect_gene_signals_single(
    gene=None,
    data_gene_persons_tissues_signals=None,
    shuffle=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues
        shuffle (list<list<int>>): matrix of shuffle indices

    raises:

    returns:
        (array<float>): gene's signals across persons

    """

    # Shuffle gene's signals.
    data_shuffle = permutation.permute_gene_signals(
        data_gene_signals=data_gene_persons_tissues_signals,
        shuffle_indices=shuffle
    )
    # Distribution
    bins = distribution.prepare_describe_distribution(
        data_gene_persons_tissues_signals=data_shuffle,
    )
    # Extract signals.
    #signals = bins["bin_aggregation"]["values"]
    data = bins["bin_aggregation"]["data_gene_persons_signals"].copy(
        deep=True
    )
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    signals = data["value"].to_numpy()

    # Return information.
    return signals


def read_permute_collect_gene_signals(
    gene=None,
    permutations=None,
    dock=None,
):
    """
    Organize sets of genes.

    arguments:
        gene (str): identifier of a gene
        permutations (str): whether to execute "single" or "all" permutations
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (list<float>): gene's pan-tissue signals across persons and
            permutations
    """

    # Read shuffles.
    source_shuffles = permutation.read_source_shuffles(dock=dock)
    # Read gene's signals.
    source_gene = permutation.read_source(
        gene=gene,
        dock=dock
    )
    # Determine whether to execute a single permutation or all permutations.
    if permutations == "single":
        # Select shuffle indices at random.
        shuffle = random.choice(source_shuffles)
        # Permute gene's signals.
        signals = permute_collect_gene_signals_single(
            gene=gene,
            data_gene_persons_tissues_signals=(
                source_gene["data_gene_persons_tissues_signals"]
            ),
            shuffle=shuffle,
        )
        pass
    elif permutations == "all":
        # Permute gene's signals.
        signals = permute_collect_gene_signals_all(
            gene=gene,
            data_gene_persons_tissues_signals=(
                source_gene["data_gene_persons_tissues_signals"]
            ),
            shuffles=source_shuffles,
        )
        pass

    # Return information.
    return signals


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
    path_probability = os.path.join(dock, "probability")
    utility.create_directory(path_probability)

    path_data_genes_discoveries = os.path.join(
        path_probability, "data_genes_discoveries.pickle"
    )
    path_data_genes_discoveries_text = os.path.join(
        path_probability, "data_genes_discoveries.tsv"
    )
    path_sets_genes_measures = os.path.join(
        path_probability, "sets_genes_measures.pickle"
    )
    path_genes_probability = os.path.join(
        path_probability, "genes_probability.pickle"
    )
    path_sets_genes_selection = os.path.join(
        path_probability, "sets_genes_selection.pickle"
    )
    path_sets_genes_unimodal = os.path.join(
        path_probability, "sets_genes_unimodal.pickle"
    )
    path_sets_genes_multimodal = os.path.join(
        path_probability, "sets_genes_multimodal.pickle"
    )
    path_gene_signals_permutation = os.path.join(
        path_probability, "gene_signals_permutation.pickle"
    )

    # Write information to file.
    information["data_genes_discoveries"].to_pickle(
        path=path_data_genes_discoveries
    )
    information["data_genes_discoveries"].to_csv(
        path_or_buf=path_data_genes_discoveries_text,
        sep="\t",
        header=True,
        index=True,
    )
    with open(path_sets_genes_measures, "wb") as file_product:
        pickle.dump(
            information["sets_genes_measures"], file_product
        )
    with open(path_genes_probability, "wb") as file_product:
        pickle.dump(
            information["genes_probability"], file_product
        )
    with open(path_sets_genes_selection, "wb") as file_product:
        pickle.dump(
            information["sets_genes_selection"], file_product
        )
    with open(path_sets_genes_unimodal, "wb") as file_product:
        pickle.dump(
            information["sets_genes_unimodal"], file_product
        )
    with open(path_sets_genes_multimodal, "wb") as file_product:
        pickle.dump(
            information["sets_genes_multimodal"], file_product
        )
    with open(path_gene_signals_permutation, "wb") as file_product:
        pickle.dump(
            information["gene_signals_permutation"], file_product
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

    # Remove previous files to avoid version or batch confusion.
    path_probability = os.path.join(dock, "probability")
    utility.remove_directory(path=path_probability)

    # Initialize paths.
    paths = initialize_paths(dock=dock)

    # Read source information from file.
    source = read_source(
        paths=paths,
        dock=dock
    )

    # Define genes for iteration.
    genes_iteration = source["genes_selection"]
    # Unimodal and multimodal genes only (290).
    #genes_iteration = copy.deepcopy(source["genes_unimodal"])
    #genes_iteration.extend(source["genes_multimodal"])
    print("count of genes: " + str(len(genes_iteration)))

    # Verify that distributions and permutations are available for all genes.
    check_genes(
        genes=genes_iteration,
        genes_distribution=source["genes_distribution"],
        genes_permutation=source["genes_permutation"],
    )

    # Collect and organize genes' scores and permutations.
    # Contain these processes within a separate function to conserve memory.
    # This procedure standardizes values of each bimodality measure.
    genes_scores_permutations = collect_organize_genes_scores_permutations(
        genes=genes_iteration,
        paths=paths,
    )

    utility.print_terminal_partition(level=2)
    print("collection, standardization, and combination done!!!")
    utility.print_terminal_partition(level=2)

    # Calculate genes' probabilities of bimodality.
    utility.print_terminal_partition(level=2)
    data_genes_probabilities = calculate_organize_probabilities_genes(
        genes_scores=genes_scores_permutations["genes_scores"],
        genes_permutations=genes_scores_permutations["genes_permutations"],
        data_gene_annotation=source["data_gene_annotation"],
    )
    print(data_genes_probabilities)

    # Calculate genes' false discovery rates.
    utility.print_terminal_partition(level=2)
    data_genes_discoveries = calculate_discoveries_genes(
        threshold=0.05,
        data_genes_probabilities=data_genes_probabilities,
    )
    print(data_genes_discoveries)

    # Select genes by false discovery rates for each measure of modality.
    sets_genes_measures = select_genes_measures_discoveries(
        data_genes_discoveries=data_genes_discoveries,
    )
    utility.print_terminal_partition(level=2)
    print("coefficient: " + str(len(sets_genes_measures["coefficient"])))
    print("mixture: " + str(len(sets_genes_measures["mixture"])))
    print("dip: " + str(len(sets_genes_measures["dip"])))

    # Select genes that are significant by any single measure of modality.
    genes_probability = utility.select_elements_by_sets(
        names=["dip", "mixture", "coefficient"],
        sets=sets_genes_measures,
        count=1,
    )
    utility.print_terminal_partition(level=2)
    print(
        "selection of genes by multiple measurements: " +
        str(len(genes_probability))
    )
    utility.print_terminal_partition(level=2)

    # Organize sets.
    bin = organize_genes_sets(
        genes_selection=source["genes_selection"],
        genes_unimodal=source["genes_unimodal"],
        genes_multimodal=source["genes_multimodal"],
        genes_probability=genes_probability,
    )

    # Collect an exemplary gene's signals from permutation.
    utility.print_terminal_partition(level=2)
    gene_signals_permutation = read_permute_collect_gene_signals(
        gene="ENSG00000183793", # NPIPA5
        permutations="single", # "single" or "all"
        dock=dock,
    )
    #print(gene_signals_permutation)
    print("Count of permuted signals: " + str(len(gene_signals_permutation)))

    # Compile information.
    information = {
        "data_genes_discoveries": data_genes_discoveries,
        "sets_genes_measures": sets_genes_measures,
        "sets_genes_selection": bin["selection"],
        "sets_genes_unimodal": bin["unimodal"],
        "sets_genes_multimodal": bin["multimodal"],
        "genes_probability": genes_probability,
        "gene_signals_permutation": gene_signals_permutation,
    }
    # Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
