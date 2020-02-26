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


# Source


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
    path_selection = os.path.join(dock, "selection", "tight")
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

    path_distribution = os.path.join(dock, "distribution")
    path_permutation = os.path.join(dock, "permutation")

    # Read information from file.
    with open(path_selection_genes, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)
    genes_distribution = utility.extract_subdirectory_names(
        path=path_distribution
    )
    genes_permutation = utility.extract_subdirectory_names(
        path=path_permutation
    )

    # Compile and return information.
    return {
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


# TODO: remove the calculation of combination scores...
# TODO: as measures do not have strong correlation, combination scores are meaningless


def collect_organize_genes_scores_permutations(
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
    collection = read_collect_genes_scores_permutations(
        genes=genes,
        path_distribution=path_distribution,
        path_permutation=path_permutation,
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
    path_distribution=None,
    path_permutation=None,
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
        path_distribution (str): path to distribution directory
        path_permutation (str): path to permutation directory

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
        # Compile information.
        genes_scores[gene] = scores
        genes_permutations[gene] = permutations

    # Compile information.
    information = dict()
    information["genes_scores"] = genes_scores
    information["genes_permutations"] = genes_permutations
    # Return information.
    return information


def read_collect_genes_patients(path=None):
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
    genes_patients = list()
    # Iterate on directories for genes.
    directories = os.listdir(path)
    for directory in directories:
        # Report progress.
        print(directory)
        # Create entry for gene.
        #genes_patients[directory] = list()
        # Specify directories and files.
        path_directory = os.path.join(path, directory)
        path_report = os.path.join(path_directory, "report_gene.pickle")
        # Read information from file.
        with open(path_report, "rb") as file_source:
            report = pickle.load(file_source)
        # Compile information.
        patients = report["data_patients_tissues"].size
        genes_patients.append(patients)
    # Return information.
    return genes_patients


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

    raises:

    returns:
        (object): Pandas data frame of genes' probabilities from permutations

    """

    # Collect information about genes.
    records = list()
    # Iterate on genes.
    for gene in genes_scores:
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





def select_genes_by_threshold(
    genes_probabilities=None,
    measure=None,
    threshold=None,
):
    """
    Selects genes by a threshold on probabilities for bimodality measures.

    arguments:
        genes_probabilities (dict<dict<float>>): information about genes'
            probabilities
        measure (str): name of a measure of bimodality
        threshold (float): value of a threshold for a measure of bimodality

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities

    """

    genes = list()
    for gene in genes_probabilities:
        if genes_probabilities[gene][measure] < threshold:
            genes.append(gene)
    return genes


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
    path_scores = os.path.join(
        path_probability, "genes_scores.pickle"
    )
    path_permutations = os.path.join(
        path_probability, "genes_permutations.pickle"
    )
    path_probabilities = os.path.join(
        path_probability, "genes_probabilities.pickle"
    )
    path_sets_genes_measures = os.path.join(
        path_probability, "sets_genes_measures.pickle"
    )
    path_genes_probability = os.path.join(
        path_probability, "genes_probability.pickle"
    )


    # Write information to file.
    with open(path_scores, "wb") as file_product:
        pickle.dump(
            information["genes_scores"], file_product
        )
    with open(path_permutations, "wb") as file_product:
        pickle.dump(
            information["genes_permutations"], file_product
        )
    with open(path_probabilities, "wb") as file_product:
        pickle.dump(
            information["genes_probabilities"], file_product
        )
    with open(path_sets_genes_measures, "wb") as file_product:
        pickle.dump(
            information["sets_genes_measures"], file_product
        )
    with open(path_genes_probability, "wb") as file_product:
        pickle.dump(
            information["genes_probability"], file_product
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

    #####*********DO THIS####
    # TODO: For a limited count of genes (of interest...) re-run the permutation
    # TODO: procedure and collect the pan-tissue signal distribution across all permutations
    # TODO: NPIP15

    # Remove previous files to avoid version or batch confusion.
    path_probability = os.path.join(dock, "probability")
    utility.remove_directory(path=path_probability)

    # Read source information from file.
    source = read_source(dock=dock)

    # Define genes for iteration.
    #genes_iteration = source["genes_selection"]
    #genes_iteration = source["genes_unimodal"]
    genes_iteration = copy.deepcopy(source["genes_unimodal"])
    genes_iteration.extend(source["genes_multimodal"])
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
    path_distribution = os.path.join(dock, "distribution", "genes")
    path_permutation = os.path.join(dock, "permutation")
    genes_scores_permutations = collect_organize_genes_scores_permutations(
        genes=genes_iteration,
        path_distribution=path_distribution,
        path_permutation=path_permutation,
    )

    utility.print_terminal_partition(level=2)
    print("collection, standardization, and combination done!!!")
    utility.print_terminal_partition(level=2)



    # TODO: I think it makes more sense to organize genes, their probabilities, and their FDRs in a dataframe
    # TODO: calculate FDRs
    # TODO: organize data frame

    # Calculate genes' probabilities of bimodality.
    utility.print_terminal_partition(level=2)
    data_genes_probabilities = calculate_organize_probabilities_genes(
        genes_scores=genes_scores_permutations["genes_scores"],
        genes_permutations=genes_scores_permutations["genes_permutations"],
    )
    print(data_genes_probabilities)

    # Calculate genes' false discovery rates.
    utility.print_terminal_partition(level=2)
    data_genes_discoveries = calculate_discoveries_genes(
        threshold=0.05,
        data_genes_probabilities=data_genes_probabilities,
    )
    print(data_genes_discoveries)

    # TODO: I still want to select genes by the set method (probably by any measure for all genes...)
    # TODO: maybe implement similar to "heritability" procedure

    if False:
        # Select genes by probabilities.
        utility.print_terminal_partition(level=2)
        print("Select genes by threshold on measures' probabilities...")
        genes_dip = select_genes_by_threshold(
            genes_probabilities=genes_probabilities,
            measure="dip",
            threshold=0.05, # 0.001 (61), 0.005 (125), 0.05 (675)
        )
        genes_mixture = select_genes_by_threshold(
            genes_probabilities=genes_probabilities,
            measure="mixture",
            threshold=0.0005, # 0.0005 (2626), 0.001 (2945), 0.005 (3895),
        )
        genes_coefficient = select_genes_by_threshold(
            genes_probabilities=genes_probabilities,
            measure="coefficient",
            threshold=0.001, # 0.0005 (1160), 0.001 (1331), 0.005 (1827)
        )
        utility.print_terminal_partition(level=2)
        print("dip: " + str(len(genes_dip)))
        print("mixture: " + str(len(genes_mixture)))
        print("coefficient: " + str(len(genes_coefficient)))
        sets_genes_measures = dict()
        sets_genes_measures["dip"] = genes_dip
        sets_genes_measures["mixture"] = genes_mixture
        sets_genes_measures["coefficient"] = genes_coefficient

        # Select genes that pass filters by multiple measures of bimodality.
        genes_probability = utility.select_elements_by_sets(
            names=["dip", "mixture", "coefficient"],
            sets=sets_genes_measures,
            count=2,
        )
        utility.print_terminal_partition(level=2)
        print(
            "selection of genes by multiple measurements: " +
            str(len(genes_probability))
        )
        utility.print_terminal_partition(level=2)

        # Compile information.
        information = {
            "genes_scores": genes_scores_permutations["genes_scores"],
            "genes_permutations": genes_scores_permutations["genes_permutations"],
            "genes_probabilities": genes_probabilities,
            "sets_genes_measures": sets_genes_measures,
            "genes_probability": genes_probability,
        }
        # Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
