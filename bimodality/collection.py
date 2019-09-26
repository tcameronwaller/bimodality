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
    path_genes = os.path.join(path_split, "genes.txt")
    path_distribution = os.path.join(dock, "distribution")
    path_permutation = os.path.join(dock, "permutation")
    # Read information from file.
    genes = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes,
    )
    genes_distribution = utility.extract_subdirectory_names(
        path=path_distribution
    )
    genes_permutation = utility.extract_subdirectory_names(
        path=path_permutation
    )

    # Compile and return information.
    return {
        "genes": genes,
        "genes_distribution": genes_distribution,
        "genes_permutation": genes_permutation,
    }


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
        "Compare lists of genes from split, distribution, and permutation " +
        "procedures."
    )

    print("count of genes from split (master): " + str(len(genes)))
    print("count of genes from distribution: " + str(len(genes_distribution)))
    print("count of genes from permutation: " + str(len(genes_permutation)))

    print(
        "Check whether all permutation list genes are in distribution list."
    )
    # Returns True if all elements in list_two are in list_one.
    match = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes_permutation,
    )
    print("permutation genes match distribution genes: " + str(match))

    utility.print_terminal_partition(level=2)

    pass


def read_collect_genes_scores_permutations(
    genes=None,
    path_distribution=None,
    path_permutation=None,
):
    """
    Collects information about genes.

    Data structure.
    - genes_scores_shuffles (dict)
    -- gene (dict)
    --- scores (dict)
    ---- coefficient (float)
    ---- dip (float)
    ---- mixture (float)
    --- permutations (dict)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
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
    genes_scores = list(genes_scores_permutations.keys())
    print("count of genes: " + str(len(genes)))
    print("count of genes with scores: " + str(len(genes_scores)))
    gene_check = random.sample(genes, 1)
    print("example gene: " + str(gene_check[0]))
    entry = genes_scores_permutations[gene_check[0]]
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
    pass


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
    utility.create_directory(path_collection)
    path_scores_permutations = os.path.join(
        path_collection, "genes_scores_permutations.pickle"
    )
    # Write information to file.
    with open(path_scores_permutations, "wb") as file_product:
        pickle.dump(
            information["genes_scores_permutations"], file_product
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

    # Verify that distributions and permutations are available for all genes.
    check_genes(
        genes=source["genes"],
        genes_distribution=source["genes_distribution"],
        genes_permutation=source["genes_permutation"],
    )

    # Collect scores from distribution procedure.
    # Collect scores from permutation procedure.
    path_distribution = os.path.join(dock, "distribution")
    path_permutation = os.path.join(dock, "permutation")
    genes_scores_permutations = read_collect_genes_scores_permutations(
        genes=source["genes_permutation"],
        path_distribution=path_distribution,
        path_permutation=path_permutation,
    )

    # Check and summarize information about genes.
    check_genes_scores_permutations(
        genes=source["genes_permutation"],
        genes_scores_permutations=genes_scores_permutations
    )

    # Compile information.
    information = {
        "genes_scores_permutations": genes_scores_permutations,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
