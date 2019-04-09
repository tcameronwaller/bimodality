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

    # Collect information about genes.
    genes_scores_distributions = dict()
    # Iterate on directories for genes.
    directories = os.listdir(path)
    for directory in directories:
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
    print("count of genes: " + str(len(genes_scores)))
    gene_entries = genes_scores_distributions["ENSG00000005483"]
    print(
        "count of shuffles: " +
        str(len(
            gene_entries["distributions"]["coefficient"]
        ))
    )
    # Check whether all genes have scores.

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

    # Compile information.
    information = {
        "genes_scores_distributions": source["genes_scores_distributions"],
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
