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
        genes=genes,
        path=path_pipe
    )
    if False:
        genes_persons = read_collect_genes_patients(
            path=path_pipe
        )
    # Compile and return information.
    return {
        "genes": genes,
        #"genes_persons": genes_persons,
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
    ---- mixture (float)
    --- distributions (dict)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)
    ---- mixture (list)
    ----- value (float)


    arguments:
        genes (list<str>): identifiers of genes
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
    # Check contents of directory.
    print("Check that directories exist for all genes.")
    directories = os.listdir(path)
    match = utility.compare_lists_by_mutual_inclusion(
        list_one=genes, list_two=directories
    )
    print("Genes and directories match: " + str(match))

    # Collect information about genes.
    genes_scores_distributions = dict()
    # Iterate on directories for genes.
    for gene in genes:
        # Report progress.
        print(gene)
        # Create entry for gene.
        genes_scores_distributions[gene] = dict()
        # Specify directories and files.
        path_gene = os.path.join(path, gene)
        path_scores = os.path.join(path_gene, "scores_real.pickle")
        path_distributions = os.path.join(
            path_gene, "scores_shuffle.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            scores = pickle.load(file_source)
        with open(path_distributions, "rb") as file_source:
            distributions = pickle.load(file_source)
        # Compile information.
        genes_scores_distributions[gene]["scores"] = scores
        genes_scores_distributions[gene]["distributions"] = distributions
    # Return information.
    return genes_scores_distributions


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
    # Check whether all genes have scores.
    # Iterate on genes.
    for gene in genes_scores:
        entry = genes_scores_distributions[gene]
        if not (
            (
                len(entry["distributions"]["coefficient"]) ==
                len(entry["distributions"]["dip"])
            ) or
            (
                len(entry["distributions"]["coefficient"]) ==
                len(entry["distributions"]["mixture"])
            )
        ):
            print(
                "***Error***: difference in counts of shuffle values for " +
                "scores!"
            )
            print(gene)
            pass
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
    path_genes_patients = os.path.join(
        path_collection, "genes_patients.pickle"
    )
    # Write information to file.
    with open(path_genes_patients, "wb") as file_product:
        pickle.dump(
            information["genes_patients"], file_product
        )
    if False:
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

    # TODO: ...
    # For each gene...
    # I need to know the real modality scores
    # I also need to know the 10,000 X shuffle modality scores...

    if False:
        # Check and summarize information about genes.
        check_genes_scores_distributions(
            genes=source["genes"],
            genes_scores_distributions=source["genes_scores_distributions"]
        )

    # Compile information.
    information = {
        #"genes_scores_distributions": source["genes_scores_distributions"],
        "genes_patients": source["genes_patients"],
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
