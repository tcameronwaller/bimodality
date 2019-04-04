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
    path_assembly = os.path.join(dock, "assembly")
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_organization = os.path.join(dock, "organization")
    path_mean = os.path.join(
        path_organization, "data_gene_signal_mean.pickle"
    )
    path_median = os.path.join(
        path_organization, "data_gene_signal_median.pickle"
    )
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    path_log = os.path.join(
        path_organization, "data_gene_signal_log.pickle"
    )
    path_standard = os.path.join(
        path_organization, "data_gene_signal_standard.pickle"
    )
    path_aggregation = os.path.join(dock, "aggregation")
    path_signal = os.path.join(
        path_aggregation, "data_gene_signal_aggregation.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_gene_signal_mean = pandas.read_pickle(path_mean)
    data_gene_signal_median = pandas.read_pickle(path_median)
    data_gene_signal_imputation = pandas.read_pickle(path_imputation)
    data_gene_signal_log = pandas.read_pickle(path_log)
    data_gene_signal_standard = pandas.read_pickle(path_standard)
    data_gene_signal_aggregation = pandas.read_pickle(path_signal)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal_mean": data_gene_signal_mean,
        "data_gene_signal_median": data_gene_signal_median,
        "data_gene_signal_imputation": data_gene_signal_imputation,
        "data_gene_signal_log": data_gene_signal_log,
        "data_gene_signal_standard": data_gene_signal_standard,
        "data_gene_signal_aggregation": data_gene_signal_aggregation,
    }


def collect_genes_names(
    data_gene_annotation=None,
    data_gene_score=None
):
    """
    Collects names of genes.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations.
        data_gene_score (object): Pandas data frame of genes' scores.

    raises:

    returns:
        (object): Pandas data frame of genes' scores.

    """

    def match_gene_name(identifier):
        return data_gene_annotation.loc[identifier, "gene_name"]
    data_gene_score.reset_index(level=["gene"], inplace=True)
    data_gene_score["name"] = (
        data_gene_score["gene"].apply(match_gene_name)
    )
    data_gene_score.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )
    return data_gene_score





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
    path_organization = os.path.join(dock, "organization")
    utility.confirm_path_directory(path_organization)
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    path_aggregation = os.path.join(
        path_organization, "data_gene_signal_aggregation.pickle"
    )
    path_log = os.path.join(
        path_organization, "data_gene_signal_log.pickle"
    )
    # Write information to file.
    with open(path_imputation, "wb") as file_product:
        pickle.dump(
            information["data_gene_signal_tissue_median"], file_product
        )
    with open(path_aggregation, "wb") as file_product:
        pickle.dump(information["data_gene_signal_aggregation"], file_product)
    with open(path_log, "wb") as file_product:
        pickle.dump(information["data_gene_signal_log"], file_product)


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

    # Each shuffle matrix should look something like this...
    # dimensions of the shuffle matrix should match the count of tissues and patients in the actual gene matrices
    # read in the gene matrices to find these counts.

    # matrix before shuffle... template
    # tissue   adipose  blood    colon
    # index
    # 1        1        1        1
    # 2        2        2        2
    # 3        3        3        3
    # 4        4        4        4
    # 5        5        5        5
    # 6        6        6        6
    # 7        7        7        7
    # 8        8        8        8
    # 9        9        9        9
    # 10       10       10       10

    # matrix after shuffle... template
    # tissue   adipose  blood    colon
    # index
    # 1        6
    # 2        3
    # 3        2
    # 4        7
    # 5        10
    # 6        4
    # 7        9
    # 8        1
    # 9        8
    # 10       5

    # when I apply the shuffle indices to the actual gene signal matrix, I might need to reset the gene signal matrix to have a numeric index
    # Consider doing this in split when setting up the gene matrices originally.
    # Consider changing the "patient" to a column instead of an index.
    # Then define some sort of mapping function to map a gene's matrix to values using the shuffle coordinates in the shuffle matrices...


    # Compile information.
    information = {}
    #Write product information to file.
    #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
