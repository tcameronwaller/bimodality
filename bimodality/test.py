
"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import pickle
import statistics

# Relevant

import numpy
import pandas
import sklearn
import sklearn.datasets
import sklearn.decomposition

# Custom

import organization
import selection
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
    gene = "ENSG00000231925" # TAPBP
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
    }



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
    data = source["data_gene_samples_signals"]
    print(data)

    gene = "ENSG00000231925" # TAPBP

    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data
    )
    data_gene_persons_tissues_signals = collection_organization["data_gene_persons_tissues_signals"]

    print(data_gene_persons_tissues_signals)
    print(data_gene_persons_tissues_signals.shape)

    # Filter genes by signal.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_row = selection.filter_rows_columns_by_threshold_proportion(
        data=data_gene_persons_tissues_signals,
        dimension="row",
        threshold=10.0,
        proportion=0.1
    )
    print(data_row)
    print(data_row.shape)

    # Filter samples by signal.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_column = selection.filter_rows_columns_by_threshold_proportion(
        data=data_gene_persons_tissues_signals,
        dimension="column",
        threshold=10.0,
        proportion=0.5
    )
    print(data_column)
    print(data_column.shape)






    pass


if (__name__ == "__main__"):
    execute_procedure()
