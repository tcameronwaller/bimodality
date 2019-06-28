
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
    path_tissue = os.path.join(dock, "tissue")
    path_gene_tissue = os.path.join(
        path_tissue, "data_gene_tissue.pickle"
    )
    path_component = os.path.join(
        path_tissue, "data_gene_tissue_component.pickle"
    )
    # Read information from file.
    data_gene_tissue = pandas.read_pickle(
        path_gene_tissue
    )
    data_gene_tissue_component = pandas.read_pickle(path_component)
    # Compile and return information.
    return {
        "data_gene_tissue": data_gene_tissue,
        "data_gene_tissue_component": data_gene_tissue_component,
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

    # Specify directories and files.
    path_tissue = os.path.join(dock, "tissue")
    utility.confirm_path_directory(path_tissue)
    path_gene_tissue = os.path.join(
        path_tissue, "data_gene_tissue.txt"
    )
    path_component = os.path.join(
        path_tissue, "data_gene_tissue_component.txt"
    )
    # Write information to file.
    source["data_gene_tissue"].to_csv(
        path_or_buf=path_gene_tissue,
        sep="\t",
        header=True,
        index=False,
    )
    source["data_gene_tissue_component"].to_csv(
        path_or_buf=path_component,
        sep="\t",
        header=True,
        index=False,
    )



    pass


if (__name__ == "__main__"):
    execute_procedure()
