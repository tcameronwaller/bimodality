
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
    path_report = os.path.join(
        path_tissue, "data_gene_sample_report.txt"
    )
    # Read information from file.
    data_gene_sample_report = pandas.read_csv(
        path_report,
        sep="\t",
        header=0,
    )
    # Compile and return information.
    return {
        "data_gene_sample_report": data_gene_sample_report,
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
    data = source["data_gene_sample_report"]
    data.set_index(
        ["sample", "person", "tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data)




    pass


if (__name__ == "__main__"):
    execute_procedure()
