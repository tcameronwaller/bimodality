
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
import statsmodels.api

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
    #source = read_source(dock=dock)

    nsample = 100
    x = numpy.linspace(0, 10, 100)
    print(x)
    X = numpy.column_stack((x, x**2))
    print(X)
    beta = numpy.array([1, 0.1, 10])
    e = numpy.random.normal(size=nsample)
    X = statsmodels.api.add_constant(X)
    y = numpy.dot(X, beta) + e

    print(X)
    print(y)








    pass


if (__name__ == "__main__"):
    execute_procedure()
