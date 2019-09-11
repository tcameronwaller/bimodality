
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
import wget
import gzip

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
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_expecto = os.path.join(dock, "expecto")
    path_columns = os.path.join(path_expecto, "colnames.final")
    path_data = os.path.join(path_expecto, "combined_snps.0.3.final.sorted")

    # Read information from file.
    if False:
        columns_list = utility.read_file_text_list(
            delimiter="\t",
            strip="\n",
            path_file=path_columns
        )

    columns = pandas.read_csv(
        path_columns,
        sep='\t'
    )
    columns_list = list(columns)
    data_expecto = pandas.read_csv(
        path_data,
        sep='\t',
        header=None,
        names=columns_list
    )

    # Compile and return information.
    return {
        "data_expecto": data_expecto,
    }


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
    path_expecto = os.path.join(dock, "expecto")
    utility.create_directory(path_expecto)
    path_data_feather = os.path.join(
        path_expecto, "data_expecto.feather"
    )
    path_data_pickle = os.path.join(
        path_expecto, "data_expecto.pickle"
    )

    # Write information to file.
    information["data_expecto"].to_feather(
        fname=path_data_feather,
    )
    pandas.to_pickle(
        information["data_expecto"],
        path=path_data_pickle,
        #compression="gzip",
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

    if False:
        # Remove previous files to avoid version or batch confusion.
        path_expecto = os.path.join(dock, "expecto")
        utility.remove_directory(path=path_expecto)

        path_remote = "http://deepsea.princeton.edu/media/code/expecto/combined_snps.0.3.zip"
        path_local = os.path.join(path_expecto, "combined_snps.0.3.zip")
        utility.remove_file(path_local)
        wget.download(path_remote, path_local)
        utility.decompress_file_gzip(path_local)

    # Read source information from file.
    source = read_source(dock=dock)

    print(source["data_expecto"])

    # Compile information.
    information = {
        "data_expecto": source["data_expecto"]
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
