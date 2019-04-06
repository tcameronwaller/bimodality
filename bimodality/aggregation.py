"""
...

"""

###############################################################################
# Notes

# TODO: move this procedure to the organization procedure... ?


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
    # Read information from file.
    data_gene_signal_mean = pandas.read_pickle(path_mean)
    data_gene_signal_median = pandas.read_pickle(path_median)
    data_gene_signal_imputation = pandas.read_pickle(path_imputation)
    data_gene_signal_log = pandas.read_pickle(path_log)
    data_gene_signal_standard = pandas.read_pickle(path_standard)
    # Compile and return information.
    return {
        "data_gene_signal_mean": data_gene_signal_mean,
        "data_gene_signal_median": data_gene_signal_median,
        "data_gene_signal_imputation": data_gene_signal_imputation,
        "data_gene_signal_log": data_gene_signal_log,
        "data_gene_signal_standard": data_gene_signal_standard,
    }


def calculate_sum_score_gene_signal_by_patient(data_gene_signal=None):
    """
    Calculates the sum of genes' signals across tissues for each patient.

    arguments:
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of sum of standard score signals for all
            genes across tissues for each patient.

    """

    # Drop the tissue designations.
    data_gene_signal.reset_index(level=["tissue"], inplace=True)
    data_gene_signal.drop(
        labels="tissue",
        axis="columns",
        inplace=True
    )
    # Define groups.
    groups = data_gene_signal.groupby(level="patient")
    # Apply calculation to groups.
    data_sum = groups.aggregate(lambda x: x.sum())
    return data_sum


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
    path_aggregation = os.path.join(dock, "aggregation")
    utility.confirm_path_directory(path_aggregation)
    path_signal = os.path.join(
        path_aggregation, "data_gene_signal_aggregation.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_gene_signal_aggregation"], path_signal
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
    # Summary.

    print(source["data_gene_signal_standard"].iloc[0:25, 0:10])

    # Calculate sum of genes' signals in standard, z-score space across all
    # tissues for each patient.
    utility.print_terminal_partition(level=2)
    print(
        "Calculation of sum of genes' signals in z-score space across " +
        "tissues for each patient."
    )
    data_gene_signal_aggregation = calculate_sum_score_gene_signal_by_patient(
        data_gene_signal=source["data_gene_signal_standard"]
    )
    print(data_gene_signal_aggregation.iloc[0:10, 0:10])

    # Compile information.
    information = {
        "data_gene_signal_aggregation": data_gene_signal_aggregation,
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
