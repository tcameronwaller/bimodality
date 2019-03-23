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
    path_access = os.path.join(dock, "access")
    path_organization = os.path.join(dock, "organization")
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    # Read information from file.
    data_gene_signal_imputation = pandas.read_pickle(path_imputation)
    # Compile and return information.
    return {
        "data_gene_signal_imputation": data_gene_signal_imputation,
    }


def calculate_logarithm_gene_signal(pseudo_count=None, data_gene_signal=None):
    """
    Calculates the base-2 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base-2 logarithm.

    arguments:
        pseudo_count (float): Pseudo count to add to gene signal before
            transformation to avoid values of zero.
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic signals for all genes
            across specific patients and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    # An alternative approach would be to set the label columns to indices.
    if False:
        data_gene_signal_index = data_gene_signal.set_index(
            ["Name", "Description"], append=True, drop=True
        )
        data_gene_signal_log = data_gene_signal_index.apply(
            lambda value: 2 * value
        )
        data_log = data_signal_index.copy()
        data_log.iloc[0:, 2:] = data_log.iloc[0:, 2:].applymap(
            lambda value: math.log((value + 1.0), 2)
        )

    data_index = data_gene_signal.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    data_log = data_index.applymap(
        lambda value: math.log((value + pseudo_count), 2)
    )
    # Reverse an index to a column.
    #dataframe["new_column"] = dataframe.index
    data = data_log.reset_index(level=["patient", "tissue"])
    return data


def calculate_standard_score_gene_signal_by_tissue(data_gene_signal=None):
    """
    Calculates the standard (z-score) of genes' signals for each patient and
    tissue.

    The standard scores are relative to tissue.
    The values of mean and standard deviation are across all patients for each
    tissue.

    arguments:
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific patients and tissues.

    raises:

    returns:
        (object): Pandas data frame of standard score signals for all genes
            across specific patients and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index(
        ["tissue", "patient"], append=False, drop=True
    )
    groups = data_gene_signal_index.groupby(level="tissue")
    # For some reason, the scipy function throws an error about illegal values.
    #groups.transform(lambda value: print(value.to_list()))
    #data_standard_index = (
    #    groups.transform(lambda value: scipy.stats.zscore(value.to_list()))
    #)
    data_standard_index = groups.transform(lambda x: (x - x.mean()) / x.std())
    data_standard = (
        data_standard_index.reset_index(level=["tissue", "patient"])
    )
    return data_standard


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
    data_gene_signal_imputation_index = (
        source["data_gene_signal_imputation"].set_index(
            ["patient", "tissue"], append=False, drop=True
        )
    )
    print(data_gene_signal_imputation_index.iloc[0:10, 0:10])
    if False:
        data_gene_signal_sort = (
            data_gene_signal.sort_index(axis="columns", ascending=False)
        )
    # Summarize tissues and patients.
    # Extract unique values of tissue.
    tissues = list(source["data_gene_signal_imputation"]["tissue"].unique())
    patients = list(source["data_gene_signal_imputation"]["patient"].unique())
    # Split data by tissue.
    data_gene_signal_imputation_index = (
        source["data_gene_signal_imputation"].set_index(
        ["tissue"], append=False, drop=True
        )
    )
    data_heart = (
        data_gene_signal_imputation_index
            .loc["Heart"].reset_index(level=["tissue"])
    )
    print(data_heart)
    #data_merger = data_skin.merge(data_heart, how="outer")


    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values. <- maybe false
    # Logarithms are not distributive.
    # To accommodate gene signals of 0.0, add a pseudo count of 2.0 to all
    # counts before calculation of the base-2 logarithm.
    # log-2 signal = log2(TPM + pseudo-count)
    # pseudo-count = 1.0
    utility.print_terminal_partition(level=2)
    print("Transformation of genes' signals to base-2 logarithmic space.")
    print(
        "To accommodate gene signals of 0.0, add a pseudo count of 2.0 to " +
        "all counts before calculation of the base-2 logarithm."
    )
    data_gene_signal_log = calculate_logarithm_gene_signal(
        pseudo_count=2.0,
        data_gene_signal=source["data_gene_signal_imputation"]
    )
    data_gene_signal_log_index = (
        data_gene_signal_log.set_index(
            ["patient", "tissue"], append=False, drop=True
        )
    )
    print(data_gene_signal_log_index.iloc[0:10, 0:10])

    # Transform signals to standard score space.
    utility.print_terminal_partition(level=2)
    print(
        "Transformation of genes' signals to standard score (z-score) " +
        "space."
    )
    data_gene_signal_standard = calculate_standard_score_gene_signal_by_tissue(
        data_gene_signal=data_gene_signal_log
    )
    data_standard_index = data_gene_signal_standard.set_index(
        ["tissue", "patient"], append=False, drop=True
    )
    print(data_standard_index.iloc[0:10, 0:10])
    # Compare summary statistics before and after transformation.
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals in base-2 logarithmic space.")
    groups = data_gene_signal_log_index.groupby(level="tissue")
    print(groups.describe())
    utility.print_terminal_partition(level=3)
    print("Summary statistics for gene signals in standard space.")
    groups_standard = data_standard_index.groupby(level="tissue")
    print(groups_standard.describe())
    print(groups_standard.mean())
    print(groups_standard.std())











    # Compile information.
    information = {}
    #Write product information to file.
    #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
