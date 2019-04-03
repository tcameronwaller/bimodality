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

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index(
        ["patient", "tissue"], append=False, drop=True
    )
    # Define groups.
    groups = data_gene_signal_index.groupby(level="patient")
    # Apply calculation to groups.
    data_sum_index = groups.aggregate(lambda x: x.sum())
    data_sum = (
        data_sum_index.reset_index(level=["patient"])
    )
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
    if False:
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

    ########################################################




    # Calculate sum of genes' signals in z-score space across all tissues for
    # each patient.
    utility.print_terminal_partition(level=2)
    print(
        "Calculation of sum of genes' signals in z-score space across " +
        "tissues for each patient."
    )
    data_gene_signal_sum = calculate_sum_score_gene_signal_by_patient(
        data_gene_signal=data_gene_signal_standard
    )
    data_sum_index = data_gene_signal_sum.set_index(
        ["patient"], append=False, drop=True
    )
    print(data_sum_index.iloc[0:10, 0:10])
    utility.print_terminal_partition(level=3)
    #print(data_sum_index.loc[:, "ENSG00000240453.1"])
    # Gene ENSG00000240453.1 has values of 0.0 for all patients.



    # Reshape data with patients as columns and genes as rows.
    if False:
        utility.print_terminal_partition(level=2)
        print(
            "Reshape of data with patients as columns and genes as rows."
        )
        data_gene_signal_axis = data_gene_signal_sum.rename_axis("gene", axis="columns")
        data_gene_signal_axis_index = data_gene_signal_axis.set_index(
            ["patient"], append=False, drop=True
        )
        print(data_gene_signal_axis_index.iloc[0:10, 0:10])
        data_gene_signal_shape = data_gene_signal_axis_index.transpose(copy=True)
        print(data_gene_signal_shape.iloc[:, :])

        print("I will want to apply the bimodality test to each column of genes... columns need to be genes...")

    # Calculate metrics of modality.
    data_bimodality = data_sum_index.aggregate(metric.calculate_bimodality_coefficient, axis="index")
    print(data_bimodality)


    data_dip = data_sum_index.aggregate(metric.calculate_dip_statistic, axis="index")
    print(data_dip)














    # Compile information.
    information = {}
    #Write product information to file.
    #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
