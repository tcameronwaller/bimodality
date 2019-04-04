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


def split_genes_signals(data_gene_signal=None):
    """
    Splits genes' signals across patients and tissues by genes.

    arguments:
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific patients and tissues.

    raises:

    returns:
        (dict<object>): Collection of matrices.

    """

    # Split genes' signals across tissues and patients by gene.
    utility.print_terminal_partition(level=2)
    print(
        "Split of genes' signals across patients and tissues by gene."
    )
    print(
        "Organize patients-tissues matrices for each gene within a " +
        "dictionary."
    )

    # Change data's structure.
    data_long = data_gene_signal.stack("gene").to_frame(name="value")
    print(data_long)
    print(data_long.shape)
    data_tissue = data_long.unstack(level="tissue")
    data_tissue.columns = data_tissue.columns.get_level_values(1)
    print(data_tissue)
    print(data_tissue.shape)
    data_tissue.reset_index(level="gene", inplace=True)
    groups = data_tissue.groupby("gene")
    # Collect matrices for each gene.
    genes_signals = dict()
    for name, group in groups:
        data = group.copy(deep=True)
        data.drop(
            labels="gene",
            axis="columns",
            inplace=True
        )
        genes_signals[name] = data
    # Print single matrix.
    print("Access to matrix for a single gene.")
    print(genes_signals["ENSG00000029363"])
    print("Access to value for single patient and tissue for a single gene.")
    print(genes_signals["ENSG00000029363"].loc["GTEX-117YW", "Blood"])

    # Return information.
    return genes_signals


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
    path_split = os.path.join(dock, "split")
    utility.confirm_path_directory(path_split)
    path_gene = os.path.join(
        path_split, "genes.txt"
    )
    path_signal = os.path.join(
        path_split, "genes_signals_patients_tissues.pickle"
    )
    # Write information to file.
    utility.write_file_text_list(
        information=information["genes"], path_file=path_gene
    )
    pandas.to_pickle(
        information["genes_signals_patients_tissues"], path_signal
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
    print(source["data_gene_signal_standard"].iloc[0:10, 0:10])

    # Split genes' signals across tissues and patients by gene.
    genes_signals_patients_tissues = split_genes_signals(
        data_gene_signal=source["data_gene_signal_standard"]
    )

    # Organize genes' identifiers.
    # Format of genes' identifiers needs to be readable by Bash as an array.
    genes = source["data_gene_signal_standard"].columns.to_list()
    print("count of genes: " + str(len(genes)))

    # Compile information.
    information = {
        "genes": genes,
        "genes_signals_patients_tissues": genes_signals_patients_tissues,
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
