"""
...

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import statistics
import pickle

# Relevant

import numpy
import pandas
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

# Custom

import measurement
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
    path_split = os.path.join(dock, "split")
    path_demo = os.path.join(
        path_split, "demonstration_gene_samples_signals.pickle"
    )
    # Read information from file.
    data_gene_samples_signals = pandas.read_pickle(path_demo)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
    }



##########
# Report.




##########
# Product.


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
    utility.create_directory(path_organization)
    path_gene_signal = os.path.join(
        path_organization, "data_gene_persons_tissues_signals.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_gene_persons_tissues_signals"], path_gene_signal
    )

    pass


###############################################################################
# Procedure


def test(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    utility.print_terminal_partition(level=1)
    print("Welcome to the test, demonstration of the organization procedure!")

    # Read source information from file.
    source = read_source(dock=dock)
    print(source["data_gene_samples_signals"])

    # Call procedure.
    information = execute_procedure(
        data_gene_samples_signals=source["data_gene_samples_signals"]
    )

    print("here is the gene dispersion report...")
    print(information["report_gene"]["data_gene_tissue_dispersion"])

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


def execute_procedure(
    data_gene_samples_signals=None
):
    """
    Function to execute module's main behavior.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): report about gene, Pandas data frame of gene's signals across
            persons and tissues

    """

    # Prepare gene report.

    # Compile information.
    information = {
        "report_gene": report_gene,
        "data_gene_persons_tissues_signals": data_gene_persons_tissues_signals,
    }
    # Return information.
    return information


if (__name__ == "__main__"):
    execute_procedure()
