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
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_selection = os.path.join(dock, "selection")
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_signal": data_gene_signal,
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
        "data frame."
    )
    print(
        "Organize patients as rows and tissues as columns."
    )
    print(
        "In subsequent analyses, patients are higher in hierarchy than " +
        "tissues."
    )

    # Change data's structure.
    data_long = data_gene_signal.stack("gene").to_frame(name="value")
    print(data_long)
    print(data_long.shape)
    data_stack = data_long.unstack(level="tissue")
    data_stack.columns = data_stack.columns.get_level_values(1)
    print(data_stack)
    print(data_stack.shape)
    data_stack.reset_index(level="gene", inplace=True)
    groups = data_stack.groupby("gene")
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
    print(source["data_gene_signal"].iloc[0:10, 0:10])

    # Before I split... I need to aggregate the tissue categories... I need the same major tissues for every patient...
    # Maybe I can go ahead and aggregate major tissue types within the tissue procedure?
    # When I do that, I should also describe the variance of each gene across minor tissue types...

    if False:

        # Consider using pandas.DataFrame.round(decimals=10)

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
