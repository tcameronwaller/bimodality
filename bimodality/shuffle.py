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
import random
import copy

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import metric
import pipe
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
    path_signal = os.path.join(
        path_split, "genes_signals_patients_tissues.pickle"
    )
    # Read information from file.
    genes_signals_patients_tissues = pandas.read_pickle(path_signal)
    # Compile and return information.
    return {
        "genes_signals_patients_tissues": genes_signals_patients_tissues,
    }


def create_matrix(
    width=None,
    length=None
):
    """
    Creates matrix of indices.

    arguments:
        width (int): Count of indices in the first level axis, columns.
        length (int): Count of indices in the second level axis, rows.

    raises:

    returns:
        (list<list<int>>): Matrix of indices.

    """

    # Base indices on zero.
    range_width = list(range(width))
    range_length = list(range(length))
    matrix = list()
    for index in range_width:
        matrix.append(copy.copy(range_length))
    return matrix


def copy_matrix_shuffle_values(
    matrix=None
):
    """
    Copies a matrix and shuffles its values.

    arguments:
        matrix (list<list<int>>): Matrix of values.

    raises:

    returns:
        (list<list<list<int>>>): Matrices of values.

    """

    # Copy values in matrix.
    matrix_copy = copy.deepcopy(matrix)
    # Shuffle values in matrix.
    for index in range(len(matrix_copy)):
        random.shuffle(matrix_copy[index])
    # Return matrix.
    return matrix_copy


def copy_matrices_shuffle_values(
    count=None,
    matrix=None
):
    """
    Copies matrices and shuffles values.

    arguments:
        count (int): Count of matrices to create.
        matrix (list<list<int>>): Template matrix of indices.

    raises:

    returns:
        (list<list<list<int>>>): Matrices of indices.

    """

    matrices = list()
    for index in range(count):
        matrix_copy = copy_matrix_shuffle_values(matrix=matrix)
        matrices.append(matrix_copy)
    return matrices


def create_shuffle_indices(
    count=None,
    width=None,
    length=None
):
    """
    Creates matrices of indices to shuffle values.

    arguments:
        count (int): Count of matrices to create.
        width (int): Count of indices in the first level axis, columns.
        length (int): Count of indices in the second level axis, rows.

    raises:

    returns:
        (list<list<list<int>>>): Matrices of indices.

    """

    # Dimensions of the shuffle matrix should match the count of tissues and
    # patients in the actual matrices for genes' signals.

    # Structure the matrices of shuffle indices as lists of lists.
    # The implicit semantics of data frames bind values within rows across
    # multiple columns.
    # Use a simple list matrix instead to make shuffles convenient.
    # The top dimension should be for tissues to allow convenient shuffles of
    # patients for each tissue.

    # Generate template matrix for permutations of patients and tissues.
    matrix = create_matrix(width=width, length=length)

    # Create and shuffle copies of the template matrix.
    matrices = copy_matrices_shuffle_values(count=count, matrix=matrix)

    # Summarize.
    print(
        "Dimensions of each matrix: " +
        str(len(matrices[0])) +
        " by " +
        str(len(matrices[0][0]))
    )
    print("Created and shuffled matrices: " + str(len(matrices)))

    # Print representation of first matrix.
    # Convert matrix to a data frame for a clear representation.
    data = pandas.DataFrame(data=matrices[0]).transpose()
    print(data)

    # Demonstrate access to values.
    # Matrix structure has tissues on first axis and patients on second axis.
    print("Access individual value from first matrix.")
    print(
        "Matrix structure of shuffle indices has tissues on first axis and " +
        "patients on second axis."
    )
    print(
        "Matrix structure of genes' signals has tissues on columns and " +
        "patients on rows."
    )
    print(
        "Notice differences in access for matrices of indices versus signals."
    )
    print(matrices[0][5][10])
    print(data.loc[10, 5])

    # Return information.
    return matrices


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
    path_shuffle = os.path.join(dock, "shuffle")
    utility.confirm_path_directory(path_shuffle)
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    # Write information to file.
    with open(path_shuffles, "wb") as file_product:
        pickle.dump(information["shuffles"], file_product)


###############################################################################
# Procedure


def execute_procedure(dock=None, count=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        count (int): count of shuffles to create and store

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(dock=dock)

    # Determine dimensions of the patient-tissue signal matrix for each gene.
    print(source["genes_signals_patients_tissues"]["ENSG00000029363"].shape)
    count_tissues = (
        source["genes_signals_patients_tissues"]["ENSG00000029363"].shape[1]
    )
    count_patients = (
        source["genes_signals_patients_tissues"]["ENSG00000029363"].shape[0]
    )

    # Report.
    utility.print_terminal_partition(level=3)
    print(
        "Creating " + str(count) + " shuffles for matrices of dimensions " +
        str(count_tissues) + " by " + str(count_patients) + "."
    )
    utility.print_terminal_partition(level=3)

    # Create shuffle indices.
    shuffles = create_shuffle_indices(
        count=count,
        width=count_tissues,
        length=count_patients
    )

    # when I apply the shuffle indices to the actual gene signal matrix, I might need to reset the gene signal matrix to have a numeric index
    # Consider doing this in split when setting up the gene matrices originally.
    # Consider changing the "patient" to a column instead of an index.
    # Then define some sort of mapping function to map a gene's matrix to values using the shuffle coordinates in the shuffle matrices...


    # Compile information.
    information = {
        "shuffles": shuffles
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
