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

import distribution
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
    path_persons = os.path.join(
        path_split, "persons.pickle"
    )
    path_tissues = os.path.join(
        path_split, "tissues.pickle"
    )
    # Read information from file.
    with open(path_persons, "rb") as file_source:
        persons = pickle.load(file_source)
    with open(path_tissues, "rb") as file_source:
        tissues = pickle.load(file_source)
    # Compile and return information.
    return {
        "persons": persons,
        "tissues": tissues,
    }


def create_matrix(
    dimension_zero=None,
    dimension_one=None
):
    """
    Creates matrix of indices.

    arguments:
        dimension_zero (int): count of indices in the first dimension or axis
        dimension_one (int): count of indices in the second dimension or axis

    raises:

    returns:
        (list<list<int>>): matrix of indices

    """

    # Base indices on zero.
    range_zero = list(range(dimension_zero))
    range_one = list(range(dimension_one))
    matrix = list()
    for index in range_zero:
        matrix.append(copy.copy(range_one))
    return matrix


def copy_matrix_shuffle_values(
    matrix=None
):
    """
    Copies a matrix and shuffles its values.

    arguments:
        matrix (list<list<int>>): matrix of values

    raises:

    returns:
        (list<list<list<int>>>): matrices of values

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
        count (int): count of matrices to create
        matrix (list<list<int>>): template matrix of indices

    raises:

    returns:
        (list<list<list<int>>>): matrices of indices

    """

    matrices = list()
    for index in range(count):
        matrix_copy = copy_matrix_shuffle_values(matrix=matrix)
        matrices.append(matrix_copy)
    return matrices


def create_shuffle_indices(
    count=None,
    dimension_zero=None,
    dimension_one=None
):
    """
    Creates matrices of indices to shuffle values.

    Dimensions of the shuffle matrix should match the count of tissues and
    persons in the actual matrices for genes' signals.

    arguments:
        count (int): count of matrices to create
        dimension_zero (int): count of indices in the first dimension or axis
        dimension_one (int): count of indices in the second dimension or axis

    raises:

    returns:
        (list<list<list<int>>>): matrices of indices

    """

    # Structure the matrices of shuffle indices as lists of lists.
    # The implicit semantics of data frames bind values within rows across
    # multiple columns.
    # Use a simple list matrix instead to make shuffles convenient.

    # Generate template matrix for permutations of persons and tissues.
    matrix = create_matrix(
        dimension_zero=dimension_zero,
        dimension_one=dimension_one
    )

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
    data = pandas.DataFrame(data=matrices[0])
    print(data)

    print(matrices[0][5][10])
    print(data.loc[5, 10])

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
    utility.create_directory(path_shuffle)
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

    # Remove previous files to avoid version or batch confusion.
    path_shuffle = os.path.join(dock, "shuffle")
    utility.remove_directory(path=path_shuffle)


    # Read source information from file.
    source = read_source(dock=dock)

    # Report.
    utility.print_terminal_partition(level=3)
    print(
        "Creating " + str(count) + " shuffles for matrices of dimension " +
        "zero: " + str(source["tissues"]) + " by dimension one: " +
        str(source["persons"]) + ". "
        "Notice that shuffles occur across dimension one (tissues for each " +
        "person)."
    )
    print(
        "Hence, values will stay matched to their respective tissues, but " +
        "they will be shuffled with respect to persons."
    )
    utility.print_terminal_partition(level=3)

    # Create shuffle indices.
    shuffles = create_shuffle_indices(
        count=count,
        dimension_zero=source["tissues"],
        dimension_one=source["persons"],
    )

    # Compile information.
    information = {
        "shuffles": shuffles
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
