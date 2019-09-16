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
import scipy
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

# Custom

import measurement
import permutation
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
    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    path_organization = os.path.join(dock, "organization")
    path_organization_signal = os.path.join(
        path_organization, "data_gene_persons_tissues_signals.pickle"
    )
    path_restriction = os.path.join(dock, "restriction")
    path_restriction_signal = os.path.join(
        path_restriction, "data_gene_persons_tissues_signals.pickle"
    )
    # Read information from file.
    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    data_organization_signals = pandas.read_pickle(path_organization_signal)
    data_restriction_signals = pandas.read_pickle(path_restriction_signal)
    # Compile and return information.
    return {
        "shuffles": shuffles,
        "data_organization_signals": data_organization_signals,
        "data_restriction_signals": data_restriction_signals,
    }



##########
# Distribution.


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
    print("Welcome to the test, demonstration of the distribution procedure!")

    # Read source information from file.
    source = read_source(dock=dock)
    print(source["data_restriction_signals"])

    # Call procedure.
    information = execute_procedure(
        data_gene_persons_tissues_signals=(
            source["data_restriction_signals"]
        ),
    )

    print("Printing distribution values...")
    print(information["values"])

    print("Printing bimodality scores...")
    print(information["scores"])

    print("Printing gene report... specifically data_gene_persons_signals")
    print(information["report_gene"]["data_gene_persons_signals"])

    ###########################################################################

    utility.print_terminal_partition(level=1)
    print("Now let's test the shuffle procedure...")

    utility.print_terminal_partition(level=2)
    print("Here's the original matrix before shuffle.")
    matrix = shuffle.create_matrix(
        dimension_zero=20,
        dimension_one=406
    )
    # Print representation of matrix.
    # Convert matrix to a data frame for a clear representation.
    data = pandas.DataFrame(data=matrix)
    print(data)

    # Transpose matrix to match data dimensions.
    matrix_array = numpy.array(matrix)
    matrix_transpose = numpy.transpose(matrix_array)
    # Print representation of matrix.
    # Convert matrix to a data frame for a clear representation.
    data = pandas.DataFrame(data=matrix_transpose)
    print(data)

    # Create and shuffle copies of the template matrix.
    utility.print_terminal_partition(level=2)
    matrices = shuffle.copy_matrices_shuffle_values(count=10, matrix=matrix)
    # Print representation of matrix.
    # Convert matrix to a data frame for a clear representation.
    matrix_shuffle = matrices[0]
    print("raw matrix after shuffle")
    print(matrix_shuffle)
    matrix_shuffle_transpose = numpy.transpose(numpy.array(matrix_shuffle))
    data = pandas.DataFrame(data=matrix_shuffle_transpose)
    print("transposition of the matrix after shuffle")
    print(data)

    ###########################################################################

    # Check that shuffles occur on proper axis.
    utility.print_terminal_partition(level=1)
    print("normalize and standardize data... then shuffle")
    print("Notice that shuffle should not change z-score mean and std.")

    # Normalize and standardize gene's signals.
    # Transform gene's signals to base-two logarithmic space.
    # Transform gene's signals to standard, z-score space.
    data_normal_standard = normalize_standardize_gene_signal(
        data_gene_persons_tissues_signals=source["data_restriction_signals"],
    )

    # Compare summary statistics before and after transformation.
    utility.print_terminal_partition(level=3)
    print("data before shuffle...")
    print(data_normal_standard.iloc[0:10, 0:10])
    data_mean = data_normal_standard.aggregate(
        lambda x: x.mean(),
        axis="index"
    )
    print("Mean")
    print(data_mean.iloc[0:10])
    data_deviation = data_normal_standard.aggregate(
        lambda x: x.std(),
        axis="index"
    )
    print("Standard deviation")
    print(data_deviation.iloc[0:10])

    # Shuffle gene's signals.
    utility.print_terminal_partition(level=3)
    data_shuffle = shuffle.shuffle_gene_signals(
        data_gene_signals=data_normal_standard,
        shuffle=matrices[0]
    )
    print("data after shuffle...")
    print(data_shuffle.iloc[0:10, 0:10])
    data_mean = data_shuffle.aggregate(
        lambda x: x.mean(),
        axis="index"
    )
    print("Mean")
    print(data_mean.iloc[0:10])
    data_deviation = data_shuffle.aggregate(
        lambda x: x.std(),
        axis="index"
    )
    print("Standard deviation")
    print(data_deviation.iloc[0:10])


    pass


if (__name__ == "__main__"):
    execute_procedure()
