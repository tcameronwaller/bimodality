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


# convert dataframe to matrix using dataframe.values
# shuffle the matrix by indices using a map function as in... https://stackoverflow.com/questions/26194389/how-to-rearrange-array-based-upon-index-array
# assign the shuffled matrix back to the same persons indices and tissues columns in a new pandas dataframe using the Pandas dataframe constructor

#sorted_matrix = np.array(list(map(
    #lambda index, sort: sort[index], index_matrix, matrix_to_sort
#)))



def shuffle_gene_signals_iterative(
    data_gene_signals=None,
    shuffle=None
):
    """
    Shuffles the association of gene's signals between tissues and persons.

    This method is iterative and inefficient.

    This method requires dimension zero to correspond to tissues and dimension
    one to correspond to persons.

    arguments:
        data_gene_signals (object): Pandas data frame of a gene's signals
            across persons and tissues
        shuffle (list<list<int>>): matrix of indices

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Determine values and counts of tissues and persons.
    tissues = data_gene_signals.columns.to_list()
    persons = data_gene_signals.index.to_list()
    count_tissues = len(tissues)
    count_persons = len(persons)

    # Collect gene's signals.
    gene_signals_shuffle = list()
    for index_person in range(count_persons):
        person = persons[index_person]
        record = dict()
        record["person"] = person
        for index_tissue in range(count_tissues):
            tissue = tissues[index_tissue]
            index_person_shuffle = shuffle[index_tissue][index_person]
            signal = (
                data_gene_signals.iloc[index_person_shuffle, index_tissue]
            )
            record[tissue] = signal
        gene_signals_shuffle.append(record)

    # Convert records to data frame.
    data = utility.convert_records_to_dataframe(
        records=gene_signals_shuffle
    )
    data.rename_axis(
        columns="tissue",
        axis="columns",
        copy=False,
        inplace=True
    )
    data.set_index(
        "person",
        drop=True,
        inplace=True,
    )
    return data


def shuffle_gene_signals(
    data_gene_signals=None,
    shuffle=None
):
    """
    Shuffles the association of gene's signals between tissues and persons.

    This method requires less memory and less process time than the iterative
    method.

    This method requires dimension zero to correspond to persons and dimension
    one to correspond to tissues.

    arguments:
        data_gene_signals (object): Pandas data frame of a gene's signals
            across persons and tissues
        shuffle (list<list<int>>): matrix of indices

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Copy data.
    data = data_gene_signals.copy(deep=True)
    # Determine values and counts of tissues and persons.
    tissues = data.columns.to_list()
    persons = data.index.to_list()
    # Convert data to matrix.
    # Transpose data matrix to match shuffle dimension.
    # Shuffle occurs on dimension one (persons).
    matrix_data = numpy.transpose(data.values)
    #print("data matrix before shuffle")
    #print(pandas.DataFrame(data=matrix_data))
    #print(matrix_data.shape)

    # Organize shuffle matrix.
    matrix_shuffle = numpy.array(shuffle)
    #print("shuffle matrix before shuffle")
    #print(pandas.DataFrame(data=matrix_shuffle))
    #print(matrix_shuffle.shape)

    # Sort data matrix's values by shuffle matrix's indices.
    matrix_data_sort = numpy.array(list(map(
        lambda index, sort: sort[index], matrix_shuffle, matrix_data
    )))
    matrix_sort_transpose = numpy.transpose(matrix_data_sort)

    # Convert records to data frame.
    data_sort = pandas.DataFrame(
        data=matrix_sort_transpose,
        index=persons,
        columns=tissues,
    )
    data_sort.rename_axis(
        columns="tissue",
        axis="columns",
        copy=False,
        inplace=True
    )
    data_sort.rename_axis(
        index="person",
        axis="index",
        copy=False,
        inplace=True
    )
    return data_sort


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
    path_permutation = os.path.join(dock, "permutation")
    utility.create_directory(path_permutation)
    path_permutations = os.path.join(
        path_permutation, "permutations.pickle"
    )
    # Write information to file.
    with open(path_permutations, "wb") as file_product:
        pickle.dump(information["permutations"], file_product)


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

    # Report.
    utility.print_terminal_partition(level=3)
    print(
        "Creating " + str(count) + " shuffles for matrices of dimension " +
        "zero: " + str(source["tissues"]) + " by dimension one: " +
        str(source["persons"]) + ". "
        "Notice that shuffles occur across dimension one."
    )
    utility.print_terminal_partition(level=3)

    # Create shuffle indices.
    permutations = create_shuffle_indices(
        count=count,
        dimension_zero=source["tissues"],
        dimension_one=source["persons"],
    )

    # Compile information.
    information = {
        "permutations": permutations
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
