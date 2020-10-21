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
import copy
import random
import itertools

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import distribution
import plot
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality



##########
# Source.


def read_source(dock=None):
    """
    Reads and organizes source information from file.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_access = os.path.join(dock, "access")
    path_data_pairs = os.path.join(
        path_access, "GSE95014_Hap1.validPairs.txt.gz"
    )

    path_structure = os.path.join(dock, "structure")
    path_data_pairs_sort = os.path.join(
        path_structure, "data_pairs_sort_raw.txt.gz"
    )

    # Read information from file.
    # As file is very large, read in chunks.
    reader_data_pairs = pandas.read_csv(
        path_data_pairs,
        compression="gzip",
        sep="\s+",
        header=None,
        names=[
            "read",
            "chromosome_one",
            "position_one",
            "direction_one",
            "chromosome_two",
            "position_two",
            "direction_two",
            "quality_one",
        ],
        dtype={
            "read": str,
            "chromosome_one": str,
            "position_one": numpy.uint64,
            "direction_one": str,
            "chromosome_two": str,
            "position_two": numpy.uint64,
            "direction_two": str,
            "quality_one": numpy.int32,
        },
        usecols=[
            #"read",
            "chromosome_one",
            "position_one",
            "direction_one",
            "chromosome_two",
            "position_two",
            "direction_two",
            #"quality_one",
        ],
        low_memory=True,
        chunksize=10000000,
        #nrows=100000,
    )
    # Compile and return information.
    return {
        "reader_data_pairs": reader_data_pairs,
    }


def translate_strand(direction=None):
    """
    Translates annotations of strand direction.

    arguments:
        direction (str): annotation of strand direction

    raises:

    returns:
        (int): annotation of strand direction

    """

    if direction == "+":
        strand = 0
    elif direction == "-":
        strand = 1
    return strand


def organize_data_pairs_chunk(
    data_chunk=None,
):
    """
    Organizes chunk of data.

    arguments:
        data_gene_annotation (object): Pandas data frame of information about
            pairs of genomic contacts

    raises:

    returns:
        (object): Pandas data frame of information about pairs of genomic
            contacts

    """

    # Organize data.
    data_chunk = data_chunk.copy(deep=True)
    # Insert dummy values for restriction fragments.
    # Juicer will ignore these values without a fragment library.
    # Insert dummy values for mapping quality.
    # Juicer will also ignore these values without a quality threshold.
    data_chunk["fragment_one"] = 0
    data_chunk["fragment_two"] = 1
    #data_chunk["quality_two"] = data_chunk["quality_one"]
    # Translate designations of strand direction.
    data_chunk["strand_one"] = data_chunk["direction_one"].apply(
        lambda value: translate_strand(direction=value)
    )
    data_chunk["strand_two"] = data_chunk["direction_two"].apply(
        lambda value: translate_strand(direction=value)
    )
    # Organize order of columns.
    data_chunk = data_chunk[[
        #"read",
        "strand_one",
        "chromosome_one",
        "position_one",
        "fragment_one",
        "strand_two",
        "chromosome_two",
        "position_two",
        "fragment_two",
        #"quality_one",
        #"quality_two",
    ]]
    # Return information.
    return data_chunk


# Write.


def write_product_chunk(dock=None, information=None):
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
    path_structure = os.path.join(dock, "structure")
    utility.create_directory(path_structure)

    path_data_pairs_text = os.path.join(
        path_structure, "data_pairs_sort_format_pandas.txt.gz"
    )

    # Write information to file.
    information["data_chunk"].to_csv(
        path_or_buf=path_data_pairs_text,
        columns=None,
        sep="\t",
        na_rep="",
        header=False,
        index=False,
        mode="a",
        compression="gzip",
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

    utility.print_terminal_partition(level=1)
    path_data_pairs_text = os.path.join(
        path_structure, "data_pairs_sort_format_pandas.txt.gz"
    )
    # Remove any previous version of file explicitly to avoid appending to
    # previous file.
    utility.remove_file(path=path_data_pairs_text)

    # Read source information from file.
    source = read_source(dock=dock)

    # Iterate on chunks of data.
    for data_chunk in source["reader_data_pairs"]:
        # Organize data chunk.
        data_chunk_novel = organize_data_pairs_chunk(
            data_chunk=data_chunk,
        )
        if False:
            utility.print_terminal_partition(level=2)
            print("data before reorganization...")
            print(data_chunk)
            utility.print_terminal_partition(level=2)
            print("data after reorganization...")
            print(data_chunk_novel)
        # Write data chunk to file.
        information = dict()
        information["data_chunk"] = data_chunk_novel
        write_product_chunk(
            dock=dock,
            information=information,
        )

    pass


if (__name__ == "__main__"):
    execute_procedure()
