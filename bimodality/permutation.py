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
import functools
import multiprocessing
import datetime
import gc
import time
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


def read_source_initial(
    dock=None
):
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
    path_selection = os.path.join(dock, "selection", "tight")
    path_selection_genes = os.path.join(
        path_selection, "genes_selection.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_unimodal = os.path.join(
        path_candidacy, "genes_unimodal.pickle"
    )
    path_genes_multimodal = os.path.join(
        path_candidacy, "genes_multimodal.pickle"
    )

    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(path_shuffle, "shuffles.pickle")
    # Read information from file.
    with open(path_selection_genes, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_genes_unimodal, "rb") as file_source:
        genes_unimodal = pickle.load(file_source)
    with open(path_genes_multimodal, "rb") as file_source:
        genes_multimodal = pickle.load(file_source)

    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes_selection": genes_selection,
        "genes_unimodal": genes_unimodal,
        "genes_multimodal": genes_multimodal,
        "shuffles": shuffles,
    }


def read_source(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution", "genes")
    path_gene = os.path.join(path_distribution, gene)
    path_data_gene_persons_tissues_signals = os.path.join(
        path_gene, "data_gene_persons_tissues_signals.pickle"
    )
    # Read information from file.
    data_gene_persons_tissues_signals = pandas.read_pickle(
        path_data_gene_persons_tissues_signals,
    )
    # Compile and return information.
    return {
        "data_gene_persons_tissues_signals": data_gene_persons_tissues_signals,
    }


# Permute genes' signals


def permute_gene_distribution(
    gene=None,
    data_gene_persons_tissues_signals=None,
    shuffles=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues
        shuffles (list<list<list<int>>>): matrices of shuffle indices

    raises:

    returns:
        (dict<list<float>>): values of modality metrics across permutations

    """

    def permute_collect_measures(
        collection=None,
        shuffle_indices=None,
        data_gene_signals=None,
    ):
        #print(collection)
        # Shuffle gene's signals.
        data_shuffle = permute_gene_signals(
            data_gene_signals=data_gene_signals,
            shuffle_indices=shuffle_indices
        )
        # Distribution
        bins = distribution.prepare_describe_distribution(
            data_gene_persons_tissues_signals=data_shuffle,
        )
        # Collect metrics.
        collection["coefficients"].append(
            bins["bin_modality"]["scores"]["coefficient"]
        )
        collection["dips"].append(
            bins["bin_modality"]["scores"]["dip"]
        )
        collection["mixtures"].append(
            bins["bin_modality"]["scores"]["mixture"]
        )
        return collection

    # Initialize collections of measures.
    collection = dict()
    collection["coefficients"] = list()
    collection["mixtures"] = list()
    collection["dips"] = list()

    # Iterate on permutations.
    information = functools.reduce(
        lambda collection, shuffle_indices: permute_collect_measures(
            collection=collection,
            shuffle_indices=shuffle_indices,
            data_gene_signals=data_gene_persons_tissues_signals,
        ),
        shuffles, # Iterable
        collection, # Initializer
    )

    # Convert lists to NumPy arrays.
    information["coefficients"] = numpy.array(
        collection["coefficients"],
        copy=False,
    )
    information["mixtures"] = numpy.array(
        collection["mixtures"],
        copy=False,
    )
    information["dips"] = numpy.array(
        collection["dips"],
        copy=False,
    )

    # Return information.
    return information


def permute_gene_distribution_iterative(
    gene=None,
    data_gene_persons_tissues_signals=None,
    shuffles=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        gene (str): identifier of gene
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues
        shuffles (list<list<list<int>>>): matrices of shuffle indices

    raises:

    returns:
        (dict<list<float>>): values of modality metrics across permutations

    """

    # Initialize collections of measures.
    collection = dict()
    collection["coefficients"] = list()
    collection["mixtures"] = list()
    collection["dips"] = list()

    # Iterate on permutations.
    for shuffle_indices in shuffles:
        # Shuffle gene's signals.
        data_shuffle = permute_gene_signals(
            data_gene_signals=data_gene_persons_tissues_signals,
            shuffle_indices=shuffle_indices
        )
        # Distribution
        bins = distribution.prepare_describe_distribution(
            data_gene_persons_tissues_signals=data_shuffle,
        )
        # Collect metrics.
        collection["coefficients"].append(
            bins["bin_modality"]["scores"]["coefficient"]
        )
        collection["dips"].append(
            bins["bin_modality"]["scores"]["dip"]
        )
        collection["mixtures"].append(
            bins["bin_modality"]["scores"]["mixture"]
        )

        pass

    # Convert lists to NumPy arrays.
    collection["coefficients"] = numpy.array(
        collection["coefficients"],
        copy=False,
    )
    collection["mixtures"] = numpy.array(
        collection["mixtures"],
        copy=False,
    )
    collection["dips"] = numpy.array(
        collection["dips"],
        copy=False,
    )

    # Return information.
    return collection


def permute_gene_signals(
    data_gene_signals=None,
    shuffle_indices=None
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
        shuffle_indices (list<list<int>>): matrix of indices

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
    matrix_shuffle = numpy.array(shuffle_indices)
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


# Shuffle


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


# Product


def write_product_gene(gene=None, dock=None, information=None):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files.
        information (object): information to write to file.

    raises:

    returns:

    """

    # Specify directories and files.
    path_permutation = os.path.join(dock, "permutation")
    utility.create_directory(path_permutation)
    path_gene = os.path.join(path_permutation, gene)
    utility.create_directory(path_gene)
    path_permutations = os.path.join(
        path_gene, "permutations.pickle"
    )

    # Write information to file.
    with open(path_permutations, "wb") as file_product:
        pickle.dump(information["permutations"], file_product)

    pass


###############################################################################
# Procedure


def execute_procedure(
    gene=None,
    data_gene_persons_tissues_signals=None,
    shuffles=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues
        shuffles (list<list<list<int>>>): matrices of shuffle indices
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Begin with gene's signals across persons and tissues before restrictive
    # selection of persons by coverage of tissues.
    # Matrices of genes' signals across persons and tissues must have same
    # dimensions.
    # Permutation hence controls for this restrictive selection of persons by
    # coverage of tissues.

    # Shuffle gene's distributions.
    permutations = permute_gene_distribution(
        gene=gene,
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals,
        shuffles=shuffles,
    )

    # Compile information.
    information = {
        "permutations": permutations
    }
    #Write product information to file.
    write_product_gene(
        gene=gene,
        dock=dock,
        information=information
    )

    pass


def execute_procedure_local(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    utility.print_terminal_partition(level=1)
    print("... Permutation procedure ...")

    # Report date and time.
    utility.print_terminal_partition(level=3)
    start = datetime.datetime.now()
    print(start)
    utility.print_terminal_partition(level=3)

    # Remove previous files to avoid version or batch confusion.
    path_permutation = os.path.join(dock, "permutation")
    utility.remove_directory(path=path_permutation)

    # Read source information from file.
    # It is an option to read directory names from "distribution"; however,
    # distribution also writes report files to its directory, and these will
    # include in the list.
    source = read_source_initial(dock=dock)

    # Define genes for iteration.
    #genes_iteration = source["genes_selection"]
    #genes_iteration = source["genes_unimodal"]
    genes_iteration = copy.deepcopy(source["genes_unimodal"])
    genes_iteration.extend(source["genes_multimodal"])
    print("count of genes: " + str(len(genes_iteration)))
    print("count of shuffles: " + str(len(source["shuffles"])))

    #report = execute_procedure_local_sub(
    #    gene="ENSG00000231925", # TAPBP
    #    shuffles=source["shuffles"],
    #    dock=dock,
    #)

    # Set up partial function for iterative execution.
    # Each iteration uses a different sequential value of the "gene" variable
    # with the same value of the "dock" variable.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
        shuffles=source["shuffles"],
        dock=dock,
    )

    # Initialize multiprocessing pool.
    # Limit parallization on halyard to avoid exceeding memory.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=8) # 7 on halyard, 200 on nrnb
    # Iterate on genes.
    report = pool.map(execute_procedure_gene, genes_iteration)

    # Pause procedure.
    time.sleep(10.0)

    # Report.

    # Report date and time.
    utility.print_terminal_partition(level=3)
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))
    utility.print_terminal_partition(level=3)

    pass


def execute_procedure_local_sub(
    gene=None,
    shuffles=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        shuffles (list<list<list<int>>>): matrices of shuffle indices
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_gene_persons_tissues_signals=(
            source["data_gene_persons_tissues_signals"]
        ),
        shuffles=shuffles,
        dock=dock
    )

    # Report progress.
    path_permutation = os.path.join(dock, "permutation")
    directories = os.listdir(path_permutation)
    count = len(directories)
    if (count % 10 == 0):
        print("complete genes: " + str(len(directories)))

    pass


def execute_procedure_remote(dock=None, gene=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Read source information from file.
    source_initial = read_source_initial(
        source_genes="split",
        dock=dock
    )

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_gene_persons_tissues_signals=(
            source["data_gene_persons_tissues_signals"]
        ),
        shuffles=source_initial["shuffles"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure()
