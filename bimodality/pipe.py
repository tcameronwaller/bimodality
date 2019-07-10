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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import shuffle
import organization
import restriction
import distribution
import metric
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source_local_initial(dock=None):
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
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    # Compile and return information.
    return {
        "genes": genes,
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
    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
        "shuffles": shuffles,
    }


##########
# Distribution


def describe_gene_signal_distribution(
    method=None,
    count=None,
    tissues=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        tissues (list<str>): specific tissues to select by "imputation" method
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Restriction
    collection_restriction = restriction.execute_procedure(
        method=method,
        count=count,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            data_gene_persons_tissues_signals
        ),
    )

    # Distribution
    collection_distribution = distribution.execute_procedure(
        data_gene_persons_tissues_signals=(
            collection_restriction["data_gene_persons_tissues_signals"]
        ),
    )

    # Compile information.
    information = {
        "report_restriction": collection_restriction["report_gene"],
        "values": collection_distribution["values"],
        "scores": collection_distribution["scores"],
        "report_distribution": collection_distribution["report_gene"],
    }

    # Return information.
    return information


##########
# Shuffle


def shuffle_gene_signal_distribution(
    method=None,
    count=None,
    tissues=None,
    shuffles=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distribution of a gene's signals across persons.

    arguments:
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        tissues (list<str>): specific tissues to select by "imputation" method
        shuffles (list<list<list<int>>>): matrices of indices
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Initialize collections of metrics.
    coefficients = list()
    dips = list()
    mixtures = list()
    # Iterate on shuffles.
    for shuffle_matrix in shuffles:
        # Shuffle gene's signals.
        data_shuffle = shuffle.shuffle_gene_signals(
            data_gene_signals=data_gene_persons_tissues_signals,
            shuffle=shuffle_matrix
        )
        # Prepare and describe distribution of shuffle of gene's signals.
        distribution_shuffle = describe_gene_signal_distribution(
            method=method,
            count=count,
            tissues=tissues,
            data_gene_persons_tissues_signals=data_shuffle
        )
        # Collect metrics.
        coefficients.append(
            distribution_shuffle["scores"]["coefficient"]
        )
        dips.append(
            distribution_shuffle["scores"]["dip"]
        )
        mixtures.append(
            distribution_shuffle["scores"]["mixture"]
        )
        pass
    # Compile information.
    information = {
        "coefficient": coefficients,
        "dip": dips,
        "mixture": mixtures,
    }

    # Return information.
    return information


##########
# Product


def write_product(gene=None, dock=None, information=None):
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
    path_pipe = os.path.join(dock, "pipe")
    utility.confirm_path_directory(path_pipe)
    path_gene = os.path.join(path_pipe, gene)
    utility.confirm_path_directory(path_gene)
    path_report_organization = os.path.join(
        path_gene, "report_organization.pickle"
    )
    path_report_restriction = os.path.join(
        path_gene, "report_restriction.pickle"
    )
    path_report_distribution = os.path.join(
        path_gene, "report_distribution.pickle"
    )
    path_scores_real = os.path.join(
        path_gene, "scores_real.pickle"
    )
    path_scores_shuffle = os.path.join(
        path_gene, "scores_shuffle.pickle"
    )
    # Write information to file.
    with open(path_report_organization, "wb") as file_product:
        pickle.dump(information["report_organization"], file_product)
    with open(path_report_restriction, "wb") as file_product:
        pickle.dump(information["report_restriction"], file_product)
    with open(path_report_distribution, "wb") as file_product:
        pickle.dump(information["report_distribution"], file_product)
    with open(path_scores_real, "wb") as file_product:
        pickle.dump(information["scores_real"], file_product)
    with open(path_scores_shuffle, "wb") as file_product:
        pickle.dump(information["scores_shuffle"], file_product)

    pass


###############################################################################
# Procedure


def execute_procedure(
    gene=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    # Report gene.

    print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Collect garbage to clear memory.
    gc.collect()

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=source["data_gene_samples_signals"]
    )

    # Prepare and describe distribution of real gene's signals.
    # Method for selection is either "availability" or "imputation".
    # Specific tissues are only relevant to "imputation" method.
    tissues = [
        "adipose", "blood", "colon", "esophagus", "heart", "muscle",
        "lung", "nerve", "skin", "thyroid",
    ]
    distribution_real = describe_gene_signal_distribution(
        method="availability",
        count=7,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Prepare and describe distributions of shuffles of gene's signals.
    scores_shuffle = shuffle_gene_signal_distribution(
        method="availability",
        count=7,
        tissues=tissues,
        shuffles=source["shuffles"],
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Compile information.
    information = {
        "report_organization": collection_organization["report_gene"],
        "report_restriction": distribution_real["report_restriction"],
        "report_distribution": distribution_real["report_distribution"],
        "scores_real": distribution_real["scores"],
        "scores_shuffle": scores_shuffle,
    }
    #Write product information to file.
    write_product(gene=gene, dock=dock, information=information)

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

    # Read source information from file.
    source = read_source_local_initial(dock=dock)

    print("count of genes: " + str(len(source["genes"])))
    #print("count of shuffles: " + str(len(source["shuffles"])))

    # Report date and time.
    print(datetime.datetime.now())

    # Remove previous files to avoid version or batch confusion.
    path_pipe = os.path.join(dock, "pipe")
    utility.remove_directory(path=path_pipe)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "shuffles", and
    # "dock" variables.
    execute_procedure_gene = functools.partial(
        execute_procedure,
        dock=dock
    )

    # Initialize multiprocessing pool.
    # Iterative shuffle procedure.
    # 4 concurrent processes require approximately 60% of 32 Gigabytes memory.
    # 5 concurrent processes require approximately 75% of 32 Gigabytes memory.
    # 6 concurrent processes require approximately 90% of 32 Gigabytes memory.
    # Index shuffle procedure.
    # 6 concurrent processes require approximately 50% of 32 Gigabytes memory.
    # 7 concurrent processes require approximately 65% of 32 Gigabytes memory.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=7)

    # Iterate on genes.
    report = pool.map(execute_procedure_gene, source["genes"][0:14])
    #report = pool.map(execute_procedure_gene, source["genes"])

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # Report date and time.
    print(datetime.datetime.now())

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

    # Execute procedure.
    execute_procedure(
        gene=gene,
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
