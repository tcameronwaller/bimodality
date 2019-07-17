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


def read_source_local_initial(
    source_genes=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        source_genes (str): name of directory from which to obtain genes list
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    if source_genes == "split":
        path_source = os.path.join(dock, "split")
    elif source_genes == "candidacy":
        path_source = os.path.join(dock, "candidacy")
    elif source_genes == "combination":
        path_source = os.path.join(dock, "combination")
    path_genes = os.path.join(
        path_source, "genes.txt"
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
    path_imputation = os.path.join(path_gene, "imputation")
    utility.confirm_path_directory(path_imputation)
    path_availability = os.path.join(path_gene, "availability")
    utility.confirm_path_directory(path_availability)

    path_report_organization = os.path.join(
        path_gene, "report_organization.pickle"
    )
    # Imputation method
    path_imputation_report_restriction = os.path.join(
        path_imputation, "report_restriction.pickle"
    )
    path_imputation_report_distribution = os.path.join(
        path_imputation, "report_distribution.pickle"
    )
    path_imputation_scores = os.path.join(
        path_imputation, "scores.pickle"
    )
    path_imputation_shuffles = os.path.join(
        path_imputation, "shuffles.pickle"
    )
    # Availability method
    path_availability_report_restriction = os.path.join(
        path_availability, "report_restriction.pickle"
    )
    path_availability_report_distribution = os.path.join(
        path_availability, "report_distribution.pickle"
    )
    path_availability_scores = os.path.join(
        path_availability, "scores.pickle"
    )
    path_availability_shuffles = os.path.join(
        path_availability, "shuffles.pickle"
    )

    # Write information to file.
    with open(path_report_organization, "wb") as file_product:
        pickle.dump(information["report_organization"], file_product)
    # Imputation method
    with open(path_imputation_report_restriction, "wb") as file_product:
        pickle.dump(information["report_restriction_imputation"], file_product)
    with open(path_imputation_report_distribution, "wb") as file_product:
        pickle.dump(
            information["report_distribution_imputation"], file_product
        )
    with open(path_imputation_scores, "wb") as file_product:
        pickle.dump(information["scores_imputation"], file_product)
    with open(path_imputation_shuffles, "wb") as file_product:
        pickle.dump(information["shuffles_imputation"], file_product)

    # Imputation method
    with open(path_availability_report_restriction, "wb") as file_product:
        pickle.dump(
            information["report_restriction_availability"], file_product
        )
    with open(path_availability_report_distribution, "wb") as file_product:
        pickle.dump(
            information["report_distribution_availability"], file_product
        )
    with open(path_availability_scores, "wb") as file_product:
        pickle.dump(information["scores_availability"], file_product)
    with open(path_availability_shuffles, "wb") as file_product:
        pickle.dump(information["shuffles_availability"], file_product)

    pass


###############################################################################
# Procedure

# Change procedure to run and collect both by "availability" and "imputation" methods...


def execute_procedure(
    gene=None,
    shuffles=None,
    data_gene_samples_signals=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        shuffles (list<list<list<int>>>): matrices of indices
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Prepare and describe distribution of real gene's signals.
    # Method for selection is either "availability" or "imputation".
    # Specific tissues are only relevant to "imputation" method.
    tissues = [
        "adipose", # 552
        #"adrenal", # 190
        "artery", # 551
        "blood", # 407
        "brain", # 254
        "colon", # 371
        "esophagus", # 513
        "heart", # 399
        #"liver", # 175
        "lung", # 427
        "muscle", # 564
        "nerve", # 414
        "pancreas", # 248
        #"pituitary", # 183
        "skin", # 583
        #"intestine", # 137
        #"salivary", # 97 # now excluded at "selection" procedure... threshold 100
        #"spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ] # 300
    distribution_imputation = describe_gene_signal_distribution(
        method="imputation",
        count=10,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )
    distribution_availability = describe_gene_signal_distribution(
        method="availability",
        count=10,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Prepare and describe distributions of shuffles of gene's signals.
    shuffles_imputation = shuffle_gene_signal_distribution(
        method="imputation",
        count=10,
        tissues=tissues,
        shuffles=shuffles,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )
    shuffles_availability = shuffle_gene_signal_distribution(
        method="availability",
        count=10,
        tissues=tissues,
        shuffles=shuffles,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Compile information.
    information = {
        "report_organization": collection_organization["report_gene"],
        "report_restriction_imputation": (
            distribution_imputation["report_restriction"]
        ),
        "report_restriction_availability": (
            distribution_availability["report_restriction"]
        ),
        "report_distribution_imputation": (
            distribution_imputation["report_distribution"]
        ),
        "report_distribution_availability": (
            distribution_availability["report_distribution"]
        ),
        "scores_imputation": distribution_imputation["scores"],
        "scores_availability": distribution_availability["scores"],
        "shuffles_imputation": shuffles_imputation,
        "shuffles_availability": shuffles_availability,
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
    source = read_source_local_initial(source_genes="candidacy", dock=dock)

    print("count of genes: " + str(len(source["genes"])))
    #print("count of shuffles: " + str(len(source["shuffles"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    # Remove previous files to avoid version or batch confusion.
    path_pipe = os.path.join(dock, "pipe")
    utility.remove_directory(path=path_pipe)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "shuffles", and
    # "dock" variables.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
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
    pool = multiprocessing.Pool(processes=8)

    # Iterate on genes.
    check_genes=[
        "ENSG00000198965",
    ]
    #report = pool.map(execute_procedure_gene, check_genes)
    report = pool.map(execute_procedure_gene, source["genes"][0:8])
    #report = pool.map(execute_procedure_gene, source["genes"])

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # Report date and time.
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))

    pass


def execute_procedure_local_sub(gene=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Enable automatic garbage collection to clear memory.
    #gc.enable()

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Collect garbage to clear memory.
    #gc.collect()

    # Execute procedure.
    execute_procedure(
        gene=gene,
        shuffles=source["shuffles"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    # Report contents of directory.
    path_pipe = os.path.join(dock, "pipe")
    directories = os.listdir(path_pipe)
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

    # Execute procedure.
    execute_procedure_remote_sub(
        gene=gene,
        dock=dock
    )

    pass


def execute_procedure_remote_sub(gene=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Enable automatic garbage collection to clear memory.
    gc.enable()

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Collect garbage to clear memory.
    gc.collect()

    # Execute procedure.
    execute_procedure(
        gene=gene,
        shuffles=source["shuffles"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
