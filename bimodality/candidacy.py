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

import pipe
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Source


def read_source_initial(
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
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Candidacy


def evaluate_gene_candidacy(
    gene=None,
    method=None,
    count=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (bool): whether the gene is a candidate

    """

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Prepare and describe distribution of real gene's signals.
    # Method for selection is either "availability" or "imputation".
    # Specific tissues are only relevant to "imputation" method.
    tissues = [
        "adipose", "blood", "colon", "esophagus", "heart", "muscle",
        "lung", "nerve", "skin", "thyroid",
    ]
    distribution = pipe.describe_gene_signal_distribution(
        method=method,
        count=count,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        )
    )

    # Access information.
    data_persons_tissues = (
        distribution["report_restriction"]["data_persons_tissues"]
    )
    data_gene_persons_signals = (
        distribution["report_distribution"]["data_gene_persons_signals"]
    )

    # Determine whether gene's aggregate scores across tissues have adequate
    # variance to apply metrics to the distribution.
    # TODO: > 10 non-missing and non-zero values...
    # TODO: standard deviation > 0.0001?


    pass


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
    source_initial = read_source_initial(
        source_genes="split",
        dock=dock
    )
    print("count of genes: " + str(len(source_initial["genes"])))

    # Intitialize collections of candidate genes.
    genes_candidacy_imputation = list()
    genes_candidacy_availability = list()

    # Initialize counter.
    counter = 0
    for gene in source["genes"]:
        # Read source information from file.
        source = read_source(gene=gene, dock=dock)

        # Evaluate gene's candidacy.
        # Imputation
        candidacy_imputation = evaluate_gene_candidacy(
            gene=gene,
            method="imputation",
            count=7,
            data_samples_tissues_persons=(
                source["data_samples_tissues_persons"]
            ),
            data_gene_samples_signals=source["data_gene_samples_signals"],
        )
        if candidacy_imputation:
            genes_candidacy_imputation.append(gene)
        # Availability
        candidacy_availability = evaluate_gene_candidacy(
            gene=gene,
            method="availability",
            count=7,
            data_samples_tissues_persons=(
                source["data_samples_tissues_persons"]
            ),
            data_gene_samples_signals=source["data_gene_samples_signals"],
        )
        if candidacy_availability:
            genes_candidacy_availability.append(gene)


        # Increment counter.
        counter += 1
        # Report progress.
        if (counter % 10 == 0):
            print("complete genes: " + str(counter))
        pass

    # Combine lists of genes... take the unique union...

    # Compile information.
    information = {
        "genes": genes_candidacy,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
