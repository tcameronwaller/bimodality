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

# Relevant

import numpy
import pandas

# Custom

import organization
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
    path_measurement = os.path.join(dock, "measurement")
    path_gene_signal = os.path.join(
        path_measurement, "data_gene_signal.pickle"
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


def select_samples(
    tissues_major=None,
    data_samples_tissues_persons=None
):
    """
    Selects samples for tissues and persons of interest.

    arguments:
        tissues_major (list<str>): major tissues of interest
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.

    raises:

    returns:
        (list<str>): samples of interest

    """

    # Select samples from persons and tissues of interest.
    utility.print_terminal_partition(level=2)
    print("Selection of samples from persons and tissues of interest.")
    print(
        "count of samples, original: " +
        str(data_samples_tissues_persons.shape[0])
    )
    data_samples = data_samples_tissues_persons.loc[
        data_samples_tissues_persons["tissue_major"].isin(tissues_major),
    ]
    print(
        "count of samples from tissues of interest: " +
        str(data_samples.shape[0])
    )
    print(data_samples)
    # Extract identifiers of samples.
    samples = data_samples.index.to_list()
    print(len(samples))
    return samples


def select_genes_detection(data_gene_signal=None):
    """
    Selects detectable genes with nonzero signals.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons.

    """

    utility.print_terminal_partition(level=2)
    print(
        "Selection of detectable genes with nonzero signals in persons " +
        "and tissues of interest."
    )
    print("genes, original: " + str(data_gene_signal.shape[1]))
    data_nonzero = (data_gene_signal != 0)
    data_signal = data_gene_signal.loc[ : , data_nonzero.any(axis="index")]
    print("genes, detection: " + str(data_signal.shape[1]))
    return data_signal


def select_samples_genes(
    samples=None,
    data_gene_signal=None
):
    """
    Selects samples and genes of interest for further analyses.

    arguments:
        samples (list<str>): samples of interest
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons

    """

    # Select samples and genes of interest.
    utility.print_terminal_partition(level=1)
    print("Selection of samples and genes of interest.")
    print(
        "count of samples, original: " +
        str(data_gene_signal.shape[1])
    )
    data_samples = data_gene_signal.loc[
        :,
        samples
    ]
    print(
        "count of samples of interest: " +
        str(data_samples.shape[1])
    )

    # Select genes with detectable, non-zero signal in tissues and persons of
    # interest.
    data_samples = select_genes_detection(
        data_gene_signal=data_samples
    )

    # Return information.
    return data_samples


def select_samples_genes_by_tissues(
    data_samples_tissues_persons=None,
    data_gene_signal=None
):
    """
    Selects samples and genes by tissues of interest.

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_signal (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons

    """

    # Define tissues of interest.
    tissues_major = [
        "adipose",
        "adrenal",
        "artery",
        "blood",
        "brain",
        "colon",
        "esophagus",
        "heart",
        "liver",
        "lung",
        "muscle",
        "nerve",
        "pancreas",
        "pituitary",
        "skin",
        "intestine",
        "spleen",
        "stomach",
        "thyroid",
    ]
    # Select samples and genes of interest.
    samples = select_samples(
        tissues_major=tissues_major,
        data_samples_tissues_persons=data_samples_tissues_persons
    )
    # Selection of samples and genes of interest.
    data_gene_tissue = select_samples_genes(
        samples=samples,
        data_gene_signal=data_gene_signal
    )
    # Return information.
    return data_gene_tissue


def select_samples_genes_by_few_tissues(
    tissues_major=None,
    data_gene_sample=None
):
    """
    Selects samples and genes by tissues of interest.

    arguments:
        tissues_major (list<str>): major tissues of interest
        data_gene_sample (object): Pandas data frame of genes' signals for all
            samples, tissues, and persons

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and persons

    """

    data_gene_sample.reset_index(
        level=["tissue_major", "tissue_minor"],
        inplace=True
    )
    data_selection = data_gene_sample.loc[
        data_gene_sample["tissue_major"].isin(tissues_major),
    ]
    data_selection.set_index(
        ["tissue_major", "tissue_minor"],
        append=False,
        drop=True,
        inplace=True
    )
    return data_selection





##########
# Product.


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
    path_selection = os.path.join(dock, "selection")
    utility.confirm_path_directory(path_selection)
    path_gene_signal = os.path.join(
        path_selection, "data_gene_signal.pickle"
    )
    path_gene_sample = os.path.join(
        path_selection, "data_gene_sample.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
    )
    pandas.to_pickle(
        information["data_gene_sample"],
        path_gene_sample
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

    # Select samples and genes by relevant major tissues.
    data_gene_signal = select_samples_genes_by_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=source["data_gene_signal"],
    )
    print(data_gene_signal)
    # Associate genes to persons and tissues for each sample.
    data_gene_sample = organization.associate_genes_samples_persons_tissues(
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_signal=data_gene_signal,
    )
    print(data_gene_sample)
    # Select data for individual tissues.
    if False:
        data_gene_selection = select_samples_genes_by_few_tissues(
            tissues_major=["brain"],
            data_gene_sample=data_gene_sample,
        )
    # Compile information.
    information = {
        "data_gene_signal": data_gene_signal,
        "data_gene_sample": data_gene_sample,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
