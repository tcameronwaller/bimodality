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
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

# Custom

import measurement
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
    path_demo = os.path.join(
        path_split, "demonstration_gene_samples_signals.pickle"
    )
    # Read information from file.
    data_gene_samples_signals = pandas.read_pickle(path_demo)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Aggregation of genes' signals by person and major tissue category


def aggregate_gene_signal_tissue(
    data_gene_samples_signals=None
):
    """
    Aggregates gene's signals within groups by person and major tissue
    category.
    Aggregates gene's signals across any minor tissue categories or redundant
    samples.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Split data by person and major tissue category.
    if False:
        data_gene_samples_signals.reset_index(
            level=["sample", "tissue_minor"],
            inplace=True
        )
        data_gene_samples_signals.drop(
            labels=["sample", "tissue_minor"],
            axis="columns",
            inplace=True
        )
    groups = data_gene_samples_signals.groupby(
        level=["person", "tissue_major"]
    )
    # Aggregate genes' signals within groups by mean.
    data_gene_signal_tissue = groups.aggregate(statistics.mean)
    # Return information.
    return data_gene_signal_tissue


##########
# Organization of gene's signals across matrices of persons by tissues.


def organize_gene_signal_persons_tissues(
    data_gene_signal_tissue=None
):
    """
    Organizes gene's signals across matrices of persons by tissues.

    arguments:
        data_gene_signal_tissue (object): Pandas data frame of a gene's signals
            across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    series = pandas.Series(
        data_gene_signal_tissue["signal"],
        index=data_gene_signal_tissue.index
    )
    data_gene_persons_tissues_signals = series.unstack(
        level="tissue_major"
    )
    # Return information.
    return data_gene_persons_tissues_signals


##########
# Report.


def calculate_gene_tissue_mean(
    data_gene_samples_signals=None
):
    """
    Calculates the mean of a gene's signals across samples within major
    tissue categories.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Calculate mean of a gene's signals across samples for each major
    # tissue category.
    # Split data by major tissue category.
    groups = data_gene_samples_signals.groupby(
        level=["tissue_major"]
    )
    # Aggregate genes' signals within groups by mean.
    data_mean = groups.aggregate(
        lambda x: x.mean()
    )
    return data_mean


def calculate_gene_tissue_variance(
    data_gene_samples_signals=None
):
    """
    Calculates the variance of a gene's signals across samples within major
    tissue categories.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Calculate variance of a gene's signals across samples for each major
    # tissue category.
    # Split data by major tissue category.
    groups = data_gene_samples_signals.groupby(
        level=["tissue_major"]
    )
    # Aggregate genes' signals within groups by mean.
    data_variance = groups.aggregate(
        lambda x: x.var()
    )
    return data_variance


def calculate_gene_tissue_deviation(
    data_gene_samples_signals=None
):
    """
    Calculates the standard deviation of a gene's signals across samples within
    major tissue categories.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Calculate standard deviation of a gene's signals across samples for each
    # major tissue category.
    # Split data by major tissue category.
    groups = data_gene_samples_signals.groupby(
        level=["tissue_major"]
    )
    # Aggregate genes' signals within groups by mean.
    data_deviation = groups.aggregate(
        lambda x: x.std()
    )
    return data_deviation


def calculate_gene_tissue_variation(
    data_gene_samples_signals=None
):
    """
    Calculates the coefficient of variation of a gene's signals across samples
    within major tissue categories.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Calculate coefficient of variation of a gene's signals across samples for
    # each major tissue category.
    # Split data by major tissue category.
    groups = data_gene_samples_signals.groupby(
        level=["tissue_major"]
    )
    # Aggregate genes' signals within groups by relative variance.
    # Some genes' signals in some groups have means of zero.
    data_variation = groups.aggregate(
        lambda x: (x.std() / x.mean()) if (x.mean() != 0) else (numpy.nan)
    )
    return data_variation


def calculate_gene_tissue_dispersion(
    data_gene_samples_signals=None
):
    """
    Calculates the dispersion of a gene's signals across samples within major
    tissue categories.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    # Calculate dispersion of a gene's signals across samples for each major
    # tissue category.
    # Split data by major tissue category.
    groups = data_gene_samples_signals.groupby(
        level=["tissue_major"]
    )
    # Aggregate genes' signals within groups by relative variance.
    # Some genes' signals in some groups have means of zero.
    data_dispersion = groups.aggregate(
        lambda x: (x.var() / x.mean()) if (x.mean() != 0) else (numpy.nan)
    )
    return data_dispersion


def prepare_report_gene(
    data_gene_samples_signals=None,
):
    """
    Prepare report for gene.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): report for gene

    """

    # Calculate mean of gene's signals across major tissue categories.
    data_gene_tissue_mean = calculate_gene_tissue_mean(
        data_gene_samples_signals=data_gene_samples_signals
    )
    # Calculate variance of gene's signals across major tissue categories.
    data_gene_tissue_variance = calculate_gene_tissue_variance(
        data_gene_samples_signals=data_gene_samples_signals
    )
    # Calculate standard deviation of gene's signals across major tissue
    # categories.
    data_gene_tissue_deviation = calculate_gene_tissue_deviation(
        data_gene_samples_signals=data_gene_samples_signals
    )
    # Calculate coefficient of variation of gene's signals across major tissue
    # categories.
    data_gene_tissue_variation = calculate_gene_tissue_variation(
        data_gene_samples_signals=data_gene_samples_signals
    )
    # Calculate dispersion in gene's signals across major tissue categories.
    data_gene_tissue_dispersion = calculate_gene_tissue_dispersion(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Compile information.
    information = {
        "data_gene_tissue_mean": data_gene_tissue_mean,
        "data_gene_tissue_variance": data_gene_tissue_variance,
        "data_gene_tissue_deviation": data_gene_tissue_deviation,
        "data_gene_tissue_variation": data_gene_tissue_variation,
        "data_gene_tissue_dispersion": data_gene_tissue_dispersion,
    }
    # Return information.
    return information


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
    path_organization = os.path.join(dock, "organization")
    utility.confirm_path_directory(path_organization)
    path_gene_signal = os.path.join(
        path_organization, "data_gene_persons_tissues_signals.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_gene_persons_tissues_signals"], path_gene_signal
    )

    pass


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
    print("Welcome to the test, demonstration of the organization procedure!")

    # Read source information from file.
    source = read_source(dock=dock)
    print(source["data_gene_samples_signals"])

    # Call procedure.
    information = execute_procedure(
        data_gene_samples_signals=source["data_gene_samples_signals"]
    )

    print("here is the gene dispersion report...")
    print(information["report_gene"]["data_gene_tissue_dispersion"])

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


def execute_procedure(
    data_gene_samples_signals=None
):
    """
    Function to execute module's main behavior.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): report about gene, Pandas data frame of gene's signals across
            persons and tissues

    """

    # Aggregate gene's signals within groups by person and major tissue
    # category.
    # Aggregate across any minor tissue categories or redundant samples.
    data_gene_signal_tissue = aggregate_gene_signal_tissue(
        data_gene_samples_signals=data_gene_samples_signals
    )
    #print(data_gene_signal_tissue)

    # Organize gene's signals across matrices of persons by tissues.
    data_gene_persons_tissues_signals = organize_gene_signal_persons_tissues(
        data_gene_signal_tissue=data_gene_signal_tissue
    )
    #print(data_gene_persons_tissues_signals)

    # Prepare gene report.
    # Describe dispersion of gene's signals across major tissue categories.
    report_gene = prepare_report_gene(
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Compile information.
    information = {
        "report_gene": report_gene,
        "data_gene_persons_tissues_signals": data_gene_persons_tissues_signals,
    }
    # Return information.
    return information


if (__name__ == "__main__"):
    execute_procedure()
