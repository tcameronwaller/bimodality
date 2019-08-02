"""
Module to provide template of common structure of modules.

Title:

    template

Imports:

    os: Package to interact with the operating system.
    sys: Package to interact with the interpreter.
    shutil: Package to perform file operations.
    importlib: Package to import packages and modules.
    argparse: Package to interpret user parameters from terminal.
    csv: Package to organize information in text.
    copy: Package to copy objects.
    pickle: Package to preserve information.
    numpy: Package to calculate with arrays of numbers.
    pandas: Package to organize collections of variables.

Classes:

    This module does not contain any classes.

Exceptions:

    This module does not contain any exceptions.

Functions:

    ...

Author:

    Thomas Cameron Waller
    tcameronwaller@gmail.com
    tcwaller@ucsd.edu
    Department of Medical Genetics
    University of California at San Diego
    Room 3A14, Biomedical Research Facility II
    9500 Gilman Drive
    La Jolla, California 92093
    United States of America

License:

    This file is part of project bimodality
    (https://github.com/tcameronwaller/bimodality/).

    Bimodality ...
    Copyright (C) 2018 Thomas Cameron Waller

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.
    If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import os
import math

# Relevant

import numpy
import pandas

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
    path_assembly = os.path.join(dock, "assembly")
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_gene_count = os.path.join(
        path_assembly, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_assembly, "data_gene_signal.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(
        path_gene_annotation
    )
    data_gene_count = pandas.read_pickle(path_gene_count)
    data_gene_signal = pandas.read_pickle(path_gene_signal)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal,
    }


def check_missing_values(data=None):
    """
    Checks data for missing values and prints reports.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for missing values in genes' signals.")
    print("shape of original data frame: " + str(data.shape))
    print("shape without missing axis 0: " + str(data.dropna(axis=0).shape))
    print("shape without missing axis 1: " + str(data.dropna(axis=1).shape))
    pass


def check_redundancy_genes(data=None):
    """
    Checks data for redundancy in genes.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for redundant genes in genes' signals.")
    print("Consider names of genes.")
    # Reset indices to consider names of genes.
    data = data.reset_index()
    print(data.iloc[0:10, 0:10])
    data_redundancy = data.duplicated(subset=None, keep="first")
    data_redundancy_list = data_redundancy.to_list()
    if any(data_redundancy_list):
        print("Redundancy in genes: Yes")
    else:
        print("Redundancy in genes: No")
    pass


def check_zero_samples(data=None):
    """
    Checks data for samples with values of 0 for all genes' signals.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for samples with values of 0 for all genes' signals.")
    print("shape of original data frame: " + str(data.shape))
    data_nonzero = (data != 0)
    print(
        "shape of data frame without zero samples: " +
        str(data.loc[ : , data_nonzero.any(axis="index")].shape)
    )
    pass


def check_zero_genes(data=None):
    """
    Checks data for genes with values of 0 for signals across all samples.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:

    """

    utility.print_terminal_partition(level=2)
    print("Check for genes with values of 0 for signals across all samples.")
    print("These genes are undetectable.")
    print("shape of original data frame: " + str(data.shape))
    data_nonzero = (data != 0)
    print(
        "shape of data frame without zero genes: " +
        str(data.loc[data_nonzero.any(axis="columns"), : ].shape)
    )
    print("Now printing a summary of data for genes with all zero signals.")
    data_zero = (data == 0)
    data_signal_zero = data.loc[data_zero.all(axis="columns"), : ]
    print(data_signal_zero.iloc[0:10, 0:10])
    #groups = data_signal_zero.groupby(level="gene")
    #print(groups.describe())
    pass


def drop_undetectable_genes(data=None):
    """
    Drops genes with values of 0 for signals across all samples.

    arguments:
        data (object): Pandas data frame of genes' signals for all samples.

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples.

    """

    utility.print_terminal_partition(level=2)
    print("Drop genes that are undetectable.")
    data_nonzero = (data != 0)
    data_signal = data.loc[data_nonzero.any(axis="columns"), : ]
    print("Data without undetectable genes.")
    print(data_signal.iloc[0:10, 0:10])
    print("data dimensions: " + str(data_signal.shape))
    return data_signal




def filter_genes_by_signal_threshold(
    data=None,
    threshold=None,
):
    """
    Filter genes to keep only those with signals beyond threshold in at least
    one sample.

    Data format should have samples across columns and genes across rows.

    arguments:
        data (object): Pandas data frame of genes' signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    utility.print_terminal_partition(level=2)
    print(
        "Filter genes to keep only those with signals beyond threshold in " +
        "at least one sample. \n" +
        "Data format should have samples across columns and genes across rows."
    )
    print("signal threshold: " + str(threshold))
    print("data dimensions before filter: " + str(data.shape))
    data_threshold = (data >= threshold)
    data_detection = data.loc[data_threshold.any(axis="columns"), : ]
    print("data dimensions after filter: " + str(data_detection.shape))
    utility.print_terminal_partition(level=3)
    return data_detection


def filter_samples_by_signal_threshold(
    data=None,
    threshold=None,
):
    """
    Filter samples to keep only those with signals beyond threshold in at least
    one gene.

    Data format should have samples across columns and genes across rows.

    arguments:
        data (object): Pandas data frame of genes' signals across samples

    raises:

    returns:
        (object): Pandas data frame of genes' signals across samples

    """

    utility.print_terminal_partition(level=2)
    print(
        "Filter samples to keep only those with signals beyond threshold in " +
        "at least one gene. \n" +
        "Data format should have samples across columns and genes across rows."
    )
    print("signal threshold: " + str(threshold))
    print("data dimensions before filter: " + str(data.shape))
    data_threshold = (data >= threshold)
    data_detection = data.loc[:, data_threshold.any(axis="index")]
    print("data dimensions after filter: " + str(data_detection.shape))
    utility.print_terminal_partition(level=3)
    return data_detection







def transform_gene_signal_log(
    data_gene_signal=None,
    pseudo_count=None,
):
    """
    Transforms values of genes' signals to base-two logarithmic space.

    arguments:
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.
        pseudo_count (float): value to add to signal before transformation

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values. <- maybe false
    # Logarithms are not distributive.
    # To accommodate gene signals of 0.0, add a pseudo count of 2.0 to all
    # counts before calculation of the base-2 logarithm.
    # log-2 signal = log2(TPM + pseudo-count)
    # pseudo-count = 1.0
    if False:
        utility.print_terminal_partition(level=2)
        print("Transformation of genes' signals to base-2 logarithmic space.")
        print(
            "To accommodate gene signals of 0.0, add a pseudo count of " +
            "1.0-5.0 to all counts before calculation of the base-2 logarithm."
        )
    data_gene_signal_log = calculate_logarithm_gene_signal(
        pseudo_count=pseudo_count,
        data_gene_signal=data_gene_signal
    )
    #print(data_gene_signal_log.iloc[0:10, 0:10])

    return data_gene_signal_log


def calculate_logarithm_gene_signal(pseudo_count=None, data_gene_signal=None):
    """
    Calculates the base-2 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base-2 logarithm.

    arguments:
        pseudo_count (float): Pseudo count to add to gene signal before
            transformation to avoid values of zero.
        data_gene_signal (object): Pandas data frame of signals for all genes
            across specific persons and tissues.

    raises:

    returns:
        (object): Pandas data frame of base-2 logarithmic signals for all genes
            across specific persons and tissues.

    """

    # lambda x: math.log((x + 1), 2)

    # An alternative approach would be to set the label columns to indices.
    if False:
        data_gene_signal_index = data_gene_signal.set_index(
            ["Name", "Description"], append=True, drop=True
        )
        data_gene_signal_log = data_gene_signal_index.apply(
            lambda value: 2 * value
        )
        data_log = data_signal_index.copy()
        data_log.iloc[0:, 2:] = data_log.iloc[0:, 2:].applymap(
            lambda value: math.log((value + 1.0), 2)
        )

    data_log = data_gene_signal.applymap(
        lambda value: math.log((value + pseudo_count), 2)
    )
    # Reverse an index to a column.
    return data_log


def select_genes_counts(
    data_gene_count=None,
    data_gene_annotation=None,
):
    """
    Selects genes' counts.

    arguments:
        data_gene_count (object): Pandas data frame of genes' counts across
            samples.
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' counts for all samples

    """

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Selection of genes' counts.")

    print("count of samples, original: " + str(data_gene_count.shape[1]))

    # Select genes that encode proteins.
    data_gene_count = select_genes_protein(
        data_gene_annotation=data_gene_annotation,
        data_gene_signal=data_gene_count
    )

    # Check for data quality.
    utility.print_terminal_partition(level=2)
    print("Check for quality of genes' counts.")

    # Check for missing values of genes' signals.
    check_missing_values(data=data_gene_count)

    # Check for redundant genes.
    check_redundancy_genes(data=data_gene_count)

    # Check for samples with values of 0 for all genes' signals.
    check_zero_samples(data=data_gene_count)

    # Check for genes with values of 0 for signals across all samples.
    check_zero_genes(data=data_gene_count)

    # Remove irrelevant signals for genes.
    utility.print_terminal_partition(level=2)
    print("Removal of signals for undetectable genes.")
    # Drop undetectable genes.
    #data_gene_count = drop_undetectable_genes(data=data_gene_count)

    # Filter genes by signal.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_gene_count = filter_genes_by_signal_threshold(
        data=data_gene_count,
        threshold=1,
    )

    # Filter samples by signal.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_gene_count = filter_samples_by_signal_threshold(
        data=data_gene_count,
        threshold=1,
    )

    utility.print_terminal_partition(level=2)

    print(data_gene_count.iloc[0:10, 0:7])

    utility.print_terminal_partition(level=2)

    print("Count of original genes: 56202")
    print("Count of protein-coding genes: 18842")
    print(
        "Count of detectable genes, count beyond 0.0 in at least 1 sample: " +
        "18813"
    )
    print(
        "Count of detectable genes, count beyond 1 in at least 1 sample: " +
        "18813"
    )
    print("count of samples, final: " + str(data_gene_count.shape[1]))

    # Return information.
    return data_gene_count


def select_genes_signals(
    data_gene_signal=None,
    data_gene_annotation=None,
):
    """
    Selects genes' signals.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals across
            samples
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' signals for all samples, tissues,
            and patients

    """

    # Organize genes' signals.
    utility.print_terminal_partition(level=1)
    print("Selection of genes' signals.")

    print("count of samples, original: " + str(data_gene_signal.shape[1]))

    # Select genes that encode proteins.
    data_gene_signal = select_genes_protein(
        data_gene_annotation=data_gene_annotation,
        data_gene_signal=data_gene_signal
    )

    # Check for data quality.
    utility.print_terminal_partition(level=2)
    print("Check for quality of genes' signals.")

    # Check for missing values of genes' signals.
    check_missing_values(data=data_gene_signal)

    # Check for redundant genes.
    check_redundancy_genes(data=data_gene_signal)

    # Check for samples with values of 0 for all genes' signals.
    check_zero_samples(data=data_gene_signal)

    # Check for genes with values of 0 for signals across all samples.
    check_zero_genes(data=data_gene_signal)

    # Remove irrelevant signals for genes.
    utility.print_terminal_partition(level=2)
    print("Removal of signals for undetectable genes.")
    # Drop undetectable genes.
    #data_gene_signal = drop_undetectable_genes(data=data_gene_signal)

    # Filter genes by signal.
    # Filter to keep only genes with signals beyond threshold in at least one
    # sample.
    data_gene_signal = filter_genes_by_signal_threshold(
        data=data_gene_signal,
        threshold=1.0,
    )

    # Filter samples by signal.
    # Filter to keep only samples with signals beyond threshold in at least one
    # gene.
    data_gene_signal = filter_samples_by_signal_threshold(
        data=data_gene_signal,
        threshold=1.0,
    )


    utility.print_terminal_partition(level=2)

    print(data_gene_signal.iloc[0:10, 0:7])

    utility.print_terminal_partition(level=2)

    print("Count of original genes: 56202")
    print("Count of protein-coding genes: 18842")
    print(
        "Count of detectable genes, signal beyond 0.0 in at least 1 sample: " +
        "18813"
    )
    print(
        "Count of detectable genes, signal beyond 1.0 in at least 1 sample: " +
        "18511"
    )
    print("count of samples, final: " + str(data_gene_signal.shape[1]))

    # Return information.
    return data_gene_signal


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
    path_measurement = os.path.join(dock, "measurement")
    utility.confirm_path_directory(path_measurement)
    path_gene_count = os.path.join(
        path_measurement, "data_gene_count.pickle"
    )
    path_gene_signal = os.path.join(
        path_measurement, "data_gene_signal.pickle"
    )

    # Write information to file.
    pandas.to_pickle(
        information["data_gene_count"],
        path_gene_count
    )
    pandas.to_pickle(
        information["data_gene_signal"],
        path_gene_signal
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

    # Organize genes' counts.
    data_gene_count = select_genes_counts(
        data_gene_count=source["data_gene_count"],
        data_gene_annotation=source["data_gene_annotation"],
    )

    # Organize genes' signals.
    data_gene_signal = select_genes_signals(
        data_gene_signal=source["data_gene_signal"],
        data_gene_annotation=source["data_gene_annotation"],
    )

    # Compile information.
    information = {
        "data_gene_count": data_gene_count,
        "data_gene_signal": data_gene_signal
    }

    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
