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
import rpy2
import rpy2.robjects
import rpy2.robjects.packages
utils = rpy2.robjects.packages.importr("utils")
#base = rpy2.robjects.packages.importr("base")
utils.install_packages("diptest")
diptest = rpy2.robjects.packages.importr("diptest")

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
    path_access = os.path.join(dock, "access")
    path_organization = os.path.join(dock, "organization")
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    # Read information from file.
    data_gene_signal_imputation = pandas.read_pickle(path_imputation)
    # Compile and return information.
    return {
        "data_gene_signal_imputation": data_gene_signal_imputation,
    }


def calculate_bimodality_coefficient(series=None):
    """
    Calculates the bimodality coefficient of a series.

    arguments:
        series (list<float>): series of values of type float

    raises:

    returns:
        (float): value of bimodality coefficient

    """

    # Calculate skewness.
    skewness = scipy.stats.skew(series, axis=None)
    # Calculate excess kurtosis.
    kurtosis = scipy.stats.kurtosis(series, axis=None, fisher=True)
    # Calculate count.
    count = len(series)
    # Calculate count factor.
    count_factor = ((count - 1) ** 2) / ((count - 2) * (count - 3))
    # Count bimodality coefficient.
    coefficient = ((skewness ** 2) + 1) / (kurtosis + (3 * count_factor))
    # Return value.
    return coefficient


def calculate_dip_statistic(series=None):
    """
    Calculates the Hartigans' dip statistic of a series.

    arguments:
        series (list<float>): series of values of type float

    raises:

    returns:
        (float): value of dip statistic

    """

    # Convert list of values to a vector in R.
    series_r = rpy2.robjects.FloatVector(series)
    # Calculate the Hartigans' dip test statistic.
    dip = diptest.dip(series_r, min_is_0=True)[0]
    #dip = diptest.dip_test(series_r, simulate_p_value=True, B=1000)
    return dip


def generate_random_values_normal(
    mean=None, deviation=None, count=None, method=None
):
    """
    Generates a series of random values from a normal, Gaussian, distribution.

    arguments:
        mean (float): value of mean, center, of normal distribution
        deviation (float): value of standard deviation of normal distribution
        count (int): count of values to generate
        method (str): method to use, either "random" or "numpy"

    raises:

    returns:
        (list<float>): series of values

    """

    if method == "random":
        values = []
        for value in range(count):
            values.append(random.gauss(mean, deviation))
    elif method == "numpy":
        values = list(
            numpy.random.normal(loc=mean, scale=deviation, size=count)
        )
    return values


def print_metric_report(series=None):
    """
    Calculates and prints reports on multiple metrics.

    arguments:
        series (list<float>): series of values of type float

    raises:

    returns:

    """

    # Calculate skewness.
    skewness = scipy.stats.skew(series, axis=None)
    print("skewness: " + str(skewness))
    # Calculate kurtosis.
    kurtosis = scipy.stats.kurtosis(series, axis=None, fisher=True)
    print("kurtosis: " + str(kurtosis))
    # Calculate bimodality coefficient.
    bimodality = calculate_bimodality_coefficient(series)
    print("bimodality coefficient: " + str(bimodality))
    # Calculate dip statistic.
    dip = calculate_dip_statistic(series=series)
    print("dip statistic: " + str(dip))

    pass


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
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    path_aggregation = os.path.join(
        path_organization, "data_gene_signal_aggregation.pickle"
    )
    path_log = os.path.join(
        path_organization, "data_gene_signal_log.pickle"
    )
    # Write information to file.
    with open(path_imputation, "wb") as file_product:
        pickle.dump(
            information["data_gene_signal_tissue_median"], file_product
        )
    with open(path_aggregation, "wb") as file_product:
        pickle.dump(information["data_gene_signal_aggregation"], file_product)
    with open(path_log, "wb") as file_product:
        pickle.dump(information["data_gene_signal_log"], file_product)


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
    #source = read_source(dock=dock)

    utility.print_terminal_partition(level=1)
    print("Test of metrics of modality.")

    # Unimodal normal distribution.

    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 10,000,000 random values with a unimodal normal " +
        "distribution."
    )
    print("Expectations for unimodal normal distribution...")
    print("skewness = 0")
    print("kurtosis = 0")
    print("bimodality coefficient < 0.55")
    print("dip statistic < 0.05")
    utility.print_terminal_partition(level=3)
    # Generate random values with a normal distribution.
    series = generate_random_values_normal(
        mean=1.0,
        deviation=3.0,
        count=1000000,
        method="random"
    )
    print_metric_report(series=series)
    utility.print_terminal_partition(level=3)

    # Bimodal normal distribution.
    # ...


    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 10,000,000 random values with a bimodal normal " +
        "distribution."
    )
    print("Expectations for bimodal normal distribution...")
    print("skewness = ?")
    print("kurtosis = ?")
    print("bimodality coefficient > 0.55")
    print("dip statistic > 0.05")
    utility.print_terminal_partition(level=3)
    # Generate random values with a normal distribution.
    series_one = generate_random_values_normal(
        mean=1.0,
        deviation=1.0,
        count=1000000,
        method="random"
    )
    series_two = generate_random_values_normal(
        mean=15.0,
        deviation=2.0,
        count=1000000,
        method="random"
    )
    #series_one.extend(series_two)
    series = series_one + series_two
    print_metric_report(series=series)
    utility.print_terminal_partition(level=3)

    # Compile information.
    information = {}
    #Write product information to file.
    #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
