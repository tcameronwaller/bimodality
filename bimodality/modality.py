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
import textwrap

# Relevant

import numpy
import pandas
import scipy.stats
if False:
    import rpy2
    import rpy2.robjects
    import rpy2.robjects.packages
    #utils = rpy2.robjects.packages.importr("utils")
    #base = rpy2.robjects.packages.importr("base")
    #utils.install_packages("diptest")
    diptest = rpy2.robjects.packages.importr("diptest")
import sklearn
import sklearn.mixture
import unidip
#import diptest

# Custom

import utility
import plot

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
    # Correct for statistical sample bias.
    skewness = scipy.stats.skew(
        series,
        axis=None,
        bias=False,
    )
    # Calculate excess kurtosis.
    # Correct for statistical sample bias.
    kurtosis = scipy.stats.kurtosis(
        series,
        axis=None,
        fisher=True,
        bias=False,
    )
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

    # Calculate the Hartigans' dip test statistic.
    # Python package "unidip".
    dip = unidip.dip.dip_fn(series, is_hist=False, just_dip=True)
    # Python package "diptest".
    if False:
        dip = diptest.dip(
            numpy.array(series),
            full_output=False,
            min_is_0=True,
            x_is_sorted=False,
            )
    #print("dip: " + str(dip))
    return dip


def calculate_dip_statistic_r(series=None):
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


def calculate_mixture_model_score(series=None):
    """
    Calculates the score from a Gaussian mixture model of a series.

    arguments:
        series (list<float>): series of values of type float

    raises:

    returns:
        (float): value of mixture model score

    """

    # Organize series.
    array = numpy.asarray(series)
    array_shape = numpy.reshape(array, (-1, 1))
    # Fit series to models.
    model_one = sklearn.mixture.GaussianMixture(
        n_components=1,
        covariance_type="full",
        n_init=3,
        max_iter=100,
    )
    model_two = sklearn.mixture.GaussianMixture(
        n_components=2,
        covariance_type="full",
        n_init=3,
        max_iter=100,
    )
    model_one.fit(array_shape)
    model_two.fit(array_shape)
    # In general, the model with the greater log likelihood is a better
    # representation of the data.
    # Determine log likelihood scores for fits to models.
    likelihood_log_one = model_one.score(array_shape)
    likelihood_log_two = model_two.score(array_shape)
    # Calculate log likelihood ratio of the bimodal model relative to the
    # unimodal model.
    # Null hypothesis is that the unimodal model (more restrictive) is better.
    # Do not use the value of the likelihood ratio test.
    # Rather use the value of the ratio itself, which is the statistic.
    # Per standard, multiply log likelihood ratio by 2 to adjust distribution.
    #likelihood_one = math.exp(likelihood_log_one)
    #likelihood_two = math.exp(likelihood_log_two)
    #ratio = -2*math.log(likelihood_one / likelihood_two)
    ratio = 2*(likelihood_log_two - likelihood_log_one)
    return ratio


# Calculate metric on basis of mixture model.


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


def prepare_metric_report_text(
    skewness=None,
    kurtosis=None,
    coefficient=None,
    dip=None,
    mixture=None,
):
    """
    Prepares a summary report on curation of metabolic sets and entities.

    arguments:
        skewness (float): value of skewness
        kurtosis (float): value of kurtosis
        coefficient (float): value of bimodality coefficient
        dip (float): value of Hardigan dip statistic
        mixture (float): value of likelihood ratio from mixture models

    returns:
        (str): report of summary information

    raises:

    """

    # Compile information.
    report = textwrap.dedent("""\

        --------------------------------------------------

        skewness: {skewness}
        kurtosis: {kurtosis}
        bimodality coefficient: {coefficient}
        dip statistic: {dip}
        mixture model: {mixture}

        --------------------------------------------------
    """).format(
        skewness=skewness,
        kurtosis=kurtosis,
        coefficient=coefficient,
        dip=dip,
        mixture=mixture,
    )
    # Return information.
    return report


def report_metrics(
    name=None,
    series=None,
    dock=None,
):
    """
    Calculates and prints reports on multiple metrics.

    arguments:
        name (str): name of series
        series (list<float>): series of values of type float
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Calculate skewness.
    skewness = scipy.stats.skew(series, axis=None)
    # Calculate kurtosis.
    kurtosis = scipy.stats.kurtosis(series, axis=None, fisher=True)
    # Calculate bimodality coefficient.
    coefficient = calculate_bimodality_coefficient(series)
    # Calculate dip statistic.
    dip = calculate_dip_statistic(series=series)
    # Calculate mixture model score.
    mixture = calculate_mixture_model_score(series=series)

    # Prepare metric report text.
    text = prepare_metric_report_text(
        skewness=skewness,
        kurtosis=kurtosis,
        coefficient=coefficient,
        dip=dip,
        mixture=mixture,
    )
    print(text)

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_distribution_histogram(
        series=series,
        name="",
        bin_method="count",
        bin_count=50,
        label_bins="Bins",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=False,
        position=0,
        text=text,
    )
    # Specify directories and files.
    path_metric = os.path.join(dock, "metric")
    path_figure = os.path.join(path_metric, "figure")
    utility.create_directory(path_figure)
    file = (name + "_distribution.svg")
    path_file = os.path.join(path_figure, file)
    # Write figure.
    plot.write_figure(
        path=path_file,
        figure=figure
    )

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
    utility.create_directory(path_organization)
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

    # Remove previous files.
    # Specify directories and files.
    path_metric = os.path.join(dock, "metric")
    utility.create_directory(path_metric)
    path_figure = os.path.join(path_metric, "figure")
    utility.remove_directory(path=path_figure)

    # Read source information from file.
    #source = read_source(dock=dock)

    utility.print_terminal_partition(level=1)
    print("Test of metrics of modality.")

    # Unimodal normal distribution.

    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 1,000,000 random values with a unimodal normal " +
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
    report_metrics(
        name="unimodality",
        series=series,
        dock=dock
    )
    utility.print_terminal_partition(level=3)

    # Bimodal normal distribution 1.
    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 1,000,000 random values with a bimodal normal " +
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
        count=500000,
        method="random"
    )
    series_two = generate_random_values_normal(
        mean=5.0,
        deviation=2.0,
        count=500000,
        method="random"
    )
    #series_one.extend(series_two)
    series = series_one + series_two
    report_metrics(
        name="bimodality_1",
        series=series,
        dock=dock
    )
    utility.print_terminal_partition(level=3)

    # Bimodal normal distribution 2.
    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 1,000,000 random values with a bimodal normal " +
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
        count=500000,
        method="random"
    )
    series_two = generate_random_values_normal(
        mean=10.0,
        deviation=2.0,
        count=500000,
        method="random"
    )
    #series_one.extend(series_two)
    series = series_one + series_two
    report_metrics(
        name="bimodality_2",
        series=series,
        dock=dock
    )
    utility.print_terminal_partition(level=3)

    # Bimodal normal distribution 3.
    utility.print_terminal_partition(level=2)
    print(
        "Simulation on 1,000,000 random values with a bimodal normal " +
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
        count=100000,
        method="random"
    )
    series_two = generate_random_values_normal(
        mean=10.0,
        deviation=2.0,
        count=900000,
        method="random"
    )
    #series_one.extend(series_two)
    series = series_one + series_two
    report_metrics(
        name="bimodality_3",
        series=series,
        dock=dock
    )
    utility.print_terminal_partition(level=3)

    # Compile information.
    information = {}
    #Write product information to file.
    #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
