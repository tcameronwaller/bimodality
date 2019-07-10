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
import scipy
#import cmapPy.pandasGEXpress.parse
#import cmapPy.pandasGEXpress.gct2gctx

# Custom

import measurement
import shuffle
import metric
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
    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    path_organization = os.path.join(dock, "organization")
    path_organization_signal = os.path.join(
        path_organization, "data_gene_persons_tissues_signals.pickle"
    )
    path_restriction = os.path.join(dock, "restriction")
    path_restriction_signal = os.path.join(
        path_restriction, "data_gene_persons_tissues_signals.pickle"
    )
    # Read information from file.
    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    data_organization_signals = pandas.read_pickle(path_organization_signal)
    data_restriction_signals = pandas.read_pickle(path_restriction_signal)
    # Compile and return information.
    return {
        "shuffles": shuffles,
        "data_organization_signals": data_organization_signals,
        "data_restriction_signals": data_restriction_signals,
    }


##########
# Transformation of gene's signals to logarithmic space.


def calculate_standard_score_gene_signal_by_tissue(
    data_gene_persons_tissues_signals=None
):
    """
    Calculates the standard (z-score) of genes' signals across each tissue.

    The standard scores are relative to tissue.
    The values of mean and standard deviation are across all samples (patients)
    for each tissue.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    data_gene_standard = data_gene_persons_tissues_signals.apply(
        lambda x: (x - x.mean()) / x.std(),
        axis="columns",
    )
    return data_gene_standard


def standardize_gene_signal(
    data_gene_persons_tissues_signals=None
):
    """
    Transforms values of genes' signals to standard or z-score space.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Transform signals to standard score space.
    data_gene_signal_standard = calculate_standard_score_gene_signal_by_tissue(
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals
    )
    if False:
        # Compare summary statistics before and after transformation.
        utility.print_terminal_partition(level=3)
        print("Summary statistics for gene signals before standardization.")
        print(data_gene_persons_tissues_signals.iloc[0:10, 0:10])
        data_mean = data_gene_persons_tissues_signals.apply(
            lambda x: x.mean(),
            axis="columns"
        )
        print("Mean")
        print(data_mean.iloc[0:10])
        data_deviation = data_gene_persons_tissues_signals.apply(
            lambda x: x.std(),
            axis="columns"
        )
        print("Standard deviation")
        print(data_deviation.iloc[0:10])
        utility.print_terminal_partition(level=3)
        print("Summary statistics for gene signals after standardization.")
        print(data_gene_signal_standard.iloc[0:10, 0:10])
        data_mean = data_gene_signal_standard.apply(
            lambda x: x.mean(),
            axis="columns"
        )
        print("Mean")
        print(data_mean.iloc[0:10])
        data_deviation = data_gene_signal_standard.apply(
            lambda x: x.std(),
            axis="columns"
        )
        print("Standard deviation")
        print(data_deviation.iloc[0:10])
    # Return information.
    return data_gene_signal_standard


def normalize_standardize_gene_signal(
    data_gene_persons_tissues_signals=None
):
    """
    Normalizes and standardizes genes' signals across persons and tissues

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    if False:
        utility.print_terminal_partition(level=2)
        print(
            "Notice that transformation to logarithmic and standard space " +
            "propagates missing values."
        )

    # Transform genes' signals to logarithmic space.
    # Transform values of genes' signals to base-2 logarithmic space.
    data_normal = measurement.transform_gene_signal_log(
        data_gene_signal=data_gene_persons_tissues_signals,
        pseudo_count=1.0,
    )
    # Transform genes' signals to standard or z-score space.
    # Calculate standard score for each gene.
    # Standard score is undefined if the series has inadequate variance.
    # Consider filtering before calculation.
    # Alternatively, filter afterwards.
    data_normal_standard = standardize_gene_signal(
        data_gene_persons_tissues_signals=data_normal
    )
    return data_normal_standard


##########
# Aggregation of gene's signals across tissues for each person.


def aggregate_gene_persons_signals(
    data_gene_persons_tissues_signals=None
):
    """
    Aggregates a gene's signals across persons.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate signals across
            persons

    """

    # Aggregate gene's signals across tissues for each person.
    # Exclude missing values from calculation of the mean.
    series_aggregate = data_gene_persons_tissues_signals.aggregate(
        lambda x: x.mean(),
        axis="columns",
    )
    data_aggregate = series_aggregate.to_frame(name="value")
    # Return information.
    return data_aggregate


##########
# Distribution.


def calculate_bimodality_metrics(
    values=None
):
    """
    Calculates bimodality metrics for a distribution.

    arguments:
        values (list): values of a distribution

    raises:

    returns:
        (dict<float>): Values of metrics of bimodality.

    """

    # Calculate metrics of bimodality.
    coefficient = metric.calculate_bimodality_coefficient(series=values)
    dip = metric.calculate_dip_statistic(series=values)
    mixture = metric.calculate_mixture_model_score(series=values)

    # Compile information.
    information = {
        "coefficient": coefficient,
        "dip": dip,
        "mixture": mixture,
    }

    # Return information.
    return information


def determine_distribution_values_scores(
    data_gene_persons_signals=None
):
    """
    Determines metrics of a distribution's modality.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (dict<float>): Values of metrics of bimodality.

    """

    # Extract values of distribution.
    # Remove missing values.
    data = data_gene_persons_signals.copy(deep=True)
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    values = data["value"].to_list()

    # Convert the scale of the values.
    if False:
        values_scale = list(scipy.stats.zscore(values))
        values_scale = numpy.apply_along_axis(
            lambda x: (x - x.mean()) / x.std(),
            0,
            numpy.array(values)
        )
        values_scale = numpy.apply_along_axis(
            lambda x: x * 10000000000,
            0,
            numpy.array(values)
        )

    # Calculate metrics of bimodality.
    scores = calculate_bimodality_metrics(values=values)

    # Prepare gene report.
    # Describe gene's aggregate signals across persons.
    report_gene = prepare_report_gene(
        data_gene_persons_signals=data,
    )

    # Compile information.
    information = {
        "report_gene": report_gene,
        "values": values,
        "scores": scores
    }
    # Return information.
    return information



##########
# Report.


def prepare_report_gene(
    data_gene_persons_signals=None,
):
    """
    Prepare report for gene.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (dict): report for gene

    """

    # Prepare gene report.
    # Describe gene's aggregate signals across persons.

    # Compile information.
    information = {
        "data_gene_persons_signals": data_gene_persons_signals,
    }
    # Return information.
    return information




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
    print("Welcome to the test, demonstration of the distribution procedure!")

    # Read source information from file.
    source = read_source(dock=dock)
    print(source["data_restriction_signals"])

    # Call procedure.
    information = execute_procedure(
        data_gene_persons_tissues_signals=(
            source["data_restriction_signals"]
        ),
    )

    print("Printing distribution values...")
    print(information["values"])

    print("Printing bimodality scores...")
    print(information["scores"])

    print("Printing gene report... specifically data_gene_persons_signals")
    print(information["report_gene"]["data_gene_persons_signals"])

    utility.print_terminal_partition(level=2)
    print("Now let's test the shuffle procedure...")

    print("data before shuffle...")
    print(source["data_organization_signals"])

    # Shuffle gene's signals.
    utility.print_terminal_partition(level=3)
    print("index procedure...")
    data_shuffle = shuffle.shuffle_gene_signals(
        data_gene_signals=source["data_organization_signals"],
        shuffle=source["shuffles"][0]
    )
    print("data after shuffle")
    print(data_shuffle)


    pass


def execute_procedure(
    data_gene_persons_tissues_signals=None
):
    """
    Function to execute module's main behavior.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): report about gene, Pandas data frame of gene's signals across
            persons and tissues

    """

    # If data come from the "availability" method of the restriction procedure,
    # they might include missing values.

    # Normalize and standardize gene's signals.
    # Transform gene's signals to base-two logarithmic space.
    # Transform gene's signals to standard, z-score space.
    data_normal_standard = normalize_standardize_gene_signal(
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals,
    )

    # Aggregate gene's signals across tissues from each person.
    # The method for aggregation across tissues is to calculate the mean.
    # This aggregation value is a mean of standard z-scores.
    # Hence the aggregation value is very near zero.
    data_gene_persons_signals = aggregate_gene_persons_signals(
        data_gene_persons_tissues_signals=data_normal_standard,
    )


    # TODO: Keep track of final distribution values for each patient's identifier...
    # TODO: construct a report for the distribution procedure...

    # Determine values and scores of distribution.
    # Rescale the values before calculating modality scores on the
    # distribution.
    collection = determine_distribution_values_scores(
        data_gene_persons_signals=data_gene_persons_signals
    )

    # Compile information.
    information = {
        "values": collection["values"],
        "scores": collection["scores"],
        "report_gene": collection["report_gene"],
    }
    # Return information.
    return information


if (__name__ == "__main__"):
    execute_procedure()