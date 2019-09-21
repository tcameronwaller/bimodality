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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import measurement
import metric
import category
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


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
    path_assembly = os.path.join(dock, "assembly")
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_selection = os.path.join(dock, "selection")
    path_families = os.path.join(path_selection, "data_families.pickle")
    if source_genes == "split":
        path_source = os.path.join(dock, "split")
    elif source_genes == "combination":
        path_source = os.path.join(dock, "combination")
    path_genes = os.path.join(
        path_source, "genes.txt"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_families = pandas.read_pickle(path_families)
    genes = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes,
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_families": data_families,
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
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Organization


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

    # Pandas sets the new columns as categorical type.
    # This categorical type is restrictive.
    # Organize data.

    data_gene_persons_tissues_signals.columns = (
        data_gene_persons_tissues_signals.columns.astype(str)
    )
    data_gene_persons_tissues_signals.rename_axis(
        columns="",
        axis="columns",
        copy=False,
        inplace=True
    )
    if False:
        data_gene_persons_tissues_signals.reset_index(
            level="person",
            drop=False,
            inplace=True
        )
        data_gene_persons_tissues_signals.set_index(
            ["person"],
            append=False,
            drop=True,
            inplace=True
        )

    # Return information.
    return data_gene_persons_tissues_signals


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


def evaluate_organization(
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


def organize_data(
    data_gene_samples_signals=None,
):
    """
    Organizes gene's signals across samples, persons, and tissues.

    arguments:
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): information from organization

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
    # Describe dispersion of gene's signals across major tissue categories.
    pack = evaluate_organization(
        data_gene_samples_signals=data_gene_samples_signals,
    )
    # Compile information.
    information = {
        "data_gene_signal_tissue": data_gene_signal_tissue,
        "data_gene_persons_tissues_signals": data_gene_persons_tissues_signals,
        "data_gene_tissue_mean": pack["data_gene_tissue_mean"],
        "data_gene_tissue_variance": pack["data_gene_tissue_variance"],
        "data_gene_tissue_deviation": pack["data_gene_tissue_deviation"],
        "data_gene_tissue_variation": pack["data_gene_tissue_variation"],
        "data_gene_tissue_dispersion": pack["data_gene_tissue_dispersion"],
    }
    # Return information.
    return information


##########
# Restriction


def remove_null_persons_tissues(
    data_gene_persons_tissues_signals=None
):
    """
    Removes persons and tissues with only missing values.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Copy data.
    data_copy = data_gene_persons_tissues_signals.copy(deep=True)

    if False:
        utility.print_terminal_partition(level=2)
        print(
            "Remove persons and tissues with only missing values."
        )

        print(
            "shape of original data frame: " +
            str(data_gene_persons_tissues_signals.shape)
        )
    data_copy.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_copy.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    if False:
        print(
            "shape of data frame: " +
            str(data_copy.shape)
        )
    return data_copy


def select_persons_tissues(
    method=None,
    count=None,
    tissues=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Select persons by availability of valid values of gene's signal for
    specific tissues.

    Data format should have tissues across columns and persons across rows.

    arguments:
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        tissues (list<str>): specific tissues to select by "imputation" method
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Copy data.
    data_original = data_gene_persons_tissues_signals.copy(deep=True)

    if method == "availability":
        data_selection = data_original.dropna(
            axis="index",
            how="all",
            thresh=count,
            inplace=False,
        )
        if False:
            data_nonzero = (data_selection != 0)
            data_selection = (
                data_selection.loc[data_nonzero.any(axis="columns"), : ]
            )

        pass
    elif method == "imputation":
        # Select tissues of interest.
        data_tissues = data_original.loc[ :, tissues]
        # Select patients with coverage of tissues.
        data_selection = data_tissues.dropna(
            axis="index",
            how="all",
            thresh=count,
            inplace=False,
        )
        if False:
            data_nonzero = (data_selection != 0)
            data_selection = (
                data_selection.loc[data_nonzero.any(axis="columns"), : ]
            )
        pass
    # Return information.
    return data_selection


def impute_gene_persons_tissues(
    data_selection=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Select persons by availability of valid values of gene's signal for
    specific tissues.

    arguments:
        data_selection (object): Pandas data frame of a gene's signals across
            selection persons and tissues
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (object): Pandas data frame of a gene's signals across persons and
            tissues

    """

    # Copy data.
    data_copy = data_selection.copy(deep=True)

    # Calculate values for imputation.
    # Calculate values for imputation from all available samples.
    imputations = data_gene_persons_tissues_signals.aggregate(
        lambda x: x.median()
    )
    #imputations = series_imputation.to_dict()

    # Insert imputations to selections of persons and tissues.
    if False:
        data_copy.apply(
            lambda x: x.fillna(
                imputations[x.name],
                inplace=True,
            ),
            axis="index",
        )
    data_imputation = data_copy.fillna(
        value=imputations,
        #axis="columns",
        inplace=False,
    )

    # Return information.
    return data_imputation


def evaluate_restriction(
    data_gene_persons_tissues_signals=None,
):
    """
    Prepare report for gene.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): report for gene

    """

    # Prepare gene report.
    # Describe mean count of tissues across persons.
    # Describe counts of tissues across persons.
    # Describe specific tissues for which persons have valid signals for
    # gene without imputation.

    # Describe specific tissues for which each person has valid signals for
    # gene without imputation.
    data_persons_tissues = data_gene_persons_tissues_signals.applymap(
        lambda x: (True) if (pandas.notna(x)) else (False)
    )

    # Describe counts of tissues for which each person has valid signals for
    # gene without imputation.
    data_persons_tissues_count = data_gene_persons_tissues_signals.aggregate(
        lambda x: x.dropna().size,
        axis="columns",
    )

    # Calculate count of persons.
    persons = data_persons_tissues.shape[0]

    # Calculate mean count of tissues per person.
    mean = data_persons_tissues_count.mean()

    # Calculate median count of tissues per person.
    median = data_persons_tissues_count.median()

    # Compile information.
    information = {
        "data_gene_persons_tissues": data_persons_tissues,
        "data_gene_persons_tissues_count": data_persons_tissues_count,
        "persons_count": persons,
        "tissues_mean": mean,
        "tissues_median": median,
    }
    # Return information.
    return information


def restrict_data(
    method=None,
    count=None,
    tissues=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Restricts gene's signals across samples, persons, and tissues.

    arguments:
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        count (int): minimal count of tissues with signal coverage for
            selection by "availability" or "imputation" methods
        tissues (list<str>): specific tissues to select by "imputation" method
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information from restriction

    """

    # Filter to select persons and tissues with valid values of gene's signals.
    data_filter = remove_null_persons_tissues(
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals
    )

    # Select persons and tissues for subsequent aggregation and analysis.
    # There are multiple methods for this selection.
    if method == "availability":
        # Select tissues and patients.
        data_selection = select_persons_tissues(
            method=method,
            count=count,
            tissues=tissues,
            data_gene_persons_tissues_signals=data_filter,
        )
        # Prepare gene report.
        # Describe mean count of tissues across persons.
        # Describe counts of tissues across persons.
        # Describe specific tissues for which persons have valid signals for
        # gene without imputation.
        pack = evaluate_restriction(
            data_gene_persons_tissues_signals=data_selection,
        )

        pass
    elif method == "imputation":
        # Select tissues and patients.
        data_temporary = select_persons_tissues(
            method=method,
            count=count,
            tissues=tissues,
            data_gene_persons_tissues_signals=data_filter,
        )
        # Impute missing values.
        data_selection = impute_gene_persons_tissues(
            data_selection=data_temporary,
            data_gene_persons_tissues_signals=data_filter,
        )

        # Prepare gene report.
        # Describe mean count of tissues across persons.
        # Describe counts of tissues across persons.
        # Describe specific tissues for which persons have valid signals for
        # gene without imputation.
        # Importantly, the gene's report is on the basis of the data before
        # imputation.
        pack = evaluate_restriction(
            data_gene_persons_tissues_signals=data_temporary,
        )

        pass

    # Compile information.
    information = {
        "data_gene_persons_tissues_signals_restriction": data_selection,
        "data_gene_persons_tissues": pack["data_gene_persons_tissues"],
        "data_gene_persons_tissues_count": (
            pack["data_gene_persons_tissues_count"]
        ),
        "persons_count_restriction": pack["persons_count"],
        "tissues_count_mean": pack["tissues_mean"],
        "tissues_count_median": pack["tissues_median"],
    }
    # Return information.
    return information


##########
# Aggregation


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

    # Calculate standard score.
    # This method inserts missing values if the standard deviation is zero.
    data_gene_standard = data_gene_persons_tissues_signals.apply(
        lambda x: (
            ((x - x.mean()) / x.std())
        ),
        axis="index",
    )
    return data_gene_standard


def standardize_gene_signal(
    data_gene_persons_tissues_signals=None
):
    """
    Transforms values of genes' signals to standard or z-score space.

    Data has persons across rows and tissues across columns.

             tissue_1 tissue_2 tissue_3 tissue_4 tissue_5
    person_1 ...      ...      ...      ...      ...
    person_2 ...      ...      ...      ...      ...
    person_3 ...      ...      ...      ...      ...
    person_4 ...      ...      ...      ...      ...
    person_5 ...      ...      ...      ...      ...

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
        data_mean = data_gene_persons_tissues_signals.aggregate(
            lambda x: x.mean(),
            axis="index"
        )
        print("Mean")
        print(data_mean.iloc[0:10])
        data_deviation = data_gene_persons_tissues_signals.aggregate(
            lambda x: x.std(),
            axis="index"
        )
        print("Standard deviation")
        print(data_deviation.iloc[0:10])
        utility.print_terminal_partition(level=3)
        print("Summary statistics for gene signals after standardization.")
        print(data_gene_signal_standard.iloc[0:10, 0:10])
        data_mean = data_gene_signal_standard.aggregate(
            lambda x: x.mean(),
            axis="index"
        )
        print("Mean")
        print(data_mean.iloc[0:10])
        data_deviation = data_gene_signal_standard.aggregate(
            lambda x: x.std(),
            axis="index"
        )
        print("Standard deviation")
        print(data_deviation.iloc[0:10])
    # Return information.
    return data_gene_signal_standard


# Transformation of gene's signals to logarithmic space.
def normalize_standardize_gene_signal(
    data_gene_persons_tissues_signals=None
):
    """
    Normalizes and standardizes genes' signals across persons and tissues.

    Data has persons across rows and tissues across columns.

             tissue_1 tissue_2 tissue_3 tissue_4 tissue_5
    person_1 ...      ...      ...      ...      ...
    person_2 ...      ...      ...      ...      ...
    person_3 ...      ...      ...      ...      ...
    person_4 ...      ...      ...      ...      ...
    person_5 ...      ...      ...      ...      ...

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


def determine_gene_persons_signals(
    data_gene_persons_signals=None
):
    """
    Determines a gene's values of tissue-aggregate scores across persons.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (dict): gene's tissue-aggregate scores across persons

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

    # Compile information.
    information = {
        "data_gene_persons_signals": data,
        "values": values,
    }
    # Return information.
    return information


def evaluate_aggregation(
    data_gene_persons_signals=None,
):
    """
    Prepare report for gene.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict): report for gene

    """

    # Prepare gene report.
    # Describe gene's aggregate signals across persons.

    # Calculate count of persons.
    persons = data_gene_persons_signals.shape[0]

    # Compile information.
    information = {
        "persons_count_aggregation": persons,
        "data_gene_persons_signals": data_gene_persons_signals,
    }
    # Return information.
    return information


def aggregate_data(
    data_gene_persons_tissues_signals=None,
):
    """
    Aggregates genes' signals within tissues across persons.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information from aggregation

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

    # Determine values of gene's tissue-aggregate signals across persons.
    # This procedure includes basic filters.
    collection = determine_gene_persons_signals(
        data_gene_persons_signals=data_gene_persons_signals
    )

    # Prepare gene report.
    # Describe gene's aggregate signals across persons.
    pack = evaluate_aggregation(
        data_gene_persons_signals=collection["data_gene_persons_signals"],
    )

    # Compile information.
    information = {
        "data_gene_persons_signals": collection["data_gene_persons_signals"],
        "values": collection["values"],
        "persons_count_aggregation": pack["persons_count_aggregation"],
    }
    # Return information.
    return information


##########
# Modality


def evaluate_gene_population(
    threshold=None,
    count=None,
    data_gene_persons_signals=None
):
    """
    Evaluates the a gene's representation across persons.

    arguments:
        threshold (float): minimal signal value
        count (int): minimal count of persons with signal beyond threshold
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene has representation across an adequate
            population of persons

    """

    # Copy data.
    data = data_gene_persons_signals.copy(
        deep=True
    )
    # Determine whether values exceed threshold.
    data_threshold = (data != threshold)
    data_pass = data.loc[data_threshold["value"], :]
    count_persons = data_pass.shape[0]
    return (count_persons >= count)


def evaluate_gene_signal(
    threshold=None,
    data_gene_persons_signals=None
):
    """
    Evaluates the gene's aggregate tissue signals across persons.

    arguments:
        threshold (float): minimal standard deviation
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene's signals are adequate for analysis

    """

    # Calculate standard deviation of values.
    deviation = numpy.std(data_gene_persons_signals["value"].values, axis=0)
    return deviation > threshold


# This function is currently excluded from evaluation of candidacy.
# TODO: this function will eventually need to accommodate additional attributes of persons...
# 1. variation in gene's signal across tissue categories could indicate problems in aggregation of minor tissues
# 2. impose a constraint on the mean count of tissues across persons
# 3.
def evaluate_gene_tissue(
    threshold_proportion=None,
    threshold_probability=None,
    data_samples_tissues_persons=None,
    data_gene_persons_tissues=None,
    data_gene_persons_signals=None,
):
    """
    Evaluates the selection of tissues for calculation of gene's aggregate
    scores.

    arguments:
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a tissue's group in analysis
        threshold_probability (float): minimal probability for candidacy
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene's selection of tissues is acceptable

    """

    # Organize information about gene's scores across persons and tissues.
    data_persons_tissues_signals = category.organize_gene_persons_signals(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=data_gene_persons_tissues,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine minimal count of persons for inclusion of a tissue's group.
    count_persons = data_gene_persons_signals.shape[0]
    threshold_count = threshold_proportion * count_persons


def calculate_bimodality_metrics(
    values=None
):
    """
    Calculates bimodality metrics for a distribution.

    arguments:
        values (list): values of a distribution

    raises:

    returns:
        (dict<float>): values of metrics of bimodality

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


def generate_null_metrics():
    """
    Generates null metrics.

    arguments:

    raises:

    returns:
        (dict<float>): values of metrics of bimodality

    """

    # Compile information.
    information = {
        "coefficient": float("nan"),
        "dip": float("nan"),
        "mixture": float("nan"),
    }

    # Return information.
    return information


# TODO: improve the candidacy checks on distributions...
def describe_distribution_modality(
    modality=None,
    values=None,
    data_gene_persons_signals=None,
):
    """
    Applies multiple metrics to describe the distribution of a gene's signals
    across persons.

    arguments:
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        values (list<float>): values of a gene's signals across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict): scores of modality metrics

    """

    # Determine distribution modality metrics.
    # Determine whether to calculate metrics for gene's distribution.
    if modality:
        population = evaluate_gene_population(
            threshold=0.0,
            count=100,
            data_gene_persons_signals=data_gene_persons_signals,
        )
        signal = evaluate_gene_signal(
            threshold=0.0,
            data_gene_persons_signals=data_gene_persons_signals,
        )
        if signal and population:
            # Calculate metrics.
            # Calculate metrics of bimodality.
            scores = calculate_bimodality_metrics(
                values=values
            )
            pass
        else:
            # Generate missing values.
            scores = generate_null_metrics()
            pass
    else:
        # Generate missing values.
        scores = generate_null_metrics()
        pass

    # Compile information.
    information = {
        "scores": scores,
    }

    # Return information.
    return information


##########
# Extraction


def extract_gene_persons_signals(
    data_gene_persons_signals=None,
    data_families=None,
):
    """
    Extracts information about a gene's distribution of pan-tissue signals
    across persons.

    This information is useful for heritability analysis in GCTA.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons
        data_families (object): Pandas data frame of families and persons

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate, pan-tissue signals
            across families and persons

    """

    # Treat the complete table of families and persons as master.
    # Introduce missing values for persons without valid signals.
    # Join
    data_families.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    data_gene_families_persons_signals = data_families.join(
        data_gene_persons_signals,
        how="left",
        on="person"
    )
    data_gene_families_persons_signals.reset_index(
        level=None,
        inplace=True
    )
    return data_gene_families_persons_signals


##########
# Report


def calculate_combination_scores(
    data_report=None,
):
    """
    Converts genes' modality metrics to standard, z-score space and then
    calculates combination metrics.

    arguments:
        data_report (object): Pandas data frame of information about genes'
            distributions

    raises:

    returns:
        (object): Pandas data frame of information about genes' distributions

    """

    # Copy data.
    data = data_report.copy(deep=True)
    # Calculate mean and standard deviation of each modality metric across
    # genes.
    coefficients = data["coefficient"].to_list()
    dips = data["dip"].to_list()
    mixtures = data["mixture"].to_list()
    mean_coefficient = statistics.mean(coefficients)
    deviation_coefficient = statistics.stdev(coefficients)
    mean_dip = statistics.mean(dips)
    deviation_dip = statistics.stdev(dips)
    mean_mixture = statistics.mean(mixtures)
    deviation_mixture = statistics.stdev(mixtures)
    # Convert gene's modality metrics to standard, z-score, space across genes.
    data["coefficient_z"] = data.apply(
        lambda row: utility.calculate_standard_score(
            value=row["coefficient"],
            mean=mean_coefficient,
            deviation=deviation_coefficient,
        ),
        axis="columns",
    )
    data["dip_z"] = data.apply(
        lambda row: utility.calculate_standard_score(
            value=row["dip"],
            mean=mean_dip,
            deviation=deviation_dip,
        ),
        axis="columns",
    )
    data["mixture_z"] = data.apply(
        lambda row: utility.calculate_standard_score(
            value=row["mixture"],
            mean=mean_mixture,
            deviation=deviation_mixture,
        ),
        axis="columns",
    )
    # Calculate combination score as mean of other modality metrics.
    data["combination"] = data.apply(
        lambda row: statistics.mean(
            [
                row["coefficient_z"],
                row["dip_z"],
                row["mixture_z"],
            ]
        ),
        axis="columns",
    )
    # Return information.
    return data


def prepare_gene_report(
    gene=None,
    persons_restriction=None,
    persons_aggregation=None,
    tissues=None,
    tissues_mean=None,
    tissues_median=None,
    coefficient=None,
    dip=None,
    mixture=None,
    data_gene_annotation=None,
):
    """
    Prepares a report about the gene's distribution of aggregate pan-tissue
    signals across persons.

    arguments:
        gene (str): identifier of gene
        persons_restriction (int): count of persons after restriction
        persons_aggregation (int): count of persons after aggregation
        tissues (int): count of tissues
        tissues_mean (float): mean count of tissues across persons
        tissues_median (float): median count of tissues across persons
        coefficient (float): bimodality coefficient
        dip (float): Hardigan's dip statistic
        mixture (float): log likelihood ratio from Gaussian mixture models for
            one or two modes
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (dict): information about gene's signals

    """

    # Name
    name = data_gene_annotation.loc[gene, "gene_name"]
    # Chromosome
    chromosome = data_gene_annotation.loc[gene, "seqname"]
    # Coordinates
    start = data_gene_annotation.loc[gene, "start"]
    end = data_gene_annotation.loc[gene, "end"]
    # Compile information.
    information = {
        "gene": gene,
        "name": name,
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "persons_restriction": persons_restriction,
        "persons_aggregation": persons_aggregation,
        "tissues": tissues,
        "tissues_mean": tissues_mean,
        "tissues_median": tissues_median,
        "coefficient": coefficient,
        "dip": dip,
        "mixture": mixture,
    }
    # Return information.
    return information


##########################################################*******************
########## More or less scrap beyond here...
# Distribution

def determine_gene_distributions(
    gene=None,
    modality=None,
    data_gene_samples_signals=None,
):
    """
    Prepares and describes the distributions of a gene's signals across
    persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Organize data for analysis.
    collection_organization = organization.execute_procedure(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Determine distributions of gene's signals by multiple methods.
    imputation = determine_gene_distribution(
        gene=gene,
        modality=modality,
        method="imputation",
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        ),
    )
    availability = determine_gene_distribution(
        gene=gene,
        modality=modality,
        method="availability",
        data_gene_persons_tissues_signals=(
            collection_organization["data_gene_persons_tissues_signals"]
        ),
    )

    # Compile information.
    information = {
        "organization": collection_organization,
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


# Use this function as a reference to extract information for a gene's
# distribution.
def extract_gene_distribution_information(
    method=None,
    observation=None,
):
    """
    Extracts information about a gene's distribution of pan-tissue signals
    across persons.

    arguments:
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

    raises:

    returns:

    """

    # Access information.
    report_organization = observation["organization"]["report_gene"]
    report_restriction = observation[method]["report_restriction"]
    report_aggregation = observation[method]["report_aggregation"]

    # Gene's signals across original, complete set of persons and tissues.
    data_gene_persons_tissues_signals_complete = (
        observation["organization"]["data_gene_persons_tissues_signals"]
    )
    # Mean of gene's signals across major tissue categories.
    data_gene_tissue_mean = report_organization["data_gene_tissue_mean"]
    # Variance of gene's signals across major tissue categories.
    data_gene_tissue_variance = (
        report_organization["data_gene_tissue_variance"]
    )
    # Standard deviation of gene's signals across major tissue categories.
    data_gene_tissue_deviation = (
        report_organization["data_gene_tissue_deviation"]
    )
    # Coefficient of variation of gene's signals across major tissue
    # categories.
    data_gene_tissue_variation = (
        report_organization["data_gene_tissue_variation"]
    )
    # Dispersion in gene's signals across major tissue categories.
    data_gene_tissue_dispersion = (
        report_organization["data_gene_tissue_dispersion"]
    )



    # Gene's signals across tissues and persons, before any imputation.
    data_gene_persons_tissues_signals = (
        report_restriction["data_gene_persons_tissues_signals"]
    )
    # Boolean availability of valid gene's signals across tissues and persons.
    data_gene_persons_tissues = (
        report_restriction["data_gene_persons_tissues"]
    )
    # Count of tissues for which gene has valid signals across persons.
    data_gene_persons_tissues_count = (
        report_restriction["data_gene_persons_tissues_count"]
    )
    # Count of persons for which gene has valid signals across adequate
    # tissues.
    restriction_persons = report_restriction["persons"]
    # Mean count of tissues per person.
    restriction_tissues_mean = report_restriction["tissues_mean"]
    # Median count of tissues per person.
    restriction_tissues_median = report_restriction["tissues_median"]



    # Gene's pan-tissue aggregate signals across persons.
    data_gene_persons_signals = (
        report_aggregation["data_gene_persons_signals"]
    )
    # Count of persons for which gene has valid signals across adequate
    # tissues.
    aggregation_persons = report_aggregation["persons"]

    pass


##########
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

    ##########
    # Specify directories and files.
    path_distribution = os.path.join(dock, "distribution")
    utility.create_directory(path_distribution)
    path_gene = os.path.join(path_distribution, gene)
    utility.create_directory(path_gene)
    # General
    path_tissues = os.path.join(
        path_gene, "tissues.txt"
    )
    path_report = os.path.join(
        path_gene, "report.pickle"
    )
    # Organization
    path_data_gene_signal_tissue = os.path.join(
        path_gene, "data_gene_signal_tissue.pickle"
    )
    path_data_gene_persons_tissues_signals = os.path.join(
        path_gene, "data_gene_persons_tissues_signals.pickle"
    )
    path_data_gene_tissue_mean = os.path.join(
        path_gene, "data_gene_tissue_mean.pickle"
    )
    path_data_gene_tissue_variance = os.path.join(
        path_gene, "data_gene_tissue_variance.pickle"
    )
    path_data_gene_tissue_deviation = os.path.join(
        path_gene, "data_gene_tissue_deviation.pickle"
    )
    path_data_gene_tissue_variation = os.path.join(
        path_gene, "data_gene_tissue_variation.pickle"
    )
    path_data_gene_tissue_dispersion = os.path.join(
        path_gene, "data_gene_tissue_dispersion.pickle"
    )
    # Restriction
    path_data_gene_persons_tissues_signals_restriction = os.path.join(
        path_gene, "data_gene_persons_tissues_signals_restriction.pickle"
    )
    path_data_gene_persons_tissues = os.path.join(
        path_gene, "data_gene_persons_tissues.pickle"
    )
    path_data_gene_persons_tissues_count = os.path.join(
        path_gene, "data_gene_persons_tissues_count.pickle"
    )
    # Aggregation
    path_data_gene_persons_signals = os.path.join(
        path_gene, "data_gene_persons_signals.pickle"
    )
    path_data_gene_families_persons_signals = os.path.join(
        path_gene, "data_gene_families_persons_signals.pickle"
    )
    path_data_gene_families_persons_signals_text = os.path.join(
        path_gene, "data_gene_families_persons_signals.tsv"
    )
    # Modality
    path_scores = os.path.join(
        path_gene, "scores.pickle"
    )

    ##########
    # Write information to file.
    # General
    utility.write_file_text_list(
        elements=information["tissues"],
        delimiter="\n",
        path_file=path_tissues
    )
    with open(path_report, "wb") as file_product:
        pickle.dump(information["report"], file_product)
    # Organization
    pandas.to_pickle(
        information["data_gene_signal_tissue"],
        path_data_gene_signal_tissue
    )
    pandas.to_pickle(
        information["data_gene_persons_tissues_signals"],
        path_data_gene_persons_tissues_signals
    )
    pandas.to_pickle(
        information["data_gene_tissue_mean"],
        path_data_gene_tissue_mean
    )
    pandas.to_pickle(
        information["data_gene_tissue_variance"],
        path_data_gene_tissue_variance
    )
    pandas.to_pickle(
        information["data_gene_tissue_deviation"],
        path_data_gene_tissue_deviation
    )
    pandas.to_pickle(
        information["data_gene_tissue_variation"],
        path_data_gene_tissue_variation
    )
    pandas.to_pickle(
        information["data_gene_tissue_dispersion"],
        path_data_gene_tissue_dispersion
    )
    # Restriction
    pandas.to_pickle(
        information["data_gene_persons_tissues_signals_restriction"],
        path_data_gene_persons_tissues_signals_restriction
    )
    pandas.to_pickle(
        information["data_gene_persons_tissues"],
        path_data_gene_persons_tissues
    )
    pandas.to_pickle(
        information["data_gene_persons_tissues_count"],
        path_data_gene_persons_tissues_count
    )
    # Aggregation
    pandas.to_pickle(
        information["data_gene_persons_signals"],
        path_data_gene_persons_signals
    )
    pandas.to_pickle(
        information["data_gene_families_persons_signals"],
        path_data_gene_families_persons_signals
    )
    information["data_gene_families_persons_signals"].to_csv(
        path_or_buf=path_data_gene_families_persons_signals_text,
        columns=["family", "person", "value"],
        sep="\t",
        na_rep="NA",
        header=False,
        index=False,
    )
    # Modality
    with open(path_scores, "wb") as file_product:
        pickle.dump(information["scores"], file_product)

    pass


def read_collect_write_gene_report(
    genes=None,
    dock=None,
):
    """
    Collects information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Check contents of directory.
    utility.print_terminal_partition(level=2)
    print("Check that directories exist for all genes.")
    path_distribution = os.path.join(dock, "distribution")
    directories = os.listdir(path_distribution)
    match = utility.compare_lists_by_mutual_inclusion(
        list_one=genes, list_two=directories
    )
    print("Genes and directories match: " + str(match))
    utility.print_terminal_partition(level=2)

    # Collect information about genes.
    records = list()
    # Iterate on directories for genes.
    for gene in directories:
        # Specify directories and files.
        path_gene = os.path.join(path_distribution, gene)
        path_report = os.path.join(
            path_gene, "report.pickle"
        )
        # Read information from file.
        with open(path_report, "rb") as file_source:
            record = pickle.load(file_source)
        # Create entry for gene.
        records.append(record)
    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records,
    )
    # Convert modality metrics to standard, z-score space.
    # Calculate combination modality metrics.
    data_report = calculate_combination_scores(
        data_report=data,
    )
    data_report.sort_values(
        by=["combination"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    utility.print_terminal_partition(level=2)
    print(data_report)
    utility.print_terminal_partition(level=2)
    # Specify directories and files.
    path_data_report = os.path.join(
        path_distribution, "data_gene_report.pickle"
    )
    path_data_report_text = os.path.join(
        path_distribution, "data_gene_report.tsv"
    )
    # Write information to file.
    data_report.to_pickle(
        path=path_data_report
    )
    data_report.to_csv(
        path_or_buf=path_data_report_text,
        sep="\t",
        header=True,
        index=False,
    )

    pass


###############################################################################
# Procedure

# TODO: export gene's aggregate pan-tissue signals across persons
# TODO: export gene's chromosome and coordinates
# TODO: heritability analysis will only consider genes on autosomes (non-sex chromosomes)

# TODO: this function should not accept permutations...

def execute_procedure(
    gene=None,
    data_gene_samples_signals=None,
    data_families=None,
    data_gene_annotation=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of gene
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        data_families (object): Pandas data frame of person's identifiers and
            families' identifiers
        data_gene_annotation (object): Pandas data frame of genes' annotations
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    print(data_gene_samples_signals)

    ##########
    # Organization
    bin_organization = organize_data(
        data_gene_samples_signals=data_gene_samples_signals
    )

    ##########
    # Restriction
    # Define tissues for inclusion.
    # This selection includes all non-sexual tissues with coverage of samples
    # from at least 100 persons.
    tissues = [
        "adipose", # 797
        "adrenal", # 258
        "artery", # 786
        "blood", # 767
        "brain", # 382
        #"breast", # 459
        "colon", # 555
        "esophagus", # 710
        "heart", # 561
        "intestine", # 187
        "liver", # 226
        "lung", # 578
        "muscle", # 803
        "nerve", # 619
        #"ovary", # 180
        "pancreas", # 328
        "pituitary", # 283
        #"prostate", # 245
        "salivary", # 162
        "skin", # 912
        "spleen", # 241
        "stomach", # 359
        #"testis", # 361
        "thyroid", # 653
        #"uterus", # 142
        #"vagina", # 156
    ]
    # 18 September 2019
    # count   persons
    # 1       948
    # 2       943
    # 3       930
    # 4       922
    # 5       911
    # 6       885
    # 7       857
    # 8       796
    # 9       737
    # 10      645
    # 11      540
    # 12      419
    # 13      297
    # 14      180
    # 15      92
    # 16      34
    # 17      8
    # 18      3
    # 19      0
    bin_restriction = restrict_data(
        method="availability", # "availability" or "imputation"
        count=11,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            bin_organization["data_gene_persons_tissues_signals"]
        ),
    )

    # Aggregation
    # Determine gene's distributions of aggregate tissue signals across
    # persons.
    # Use the final "data_gene_persons_tissues_signals" from restriction
    # procedure after imputation.
    bin_aggregation = aggregate_data(
        data_gene_persons_tissues_signals=(
            bin_restriction["data_gene_persons_tissues_signals_restriction"]
        ),
    )

    # Modality
    bin_modality = describe_distribution_modality(
        modality=True,
        values=bin_aggregation["values"],
        data_gene_persons_signals=bin_aggregation["data_gene_persons_signals"],
    )

    # Extraction
    data_gene_families_persons_signals = extract_gene_persons_signals(
        data_gene_persons_signals=bin_aggregation["data_gene_persons_signals"],
        data_families=data_families,
    )

    print("testing now...")

    print(data_families)
    print(bin_aggregation["data_gene_persons_signals"])
    print(data_gene_families_persons_signals)

    # Report
    report = prepare_gene_report(
        gene=gene,
        persons_restriction=bin_restriction["persons_count_restriction"],
        persons_aggregation=bin_aggregation["persons_count_aggregation"],
        tissues=len(tissues),
        tissues_mean=bin_restriction["tissues_count_mean"],
        tissues_median=bin_restriction["tissues_count_median"],
        coefficient=bin_modality["scores"]["coefficient"],
        dip=bin_modality["scores"]["dip"],
        mixture=bin_modality["scores"]["mixture"],
        data_gene_annotation=data_gene_annotation,
    )

    # Compile information.
    information = bin_organization.copy()
    information.update(bin_restriction)
    information.update(bin_aggregation)
    information.update(bin_modality)
    information["gene"] = gene
    information["tissues"] = tissues
    information["data_gene_families_persons_signals"] = (
        data_gene_families_persons_signals
    )
    information["report"] = report
    # Write product information to file.
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
    print("... Distribution procedure ...")

    # Report date and time.
    utility.print_terminal_partition(level=3)
    start = datetime.datetime.now()
    print(start)
    utility.print_terminal_partition(level=3)

    # Remove previous files to avoid version or batch confusion.
    path_distribution = os.path.join(dock, "distribution")
    utility.remove_directory(path=path_distribution)

    # Read source information from file.
    source = read_source_initial(
        source_genes="split",
        dock=dock
    )
    print("count of genes: " + str(len(source["genes"])))

    if True:
        report = execute_procedure_local_sub(
            gene="ENSG00000231925", # TAPBP
            data_families=source["data_families"],
            data_gene_annotation=source["data_gene_annotation"],
            dock=dock,
        )
        report = execute_procedure_local_sub(
            gene="ENSG00000134184", # TAPBP
            data_families=source["data_families"],
            data_gene_annotation=source["data_gene_annotation"],
            dock=dock,
        )
    if False:
        # Set up partial function for iterative execution.
        # Each iteration uses a different sequential value of the "gene" variable
        # with the same value of the "dock" variable.
        execute_procedure_gene = functools.partial(
            execute_procedure_local_sub,
            data_families=source["data_families"],
            data_gene_annotation=source["data_gene_annotation"],
            dock=dock,
        )
        # Initialize multiprocessing pool.
        #pool = multiprocessing.Pool(processes=os.cpu_count())
        pool = multiprocessing.Pool(processes=8)
        # Iterate on genes.
        check_genes=[
            "ENSG00000231925", # TAPBP
        ]
        #report = pool.map(execute_procedure_gene, check_genes)
        report = pool.map(execute_procedure_gene, source["genes"][0:1000])
        #report = pool.map(execute_procedure_gene, source["genes"])

    # Pause procedure.
    time.sleep(10.0)

    # Report.
    read_collect_write_gene_report(
        genes=source["genes"],
        dock=dock,
    )
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    # TODO: execute a collection procedure to assemble genes and a report...

    # Report date and time.
    utility.print_terminal_partition(level=3)
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))
    utility.print_terminal_partition(level=3)

    pass


def execute_procedure_local_sub(
    gene=None,
    data_families=None,
    data_gene_annotation=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        data_families (object): Pandas data frame of person's identifiers and
            families' identifiers
        data_gene_annotation (object): Pandas data frame of genes' annotations
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
        data_gene_samples_signals=source["data_gene_samples_signals"],
        data_families=data_families,
        data_gene_annotation=data_gene_annotation,
        dock=dock
    )

    # Report progress.
    path_distribution = os.path.join(dock, "distribution")
    directories = os.listdir(path_distribution)
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
        permutations=source["permutations"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
