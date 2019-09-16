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

import metric
import candidacy
import category
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
    path_permutation = os.path.join(dock, "permutation")
    path_permutations = os.path.join(
        path_permutation, "permutations.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    with open(path_permutations, "rb") as file_source:
        permutations = pickle.load(file_source)
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_gene_samples_signals": data_gene_samples_signals,
        "permutations": permutations,
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
        "persons": persons,
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
        "data_gene_persons_tissues": pack["data_persons_tissues"],
        "data_gene_persons_tissues_count": pack["data_persons_tissues_count"],
        "persons_restriction": pack["persons"],
        "tissues_mean": pack["tissues_mean"],
        "tissues_median": pack["tissues_median"],
    }
    # Return information.
    return information


##########
# Aggregation


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
        "persons": persons,
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
        data_gene_persons_signals=collection["data"],
    )

    # Compile information.
    information = {
        "data": collection["data"],
        "values": collection["values"],
        "report_gene": report_gene,
    }
    # Return information.
    return information


##########
# Distribution


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


# This function includes important parameters.
# Tissues for selection in restriction procedure.
# Count of tissues to require coverage in restriction procedure.
def determine_gene_distribution(
    gene=None,
    modality=None,
    method=None,
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distributions of a gene's signals across
    persons.

    arguments:
        gene (str): identifier of gene
        modality (bool): whether to calculate metrics for the modality of
            gene's distribution
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information about the distribution of a gene's signals across
            persons

    """

    # Specify tissues for the imputation method of the restriction procedure.
    # This selection includes all non-sexual tissues with coverage of samples
    # from at least 200 persons.
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
        #"spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ]

    # Restriction
    collection_restriction = restriction.execute_procedure(
        method=method,
        count=10,
        tissues=tissues,
        data_gene_persons_tissues_signals=(
            data_gene_persons_tissues_signals
        ),
    )

    # Aggregation
    # Use the final "data_gene_persons_tissues_signals" from restriction
    # procedure after imputation.
    collection_aggregation = aggregation.execute_procedure(
        data_gene_persons_tissues_signals=(
            collection_restriction["data_gene_persons_tissues_signals"]
        ),
    )

    # Determine distribution modality metrics.
    # Determine whether to calculate metrics for gene's distribution.
    if modality:
        population = candidacy.evaluate_gene_population(
            threshold=100,
            data_gene_persons_signals=collection_aggregation["data"],
        )
        signal = candidacy.evaluate_gene_signal(
            threshold=0.0,
            data_gene_persons_signals=collection_aggregation["data"],
        )
        if signal and population:
            # Calculate metrics.
            # Calculate metrics of bimodality.
            scores = calculate_bimodality_metrics(
                values=collection_aggregation["values"]
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
        "report_restriction": collection_restriction["report_gene"],
        "values": collection_aggregation["values"],
        "scores": scores,
        "report_aggregation": collection_aggregation["report_gene"],
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
    path_distribution = os.path.join(dock, "distribution")
    utility.create_directory(path_distribution)
    path_gene = os.path.join(path_distribution, gene)
    utility.create_directory(path_gene)
    path_imputation = os.path.join(path_gene, "imputation")
    utility.create_directory(path_imputation)
    path_availability = os.path.join(path_gene, "availability")
    utility.create_directory(path_availability)

    path_report_organization = os.path.join(
        path_gene, "report_organization.pickle"
    )
    # Imputation method
    path_imputation_report_restriction = os.path.join(
        path_imputation, "report_restriction.pickle"
    )
    path_imputation_report_aggregation = os.path.join(
        path_imputation, "report_aggregation.pickle"
    )
    path_imputation_scores = os.path.join(
        path_imputation, "scores.pickle"
    )
    path_imputation_permutations = os.path.join(
        path_imputation, "permutations.pickle"
    )
    # Availability method
    path_availability_report_restriction = os.path.join(
        path_availability, "report_restriction.pickle"
    )
    path_availability_report_aggregation = os.path.join(
        path_availability, "report_aggregation.pickle"
    )
    path_availability_scores = os.path.join(
        path_availability, "scores.pickle"
    )
    path_availability_permutations = os.path.join(
        path_availability, "permutations.pickle"
    )

    # Write information to file.
    with open(path_report_organization, "wb") as file_product:
        pickle.dump(information["report_organization"], file_product)

    # Imputation method
    with open(path_imputation_report_restriction, "wb") as file_product:
        pickle.dump(information["imputation"]["report_restriction"], file_product)
    with open(path_imputation_report_aggregation, "wb") as file_product:
        pickle.dump(
            information["imputation"]["report_aggregation"], file_product
        )
    with open(path_imputation_scores, "wb") as file_product:
        pickle.dump(information["imputation"]["scores"], file_product)
    with open(path_imputation_permutations, "wb") as file_product:
        pickle.dump(information["imputation"]["permutations"], file_product)

    # Imputation method
    with open(path_availability_report_restriction, "wb") as file_product:
        pickle.dump(
            information["availability"]["report_restriction"], file_product
        )
    with open(path_availability_report_aggregation, "wb") as file_product:
        pickle.dump(
            information["availability"]["report_aggregation"], file_product
        )
    with open(path_availability_scores, "wb") as file_product:
        pickle.dump(information["availability"]["scores"], file_product)
    with open(path_availability_permutations, "wb") as file_product:
        pickle.dump(information["availability"]["permutations"], file_product)

    pass


###############################################################################
# Procedure

# Change procedure to run and collect both by "availability" and "imputation" methods...


def execute_procedure(
    gene=None,
    permutations=None,
    data_gene_samples_signals=None,
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of gene
        permutations (list<list<list<int>>>): matrices of permutation indices
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    ##########
    # Organization
    bin_organization = organize_data(
        data_gene_samples_signals=data_gene_samples_signals
    )

    ##########
    # Restriction
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
        #"spleen", # 162
        "stomach", # 262
        "thyroid", # 446
    ]
    bin_restriction = restrict_data(
        method="availability", # "availability" or "imputation"
        count=10,
        tissues=tissues,
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Aggregation
    # Use the final "data_gene_persons_tissues_signals" from restriction
    # procedure after imputation.
    bin_aggregation = aggregate_data(
        data_gene_persons_tissues_signals=(
            bin_restriction["data_gene_persons_tissues_signals_restriction"]
        ),
    )



    # Modality


    # TODO: prepare a table with all original persons with missing values if they were'nt used in final distribution


    # Compilation


    ##################Old Stuff...###########################################3

    # Prepare and describe distribution of real gene's signals.

    # Determine gene's distributions of aggregate tissue signals across
    # persons.
    observation = determine_gene_distributions(
        gene=gene,
        modality=True,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Compile information.
    information_imputation = {
        "report_restriction": (
            observation["imputation"]["report_restriction"]
        ),
        "report_aggregation": (
            observation["imputation"]["report_aggregation"]
        ),
        "scores": observation["imputation"]["scores"],
        "permutations": shuffle["imputation"],
    }
    information_availability = {
        "report_restriction": (
            observation["availability"]["report_restriction"]
        ),
        "report_aggregation": (
            observation["availability"]["report_aggregation"]
        ),
        "scores": observation["availability"]["scores"],
        "permutations": shuffle["availability"],
    }
    information = {
        "report_organization": observation["organization"]["report_gene"],
        "imputation": information_imputation,
        "availability": information_availability,
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
    source = read_source_local_initial(source_genes="split", dock=dock)

    print("count of genes: " + str(len(source["genes"])))
    #print("count of permutations: " + str(len(source["permutations"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    # Remove previous files to avoid version or batch confusion.
    path_distribution = os.path.join(dock, "distribution")
    utility.remove_directory(path=path_distribution)

    # Set up partial function for iterative execution.
    # Each iteration uses the same values for "genes_signals", "permutations", and
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
    pool = multiprocessing.Pool(processes=100)

    # Iterate on genes.
    check_genes=[
        "ENSG00000198965",
    ]
    #report = pool.map(execute_procedure_gene, check_genes)
    #report = pool.map(execute_procedure_gene, source["genes"][0:32])
    report = pool.map(execute_procedure_gene, source["genes"])

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
        permutations=source["permutations"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    # Report contents of directory.
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
