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
import random
import copy

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import modality
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(dock=None):
    """
    Initialize directories for procedure's product files.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = dock
    paths["distribution"] = os.path.join(paths["dock"], "distribution")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["distribution"])
    utility.create_directory(path=paths["distribution"])
    # Define paths for groups of persons.
    groups = list()
    groups.append("selection")
    groups.append("ventilation")
    for group in groups:
        paths[group] = dict()
        paths[group]["genes"] = os.path.join(
            paths["distribution"], group, "genes"
        )
        paths[group]["collection"] = os.path.join(
            paths["distribution"], group, "collection"
        )
        paths[group]["trait"] = os.path.join(
            paths["distribution"], group, "trait"
        )
        # Initialize directories.
        utility.create_directories(path=paths[group]["genes"])
        utility.create_directories(path=paths[group]["collection"])
        utility.create_directories(path=paths[group]["trait"])
    # Return information.
    return paths


##########
# Source


def read_source_initial(
    group=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        group (str): group of persons, either selection or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )
    path_persons_sets = os.path.join(
        dock, "selection", "tight", "persons_properties",
        "persons_sets.pickle"
    )
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", group,
        "data_persons_properties.pickle"
    )
    path_data_families_persons = os.path.join(
        dock, "selection", "tight", "persons_properties", group,
        "heritability", "data_families_persons.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    with open(path_persons_sets, "rb") as file_source:
        persons_sets = pickle.load(file_source)
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_families_persons = pandas.read_pickle(path_data_families_persons)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "persons": persons_sets[group],
        "data_persons_properties": data_persons_properties,
        "data_families_persons": data_families_persons,
    }


# TODO: I think I should make this read separate collections for different person selections...
def read_source(
    gene=None,
    group=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process
        group (str): group of persons, either selection or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(
        dock, "split", group, "collection", (gene + ".pickle")
    )
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


def calculate_gene_signal_mean_by_persons_tissues_scrap(
    persons_tissues_samples=None,
    data_gene_signal=None
):
    """
    Calculates the mean values of genes' signals in each tissue of each
    person.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # bob     skin    5.0    10.0   7.0
    # bob     liver   5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # april   liver   2.0    3.0    4.0
    # april   skin    2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal (object): Pandas data frame of gene's signals for all
            samples.

    raises:

    returns:
        (object): Pandas data frame of mean values of genes' signals for all
            tissues of all persons.

    """

    # Create index for convenient access.
    data_gene_signal_index = data_gene_signal.set_index("Name")
    # Determine list of all genes with signals.
    genes = data_gene_signal["Name"].to_list()
    # Collect mean signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons.
    for person in persons_tissues_samples:
        # Iterate on tissues for which person has samples.
        for tissue in persons_tissues_samples[person]:
            # Collect signals for all genes in each sample.
            signals_genes = dict()
            for sample in persons_tissues_samples[person][tissue]:
                for gene in genes:
                    signal = data_gene_signal_index.at[gene, sample]
                    # Determine whether an entry already exists for the gene.
                    if gene in signals_genes:
                        signals_genes[gene].append(signal)
                    else:
                        signals_genes[gene] = list([signal])
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Calculate mean signal for each gene in each person and tissue.
            for gene in genes:
                signal_mean = statistics.mean(signals_genes[gene])
                record[gene] = signal_mean
            # Include record.
            summary.append(record)
    # Organize data.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary


def calculate_gene_signal_median_by_tissues_scrap(
    tissues=None,
    persons=None,
    tissues_persons_samples=None,
    data_gene_signal_mean=None
):
    """
    Calculates the median values of signals for all genes across specific
    persons and tissues.

    Redundant samples exist for many persons and tissues.
    Calculate the mean signal for each gene from all of these redundant
    samples.

    Product data format.
    Index by gene allows convenient access.
    ##################################################
    # gene   brain colon liver heart kidney
    # abc1   7     3     2     6     3
    # mpc1   5     3     2     6     3
    # rip1   3     3     2     6     3
    # rnt1   2     3     2     6     3
    # rnx3   9     3     2     6     3
    # sst5   3     3     2     6     3
    # sph6   2     3     2     6     3
    # xtr4   8     3     2     6     3
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        tissues_persons_samples (dict<dict<list<str>>): Samples for each
            person of each tissue.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of median values of genes' signals for all
            tissues.

    """

    # Create index for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect median signals of all genes for each tissue.
    summary = list()
    # Iterate on genes.
    for gene in genes:
        # Compile record.
        record = dict()
        record["gene"] = gene
        # Iterate on tissues.
        for tissue in tissues:
            # Collect gene's signals in the tissue across all persons.
            signals_tissue = list()
            # Iterate on persons.
            for person in tissues_persons_samples[tissue]:
                signal_tissue_person = data_gene_signal_mean_index.loc[
                    (person, tissue), gene
                ]
                signals_tissue.append(signal_tissue_person)
            # Calculate median value of gene's signal in tissue across all
            # persons.
            signal_tissue = statistics.median(signals_tissue)
            record[tissue] = signal_tissue
        # Include record.
        summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    return data_summary


def collect_gene_signal_by_persons_tissues_scrap(
    tissues=None,
    persons=None,
    persons_tissues_samples=None,
    data_gene_signal_median_tissue=None,
    data_gene_signal_mean=None
):
    """
    Collects the signals for all genes in each tissue of each person.

    Product data format.
    Indices by person and tissue allow convenient access.
    ##################################################
    # index   index   gene_1 gene_2 gene_3
    # person tissue
    # bob     kidney  5.0    10.0   7.0
    # april   brain   2.0    3.0    4.0
    # martin  liver   1.0    6.0    9.0
    # sue     heart   3.0    2.0    5.0
    ##################################################

    arguments:
        tissues (list<str>): Tissues of interest.
        persons (list<str>): persons with signal for tissues of interest.
        persons_tissues_samples (dict<dict<list<str>>): Samples for each
            tissue of each person.
        data_gene_signal_median_tissue (object): Pandas data frame of median
            values of genes' signals for all tissues.
        data_gene_signal_mean (object): Pandas data frame of mean values of
            genes' signals for all tissues of all persons.

    raises:

    returns:
        (object): Pandas data frame of signals for all genes across
            specific persons and tissues.

    """

    # Create indices for convenient access.
    data_gene_signal_mean_index = data_gene_signal_mean.set_index(
        ["person", "tissue"], append=False, drop=True
    )
    data_gene_signal_median_tissue_index = (
        data_gene_signal_median_tissue.set_index("gene")
    )
    # Extract genes from names of columns.
    genes = data_gene_signal_mean_index.columns.to_list()
    # Collect signals of all genes for each person and tissue.
    summary = list()
    # Iterate on persons of interest.
    for person in persons:
        # Iterate on tissues of interest.
        for tissue in tissues:
            # Compile record.
            record = dict()
            record["person"] = person
            record["tissue"] = tissue
            # Determine whether the person has samples for the tissue.
            if tissue in persons_tissues_samples[person]:
                # Collect mean values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_mean_index.loc[
                        (person, tissue), gene
                    ]
                    record[gene] = signal
            else:
                # Collect imputation values for genes' signals.
                # Iterate on genes.
                for gene in genes:
                    signal = data_gene_signal_median_tissue_index.loc[
                        gene, tissue
                    ]
                    record[gene] = signal
            # Include record.
            summary.append(record)
    # Organize information within a data frame.
    data_summary = utility.convert_records_to_dataframe(
        records=summary
    )
    #data_summary_index = data_summary.set_index(
    #    ["person", "tissue"], append=False, drop=False
    #)
    return data_summary


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
        lambda series: scipy.stats.zscore(
            series.to_numpy(),
            axis=0,
            ddof=1, # Sample standard deviation.
            nan_policy="omit",
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


def transform_gene_signal_pseudo_logarithm(
    data_gene_signal=None,
):
    """
    Calculates the base 2.0 logarithm of genes' signals in each sample.

    Original gene signals are in transcript counts per million (TPM).

    To accommodate gene signals of 0.0, add a pseudo count of 1.0 to all counts
    before calculation of the base 2.0 logarithm.

    arguments:
        data_gene_signal (object): Pandas data frame of genes' signals

    raises:

    returns:
        (object): Pandas data frame of genes' signals in logarithmic space

    """

    # Transform genes' signals to base-2 logarithmic space.
    # Transform before calculation of any median or mean values. <- maybe false
    # Logarithms are not distributive.
    # To accommodate gene signals of 0.0, add a pseudo count of 2.0 to all
    # counts before calculation of the base-2 logarithm.
    # log-2 signal = log2(TPM + pseudo-count)
    # pseudo-count = 1.0
    #utility.print_terminal_partition(level=2)
    #print("Transformation of genes' signals to base-2 logarithmic space.")
    #print(
    #    "To accommodate gene signals of 0.0, add a pseudo count of " +
    #    "1.0-5.0 to all counts before calculation of the base-2 logarithm."
    #)
    data_gene_signal_log = utility.calculate_pseudo_logarithm_signals(
        pseudo_count=1.0,
        base=2.0,
        data=data_gene_signal,
    )
    return data_gene_signal_log


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
    data_normal = transform_gene_signal_pseudo_logarithm(
        data_gene_signal=data_gene_persons_tissues_signals,
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

    data = data_gene_persons_tissues_signals.copy(deep=True)

    # If data come from the "availability" method of the restriction procedure,
    # they might include missing values.

    # Normalize and standardize gene's signals.
    # Transform gene's signals to base-two logarithmic space.
    # Transform gene's signals to standard, z-score space.
    data_normal_standard = normalize_standardize_gene_signal(
        data_gene_persons_tissues_signals=data,
    )

    # Aggregate gene's signals across tissues from each person.
    # The method for aggregation across tissues is to calculate the mean.
    # This aggregation value is a mean of standard z-scores.
    # Hence the aggregation value is very near zero.
    data_gene_persons_signals = aggregate_gene_persons_signals(
        data_gene_persons_tissues_signals=data_normal_standard,
    )

    # Determine values of gene's tissue-aggregate signals across persons.
    # This procedure includes basic filters to remove missing values.
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
        "data_gene_signals_tissues_persons_normal_standard": (
            data_normal_standard
        ),
        "data_gene_persons_signals": data_gene_persons_signals, # These data include signals across all persons.
        "values": collection["values"], # These values include only nonmissing signals across persons.
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


def measure_bimodality(
    values=None
):
    """
    Calculates bimodality measures for a distribution.

    arguments:
        values (list): values of a distribution

    raises:

    returns:
        (dict<float>): values of measures of bimodality

    """

    # Calculate measures of bimodality.
    coefficient = modality.calculate_bimodality_coefficient(series=values)
    dip = modality.calculate_dip_statistic(series=values)
    mixture = modality.calculate_mixture_model_score(series=values)
    # Compile information.
    information = {
        "coefficient": coefficient,
        "dip": dip,
        "mixture": mixture,
    }

    # Return information.
    return information


def generate_null_measures():
    """
    Generates null measures.

    arguments:

    raises:

    returns:
        (dict<float>): values of measures of bimodality

    """

    # Compile information.
    information = {
        "coefficient": float("nan"),
        "dip": float("nan"),
        "mixture": float("nan"),
    }

    # Return information.
    return information


def evaluate_gene_candidacy(
    data_gene_persons_signals=None,
):
    """
    Evaluates a gene's candidacy for further bimodal analysis by multiple
    criteria.

    1. candidacy by population
    - gene must have tissue-aggregate signals across at least 100 persons
    2. candidacy by signal
    - gene's tissue-aggregate signals across persons must have a standard
        deviation greater than 0
    3. candidacy by tissue independence
    - consider tissues with representation in at least 5% of persons

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict<bool>): information about gene's candidacy

    """

    # Determine whether the gene's distribution is adequate.
    # Distribution procedure applies filters to the gene's tissue-aggregate
    # scores across persons.

    # Determine whether gene has representation across adequate persons.
    population = evaluate_gene_population(
        threshold=0.0,
        count=100,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine whether gene's tissue-aggregate scores across persons have
    # adequate variance to apply measures to the distribution.
    signal = evaluate_gene_signal(
        threshold=0.0,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    return (population and signal)


# TODO: improve the candidacy checks on distributions...
def describe_distribution_modality(
    modality=None,
    values=None,
    data_gene_persons_signals=None,
):
    """
    Applies multiple measures to describe the distribution of a gene's signals
    across persons.

    arguments:
        modality (bool): whether to calculate measures for the modality of
            gene's distribution
        values (list<float>): values of a gene's signals across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict): scores of modality measures

    """

    # Determine distribution modality measures.
    # Determine whether to calculate measures for gene's distribution.
    if modality:
        candidacy = evaluate_gene_candidacy(
            data_gene_persons_signals=data_gene_persons_signals,
        )
        if candidacy:
            # Calculate measures.
            # Calculate measures of bimodality.
            scores = measure_bimodality(
                values=values
            )
            pass
        else:
            # Generate missing values.
            scores = generate_null_measures()
            pass
    else:
        # Generate missing values.
        scores = generate_null_measures()
        pass

    # Compile information.
    information = {
        "scores": scores,
    }

    # Return information.
    return information


##########
# Distribution
# The distribution subprocedure combines restriction, aggregation, and modality
# functions.
# This subprocedure is also useful in the permutation procedure.


def prepare_describe_distribution(
    data_gene_persons_tissues_signals=None,
):
    """
    Prepares and describes the distribution of a gene's signals across samples,
    persons, and tissues.

    arguments:
        data_gene_persons_tissues_signals (object): Pandas data frame of a
            gene's signals across persons and tissues

    raises:

    returns:
        (dict): information from restriction, aggregation, and modality

    """

    ##########
    # Restriction
    # Define tissues for inclusion.
    # This selection includes all non-sexual tissues with coverage of samples
    # from at least 100 persons.
    method = "availability" # "availability" or "imputation"
    count = 10
    tissues = [
        "adipose", "adrenal", "artery", "blood", "brain", "colon",
        "esophagus", "heart", "intestine", "liver", "lung", "muscle",
        "nerve", "pancreas", "pituitary", "salivary", "skin", "spleen",
        "stomach", "thyroid",
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
        method=method,
        count=count,
        tissues=tissues,
        data_gene_persons_tissues_signals=data_gene_persons_tissues_signals,
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

    # Compile information.
    information = {
        "method": method,
        "count": count,
        "tissues": tissues,
        "bin_restriction": bin_restriction,
        "bin_aggregation": bin_aggregation,
        "bin_modality": bin_modality,
    }
    # Return information.
    return information


##########
# Format for heritability in GCTA


def extract_gene_persons_signals(
    data_gene_persons_signals=None,
    data_families_persons=None,
):
    """
    Extracts information about a gene's distribution of pan-tissue signals
    across persons.

    This information is useful for heritability analysis in GCTA.

    arguments:
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons
        data_families_persons (object): Pandas data frame of families and persons

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate, pan-tissue signals
            across families and persons

    """

    # Copy data before modification.
    data_families_persons = data_families_persons.copy(deep=True)

    # Treat the complete table of families and persons as master.
    # Introduce missing values for persons without valid signals.
    # Join
    data_gene_families_persons_signals = data_families_persons.join(
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
    Converts genes' modality measures to standard, z-score space and then
    calculates combination measures.

    arguments:
        data_report (object): Pandas data frame of information about genes'
            distributions

    raises:

    returns:
        (object): Pandas data frame of information about genes' distributions

    """

    # Copy data.
    data = data_report.copy(deep=True)
    # Calculate mean and standard deviation of each modality measure across
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
    # Convert gene's modality measures to standard, z-score, space across genes.
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
    # Calculate combination score as mean of other modality measures.
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
    method=None,
    count=None,
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
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        count (int): minimal count of tissues with signal coverage for
            selection by "availability" or "imputation" methods
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
        "identifier": gene,
        "name": name,
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "persons_restriction": persons_restriction,
        "persons_aggregation": persons_aggregation,
        "method": method,
        "count": count,
        "tissues": tissues,
        "tissues_mean": tissues_mean,
        "tissues_median": tissues_median,
        "coefficient": coefficient,
        "dip": dip,
        "mixture": mixture,
    }
    # Return information.
    return information


##########
# Product


def write_product_gene(
    gene=None,
    group=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        group (str): group of persons, either selection or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    ##########
    # Specify directories and files.
    path_gene = os.path.join(paths[group]["genes"], gene)
    utility.create_directory(path_gene)

    # General
    path_tissues = os.path.join(
        path_gene, "tissues.txt"
    )
    path_report = os.path.join(
        path_gene, "report.pickle"
    )
    path_chromosome = os.path.join(
        path_gene, "chromosome.txt"
    )
    # Organization
    path_data_gene_signal_tissue = os.path.join(
        path_gene, "data_gene_signal_tissue.pickle"
    )

    # This data frame includes a gene's signals across persons (rows) and
    # major tissues (columns).
    # The organization segment of the distribution procedure has already
    # aggregated the signals of any minor tissues by mean.
    # The permutation procedure reads and permutes this data frame for each
    # gene.
    # The distribution procedure transforms this data frame further.
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
    path_data_gene_signals_tissues_persons_normal_standard = os.path.join(
        path_gene, "data_gene_signals_tissues_persons_normal_standard.pickle"
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
    utility.write_file_text_list(
        elements=[information["chromosome"]],
        delimiter="\n",
        path_file=path_chromosome
    )
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
    pandas.to_pickle(
        information["data_gene_signals_tissues_persons_normal_standard"],
        path_data_gene_signals_tissues_persons_normal_standard
    )

    # Modality
    with open(path_scores, "wb") as file_product:
        pickle.dump(information["scores"], file_product)

    pass


def write_product_collection(
    group=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        group (str): group of persons, either selection or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_data_signals_genes_persons = os.path.join(
        paths[group]["collection"], "data_signals_genes_persons.pickle"
    )
    path_genes_scores = os.path.join(
        paths[group]["collection"], "genes_scores.pickle"
    )
    path_scores = os.path.join(
        paths[group]["collection"], "scores.pickle"
    )
    path_data_report = os.path.join(
        paths[group]["collection"], "data_gene_report.pickle"
    )
    path_data_report_text = os.path.join(
        paths[group]["collection"], "data_gene_report.tsv"
    )
    path_data_signals_genes_persons_trait = os.path.join(
        paths[group]["trait"], "data_signals_genes_persons_trait.tsv"
    )

    # Write information to file.
    information["data_signals_genes_persons"].to_pickle(
        path=path_data_signals_genes_persons,
    )
    information["data_report"].to_pickle(
        path=path_data_report
    )
    information["data_report"].to_csv(
        path_or_buf=path_data_report_text,
        sep="\t",
        header=True,
        index=False,
    )
    with open(path_genes_scores, "wb") as file_product:
        pickle.dump(information["genes_scores"], file_product)
    with open(path_scores, "wb") as file_product:
        pickle.dump(information["scores"], file_product)
    information["data_signals_genes_persons_trait"].to_csv(
        path_or_buf=path_data_signals_genes_persons_trait,
        sep="\t",
        na_rep="NA",
        header=True,
        index=False,
    )

    pass


##########
# Collection


def read_collection_signals_gene_persons(
    group=None,
    gene=None,
    paths=None
):
    """
    Reads and organizes source information from file

    arguments:
        group (str): group of persons, either selection or ventilation
        gene (str): identifier of single gene for which to execute the process.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_gene = os.path.join(paths[group]["genes"], gene)
    path_data_gene_persons_signals = os.path.join(
        path_gene, "data_gene_persons_signals.pickle"
    )

    # Read information from file.
    data_gene_persons_signals = (
        pandas.read_pickle(path_data_gene_persons_signals)
    )
    # Compile and return information.
    return {
        "data_gene_persons_signals": data_gene_persons_signals,
    }


def read_collect_signals_genes_persons(
    group=None,
    genes=None,
    paths=None,
):
    """
    Collects genes' pantissue signals across persons.

    arguments:
        group (str): group of persons, either selection or ventilation
        genes (list<str>): identifiers of genes
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (object): Pandas data frame of genes' pantissue signals across persons

    """

    """
    Collects genes' aggregate, pantissue signals across persons.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        genes (list<str>): identifiers of genes

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' aggregate,
            pantissue signals across persons

    """

    # Collect information for genes.
    collection = dict()
    # Iterate on genes.
    for gene in genes:
        # Read source information from file.
        source_gene = read_collection_signals_gene_persons(
            group=group,
            gene=gene,
            paths=paths
        )
        # Organize data.
        collection[gene] = source_gene["data_gene_persons_signals"]
        pass
    # Return information.
    return collection


def read_collect_organize_signals_genes_persons(
    group=None,
    genes=None,
    data_persons_properties=None,
    paths=None,
):
    """
    Collects genes' pantissue signals across persons.

    arguments:
        group (str): group of persons, either selection or ventilation
        genes (list<str>): identifiers of genes
        data_persons_properties (object): Pandas data frame of persons and
            their properties
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (object): Pandas data frame of genes' pan-tissue signals across persons

    """

    # Organize template of data for persons.
    data_collection = data_persons_properties.copy(deep=True)
    data_collection.reset_index(
        level=None,
        inplace=True
    )
    data_collection = data_collection.loc[
        :, data_collection.columns.isin([
            "person", "sex_text"
        ])
    ]
    data_collection.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    #data_collection.reindex()
    data_collection.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    # Collect aggregate, pantissue signals across persons for all priority
    # genes.
    collection_signals = read_collect_signals_genes_persons(
        group=group,
        genes=genes,
        paths=paths,
    )

    # Collect information about genes.
    # Iterate on genes.
    for gene in collection_signals.keys():
        data_gene = collection_signals[gene]
        # Rename the aggregate, pantissue signals.
        data_gene.rename(
            columns={
                "value": gene,
            },
            inplace=True,
        )
        # Introduce aggregate, pantissue signals for each gene to person data.
        # Join
        data_collection = data_collection.join(
            data_gene,
            how="left",
            on="person"
        )
        pass

    # Organize data.
    data_collection.drop(
        labels="sex_text",
        axis="columns",
        inplace=True
    )
    data_collection.rename_axis(
        columns="gene",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Remove persons without signal for any genes.
    data_collection.dropna(
        axis="index",
        how="all",
        thresh=1,
        inplace=True,
    )

    # Report
    utility.print_terminal_partition(level=2)
    print("Here're genes' pantissue signals across persons...")
    print(data_collection)
    print("count of genes: " + str(data_collection.shape[1]))

    # Return information.
    return data_collection


def read_collect_genes_scores(
    group=None,
    genes=None,
    paths=None,
):
    """
    Collects information about genes.

    Data structure.

    - genes_scores (dict)
    -- gene (dict)
    --- dip (float)
    --- mixture (float)
    --- coefficient (float)

    - scores (dict)
    -- dip (dict)
    -- mixture (dict)
    -- coefficient (dict)

    arguments:
        group (str): group of persons, either selection or ventilation
        genes (list<str>): identifiers of genes
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (dict): information about genes' scores and permutations

    """

    # Collect genes' scores and permutations.
    utility.print_terminal_partition(level=3)
    print("Collecting genes' scores.")
    # Collect information about genes.
    genes_scores = dict()
    # Collect information about bimodality measures.
    scores = dict()
    measures = ["dip", "mixture", "coefficient"]
    for measure in measures:
        scores[measure] = dict()
        scores[measure]["values"] = list()
    # Iterate on genes.
    for gene in genes:
        # Specify directories and files.
        path_gene = os.path.join(paths[group]["genes"], gene)
        path_scores = os.path.join(
            path_gene, "scores.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            gene_scores = pickle.load(file_source)
        # Create entry for gene.
        # Compile information.
        genes_scores[gene] = gene_scores
        # Iterate on bimodality measures.
        for measure in measures:
            scores[measure]["values"].append(gene_scores[measure])

    # Compile information about scores.
    # Iterate on bimodality measures.
    for measure in measures:
        scores[measure]["mean"] = statistics.mean(scores[measure]["values"])
        scores[measure]["deviation"] = (
            statistics.stdev(scores[measure]["values"])
        )
    # Compile information.
    information = dict()
    information["genes_scores"] = genes_scores
    information["scores"] = scores
    # Return information.
    return information


def read_collect_gene_report(
    group=None,
    genes=None,
    paths=None,
):
    """
    Collects information about genes.

    arguments:
        group (str): group of persons, either selection or ventilation
        genes (list<str>): identifiers of genes
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:
        (object): Pandas data frame of information about genes' distributions

    """

    # Collect information about genes.
    records = list()
    # Iterate on directories for genes.
    for gene in genes:
        # Specify directories and files.
        path_gene = os.path.join(paths[group]["genes"], gene)
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
    data = data[[
        "identifier",
        "name",
        "chromosome",
        "start",
        "end",
        "persons_restriction",
        "persons_aggregation",
        "method",
        "count",
        "tissues",
        "tissues_mean",
        "tissues_median",
        "dip",
        "mixture",
        "coefficient",
    ]]
    data.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True
    )

    # Convert modality measures to standard, z-score space.
    # Calculate combination modality measures.
    if False:
        data_report = calculate_combination_scores(
            data_report=data,
        )
        data_report.sort_values(
            by=["combination"],
            axis="index",
            ascending=False,
            inplace=True,
        )
        #data_report.reset_index(
        #    level=None,
        #    inplace=True
        #)

    else:
        data_report = data

    utility.print_terminal_partition(level=2)
    print("Here's the data report...")
    print(data_report)
    print("count of genes: " + str(data_report.shape[0]))
    data_valid = data_report.loc[
        :, ["gene", "name", "coefficient", "dip", "mixture"]
    ]
    data_valid.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    print("count of genes with valid modalities: " + str(data_valid.shape[0]))
    utility.print_terminal_partition(level=2)

    return data_report


def organize_persons_genes_signals_quantitative_trait_loci(
    data_signals_genes_persons=None,
    data_gene_annotation=None,
):
    """
    Organizes genes' pan-tissue signals across persons in format for
    Quantitative Trait Loci (QTL) analysis in QTLtools.

    Format (text, tab and new line delimiters)

    #Chr  start    end      pid  gid        strand GTEX-111CU GTEX-111FC
    chr11 62337448 62393412 gene ENSG162174 +      -0.190985  -0.328112
    chr16 19856691 19886167 gene ENSG167191 +      -0.292748  -0.083921
    chr21 39867438 39929397 gene ENSG183036 -      -0.192331   0.886126
    chr6  10492223 10629368 gene ENSG111846 +      -0.543567   0.366667
    chr19 35248330 35267964 gene ENSG105699 -       0.152548  -0.275815


    arguments:
        data_sginals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' pan-tissue signals across persons
            in format for FastQTL
    """

    # Copy data.
    data_gene_annotation = data_gene_annotation.copy(deep=True)
    data_signals_genes_persons = data_signals_genes_persons.copy(deep=True)

    # Access identifiers of persons.
    persons = data_signals_genes_persons.index.to_numpy().tolist()

    # Transpose columns to rows and rows to columns.
    # Organize genes across rows and persons across columns.
    data_transposition = data_signals_genes_persons.transpose(copy=True)
    data_transposition.reset_index(
        level=None,
        inplace=True
    )
    # Introduce chromosomal coordinates.
    data_transposition["#Chr"] = data_transposition["gene"].apply(
        lambda gene:
            data_gene_annotation.loc[gene, "seqname"],
    )
    data_transposition["start"] = data_transposition["gene"].apply(
        lambda gene:
            data_gene_annotation.loc[gene, "start"],
    )
    data_transposition["end"] = data_transposition["gene"].apply(
        lambda gene:
            data_gene_annotation.loc[gene, "end"],
    )
    data_transposition["pid"] = "gene"
    data_transposition["strand"] = data_transposition["gene"].apply(
        lambda gene:
            data_gene_annotation.loc[gene, "strand"],
    )
    # Rename column.
    data_transposition.rename(
        columns={
            "gene": "gid",
        },
        inplace=True,
    )
    data_transposition.rename_axis(
        columns="",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Sort columns according to format requirements.
    columns = copy.deepcopy(persons)
    columns.insert(0, "strand")
    columns.insert(0, "gid")
    columns.insert(0, "pid")
    columns.insert(0, "end")
    columns.insert(0, "start")
    columns.insert(0, "#Chr")
    data_sort = data_transposition[[*columns]]

    # Return information.
    return data_sort


def read_collect_write_iterations(
    group=None,
    genes=None,
    data_gene_annotation=None,
    data_persons_properties=None,
    paths=None,
):
    """
    Reads, collects, and writes information from iterations.

    arguments:
        group (str): group of persons, either selection or ventilation
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report process.
    utility.print_terminal_partition(level=1)
    print(
        "Reading, collecting, and writing information about genes' " +
        "distributions."
    )
    utility.print_terminal_partition(level=2)

    # Check contents of directory.
    utility.print_terminal_partition(level=3)
    print("Check that directories exist for all genes to collect.")
    genes_distribution = utility.extract_subdirectory_names(
        path=paths[group]["genes"]
    )
    match_distribution = utility.compare_lists_by_inclusion(
        list_one=genes_distribution,
        list_two=genes
    )
    print(
        "Genes and distribution directories match: " +
        str(match_distribution)
    )
    utility.print_terminal_partition(level=2)

    # Genes' pantissue signals across persons.
    data_signals_genes_persons = read_collect_organize_signals_genes_persons(
        group=group,
        genes=genes_distribution,
        data_persons_properties=data_persons_properties,
        paths=paths,
    )

    # Genes' scores.
    genes_scores = read_collect_genes_scores(
        group=group,
        genes=genes_distribution,
        paths=paths,
    )
    # Report of genes' distributions.
    data_report = read_collect_gene_report(
        group=group,
        genes=genes_distribution,
        paths=paths,
    )

    # Format for Quantitative Trait Loci (QTL) in FastQTL
    # Extract distribution of gene's signals across persons for heritability
    # analysis.
    data_signals_genes_persons_trait = (
        organize_persons_genes_signals_quantitative_trait_loci(
            data_signals_genes_persons=data_signals_genes_persons,
            data_gene_annotation=data_gene_annotation,
    ))

    # Compile information.
    information = dict()
    information["data_signals_genes_persons"] = data_signals_genes_persons
    information["genes_scores"] = genes_scores["genes_scores"]
    information["scores"] = genes_scores["scores"]
    information["data_report"] = data_report
    information["data_signals_genes_persons_trait"] = (
        data_signals_genes_persons_trait
    )
    # Write product information to file.
    write_product_collection(
        group=group,
        information=information,
        paths=paths,
    )
    pass


###############################################################################
# Procedure


def execute_procedure(
    gene=None,
    group=None,
    data_gene_samples_signals=None,
    data_families_persons=None,
    data_gene_annotation=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of gene
        group (str): group of persons, either selection or ventilation
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        data_families_persons (object): Pandas data frame of person's
            identifiers and families' identifiers
        data_gene_annotation (object): Pandas data frame of genes' annotations
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Organization
    bin_organization = organize_data(
        data_gene_samples_signals=data_gene_samples_signals
    )

    # Distribution
    # The "aggregate_data" function prepares a data frame with the name
    # "data_gene_signals_tissues_persons_normal_standard".
    # This data frame gives gene's signals across tissues and persons after
    # normalization and standardization.
    bins = prepare_describe_distribution(
        data_gene_persons_tissues_signals=(
            bin_organization["data_gene_persons_tissues_signals"]
        ),
    )

    # Format for heritability in GCTA
    # Extract distribution of gene's signals across persons for heritability
    # analysis.
    data_gene_families_persons_signals = extract_gene_persons_signals(
        data_gene_persons_signals=(
            bins["bin_aggregation"]["data_gene_persons_signals"]
        ),
        data_families_persons=data_families_persons,
    )

    # Report
    report = prepare_gene_report(
        gene=gene,
        persons_restriction=(
            bins["bin_restriction"]["persons_count_restriction"]
        ),
        persons_aggregation=(
            bins["bin_aggregation"]["persons_count_aggregation"]
        ),
        method=bins["method"],
        count=bins["count"],
        tissues=len(bins["tissues"]),
        tissues_mean=bins["bin_restriction"]["tissues_count_mean"],
        tissues_median=bins["bin_restriction"]["tissues_count_median"],
        coefficient=bins["bin_modality"]["scores"]["coefficient"],
        dip=bins["bin_modality"]["scores"]["dip"],
        mixture=bins["bin_modality"]["scores"]["mixture"],
        data_gene_annotation=data_gene_annotation,
    )

    # Compile information.
    information = bin_organization.copy()
    information.update(bins["bin_restriction"])
    information.update(bins["bin_aggregation"])
    information.update(bins["bin_modality"])
    information["gene"] = gene
    information["tissues"] = bins["tissues"]
    information["data_gene_families_persons_signals"] = (
        data_gene_families_persons_signals
    )
    information["report"] = report
    information["chromosome"] = report["chromosome"].replace("chr", "")
    # Write product information to file.
    write_product_gene(
        gene=gene,
        group=group,
        information=information,
        paths=paths,
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
    # Initialize directories.
    paths = initialize_directories(dock=dock)

    # Execute procedure for each group of persons.
    execute_procedure_local_group(
        group="selection",
        paths=paths,
    )
    execute_procedure_local_group(
        group="ventilation",
        paths=paths,
    )
    # Report date and time.
    utility.print_terminal_partition(level=3)
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))
    utility.print_terminal_partition(level=3)

    pass


def execute_procedure_local_group(
    group=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        group (str): group of persons, either selection or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("... Distribution procedure for: " + str(group) + " persons...")

    # Read source information from file.
    source = read_source_initial(
        group=group,
        dock=paths["dock"],
    )
    print("count of genes: " + str(len(source["genes_selection"])))
    # Specify genes on which to iterate.
    #genes_iteration = random.sample(source["genes_selection"], 250)
    genes_iteration = source["genes_selection"]#[0:100]

    # Set up partial function for iterative execution.
    # Each iteration uses a different sequential value of the "gene" variable
    # with the same value of the "dock" variable.
    execute_procedure_gene = functools.partial(
        execute_procedure_local_sub,
        group=group,
        data_families_persons=source["data_families_persons"],
        data_gene_annotation=source["data_gene_annotation"],
        paths=paths,
    )
    # Initialize multiprocessing pool.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=8)
    # Iterate on genes.
    report = pool.map(
        execute_procedure_gene,
        genes_iteration,
    )
    # Pause procedure.
    time.sleep(10.0)
    # Read, collect, and write information from iterations.
    read_collect_write_iterations(
        group=group,
        genes=genes_iteration,
        data_persons_properties=source["data_persons_properties"],
        data_gene_annotation=source["data_gene_annotation"],
        paths=paths,
    )
    pass


def execute_procedure_local_sub(
    gene=None,
    group=None,
    data_families_persons=None,
    data_gene_annotation=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        group (str): group of persons, either selection or ventilation
        data_families_persons (object): Pandas data frame of persons'
            identifiers and families' identifiers
        data_gene_annotation (object): Pandas data frame of genes' annotations
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(
        gene=gene,
        group=group,
        dock=paths["dock"],
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        group=group,
        data_gene_samples_signals=source["data_gene_samples_signals"],
        data_families_persons=data_families_persons,
        data_gene_annotation=data_gene_annotation,
        paths=paths,
    )

    # Report progress.
    directories = os.listdir(paths[group]["genes"])
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
