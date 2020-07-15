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

import utility
import assembly

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
    paths["candidacy"] = os.path.join(paths["dock"], "candidacy")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["candidacy"])
    utility.create_directory(path=paths["candidacy"])
    # Define paths for cohorts of persons.
    cohorts = list()
    cohorts.append("selection")
    cohorts.append("respiration")
    cohorts.append("ventilation")
    for cohort in cohorts:
        paths[cohort] = dict()
        paths[cohort]["threshold"] = os.path.join(
            paths["candidacy"], cohort, "threshold"
        )
        # Define paths for groups of genes by their distributions.
        paths[cohort]["distribution"] = dict()
        paths[cohort]["distribution"]["multimodal"] = os.path.join(
            paths["candidacy"], cohort, "distribution", "multimodal"
        )
        paths[cohort]["distribution"]["unimodal"] = os.path.join(
            paths["candidacy"], cohort, "distribution", "unimodal"
        )
        paths[cohort]["distribution"]["nonmultimodal"] = os.path.join(
            paths["candidacy"], cohort, "distribution", "nonmultimodal"
        )
        paths[cohort]["distribution"]["any"] = os.path.join(
            paths["candidacy"], cohort, "distribution", "any"
        )
        # Initialize directories.
        utility.create_directories(path=paths[cohort]["threshold"])
        utility.create_directories(
            path=paths[cohort]["distribution"]["any"]
        )
        utility.create_directories(
            path=paths[cohort]["distribution"]["multimodal"]
        )
        utility.create_directories(
            path=paths[cohort]["distribution"]["unimodal"]
        )
        utility.create_directories(
            path=paths[cohort]["distribution"]["nonmultimodal"]
        )
    # Return information.
    return paths


def read_source(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_data_gene_annotation = os.path.join(
        dock, "selection", "tight", "gene_annotation",
        "data_gene_annotation_gencode.pickle"
    )
    path_genes_selection = os.path.join(
        dock, "selection", "tight", "samples_genes_signals",
        "genes.pickle"
    )

    # Use genes' bimodality measures from distribution procedure (not those
    # from the probability procedure) because here we want the raw values
    # before standardization.
    path_distribution_genes = os.path.join(
        dock, "distribution", cohort, "genes"
    )
    path_genes_scores = os.path.join(
        dock, "distribution", cohort, "collection", "genes_scores.pickle"
    )
    path_scores = os.path.join(
        dock, "distribution", cohort, "collection", "scores.pickle"
    )
    path_data_distribution_report = os.path.join(
        dock, "distribution", cohort, "collection", "data_gene_report.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )


    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_data_gene_annotation)
    with open(path_genes_selection, "rb") as file_source:
        genes_selection = pickle.load(file_source)
    genes_distribution = utility.extract_subdirectory_names(
        path=path_distribution_genes
    )
    with open(path_genes_scores, "rb") as file_source:
        genes_scores = pickle.load(file_source)
    with open(path_scores, "rb") as file_source:
        scores = pickle.load(file_source)
    data_distribution_report = pandas.read_pickle(
        path_data_distribution_report
    )
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes_selection": genes_selection,
        "genes_distribution": genes_distribution,
        "genes_scores": genes_scores,
        "scores": scores,
        "data_distribution_report": data_distribution_report,
        "data_signals_genes_persons": data_signals_genes_persons,
    }


def determine_selection_distribution_genes_valid_modalities(
    genes_selection=None,
    genes_distribution=None,
    data_distribution_report=None,
    report=None,
):
    """
    Determines selection genes with valid modalities from distribution
    procedure.

    Only genes with adequate signal coverage across tissues and persons have
    valid modalities from distribution procedure.

    arguments:
        genes_selection (list<str>): identifiers of genes from selection
        genes_distribution (list<str>): identifiers of genes from distribution
            procedure
        data_distribution_report (object): Pandas data frame of information
            about genes and their measures of modality
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information

    """

    # Organize data.
    genes_selection = copy.deepcopy(genes_selection)
    data_report = data_distribution_report.copy(deep=True)
    # Select genes with valid distribution modalities.
    data_valid = data_report.loc[
        :, data_report.columns.isin(["gene", "coefficient", "mixture", "dip"])
    ]
    data_valid.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    genes_distribution_valid = utility.collect_unique_elements(
        elements_original=data_valid.index.to_list()
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "count of all distribution genes: " +
            str(len(genes_distribution))
        )
        utility.print_terminal_partition(level=4)
        print(
            "count of all distribution genes with valid modalities: " +
            str(len(genes_distribution_valid))
        )
        utility.print_terminal_partition(level=2)
        pass
    # Select genes from selection procedure with valid distribution modalities.
    genes_selection_distribution_valid = utility.filter_common_elements(
        list_one=genes_selection,
        list_two=genes_distribution_valid,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "count of all selection genes: " +
            str(len(genes_selection))
        )
        utility.print_terminal_partition(level=4)
        print(
            "count of all selection genes with valid modalities: " +
            str(len(genes_selection_distribution_valid))
        )
        utility.print_terminal_partition(level=2)
        pass
    # Compile information.
    bin = dict()
    bin["data_distribution_report"] = data_valid
    bin["genes_selection_distribution"] = genes_selection_distribution_valid
    # Return information.
    return bin


def report_modality_measure_correlations(
    genes=None,
    measures=None,
    data_distribution_report=None,
):
    """
    Reports the pairwise correlation coefficients between measures of modality.

    arguments:
        genes (list<str>): identifiers of genes for which to report
            correlations
        measures (list<str>): measures of modality
        data_distribution_report (object): Pandas data frame of information
            about genes and their measures of modality

    raises:

    returns:

    """
    # Organize data.
    genes = copy.deepcopy(genes)
    data_report = data_distribution_report.copy(deep=True)
    data_report_genes = data_report.loc[
        data_report.index.isin(genes), :
    ]
    utility.print_terminal_partition(level=2)
    print(
        "Report modality measure correlations across genes: " +
        str(len(data_report_genes.index.to_list()))
    )
    # Calculate correlations between gene pairs of their pantissue signals
    # across persons.
    correlations = utility.collect_pairwise_correlations_pobabilities(
        method="spearman",
        features=measures,
        data=data_report_genes,
    )
    # Organize data matrix for correlation coefficients.
    data_correlations = utility.organize_symmetric_adjacency_matrix(
        features=measures,
        collection=correlations,
        key="correlation",
        threshold=False,
        fill=0.0,
    )
    utility.print_terminal_partition(level=2)
    print("Spearman correlation coefficients between pairs of measures.")
    utility.print_terminal_partition(level=3)
    print(data_correlations)

    # Organize data matrix for correlation probabilities.
    data_probabilities = utility.organize_symmetric_adjacency_matrix(
        features=measures,
        collection=correlations,
        key="probability",
        threshold=False,
        fill=0.0,
    )
    utility.print_terminal_partition(level=2)
    print("Probabilities of correlations between pairs of measures.")
    utility.print_terminal_partition(level=3)
    print(data_probabilities)

    pass


def select_genes_by_modality_measures_ranks(
    genes=None,
    proportion_least=None,
    proportion_greatest=None,
    measures=None,
    data_distribution_report=None,
    report=None,
):
    """
    Selects genes with least and greatest values of measures of modality.

    arguments:
        genes (list<str>): identifiers of genes for which to consider
            modalities
        proportion_least (float): proportion of genes to select from those with
            least values of modality measures
        proportion_greatest (float): proportion of genes to select from those
            with greatest values of modality measures
        measures (list<str>): measures of modality
        data_distribution_report (object): Pandas data frame of information
            about genes and their measures of modality
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information about selection of genes

    """

    # Organize data.
    genes = copy.deepcopy(genes)
    data_report = data_distribution_report.copy(deep=True)
    data_report_genes = data_report.loc[data_report.index.isin(genes), :]
    # Calculate count of genes to select from least and greatest extremes.
    count_total = len(genes)
    count_least = round(proportion_least * count_total)
    count_greatest = round(proportion_greatest * count_total)
    # Report.
    if report:
        print(
            "selection percentage least: " +
            str(round((proportion_least * 100), 2))
        )
        print("selection count least: " + str(count_least))
        utility.print_terminal_partition(level=3)
        print(
            "selection percentage greatest: " +
            str(round((proportion_greatest * 100), 2))
        )
        print("selection count greatest: " + str(count_greatest))
        pass
    # Iterate on measures of modality.
    bin = dict()
    for measure in measures:
        # Copy data.
        data_measure = data_report_genes.copy(deep=True)
        data_measure = data_measure.loc[:, ["name", measure]]
        # Sort by values of the measure.
        data_measure.sort_values(
            by=[measure],
            axis="index",
            ascending=True,
            inplace=True,
        )
        # Select least and greatest genes.
        # Pay attention to index values.
        # I validated the selection of threshold values.
        threshold_least = data_measure.iat[(count_least - 1), 1]
        data_least = data_measure.iloc[:count_least]
        genes_least = utility.collect_unique_elements(
            elements_original=data_least.index.to_list()
        )
        threshold_greatest = (
            data_measure.iat[(count_total - (count_greatest)), 1]
        )
        data_greatest = data_measure.iloc[(count_total - count_greatest):]
        genes_greatest = data_greatest.index.to_list()
        genes_greatest = utility.collect_unique_elements(
            elements_original=data_greatest.index.to_list()
        )
        # Collect information.
        bin[measure] = dict()
        bin[measure]["least"] = dict()
        bin[measure]["least"]["threshold"] = threshold_least
        bin[measure]["least"]["genes"] = genes_least
        bin[measure]["greatest"] = dict()
        bin[measure]["greatest"]["threshold"] = threshold_greatest
        bin[measure]["greatest"]["genes"] = genes_greatest
        pass
    # Organize measures' thresholds for plot.
    bin["measures_thresholds"] = dict()
    for measure in measures:
        bin["measures_thresholds"][measure] = (
            bin[measure]["greatest"]["threshold"]
        )
    # Return information.
    return bin


def validate_report_selection_thresholds(
    measures=None,
    selection=None,
    genes_scores=None,
):
    """
    Validates thresholds from selection of genes with least and greatest values
        of measures of modality.

    arguments:
        measures (list<str>): measures of modality
        selection (dict): selections of genes
        genes_scores (dict): information about genes' measures of modality

    raises:

    returns:
        (dict): information about selection of genes

    """

    utility.print_terminal_partition(level=2)
    print(
        "Validation of thresholds for selection of unimodal and multimodal " +
        "genes."
    )
    # Iterate on measures of modality.
    for measure in measures:
        for direction in ["least", "greatest"]:
            # Collect values of measure for selection of genes.
            values = list()
            for gene in selection[measure][direction]["genes"]:
                value = genes_scores[gene][measure]
                values.append(value)
            if direction == "least":
                validation = max(values)
            elif direction == "greatest":
                validation = min(values)
            selection[measure][direction]["threshold_validation"] = validation
            threshold = selection[measure][direction]["threshold"]
            utility.print_terminal_partition(level=3)
            print("measure: " + measure)
            print("direction: " + direction)
            print("threshold: " + str(round(threshold, 5)))
            print("validation: " + str(round(validation, 5)))
    # Return information.
    return selection


def calculate_report_threshold_deviation_factors(
    measures=None,
    selection=None,
    genes=None,
    data_distribution_report=None,
):
    """
    Calculates how many standard deviations by which each threshold diverges
    from mean.

    arguments:
        measures (list<str>): measures of modality
        selection (dict): selections of genes
        genes (list<str>): identifiers of genes for which to consider
            modalities
        data_distribution_report (object): Pandas data frame of information
            about genes and their measures of modality

    raises:

    returns:
        (dict): information about selection of genes

    """

    # Organize data.
    genes = copy.deepcopy(genes)
    data_report = data_distribution_report.copy(deep=True)
    data_report_measures = data_report.loc[
        :, data_report.columns.isin(["gene", "coefficient", "mixture", "dip"])
    ]
    data_report_genes = data_report_measures.loc[
        data_report_measures.index.isin(genes), :
    ]

    utility.print_terminal_partition(level=2)
    print(
        "By how many standard deviations do these thresholds differ from " +
        "means?"
    )
    # Iterate on measures of modality.
    factors = dict()
    for measure in measures:
        factors[measure] = dict()
        for direction in ["least", "greatest"]:
            threshold = selection[measure][direction]["threshold"]
            mean = numpy.nanmean(data_report_genes[measure].to_numpy())
            deviation = numpy.nanstd(
                data_report_genes[measure].to_numpy(),
                axis=0,
                ddof=1, # sample standard deviation
            )
            # Calculate factor.
            if direction == "least":
                # Least:
                # threshold = mean - (factor * deviation)
                # (factor * deviation) = (mean - threshold)
                # factor = (mean - threshold) / deviation
                factor = (mean - threshold) / deviation
                factors[measure][direction] = factor
            elif direction == "greatest":
                # Greatest:
                # threshold = mean + (factor * deviation)
                # (factor * deviation) = (threshold - mean)
                # factor = (threshold - mean) / deviation
                factor = (threshold - mean) / deviation
                factors[measure][direction] = factor
            utility.print_terminal_partition(level=3)
            print("measure: " + measure)
            print("direction: " + direction)
            print("mean: " + str(round(mean, 5)))
            print("deviation: " + str(round(deviation, 5)))
            print("threshold: " + str(round(threshold, 5)))
            print("factor: " + str(round(factor, 2)))
    # Return information.
    return factors


def extract_genes_modality_sets(
    direction=None,
    measures=None,
    selection=None,
):
    """
    Extracts identifiers of unique genes from selection by modality measures.

    arguments:
        direction (str): direction of distribution from which to select, lesser
            or greater
        measures (list<str>): measures of modality
        selection (dict): selections of genes

    raises:

    returns:
        (dict<list<str>>): identifiers of genes

    """

    # Organize sets of genes.
    sets = dict()
    for measure in measures:
        sets[measure] = selection[measure][direction]["genes"]

    # Select genes that pass filters by multiple measures of bimodality.
    genes_1 = utility.select_elements_by_sets(
        names=measures,
        sets=sets,
        count=1,
    )
    genes_2 = utility.select_elements_by_sets(
        names=measures,
        sets=sets,
        count=2,
    )
    genes_3 = utility.select_elements_by_sets(
        names=measures,
        sets=sets,
        count=3,
    )

    # Summarize information.
    utility.print_terminal_partition(level=2)
    print("Selection of genes by: " + direction)
    print("... any 1 sets: " + str(len(genes_1)))
    print("... any 2 sets: " + str(len(genes_2)))
    print("... any 3 sets: " + str(len(genes_3)))

    # Collect information.
    bin = dict()
    bin["measures_1"] = genes_1
    bin["measures_2"] = genes_2
    bin["measures_3"] = genes_3
    bin["sets_genes_measures"] = sets
    # Return information.
    return bin


def prepare_gene_report(
    genes=None,
    direction=None,
    measures=None,
    selection=None,
    data_gene_annotation=None,
):
    """
    Prepares a report of genes that have unimodal or multimodal distributions.

    arguments:
        genes (list<str>): identifiers of genes
        direction (str): direction of distribution from which to select, lesser
            or greater
        measures (list<str>): measures of modality
        selection (dict): selections of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of information about genes

    """

    # Iterate on genes.
    # Collect information.
    records = list()
    for gene in genes:
        record = dict()

        # Organize gene's properties.
        record["identifier"] = gene
        record["name"] = assembly.access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
        record["chromosome"] = data_gene_annotation.loc[gene, "seqname"]
        record["start"] = data_gene_annotation.loc[gene, "start"]
        record["end"] = data_gene_annotation.loc[gene, "end"]

        # Iterate on measures of modality.
        # Collect information about each measure.
        matches = 0
        for measure in measures:
            # Determine whether gene passes threshold by measure of modality.
            name = (measure + "_pass")
            match = (gene in selection[measure][direction]["genes"])
            record[name] = match
            if match:
                matches += 1
        record["passes"] = matches
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
        "dip_pass",
        "mixture_pass",
        "coefficient_pass",
        "passes",
    ]]
    # Sort by values of the measure.
    data.sort_values(
        by=["passes"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    data.reset_index(
        drop=True,
        inplace=True,
    )

    # Return information.
    return data


# Previous method
# Selection of genes strictly by mean and standard deviations


def determine_measures_thresholds(
    measures=None,
    scores=None,
    factors=None,
):
    """
    Determines values of thresholds for genes' measures of bimodality.

    arguments:
        measures (list<str>): measures of bimodality
        scores (dict<dict>): information about genes' measures of bimodality
        factors (dict<float>): factors for each measure of bimodality

    raises:

    returns:
        (dict<float>): values of thresholds for genes' measures of bimodality

    """

    # Collect thresholds for measures of bimodality.
    thresholds = dict()
    thresholds["lesser"] = dict()
    thresholds["greater"] = dict()
    for measure in measures:
        # Calculate selection thresholds.
        thresholds["lesser"][measure] = (
            scores[measure]["mean"] -
            (factors[measure] * scores[measure]["deviation"])
        )
        thresholds["greater"][measure] = (
            scores[measure]["mean"] +
            (factors[measure] * scores[measure]["deviation"])
        )
    return thresholds


def summarize_measures_thresholds(
    measures=None,
    scores=None,
    thresholds=None,
):
    """
    Summarizes values of thresholds for genes' measures of bimodality.

    arguments:
        measures (list<str>): measures of bimodality
        scores (dict<dict>): information about genes' measures of bimodality
        thresholds (dict<float>): values of thresholds for genes' measures of
            bimodality

    raises:

    returns:


    """

    utility.print_terminal_partition(level=2)
    for measure in measures:
        utility.print_terminal_partition(level=3)
        print("measure: " + str(measure))
        print("mean: " + str(scores[measure]["mean"]))
        print("deviation: " + str(scores[measure]["deviation"]))
        print("threshold lesser: " + str(thresholds["lesser"][measure]))
        print("threshold greater: " + str(thresholds["greater"][measure]))
    pass
    utility.print_terminal_partition(level=2)


def copy_split_minimal_gene_data(
    measure=None,
    data_genes_distributions=None,
):
    """
    Copy and split information about genes.

    arguments:
        measure (str): name of a measure of bimodality
        data_genes_distributions (object): Pandas data frame of information
            about genes and their measures of bimodality

    raises:

    returns:
        (object): Pandas data frame of genes' identifiers and values of a
            measure of bimodality

    """

    # Select relevant information.
    columns = list()
    columns.append("identifier")
    columns.append(measure)
    data_copy = data_genes_distributions.copy(deep=True)
    data = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    # Return information.
    return data


def filter_genes_by_bimodality_threshold(
    data=None,
    measure=None,
    threshold=None,
    direction=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality measures with permutation probabilities below threshold.

    arguments:
        data (object): Pandas data frame of information about genes and their
            measures of bimodality
        measure (str): name of a measure of bimodality
        threshold (float): value of a threshold for a measure of bimodality
        direction (str): direction of distribution from which to select, lesser
            or greater

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities

    """

    # Copy data.
    data = data.copy(deep=True)
    # Determine whether values pass threshold.
    # Only consider the original modality measures for this threshold.
    if direction == "greater":
        data_threshold = data.loc[(data[measure] > threshold), :]
    elif direction == "lesser":
        data_threshold = data.loc[(data[measure] < threshold), :]
    # Return information.
    return data_threshold


def filter_genes_probabilities_threshold_old(
    data=None,
    threshold=None,
    count=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality measures with permutation probabilities below threshold.

    arguments:
        data (object): Pandas data frame of probabilities from genes' scores
            and permutations
        threshold (float): maximal probability
        count (int): minimal count of probabilities that must pass threshold

    raises:

    returns:
        (object): Pandas data frame of probabilities from genes' scores and
            distributions

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Count how many of the gene's modality probabilities are below threshold
    # Keep genes with at least 2 measures below threshold

    # Determine whether values pass threshold.
    # Only consider the original modality measures for this threshold.
    data_selection = data.loc[ :, ["coefficient", "dip", "mixture"]]
    data_threshold = (data_selection <= threshold)
    # This aggregation operation produces a series.
    # Aggregate columns to determine which rows to keep in the next step.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=2),
        axis="columns",
    )

    # Select rows and columns with appropriate values.
    data_pass = data.loc[data_count, : ]
    # Return information.
    return data_pass


def filter_genes_by_bimodality_thresholds(
    measures=None,
    thresholds=None,
    data_genes_distributions=None,
    direction=None,
):
    """
    Copy and split information about genes.

    arguments:
        measures (list<str>): measures of bimodality
        thresholds (dict<float>): values of thresholds for measures of
            bimodality
        data_genes_distributions (object): Pandas data frame of information
            about genes and their measures of bimodality
        direction (str): direction of distribution from which to select, lesser
            or greater

    raises:

    returns:
        (dict<list<str>>): identifiers of genes that pass filtration by
            thresholds on each measure of bimodality

    """

    utility.print_terminal_partition(level=1)
    print(
        "count of genes filtered by probabilities of each bimodality " +
        "measurement"
    )

    # Collect genes from filtration by each measurement of bimodality.
    entries = dict()
    for measure in measures:
        # Copy minimal genes' data for each measure of bimodality.
        data_measure = copy_split_minimal_gene_data(
            measure=measure,
            data_genes_distributions=data_genes_distributions,
        )

        # Filter genes by threshold on each measure's probabilities.
        data_filter = filter_genes_by_bimodality_threshold(
            data=data_measure,
            measure=measure,
            threshold=thresholds[direction][measure],
            direction=direction,
        )

        # Extract genes' identifiers.
        genes = data_filter["identifier"].tolist()
        utility.print_terminal_partition(level=3)
        print(measure + ": " + str(len(genes)))

        # Compile information.
        entries[measure] = genes

    # Return information.
    return entries


def old_execution_method():
    """
    Determines values of thresholds for genes' measures of bimodality.

    raises:

    returns:

    """

    # Set measures of modality.
    measures = ["dip", "mixture", "coefficient"]
    # Set factors for each measure of bimodality.
    # Set factors to select 200-300 genes by each measure of bimodality.
    # Values of bimodality coefficient greater than 0.555 are multimodal.
    factors = dict()
    factors["dip"] = 2 # (mean + (2 * sigma) = 0.017): 169
    factors["mixture"] = 3 # (mean + (3 * sigma) = 0.489): 213
    factors["coefficient"] = 3 # (mean + (3 * sigma) = 0.556): 226
    # Calculate thresholds for each measure of bimodality.
    thresholds = determine_measures_thresholds(
        measures=measures,
        scores=source["scores"],
        factors=factors,
    )
    utility.print_terminal_partition(level=2)
    print("thresholds for unimodal and multimodal genes")
    summarize_measures_thresholds(
        measures=measures,
        scores=source["scores"],
        thresholds=thresholds,
    )

    # Filter genes by thresholds on measures of bimodality.
    genes_measures_unimodal = filter_genes_by_bimodality_thresholds(
        measures=measures,
        thresholds=thresholds,
        data_genes_distributions=source["data_distribution_report"],
        direction="lesser",
    )
    genes_measures_multimodal = filter_genes_by_bimodality_thresholds(
        measures=measures,
        thresholds=thresholds,
        data_genes_distributions=source["data_distribution_report"],
        direction="greater",
    )

    #############################
    # TODO:
    # select nonunimodal genes by 1, 2, or 3 measures
    # select unimodal genes

    ############################

    utility.print_terminal_partition(level=2)
    print(
        "genes_unimodal_1: " + str(len(genes_unimodal_1))
    )
    print(
        "genes_unimodal_2: " + str(len(genes_unimodal_2))
    )
    print(
        "genes_unimodal_3: " + str(len(genes_unimodal_3))
    )
    print(
        "genes_multimodal_1: " + str(len(genes_multimodal_1))
    )
    print(
        "genes_multimodal_2: " + str(len(genes_multimodal_2))
    )
    print(
        "genes_multimodal_3: " + str(len(genes_multimodal_3))
    )

    pass


##########
# Product


def write_product(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file.
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_measures_thresholds = os.path.join(
        paths[cohort]["threshold"], "measures_thresholds.pickle"
    )
    path_sets_unimodal = os.path.join(
        paths[cohort]["distribution"]["unimodal"], "sets_unimodal.pickle"
    )
    path_sets_multimodal = os.path.join(
        paths[cohort]["distribution"]["multimodal"], "sets_multimodal.pickle"
    )

    path_genes_any = os.path.join(
        paths[cohort]["distribution"]["any"], "genes.pickle"
    )
    path_genes_any_text = os.path.join(
        paths[cohort]["distribution"]["any"], "genes.txt"
    )
    path_genes_multimodal = os.path.join(
        paths[cohort]["distribution"]["multimodal"], "genes.pickle"
    )
    path_genes_multimodal_text = os.path.join(
        paths[cohort]["distribution"]["multimodal"], "genes.txt"
    )
    path_genes_nonmultimodal = os.path.join(
        paths[cohort]["distribution"]["nonmultimodal"], "genes.pickle"
    )
    path_genes_nonmultimodal_text = os.path.join(
        paths[cohort]["distribution"]["nonmultimodal"], "genes.txt"
    )
    path_genes_unimodal = os.path.join(
        paths[cohort]["distribution"]["unimodal"], "genes.pickle"
    )
    path_genes_unimodal_text = os.path.join(
        paths[cohort]["distribution"]["unimodal"], "genes.txt"
    )

    path_data_genes_unimodal = os.path.join(
        paths[cohort]["distribution"]["unimodal"], "data_genes_unimodal.pickle"
    )
    path_data_genes_unimodal_text = os.path.join(
        paths[cohort]["distribution"]["unimodal"], "data_genes_unimodal.tsv"
    )
    path_data_genes_multimodal = os.path.join(
        paths[cohort]["distribution"]["multimodal"],
        "data_genes_multimodal.pickle"
    )
    path_data_genes_multimodal_text = os.path.join(
        paths[cohort]["distribution"]["multimodal"],
        "data_genes_multimodal.tsv"
    )

    # Write information to file.
    with open(path_measures_thresholds, "wb") as file_product:
        pickle.dump(information["measures_thresholds"], file_product)
    with open(path_sets_unimodal, "wb") as file_product:
        pickle.dump(information["sets_unimodal"], file_product)
    with open(path_sets_multimodal, "wb") as file_product:
        pickle.dump(information["sets_multimodal"], file_product)

    with open(path_genes_any, "wb") as file_product:
        pickle.dump(
            information["genes_any"], file_product
        )
    utility.write_file_text_list(
        elements=information["genes_any"],
        delimiter="\n",
        path_file=path_genes_any_text
    )
    with open(path_genes_multimodal, "wb") as file_product:
        pickle.dump(
            information["genes_multimodal"], file_product
        )
    utility.write_file_text_list(
        elements=information["genes_multimodal"],
        delimiter="\n",
        path_file=path_genes_multimodal_text
    )
    with open(path_genes_nonmultimodal, "wb") as file_product:
        pickle.dump(
            information["genes_nonmultimodal"], file_product
        )
    utility.write_file_text_list(
        elements=information["genes_nonmultimodal"],
        delimiter="\n",
        path_file=path_genes_nonmultimodal_text
    )
    with open(path_genes_unimodal, "wb") as file_product:
        pickle.dump(
            information["genes_unimodal"], file_product
        )
    utility.write_file_text_list(
        elements=information["genes_unimodal"],
        delimiter="\n",
        path_file=path_genes_unimodal_text
    )

    information["data_genes_unimodal"].to_pickle(
        path=path_data_genes_unimodal
    )
    information["data_genes_unimodal"].to_csv(
        path_or_buf=path_data_genes_unimodal_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_genes_multimodal"].to_pickle(
        path=path_data_genes_multimodal
    )
    information["data_genes_multimodal"].to_csv(
        path_or_buf=path_data_genes_multimodal_text,
        sep="\t",
        header=True,
        index=True,
    )

    pass


###############################################################################
# Procedure


def select_report_write_candidate_modality_genes_cohort(
    cohort=None,
    paths=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=1)
    print("... Candidacy procedure for: " + str(cohort) + " persons...")
    utility.print_terminal_partition(level=2)

    # Read source information from file.
    source = read_source(
        cohort=cohort,
        dock=paths["dock"],
    )
    # Determine genes with valid modalities from distribution procedure.
    bin_valid = (
        determine_selection_distribution_genes_valid_modalities(
            genes_selection=source["genes_selection"],
            genes_distribution=source["genes_distribution"],
            data_distribution_report=source["data_distribution_report"],
            report=False,
    ))

    # Set measures of modality.
    measures = list(source["scores"].keys())
    # Describe correlations between pairs of modality measures.
    if False:
        report_modality_measure_correlations(
            genes=bin_valid["genes_selection_distribution"],
            measures=measures,
            data_distribution_report=source["data_distribution_report"],
        )
        pass

    # Specify thresholds for each cohort.
    # Aim for <= 785 genes in each group.
    if cohort == "selection":
        proportion_least = 0.2915 # 784
        proportion_greatest = 0.02155 # 785
    elif cohort == "respiration":
        proportion_least = 0.3065 # 785
        proportion_greatest = 0.0213825 # 784
    elif cohort == "ventilation":
        proportion_least = 0.2923 # 784
        proportion_greatest = 0.022 # 785

    # Select genes with least and greatest values of each measure of modality.
    # Use less stringency for unimodal genes in order to select those that are
    # not multimodal by any measures.
    # Select same count of unimodal genes and multimodal genes.
    utility.print_terminal_partition(level=2)
    selection = select_genes_by_modality_measures_ranks(
        genes=bin_valid["genes_selection_distribution"],
        proportion_least=proportion_least,
        proportion_greatest=proportion_greatest,
        measures=measures,
        data_distribution_report=source["data_distribution_report"],
        report=False,
    )
    # Validate genes and thresholds.
    if True:
        validation = validate_report_selection_thresholds(
            measures=measures,
            selection=selection,
            genes_scores=source["genes_scores"],
        )
    # Calculate how many standard deviations by which each threshold diverges
    # from mean.
    if True:
        factors = calculate_report_threshold_deviation_factors(
            measures=measures,
            selection=selection,
            genes=bin_valid["genes_selection_distribution"],
            data_distribution_report=source["data_distribution_report"],
        )

    # Extract genes that pass thresholds by measures of modality.
    utility.print_terminal_partition(level=2)
    print(
        "Extract identifiers of genes by counts of thresholds they pass."
    )
    bin_genes_unimodal = extract_genes_modality_sets(
        direction="least",
        measures=measures,
        selection=selection,
    )
    bin_genes_multimodal = extract_genes_modality_sets(
        direction="greatest",
        measures=measures,
        selection=selection,
    )

    # Determine all genes that are not in the multimodality set.
    genes_nonmultimodal = utility.filter_unique_exclusion_elements(
        elements_exclusion=bin_genes_multimodal["measures_1"],
        elements_total=bin_valid["genes_selection_distribution"],
    )
    utility.print_terminal_partition(level=2)
    print(
        "Count of all selection genes with valid modalities: " +
        str(len(bin_valid["genes_selection_distribution"]))
    )
    print(
        "Count of multimodal genes: " +
        str(len(bin_genes_multimodal["measures_1"]))
    )
    print(
        "Count of unimodal genes: " +
        str(len(bin_genes_unimodal["measures_3"]))
    )
    print(
        "Count of all genes not multimodal (other): " +
        str(len(genes_nonmultimodal))
    )

    # Rank genes by the counts of measures by which they pass thresholds.
    # Unimodal genes: Consider genes that have extreme least values by all
    # measures of modality.
    # Multimodal genes: Consider genes that have extreme greatest values by any
    # measure of modality.
    data_genes_unimodal = prepare_gene_report(
        genes=bin_genes_unimodal["measures_3"],
        direction="least",
        measures=measures,
        selection=selection,
        data_gene_annotation=source["data_gene_annotation"],
    )
    data_genes_multimodal = prepare_gene_report(
        genes=bin_genes_multimodal["measures_1"],
        direction="greatest",
        measures=measures,
        selection=selection,
        data_gene_annotation=source["data_gene_annotation"],
    )
    if False:
        utility.print_terminal_partition(level=2)
        print("Summary of unimodal genes.")
        print(data_genes_unimodal)
        utility.print_terminal_partition(level=2)
        print("Summary of multimodal genes.")
        print(data_genes_multimodal)

    # Compile information.
    information = {
        "genes_any": bin_valid["genes_selection_distribution"],
        "genes_unimodal": bin_genes_unimodal["measures_3"],
        "genes_multimodal": bin_genes_multimodal["measures_1"],
        "genes_nonmultimodal": genes_nonmultimodal,
        "data_genes_unimodal": data_genes_unimodal,
        "data_genes_multimodal": data_genes_multimodal,
        "measures_thresholds": selection["measures_thresholds"],
        "sets_unimodal": bin_genes_unimodal["sets_genes_measures"],
        "sets_multimodal": bin_genes_multimodal["sets_genes_measures"],
    }
    # Write product information to file.
    write_product(
        cohort=cohort,
        information=information,
        paths=paths,
    )
    pass


def execute_procedure(
    dock=None
):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    paths = initialize_directories(dock=dock)

    # Execute procedure for each cohort of persons.
    cohorts = [
        "selection",
        "respiration",
        "ventilation",
    ]
    for cohort in cohorts:
        select_report_write_candidate_modality_genes_cohort(
            cohort=cohort,
            paths=paths,
        )
    pass


if (__name__ == "__main__"):
    execute_procedure_local()
