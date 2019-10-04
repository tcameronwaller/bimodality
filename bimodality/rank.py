"""
...

"""

###############################################################################
# Notes

# If the count of bootstrap shuffles is too small, then it is likely not to
# have any values equal to or greater than the real value, in which case the
# probability is zero.

###############################################################################
# Installation and importation

# Standard

import os
import math
import statistics
import pickle
import copy
import random
import itertools

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import distribution
import plot
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
    path_genes = os.path.join(path_split, "genes.txt")
    path_selection = os.path.join(dock, "selection")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation.pickle"
    )
    path_distribution = os.path.join(dock, "distribution")
    path_distribution_report = os.path.join(
        path_distribution, "data_gene_report.pickle"
    )
    path_permutation = os.path.join(dock, "permutation")
    path_collection = os.path.join(dock, "collection")
    path_scores = os.path.join(
        path_collection, "genes_scores.pickle"
    )
    path_permutations = os.path.join(
        path_collection, "genes_permutations.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes,
    )
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_gene_distribution_report = pandas.read_pickle(
        path_distribution_report
    )
    with open(path_scores, "rb") as file_source:
        genes_scores = pickle.load(file_source)
    with open(path_permutations, "rb") as file_source:
        genes_permutations = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "data_gene_annotation": data_gene_annotation,
        "data_gene_distribution_report": data_gene_distribution_report,
        "genes_scores": genes_scores,
        "genes_permutations": genes_permutations,
    }


# Probability


def calculate_probability_equal_greater(
    value=None,
    distribution=None,
):
    """
    Calculates from a distribution the probability of obtaining a value equal
    to or greater than a specific value.

    arguments:
        value (float): value
        distribution (list<float>): distribution of values

    raises:

    returns:
        (float): probability of obtaining from the distribution a value equal
            to or greater than the specific value

    """

    count_total = len(distribution)
    count_matches = 0
    for value_distribution in distribution:
        if value_distribution >= value:
            count_matches += 1
    probability = count_matches / count_total
    return probability


def calculate_probabilities_genes(
    genes_scores=None,
    genes_permutations=None,
):
    """
    Calculates probabilities (p-values) from genes' scores and random
    distributions.

    Data structure.
    - genes_p_values (dict)
    -- gene (dict)
    ---- coefficient (float)
    ---- dip (float)
    ---- mixture (float)

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations

    raises:

    returns:
        (dict<dict<float>>): information about genes' probabilities

    """

    # Collect information about genes.
    entries = dict()
    # Iterate on genes.
    for gene in genes_scores:
        # Access information about gene's scores and distributions.
        scores = genes_scores[gene]
        permutations = genes_permutations[gene]
        # Scores
        dip = scores["dip"]
        mixture = scores["mixture"]
        coefficient = scores["coefficient"]
        combination = scores["combination"]
        # Permutations
        dips = permutations["dip"].tolist()
        mixtures = permutations["mixture"].tolist()
        coefficients = permutations["coefficient"].tolist()
        combinations = permutations["combination"].tolist()
        # Calculate p-values.
        # These scores of bimodality indicate greater bimodality as values
        # increase.
        # The scores are unidirectional, so the hypothesis is unidirectional.
        # The p-value is the probability of obtaining by random chance a value
        # equal to or greater than the actual score.
        probability_dip = calculate_probability_equal_greater(
            value=dip,
            distribution=dips,
        )
        probability_mixture = calculate_probability_equal_greater(
            value=mixture,
            distribution=mixtures,
        )
        probability_coefficient = calculate_probability_equal_greater(
            value=coefficient,
            distribution=coefficients,
        )
        probability_combination = calculate_probability_equal_greater(
            value=combination,
            distribution=combinations,
        )
        # Calculate combination of p-values.
        if False:
            probability_combination = scipy.stats.combine_pvalues(
                [probability_coefficient, probability_dip, probability_mixture],
                method="stouffer",
                weights=None,
            )[1]
        # Compile information.
        entry = dict()
        entry["gene"] = gene
        entry["dip"] = probability_dip
        entry["mixture"] = probability_mixture
        entry["coefficient"] = probability_coefficient
        entry["combination"] = probability_combination
        entries[gene] = entry

    # Return information.
    return entries


# Summary


def organize_genes_scores_probabilities(
    genes_scores=None,
    genes_permutations=None,
    genes_probabilities=None,
    data_gene_distribution_report=None,
    data_gene_annotation=None,
):
    """
    Collects tissues, patients, and genes' signals for each sample.

    Data structure
    identifier        name    coefficient  p_coefficient  dip     p_dip   ...
     ENSG00000078808   TK421   0.75         0.003          0.12    0.005  ...
     ENSG00000029363   R2D2    0.55         0.975          0.23    0.732  ...

    arguments:
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_permutations (dict<dict<list<float>>>): information about genes'
            permutations
        genes_probabilities (dict<dict<float>>): information about genes'
            probabilities
        data_gene_distribution_report (object): Pandas data frame of
            information about genes' distributions
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality scores,
            and probabilities

    """

    # Organize data.
    data_gene_distribution_report.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )

    # Collect information about genes.
    records = list()
    # Iterate on genes.
    for gene in genes_scores:
        # Access information about gene.

        name = data_gene_annotation.loc[gene, "gene_name"]
        chromosome = data_gene_annotation.loc[gene, "seqname"]
        start = data_gene_annotation.loc[gene, "start"]
        end = data_gene_annotation.loc[gene, "end"]
        persons = (
            data_gene_distribution_report.loc[gene, "persons_aggregation"]
        )
        tissues = data_gene_distribution_report.loc[gene, "tissues"]
        method = data_gene_distribution_report.loc[gene, "method"]
        count = data_gene_distribution_report.loc[gene, "count"]
        tissues_mean = data_gene_distribution_report.loc[gene, "tissues_mean"]
        tissues_median = (
            data_gene_distribution_report.loc[gene, "tissues_median"]
        )

        scores = genes_scores[gene]
        score_dip = scores["dip"]
        score_mixture = scores["mixture"]
        score_coefficient = scores["coefficient"]
        score_combination = scores["combination"]

        probabilities = genes_probabilities[gene]
        probability_dip = probabilities["dip"]
        probability_mixture = probabilities["mixture"]
        probability_coefficient = probabilities["coefficient"]
        probability_combination = probabilities["combination"]

        # Create record for gene.
        record = dict()
        # Compile information.

        record["identifier"] = gene
        record["name"] = name
        record["chromosome"] = chromosome
        record["start"] = start
        record["end"] = end
        record["persons"] = persons
        record["tissues"] = tissues
        record["method"] = method
        record["count"] = count
        record["tissues_mean"] = tissues_mean
        record["tissues_median"] = tissues_median

        record["dip"] = score_dip
        record["dip_probability"] = probability_dip
        record["mixture"] = score_mixture
        record["mixture_probability"] = probability_mixture
        record["coefficient"] = score_coefficient
        record["coefficient_probability"] = probability_coefficient
        record["combination"] = score_combination
        record["combination_probability"] = probability_combination

        records.append(record)

    # Organize data.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    data = data[[
        "identifier",
        "name",
        "chromosome",
        "start",
        "end",
        "persons",
        "tissues",
        "method",
        "count",
        "tissues_mean",
        "tissues_median",
        "dip",
        "dip_probability",
        "mixture",
        "mixture_probability",
        "coefficient",
        "coefficient_probability",
        "combination",
        "combination_probability",
    ]]
    data["dip"] = data["dip"].astype("float32")
    data["dip_probability"] = data["dip_probability"].astype("float32")
    data["mixture"] = data["mixture"].astype("float32")
    data["mixture_probability"] = data["mixture_probability"].astype("float32")
    data["coefficient"] = data["coefficient"].astype("float32")
    data["coefficient_probability"] = (
        data["coefficient_probability"].astype("float32")
    )
    data["combination"] = data["combination"].astype("float32")
    data["combination_probability"] = (
        data["combination_probability"].astype("float32")
    )
    data.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="identifier",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="properties",
        axis="columns",
        copy=False,
        inplace=True
    )
    data.sort_values(
        by=["combination"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Return information.
    return data


# Filtration


def copy_split_minimal_gene_data(
    measurement=None,
    probability=None,
    data_genes_scores_probabilities=None,
):
    """
    Copy and split information about genes.

    arguments:
        measurement (str): name of a measurement of bimodality
        probability (str): name of a measurement's probability
        data_genes_scores_probabilities (object): Pandas data frame of genes'
            properties, bimodality scores, and probabilities

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities

    """

    # Select relevant information.
    columns = list()
    columns.append(measurement)
    columns.append(probability)
    data_copy = data_genes_scores_probabilities.copy(deep=True)
    data = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    # Return information.
    return data


def filter_genes_probabilities_threshold(
    data=None,
    probability=None,
    threshold=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality measures with permutation probabilities below threshold.

    arguments:
        data (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities
        probability (str): name of a measurement's probability
        threshold (float): maximal probability

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities

    """

    # Determine whether values pass threshold.
    # Only consider the original modality measures for this threshold.
    data_threshold = data.loc[(data[probability] < threshold), :]
    # Return information.
    return data_threshold


# TODO: this function is useful...
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


def filter_genes_by_probabilities(
    measurements=None,
    threshold=None,
    data_genes_scores_probabilities=None,
):
    """
    Copy and split information about genes.

    arguments:
        measurements (list<str>): names of measurements of bimodality
        threshold (float): maximal probability
        data_genes_scores_probabilities (object): Pandas data frame of genes'
            properties, bimodality scores, and probabilities

    raises:

    returns:
        (dict<list<str>>): identifiers of genes that pass filtration by
            probabilities of each bimodality measurement

    """

    utility.print_terminal_partition(level=1)
    print(
        "count of genes filtered by probabilities of each bimodality " +
        "measurement"
    )

    # Collect genes from filtration by each measurement of bimodality.
    entries = dict()
    for measurement in measurements:
        probability = (measurement + "_probability")
        # Copy minimal genes' data for each modality measure.
        data_measurement = copy_split_minimal_gene_data(
            measurement=measurement,
            probability=probability,
            data_genes_scores_probabilities=data_genes_scores_probabilities,
        )

        # Filter genes by threshold on each measure's probabilities.
        data_filter = filter_genes_probabilities_threshold(
            data=data_measurement,
            probability=probability,
            threshold=threshold,
        )

        # Extract genes' identifiers.
        genes = data_filter.index.tolist()
        utility.print_terminal_partition(level=3)
        print(measurement + ": " + str(len(genes)))

        # Compile information.
        entries[measurement] = genes

    # Return information.
    return entries


# Rank


def select_rank_genes(
    data_genes_scores_probabilities=None,
    selection=None,
    rank=None,
):
    """
    Prepares ranks of genes for subsequent functional analyses.

    arguments:
        data_genes_scores_probabilities (object): Pandas data frame of genes'
            properties, bimodality scores, and probabilities
        selection (list<str>): genes for selection
        rank (str): property to use for ranks

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality scores,
            and probabilities

    """

    # Copy data.
    data_copy = data_genes_scores_probabilities.copy(deep=True)
    # Select data for genes of interest.
    data = data_copy.loc[
        data_copy.index.isin(selection), :
    ]
    # Rank genes.
    data.sort_values(
        by=[rank],
        axis="index",
        ascending=False,
        inplace=True,
    )
    # Organize information.
    # Return information.
    return data


# Exportation


def export_genes_function(
    data_genes=None,
    identifier=None,
    rank=None,
):
    """
    Prepares table of genes for export to functional analysis.

    arguments:
        data_genes (object): Pandas data frame of genes' properties, bimodality
            scores, and probabilities
        identifier (str): identifier to use, "ensembl" or "hugo"
        rank (bool): whether to include gene rankings in the table

    raises:

    returns:
        (object): Pandas data frame of genes' properties, bimodality scores,
            and probabilities

    """

    # Copy data.
    data_copy = data_genes.copy(deep=True)
    # Organize data.
    data_copy.reset_index(
        level=None,
        inplace=True
    )
    # Select columns.
    columns = list()
    if identifier == "ensembl":
        columns.append("identifier")
    elif identifier == "hugo":
        columns.append("name")
    if rank:
        columns.append("combination")
    data = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]

    # Return information.
    return data


# Write.


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
    path_rank = os.path.join(dock, "rank")
    utility.create_directory(path_rank)
    path_genes_scores_probabilities = os.path.join(
        path_rank, "data_genes_scores_probabilites.pickle"
    )
    path_genes_scores_probabilities_text = os.path.join(
        path_rank, "data_genes_scores_probabilites.tsv"
    )
    path_genes_selection_rank = os.path.join(
        path_rank, "data_genes_selection_rank.pickle"
    )
    path_genes_selection_rank_text = os.path.join(
        path_rank, "data_genes_selection_rank.tsv"
    )
    path_export_genes_selection = os.path.join(
        path_rank, "export_genes_selection.tsv"
    )
    path_export_genes_total = os.path.join(
        path_rank, "export_genes_total.tsv"
    )
    path_export_genes_rank = os.path.join(
        path_rank, "export_genes_rank.tsv"
    )

    # Write information to file.
    information["data_genes_scores_probabilities"].to_pickle(
        path_genes_scores_probabilities
    )
    information["data_genes_selection_rank"].to_pickle(
        path_genes_selection_rank
    )
    information["data_genes_scores_probabilities"].to_csv(
        path_or_buf=path_genes_scores_probabilities_text,
        columns=None,
        sep="\t",
        na_rep="",
        header=True,
        index=True,
    )
    information["data_genes_selection_rank"].to_csv(
        path_or_buf=path_genes_selection_rank_text,
        columns=None,
        sep="\t",
        na_rep="",
        header=True,
        index=True,
    )

    # Exportation
    information["data_export_genes_selection"].to_csv(
        path_or_buf=path_export_genes_selection,
        columns=None,
        sep="\t",
        na_rep="",
        header=False,
        index=False,
    )
    information["data_export_genes_total"].to_csv(
        path_or_buf=path_export_genes_total,
        columns=None,
        sep="\t",
        na_rep="",
        header=False,
        index=False,
    )
    information["data_export_genes_rank"].to_csv(
        path_or_buf=path_export_genes_rank,
        columns=None,
        sep="\t",
        na_rep="",
        header=False,
        index=False,
    )

    pass





###########################
# TODO: stuff below here needs work...


# Thorough summary.

# TODO: I think this is the most useful of the functions below...
# TODO: ***for each measure separately***, split the gene ranks into low, middle, and high groups
# TODO: I want to compare the scores and distributions of genes in these groups for each measure as a quality control
def define_report_genes(
    data_summary_genes=None,
    rank=None
):
    """
    Defines genes of interest for thorough summary in reports.

    arguments:
        data_summary_genes (object): Pandas data frame of genes' identifiers,
            names, and probability values by measures of modality
        rank (str): key of data column to use for ranks

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Collect genes.
    genes = list()
    # Collect sample genes from different ranges of modality scores.
    # Select ranges.
    data_one = data_summary_genes.loc[(data_summary_genes[rank] < 0.01)]
    data_two = (data_summary_genes.loc[
        (data_summary_genes[rank] > 0.01) & (data_summary_genes[rank] < 0.05)
    ])
    data_three = (data_summary_genes.loc[
        (data_summary_genes[rank] > 0.05) & (data_summary_genes[rank] < 0.1)
    ])
    data_four = data_summary_genes.loc[(data_summary_genes[rank] > 0.1)]
    # Extract genes' identifiers.
    genes_one = data_one.index.to_list()
    genes_two = data_two.index.to_list()
    genes_three = data_three.index.to_list()
    genes_four = data_four.index.to_list()
    # Select sample from each range of genes.
    gene_one = random.choice(genes_one)
    gene_two = random.choice(genes_two)
    gene_three = random.choice(genes_three)
    gene_four = random.choice(genes_four)
    genes.append(gene_one)
    genes.append(gene_two)
    genes.append(gene_three)
    genes.append(gene_four)
    # Include custom genes.
    # UTY
    genes.append("ENSG00000183878")
    # TAPBP
    genes.append("ENSG00000231925")
    # XBP1
    genes.append("ENSG00000100219")
    # PHGDH
    genes.append("ENSG00000092621")
    return genes


def report_gene_abundance_distribution_real(
    data_gene_signals=None,
):
    """
    Reports the real distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.

    raises:

    returns:
        (list<float>): values of gene abundances

    """

    # Analyze gene's signals.
    information = pipe.analyze_gene_signal_distribution(
        data_gene_signals=data_gene_signals
    )
    values = information["values"]
    return values


def report_gene_abundance_distribution_shuffle(
    data_gene_signals=None,
    shuffles=None,
):
    """
    Reports the shuffle distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues
        shuffles (list<list<list<int>>>): Matrices of indices

    raises:

    returns:
        (list<float>): values of gene abundances

    """

    # Collect shuffle distribution of values.
    values = pipe.collect_shuffle_distribution_value(
        data_gene_signals=data_gene_signals,
        shuffles=shuffles
    )
    return values


def report_gene_modality_scores_distributions(
    gene_scores_distributions=None,
):
    """
    Reports the real distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        name (str): name of gene
        gene_scores_distributions (dict<dict>): information about a gene's
            scores and distributions for modality
        path (str): path to a directory


    raises:

    returns:
        (dict<dict>): scores and distributions by multiple measures

    """

    # Report modality scores.
    # Collect scores and distributions.
    information = dict()
    for type in ["coefficient", "dip", "mixture", "combination"]:
        # Access information.
        score = gene_scores_distributions["scores"][type]
        values = gene_scores_distributions["distributions"][type]
        # Compile information.
        information[type] = dict()
        information[type]["score"] = score
        information[type]["distribution"] = values
    # Return information.
    return information


def prepare_reports_genes(
    genes=None,
    data_gene_annotation=None,
    genes_signals_patients_tissues=None,
    shuffles=None,
    genes_scores_distributions=None,
):
    """
    Prepares thorough reports about analyses for a few genes of interest.

    arguments:
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        genes_signals_patients_tissues (dict<object>): Collection of matrices.
        shuffles (list<list<list<int>>>): Matrices of indices.
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:
        (dict): information for reports

    """

    # Collect reports.
    reports = dict()
    # Iterate on genes.
    for gene in genes:
        report = prepare_report_gene(
            gene=gene,
            data_gene_annotation=data_gene_annotation,
            genes_signals_patients_tissues=genes_signals_patients_tissues,
            shuffles=shuffles,
            genes_scores_distributions=genes_scores_distributions,
        )
        # Collect report.
        reports[gene] = report
    return reports


def prepare_report_gene(
    gene=None,
    data_gene_annotation=None,
    genes_signals_patients_tissues=None,
    shuffles=None,
    genes_scores_distributions=None,
):
    """
    Prepares a thorough report about analyses for a single gene of interest.

    arguments:
        gene (str): identifier of a single gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        genes_signals_patients_tissues (dict<object>): Collection of matrices.
        shuffles (list<list<list<int>>>): Matrices of indices.
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:
        (dict): information for report

    """

    # Access information.
    name = data_gene_annotation.loc[gene, "gene_name"]
    data_gene_signals = genes_signals_patients_tissues[gene].copy(deep=True)
    gene_scores_distributions = genes_scores_distributions[gene]

    # Report real distribution across patients of a gene's aggregate
    # abundance across tissues.
    abundances_real = report_gene_abundance_distribution_real(
        data_gene_signals=data_gene_signals,
    )
    # Report shuffle distribution across patients of a gene's aggregate
    # abundance across tissues.
    abundances_shuffle = report_gene_abundance_distribution_shuffle(
        data_gene_signals=data_gene_signals,
        shuffles=shuffles,
    )
    # Report distribution across shuffles of a gene's modality score.
    scores_distributions = report_gene_modality_scores_distributions(
        gene_scores_distributions=gene_scores_distributions,
    )
    # Compile information.
    information = dict()
    information["name"] = name
    information["abundances_real"] = abundances_real
    information["abundances_shuffle"] = abundances_shuffle
    information["scores_distributions"] = scores_distributions
    # Return information.
    return information


def write_product_integration(dock=None, information=None):
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
    path_analysis = os.path.join(dock, "analysis")
    utility.create_directory(path_analysis)
    path_probabilities = os.path.join(
        path_analysis, "data_genes_probabilities.pickle"
    )
    path_summary = os.path.join(
        path_analysis, "data_summary_genes.pickle"
    )
    path_summary_text = os.path.join(
        path_analysis, "data_summary_genes.txt"
    )
    path_summary_text_alternative = os.path.join(
        path_analysis, "data_summary_genes_alternative.txt"
    )
    path_data_rank = os.path.join(
        path_analysis, "data_rank_genes.txt"
    )
    path_rank_ensembl = os.path.join(
        path_analysis, "genes_ranks_ensembl.txt"
    )
    path_rank_hugo = os.path.join(
        path_analysis, "genes_ranks_hugo.txt"
    )
    path_genes_ensembl = os.path.join(
        path_analysis, "genes_ensembl.txt"
    )
    path_genes_hugo = os.path.join(
        path_analysis, "genes_hugo.txt"
    )
    path_genes_report = os.path.join(
        path_analysis, "genes_report.txt"
    )
    path_reports = os.path.join(
        path_analysis, "reports.pickle"
    )
    # Write information to file.
    information["data_genes_probabilities"].to_pickle(path_summary)
    information["data_summary_genes"].to_pickle(path_summary)
    information["data_summary_genes"].to_csv(
        path_or_buf=path_summary_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_summary_genes"].reset_index(
        level="identifier", inplace=True
    )
    summary_genes = utility.convert_dataframe_to_records(
        data=information["data_summary_genes"]
    )
    utility.write_file_text_table(
        information=summary_genes,
        path_file=path_summary_text_alternative,
        names=summary_genes[0].keys(),
        delimiter="\t",
        header=True,
    )
    information["data_rank_genes"].to_csv(
        path_or_buf=path_data_rank,
        sep="\t",
        header=True,
        index=True,
    )
    information["data_rank_genes_ensembl"].to_csv(
        path_or_buf=path_rank_ensembl,
        sep="\t",
        header=False,
        index=False,
    )
    information["data_rank_genes_hugo"].to_csv(
        path_or_buf=path_rank_hugo,
        sep="\t",
        header=False,
        index=False,
    )
    utility.write_file_text_list(
        elements=information["genes_ensembl"],
        delimiter="\n",
        path_file=path_genes_ensembl
    )
    utility.write_file_text_list(
        elements=information["genes_hugo"],
        delimiter="\n",
        path_file=path_genes_hugo
    )
    utility.write_file_text_list(
        elements=information["genes_report"],
        delimiter="\n",
        path_file=path_genes_report
    )
    with open(path_reports, "wb") as file_product:
        pickle.dump(information["reports"], file_product)

    pass

# TODO: obsolete?
def write_product_reports(dock=None, information=None):
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
    path_analysis = os.path.join(dock, "analysis")
    utility.create_directory(path_analysis)
    path_reports = os.path.join(path_analysis, "reports")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_reports)
    utility.create_directory(path_reports)
    # Iterate on reports.
    for gene in information["genes_report"]:
        # Access information.
        name = information["reports"][gene]["name"]
        # Specify directories and files.
        path_report = os.path.join(
            path_reports, (name + "_reports.pickle")
        )
        # Write information to file.

        pass

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

    # Remove previous files to avoid version or batch confusion.
    path_rank = os.path.join(dock, "rank")
    utility.remove_directory(path=path_rank)

    # Read source information from file.
    source = read_source(dock=dock)

    # Calculate genes' probabilities of bimodality.
    genes_probabilities = calculate_probabilities_genes(
        genes_scores=source["genes_scores"],
        genes_permutations=source["genes_permutations"],
    )

    # Organize gene summary.
    data_genes_scores_probabilities = organize_genes_scores_probabilities(
        genes_scores=source["genes_scores"],
        genes_permutations=source["genes_permutations"],
        genes_probabilities=genes_probabilities,
        data_gene_distribution_report=source["data_gene_distribution_report"],
        data_gene_annotation=source["data_gene_annotation"],
    )
    print(data_genes_scores_probabilities)
    print(data_genes_scores_probabilities.iloc[0:10, 0:13])

    # Filter genes by probabilities of each bimodality measure.
    genes_measurements = filter_genes_by_probabilities(
        measurements=["dip", "mixture", "coefficient", "combination"],
        threshold=0.01, # set threshold at 1/1000 to give 10-fold grace
        data_genes_scores_probabilities=data_genes_scores_probabilities,
    )

    # Select genes that pass filters by at least two of the three primary
    # bimodality measurements.
    # Multiple hypothesis threshold:
    # 0.05 (probability) / 16000 (genes) = 3E-6
    # Overall probability is a multiple of each measurement.
    # (0.001)^3 = 1E-9
    # That calculation of overall probability is inaccurate as bimodality
    # measurements are not independent.
    genes_selection = utility.select_elements_by_sets(
        names=["dip", "mixture", "coefficient"],
        sets=genes_measurements,
        count=3,
    )
    utility.print_terminal_partition(level=2)
    print(
        "selection of genes by multiple measurements: " +
        str(len(genes_selection))
    )
    utility.print_terminal_partition(level=2)

    # Select and rank genes by probabilities.
    data_genes_selection_rank = select_rank_genes(
        data_genes_scores_probabilities=data_genes_scores_probabilities,
        selection=genes_selection,
        rank="combination",
    )
    print(data_genes_selection_rank)

    # Exportation of genes for functional analysis.
    utility.print_terminal_partition(level=1)
    data_export_genes_selection = export_genes_function(
        data_genes=data_genes_selection_rank,
        identifier="ensembl",
        rank=False,
    )
    print(data_export_genes_selection)
    data_export_genes_total = export_genes_function(
        data_genes=data_genes_scores_probabilities,
        identifier="ensembl",
        rank=False,
    )
    print(data_export_genes_total)

    data_export_genes_rank = export_genes_function(
        data_genes=data_genes_scores_probabilities,
        identifier="ensembl",
        rank=True,
    )
    print(data_export_genes_rank)

    # Compile information.
    information = {
        "data_genes_scores_probabilities": data_genes_scores_probabilities,
        "data_genes_selection_rank": data_genes_selection_rank,
        "data_export_genes_selection": data_export_genes_selection,
        "data_export_genes_total": data_export_genes_total,
        "data_export_genes_rank": data_export_genes_rank,
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    if False:

        # Define genes of interest.
        # Genes of interest are few for thorough summary.
        genes_report = define_report_genes(
            data_summary_genes=data_summary_genes,
            rank="combination"
        )

        # Organize thorough summaries for a few genes of interest.
        reports = prepare_reports_genes(
            genes=genes_report,
            data_gene_annotation=source["data_gene_annotation"],
            genes_signals_patients_tissues=(
                source["genes_signals_patients_tissues"]
            ),
            shuffles=source["shuffles"][0:500],
            genes_scores_distributions=source["genes_scores_distributions"],
        )

    pass


if (__name__ == "__main__"):
    execute_procedure()
