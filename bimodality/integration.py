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
    path_selection = os.path.join(dock, "selection")
    path_gene_annotation = os.path.join(
        path_selection, "data_gene_annotation.pickle"
    )
    path_samples_properties = os.path.join(
        path_selection, "data_samples_tissues_persons.pickle"
    )
    path_persons_properties = os.path.join(
        path_selection, "data_persons_properties.pickle"
    )

    path_split = os.path.join(dock, "split")
    path_genes_split = os.path.join(path_split, "genes.txt")

    path_distribution = os.path.join(dock, "distribution")
    path_collection = os.path.join(path_distribution, "collection")
    path_distribution_report = os.path.join(
        path_collection, "data_gene_report.pickle"
    )
    path_data_signals_genes_persons = os.path.join(
        path_collection, "data_signals_genes_persons.pickle"
    )

    path_candidacy = os.path.join(dock, "candidacy")
    path_genes_candidacy = os.path.join(
        path_candidacy, "genes_candidacy.pickle"
    )

    path_probability = os.path.join(dock, "probability")
    path_genes_scores = os.path.join(
        path_probability, "genes_scores.pickle"
    )
    path_genes_probabilities = os.path.join(
        path_probability, "genes_probabilities.pickle"
    )
    path_genes_probability = os.path.join(
        path_probability, "genes_probability.pickle"
    )

    path_heritability = os.path.join(dock, "heritability")
    path_collection = os.path.join(path_heritability, "collection")
    path_genes_heritabilities_simple = os.path.join(
        path_collection, "genes_heritabilities_simple.pickle"
    )
    path_genes_heritabilities_complex = os.path.join(
        path_collection, "genes_heritabilities_complex.pickle"
    )
    path_genes_heritability = os.path.join(
        path_collection, "genes_heritability.pickle"
    )

    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_samples_tissues_persons = pandas.read_pickle(path_samples_properties)
    data_persons_properties = pandas.read_pickle(path_persons_properties)
    genes_split = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes_split,
    )

    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    data_gene_distribution_report = pandas.read_pickle(
        path_distribution_report
    )
    with open(path_genes_candidacy, "rb") as file_source:
        genes_candidacy = pickle.load(file_source)

    with open(path_genes_scores, "rb") as file_source:
        genes_scores = pickle.load(file_source)
    with open(path_genes_probabilities, "rb") as file_source:
        genes_probabilities = pickle.load(file_source)
    with open(path_genes_probability, "rb") as file_source:
        genes_probability = pickle.load(file_source)

    with open(path_genes_heritabilities_simple, "rb") as file_source:
        genes_heritabilities_simple = pickle.load(file_source)
    with open(path_genes_heritabilities_complex, "rb") as file_source:
        genes_heritabilities_complex = pickle.load(file_source)
    with open(path_genes_heritability, "rb") as file_source:
        genes_heritability = pickle.load(file_source)

    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_persons_properties": data_persons_properties,
        "genes_split": genes_split,
        "data_signals_genes_persons": data_signals_genes_persons,
        "data_gene_distribution_report": data_gene_distribution_report,
        "genes_candidacy": genes_candidacy,

        "genes_scores": genes_scores,
        "genes_probabilities": genes_probabilities,
        "genes_probability": genes_probability,

        "genes_heritabilities_simple": genes_heritabilities_simple,
        "genes_heritabilities_complex": genes_heritabilities_complex,
        "genes_heritability": genes_heritability,
    }


# Summary


def access_gene_name(
    identifier=None,
    data_gene_annotation=None,
):
    """
    Combines elements in ordered pairs.

    arguments:
        identifier (str): identifier of gene
        data_gene_annotation (object): Pandas data frame of genes' annotations


    returns:
        (str): name of gene

    raises:

    """

    name = data_gene_annotation.loc[identifier, "gene_name"].replace(".", "-")
    return name


def organize_genes_integration(
    genes=None,
    data_gene_annotation=None,
    data_samples_tissues_persons=None,
    data_gene_distribution_report=None,
    genes_scores=None,
    genes_probabilities=None,
    genes_heritabilities_simple=None,
    genes_heritabilities_complex=None,
):
    """
    Organizes information about genes.

    Data structure
    identifier        name    coefficient  p_coefficient  dip     p_dip   ...
     ENSG00000078808   TK421   0.75         0.003          0.12    0.005  ...
     ENSG00000029363   R2D2    0.55         0.975          0.23    0.732  ...

    arguments:
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_distribution_report (object): Pandas data frame of
            information about genes' distributions
        genes_scores (dict<dict<float>>): information about genes' scores
        genes_probabilities (dict<dict<float>>): information about genes'
            probabilities
        genes_heritabilities_simple (dict<dict<float>>): information about
            genes' heritabilities from a model without covariates
        genes_heritabilities_complex (dict<dict<float>>): information about
            genes' heritabilities from a model with covariates

    raises:

    returns:
        (object): Pandas data frame of genes' properties

    """

    # Organize data.
    data_gene_distribution_report.set_index(
        ["identifier"],
        append=False,
        drop=True,
        inplace=True
    )

    # Collect information about genes.
    records = list()
    # Iterate on genes.
    for gene in genes:
        # Access information about gene.

        name = access_gene_name(
            identifier=gene,
            data_gene_annotation=data_gene_annotation,
        )
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

        heritabilities = genes_heritabilities_complex[gene]
        heritability_proportion = heritabilities["proportion"]
        heritability_probability = heritabilities["probability"]

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

        record["heritability_proportion"] = heritability_proportion
        record["heritability_probability"] = heritability_probability

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
        "heritability_proportion",
        "heritability_probability",
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
    data["heritability_proportion"] = (
        data["heritability_proportion"].astype("float32")
    )
    data["heritability_probability"] = (
        data["heritability_probability"].astype("float32")
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


# Correlations


def calculate_pairwise_signal_correlations(
    entities=None,
    data_signals=None,
):
    """
    Collects genes' pantissue signals across persons.

    arguments:
        entities (list<str>): identifiers of entities
        data_signals (object): Pandas data frame of entities' signals across
            observations

    raises:

    returns:
        (dict): correlation coefficients for pairs of entities

    """

    # Determine orderless, pairwise combinations of genes.
    pairs = utility.combine_unique_elements_pairwise_order(
        elements=entities,
    )
    print("count of orderless pairs of unique entities:" + str(len(pairs)))
    utility.print_terminal_partition(level=1)

    # Collect counts and correlations across pairs of features.
    correlations = dict()
    # Iterate on pairs of features.
    for pair in pairs:
        # Select data for pair of features.
        data_pair = data_signals.loc[:, list(pair)]
        # Remove observations with missing values for either feature.
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        # Determine observations for which pair of features has matching
        # values.
        count = data_pair.shape[0]
        # Calculate correlation.
        if count > 1:
            correlation = scipy.stats.pearsonr(
                data_pair[pair[0]].values,
                data_pair[pair[1]].values,
            )[0]
            pass
        else:
            correlation = float("nan")
            pass

        # Collect information.
        if not pair[0] in correlations:
            correlations[pair[0]] = dict()
            pass
        correlations[pair[0]][pair[1]] = dict()
        correlations[pair[0]][pair[1]]["count"] = count
        correlations[pair[0]][pair[1]]["correlation"] = correlation

        pass

    # Return information.
    return correlations


def organize_correlations_matrix(
    entities=None,
    correlations=None,
):
    """
    Collects genes' pantissue signals across persons.

    arguments:
        entities (list<str>): identifiers of entities
        (dict): correlation coefficients for pairs of entities

    raises:

    returns:
        (object): Pandas data frame of correlations between entities

    """

    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on dimension one of entities.
    for entity_one in entities:
        # Collect information.
        record = dict()
        record["gene"] = entity_one

        # Iterate on dimension two of entities.
        for entity_two in entities:
            # Access value.
            if entity_two in correlations[entity_one]:
                correlation = (
                    correlations[entity_one][entity_two]["correlation"]
                )
            else:
                correlation = float("nan")
            pass
            # Collect information.
            record[entity_two] = correlation

        # Collect information.
        records.append(record)
        pass

    # Organize data.
    data = utility.convert_records_to_dataframe(records=records)
    data.sort_values(
        by=["gene"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    data.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )
    data.rename_axis(
        columns="genes",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data


# Rank

def select_rank_genes(
    data_genes_integration=None,
    selection=None,
    rank=None,
):
    """
    Prepares ranks of genes for subsequent functional analyses.

    arguments:
        data_genes_integration (object): Pandas data frame of genes'
            properties, bimodality scores, and probabilities
        selection (list<str>): genes for selection
        rank (str): property to use for ranks

    raises:

    returns:
        (object): Pandas data frame of genes' properties

    """

    # Copy data.
    data_copy = data_genes_integration.copy(deep=True)
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
    path_integration = os.path.join(dock, "integration")
    utility.create_directory(path_integration)

    path_genes_integration = os.path.join(
        path_integration, "genes_integration.pickle"
    )
    path_data_genes_integration = os.path.join(
        path_integration, "data_genes_integration.pickle"
    )
    path_data_genes_correlations = os.path.join(
        path_integration, "data_genes_correlations.pickle"
    )

    path_data_genes_selection = os.path.join(
        path_integration, "data_genes_selection.pickle"
    )
    path_data_genes_selection_text = os.path.join(
        path_integration, "data_genes_selection.tsv"
    )

    path_export_genes_selection = os.path.join(
        path_integration, "export_genes_selection.tsv"
    )
    path_export_genes_total = os.path.join(
        path_integration, "export_genes_total.tsv"
    )

    # Write information to file.
    with open(path_genes_integration, "wb") as file_product:
        pickle.dump(
            information["genes_integration"], file_product
        )
    information["data_genes_correlations"].to_pickle(
        path_data_genes_correlations
    )

    if False:

        information["data_genes_integration"].to_pickle(
            path_data_genes_integration
        )
        information["data_genes_selection"].to_pickle(
            path_data_genes_selection
        )
        information["data_genes_selection"].to_csv(
            path_or_buf=path_data_genes_selection_text,
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

    pass



##########################Scrap#########################3




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
    path_integration = os.path.join(dock, "integration")
    utility.remove_directory(path=path_integration)

    # Read source information from file.
    source = read_source(dock=dock)

    # Select genes integration of methods from candidacy, probability, and
    # heritability procedures.
    sets_genes_integration = dict()
    sets_genes_integration["candidacy"] = source["genes_candidacy"]
    sets_genes_integration["probability"] = source["genes_probability"]
    sets_genes_integration["heritability"] = source["genes_heritability"]
    genes_integration = utility.select_elements_by_sets(
        names=["candidacy", "probability", "heritability"],
        sets=sets_genes_integration,
        count=3,
    )
    utility.print_terminal_partition(level=2)
    print(
        "selection of genes by candidacy, probability, and heritability: " +
        str(len(genes_integration))
    )
    utility.print_terminal_partition(level=2)

    # Calculate correlations between gene pairs of their pantissue signals
    # across persons.
    correlations_genes = calculate_pairwise_signal_correlations(
        entities=genes_integration,
        data_signals=source["data_signals_genes_persons"],
    )
    print("count of genes: " + str(len(correlations_genes.keys())))
    # Organize data matrix.
    data_genes_correlations = organize_correlations_matrix(
        entities=genes_integration,
        correlations=correlations_genes,
    )
    print(data_genes_correlations)

    if False:

        # Integrate and organize information about all genes.
        data_genes_integration = organize_genes_integration(
            genes=source["genes_split"],
            data_gene_annotation=source["data_gene_annotation"],
            data_samples_tissues_persons=source["data_samples_tissues_persons"],
            data_gene_distribution_report=source["data_gene_distribution_report"],
            genes_scores=source["genes_scores"],
            genes_probabilities=source["genes_probabilities"],
            genes_heritabilities_simple=source["genes_heritabilities_simple"],
            genes_heritabilities_complex=source["genes_heritabilities_complex"],
        )

        # Select genes of interest.
        data_genes_selection = select_rank_genes(
            data_genes_integration=data_genes_integration,
            selection=genes_integration,
            rank="combination",
        )
        print(data_genes_selection)

        # Exportation of genes for functional analysis.
        utility.print_terminal_partition(level=1)
        data_export_genes_selection = export_genes_function(
            data_genes=data_genes_selection,
            identifier="ensembl",
            rank=False,
        )
        print(data_export_genes_selection)
        data_export_genes_total = export_genes_function(
            data_genes=data_genes_integration,
            identifier="ensembl",
            rank=False,
        )
        print(data_export_genes_total)

    # Compile information.
    information = {
        "genes_integration": genes_integration,
        #"data_genes_integration": data_genes_integration,
        "data_genes_correlations": data_genes_correlations,
        #"data_genes_selection": data_genes_selection,
        #"data_export_genes_selection": data_export_genes_selection,
        #"data_export_genes_total": data_export_genes_total,
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
