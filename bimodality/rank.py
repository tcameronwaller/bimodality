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
import metric
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
    path_permutation = os.path.join(dock, "permutation")
    path_collection = os.path.join(dock, "collection")
    path_scores_permutations = os.path.join(
        path_collection, "genes_scores_permutations.pickle"
    )
    # Read information from file.
    genes = utility.read_file_text_list(
        delimiter="\n",
        path_file=path_genes,
    )
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_scores_permutations, "rb") as file_source:
        genes_scores_permutations = pickle.load(file_source)
    # Compile and return information.
    return {
        "genes": genes,
        "data_gene_annotation": data_gene_annotation,
        "genes_scores_permutations": genes_scores_permutations,
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
    genes_scores_permutations=None
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
        genes_scores_permutations (dict<dict<dict>>): information about genes

    raises:

    returns:
        (object): Pandas data frame of probabilities from genes' scores and
            permutations

    """

    # Collect information about genes.
    genes_probabilities = list()
    # Extract identifiers of genes with scores.
    genes = list(genes_scores_permutations.keys())
    # Iterate on genes.
    for gene in genes:
        # Access information about gene's scores and distributions.
        entry = genes_scores_permutations[gene]
        score_coefficient = entry["scores"]["coefficient"]
        score_dip = entry["scores"]["dip"]
        score_mixture = entry["scores"]["mixture"]
        #score_combination = entry["scores"]["combination"]
        permutations_coefficient = entry["permutations"]["coefficient"]
        permutations_dip = entry["permutations"]["dip"]
        permutations_mixture = entry["permutations"]["mixture"]
        #permutations_combination = entry["permutations"]["combination"]
        # Calculate p-values.
        # These scores of bimodality indicate greater bimodality as values
        # increase.
        # The scores are unidirectional, so the hypothesis is unidirectional.
        # The p-value is the probability of obtaining by random chance a value
        # equal to or greater than the actual score.
        probability_coefficient = calculate_probability_equal_greater(
            value=score_coefficient,
            distribution=permutations_coefficient
        )
        probability_dip = calculate_probability_equal_greater(
            value=score_dip,
            distribution=permutations_dip
        )
        probability_mixture = calculate_probability_equal_greater(
            value=score_mixture,
            distribution=permutations_mixture
        )
        #probability_combination = calculate_probability_equal_greater(
        #    value=score_combination,
        #    distribution=permutations_combination
        #)
        # Calculate combination of p-values.
        if False:
            probability_combination = scipy.stats.combine_pvalues(
                [probability_coefficient, probability_dip, probability_mixture],
                method="stouffer",
                weights=None,
            )[1]
        # Compile information.
        record = dict()
        record["gene"] = gene
        record["coefficient"] = probability_coefficient
        record["dip"] = probability_dip
        record["mixture"] = probability_mixture
        #record["combination"] = probability_combination
        genes_probabilities.append(record)
    # Organize information.
    data = utility.convert_records_to_dataframe(
        records=genes_probabilities
    )
    data["coefficient"] = data["coefficient"].astype("float32")
    data["dip"] = data["dip"].astype("float32")
    data["mixture"] = data["mixture"].astype("float32")
    #data["combination"] = data["combination"].astype("float32")
    data.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    data.rename_axis(
        index="gene",
        axis="index",
        copy=False,
        inplace=True,
    )
    data.rename_axis(
        columns="probabilities",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data


def filter_genes_probabilities_threshold(
    data=None,
    threshold=None,
    count=None,
):
    """
    Filters genes to keep only those with scores from at least two of the three
    modality metrics with permutation probabilities below threshold.

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
    # Keep genes with at least 2 metrics below threshold

    # Determine whether values pass threshold.
    # Only consider the original modality metrics for this threshold.
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




###########################
# TODO: stuff below here needs work...

# Rank


def rank_genes(
    data_genes_probabilities=None,
    rank=None,
    method=None,
    threshold=None,
    count=None,
):
    """
    Prepares ranks of genes for subsequent functional analyses.

    arguments:
        data_genes_probabilities (object): Pandas data frame of probabilities
            from genes' scores and distributions
        rank (str): key of data column to use for ranks
        method (str): method of filter to apply, either threshold or count
        threshold (float): value for threshold
        count (int): Count of genes with greatest ranks to select

    raises:

    returns:
        (dict): collection of representations of ranks of genes

    """

    # Copy data.
    data_genes_ranks = data_genes_probabilities.copy(deep=True)
    data_genes_ranks.reset_index(level="gene", inplace=True)
    data_genes_ranks.sort_values(
        by=[rank],
        axis="index",
        ascending=True,
        inplace=True,
    )
    data_genes_ranks.reset_index(
        inplace=True,
        drop=True,
    )
    # Filter genes.
    if method == "threshold":
        data_genes_ranks_filter = (
            data_genes_ranks.loc[(data_genes_ranks[rank] < threshold)]
        )
    elif method == "count":
        data_genes_ranks_filter = (
            data_genes_ranks.iloc[0:count]
        )
    # Organize information.
    data_genes_ranks_filter.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True,
    )
    data_genes_ranks_filter.rename_axis(
        index="gene",
        axis="index",
        copy=False,
        inplace=True,
    )
    data_genes_ranks_filter.rename_axis(
        columns="probabilities",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data_genes_ranks_filter


# Summary


def organize_summary_genes(
    genes_scores_distributions=None,
    data_genes_probabilities=None,
    data_gene_annotation=None,
):
    """
    Collects tissues, patients, and genes' signals for each sample.

    Data structure
    identifier        name    coefficient  p_coefficient  dip     p_dip   ...
     ENSG00000078808   TK421   0.75         0.003          0.12    0.005  ...
     ENSG00000029363   R2D2    0.55         0.975          0.23    0.732  ...

    arguments:
        genes_scores_distributions (dict<dict<dict>>): information about genes
        data_genes_probabilities (object): Pandas data frame of probabilities
            from genes' scores and distributions
        data_gene_annotation (object): Pandas data frame of genes' annotations

    raises:

    returns:
        (object): Pandas data frame of genes' identifiers, names, and
            probability values by bimodality coefficient and dip statistic

    """

    def match_gene_name(identifier):
        name = data_gene_annotation.loc[identifier, "gene_name"]
        if len(name) < 1:
            print("*************Error?****************")
        return name

    def calculate_mean_probability(row):
        return statistics.mean([
            row["p_coefficient"],
            row["p_dip"],
            row["p_mixture"],
            row["p_combination"]
        ])

    # Collect information about genes.
    records = list()
    # Extract identifiers of genes with scores.
    genes = list(genes_scores_distributions.keys())
    # Iterate on genes.
    for gene in genes:
        # Create record for gene.
        record = dict()
        # Access information about gene.
        name = match_gene_name(gene)
        gene_scores = genes_scores_distributions[gene]
        score_coefficient = gene_scores["scores"]["coefficient"]
        score_dip = gene_scores["scores"]["dip"]
        score_mixture = gene_scores["scores"]["mixture"]
        score_combination = gene_scores["scores"]["combination"]
        probability_coefficient = (
            data_genes_probabilities.loc[gene, "coefficient"]
        )
        probability_dip = data_genes_probabilities.loc[gene, "dip"]
        probability_mixture = data_genes_probabilities.loc[gene, "mixture"]
        probability_combination = (
            data_genes_probabilities.loc[gene, "combination"]
        )
        # Compile information.
        record["identifier"] = gene
        record["name"] = name
        record["coefficient"] = score_coefficient
        record["p_coefficient"] = probability_coefficient
        record["dip"] = score_dip
        record["p_dip"] = probability_dip
        record["mixture"] = score_mixture
        record["p_mixture"] = probability_mixture
        record["combination"] = score_combination
        record["p_combination"] = probability_combination
        records.append(record)
    # Organize information.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    # Calculate mean probabilities from all scores.
    # This mean probability will facilitate ranking.
    data["p_mean"] = data.apply(
        calculate_mean_probability,
        axis="columns",
    )
    data["coefficient"] = data["coefficient"].astype("float32")
    data["p_coefficient"] = data["p_coefficient"].astype("float32")
    data["dip"] = data["dip"].astype("float32")
    data["p_dip"] = data["p_dip"].astype("float32")
    data["mixture"] = data["mixture"].astype("float32")
    data["p_mixture"] = data["p_mixture"].astype("float32")
    data["combination"] = data["combination"].astype("float32")
    data["p_combination"] = data["p_combination"].astype("float32")
    data["p_mean"] = data["p_mean"].astype("float32")
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
        by=["p_mean"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Return information.
    return data


def collect_genes_names(
    data_gene_annotation=None,
    data_gene_score=None
):
    """
    Collects names of genes.

    arguments:
        data_gene_annotation (object): Pandas data frame of genes' annotations.
        data_gene_score (object): Pandas data frame of genes' scores.

    raises:

    returns:
        (object): Pandas data frame of genes' scores.

    """

    def match_gene_name(identifier):
        return data_gene_annotation.loc[identifier, "gene_name"]
    data_gene_score.reset_index(level=["gene"], inplace=True)
    data_gene_score["name"] = (
        data_gene_score["gene"].apply(match_gene_name)
    )
    data_gene_score.set_index(
        ["gene"],
        append=False,
        drop=True,
        inplace=True
    )
    return data_gene_score


# Thorough summary.

# TODO: I think this is the most useful of the functions below...
# TODO: ***for each metric separately***, split the gene ranks into low, middle, and high groups
# TODO: I want to compare the scores and distributions of genes in these groups for each metric as a quality control
def define_report_genes(
    data_summary_genes=None,
    rank=None
):
    """
    Defines genes of interest for thorough summary in reports.

    arguments:
        data_summary_genes (object): Pandas data frame of genes' identifiers,
            names, and probability values by metrics of modality
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
        (dict<dict>): scores and distributions by multiple metrics

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


# Write.


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
    path_rank_availability = os.path.join(
        path_rank, "data_gene_rank_availability.pickle"
    )
    path_rank_imputation = os.path.join(
        path_rank, "data_gene_rank_imputation.pickle"
    )
    path_genes_availability = os.path.join(
        path_rank, "genes_availability.pickle"
    )
    path_genes_imputation = os.path.join(
        path_rank, "genes_imputation.pickle"
    )
    path_genes_consensus = os.path.join(
        path_rank, "genes_consensus.pickle"
    )
    path_genes_consensus_text = os.path.join(
        path_rank, "genes_consensus.txt"
    )
    path_genes_null_text = os.path.join(
        path_rank, "genes_null.txt"
    )

    # Write information to file.
    information["data_gene_rank_availability"].to_pickle(
        path_rank_availability
    )
    information["data_gene_rank_imputation"].to_pickle(
        path_rank_imputation
    )
    with open(path_genes_availability, "wb") as file_product:
        pickle.dump(information["genes_availability"], file_product)
    with open(path_genes_imputation, "wb") as file_product:
        pickle.dump(information["genes_imputation"], file_product)
    with open(path_genes_consensus, "wb") as file_product:
        pickle.dump(information["genes_consensus"], file_product)

    utility.write_file_text_list(
        elements=information["genes_consensus"],
        delimiter="\n",
        path_file=path_genes_consensus_text
    )
    utility.write_file_text_list(
        elements=information["genes_null"],
        delimiter="\n",
        path_file=path_genes_null_text
    )


    pass


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

    #############################
    # TODO:
    # 1. calculate probabilities
    # 2. split genes into 3 copies for each metric --> dip, mixture, coefficient
    # 3. filter each set of genes by probability threshold
    # 4. rank each set of genes by bimodality metric
    # 5. filter each set of genes for top 10% of genes by bimodal metric rank
    # 6. keep genes that are in the final ranked and filtered lists by at least 2 of the 3 metrics
    # 7. rank final genes by combination score
    # 8. plot examples?
    ##############################

    # Read source information from file.
    source = read_source(dock=dock)

    # Calculate genes' probabilities of bimodality.
    data_genes_probabilities = calculate_probabilities_genes(
        genes_scores_permutations=source["genes_scores_permutations"],
    )
    print(data_genes_probabilities)

    # Filter genes by probabilities.
    # Keep genes only if scores from at least two of the three modality metrics
    # have permutation probabilities below threshold.
    data_gene_probabilites_threshold = filter_genes_probabilities_threshold(
        data=data_genes_probabilities,
        threshold=0.05,
        count=3,
    )
    print(data_gene_probabilites_threshold)


    # TODO: after ranking genes...
    # TODO: Prepare a gene report, similar to the report from the distribution procedure

    if False:
        # Rank genes by probabilities.
        # TODO: rank genes separately by "coefficient", "dip", and "mixture"...
        # TODO: Then take the union of the top genes from each ranking...
        data_gene_rank = rank_genes(
            data_genes_probabilities=data_genes_probabilities,
            rank="combination",
            method="threshold",
            threshold=0.001,
            count=1000,
        )
        print(data_gene_rank)

        # Extract identifiers of genes.
        genes_availability = data_gene_rank_availability.index.to_list()
        genes_imputation = data_gene_rank_imputation.index.to_list()

        # Combine consensus genes.
        genes_consensus = utility.filter_unique_union_elements(
            list_one=genes_availability,
            list_two=genes_imputation,
        )
        print("Count of availability genes: " + str(len(genes_availability)))
        print("Count of imputation genes: " + str(len(genes_imputation)))
        print("Count of consensus genes: " + str(len(genes_consensus)))

        # Determine genes not in the consensus list.
        genes_null = utility.filter_unique_exclusion_elements(
            list_exclusion=genes_consensus,
            list_total=source["genes"],
        )
        print("Count of null genes: " + str(len(genes_null)))
        print("Count of total genes: " + str(len(source["genes"])))

        # Organize gene summary.
        data_summary_genes = organize_summary_genes(
            genes_scores_distributions=source["genes_scores_distributions"],
            data_genes_probabilities=data_genes_probabilities,
            data_gene_annotation=source["data_gene_annotation"],
        )
        print(data_summary_genes.iloc[0:25, 0:10])

        # Define genes of interest.
        # Genes of interest are few for thorough summary.
        genes_report = define_report_genes(
            data_summary_genes=data_summary_genes,
            rank="p_mean"
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

        # Compile information.
        information = {
            "data_genes_probabilities": data_genes_probabilities,
            "data_summary_genes": data_summary_genes,
            "data_rank_genes": ranks["data_rank_genes"],
            "data_rank_genes_ensembl": ranks["data_rank_genes_ensembl"],
            "data_rank_genes_hugo": ranks["data_rank_genes_hugo"],
            "genes_ensembl": ranks["genes_ensembl"],
            "genes_hugo": ranks["genes_hugo"],
            "genes_report": genes_report,
            "reports": reports,
        }

        # Compile information.
        information = {
            "data_gene_rank_availability": data_gene_rank_availability,
            "data_gene_rank_imputation": data_gene_rank_imputation,
            "genes_availability": genes_availability,
            "genes_imputation": genes_imputation,
            "genes_consensus": genes_consensus,
            "genes_null": genes_null,
        }
        #Write product information to file.
        write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
