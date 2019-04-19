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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import pipe
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
    path_assembly = os.path.join(dock, "assembly")
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_signal = os.path.join(
        path_split, "genes_signals_patients_tissues.pickle"
    )
    path_shuffle = os.path.join(dock, "shuffle")
    path_shuffles = os.path.join(
        path_shuffle, "shuffles.pickle"
    )
    path_collection = os.path.join(dock, "collection")
    path_distributions = os.path.join(
        path_collection, "genes_scores_distributions.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    with open(path_signal, "rb") as file_source:
        genes_signals_patients_tissues = pickle.load(file_source)
    with open(path_shuffles, "rb") as file_source:
        shuffles = pickle.load(file_source)
    genes = utility.read_file_text_list(path_genes)
    with open(path_distributions, "rb") as file_source:
        genes_scores_distributions = pickle.load(file_source)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes": genes,
        "genes_signals_patients_tissues": genes_signals_patients_tissues,
        "shuffles": shuffles,
        "genes_scores_distributions": genes_scores_distributions
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
    genes_scores_distributions=None
):
    """
    Calculates probabilities (p-values) from genes' scores and random
    distributions.

    Data structure.
    - genes_p_values (dict)
    -- gene (dict)
    ---- combination (float)
    ---- coefficient (float)
    ---- dip (float)

    arguments:
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:
        (dict<dict<float>>): p-values from genes' scores and distributions


    """

    # Collect information about genes.
    genes_p_values = dict()
    # Extract identifiers of genes with scores.
    genes = list(genes_scores_distributions.keys())
    # Iterate on genes.
    for gene in genes:
        # Create entry for gene.
        genes_p_values[gene] = dict()
        # Access information about gene's scores and distributions.
        gene_scores = genes_scores_distributions[gene]
        score_combination = gene_scores["scores"]["combination"]
        score_coefficient = gene_scores["scores"]["coefficient"]
        score_dip = gene_scores["scores"]["dip"]
        distribution_combination = gene_scores["distributions"]["combination"]
        distribution_coefficient = gene_scores["distributions"]["coefficient"]
        distribution_dip = gene_scores["distributions"]["dip"]
        # Calculate p-values.
        # These scores of bimodality indicate greater bimodality as values
        # increase.
        # The scores are unidirectional, so the hypothesis is unidirectional.
        # The p-value is the probability of obtaining by random chance a value
        # equal to or greater than the actual score.
        probability_combination = calculate_probability_equal_greater(
            value=score_combination,
            distribution=distribution_combination
        )
        probability_coefficient = calculate_probability_equal_greater(
            value=score_coefficient,
            distribution=distribution_coefficient
        )
        probability_dip = calculate_probability_equal_greater(
            value=score_dip,
            distribution=distribution_dip
        )
        # Compile information.
        genes_p_values[gene]["combination"] = probability_combination
        genes_p_values[gene]["coefficient"] = probability_coefficient
        genes_p_values[gene]["dip"] = probability_dip
    # Return information.
    return genes_p_values


def organize_summary_genes(
    genes_scores_distributions=None,
    genes_probabilities=None,
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
        genes_probabilities (dict<dict<float>>): p-values from genes' scores
            and distributions
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
        score_combination = gene_scores["scores"]["combination"]
        score_coefficient = gene_scores["scores"]["coefficient"]
        score_dip = gene_scores["scores"]["dip"]
        probability_combination = genes_probabilities[gene]["combination"]
        probability_coefficient = genes_probabilities[gene]["coefficient"]
        probability_dip = genes_probabilities[gene]["dip"]
        # Compile information.
        record["identifier"] = gene
        record["name"] = name
        record["combination"] = score_combination
        record["p_combination"] = probability_combination
        record["coefficient"] = score_coefficient
        record["p_coefficient"] = probability_coefficient
        record["dip"] = score_dip
        record["p_dip"] = probability_dip
        records.append(record)
    # Organize information.
    data = utility.convert_records_to_dataframe(
        records=records
    )
    if False:
        data["coefficient"] = data["coefficient"].astype("float32")
        data["p_coefficient"] = data["p_coefficient"].astype("float32")
        data["dip"] = data["dip"].astype("float32")
        data["p_dip"] = data["p_dip"].astype("float32")
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
        by=["p_combination", "p_coefficient", "p_dip"],
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


def extract_priority_genes(
    data_summary_genes=None, count=None
):
    """
    Extracts identifiers of genes with greatest ranks.

    arguments:
        data_summary_genes (object): Pandas data frame of genes' identifiers,
            names, and probability values by bimodality coefficient and dip
            statistic
        count (int): Count of genes with greatest ranks to select

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    # Extract genes' identifiers.
    data_summary_genes.reset_index(level="identifier", inplace=True)
    genes = data_summary_genes["identifier"].to_list()
    print("count of genes: " + str(len(genes)))
    return genes[0:count]


# Thorough summary.


def define_genes_interest():
    """
    Defines genes of interest for thorough summary.

    arguments:

    raises:

    returns:
        (list<str>): identifiers of genes

    """

    genes = list()
    # TAPBP
    genes.append("ENSG00000231925")
    # XBP1
    genes.append("ENSG00000100219")
    # FOXP4
    genes.append("ENSG00000137166")
    # ADAT3
    genes.append("ENSG00000213638")
    # LDHA
    genes.append("ENSG00000134333")
    # LDHB
    genes.append("ENSG00000111716")
    return genes


def report_gene_abundance_distribution_real(
    name=None,
    data_gene_signals=None,
    dock=None,
):
    """
    Reports the real distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        name (str): name of gene
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues.
        dock (str): path to root or dock directory for source and product
            directories and files.


    raises:

    returns:

    """

    # Analyze gene's signals.
    information = pipe.analyze_gene_signal_distribution(
        data_gene_signals=data_gene_signals
    )
    values = information["values"]

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_distribution_histogram(
        series=values,
        name="",
        bin_method="count",
        bin_count=50,
        label_bins="Bins",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=False,
        position=0
    )
    # Specify directories and files.
    path_analysis = os.path.join(dock, "analysis")
    path_figure = os.path.join(path_analysis, "figure")
    utility.confirm_path_directory(path_figure)
    file = (name + "_abundance_distribution_real.svg")
    path_file = os.path.join(path_figure, file)
    # Write figure.
    plot.write_figure(
        path=path_file,
        figure=figure
    )

    pass


def report_gene_abundance_distribution_shuffle(
    name=None,
    data_gene_signals=None,
    shuffles=None,
    dock=None,
):
    """
    Reports the shuffle distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        name (str): name of gene
        data_gene_signals (object): Pandas data frame of a single gene's
            signals across specific patients and tissues
        shuffles (list<list<list<int>>>): Matrices of indices
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Collect shuffle distribution of values.
    values = pipe.collect_shuffle_distribution_value(
        data_gene_signals=data_gene_signals,
        shuffles=shuffles
    )

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Create figure.
    figure = plot.plot_distribution_histogram(
        series=values,
        name="",
        bin_method="count",
        bin_count=50,
        label_bins="Bins",
        label_counts="Counts",
        fonts=fonts,
        colors=colors,
        line=False,
        position=0
    )
    # Specify directories and files.
    path_analysis = os.path.join(dock, "analysis")
    path_figure = os.path.join(path_analysis, "figure")
    utility.confirm_path_directory(path_figure)
    file = (name + "_abundance_distribution_shuffle.svg")
    path_file = os.path.join(path_figure, file)
    # Write figure.
    plot.write_figure(
        path=path_file,
        figure=figure
    )

    pass


def report_gene_modality_scores_distributions(
    name=None,
    gene_scores_distributions=None,
    dock=None,
):
    """
    Reports the real distribution across patients of a gene's aggregate
    abundance across tissues.

    arguments:
        name (str): name of gene
        gene_scores_distributions (dict<dict>): information about a gene's
            scores and distributions for modality
        dock (str): path to root or dock directory for source and product
            directories and files.


    raises:

    returns:

    """

    # Define fonts.
    fonts = plot.define_font_properties()
    # Define colors.
    colors = plot.define_color_properties()

    # Report modality scores.
    for type in ["combination", "coefficient", "dip"]:
        # Access information.
        score = gene_scores_distributions["scores"][type]
        values = gene_scores_distributions["distributions"][type]

        # Create figure.
        figure = plot.plot_distribution_histogram(
            series=values,
            name="",
            bin_method="count",
            bin_count=50,
            label_bins="Bins",
            label_counts="Counts",
            fonts=fonts,
            colors=colors,
            line=True,
            position=score,
        )
        # Specify directories and files.
        path_analysis = os.path.join(dock, "analysis")
        path_figure = os.path.join(path_analysis, "figure")
        utility.confirm_path_directory(path_figure)
        file = (name + "_score_distribution_" + type + ".svg")
        path_file = os.path.join(path_figure, file)
        # Write figure.
        plot.write_figure(
            path=path_file,
            figure=figure
        )

    pass


def prepare_reports_genes(
    genes=None,
    data_gene_annotation=None,
    genes_signals_patients_tissues=None,
    shuffles=None,
    genes_scores_distributions=None,
    dock=None,
):
    """
    Prepares thorough reports about analyses for a few genes of interest.

    arguments:
        genes (list<str>): identifiers of genes
        data_gene_annotation (object): Pandas data frame of genes' annotations
        genes_signals_patients_tissues (dict<object>): Collection of matrices.
        shuffles (list<list<list<int>>>): Matrices of indices.
        genes_scores_distributions (dict<dict<dict>>): information about genes
        dock (str): path to root or dock directory for source and product
            directories and files.


    raises:

    returns:

    """

    for gene in genes:
        prepare_report_gene(
            gene=gene,
            data_gene_annotation=data_gene_annotation,
            genes_signals_patients_tissues=genes_signals_patients_tissues,
            shuffles=shuffles,
            genes_scores_distributions=genes_scores_distributions,
            dock=dock,
        )
    pass


def prepare_report_gene(
    gene=None,
    data_gene_annotation=None,
    genes_signals_patients_tissues=None,
    shuffles=None,
    genes_scores_distributions=None,
    dock=None,
):
    """
    Prepares a thorough report about analyses for a single gene of interest.

    arguments:
        gene (str): identifier of a single gene
        data_gene_annotation (object): Pandas data frame of genes' annotations
        genes_signals_patients_tissues (dict<object>): Collection of matrices.
        shuffles (list<list<list<int>>>): Matrices of indices.
        genes_scores_distributions (dict<dict<dict>>): information about genes
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Access information.
    name = data_gene_annotation.loc[gene, "gene_name"]
    data_gene_signals = genes_signals_patients_tissues[gene].copy(deep=True)
    gene_scores_distributions = genes_scores_distributions[gene]

    # Summarize real distribution across patients of a gene's aggregate
    # abundance across tissues.
    report_gene_abundance_distribution_real(
        name=name,
        data_gene_signals=data_gene_signals,
        dock=dock,
    )

    # Summarize shuffle distribution across patients of a gene's aggregate
    # abundance across tissues.
    report_gene_abundance_distribution_shuffle(
        name=name,
        data_gene_signals=data_gene_signals,
        shuffles=shuffles,
        dock=dock,
    )

    # Summarize distribution across shuffles of a gene's modality score.
    report_gene_modality_scores_distributions(
        name=name,
        gene_scores_distributions=gene_scores_distributions,
        dock=dock,
    )

    pass


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
    path_analysis = os.path.join(dock, "analysis")
    utility.confirm_path_directory(path_analysis)
    path_probabilities = os.path.join(
        path_analysis, "genes_probabilities.pickle"
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
    path_genes = os.path.join(
        path_analysis, "genes.txt"
    )
    # Write information to file.
    with open(path_probabilities, "wb") as file_product:
        pickle.dump(information["genes_probabilities"], file_product)
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
        delimiter="\t"
    )
    utility.write_file_text_list(
        information=information["genes"], path_file=path_genes
    )
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

    # Read source information from file.
    source = read_source(dock=dock)

    # Calculate probability of each gene's scores by chance.
    genes_probabilities = calculate_probabilities_genes(
        genes_scores_distributions=source["genes_scores_distributions"]
    )
    #print(source["genes_scores_distributions"]["ENSG00000005483"])
    #print(genes_probabilities["ENSG00000005483"])

    # Organize gene summary.
    data_summary_genes = organize_summary_genes(
        genes_scores_distributions=source["genes_scores_distributions"],
        genes_probabilities=genes_probabilities,
        data_gene_annotation=source["data_gene_annotation"],
    )

    print(data_summary_genes.iloc[1700: , 0:10])

    # Extract identifiers of genes of interest.
    # These genes are for further subsequent analysis.
    genes_priority = extract_priority_genes(
        data_summary_genes=data_summary_genes.copy(deep=True),
        count=250
    )

    # Organize thorough summaries for genes of interest.
    prepare_reports_genes(
        genes=define_genes_interest(),
        data_gene_annotation=source["data_gene_annotation"],
        genes_signals_patients_tissues=(
            source["genes_signals_patients_tissues"]
        ),
        shuffles=source["shuffles"][0:500],
        genes_scores_distributions=source["genes_scores_distributions"],
        dock=dock,
    )

    # Compile information.
    information = {
        "genes_probabilities": genes_probabilities,
        "data_summary_genes": data_summary_genes,
        "genes": genes_priority
    }
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
