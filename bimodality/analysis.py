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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

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
    path_assembly = os.path.join(dock, "assembly")
    path_gene_annotation = os.path.join(
        path_assembly, "data_gene_annotation.pickle"
    )
    path_split = os.path.join(dock, "split")
    path_genes = os.path.join(
        path_split, "genes.txt"
    )
    path_pipe = os.path.join(dock, "pipe")
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    genes = utility.read_file_text_list(path_genes)
    genes_scores_distributions = read_collect_genes_scores_distributions(
        path=path_pipe
    )
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "genes": genes,
        "genes_scores_distributions": genes_scores_distributions
    }


def read_collect_genes_scores_distributions(path=None):
    """
    Collects information about genes.

    Data structure.
    - genes_scores_distributions (dict)
    -- gene (dict)
    --- scores (dict)
    ---- coefficient (float)
    ---- dip (float)
    --- distributions (dict)
    ---- coefficient (list)
    ----- value (float)
    ---- dip (list)
    ----- value (float)

    arguments:
        path (str): path to directory

    raises:

    returns:
        (dict<dict<dict>>): information about genes

    """

    # Collect information about genes.
    genes_scores_distributions = dict()
    # Iterate on directories for genes.
    directories = os.listdir(path)
    for directory in directories:
        # Create entry for gene.
        genes_scores_distributions[directory] = dict()
        # Specify directories and files.
        path_directory = os.path.join(path, directory)
        path_scores = os.path.join(path_directory, "scores.pickle")
        path_distributions = os.path.join(
            path_directory, "distributions.pickle"
        )
        # Read information from file.
        with open(path_scores, "rb") as file_source:
            scores = pickle.load(file_source)
        with open(path_distributions, "rb") as file_source:
            distributions = pickle.load(file_source)
        # Compile information.
        genes_scores_distributions[directory]["scores"] = scores
        genes_scores_distributions[directory]["distributions"] = distributions
    # Return information.
    return genes_scores_distributions


def check_genes_scores_distributions(
    genes=None,
    genes_scores_distributions=None
):
    """
    Checks and summarizes information about genes.

    arguments:
        genes (list<str>): identifiers of genes
        genes_scores_distributions (dict<dict<dict>>): information about genes

    raises:

    returns:


    """

    #print(genes_scores_distributions)
    #print(genes_scores_distributions["ENSG00000005483"])

    # Extract identifiers of genes with scores.
    genes_scores = list(genes_scores_distributions.keys())
    # Check whether all genes have scores.

    pass


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


def calculate_p_values_genes(
    genes_scores_distributions=None
):
    """
    Calculates p-values from genes' scores and random distributions.

    Data structure.
    - genes_p_values (dict)
    -- gene (dict)
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
        score_coefficient = gene_scores["scores"]["coefficient"]
        score_dip = gene_scores["scores"]["dip"]
        distribution_coefficient = gene_scores["distributions"]["coefficient"]
        distribution_dip = gene_scores["distributions"]["dip"]
        # Calculate p-values.
        # These scores of bimodality indicate greater bimodality as values
        # increase.
        # The scores are unidirectional, so the hypothesis is unidirectional.
        # The p-value is the probability of obtaining by random chance a value
        # equal to or greater than the actual score.
        probability_coefficient = calculate_probability_equal_greater(
            value=score_coefficient,
            distribution=distribution_coefficient
        )
        probability_dip = calculate_probability_equal_greater(
            value=score_dip,
            distribution=distribution_dip
        )
        # Compile information.
        genes_p_values[gene]["coefficient"] = probability_coefficient
        genes_p_values[gene]["dip"] = probability_dip
    # Return information.
    return genes_p_values



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
    utility.confirm_path_directory(path_organization)
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

    # Read source information from file.
    source = read_source(dock=dock)

    # Check and summarize information about genes.
    check_genes_scores_distributions(
        genes=source["genes"],
        genes_scores_distributions=source["genes_scores_distributions"]
    )

    # Calculate probability of each gene's scores by chance.
    genes_p_values = calculate_p_values_genes(
        genes_scores_distributions=source["genes_scores_distributions"]
    )

    print(source["genes_scores_distributions"]["ENSG00000005483"])
    print(genes_p_values["ENSG00000005483"])

    # Collect names of genes and organize a data frame of genes.

    if False:

        print(source["data_gene_annotation"].iloc[0:10, 0:10])

        print(source["data_gene_signal_aggregation"].iloc[0:10, 0:10])
        #print(data_sum_index.loc[:, "ENSG00000240453.1"])
        # Gene ENSG00000240453.1 has values of 0.0 for all patients.

        # Calculate metrics of modality.
        series_bimodality = source["data_gene_signal_aggregation"].aggregate(metric.calculate_bimodality_coefficient, axis="index")
        data_bimodality_sparse = series_bimodality.to_frame(name="value").sort_values("value", axis="index", ascending=False)
        data_bimodality = collect_genes_names(
            data_gene_annotation=source["data_gene_annotation"],
            data_gene_score=data_bimodality_sparse
        )
        print(data_bimodality.iloc[0:50, : ])
        print(data_bimodality.shape)
        utility.print_terminal_partition(level=3)
        print("Tapasin... TAPBP...")
        print(data_bimodality.loc["ENSG00000231925"])


        #ENSG00000221947
        print(data_bimodality.loc["ENSG00000221947"])


        # Summarize table of patients' attributes.
        utility.print_terminal_partition(level=2)


        series_dip = source["data_gene_signal_aggregation"].aggregate(metric.calculate_dip_statistic, axis="index")
        data_dip_sparse = series_dip.to_frame(name="value").sort_values("value", axis="index", ascending=False)
        data_dip = collect_genes_names(
            data_gene_annotation=source["data_gene_annotation"],
            data_gene_score=data_dip_sparse
        )
        print(data_dip.iloc[0:50, : ])
        print(data_dip.shape)
        utility.print_terminal_partition(level=3)
        print("Tapasin... TAPBP...")
        print(data_dip.loc["ENSG00000231925"])





    if False:

        data_gene_signal_imputation_index = (
            source["data_gene_signal_imputation"].set_index(
                ["patient", "tissue"], append=False, drop=True
            )
        )
        print(data_gene_signal_imputation_index.iloc[0:10, 0:10])
        if False:
            data_gene_signal_sort = (
                data_gene_signal.sort_index(axis="columns", ascending=False)
            )

        # Split data by tissue.
        if False:
            data_gene_signal_imputation_index = (
                source["data_gene_signal_imputation"].set_index(
                ["tissue"], append=False, drop=True
                )
            )
            data_heart = (
                data_gene_signal_imputation_index
                    .loc["Heart"].reset_index(level=["tissue"])
            )
            print(data_heart)
            #data_merger = data_skin.merge(data_heart, how="outer")

        ########################################################




        # Reshape data with patients as columns and genes as rows.
        if False:
            utility.print_terminal_partition(level=2)
            print(
                "Reshape of data with patients as columns and genes as rows."
            )
            data_gene_signal_axis = data_gene_signal_sum.rename_axis("gene", axis="columns")
            data_gene_signal_axis_index = data_gene_signal_axis.set_index(
                ["patient"], append=False, drop=True
            )
            print(data_gene_signal_axis_index.iloc[0:10, 0:10])
            data_gene_signal_shape = data_gene_signal_axis_index.transpose(copy=True)
            print(data_gene_signal_shape.iloc[:, :])

            print("I will want to apply the bimodality test to each column of genes... columns need to be genes...")



        # Compile information.
        information = {}
        #Write product information to file.
        #write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
