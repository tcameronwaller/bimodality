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
    path_organization = os.path.join(dock, "organization")
    path_mean = os.path.join(
        path_organization, "data_gene_signal_mean.pickle"
    )
    path_median = os.path.join(
        path_organization, "data_gene_signal_median.pickle"
    )
    path_imputation = os.path.join(
        path_organization, "data_gene_signal_imputation.pickle"
    )
    path_log = os.path.join(
        path_organization, "data_gene_signal_log.pickle"
    )
    path_standard = os.path.join(
        path_organization, "data_gene_signal_standard.pickle"
    )
    path_aggregation = os.path.join(dock, "aggregation")
    path_signal = os.path.join(
        path_aggregation, "data_gene_signal_aggregation.pickle"
    )
    # Read information from file.
    data_gene_annotation = pandas.read_pickle(path_gene_annotation)
    data_gene_signal_mean = pandas.read_pickle(path_mean)
    data_gene_signal_median = pandas.read_pickle(path_median)
    data_gene_signal_imputation = pandas.read_pickle(path_imputation)
    data_gene_signal_log = pandas.read_pickle(path_log)
    data_gene_signal_standard = pandas.read_pickle(path_standard)
    data_gene_signal_aggregation = pandas.read_pickle(path_signal)
    # Compile and return information.
    return {
        "data_gene_annotation": data_gene_annotation,
        "data_gene_signal_mean": data_gene_signal_mean,
        "data_gene_signal_median": data_gene_signal_median,
        "data_gene_signal_imputation": data_gene_signal_imputation,
        "data_gene_signal_log": data_gene_signal_log,
        "data_gene_signal_standard": data_gene_signal_standard,
        "data_gene_signal_aggregation": data_gene_signal_aggregation,
    }


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
    # Summary.

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
