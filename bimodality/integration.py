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

# TODO: "integration" procedure should be about integrating ranks of genes from modality with heritability analysis


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
    path_combination = os.path.join(dock, "combination")
    path_distributions = os.path.join(
        path_combination, "genes_scores_distributions.pickle"
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



# TODO: write_product_reports()


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


    #############################################################

    # Calculate genes' probabilities of bimodality.
    data_genes_probabilities_availability = rank.calculate_probabilities_genes(
        genes_scores_permutations=(
            #source["genes_scores_permutations_availability"]
            genes_scores_permutations_availability
        ),
    )
    data_genes_probabilities_imputation = rank.calculate_probabilities_genes(
        genes_scores_permutations=(
            #source["genes_scores_permutations_imputation"]
            genes_scores_permutations_imputation
        ),
    )
    print(data_genes_probabilities_availability)
    print(data_genes_probabilities_imputation)

    # Rank genes by probabilities.
    data_ranks_availability = rank.rank_genes(
        data_genes_probabilities=data_genes_probabilities_availability,
        rank="combination",
        method="threshold",
        threshold=0.001,
        count=1000,
    )
    data_ranks_imputation = rank.rank_genes(
        data_genes_probabilities=data_genes_probabilities_imputation,
        rank="combination",
        method="threshold",
        threshold=0.001,
        count=1000,
    )
    print(data_ranks_availability)
    print(data_ranks_imputation)

    # Extract identifiers of genes.
    genes_availability = data_ranks_availability.index.to_list()
    genes_imputation = data_ranks_imputation.index.to_list()

    # Combine consensus genes.
    genes_combination = utility.filter_unique_union_elements(
        list_one=genes_availability,
        list_two=genes_imputation,
    )
    print("Count of consensus genes: " + str(len(genes_combination)))

    ############################################################

    # Calculate probability of each gene's scores by chance.
    data_genes_probabilities = calculate_probabilities_genes(
        genes_scores_distributions=source["genes_scores_distributions"]
    )
    #print(source["genes_scores_distributions"]["ENSG00000005483"])
    #print(genes_probabilities["ENSG00000005483"])

    # Rank genes for subsequent functional analyses.
    ranks = rank_genes(
        data_genes_probabilities=data_genes_probabilities,
        rank="combination",
        method="threshold",
        threshold=0.05,
        count=500,
    )
    print(ranks["data_rank_genes"])



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
    #Write product information to file.
    write_product(dock=dock, information=information)

    pass


if (__name__ == "__main__"):
    execute_procedure()
