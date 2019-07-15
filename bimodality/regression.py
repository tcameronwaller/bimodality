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
from sklearn.linear_model import LinearRegression

# Custom

import pipe
import metric
import plot
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source(
    gene=None,
    method=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of a gene
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_pipe = os.path.join(dock, "pipe")
    path_gene = os.path.join(path_pipe, gene)
    path_method = os.path.join(path_gene, method)

    path_report_restriction = os.path.join(
        path_method, "report_restriction.pickle"
    )
    path_report_distribution = os.path.join(
        path_method, "report_distribution.pickle"
    )

    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    with open(path_report_restriction, "rb") as file_source:
        report_restriction = pickle.load(file_source)
    with open(path_report_distribution, "rb") as file_source:
        report_distribution = pickle.load(file_source)

    # Access information.
    data_persons_tissues = report_restriction["data_persons_tissues"]
    data_gene_persons_signals = (
        report_distribution["data_gene_persons_signals"]
    )

    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_persons_tissues": data_persons_tissues,
        "data_gene_persons_signals": data_gene_persons_signals,
    }



def organize_data(
    data_samples_tissues_persons=None,
    data_persons_tissues=None,
    data_gene_persons_signals=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of a gene
        method (str): method for selection of tissues and patients, either
            "availability" for selection by minimal count of tissues, or
            "imputation" for selection by same tissues with imputation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    print(source["data_samples_tissues_persons"])
    print(source["data_persons_tissues"])
    print(source["data_gene_persons_signals"])


    # Specify directories and files.
    path_assembly = os.path.join(dock, "assembly")
    path_samples_tissues_persons = os.path.join(
        path_assembly, "data_samples_tissues_persons.pickle"
    )
    path_pipe = os.path.join(dock, "pipe")
    path_gene = os.path.join(path_pipe, gene)
    path_method = os.path.join(path_gene, method)

    path_report_restriction = os.path.join(
        path_method, "report_restriction.pickle"
    )
    path_report_distribution = os.path.join(
        path_method, "report_distribution.pickle"
    )

    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    with open(path_report_restriction, "rb") as file_source:
        report_restriction = pickle.load(file_source)
    with open(path_report_distribution, "rb") as file_source:
        report_distribution = pickle.load(file_source)

    # Access information.
    data_persons_tissues = report_restriction["data_persons_tissues"]
    data_gene_persons_signals = (
        report_distribution["data_gene_persons_signals"]
    )

    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_persons_tissues": data_persons_tissues,
        "data_gene_persons_signals": data_gene_persons_signals,
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
    utility.confirm_path_directory(path_analysis)
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
        delimiter="\t"
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
        information=information["genes_ensembl"],
        path_file=path_genes_ensembl
    )
    utility.write_file_text_list(
        information=information["genes_hugo"],
        path_file=path_genes_hugo
    )
    utility.write_file_text_list(
        information=information["genes_report"],
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
    utility.confirm_path_directory(path_analysis)
    path_reports = os.path.join(path_analysis, "reports")
    # Remove previous files since they change from run to run.
    utility.remove_directory(path=path_reports)
    utility.confirm_path_directory(path_reports)
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

    # Specify genes of interest.
    gene = "ENSG00000231925" # TAPBP
    #gene = "ENSG00000183878" # UTY ... Y-linked analog to KDM6A
    #gene = "ENSG00000147050" # KDM6A ... X-linked analog to UTY

    # Read source information from file.
    source = read_source(
        gene=gene,
        method="availability",
        dock=dock
    )

    # Organize persons' signals and properties.

    print(source["data_samples_tissues_persons"])
    print(source["data_persons_tissues"])
    print(source["data_gene_persons_signals"])

    persons = source["data_gene_persons_signals"].index.to_list()
    print("persons")
    data_samples_tissues_persons = source["data_samples_tissues_persons"].loc[
        source["data_samples_tissues_persons"]["person"].isin(persons),
    ]
    print(data_samples_tissues_persons)
    data_samples_tissues_persons.reset_index(
        level="sample",
        inplace=True
    )
    data_samples_tissues_persons.drop(
        labels=["sample", "tissue_major", "tissue_minor"],
        axis="columns",
        inplace=True
    )
    data_samples_tissues_persons.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_samples_tissues_persons.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )
    print(data_samples_tissues_persons)

    # Translate tissue
    def translate_tissue(value):
        if value:
            return 1
        else:
            return 0
    data_persons_tissues = source["data_persons_tissues"].applymap(
        translate_tissue,
    )
    print(data_persons_tissues)




    # Join
    data_persons_signals = source["data_gene_persons_signals"].join(
        data_samples_tissues_persons,
        how="left",
        on="person"
    )
    print(data_persons_signals)

    data_persons_tissues_signals = data_persons_signals.join(
        data_persons_tissues,
        how="left",
        on="person"
    )
    print(data_persons_tissues_signals)

    # Translate sex.
    def translate_sex(sex):
        if sex == "female":
            return 2
        elif sex == "male":
            return 1
    data_persons_tissues_signals["sex"] = (
        data_persons_tissues_signals["sex"].apply(translate_sex)
    )
    # Translate age.
    def translate_age(age):
        if age == "20-29":
            return 25
        elif age == "30-39":
            return 35
        elif age == "40-49":
            return 45
        elif age == "50-59":
            return 55
        elif age == "60-69":
            return 65
        elif age == "70-79":
            return 75
        elif age == "80-89":
            return 85
        elif age == "90-99":
            return 95
        elif age == "100-109":
            return 105
    data_persons_tissues_signals["age"] = (
        data_persons_tissues_signals["age"].apply(translate_age)
    )
    print(data_persons_tissues_signals)





    if True:
        data_independence = data_persons_tissues_signals.drop(
            labels=["value", "age", "sex"],
            axis="columns",
            inplace=False
        )
        pass
    if False:
        data_independence = data_persons_tissues_signals.loc[
            :, "age"
        ]
        pass
    #print(data_independence.values.reshape(-1, 1))
    data_dependence = data_persons_tissues_signals.loc[
        :, "value"
    ]
    #print(data_dependence.values.reshape(-1, 1))

    # Regression...
    regression = LinearRegression(
        fit_intercept=True,
        normalize=False,
        copy_X=True,
        n_jobs=None
    ).fit(
        data_independence.values,#.reshape(-1, 1),
        data_dependence.values.reshape(-1, 1),
        sample_weight=None
    )
    score = regression.score(
        data_independence.values,#.reshape(-1, 1),
        data_dependence.values.reshape(-1, 1),
        sample_weight=None
    )
    print(score)




    if False:

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
