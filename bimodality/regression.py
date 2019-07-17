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


##########
# Gene's persons and signals


def organize_persons_attributes(
    persons=None,
    data_samples_tissues_persons=None,
):
    """
    Organizes persons' attributes.

    arguments:
        persons (list<str>): list of persons of interest
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples

    raises:

    returns:
        (object): Pandas data frame of persons and tissues for all samples

    """

    # Define functions
    def translate_sex(sex):
        if sex == "female":
            return 1
        elif sex == "male":
            return 0

    def translate_age(age):
        ages_text = [
            "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89",
            "90-99", "100-109",
        ]
        ages_number = [
            25, 35, 45, 55, 65, 75, 85, 95, 105
        ]
        # Match indices.
        if age in ages_text:
            return ages_number[ages_text.index(age)]
        else:
            return float("nan")

    # Copy data.
    data_persons = data_samples_tissues_persons.copy(deep=True)

    # Select persons and attributes of interest.
    data_persons.reset_index(
        level="sample",
        inplace=True
    )
    data_persons = data_persons.loc[
        data_persons["person"].isin(persons),
    ]
    data_persons.drop(
        labels=["sample", "tissue_major", "tissue_minor"],
        axis="columns",
        inplace=True
    )
    data_persons.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_persons.set_index(
        ["person"],
        append=False,
        drop=True,
        inplace=True
    )

    # Organize persons' sex and age.
    data_persons["female"] = (
        data_persons["sex"].apply(translate_sex)
    )
    data_persons["age"] = (
        data_persons["age_range"].apply(translate_age)
    )
    data_persons.drop(
        labels=["sex", "age_range"],
        axis="columns",
        inplace=True
    )
    # Return data.
    return data_persons


def organize_gene_persons_signals(
    data_samples_tissues_persons=None,
    data_gene_persons_tissues=None,
    data_gene_persons_signals=None
):
    """
    Reads and organizes source information from file

    arguments:
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate signals across persons

    raises:

    returns:
        (object): Pandas data frame of a gene's aggregate signals across
            persons with persons's attributes

    """

    # Select persons of interest.
    persons = data_gene_persons_signals.index.to_list()

    # Organize information about persons.
    data_persons_attributes = organize_persons_attributes(
        persons=persons,
        data_samples_tissues_persons=data_samples_tissues_persons,
    )

    # Organize information about persons' tissues.

    # Join
    data_persons_signals = data_gene_persons_signals.join(
        data_persons_attributes,
        how="left",
        on="person"
    )
    data_persons_tissues_signals = data_persons_signals.join(
        data_gene_persons_tissues,
        how="left",
        on="person"
    )

    # Return information
    return data_persons_tissues_signals


def select_values_by_binary_category(
    category=None,
    value=None,
    data=None
):
    """
    Selects a gene's aggregate scores by binary category.

    arguments:
        category (str): name of a category
        value (int): value of a category
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (list<float>): values for a group

    """

    # Select values of gene's scores by attribute category.
    data_selection = data.loc[
        data[category] == value, :
    ]
    return data_selection["value"].to_list()


def evaluate_variance_by_binary_groups(
    groups=None,
    data=None
):
    """
    Evaluates the variance between values in binary groups by one-way analysis
    of variance.

    The probability is sometimes a missing value.
    This missing value does not seem to indicate missing values in the groups.
    Potential problems
        1. a specific group lacks values
        2. multiple groups have the same values

    arguments:
        groups (list<str>): names of binary groups
        data (object): Pandas data frame of values across binary groups

    raises:

    returns:
        (dict<float>): values of results from the analysis of variance

    """

    # Collect values in each group.
    groups_values = list()
    for group in groups:
        values = select_values_by_binary_category(
            category=group,
            value=1,
            data=data,
        )
        if len(values) > 2:
            groups_values.append(values)
    # Analyze variance.
    # Use splat operator to extract arguments from list.
    # Alternative: scipy.stats.f_oneway(*[df[col] for col in df])
    # The probability is sometimes a missing value.
    # Even when this probability is missing, there do not seem to be
    # irregularities in the values.
    if (len(groups_values) > 1):
        report = scipy.stats.f_oneway(*groups_values)
        # Compile information.
        information = {
            "statistic": report[0],
            "probability": report[1],
        }
    else:
        # Compile information.
        information = {
            "statistic": float("nan"),
            "probability": float("nan"),
        }
    # Return information.
    return information




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
