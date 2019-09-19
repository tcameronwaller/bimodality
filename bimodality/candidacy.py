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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom

import distribution
import organization
import category
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Source



##########
# Candidacy


def evaluate_gene_candidacy_method(
    gene=None,
    method=None,
    data_samples_tissues_persons=None,
    data_gene_persons_signals=None,
):
    """
    Evaluates a gene's candidacy for further bimodal analysis by multiple
    criteria.

    1. candidacy by population
    - gene must have tissue-aggregate signals across at least 100 persons
    2. candidacy by signal
    - gene's tissue-aggregate signals across persons must have a standard
        deviation greater than 0
    3. candidacy by tissue independence
    - consider tissues with representation in at least 5% of persons

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict<bool>): information about gene's candidacy

    """

    # Determine whether the gene's distribution is adequate.
    # Distribution procedure applies filters to the gene's tissue-aggregate
    # scores across persons.

    # Determine whether gene has representation across adequate persons.
    population = evaluate_gene_population(
        threshold=100,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine whether gene's tissue-aggregate scores across persons have
    # adequate variance to apply metrics to the distribution.
    signal = evaluate_gene_signal(
        threshold=0.0,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine whether the gene's distribution is dependent on tissue
    # selection.
    if False:
        tissue = evaluate_gene_tissue(
            threshold_proportion=0.1,
            threshold_probability=0.05,
            data_samples_tissues_persons=data_samples_tissues_persons,
            data_gene_persons_tissues=data_gene_persons_tissues,
            data_gene_persons_signals=data_gene_persons_signals,
        )

    # Compile information.
    information = {
        "population": population,
        "signal": signal,
        "tissue": True,
    }

    # Return information.
    return information


def evaluate_gene_candidacy(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_persons_signals=None,
):
    """
    Evaluates a gene's candidacy for further bimodal analysis by multiple
    criteria.

    1. candidacy by population
    - gene must have tissue-aggregate signals across at least 100 persons
    2. candidacy by signal
    - gene's tissue-aggregate signals across persons must have a standard
        deviation greater than 0
    3. candidacy by tissue independence
    - consider tissues with representation in at least 5% of persons

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (dict<dict<bool>>): information about gene's candidacy

    """

    imputation = evaluate_gene_candidacy_method(
        gene=gene,
        method="imputation",
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    availability = evaluate_gene_candidacy_method(
        gene=gene,
        method="availability",
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Compile information.
    information = {
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


##########
# Report


def prepare_gene_report_method(
    gene=None,
    method=None,
    data_samples_tissues_persons=None,
    observation=None,
):
    """
    Prepares a report about a gene's multiple factors.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): method for selection of tissues and persons in
            restriction procedure, either "imputation" for selection by
            specific tissues with imputation or "availability" for selection by
            minimal count of tissues
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

    raises:

    returns:
        (dict): information about gene's candidacy

    """

    # Access information.
    report_organization = observation["organization"]["report_gene"]
    report_restriction = observation[method]["report_restriction"]
    report_distribution = observation[method]["report_distribution"]
    data_gene_persons_tissues = (
        observation[method]["report_restriction"]["data_gene_persons_tissues"]
    )
    data_gene_persons_signals = (
        observation[method]["report_distribution"]["data_gene_persons_signals"]
    )

    # Compile report.
    report = dict()
    report["gene"] = gene
    #report["persons_restriction"] = report_restriction["persons"]
    report["persons"] = report_distribution["persons"]
    report["tissues_mean"] = report_restriction["tissues_mean"]
    report["tissues_median"] = report_restriction["tissues_median"]
    # Compile information.
    information = {
        "report": report,
    }

    # Return information.
    return information


def prepare_gene_report(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
    candidacy=None,
):
    """
    Prepares a report about a gene's multiple factors.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        candidacy (dict<dict<bool>>): information about gene's candidacy

    raises:

    returns:
        (dict): information about gene's candidacy

    """

    imputation = prepare_gene_report_method(
        gene=gene,
        method="imputation",
        data_samples_tissues_persons=data_samples_tissues_persons,
        observation=observation,
    )

    availability = prepare_gene_report_method(
        gene=gene,
        method="availability",
        data_samples_tissues_persons=data_samples_tissues_persons,
        observation=observation,
    )

    # Compile information.
    information = {
        "imputation": imputation,
        "availability": availability,
    }

    # Return information.
    return information


##########
# Product



########################Scrap#################################
# TODO: current "report_gene_candidacy" isn't very useful...
# TODO: Compose new report for each gene...
# Gene ... count_persons ... mean_tissues_per_person ... tissue_ANOVA_p-val ... tissue_correlation_regression
# Save a dictionary for each gene... then compile later...
# ...
def report_gene_candidacy(
    gene=None,
    candidacy=None,
    dock=None,
):
    """
    Reports a gene's candidacy.

    arguments:
        gene (str): identifier of a gene
        candidacy (dict<dict<bool>>): information about gene's candidacy
            and tissues for all samples.
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:


    """

    # Write gene to dummy list to keep track.
    write_product(
        gene=gene,
        method="count",
        type="null",
        dock=dock,
    )

    # All
    if (
        candidacy["imputation"]["population"] and
        candidacy["imputation"]["signal"] and
        candidacy["imputation"]["tissue"] and
        candidacy["availability"]["population"] and
        candidacy["availability"]["signal"] and
        candidacy["availability"]["tissue"]
    ):
        write_product(
            gene=gene,
            method="all",
            type="null",
            dock=dock,
        )

    # Imputation
    if (
        candidacy["imputation"]["population"] and
        candidacy["imputation"]["signal"] and
        candidacy["imputation"]["tissue"]
    ):
        write_product(
            gene=gene,
            method="imputation",
            type="all",
            dock=dock,
        )
    # Population
    if candidacy["imputation"]["population"]:
        write_product(
            gene=gene,
            method="imputation",
            type="population",
            dock=dock,
        )
    # Signal
    if candidacy["imputation"]["signal"]:
        write_product(
            gene=gene,
            method="imputation",
            type="signal",
            dock=dock,
        )
    # Tissue
    if candidacy["imputation"]["tissue"]:
        write_product(
            gene=gene,
            method="imputation",
            type="tissue",
            dock=dock,
        )

    # Availability
    if (
        candidacy["availability"]["population"] and
        candidacy["availability"]["signal"] and
        candidacy["availability"]["tissue"]
    ):
        write_product(
            gene=gene,
            method="availability",
            type="all",
            dock=dock,
        )
    # Population
    if candidacy["availability"]["population"]:
        write_product(
            gene=gene,
            method="availability",
            type="population",
            dock=dock,
        )
    # Signal
    if candidacy["availability"]["signal"]:
        write_product(
            gene=gene,
            method="availability",
            type="signal",
            dock=dock,
        )
    # Tissue
    if candidacy["availability"]["tissue"]:
        write_product(
            gene=gene,
            method="availability",
            type="tissue",
            dock=dock,
        )

    pass

def write_product_scrap(
    gene=None,
    method=None,
    type=None,
    dock=None,
):
    """
    Writes product information to file.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        method (str): list to which to append the gene, "all", "imputation", or
            "availability"
        type (str): type of criterion for candidacy
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Specify directories and files.
    path_candidacy = os.path.join(dock, "candidacy")
    #utility.create_directory(path_candidacy)
    path_method = os.path.join(path_candidacy, method)
    #utility.create_directory(path_method)

    if (method == "count") or (method == "all"):
        path_gene = os.path.join(path_method, (gene))
        pass
    else:
        path_type = os.path.join(path_method, type)
        #utility.create_directory(path_type)
        path_gene = os.path.join(path_type, (gene))
        pass

    # Write information to file.
    utility.write_file_text_list(
        elements=["null"],
        delimiter="\n",
        path_file=path_gene
    )

    pass

def create_directories(
    dock=None,
):
    """
    Creates directories.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:

    """

    # Create directories.
    path_candidacy = os.path.join(dock, "candidacy")
    path_count = os.path.join(path_candidacy, "count")
    path_all = os.path.join(path_candidacy, "all")

    path_imputation = os.path.join(path_candidacy, "imputation")
    path_imputation_all = os.path.join(path_imputation, "all")
    path_imputation_population = os.path.join(path_imputation, "population")
    path_imputation_signal = os.path.join(path_imputation, "signal")
    path_imputation_tissue = os.path.join(path_imputation, "tissue")

    path_availability = os.path.join(path_candidacy, "availability")
    path_availability_all = os.path.join(path_availability, "all")
    path_availability_population = os.path.join(path_availability, "population")
    path_availability_signal = os.path.join(path_availability, "signal")
    path_availability_tissue = os.path.join(path_availability, "tissue")

    utility.create_directory(path_candidacy)
    utility.create_directory(path_count)
    utility.create_directory(path_all)

    utility.create_directory(path_imputation)
    utility.create_directory(path_imputation_all)
    utility.create_directory(path_imputation_population)
    utility.create_directory(path_imputation_signal)
    utility.create_directory(path_imputation_tissue)

    utility.create_directory(path_availability)
    utility.create_directory(path_availability_all)
    utility.create_directory(path_availability_population)
    utility.create_directory(path_availability_signal)
    utility.create_directory(path_availability_tissue)

    pass
######################Scrap##################################

###############################################################################
# Procedure

def execute_procedure(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
    data_families=None,
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        data_families (object): Pandas data frame of families and persons
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # TODO: I need an export table of each gene's pan-tissue signals across persons
    # 3. evaluate gene's candidacy
    # 3.1. report gene's summary info (including candidacy) within a table....

    # Determine gene's distributions of aggregate tissue signals across
    # persons.
    observation = distribution.determine_gene_distributions(
        gene=gene,
        modality=False,
        data_gene_samples_signals=data_gene_samples_signals,
    )
    report_restriction = observation["availability"]["report_restriction"]
    report_aggregation = observation["availability"]["report_aggregation"]
    data_gene_persons_signals = report_aggregation["data_gene_persons_signals"]

    # Export gene's pan-tissue signals across persons.
    # Design format compatible with Genome-wide Complex Trait Analysis (GCTA).
    data_gene_families_persons_signals = export_gene_persons_signals(
        data_gene_persons_signals=data_gene_persons_signals,
        data_families=data_families,
    )

    # Evaluate gene's candidacy.
    candidacy = evaluate_gene_candidacy(
        gene=gene,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Compile report about gene's distribution.
    # pull in information from candidacy, categories, organization, restriction, distribution, bimodality,
    # TODO: consider including name and chromosome in this report...
    if False:
        report = prepare_gene_report(
            gene=gene,
            candidacy=candidacy,
            persons=report_aggregation["persons"],
            tissues_mean=report_restriction["tissues_mean"],
            tissues_median=report_restriction["tissues_median"],
        )

    # Compile information.
    information = {
        #"report": report,
        "data_gene_families_persons_signals": (
            data_gene_families_persons_signals
        ),
    }
    #Write product information to file.
    write_product_gene(dock=dock, gene=gene, information=information)

    pass


def execute_procedure_local(dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Collect genes.
    #read_collect_write_genes(dock=dock)

    pass


def execute_procedure_local_sub(gene=None, dock=None):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files
        gene (str): identifier of single gene for which to execute the process.

    raises:

    returns:

    """

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        data_families=source["data_families"],
        dock=dock
    )

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
