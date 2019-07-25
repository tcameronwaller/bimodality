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

import organization
import pipe
import category
import utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Source


def read_source_initial(
    source_genes=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        source_genes (str): name of directory from which to obtain genes list
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    if source_genes == "split":
        path_source = os.path.join(dock, "split")
    elif source_genes == "combination":
        path_source = os.path.join(dock, "combination")
    path_genes = os.path.join(
        path_source, "genes.txt"
    )
    # Read information from file.
    genes = utility.read_file_text_list(path_genes)
    # Compile and return information.
    return {
        "genes": genes,
    }


def read_source(
    gene=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        gene (str): identifier of single gene for which to execute the process.
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
    path_split = os.path.join(dock, "split")
    path_collection = os.path.join(path_split, "collection")
    path_gene = os.path.join(path_collection, (gene + ".pickle"))
    # Read information from file.
    data_samples_tissues_persons = pandas.read_pickle(
        path_samples_tissues_persons
    )
    data_gene_samples_signals = pandas.read_pickle(path_gene)
    # Compile and return information.
    return {
        "data_samples_tissues_persons": data_samples_tissues_persons,
        "data_gene_samples_signals": data_gene_samples_signals,
    }


##########
# Candidacy


def evaluate_gene_population(
    threshold=None,
    data_gene_persons_signals=None
):
    """
    Evaluates the a gene's representation across persons.

    arguments:
        threshold (int): minimal count of persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene has representation across an adequate
            population of persons

    """

    count_persons = data_gene_persons_signals.shape[0]
    return (count_persons > threshold)


def evaluate_gene_signal(
    threshold=None,
    data_gene_persons_signals=None
):
    """
    Evaluates the gene's aggregate tissue signals across persons.

    arguments:
        threshold (float): minimal standard deviation
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene's signals are adequate for analysis

    """

    # Calculate standard deviation of values.
    deviation = numpy.std(data_gene_persons_signals["value"].values, axis=0)
    return deviation > threshold


# TODO: this function will eventually need to accommodate additional attributes of persons...

def evaluate_gene_tissue(
    threshold_proportion=None,
    threshold_probability=None,
    data_samples_tissues_persons=None,
    data_gene_persons_tissues=None,
    data_gene_persons_signals=None,
):
    """
    Evaluates the selection of tissues for calculation of gene's aggregate
    scores.

    arguments:
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a tissue's group in analysis
        threshold_probability (float): minimal probability for candidacy
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples
        data_gene_persons_tissues (object): Pandas data frame of a gene's
            selection of tissues across persons
        data_gene_persons_signals (object): Pandas data frame of a gene's
            aggregate, pan-tissue signals across persons

    raises:

    returns:
        (bool): whether the gene's selection of tissues is acceptable

    """

    # Organize information about gene's scores across persons and tissues.
    data_persons_tissues_signals = category.organize_gene_persons_signals(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=data_gene_persons_tissues,
        data_gene_persons_signals=data_gene_persons_signals,
    )

    # Determine minimal count of persons for inclusion of a tissue's group.
    count_persons = data_gene_persons_signals.shape[0]
    threshold_count = threshold_proportion * count_persons

    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    #print(data_persons_tissues_signals)
    tissues = data_persons_tissues_signals.columns.to_list()
    tissues.remove("value")
    tissues.remove("female")
    tissues.remove("age")
    report = category.evaluate_variance_by_binary_groups(
        groups=tissues,
        threshold=threshold_count,
        data=data_persons_tissues_signals,
    )

    # Evaluate whether gene's scores differ by tissue selection.
    if (
        (math.isnan(report["probability"])) or
        (report["probability"] >= threshold_probability)
    ):
        return True
    else:
        return False


def evaluate_gene_candidacy_method(
    gene=None,
    method=None,
    data_samples_tissues_persons=None,
    observation=None,
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
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

    raises:

    returns:
        (dict<bool>): information about gene's candidacy

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
    data_gene_samples_signals=None,
    observation=None,
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
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

    raises:

    returns:
        (dict<dict<bool>>): information about gene's candidacy

    """

    imputation = evaluate_gene_candidacy_method(
        gene=gene,
        method="imputation",
        data_samples_tissues_persons=data_samples_tissues_persons,
        observation=observation,
    )

    availability = evaluate_gene_candidacy_method(
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


def evaluate_gene_categories(
    gene=None,
    threshold_proportion=None,
    method=None,
    data_samples_tissues_persons=None,
    observation=None,
):
    """
    Evaluates the selection of tissues for calculation of gene's aggregate
    scores.

    arguments:
        gene (str): identifier of single gene for which to execute the process
        threshold_proportion (float): minimal proportion of persons for
            inclusion of a group in analyses
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
        (bool): whether the gene's selection of tissues is acceptable

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

    # Calculate tissue-tissue correlations across persons.
    # Use the data before imputation to grasp correlations across real values.
    data_gene_mean_correlations = (
        category.calculate_mean_tissue_pairwise_correlations(
            data_gene_persons_tissues_signals=(
                report_restriction["data_gene_persons_tissues_signals"]
            )
        )
    )

    # Organize information about gene's scores across persons and tissues.
    data_persons_signals_groups = category.organize_persons_signals_groups(
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_persons_tissues=data_gene_persons_tissues,
        data_gene_persons_signals=data_gene_persons_signals,
        data_gene_mean_correlations=data_gene_mean_correlations,
    )

    #print(data_persons_signals_groups)

    # Determine minimal count of persons for inclusion of a group in
    # analyses.
    count_persons = data_persons_signals_groups.shape[0]
    threshold_count = threshold_proportion * count_persons

    # Analysis of Variance (ANOVA)

    # Signals versus tissues
    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    tissues = data_gene_persons_tissues.columns.values.tolist()
    report = category.evaluate_variance_by_binary_groups(
        groups=tissues,
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )
    anova_tissue_p = report["probability"]

    # Signals versus sex
    # Determine whether tissue selection confounds gene's aggregate scores
    # across persons.
    sexes = ["male", "female"]
    report = category.evaluate_variance_by_binary_groups(
        groups=sexes,
        threshold=threshold_count,
        data=data_persons_signals_groups,
    )
    anova_sex_p = report["probability"]

    # Regression

    # Signals versus sex
    regression_sex_r = category.evaluate_correlation_by_linearity(
        dependence="value",
        independence=["female"],
        data=data_persons_signals_groups,
    )

    # Signals versus age
    regression_age_r = category.evaluate_correlation_by_linearity(
        dependence="value",
        independence=["age"],
        data=data_persons_signals_groups,
    )

    # Signals versus tissues
    regression_tissue_r = category.evaluate_correlation_by_linearity(
        dependence="value",
        independence=tissues,
        data=data_persons_signals_groups,
    )

    # Signals versus mean tissue-tissue correlation
    regression_correlation_r = category.evaluate_correlation_by_linearity(
        dependence="value",
        independence=["correlation"],
        data=data_persons_signals_groups,
    )

    # Compile information.
    information = {
        "anova_tissue_p": anova_tissue_p,
        "anova_sex_p": anova_sex_p,
        "regression_sex_r": regression_sex_r,
        "regression_age_r": regression_age_r,
        "regression_tissue_r": regression_tissue_r,
        "regression_correlation_r": regression_correlation_r,
    }

    # Return information.
    return information


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

    # Evaluate gene's distribution between categories.
    categories = evaluate_gene_categories(
        gene=gene,
        threshold_proportion=0.05,
        data_samples_tissues_persons=data_samples_tissues_persons,
        method=method,
        observation=observation,
    )

    # Compile report.
    # sex anova p-val
    # tissue anova p-val
    # tissue regression metric?
    # age regression metric?
    # tissue-specific mean, variance, dispersion, etc...
    report = dict()
    report["gene"] = gene
    #report["persons_restriction"] = report_restriction["persons"]
    report["persons"] = report_distribution["persons"]
    report["tissues_mean"] = report_restriction["tissues_mean"]
    report["tissues_median"] = report_restriction["tissues_median"]

    report["anova_sex_p"] = categories["anova_sex_p"]
    report["anova_tissue_p"] = categories["anova_tissue_p"]

    report["regression_sex_r"] = categories["regression_sex_r"]
    report["regression_age_r"] = categories["regression_age_r"]
    report["regression_tissue_r"] = categories["regression_tissue_r"]
    report["regression_correlation_r"] = categories["regression_correlation_r"]

    # TODO: include metrics from the organization report...
    # TODO: include maximum variation and dispersion across tissues...
    # variation_max
    # dispersion_max

    print(report)


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
    observation=None,
):
    """
    Prepares a report about a gene's multiple factors.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples
        observation (dict): information about a gene's actual distribution of
            signals across persons and tissues

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


def evaluate_gene(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
):
    """
    Evaluates a gene by multiple factors.

    arguments:
        gene (str): identifier of a gene
        data_samples_tissues_persons (object): Pandas data frame of persons
            and tissues for all samples.
        data_gene_samples_signals (object): Pandas data frame of a gene's
            signals across samples

    raises:

    returns:
        (dict<dict<bool>>): information about gene's candidacy

    """

    # Determine gene's distributions of aggregate tissue scores across persons.
    observation = pipe.determine_gene_distributions(
        gene=gene,
        metric=False,
        correlation=False,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Report.
    report = prepare_gene_report(
        gene=gene,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
        observation=observation,
    )

    # Candidacy.
    candidacy = evaluate_gene_candidacy(
        gene=gene,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
        observation=observation,
    )

    # Compile information.
    information = {
        "report": report,
        "candidacy": candidacy,
    }

    # Return information.
    return information



# TODO: Compose new report for each gene...
# Gene ... count_persons ... mean_tissues_per_person ... tissue_ANOVA_p-val ... tissue_correlation_regression
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


##########
# Product


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

    utility.confirm_path_directory(path_candidacy)
    utility.confirm_path_directory(path_count)
    utility.confirm_path_directory(path_all)

    utility.confirm_path_directory(path_imputation)
    utility.confirm_path_directory(path_imputation_all)
    utility.confirm_path_directory(path_imputation_population)
    utility.confirm_path_directory(path_imputation_signal)
    utility.confirm_path_directory(path_imputation_tissue)

    utility.confirm_path_directory(path_availability)
    utility.confirm_path_directory(path_availability_all)
    utility.confirm_path_directory(path_availability_population)
    utility.confirm_path_directory(path_availability_signal)
    utility.confirm_path_directory(path_availability_tissue)

    pass


def write_product(
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
    #utility.confirm_path_directory(path_candidacy)
    path_method = os.path.join(path_candidacy, method)
    #utility.confirm_path_directory(path_method)

    if (method == "count") or (method == "all"):
        path_gene = os.path.join(path_method, (gene))
        pass
    else:
        path_type = os.path.join(path_method, type)
        #utility.confirm_path_directory(path_type)
        path_gene = os.path.join(path_type, (gene))
        pass

    # Write information to file.
    utility.write_file_text_list(
        information=["null"], path_file=path_gene
    )

    pass


def read_collect_write_genes(
    dock=None,
):
    """
    Collects information about genes.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files.

    raises:

    returns:
        (list<str>): information about genes

    """

    # Specify directories and files.
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

    path_genes_all = os.path.join(path_candidacy, "genes_all.txt")
    path_genes_imputation_all = os.path.join(
        path_candidacy, "genes_imputation_all.txt"
    )
    path_genes_availability_all = os.path.join(
        path_candidacy, "genes_availability_all.txt"
    )

    # Read genes.
    genes_all = os.listdir(path_all)

    genes_imputation_all = os.listdir(path_imputation_all)
    genes_imputation_population = os.listdir(path_imputation_population)
    genes_imputation_signal = os.listdir(path_imputation_signal)
    genes_imputation_tissue = os.listdir(path_imputation_tissue)

    genes_availability_all = os.listdir(path_availability_all)
    genes_availability_population = os.listdir(path_availability_population)
    genes_availability_signal = os.listdir(path_availability_signal)
    genes_availability_tissue = os.listdir(path_availability_tissue)

    utility.remove_directory(path=path_all)
    utility.remove_directory(path=path_count)
    utility.remove_directory(path=path_imputation)
    utility.remove_directory(path=path_availability)

    utility.print_terminal_partition(level=1)
    print("Counts of genes in candidacy categories...")

    print(
        "all (both imputation and availability): " +
        str(len(genes_all))
    )

    print(
        "imputation all (population, signal, and tissue): " +
        str(len(genes_imputation_all))
    )
    print(
        "imputation population: " +
        str(len(genes_imputation_population))
    )
    print(
        "imputation signal: " +
        str(len(genes_imputation_signal))
    )
    print(
        "imputation tissue: " +
        str(len(genes_imputation_tissue))
    )


    print(
        "availability all (population, signal, and tissue): " +
        str(len(genes_availability_all))
    )
    print(
        "availability population: " +
        str(len(genes_availability_population))
    )
    print(
        "availability signal: " +
        str(len(genes_availability_signal))
    )
    print(
        "availability tissue: " +
        str(len(genes_availability_tissue))
    )


    # Write information to file.
    utility.write_file_text_list(
        information=genes_all, path_file=path_genes_all
    )
    utility.write_file_text_list(
        information=genes_imputation_all, path_file=path_genes_imputation_all
    )
    utility.write_file_text_list(
        information=genes_availability_all,
        path_file=path_genes_availability_all,
    )

    pass


# TODO: need a new function to compile

def report_gene_tissue_persons(
    gene=None,
    count=None,
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        gene (str): identifier of a gene
        count (int): either minimal count of tissues for selection by
            "availability" or maximal count of imputation for selection by
            "imputation"
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )
    # Determine gene's distribution of aggregate tissue scores across persons.
    imputation = determine_gene_signal_distribution(
        gene=gene,
        method="imputation",
        count=count,
        data_gene_samples_signals=source["data_gene_samples_signals"],
    )
    availability = determine_gene_signal_distribution(
        gene=gene,
        method="availability",
        count=count,
        data_gene_samples_signals=source["data_gene_samples_signals"],
    )
    utility.print_terminal_partition(level=2)
    print("Persons by imputation and availability methods...")
    print(imputation["data_gene_persons_signals"].size)
    print(availability["data_gene_persons_signals"].size)

    pass



###############################################################################
# Procedure

def execute_procedure(
    gene=None,
    data_samples_tissues_persons=None,
    data_gene_samples_signals=None,
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
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # TODO: Make this a 2-part procedure...
    # 1. describe the gene
    # 2. assign the gene candidacy...

    # Evaluate gene's candidacy.
    collection = evaluate_gene(
        gene=gene,
        data_samples_tissues_persons=data_samples_tissues_persons,
        data_gene_samples_signals=data_gene_samples_signals,
    )

    # Determine gene's candidacy and write reports.
    report_gene_candidacy(
        gene=gene,
        candidacy=collection["candidacy"],
        dock=dock,
    )

    # Compile a report about the gene.
    # pull in information from candidacy, categories, organization, restriction, distribution, bimodality,


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

    # Remove previous files to avoid version or batch confusion.
    path_candidacy = os.path.join(dock, "candidacy")
    utility.remove_directory(path=path_candidacy)

    # Create directories.
    create_directories(dock=dock)

    # Read source information from file.
    source = read_source_initial(
        source_genes="split",
        dock=dock
    )
    print("count of genes: " + str(len(source["genes"])))

    # Report date and time.
    start = datetime.datetime.now()
    print(start)

    if True:
        execute_procedure_local_sub(
            gene="ENSG00000231925", # TAPBP
            dock=dock,
        )

    if False:
        # Set up partial function for iterative execution.
        # Each iteration uses the same values for "genes_signals", "shuffles", and
        # "dock" variables.
        execute_procedure_gene = functools.partial(
            execute_procedure_local_sub,
            dock=dock
        )

        # Initialize multiprocessing pool.
        #pool = multiprocessing.Pool(processes=os.cpu_count())
        pool = multiprocessing.Pool(processes=8)

        # Iterate on genes.
        check_genes=[
            "ENSG00000231925", # TAPBP
        ]
        report = pool.map(execute_procedure_gene, check_genes)
        #report = pool.map(execute_procedure_gene, source["genes"][0:8])
        #report = pool.map(execute_procedure_gene, source["genes"])

    # Pause procedure.
    time.sleep(5.0)

    # Report.
    #print("Process complete for the following genes...")
    #print(str(len(report)))

    if False:
        report_gene_tissue_persons(
            gene="ENSG00000231925", # TAPBP
            count=10,
            dock=dock,
        )

    # Collect genes.
    read_collect_write_genes(dock=dock)

    # Report date and time.
    end = datetime.datetime.now()
    print(end)
    print("duration: " + str(end - start))

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

    # Report gene.
    #print("gene: " + gene)

    # Read source information from file.
    source = read_source(
        gene=gene,
        dock=dock
    )

    # Execute procedure.
    execute_procedure(
        gene=gene,
        data_samples_tissues_persons=source["data_samples_tissues_persons"],
        data_gene_samples_signals=source["data_gene_samples_signals"],
        dock=dock
    )

    # Report contents of directory.
    path_candidacy = os.path.join(dock, "candidacy")
    path_count = os.path.join(path_candidacy, "count")
    files = os.listdir(path_count)
    count = len(files)
    if (count % 10 == 0):
        print("complete genes: " + str(len(files)))

    pass


if (__name__ == "__main__"):
    execute_procedure_local()
