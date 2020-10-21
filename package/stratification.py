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
import copy
import random
import itertools
import functools
import multiprocessing
import datetime
import time


# Relevant

import numpy
import pandas
import scipy.stats
#from sklearn.linear_model import LinearRegression
import statsmodels.api
import statsmodels.stats.outliers_influence

# Custom

import integration
import selection
import distribution
import modality
import prediction
import utility
import assembly

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(dock=None):
    """
    Initialize directories for procedure's product files.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = dock
    paths["stratification"] = os.path.join(paths["dock"], "stratification")
    # Remove previous files to avoid version or batch confusion.
    utility.remove_directory(path=paths["stratification"])
    utility.create_directory(path=paths["stratification"])

    # Define paths for cohorts of persons.
    cohorts = list()
    cohorts.append("selection")
    cohorts.append("respiration")
    cohorts.append("ventilation")
    for cohort in cohorts:
        paths[cohort] = dict()
        paths[cohort]["component"] = os.path.join(
            paths["stratification"], cohort, "component"
        )
        paths[cohort]["regression"] = os.path.join(
            paths["stratification"], cohort, "regression"
        )
        paths[cohort]["summary"] = os.path.join(
            paths["stratification"], cohort, "summary"
        )
        # Initialize directories.
        utility.create_directories(path=paths[cohort]["component"])
        utility.create_directories(path=paths[cohort]["regression"])
        utility.create_directories(path=paths[cohort]["summary"])
    # Return information.
    return paths


##########
# Components by gene expression
# 1. selection cohort
# 2. respiration cohort
# 3. ventilation cohort


def read_source_cohort_gene_components(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information
    print(source["genes_function"].keys())

    """

    # Specify directories and files.
    path_data_signals_genes_persons = os.path.join(
        dock, "distribution", cohort, "collection",
        "data_signals_genes_persons.pickle"
    )
    # Read information from file.
    data_signals_genes_persons = pandas.read_pickle(
        path_data_signals_genes_persons
    )
    # Read genes sets.
    genes_candidacy = integration.read_source_genes_sets_candidacy(
        cohort=cohort,
        dock=dock,
    )
    # Return information.
    return {
        "data_signals_genes_persons": data_signals_genes_persons,
        "genes_candidacy": genes_candidacy,
    }


def organize_data_cohort_multimodal_genes_signals(
    data_signals_genes_persons=None,
    report=None,
):
    """
    Organizes genes' pan-tissue signals across persons for principal component
    analysis.

    Data has persons across rows and genes across columns.

             gene_1 gene_2 gene_3 gene_4 gene_5
    person_1 ...    ...    ...    ...    ...
    person_2 ...    ...    ...    ...    ...
    person_3 ...    ...    ...    ...    ...
    person_4 ...    ...    ...    ...    ...
    person_5 ...    ...    ...    ...    ...

    1. normalize genes' pan-tissue signals by logarithmic transformation
    2. standardize genes' pan-tissue signals by z-score transformation

    arguments:
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Copy data.
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Normalize genes' pan-tissue signals across persons.
    # Transform genes' signals to logarithmic space.
    # Shift values arithmetically to avoid negative values.
    data_normal = utility.calculate_pseudo_logarithm_signals_negative(
        pseudo_count=1.0,
        base=2,
        axis="index", # apply function to each column of data frame
        data=data_signals,
    )
    # Standardize genes' pan-tissue signals across persons.
    # This method inserts missing values if the standard deviation is zero.
    data_standard = data_normal.apply(
        lambda series: scipy.stats.zscore(
            series.to_numpy(),
            axis=0,
            ddof=1, # Sample standard deviation.
            nan_policy="omit",
        ),
        axis="index", # apply function to each column of data frame
    )

    # Compare summary statistics before and after transformation.
    utility.print_terminal_partition(level=3)
    data_mean = data_standard.aggregate(
        lambda x: x.mean(),
        axis="index", # apply function to each column of data frame
    )
    print("Mean")
    print(data_mean.iloc[0:10])
    data_deviation = data_standard.aggregate(
        lambda x: x.std(),
        axis="index", # apply function to each column of data frame
    )
    print("Standard deviation")
    print(data_deviation.iloc[0:10])
    utility.print_terminal_partition(level=3)

    if False:

        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("Selection of genes with pan-tissue signals across persons.")
            utility.print_terminal_partition(level=3)
            print(data_selection)
        # Return information.
        return result
    pass



def organize_multimodal_genes_signals_persons_components(
    genes=None,
    data_signals_genes_persons=None,
    report=None,
):
    """
    Organizes a principal components analysis on genes' pan-tissue signals as
    features across persons as instances.

    arguments:
        genes (list<str>): identifiers of genes
        data_signals_genes_persons (object): Pandas data frame of genes'
            pan-tissue signals across persons
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames of genes' pairwise
            correlations

    """

    # Copy data.
    data_signals = data_signals_genes_persons.copy(deep=True)
    # Select genes of interest.
    data_selection = data_signals.loc[
        :, data_signals.columns.isin(genes)
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Selection of genes with pan-tissue signals across persons.")
        utility.print_terminal_partition(level=3)
        print(data_selection)
    # Reduce dimensionality.
    components = min(int(len(genes)), int(data_selection.shape[0]))
    result = utility.calculate_principal_components(
        data=data_selection,
        components=components,
        report=report,
    )
    # Return information.
    return result


def write_product_cohort_gene_components(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_data_components = os.path.join(
        paths[cohort]["component"], "data_persons_genes_components.pickle"
    )
    path_data_variances = os.path.join(
        paths[cohort]["component"], "data_persons_genes_variances.pickle"
    )
    # Write information to file.
    pandas.to_pickle(
        information["data_persons_genes_components"],
        path_data_components
    )
    pandas.to_pickle(
        information["data_persons_genes_variances"],
        path_data_variances
    )
    pass


# TODO:
# normalize signals for multimodal genes and save them to file in order to plot them later
#


def organize_cohort_gene_components(
    cohort=None,
    paths=None,
    report=None,
):
    """
    Organizes evaluation of subpopulation structure on the basis of pan-tissue
    expression of genes of interest.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("cohort: " + cohort)

    # Read source information from file.
    source = read_source_cohort_gene_components(
        cohort=cohort,
        dock=paths["dock"],
    )
    # Organize data for principal component analysis.
    data_signals_genes_persons = organize_data_cohort_multimodal_genes_signals(
        data_signals_genes_persons=source["data_signals_genes_persons"],
        report=report,
    )

    if False:
        # Calculate principal components on genes across persons.
        bin = calculate_multimodal_genes_signals_persons_components(
            genes=source["genes_candidacy"]["multimodal"],
            data_signals_genes_persons=data_signals_genes_persons,
            report=report,
        )
        # Compile information.
        information = dict()
        information["data_persons_genes_components"] = bin[
            "data_observations_components"
        ]
        information["data_persons_genes_variances"] = bin[
            "data_components_variances"
        ]
        # Write information to file.
        write_product_cohort_gene_components(
            cohort=cohort,
            information=information,
            paths=paths,
        )
    pass


##########
# Regressions on components


def read_source_cohort_components_regressions(
    cohort=None,
    dock=None
):
    """
    Reads and organizes source information from file

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:
        (object): source information
    """

    # Specify directories and files.
    path_data_persons_properties = os.path.join(
        dock, "selection", "tight", "persons_properties", cohort,
        "data_persons_properties.pickle"
    )
    path_data_persons_genes_components = os.path.join(
        dock, "stratification", cohort, "component",
        "data_persons_genes_components.pickle"
    )
    path_data_persons_genes_variances = os.path.join(
        dock, "stratification", cohort, "component",
        "data_persons_genes_variances.pickle"
    )
    # Read information from file.
    data_persons_properties = pandas.read_pickle(path_data_persons_properties)
    data_persons_genes_components = pandas.read_pickle(
        path_data_persons_genes_components
    )
    data_persons_genes_variances = pandas.read_pickle(
        path_data_persons_genes_variances
    )
    # Return information.
    return {
        "data_persons_properties": data_persons_properties,
        "data_persons_genes_components": data_persons_genes_components,
        "data_persons_genes_variances": data_persons_genes_variances,
    }


def organize_dependent_independent_variables(
    variables=None,
    threshold_variance=None,
    data_persons_properties=None,
    data_persons_genes_components=None,
    data_persons_genes_variances=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        variables (list<str>): names of independent variables for regression
        threshold_variance (float): minimal value of variance
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        data_persons_genes_components (object): Pandas data frame of components
            on multimodal genes across persons
        data_persons_genes_variances (object): Pandas data frame of components'
            variances

    raises:

    returns:
        (object): Pandas data frame of dependent and independent variables

    """

    # Copy data.
    data_properties = data_persons_properties.copy(deep=True)
    data_components = data_persons_genes_components.copy(
        deep=True
    )
    data_variances = data_persons_genes_variances.copy(
        deep=True
    )
    # Select components for regression.
    data_variances_threshold = data_variances.loc[
        data_variances["variance"] > threshold_variance, :
    ]
    components_indices = data_variances_threshold["component"].tolist()
    components = list()
    for index in components_indices:
        component = str("component_" + str(index))
        components.append(component)
        pass
    # Organize data.
    data_properties = data_properties.loc[
        :, data_properties.columns.isin(variables)
    ]
    data_components = data_components.loc[
        :, data_components.columns.isin(components)
    ]
    # Join signals on persons' properties.
    data_variables = data_properties.join(
        data_components,
        how="left",
        on="person"
    )
    # Compile information.
    bin = dict()
    bin["components"] = components
    bin["data_variables"] = data_variables
    # Return information.
    return bin


def regress_cases(
    cases=None,
    variables=None,
    data_variables=None,
    report=None,
):
    """
    Drives iterative regression across cases.

    arguments:
        cases (list<str>): names of independent cases for which to regress
            variables
        variables (list<str>): names of independent variables for regression
        data_variables (object): Pandas data frame of dependent and independent
            variables
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of residuals for regression on each case and data
            summary on regressions

    """

    # Regress.
    count_total = len(cases)
    records = list()
    residuals_cases = dict()
    counter = 0
    for case in cases:
        # Report.
        counter += 1
        if report:
            #print("count " + str(counter) + ": " + str(case))
            percentage = (counter / count_total) * 100
            if (percentage % 10 == 0):
                print("complete cases: " + str(round(percentage)) + "%")
        bin_case = prediction.regress_signal_ordinary_residuals(
            dependence=case,
            independence=variables,
            proportion=0.5,
            data=data_variables,
        )
        # Collect residuals.
        residuals_cases[case] = bin_case["residuals"]
        # Collect reports.
        report = bin_case["report"]
        report["case"] = case
        records.append(report)
        pass
    # Organize data.
    data_regression = utility.convert_records_to_dataframe(
        records=records
    )
    data_regression.set_index(
        ["case"],
        append=False,
        drop=True,
        inplace=True
    )
    columns = list()
    columns.append("observations")
    columns.append("freedom")
    columns.append("r_square")
    columns.append("r_square_adjust")
    columns.append("log_likelihood")
    columns.append("akaike")
    columns.append("bayes")
    columns.append("condition")
    columns.append("intercept_parameter")
    columns.append("intercept_probability")
    for variable in variables:
        columns.append(str(variable + ("_parameter")))
        columns.append(str(variable + ("_probability")))
        columns.append(str(variable + ("_inflation")))
        pass
    data_regression = data_regression[[*columns]]
    # Compile information.
    bin = dict()
    bin["residuals_cases"] = residuals_cases
    bin["data_regression"] = data_regression
    # Return information.
    return bin


def organize_data_regress_cases_report(
    variables_regression=None,
    data_persons_properties=None,
    data_persons_genes_components=None,
    data_persons_genes_variances=None,
    threshold_discovery=None,
    discovery=None,
    report=None,
):
    """
    Organizes and combines information about dependent and independent
    variables for regression.

    arguments:
        variables_regression (list<str>): names of independent variables for
            regression
        data_persons_properties (object): Pandas data frame of persons and
            their relevant properties
        data_persons_genes_components (object): Pandas data frame of components
            on multimodal genes across persons
        data_persons_genes_variances (object): Pandas data frame of components'
            variances
        threshold_discovery (float): threshold by false discovery rate
        discovery (bool): whether to calculate false discovery rates across all
            genes
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Organize dependent and independent variables for regression analysis.
    bin_variables = organize_dependent_independent_variables(
        variables=variables_regression,
        threshold_variance=0.01,
        data_persons_properties=data_persons_properties,
        data_persons_genes_components=data_persons_genes_components,
        data_persons_genes_variances=data_persons_genes_variances,
    )
    # Regress each gene's signal across persons.
    # Iterate on genes.
    bin_regression = regress_cases(
        cases=bin_variables["components"],
        variables=variables_regression,
        data_variables=bin_variables["data_variables"],
        report=True,
    )
    # Summarize regression quality.
    # Report means of statistics across independent models.
    if report:
        prediction.report_regression_models_quality(
            variables=variables_regression,
            data_regression_models=bin_regression["data_regression"],
        )
        pass
    # Return information.
    return bin_regression


def write_product_cohort_components_regressions(
    cohort=None,
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Specify directories and files.
    path_residuals = os.path.join(
        paths[cohort]["regression"], "residuals_cases.pickle"
    )
    path_data_regression = os.path.join(
        paths[cohort]["regression"], "data_regression_components.pickle"
    )
    path_data_regression_text = os.path.join(
        paths[cohort]["regression"], "data_regression_components.tsv"
    )
    # Write information to file.
    with open(path_residuals, "wb") as file_product:
        pickle.dump(information["residuals_cases"], file_product)
    pandas.to_pickle(
        information["data_regression"],
        path_data_regression
    )
    information["data_regression"].to_csv(
        path_or_buf=path_data_regression_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


# TODO: need to introduce different models for each cohort... probably just "main" for each cohort


def organize_cohort_components_regressions(
    cohort=None,
    paths=None,
    report=None,
):
    """
    Organizes evaluation of subpopulation structure on the basis of pan-tissue
    expression of genes of interest.

    arguments:
        cohort (str): cohort of persons--selection, respiration, or ventilation
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("cohort: " + cohort)
    # Read source information from file.
    source = read_source_cohort_components_regressions(
        cohort=cohort,
        dock=paths["dock"],
    )
    # Define variables for regression models.
    variables = selection.define_variables()
    # Organize data and regress across components.
    bin_regression = organize_data_regress_cases_report(
        variables_regression=(
            variables[cohort]["model_hypothesis"]
        ),
        data_persons_properties=source["data_persons_properties"],
        data_persons_genes_components=source["data_persons_genes_components"],
        data_persons_genes_variances=source["data_persons_genes_variances"],
        threshold_discovery=0.05,
        discovery=False,
        report=True,
    )
    # Write information to file.
    write_product_cohort_components_regressions(
        cohort=cohort,
        information=bin_regression,
        paths=paths,
    )
    pass


###############################################################################
# Procedure


def execute_procedure(
    dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        dock (str): path to root or dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Initialize directories.
    paths = initialize_directories(dock=dock)

    # Organize data and regress cases.
    cohorts = [
        "selection",
        "respiration",
        "ventilation",
    ]
    for cohort in cohorts:
        organize_cohort_gene_components(
            cohort=cohort,
            paths=paths,
            report=True,
        )
        if False:
            organize_cohort_components_regressions(
                cohort=cohort,
                paths=paths,
                report=True,
            )
    pass


if (__name__ == "__main__"):
    execute_procedure()
